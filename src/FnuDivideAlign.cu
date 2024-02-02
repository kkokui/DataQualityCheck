#include "FnuDivideAlign.h"

#include <stdio.h>
#include <EdbPattern.h>
#include <TFile.h>

#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <thrust/sort.h>
#include <thrust/execution_policy.h>

struct cudaSegment
{
	int flag, pid;
	float x, y, z;
};

struct cudaTrack
{
	float x, y, z, tx, ty, tx_first8, ty_first8;
	int nseg;
	cudaSegment segments[NPIDMAX];
};

float robustFactor = 1.0;
int ncall;
TObjArray *gTracks;

// Data buffer for the GPU process
double *h_params;
double *d_params;
cudaTrack *h_tracks;
cudaTrack *d_tracks;
float *h_chi2;
float *d_chi2;

FnuDivideAlign::FnuDivideAlign()
	: binWidth(2000), rangeXY(8500)
{
}

FnuDivideAlign::~FnuDivideAlign()
{
}

void FnuDivideAlign::SetBinWidth(double bwidth)
{
	binWidth = bwidth;
}

void FnuDivideAlign::SetRobustFactor(float rfactor)
{
	robustFactor = rfactor;
}

double FnuDivideAlign::GetBinWidth()
{
	return binWidth;
}

float FnuDivideAlign::GetRobustFactor()
{
	return robustFactor;
}

__global__ void lsm_kernel(int n, cudaTrack *d_trk, double *d_param)
{
	// Calculation of least square method

	// access thread id
	const unsigned int tid = threadIdx.x;
	const unsigned int tsize = blockDim.x;
	const unsigned int bid = blockIdx.x;
	// x = a0 + a1*z
	int pos = tid + tsize * bid;
	if (pos < n)
	{
		cudaTrack *t = &d_trk[pos];
		int i;
		double A00 = 0, A01 = 0, A02 = 0, A11 = 0, A12 = 0;
		double B00 = 0, B01 = 0, B02 = 0, B11 = 0, B12 = 0;
		for (i = 0; i < NPIDMAX; i++)
		{
			cudaSegment *s = &t->segments[i];
			float x = s->x + d_param[s->pid * 2];
			float y = s->y + d_param[s->pid * 2 + 1];
			float z = s->z;
			if (s->flag)
			{
				A00 += 1.0;
				A01 += z;
				A02 += x;
				A11 += z * z;
				A12 += z * x;
				B00 += 1.0;
				B01 += z;
				B02 += y;
				B11 += z * z;
				B12 += z * y;
			}
		}

		t->x = (A02 * A11 - A01 * A12) / (A00 * A11 - A01 * A01);
		t->tx = (A00 * A12 - A01 * A02) / (A00 * A11 - A01 * A01);

		t->y = (B02 * B11 - B01 * B12) / (B00 * B11 - B01 * B01);
		t->ty = (B00 * B12 - B01 * B02) / (B00 * B11 - B01 * B01);
		t->z = 0;
	}
	__syncthreads();
}

__global__ void calc_chi2_kernel(int n, cudaTrack *d_trk, double *p, float *d_chi2)
{
	// Calculate chi2 of a track

	// access thread id
	const unsigned int tid = threadIdx.x;
	const unsigned int tsize = blockDim.x;
	const unsigned int bid = blockIdx.x;

	int pos = tid + tsize * bid;
	float chi2 = 0;
	if (pos < n)
	{
		float sigmaPos2 = 0.36; // 0.6 * 0.6;
		float sigmaAng2 = 4e-6; // 0.002 * 0.002;

		cudaTrack *t = &d_trk[pos];

		int nseg = 0;
		for (int i = 0; i < NPIDMAX; i++)
		{
			cudaSegment *s = &t->segments[i];
			if (s->flag == 0)
				continue;
			float x = s->x + p[i * 2];
			float y = s->y + p[i * 2 + 1];
			float z = s->z;
			float dx = x - (t->x + t->tx * (z - t->z));
			float dy = y - (t->y + t->ty * (z - t->z));
			chi2 += dx * dx + dy * dy;
			nseg++;
		}
		chi2 /= sigmaPos2 * nseg;
		float txFit = t->tx;
		float tyFit = t->ty;
		float tx = t->tx_first8;
		float ty = t->ty_first8;
		chi2 += ((tx - txFit) * (tx - txFit) + (ty - tyFit) * (ty - tyFit)) / sigmaAng2;
		d_chi2[pos] = chi2;
	}
	__syncthreads();
}

void fitfuncRobust(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *p, Int_t iflag)
{
	// Fit function for TMinuit

	double delta2 = 0.0;

	for (int i = 0; i < NPIDMAX * 2; i++)
	{
		h_params[i] = p[i];
	}

	checkCudaErrors(cudaMemcpy(d_params, h_params, sizeof(double) * NPIDMAX * 2, cudaMemcpyHostToDevice));

	int ntrk = gTracks->GetEntriesFast();
	int numthread = 512;
	int numblock = (ntrk + numthread - 1) / numthread;
	dim3 threads(numthread, 1, 1);
	dim3 blocks(numblock, 1, 1);

	lsm_kernel<<<blocks, threads>>>(ntrk, d_tracks, d_params);
	cudaDeviceSynchronize();
	// check if kernel execution generated and error
	getLastCudaError("lsm Kernel execution failed");

	calc_chi2_kernel<<<blocks, threads>>>(ntrk, d_tracks, d_params, d_chi2);
	cudaDeviceSynchronize();
	getLastCudaError("calc chi2 Kernel execution failed");

	// thrust, sort on GPU
	thrust::sort(thrust::device, d_chi2, d_chi2 + ntrk);
	checkCudaErrors(cudaMemcpy(h_chi2, d_chi2, sizeof(float) * ntrk, cudaMemcpyDeviceToHost));
	int nrobust = ntrk * robustFactor; // For example, if robustFactor = 0.5, only 50% of tracks will be used.
	for (int i = 0; i < nrobust; i++)
	{
		delta2 += h_chi2[i];
	}
	// delta2/=nrobust;

	// regularization
	//  double lambda = 0.1;
	//  for(int i=0;i<nPID*2)
	//  {
	//  	delta2+=lambda*p[i]*p[i]; //L2 regularization
	//  	// delta2+=lambda*abs(p[i]); //L1 regularization
	//  }

	fval = delta2;
	if (ncall % 1000 == 0)
	{
		printf("ncall=%d fval = %lf ", ncall, fval);
		for (int i = 0; i < 10; i++)
		{
			printf("%4.1lf ", p[i]);
		}
		printf("\n");
	}
	ncall++;
}

void FnuDivideAlign::CalcAlignPar(TObjArray *tracks, double iX, double iY, int fixflag)
{
	// Calculate alignment parameters in a divided area

	gTracks = tracks;
	int ntrk = gTracks->GetEntriesFast();
	// Cuda data buffers
	checkCudaErrors(cudaMallocHost((void **)&h_tracks, sizeof(cudaTrack) * ntrk));
	checkCudaErrors(cudaMallocHost((void **)&h_chi2, sizeof(float) * ntrk));
	checkCudaErrors(cudaMallocHost((void **)&h_params, sizeof(double) * NPIDMAX * 2));

	checkCudaErrors(cudaMalloc((void **)&d_tracks, sizeof(cudaTrack) * ntrk));
	checkCudaErrors(cudaMalloc((void **)&d_chi2, sizeof(float) * ntrk));
	checkCudaErrors(cudaMalloc((void **)&d_params, sizeof(double) * NPIDMAX * 2));

	for (int i = 0; i < ntrk; i++)
	{
		// Setup structures for tracks
		EdbTrackP *t = (EdbTrackP *)gTracks->At(i);
		cudaTrack *ct = &h_tracks[i];
		ct->tx_first8 = t->TX();
		ct->ty_first8 = t->TY();
		ct->nseg = t->N();
		for (int ipid = 0; ipid < NPIDMAX; ipid++)
		{
			ct->segments[ipid].flag = 0;
		} // Clear initial values
		for (int iseg = 0; iseg < t->N(); iseg++)
		{
			EdbSegP *s = t->GetSegment(iseg);
			if (fabs(s->X() - iX) < binWidth / 2 && fabs(s->Y() - iY) < binWidth / 2)
			{
				int pid = s->PID();
				ct->segments[pid].flag = 1;
				ct->segments[pid].x = s->X();
				ct->segments[pid].y = s->Y();
				ct->segments[pid].z = s->Z();
			}
		}
	}
	checkCudaErrors(cudaMemcpy(d_tracks, h_tracks, sizeof(cudaTrack) * ntrk, cudaMemcpyHostToDevice));
	// The default minimizer is Minuit, you can also try Minuit2
	TVirtualFitter::SetDefaultFitter("Minuit");
	// minuit->BuildArrays(30);
	// Int_t SetParameter(Int_t ipar, const char* parname, Double_t value, Double_t verr, Double_t vlow, Double_t vhigh)

	// Set parameters of the initial plate
	int pid = 0; // initial plate
	minuit->SetParameter(2 * pid, Form("dx%d", pid), 0, 0, 0, 0);
	minuit->SetParameter(2 * pid + 1, Form("dy%d", pid), 0, 0, 0, 0);
	// Set parameters of the middle plate
	for (pid = 1; pid < nPID - 1; pid++)
	{
		minuit->SetParameter(2 * pid, Form("dx%d", pid), 0, 0.1, 0, 0);
		minuit->SetParameter(2 * pid + 1, Form("dy%d", pid), 0, 0.1, 0, 0);
	}
	// Set parameters of the last plate
	pid = nPID - 1; // last plate
	if (fixflag == 1)
	{
		minuit->SetParameter(2 * pid, Form("dx%d", pid), 0, 0, 0,0);
		minuit->SetParameter(2 * pid + 1, Form("dy%d", pid), 0, 0, 0,0);
	}
	else
	{
		minuit->SetParameter(2 * pid, Form("dx%d", pid), 0, 0.1, 0, 0);
		minuit->SetParameter(2 * pid + 1, Form("dy%d", pid), 0, 0.1, 0, 0);
	}

	minuit->SetFCN(fitfuncRobust);

	double arglist[200];
	arglist[0] = 0;
	// set print level. arglist[0]==0 is minimum print.
	minuit->ExecuteCommand("SET PRIntout", arglist, 1);

	// minimize
	arglist[0] = 100000; // number of function calls
	arglist[1] = 50;	 // tolerance of estimated verical distance to minimum
	// minuit->SetMaxIterations(10000);
	ncall = 0;
	printf("Aligning iX = %.0f, iY = %.0f, ntrk = %d, robustFactor = %.1f\n", iX, iY, ntrk, robustFactor);
	minuit->ExecuteCommand("MIGRAD2", arglist, 2);

	// get result
	for (int i = 0; i < nPID * 2; ++i)
	{
		p[i] = minuit->GetParameter(i);
		// parErrors[i] = minuit->GetParError(i);
	}

	checkCudaErrors(cudaFreeHost(h_tracks));
	checkCudaErrors(cudaFreeHost(h_chi2));
	checkCudaErrors(cudaFreeHost(h_params));

	checkCudaErrors(cudaFree(d_tracks));
	checkCudaErrors(cudaFree(d_chi2));
	checkCudaErrors(cudaFree(d_params));
}

int FnuDivideAlign::CountPassedSeg(EdbTrackP *t, double iX, double iY)
{
	// Count a number of segments in one track which passed a divided area
	int count = 0;
	for (int iseg = 0; iseg < t->N(); iseg++)
	{
		EdbSegP *s = t->GetSegment(iseg);
		if (fabs(s->X() - iX) < binWidth / 2 && fabs(s->Y() - iY) < binWidth / 2)
			count++;
	}
	return count;
}

void FnuDivideAlign::ApplyAlign(EdbTrackP *t, double iX, double iY)
{
	// Apply alignment to segments which passed a divided area
	for (int iseg = 0; iseg < t->N(); iseg++)
	{
		EdbSegP *s = t->GetSegment(iseg);
		if (fabs(s->X() - iX) < binWidth / 2 && fabs(s->Y() - iY) < binWidth / 2)
		{
			int pid = s->PID();
			s->SetX(s->X() + p[pid * 2]);
			s->SetY(s->Y() + p[pid * 2 + 1]);
		}
	}
}
void FnuDivideAlign::ApplyAlignBicubic(EdbSegP *s,double Xcenter,double Ycenter)
{
	// Apply alignment to segments which passed a divided area

	int pid = s->PID();

	// get shift value calculated by divide align
	// double 

	// get x and y values of the segment
	double segmentPositionX = s->X();
	double segmentPositionY = s->Y();
	// detect nearest middle position of division
	int nearestMiddlePointNumberX = (segmentPositionX - (Xcenter - rangeXY) - binWidth/2) / binWidth;
	int nearestMiddlePointNumberY = (segmentPositionY - (Ycenter - rangeXY) - binWidth/2) / binWidth;
	double nearestMiddlePointPositionX = Xcenter - rangeXY + nearestMiddlePointNumberX*binWidth;
	double nearestMiddlePointPositionY = Ycenter - rangeXY + nearestMiddlePointNumberY*binWidth;
	// detect 16 reference values
	double referencePositionX[4][4];
	double referencePositionY[4][4];

	double referenceShiftX[4][4];
	double referenceShiftY[4][4];
	// if(segmentPositionX<nearestBinPositionX)
	// {

	// 	referenceShiftX[0][0] = 
	// }
	// calculate alignment parameter
	// apply alignment

}
int FnuDivideAlign::Align(TObjArray *tracks, double Xcenter, double Ycenter, int nPatterns)
{
	// Divide the area, Calculate alignment parameters and apply alignment
	nPID = nPatterns;
	minuit = TVirtualFitter::Fitter(0, nPID * 2);
	alignPar = new TTree("alignPar", "alignPar");
	alignPar->Branch("iX", &iXBranchValue);
	alignPar->Branch("iY", &iYBranchValue);
	alignPar->Branch("shiftX", &shiftXBranchValue);
	alignPar->Branch("shiftY", &shiftYBranchValue);
	alignPar->Branch("pid", &pidBranchValue);

	// int nDivisionX = rangeXY*2/binWidth +1;
	// int nDivisionY = rangeXY*2/binWidth +1;
	// double ***shiftXEachDivision = new double**[nPatterns];
	// double ***shiftYEachDivision = new double**[nPatterns];
	// for(int pid=0;pid<nPatterns;pid++)
	// {
	// 	shiftXEachDivision[pid] = new double*[nDivisionX];
	// 	shiftYEachDivision[pid] = new double*[nDivisionX];
	// 	for(int iDivisionX=0;iDivisionX<nDivisionX;iDivisionX++)
	// 	{
	// 		shiftXEachDivision[pid][iDivisionX] = new double[nDivisionY];
	// 		shiftYEachDivision[pid][iDivisionX] = new double[nDivisionY];
	// 	}
	// }

	int ntrk = tracks->GetEntriesFast();

	double angleXSum = 0;
	double angleYSum = 0;
	for (int itrk = 0; itrk < ntrk; itrk++)
	{
		EdbTrackP *t = (EdbTrackP*)tracks->At(itrk);
		angleXSum += t->TX();
		angleYSum += t->TY();
	}
	double angleXMean = angleXSum/ntrk;
	double angleYMean = angleYSum/ntrk;

	// Divide the area into binWidth*binWidth um^2 areas
	for (iYBranchValue = Ycenter - rangeXY + binWidth / 2; iYBranchValue <= Ycenter + rangeXY; iYBranchValue += binWidth)
	// max of iY is not correct sometimes...
	{
		for (iXBranchValue = Xcenter - rangeXY + binWidth / 2; iXBranchValue <= Xcenter + rangeXY; iXBranchValue += binWidth)
		{
			TObjArray *tracks2 = new TObjArray;

			for (int itrk = 0; itrk < ntrk; itrk++)
			{
				EdbTrackP *t = (EdbTrackP *)tracks->At(itrk);
				if (t->N() < 10 || abs(t->TX() -angleXMean) >= 0.01 || abs(t->TY() - angleYMean) >= 0.01)
				// if (t->N() < 10 )
				{
					continue;
				}
				if (10 <= CountPassedSeg(t, iXBranchValue, iYBranchValue)) //  check if the track passes the area
				{
					tracks2->Add(t);
				}
			}
			if (tracks2->GetEntries() < 20)
			{
				continue;
			}

			// calculate the alignment parameters several times.
			for (int j = 0; j < 1; j++)
			{
				CalcAlignPar(tracks2, iXBranchValue, iYBranchValue, 0);
			}

			int iXID = (iXBranchValue - (Xcenter - rangeXY + binWidth / 2))/binWidth;
			int iYID = (iYBranchValue - (Ycenter - rangeXY + binWidth / 2))/binWidth;
			for (pidBranchValue = 0; pidBranchValue < nPID; pidBranchValue++)
			{
				shiftXBranchValue = p[pidBranchValue * 2];
				shiftYBranchValue = p[pidBranchValue * 2 + 1];
				alignPar->Fill();
			}
			for (int itrk = 0; itrk < ntrk; itrk++)
			{
				EdbTrackP *t = (EdbTrackP *)tracks->At(itrk);
				ApplyAlign(t, iXBranchValue, iYBranchValue);
			}
			delete tracks2;
		}
	}
	return 0;
}

void FnuDivideAlign::WriteAlignPar(TString filename)
{
	// Write TTree for Shifts of alignment
	TFile fout1(filename, "recreate");
	alignPar->Write();
	fout1.Close();
}
