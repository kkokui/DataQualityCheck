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

float robustFactor;
int ncall;
TObjArray *gTracks;

// Data buffer for the GPU process
double *h_params;
double *d_params;
cudaTrack *h_tracks;
cudaTrack* d_tracks;
float* h_chi2;
float* d_chi2;

FnuDivideAlign::FnuDivideAlign()
{
	binWidth = 2000;
	robustFactor = 1.0;
	rangeXY = 8500;
}

FnuDivideAlign::~FnuDivideAlign()
{}

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

__global__ void lsm_kernel(int n, cudaTrack* d_trk, double *d_param)
{
	// Calculation of least square method

	// access thread id
	const unsigned int tid = threadIdx.x;
	const unsigned int tsize = blockDim.x;
	const unsigned int bid = blockIdx.x;
	// x = a0 + a1*z
	int pos = tid + tsize*bid;
	if (pos < n) {
		cudaTrack *t = &d_trk[pos];
        int i;
        double A00=0 ,A01=0, A02=0, A11=0, A12=0;
        double B00=0 ,B01=0, B02=0, B11=0, B12=0;
        for (i=0;i<NPIDMAX;i++) {
			cudaSegment *s = &t->segments[i];
	 		float x = s->x + d_param[s->pid*2];
	 		float y = s->y + d_param[s->pid*2+1];
	 		float z = s->z;
			if(s->flag){
                A00+=1.0;
                A01+=z;
                A02+=x;
                A11+=z*z;
                A12+=z*x;
                B00+=1.0;
                B01+=z;
                B02+=y;
                B11+=z*z;
                B12+=z*y;
            }
        }
 
        t->x = (A02*A11-A01*A12) / (A00*A11-A01*A01);
        t->tx = (A00*A12-A01*A02) / (A00*A11-A01*A01);

        t->y = (B02*B11-B01*B12) / (B00*B11-B01*B01);
        t->ty = (B00*B12-B01*B02) / (B00*B11-B01*B01);
        t->z = 0;
	}
	__syncthreads();
}

__global__ void calc_chi2_kernel(int n, cudaTrack* d_trk, double *p, float *d_chi2)
{
	// Calculate chi2 of a track

	// access thread id
	const unsigned int tid = threadIdx.x;
	const unsigned int tsize = blockDim.x;
	const unsigned int bid = blockIdx.x;
	
	int pos = tid + tsize*bid;
	float chi2 = 0;
	if (pos < n) {
		float sigmaPos2 = 0.36;//0.6 * 0.6;
		float sigmaAng2 = 4e-6;//0.002 * 0.002;
		
		cudaTrack *t = &d_trk[pos];
		
		int nseg=0;
		for(int i=0; i<NPIDMAX; i++)
		{
			cudaSegment *s = &t->segments[i];
			if(s->flag==0) continue;
			float x = s->x + p[i*2];
			float y = s->y + p[i*2+1];
			float z = s->z;
			float dx = x - (t->x+t->tx*(z-t->z));
			float dy = y - (t->y+t->ty*(z-t->z));
			chi2 += dx * dx + dy * dy;
			nseg++;
		}
		chi2 /= sigmaPos2*nseg;
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
	
	for(int i=0; i<NPIDMAX*2; i++){h_params[i]=p[i];}
	
	checkCudaErrors( cudaMemcpy( d_params, h_params, sizeof(double)*NPIDMAX*2, cudaMemcpyHostToDevice) );
	
	int ntrk = gTracks->GetEntriesFast();
	int numthread = 512;
	int numblock = (ntrk + numthread -1)/numthread;
	dim3 threads(numthread, 1, 1);
	dim3 blocks(numblock, 1, 1);
	
	lsm_kernel <<< blocks, threads >>> (ntrk, d_tracks, d_params);
	cudaDeviceSynchronize();
	// check if kernel execution generated and error
	getLastCudaError("lsm Kernel execution failed");
	
	calc_chi2_kernel<<< blocks, threads >>> (ntrk, d_tracks, d_params, d_chi2);
	cudaDeviceSynchronize();
	getLastCudaError("calc chi2 Kernel execution failed");
	
	// thrust, sort on GPU
	thrust::sort(thrust::device, d_chi2, d_chi2+ntrk);
	checkCudaErrors( cudaMemcpy( h_chi2, d_chi2, sizeof(float)*ntrk, cudaMemcpyDeviceToHost) );
	int nrobust = ntrk * robustFactor; // For example, if robustFactor = 0.5, only 50% of tracks will be used.
	for (int i = 0; i < nrobust; i++) 
	{
		delta2 += h_chi2[i];
	}
	
	//regularization
	// double lambda = 0.1;
	// for(int i=0;i<nPID*2)
	// {
	// 	delta2+=lambda*p[i]*p[i]; //L2 regularization
	// 	// delta2+=lambda*abs(p[i]); //L1 regularization
	// }

	fval = delta2;
	if(ncall%1000==0) {
		printf("ncall=%d fval = %lf ", ncall, fval);
		for(int i=0; i<10; i++){printf("%4.1lf ", p[i]);}
		printf("\n");
	}
	ncall++;
}

void FnuDivideAlign::CalcAlignPar(TObjArray *tracks,double iX, double iY, int fixflag)
{
	// Calculate alignment parameters in a divided area

	gTracks = tracks;
	int ntrk = gTracks->GetEntriesFast();
	// Cuda data buffers
	checkCudaErrors( cudaMallocHost( (void**) &h_tracks, sizeof(cudaTrack)*ntrk) );
	checkCudaErrors( cudaMallocHost( (void**) &h_chi2, sizeof(float)*ntrk) );
	checkCudaErrors( cudaMallocHost( (void**) &h_params, sizeof(double)*NPIDMAX*2) );
	
	checkCudaErrors( cudaMalloc( (void**) &d_tracks, sizeof(cudaTrack)*ntrk) );
	checkCudaErrors( cudaMalloc( (void**) &d_chi2, sizeof(float)*ntrk) );
	checkCudaErrors( cudaMalloc( (void**) &d_params, sizeof(double)*NPIDMAX*2) );
	
	for (int i = 0; i < ntrk; i++)
	{
		// Setup structures for tracks
		EdbTrackP *t = (EdbTrackP *)gTracks->At(i);
		cudaTrack *ct = &h_tracks[i];
		ct->tx_first8 = t->TX();
		ct->ty_first8 = t->TY();
		ct->nseg = t->N();
		for(int ipid=0; ipid<NPIDMAX; ipid++){ ct->segments[ipid].flag=0;} // Clear initial values
		for (int iseg = 0; iseg < t->N(); iseg++)
		{
			EdbSegP *s = t->GetSegment(iseg);
			if (fabs(s->X() - iX) < binWidth / 2 && fabs(s->Y() - iY) < binWidth / 2)
			{
				int pid = s->PID();
				ct->segments[pid].flag=1;
				ct->segments[pid].x=s->X();
				ct->segments[pid].y=s->Y();
				ct->segments[pid].z=s->Z();
			}
		}
	}
	checkCudaErrors( cudaMemcpy( d_tracks, h_tracks, sizeof(cudaTrack)*ntrk, cudaMemcpyHostToDevice) );
	// The default minimizer is Minuit, you can also try Minuit2
	TVirtualFitter::SetDefaultFitter("Minuit");
	// minuit->BuildArrays(30);
	// Int_t SetParameter(Int_t ipar, const char* parname, Double_t value, Double_t verr, Double_t vlow, Double_t vhigh)

	// Set parameters of the initial plate
	int pid=0; //initial plate
	minuit->SetParameter(2*pid,Form("dx%d",pid),0,0,0,0);
	minuit->SetParameter(2*pid+1,Form("dy%d",pid),0,0,0,0);
	// Set parameters of the middle plate
	for(pid=1;pid<nPID-1;pid++)
	{
		minuit->SetParameter(2*pid,Form("dx%d",pid),0,0.1,-30,30);
		minuit->SetParameter(2*pid+1,Form("dy%d",pid),0,0.1,-30,30);
	}
	// Set parameters of the last plate
	pid=nPID-1; // last plate
	if(fixflag==1)
	{
		minuit->SetParameter(2*pid,Form("dx%d",pid),0,0,0,0);
		minuit->SetParameter(2*pid+1,Form("dy%d",pid),0,0,0,0);
	}else{
		minuit->SetParameter(2*pid,Form("dx%d",pid),0,0.1,-30,30);
		minuit->SetParameter(2*pid+1,Form("dy%d",pid),0,0.1,-30,30);
	}
	
	minuit->SetFCN(fitfuncRobust);
	
	double arglist[200];
	arglist[0] = 0;
	// set print level. arglist[0]==0 is minimum print.
	minuit->ExecuteCommand("SET PRIntout",arglist,1);

	// minimize
	arglist[0] = 200000; // number of function calls
	arglist[1] = 0.001; // tolerance
	minuit->SetMaxIterations( 10000 );
	ncall =0;
	printf("Aligning iX = %.0f, iY = %.0f, ntrk = %d, robustFactor = %.1f\n",iX,iY,ntrk,robustFactor);
	minuit->ExecuteCommand("MIGRAD2", arglist, 2);
	
	// get result
	for (int i = 0; i < nPID*2; ++i)
	{
		p[i] = minuit->GetParameter(i);
		// parErrors[i] = minuit->GetParError(i);
	}
	
	checkCudaErrors( cudaFreeHost( h_tracks) );
	checkCudaErrors( cudaFreeHost( h_chi2) );
	checkCudaErrors( cudaFreeHost( h_params) );

	checkCudaErrors( cudaFree( d_tracks) );
	checkCudaErrors( cudaFree( d_chi2) );
	checkCudaErrors( cudaFree( d_params) );
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

int FnuDivideAlign::Align(TObjArray *tracks,double Xcenter, double Ycenter,int nPatterns)
{
	// Divide the area, Calculate alignment parameters and apply alignment
	nPID = nPatterns;
	minuit = TVirtualFitter::Fitter(0, nPID*2);
	alignPar = new TTree("alignPar","alignPar");
	double iX, iY, shiftX, shiftY;
	int pid;
	alignPar->Branch("iX",&iX);
	alignPar->Branch("iY",&iY);
	alignPar->Branch("shiftX",&shiftX);
	alignPar->Branch("shiftY",&shiftY);
	alignPar->Branch("pid",&pid);
	int ntrk = tracks->GetEntriesFast();

	// Divide the area into binWidth*binWidth mm^2 areas
	for (iY = Ycenter - rangeXY + binWidth / 2; iY <= Ycenter + rangeXY; iY += binWidth)
	{
		for (iX = Xcenter - rangeXY + binWidth / 2; iX <= Xcenter + rangeXY; iX += binWidth)
		{
			TObjArray *tracks2 = new TObjArray;

			for (int itrk = 0; itrk < ntrk; itrk++)
			{
				EdbTrackP *t = (EdbTrackP *)tracks->At(itrk);
				if (t->N() < 10 || abs(t->TX() + 0.01) >= 0.01 || abs(t->TY() - 0.004) >= 0.01)
				{
					continue;
				}
				if (10 <= CountPassedSeg(t, iX, iY)) //  check if the track passes the area
				{
					tracks2->Add(t);
				}
			}
			if (tracks2->GetEntries() == 0)
			{
				continue;
			}
			
			// calculate the alignment parameters several times.
			for (int j = 0; j < 1; j++)
			{
				CalcAlignPar(tracks2, iX, iY, 0);
			}
			
			for(pid=0;pid<nPID;pid++)
			{
				shiftX = p[pid*2];
				shiftY = p[pid*2+1];
				alignPar->Fill();
			}
			for (int itrk = 0; itrk < ntrk; itrk++)
			{
				EdbTrackP *t = (EdbTrackP *)tracks->At(itrk);
				ApplyAlign(t, iX, iY);
			}
			delete tracks2;
		}
	}
	return 0;
}

void FnuDivideAlign::WriteAlignPar(TString filename)
{
	//Write TTree for Shifts of alignment
	TFile fout1(filename, "recreate");
	alignPar->Write();
	fout1.Close();
}
