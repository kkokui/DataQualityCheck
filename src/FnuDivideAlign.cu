#include "FnuDivideAlign.h"

#include <EdbDataSet.h>
#include <TObjArray.h>
#include <TVirtualFitter.h>
#include <TMath.h>

#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <thrust/sort.h>
#include <thrust/execution_policy.h>


float robustFactor;
int ncall = 0;
TObjArray *gTracks;

double *h_params;
double *d_params;
cudaTrack* d_tracks;
float* d_chi2;
float* h_chi2;

FnuDivideAlign::FnuDivideAlign()
{
	binWidth = 2000;
	robustFactor = 1.0;
	XYrange = 8500;
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

// 最小二乗法の計算 


__global__ void lsm_kernel(int n, cudaTrack* d_trk, double *d_param)
{
	// access thread id
	const unsigned int tid = threadIdx.x;
	const unsigned int tsize = blockDim.x;
	const unsigned int bid = blockIdx.x;
	// y = a0 + a1*x
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
	
	/*if(ncall==0){
		checkCudaErrors( cudaMemcpy( h_tracks, d_tracks, sizeof(cudaTrack)*ntrk, cudaMemcpyDeviceToHost) );
		printf("h_track %f %f %f %f %f\n", h_tracks[0].x, h_tracks[0].y, h_tracks[0].z, h_tracks[0].tx,  h_tracks[0].ty);
	}*/
	
	calc_chi2_kernel<<< blocks, threads >>> (ntrk, d_tracks, d_params, d_chi2);
	cudaDeviceSynchronize();
	getLastCudaError("calc chi2 Kernel execution failed");
	
	// thrust, sort on GPU
	thrust::sort(thrust::device, d_chi2, d_chi2+ntrk);
	checkCudaErrors( cudaMemcpy( h_chi2, d_chi2, sizeof(float)*ntrk, cudaMemcpyDeviceToHost) );
	// float robustFactor = 0.5; // only % of segments will be used.
	int nrobust = ntrk*robustFactor;
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

void FnuDivideAlign::calc_align_par(TObjArray *tracks,double iX, double iY, int fixflag)
{
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
		EdbTrackP *t = (EdbTrackP *)gTracks->At(i);
		cudaTrack *ct = &h_tracks[i];
		ct->tx_first8 = t->TX();
		ct->ty_first8 = t->TY();
		ct->nseg = t->N();
		for(int ipid=0; ipid<NPIDMAX; ipid++){ ct->segments[ipid].flag=0;} //初期値クリア
		for (int iseg = 0; iseg < t->N(); iseg++)
		{
			EdbSegP *s = t->GetSegment(iseg);
			// if (1)
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
	// gMinuit->BuildArrays(30);
	// Int_t SetParameter(Int_t ipar, const char* parname, Double_t value, Double_t verr, Double_t vlow, Double_t vhigh)

	int pid=0; //最初のプレート
	gMinuit->SetParameter(2*pid,Form("dx%d",pid),0,0,0,0);
	gMinuit->SetParameter(2*pid+1,Form("dy%d",pid),0,0,0,0);
	for(pid=1;pid<nPID-1;pid++)
	{
		gMinuit->SetParameter(2*pid,Form("dx%d",pid),0,0.1,-30,30);
		gMinuit->SetParameter(2*pid+1,Form("dy%d",pid),0,0.1,-30,30);
	}
	pid=nPID-1; // 最後のプレート
	if(fixflag==1)
	{
		gMinuit->SetParameter(2*pid,Form("dx%d",pid),0,0,0,0);
		gMinuit->SetParameter(2*pid+1,Form("dy%d",pid),0,0,0,0);
	}else{
		gMinuit->SetParameter(2*pid,Form("dx%d",pid),0,0.1,-30,30);
		gMinuit->SetParameter(2*pid+1,Form("dy%d",pid),0,0.1,-30,30);
	}
	
	
	gMinuit->SetFCN(fitfuncRobust);
	
	
	double arglist[200];

	arglist[0] = 0;
	// set print level. arglist[0]==0 is minimum print.
	gMinuit->ExecuteCommand("SET PRIntout",arglist,1);

	// minimize
	arglist[0] = 200000; // number of function calls
	arglist[1] = 0.001; // tolerance
	gMinuit->SetMaxIterations( 10000 );
	ncall =0;
	printf("Aligning iX = %.0f, iY = %.0f, ntrk = %d, robustFactor = %.1f\n",iX,iY,ntrk,robustFactor);
	gMinuit->ExecuteCommand("MIGRAD2", arglist, 2);
	/*
	gMinuit->SetFCN(fitfuncRobust3);
	gMinuit->ExecuteCommand("MIGRAD",arglist,2);
	*/
	
	/*
	double p[3];
	double parErrors[3];
	*/
	// get result
	for (int i = 0; i < nPID*2; ++i)
	{
		p[i] = gMinuit->GetParameter(i);

		// parErrors[i] = minuit->GetParError(i);
	}
	
	
	checkCudaErrors( cudaFreeHost( h_tracks) );
	checkCudaErrors( cudaFreeHost( h_chi2) );
	checkCudaErrors( cudaFreeHost( h_params) );

	checkCudaErrors( cudaFree( d_tracks) );
	checkCudaErrors( cudaFree( d_chi2) );
	checkCudaErrors( cudaFree( d_params) );
	

}


int FnuDivideAlign::count_passed_seg(EdbTrackP *t, double iX, double iY)
{
	int count = 0;
	for (int iseg = 0; iseg < t->N(); iseg++)
	{
		EdbSegP *s = t->GetSegment(iseg);
		if (fabs(s->X() - iX) < binWidth / 2 && fabs(s->Y() - iY) < binWidth / 2)
			count++;
	}
	return count;
}
void FnuDivideAlign::apply_align(EdbTrackP *t, double iX, double iY)
{
	for (int iseg = 0; iseg < t->N(); iseg++)
	{
		EdbSegP *s = t->GetSegment(iseg);
		// if (1)
		if (fabs(s->X() - iX) < binWidth / 2 && fabs(s->Y() - iY) < binWidth / 2)
		{
			int pid = s->PID();
			s->SetX(s->X() + p[pid * 2]);
			s->SetY(s->Y() + p[pid * 2 + 1]);
		}
	}
}

int FnuDivideAlign::dedicated_align(EdbPVRec *pvr,double Xcenter, double Ycenter)
{
	TObjArray *tracks = pvr->GetTracks();
	nPID = pvr->Npatterns();
	int plMin = pvr->GetPattern(0)->Plate();
	int plMax = pvr->GetPattern(nPID-1)->Plate();
	gMinuit = TVirtualFitter::Fitter(0, 300);
	alignPar = new TTree("alignPar","alignPar");
	double iX, iY, shiftX, shiftY;
	int pid;
	alignPar->Branch("iX",&iX);
	alignPar->Branch("iY",&iY);
	alignPar->Branch("shiftX",&shiftX);
	alignPar->Branch("shiftY",&shiftY);
	alignPar->Branch("pid",&pid);
	int ntrk = tracks->GetEntriesFast();

	for (iY = Ycenter - XYrange + binWidth / 2; iY <= Ycenter + XYrange; iY += binWidth) // Divide the area into 2*2 mm^2 areas
	{
		for (iX = Xcenter - XYrange + binWidth / 2; iX <= Xcenter + XYrange; iX += binWidth)
		{
			TObjArray *tracks2 = new TObjArray;

			for (int itrk = 0; itrk < ntrk; itrk++)
			{
				EdbTrackP *t = (EdbTrackP *)tracks->At(itrk);
				// if (t->N()<10|| abs(t->TX() + 0.01) >= 0.01 || abs(t->TY()-0.004) >= 0.01||t->GetSegment(0)->PID()>=10)
				if (t->N() < 10 || abs(t->TX() + 0.01) >= 0.01 || abs(t->TY() - 0.004) >= 0.01)
					continue;
				// if (fabs(t->X() - iX) < binWidth / 2 && fabs(t->Y() - iY) < binWidth / 2)
				if (10 <= count_passed_seg(t, iX, iY)) //  check if the track passes the area
					tracks2->Add(t);
			}
			// printf("iX = %.0f, iY = %.0f, ntrk = %d\n", iX, iY, tracks2->GetEntries());
			if (tracks2->GetEntries() == 0)
				continue;
			// calculate the alignment parameters several times.
			for (int j = 0; j < 1; j++)
			{
				calc_align_par(tracks2,iX,iY,0); //4th is fixflag
			}
			
			for(pid=0;pid<nPID;pid++)
			{
				shiftX = p[pid*2];
				shiftY = p[pid*2+1];
				alignPar->Fill();
			}
			// Apply alignment parameter.
			for (int itrk = 0; itrk < ntrk; itrk++)
			{
				EdbTrackP *t = (EdbTrackP *)tracks->At(itrk);
				// if (fabs(t->X() - iX) < binWidth / 2 && fabs(t->Y() - iY) < binWidth / 2)
				apply_align(t, iX, iY);
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
