#include <stdio.h>

#include <EdbDataSet.h>
#include <TGraph.h>
#include <TH2.h>
#include <TNtupleD.h>
#include <TObjArray.h>
#include <EdbPattern.h>
#include <TFitter.h>
#include <TVirtualFitter.h>
#include <TMath.h>

#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <thrust/sort.h>
#include <thrust/execution_policy.h>

const int NPIDMAX=800;
TObjArray *gTracks;
TVirtualFitter *gMinuit;
int gpidMax;
int gpidMin;
int gpidcutntrk;
int gflag_chi2dis;
float robustfactor;
double p[NPIDMAX*2] ={};
int nPID;
int ncall = 0;
const double XYrange = 8500;

struct cudaSegment{
	int flag, pid;
	float x,y,z;
};


struct cudaTrack{
	float x,y,z,tx,ty,tx_first8,ty_first8;
	int nseg;
	cudaSegment segments[NPIDMAX];
};

// Data buffer for the GPU process
cudaTrack* h_tracks;
double *h_params;
float* h_chi2;
cudaTrack* d_tracks;
float* d_chi2;
double *d_params;
int* indexArray;





// 最小二乗法の計算 
__global__ void lsm_kernel(int n, cudaTrack* d_trk, double *d_param) {
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

void lsm(double x[],double y[], int N, double &a0, double &a1)
{
	// y = a0 + a1*x
    int i;
    double A00=0 ,A01=0, A02=0, A11=0, A12=0;

    for (i=0;i<N;i++) {
        A00+=1.0;
        A01+=x[i];
        A02+=y[i];
        A11+=x[i]*x[i];
        A12+=x[i]*y[i];
    }
    a0 = (A02*A11-A01*A12) / (A00*A11-A01*A01);
    a1 = (A00*A12-A01*A02) / (A00*A11-A01*A01);
}

	

__global__ void calc_chi2_kernel(int n, cudaTrack* d_trk, double *p, float *d_chi2) {
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
	// float robustfactor = 0.5; // only % of segments will be used.
	int nrobust = ntrk*robustfactor;
	for (int i = 0; i < nrobust; i++) 
	{
		delta2 += h_chi2[i];
	}
	double lambda = 0.1;
	
	//regularization
	/*for(int i=0;i<nPID*2)
	{
		delta2+=lambda*p[i]*p[i]; //L2 regularization
		// delta2+=lambda*abs(p[i]); //L1 regularization
	}*/
	fval = delta2;
	if(ncall%1000==0) {
		printf("ncall=%d fval = %lf ", ncall, fval);
		for(int i=0; i<10; i++){printf("%4.1lf ", p[i]);}
		printf("\n");
	}
	ncall++;
}


void calc_align_par(TObjArray *tracks,double iX, double iY, double bin_width, int fixflag, TNtupleD *sta)
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
			if (fabs(s->X() - iX) < bin_width / 2 && fabs(s->Y() - iY) < bin_width / 2)
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
	gpidcutntrk = 0;
	gflag_chi2dis = 0;
	gMinuit->SetMaxIterations( 10000 );
	ncall =0;
	printf("Aligning iX = %.0f, iY = %.0f, ntrk = %d, robustfactor = %.1f\n",iX,iY,ntrk,robustfactor);
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
	
	for(int pid=0;pid<nPID;pid++)
	{
		sta->Fill(iX,iY,p[pid*2],p[pid*2+1],pid);
	}
	
	checkCudaErrors( cudaFreeHost( h_tracks) );
	checkCudaErrors( cudaFreeHost( h_chi2) );
	checkCudaErrors( cudaFreeHost( h_params) );

	checkCudaErrors( cudaFree( d_tracks) );
	checkCudaErrors( cudaFree( d_chi2) );
	checkCudaErrors( cudaFree( d_params) );
	

}

bool pass(EdbTrackP *t, double iX, double iY,double bin_width)
{
	for (int iseg = 0; iseg < t->N(); iseg++)
	{
		EdbSegP *s = t->GetSegment(iseg);
		if (fabs(s->X() - iX) < bin_width / 2 && fabs(s->Y() - iY) < bin_width / 2)
		{
			return true;
		}
	}
	return false;
}
int count_passed_seg(EdbTrackP *t, double iX, double iY, double bin_width)
{
	int count = 0;
	for (int iseg = 0; iseg < t->N(); iseg++)
	{
		EdbSegP *s = t->GetSegment(iseg);
		if (fabs(s->X() - iX) < bin_width / 2 && fabs(s->Y() - iY) < bin_width / 2)
			count++;
	}
	return count;
}
void apply_align(EdbTrackP *t, double iX, double iY,double bin_width)
{
	for (int iseg = 0; iseg < t->N(); iseg++)
	{
		EdbSegP *s = t->GetSegment(iseg);
		// if (1)
		if (fabs(s->X() - iX) < bin_width / 2 && fabs(s->Y() - iY) < bin_width / 2)
		{
			int pid = s->PID();
			s->SetX(s->X() + p[pid * 2]);
			s->SetY(s->Y() + p[pid * 2 + 1]);
		}
	}
}
int calc_dxy(EdbPVRec *pvr,TTree *tree, int ntrk, double Xcenter, double Ycenter, double bin_width)
{
	
	double deltaX, deltaY, tx, ty, deltaTX, deltaTY, x_t, y_t, slopeX, slopeY;
	int plate, cross_the_line,trid,nseg;
	tree->Branch("deltaX", &deltaX);
	tree->Branch("deltaY", &deltaY);
	tree->Branch("deltaTX", &deltaTX);
	tree->Branch("deltaTY", &deltaTY);
	tree->Branch("x", &x_t);
	tree->Branch("y", &y_t);
	tree->Branch("slopeX", &slopeX);
	tree->Branch("slopeY", &slopeY);
	tree->Branch("pl", &plate);
	tree->Branch("cross_the_line", &cross_the_line);
	tree->Branch("trid", &trid);
	tree->Branch("nseg", &nseg);
	double tx3;
	double stx;
	double ty3;
	double sty;
	// Loop over the tracks
	for (int itrk = 0; itrk < ntrk; itrk++)
	{

		stx = 0.;
		sty = 0.;
		EdbTrackP *t = pvr->GetTrack(itrk);
		if (abs(t->TX() + 0.01) >= 0.01 || abs(t->TY() - 0.004) >= 0.01 || t->N() < 5)
			continue;

		trid=t->ID();
		nseg = t->N();
		TGraph grX;
		TGraph grY;
		int pl;
		int plold = 10000;
		int consecutive_seg = 0;

		std::vector<double> vx(nseg);
		std::vector<double> vy(nseg);
		std::vector<double> vz(nseg);
		std::vector<double> vtx(nseg);
		std::vector<double> vty(nseg);

		double x3, y3, z3;
		for(int iPID=2;iPID<nPID-2;iPID++)
		{
			int count=0;
			double x[5];
			double y[5];
			double z[5];
			for (int iseg = 0; iseg < nseg; iseg++)
			{
				EdbSegP *s = t->GetSegment(iseg);
				for(int ipoint=0;ipoint<5;ipoint++)
				{
					if(s->PID()==iPID-2+ipoint)
					{
						x[ipoint]=s->X();
						y[ipoint]=s->Y();
						z[ipoint]=s->Z();
						count++;
					}
				}
				if(s->PID()==iPID)
				{
					tx3=s->TX();
					ty3=s->TY();
				}
			}
			if(count!=5) continue;
			// printf("%d\n",grX.GetN());
			int areaX[5];
			int areaY[5];
			for (int ipoint = 0; ipoint < 5; ipoint++)
			{
				areaX[ipoint] = (x[ipoint] - (Xcenter - XYrange)) / bin_width;
				areaY[ipoint] = (y[ipoint] - (Ycenter - XYrange)) / bin_width;
			}
			cross_the_line = 0;
			for (int ipoint = 0; ipoint < 5 - 1; ipoint++)
			{
				if (areaX[ipoint] != areaX[ipoint + 1] || areaY[ipoint] != areaY[ipoint + 1])
				{
					cross_the_line = 1;
					break;
				}
			}
			double x_updown[4],y_updown[4],z_updown[4];
			for(int i=0;i<2;i++)
			{
				x_updown[i] = x[i];
				y_updown[i] = y[i];
				z_updown[i] = z[i];
			}
			for(int i=2;i<4;i++)
			{
				x_updown[i] = x[i+1];
				y_updown[i] = y[i+1];
				z_updown[i] = z[i+1];
			}
			double a0;
			lsm(z_updown,x_updown,4,a0,slopeX);
			double x3fit = a0+slopeX*z[2];
			lsm(z_updown,y_updown,4,a0,slopeY);
			double y3fit = a0+slopeY*z[2];

			// Calculate delta X and delta Y.
			deltaX = x[2] - x3fit;
			deltaY = y[2] - y3fit;
			deltaTX = tx3 - slopeX;
			deltaTY = ty3 - slopeY;
			x_t = t->X();
			y_t = t->Y();
			plate = pvr->GetPattern(iPID)->Plate();
			tree->Fill();
		}
	}
	return 0;
}
int dedicated_align(TObjArray *tracks,double Xcenter, double Ycenter, double bin_width)
{
	gMinuit = TVirtualFitter::Fitter(0, 300);
	TNtupleD *sta = new TNtupleD("sta", "ShiftTAlign", "iX:iY:shiftX:shiftY:pid");
	int ntrk = tracks->GetEntriesFast();

	for (double iY = Ycenter - XYrange + bin_width / 2; iY <= Ycenter + XYrange; iY += bin_width) // Divide the area into 2*2 mm^2 areas
	{
		for (double iX = Xcenter - XYrange + bin_width / 2; iX <= Xcenter + XYrange; iX += bin_width)
		{
			TObjArray *tracks2 = new TObjArray;

			for (int itrk = 0; itrk < ntrk; itrk++)
			{
				EdbTrackP *t = (EdbTrackP *)tracks->At(itrk);
				// if (t->N()<10|| abs(t->TX() + 0.01) >= 0.01 || abs(t->TY()-0.004) >= 0.01||t->GetSegment(0)->PID()>=10)
				if (t->N() < 10 || abs(t->TX() + 0.01) >= 0.01 || abs(t->TY() - 0.004) >= 0.01)
					continue;
				// if (fabs(t->X() - iX) < bin_width / 2 && fabs(t->Y() - iY) < bin_width / 2)
				if (10 <= count_passed_seg(t, iX, iY, bin_width)) //  check if the track passes the area
					tracks2->Add(t);
			}
			// printf("iX = %.0f, iY = %.0f, ntrk = %d\n", iX, iY, tracks2->GetEntries());
			if (tracks2->GetEntries() == 0)
				continue;
			// calculate the alignment parameters several times.
			for (int j = 0; j < 1; j++)
			{
				calc_align_par(tracks2,iX,iY,bin_width,0,sta); //4th is fixflag
			}
			// Apply alignment parameter.
			for (int itrk = 0; itrk < ntrk; itrk++)
			{
				EdbTrackP *t = (EdbTrackP *)tracks->At(itrk);
				// if (fabs(t->X() - iX) < bin_width / 2 && fabs(t->Y() - iY) < bin_width / 2)
				apply_align(t, iX, iY, bin_width);
			}
			delete tracks2;
		}
	}
	// Ntuple for Shifts of TAlign
	// TFile fout1(Form("ShiftPar/sta_%.0fbinwidth_t1_NoFix_Robust%.1f_NoStraySeg_over10seg.root", bin_width, robustfactor), "recreate");
	// sta->Write();
	// fout1.Close();
	return 0;
}
int main(int argc, char *argv[])
{
	bool alignFlag = true;
	if(argc>5)
	{
		alignFlag = strcmp(argv[5], "--NoAlign") != 0;
	}else{
		printf("Usage: ./calc_dxy linked_tracks.root title Xcenter Ycenter bin_width robustfactor\n");
		printf("Usage: ./calc_dxy linked_tracks.root title Xcenter Ycenter --NoAlign\n");
		return 1;
	}
	if ( argc<7&&alignFlag)
	{
		printf("Usage: ./calc_dxy linked_tracks.root title Xcenter Ycenter bin_width robustfactor\n");
		printf("Usage: ./calc_dxy linked_tracks.root title Xcenter Ycenter --NoAlign\n");
		return 1;
	}

	TString filename_linked_tracks = argv[1];
	TString title = argv[2];
	double Xcenter, Ycenter, bin_width;
	bin_width=20000;
	sscanf(argv[3], "%lf", &Xcenter);
	sscanf(argv[4], "%lf", &Ycenter);
	TString Align_or_NoAlign = "";
	if(!alignFlag) Align_or_NoAlign = "_BeforeAlign";
	
	EdbDataProc *dproc = new EdbDataProc;
	EdbPVRec *pvr = new EdbPVRec;
	dproc->ReadTracksTree(*pvr, filename_linked_tracks, "1");

	nPID = pvr->Npatterns();
	int plMin = pvr->GetPattern(0)->Plate();
	int plMax = pvr->GetPattern(nPID-1)->Plate();

	TObjArray *tracks = pvr->GetTracks();
	int ntrk = tracks->GetEntriesFast();
	
	if (ntrk == 0)
	{
		printf("ntrk==0\n");
		return 0;
	}
	if (alignFlag) //if perform alignment
	{
		sscanf(argv[5], "%lf", &bin_width);
		sscanf(argv[6], "%f", &robustfactor);
		dedicated_align(tracks,Xcenter,Ycenter,bin_width);
	}

	TTree *tree = new TTree("tree", "deltaXY");
	TObjString *info = new TObjString(Form("plMin=%d, plMax=%d, Xcenter=%f, Ycenter=%f", plMin, plMax, Xcenter, Ycenter));
	tree->GetUserInfo()->Add(info);

	calc_dxy(pvr,tree,ntrk,Xcenter,Ycenter,bin_width);
	// Ntuple for deltaXY
	// TFile fout(Form("deltaXY%s_nt_reconnected_%.0fbinwidth_func2_t1.root", module,bin_width), "recreate");
	// TFile fout(Form("deltaXY/nt_aligntfd_NAlign_%.0f_%.0f.root",Xcenter,Ycenter), "recreate");
	// TFile fout(Form("deltaXY/nt_reconnected_%.0fbinwidth_t1_Fix_0to94_Robust%.1f_%.0f_%.0f.root",bin_width,robustfactor,Xcenter,Ycenter), "recreate");
	// TFile fout(Form("deltaXY/nt_aligntfd_reconnected_%.0fbinwidth_NAlign_%.0f_%.0f.root",bin_width,Xcenter,Ycenter), "recreate");
	// TFile fout(Form("deltaXY%s_nt_reconnected_%.0fbinwidth_func2_ttest.root", module,bin_width), "recreate");
	// TFile fout(Form("deltaXY/nt_%.0fbinwidth_t1_NoFix_Robust%.1f_NoStraySeg_over10seg.root",bin_width,robustfactor), "recreate");
	TFile fout("deltaXY/tree_" + title + Align_or_NoAlign + ".root", "recreate");
	tree->Write();
	fout.Close();
	// トラックアライメント後のlinked_tracksを作る
	//  TObjArray *tracks_t = new TObjArray;
	//  for(int itrk=0;itrk<ntrk;itrk++)
	//  {
	//  	EdbTrackP *t = (EdbTrackP *)tracks->At(itrk);
	//  	tracks_t->Add(t);
	//  }
	//  dproc->MakeTracksTree(*tracks_t,0,0,Form("/data/Users/kokui/FASERnu/F222/test/F222_zone3_vertex003_test2/reco43_095000_065000/v13/linked_tracks_AfterAlign_Robust%.1f_NoStraySeg_over10seg_StartAtFirst10.root",robustfactor));

	return 0;
	
}