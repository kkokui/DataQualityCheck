#include "FnuDivideAlign.h"

int main(int argc, char *argv[])
{
	if(argc<7)
	{
		printf("Usage: ./calc_dxy linked_tracks.root title Xcenter Ycenter bin_width robustfactor\n");
		return 1;
	}

	TString filename_linked_tracks = argv[1];
	TString title = argv[2];
	double Xcenter, Ycenter, bin_width;
    float robustfactor;
	sscanf(argv[3], "%lf", &Xcenter);
	sscanf(argv[4], "%lf", &Ycenter);
	sscanf(argv[5], "%lf", &bin_width);
	sscanf(argv[6], "%f", &robustfactor);

	EdbDataProc *dproc = new EdbDataProc;
	EdbPVRec *pvr = new EdbPVRec;
	dproc->ReadTracksTree(*pvr, filename_linked_tracks, "1");

	TObjArray *tracks = pvr->GetTracks();
	double test_x1 = pvr->GetTrack(10)->GetSegment(1)->X();
	int ntrk = tracks->GetEntriesFast();
	
	if (ntrk == 0)
	{
		printf("ntrk==0\n");
		return 0;
	}
	
	FnuDivideAlign da;
	da.SetRobustFactor(robustfactor);
	da.dedicated_align(pvr, Xcenter, Ycenter);
	double test_x2 = pvr->GetTrack(10)->GetSegment(1)->X();
	printf("%f %f\n",test_x1,test_x2);
	da.WriteAlignPar("alignPar_test.root");
    return 0;
}