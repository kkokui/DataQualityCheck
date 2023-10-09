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
	int ntrk = tracks->GetEntriesFast();
	
	if (ntrk == 0)
	{
		printf("ntrk==0\n");
		return 0;
	}
	
	FnuDivideAlign align;
	align.SetRobustFactor(robustfactor);
	align.SetBinWidth(bin_width);
	align.dedicated_align(pvr, Xcenter, Ycenter);
	align.WriteAlignPar("alignPar_"+title + ".root");
    return 0;
}