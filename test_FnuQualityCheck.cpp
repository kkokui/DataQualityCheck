#include "FnuQualityCheck.h"

#include <stdio.h>

#include <EdbDataSet.h>


int main(int argc, char *argv[])
{
	if ( argc<5)
	{
		printf("Usage: ./calc_dxy linked_tracks.root title Xcenter Ycenter\n");
		return 1;
	}

	TString filename_linked_tracks = argv[1];
	TString title = argv[2];
	double Xcenter, Ycenter, bin_width;
	bin_width=20000;
	sscanf(argv[3], "%lf", &Xcenter);
	sscanf(argv[4], "%lf", &Ycenter);
	
	EdbDataProc *dproc = new EdbDataProc;
	EdbPVRec *pvr = new EdbPVRec;
	
	// dproc->ReadTracksTree(*pvr, filename_linked_tracks, "nseg>=5");
	dproc->ReadTracksTree(*pvr, filename_linked_tracks, "Entry$<5000");

	TObjArray *tracks = pvr->GetTracks();
	int ntrk = tracks->GetEntriesFast();
	
	if (ntrk == 0)
	{
		printf("ntrk==0\n");
		return 0;
	}
    FnuQualityCheck qc(title);
	qc.CalcDeltaXY(pvr,ntrk,Xcenter,Ycenter,bin_width);
    qc.FitDeltaXY();
    qc.WritePosResPar();
    qc.PlotPosRes();
	
	return 0;
	
}