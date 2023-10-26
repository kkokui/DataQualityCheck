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
	
	dproc->ReadTracksTree(*pvr, filename_linked_tracks, "nseg>=4");
	// dproc->ReadTracksTree(*pvr, filename_linked_tracks, "Entry$<5000");

	TObjArray *tracks = pvr->GetTracks();
	int ntrk = tracks->GetEntriesFast();
	
	if (ntrk == 0)
	{
		printf("ntrk==0\n");
		return 0;
	}
    FnuQualityCheck qc(pvr, title);
	// qc.CalcDeltaXY(Xcenter,Ycenter,bin_width);
	// qc.FitDeltaXY();
	// qc.PlotPosRes("pos_res/sigmaPar_" + title + ".pdf");
	// qc.PrintDeltaXYHist("pos_res/deltaxy_" + title + ".pdf");
	// qc.WritePosResPar("pos_res/sigmaPar_" + title + ".root");
	// qc.WriteDeltaXY("deltaXY/tree_" + title + ".root");
	qc.CalcEfficiency();

	return 0;
	
}