#include "FnuDivideAlign.h"
#include <EdbDataSet.h>
int main(int argc, char *argv[])
{
	if (argc < 6)
	{
		printf("Usage: ./calc_dxy linked_tracks.root title reco binWidth robustFactor\n");
		return 1;
	}

	TString filename_linked_tracks = argv[1];
	TString title = argv[2];
	double binWidth;
	float robustFactor;
	int reco;
	sscanf(argv[3], "%d", &reco);
	double Xcenter = (reco - 1) % 9 * 15000 + 5000;
	double Ycenter = (reco - 1) / 9 * 15000 + 5000;
	sscanf(argv[4], "%lf", &binWidth);
	sscanf(argv[5], "%f", &robustFactor);

	EdbDataProc *dproc = new EdbDataProc;
	EdbPVRec *pvr = new EdbPVRec;
	dproc->ReadTracksTree(*pvr, filename_linked_tracks, "1");

	TObjArray *tracks = pvr->GetTracks();
	int nPatterns = pvr->Npatterns();
	int ntrk = tracks->GetEntriesFast();

	if (ntrk == 0)
	{
		printf("ntrk==0\n");
		return 0;
	}

	FnuDivideAlign align;
	align.SetRobustFactor(robustFactor);
	align.SetBinWidth(binWidth);
	align.Align(tracks, Xcenter, Ycenter, nPatterns);
	align.WriteAlignPar("align_output/alignPar_" + title + ".root");

	TObjArray *tracks_t = new TObjArray;
	for (int itrk = 0; itrk < ntrk; itrk++)
	{
		EdbTrackP *t = (EdbTrackP *)tracks->At(itrk);
		tracks_t->Add(t);
	}
	dproc->MakeTracksTree(*tracks_t, 0, 0, Form("/data/Users/kokui/FASERnu/F222/zone4/temp/TFD/vert32063_pl053_167_new/reco32_065000_050000/v15/linked_tracks_after_align_binWidth%.0f_robustFactor%.1f.root", binWidth, robustFactor));

	return 0;
}