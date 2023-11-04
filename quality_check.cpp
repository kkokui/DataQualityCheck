#include <FnuQualityCheck.h>

#include <stdio.h>

#include <EdbDataSet.h>

int main(int argc, char *argv[])
{
	if (argc < 6)
	{
		printf("Usage: ./test_FnuQualityCheck linked_tracks.root title Xcenter Ycenter binWidth\n");
		return 1;
	}

	TString filename_linked_tracks = argv[1];
	TString title = argv[2];
	double Xcenter, Ycenter, bin_width;
	bin_width = 20000;
	sscanf(argv[3], "%lf", &Xcenter);
	sscanf(argv[4], "%lf", &Ycenter);
	sscanf(argv[5], "%lf", &bin_width);

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
	FnuQualityCheck qc(pvr, title);
	qc.CalcDeltaXY(Xcenter, Ycenter, bin_width);
	qc.FitDeltaXY();
	qc.MakePosResGraphHist();
	// TString outputDir = "/data/Users/kokui/FASERnu/F222/zone4/temp/TFD/vert32063_pl053_167_new/reco32_065000_050000/v15/";
	// qc.PrintPosResGraphHist(outputDir + "pos_res/sigma_par_" + title + ".pdf");
	// qc.WritePosResGraphHist(outputDir + "pos_res/graph_hist_" + title + ".root");
	// qc.PrintDeltaXYHist(outputDir + "pos_res/deltaxy_hist" + title + ".pdf");
	// qc.WritePosResPar(outputDir + "pos_res/sigma_par_" + title + ".root");
	// qc.WriteDeltaXY(outputDir + "deltaXY/tree_" + title + ".root");
	qc.CalcEfficiency();
	// qc.PlotEfficiency("efficiency_output/hist_efficiency_" + title + ".pdf");
	// qc.WriteEfficiency("efficiency_output/efficiency_" + title + ".root");
	// qc.WriteEfficiencyTree(Form("efficiency_output/effinfo_%s.root", title.Data()));
	qc.MakePositionHist();
	// qc.PrintPositionHist("position_distribution_"+title+".pdf");
	// qc.WritePositionHist("position_distribution_"+title+".root");
	qc.MakeAngleHist();
	// qc.PrintAngleHist("angle_distribution_"+title+".pdf");
	// qc.WriteAngleHist("angle_distribution_"+title+".root");
	qc.PrintSummaryPlot();
	return 0;
}