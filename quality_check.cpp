#include <FnuQualityCheck.h>

#include <stdio.h>
#include <chrono>

#include <EdbDataSet.h>

int main(int argc, char *argv[])
{
	auto start = std::chrono::system_clock::now();
	if (argc < 5)
	{
		printf("Usage: ./quality_check linked_tracks.root title Xcenter Ycenter\n");
		return 1;
	}

	TString filename_linked_tracks = argv[1];
	TString title = argv[2];
	double Xcenter, Ycenter, bin_width;
	bin_width = 20000;
	sscanf(argv[3], "%lf", &Xcenter);
	sscanf(argv[4], "%lf", &Ycenter);
	// sscanf(argv[5], "%lf", &bin_width);

	EdbDataProc *dproc = new EdbDataProc;
	EdbPVRec *pvr = new EdbPVRec;

	dproc->ReadTracksTree(*pvr, filename_linked_tracks, "nseg>=4");
	// dproc->ReadTracksTree(*pvr, filename_linked_tracks, "Entry$%20000==0");
	// dproc->ReadTracksTree(*pvr, filename_linked_tracks, "Entry$<5000");

	TObjArray *tracks = pvr->GetTracks();
	int ntrk = tracks->GetEntriesFast();

	if (ntrk == 0)
	{
		printf("ntrk==0\n");
		return 0;
	}
	FnuQualityCheck qc(pvr, title);
	// qc.CalcDeltaXY(Xcenter, Ycenter, bin_width);
	// qc.CalcMeanDeltaXY(Xcenter, Ycenter, 500);
	// qc.PrintMeanDeltaXYArrowPlot("deltaXY_XYdis/arrow_deltaXY_" + title + ".pdf");
	// qc.WriteMeanDeltaXY("deltaXY_XYdis/mean_deltaXY_" + title + ".root"); // Writing should be done after all methods that use a related tree.
	// qc.FitDeltaXY();
	// TFile *fileDeltaXY = new TFile("/data/Users/kokui/FASERnu/F222/TFD_volumes/btfiltering/deltaXY/tree_zone3_rearranged_006_vert31025_pl013-127.root");
	// TTree *deltaXY = (TTree*)gDirectory->Get("deltaXY");
	// TFile *filetest = new TFile("test.root","recreate");
	// TTree *treeHistDetlaTXY = qc.MakeHistDeltaTXY(deltaXY);
	// TTree *angleResolutionPar = qc.FitHistDeltaTXY(treeHistDetlaTXY);
	// qc.MakeGraphHistAngleResolution(angleResolutionPar);
	// qc.PrintGraphHistAngleResolution("test_angle_resolution.pdf");
	// angleResolutionPar->Write();
	// qc.PrintDeltaTXYHist(treeHistDetlaTXY,"test_histDeltaTXY.pdf");
	// qc.MakePosResGraphHist();
	// TString outputDir = "/data/Users/kokui/FASERnu/F222/zone4/temp/TFD/vert32063_pl053_167_new/reco32_065000_050000/v15/";
	// qc.PrintPosResGraphHist(outputDir + "pos_res/sigma_par_" + title + ".pdf");
	// qc.WritePosResGraphHist(outputDir + "pos_res/graph_hist_" + title + ".root");
	// qc.PrintDeltaXYHist(outputDir + "pos_res/deltaxy_hist_" + title + ".pdf");
	// qc.WriteDeltaXY(outputDir + "deltaXY/tree_" + title + ".root");
	// qc.WritePosResPar(outputDir + "pos_res/sigma_par_" + title + ".root");
	// qc.CalcEfficiency();
	// qc.PrintEfficiency("efficiency_output/hist_efficiency_" + title + ".pdf");
	// qc.WriteEfficiency("efficiency_output/efficiency_" + title + ".root");
	// qc.WriteEfficiencyTree(Form("efficiency_output/effinfo_%s.root", title.Data()));
	// qc.MakePositionHist();
	// qc.PrintPositionHist("position_distribution_"+title+".pdf");
	// qc.WritePositionHist("position_distribution_"+title+".root");
	// qc.MakeAngleHist();
	// qc.PrintAngleHist("angle_distribution_"+title+".pdf");
	// qc.WriteAngleHist("angle_distribution_"+title+".root");
	qc.MakeNsegHist();
	// qc.PrintNsegHist("nseg_"+title+".pdf");
	// qc.WriteNsegHist("nseg_"+title+".root");
	qc.MakeNplHist();
	// qc.PrintNplHist("npl_" + title + ".pdf");
	// qc.WriteNplHist("npl_" + title + ".root");
	qc.MakeFirstLastPlateHist();
	// qc.PrintFirstLastPlateHist("first_last_plate_" + title + ".pdf");
	// qc.WriteFirstLastPlateHist("first_last_plate_" + title + ".root");

	// TTree *secondDifferenceTree = qc.CalcSecondDifference(32);
	
	// TFile fout(title+".root","recreate");
	// secondDifferenceTree->Write();
	// fout.Close();
	// TFile fin("test.root");
	// TTree *secondDifferenceTree = (TTree*)gDirectory->Get("secondDifferenceTree");
	// qc.MakeSecondDifferenceHist(secondDifferenceTree, 32);

	qc.Summarize(Xcenter, Ycenter, bin_width);

	auto end = std::chrono::system_clock::now();
	auto dur = end - start;
	auto secondsPassed = std::chrono::duration_cast<std::chrono::seconds>(dur).count();
	std::cout << secondsPassed << "seconds\n";
	return 0;
}
