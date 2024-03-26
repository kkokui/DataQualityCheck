#include "FnuQualityCheck.h"
#include "FnuResolution.h"
#include "FnuEfficiency.h"
#include "FnuPositionDistribution.h"
#include "FnuAngleDistribution.h"

#include <stdio.h>
#include <numeric>

#include <EdbDataSet.h>
#include <TGraph.h>
#include <TH2.h>
#include <TObjArray.h>
#include <EdbPattern.h>
#include <TMath.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TGaxis.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TEventList.h>
#include <TGraphAsymmErrors.h>
#include <TPaletteAxis.h>
#include <TPaveStats.h>
#include <TArrow.h>
#include <TLatex.h>

FnuQualityCheck::FnuQualityCheck(EdbPVRec *pvr, TString title)
	: pvr(pvr),
	  title(title),
	  nPID(pvr->Npatterns()),
	  plMin(pvr->GetPattern(0)->Plate()),
	  plMax(pvr->GetPattern(nPID - 1)->Plate()),
	  ntrk(pvr->Ntracks()),
	  deltaTXV(), deltaXV(), deltaYV(), deltaTYV(), xV(), yV(), slopeXV(), slopeYV(),
	  crossTheLineV(), tridV(), nsegV(),
	  angleCut(0.01),
	  XYrange(8500)
{
}

FnuQualityCheck::~FnuQualityCheck()
{
}

TH1I *FnuQualityCheck::MakePHHist(int nsegMin)
{
	TH1I *histPH = new TH1I(Form("histPHNsegMin%d",nsegMin), Form("PH (nseg>=%d)", nsegMin), 20, 13 - 0.5, 32 + 0.5);
	for (int itrk = 0; itrk < ntrk; itrk++)
	{
		EdbTrackP *t = pvr->GetTrack(itrk);
		int nseg = t->N();
		if (nseg < nsegMin)
			continue;
		for (int iseg = 0; iseg < nseg; iseg++)
		{
			int PH = t->GetSegment(iseg)->W();
			histPH->Fill(PH);
		}
	}
	return histPH;
}
TH1D *FnuQualityCheck::MakePHMeanHist(int nsegMin)
{
	TH1D *histPHMean = new TH1D("histPHMean", Form("PH mean (nseg>=%d)", nsegMin), 20, 13 - 0.5, 32 + 0.5);
	for (int itrk = 0; itrk < ntrk; itrk++)
	{
		EdbTrackP *t = pvr->GetTrack(itrk);
		int nseg = t->N();
		if (nseg < nsegMin)
			continue;
		double PHMean = 0;
		for (int iseg = 0; iseg < nseg; iseg++)
		{
			PHMean += t->GetSegment(iseg)->W();
		}
		PHMean /= nseg;
		histPHMean->Fill(PHMean);
	}
	return histPHMean;
}

void FnuQualityCheck::MakeNsegHist()
{
	nsegHist = new TH1I("nsegHist", "nseg (" + title + ");nseg;N tracks", nPID, 0.5, nPID + 0.5);
	for (int itrk = 0; itrk < ntrk; itrk++)
	{
		int nseg = pvr->GetTrack(itrk)->N();
		nsegHist->Fill(nseg);
	}
}
void FnuQualityCheck::PrintNsegHist(TString filename)
{
	TCanvas ctemp;
	gPad->SetLogy();
	nsegHist->Draw();
	ctemp.Print(filename);
}
void FnuQualityCheck::WriteNsegHist(TString filename)
{
	TFile fout(filename, "recreate");
	nsegHist->Write();
	fout.Close();
}
void FnuQualityCheck::MakeNplHist()
{
	nplHist = new TH1I("nplHist", "npl (" + title + ");npl;N tracks", plMax - plMin + 1, 0.5, plMax - plMin + 1.5);
	for (int itrk = 0; itrk < ntrk; itrk++)
	{
		int npl = pvr->GetTrack(itrk)->Npl();
		nplHist->Fill(npl);
	}
}
void FnuQualityCheck::PrintNplHist(TString filename)
{
	TCanvas ctemp;
	gPad->SetLogy();
	nplHist->Draw();
	ctemp.Print(filename);
}
void FnuQualityCheck::WriteNplHist(TString filename)
{
	TFile fout(filename, "recreate");
	nplHist->Write();
	fout.Close();
}
void FnuQualityCheck::MakeFirstLastPlateHist()
{
	firstPlateHist = new TH1I("firstPlateHist", "first plate (" + title + ");plate;N tracks", plMax - plMin + 1, plMin - 0.5, plMax + 0.5);
	lastPlateHist = new TH1I("lastPlateHist", "last plate (" + title + ");plate;N tracks", plMax - plMin + 1, plMin - 0.5, plMax + 0.5);
	for (int itrk = 0; itrk < ntrk; itrk++)
	{
		EdbTrackP *t = pvr->GetTrack(itrk);
		firstPlateHist->Fill(t->GetSegmentFirst()->Plate());
		lastPlateHist->Fill(t->GetSegmentLast()->Plate());
	}
}
void FnuQualityCheck::PrintFirstLastPlateHist(TString filename)
{
	TString originalHistTitle = firstPlateHist->GetTitle();
	firstPlateHist->SetTitle("start and end plate");
	gStyle->SetPadTopMargin(0.13);
	TCanvas c;
	gStyle->SetOptStat("e");
	gStyle->SetStatY(0.87);
	firstPlateHist->Draw();
	lastPlateHist->SetLineColor(kRed);
	lastPlateHist->Draw("same");
	TLegend *legFirstLast = new TLegend(0.2, 0.87, 0.9, 0.93);
	legFirstLast->AddEntry(firstPlateHist, "first plate", "l");
	legFirstLast->AddEntry(lastPlateHist, "last plate", "l");
	legFirstLast->SetNColumns(2);
	legFirstLast->Draw();
	c.Print(filename);
	firstPlateHist->SetTitle(originalHistTitle);
	gROOT->SetStyle("Modern");
	gStyle->SetOptStat();
	gStyle->SetStatY(0.935);
	firstPlateHist->UseCurrentStyle();
	lastPlateHist->UseCurrentStyle();
}
void FnuQualityCheck::WriteFirstLastPlateHist(TString filename)
{
	TFile fout(filename, "recreate");
	firstPlateHist->Write();
	lastPlateHist->Write();
	fout.Close();
}

TTree *FnuQualityCheck::CalcSecondDifference(int cellLength)
{
	// Calculate second difference for checking long range position displacement.

	TTree *secondDifferenceTree = new TTree("secondDifferenceTree", "Tree of 2nd difference");
	secondDifferenceX = new double[ntrk];
	secondDifferenceY = new double[ntrk];
	slopeX = new double[ntrk];
	slopeY = new double[ntrk];
	secondDifferenceTree->Branch("plate", &plate);
	secondDifferenceTree->Branch("entriesParPlate", &entriesParPlate);
	secondDifferenceTree->Branch("secondDiffX", secondDifferenceX, "secondDiffX[entriesParPlate]/D");
	secondDifferenceTree->Branch("secondDiffY", secondDifferenceY, "secondDiffY[entriesParPlate]/D");
	secondDifferenceTree->Branch("slopeX", slopeX, "slopeX[entriesParPlate]/D");
	secondDifferenceTree->Branch("slopeY", slopeY, "slopeY[entriesParPlate]/D");
	for (int iPID = 0; iPID < nPID; iPID++)
	{
		entriesParPlate = 0;
		EdbPattern *pattern = pvr->GetPattern(iPID);
		if (pattern == NULL)
			continue;
		int plate0 = pattern->Plate();
		int plate1 = plate0 + cellLength;
		int plate2 = plate1 + cellLength;
		// loop over the tracks
		for (int itrk = 0; itrk < ntrk; itrk++)
		{
			EdbTrackP *t = pvr->GetTrack(itrk);
			const int npl = t->Npl();
			if (npl <= cellLength * 2)
			{
				continue;
			}
			// printf("this tracks has enough npl\n");
			int nseg = t->N();
			int count = 0;
			int iseg0, iseg1, iseg2;
			// loop over the segments
			for (int iseg = 0; iseg < nseg; iseg++)
			{
				EdbSegP *s = t->GetSegment(iseg);
				int plateOfTheSegment = s->Plate();
				if (plateOfTheSegment > plate2)
				{
					break;
				}
				if (plateOfTheSegment == plate0)
				{
					count++;
					iseg0 = iseg;
				}
				if (plateOfTheSegment == plate1)
				{
					count++;
					iseg1 = iseg;
				}
				if (plateOfTheSegment == plate2)
				{
					count++;
					iseg2 = iseg;
				}
			}
			if (3 != count)
			{
				continue;
			}
			EdbSegP *s0 = t->GetSegment(iseg0);
			EdbSegP *s1 = t->GetSegment(iseg1);
			EdbSegP *s2 = t->GetSegment(iseg2);

			double X0 = s0->X();
			double Y0 = s0->Y();
			double Z0 = s0->Z();

			double X1 = s1->X();
			double Y1 = s1->Y();
			double Z1 = s1->Z();

			double X2 = s2->X();
			double Y2 = s2->Y();
			double Z2 = s2->Z();

			double X2prediction = X1 + (X1 - X0) / (Z1 - Z0) * (Z2 - Z1); // prediction by extrapolating with z correction
			double Y2prediction = Y1 + (Y1 - Y0) / (Z1 - Z0) * (Z2 - Z1); // prediction by extrapolating with z correction
			secondDifferenceX[entriesParPlate] = X2 - X2prediction;
			secondDifferenceY[entriesParPlate] = Y2 - Y2prediction;
			slopeX[entriesParPlate] = (X1 - X0) / (Z1 - Z0);
			slopeY[entriesParPlate] = (Y1 - Y0) / (Z1 - Z0);
			entriesParPlate++;
		}
		plate = plate2;
		secondDifferenceTree->Fill();
	}
	return secondDifferenceTree;
}
void FnuQualityCheck::MakeSecondDifferenceHist(TTree *secondDifferenceTree, int cellLength)
{
	grErrorsSecondDifferenceMeanX = new TGraphErrors();
	grErrorsSecondDifferenceMeanY = new TGraphErrors();
	grErrorsSecondDifferenceMeanX->SetTitle(Form("mean of second difference in x (cell length %d);plate;mean of second difference (#mum)", cellLength));
	grErrorsSecondDifferenceMeanY->SetTitle(Form("mean of second difference in y (cell length %d);plate;mean of second difference (#mum)", cellLength));
	secondDifferenceHist = new TH1D("secondDifferenceHist", Form("second difference (cell length %d);second differecne (#mum);entries", cellLength), 100, -30, 30);
	int ipoint = 0;
	for (int ipl = plMin + cellLength * 2; ipl <= plMax; ipl++)
	{
		secondDifferenceTree->Draw("secondDiffX>>secondDifferenceHist", Form("plate==%d", ipl));
		grErrorsSecondDifferenceMeanX->SetPoint(ipoint, ipl, secondDifferenceHist->GetMean());
		grErrorsSecondDifferenceMeanX->SetPointError(ipoint, 0, secondDifferenceHist->GetMeanError());
		secondDifferenceTree->Draw("secondDiffY>>secondDifferenceHist", Form("plate==%d", ipl));
		grErrorsSecondDifferenceMeanY->SetPoint(ipoint, ipl, secondDifferenceHist->GetMean());
		grErrorsSecondDifferenceMeanY->SetPointError(ipoint, 0, secondDifferenceHist->GetMeanError());
		ipoint++;
		// secondDifferenceHist->Reset();
	}
}
void FnuQualityCheck::Summarize(double Xcenter, double Ycenter, double bin_width)
{
	gStyle->SetPadLeftMargin(0.14);
	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadTopMargin(0.07);
	gStyle->SetTitleOffset(1.1, "Y");
	gStyle->SetTitleOffset(1.0, "Z");
	gStyle->SetTitleSize(0.06, "XYZ");
	gStyle->SetLabelSize(0.05, "XYZ");
	gStyle->SetNdivisions(505, "xyz");
	TGaxis::SetMaxDigits(4);
	TCanvas c("c", "summary plot", 1600, 900);
	c.Divide(4, 3);

	// position distribution
	c.cd(1);
	// gPad->SetRightMargin(0.4);
	gStyle->SetPadRightMargin(0.25);
	gStyle->SetOptStat("emr");
	gStyle->SetStatH(0.25);
	gStyle->SetStatW(0.25);
	gStyle->SetStatX(1);
	gStyle->SetStatY(1);
	FnuPositionDistribution pos(pvr, title);
	positionHist = pos.MakePositionHist();
	TString positionHistOriginalTitle = positionHist->GetTitle();
	positionHist->SetTitle("position distribution");
	positionHist->Draw("colz");
	gPad->UseCurrentStyle();
	gPad->Update();
	// gStyle->SetPadLeftMargin(0.15);
	// gStyle->SetTitleOffset(1.5, "Y");
	// ((TGaxis *)positionHist->GetXaxis())->SetMaxDigits(2);
	// ((TGaxis *)positionHist->GetYaxis())->SetMaxDigits(2);
	TPaletteAxis *palette = (TPaletteAxis *)positionHist->GetListOfFunctions()->FindObject("palette");
	palette->SetX1NDC(0.845);
	palette->SetX2NDC(0.88);
	palette->SetY2NDC(0.68);
	// gPad->RedrawAxis();
	gPad->UseCurrentStyle();

	// angle distribution
	FnuAngleDistribution angleDistribution(pvr, title);
	// angleDistribution.CalcAngle();
	// angle distribution wide
	c.cd(2);
	gStyle->SetPadRightMargin(0.25);
	angleHistWide = angleDistribution.MakeAngleHistWide();
	TString angleHistWideOriginalTitle = angleHistWide->GetTitle();
	angleHistWide->SetTitle("angle distribution wide");
	angleHistWide->Draw("colz");
	gPad->UseCurrentStyle();
	gPad->Update();
	palette = (TPaletteAxis *)angleHistWide->GetListOfFunctions()->FindObject("palette");
	palette->SetX1NDC(0.845);
	palette->SetX2NDC(0.88);
	palette->SetY2NDC(0.68);
	gPad->UseCurrentStyle();
	gPad->SetLogz();

	// angle distribution narrow
	c.cd(3);
	gStyle->SetPadRightMargin(0.25);
	angleHistNarrow = angleDistribution.MakeAngleHistNarrow();
	TString angleHistNarrowOriginalTitle = angleHistNarrow->GetTitle();
	angleHistNarrow->SetTitle("angle distribution narrow");
	angleHistNarrow->Draw("colz");
	gPad->UseCurrentStyle();
	gPad->Update();
	palette = (TPaletteAxis *)angleHistNarrow->GetListOfFunctions()->FindObject("palette");
	palette->SetX1NDC(0.845);
	palette->SetX2NDC(0.88);
	palette->SetY2NDC(0.68);
	gPad->UseCurrentStyle();
	gPad->SetLogz();

	// pulse height
	c.cd(4);
	gStyle->SetPadRightMargin(0.0);
	gStyle->SetPadTopMargin(0.13);
	gStyle->SetOptStat("e");
	gStyle->SetStatH(0.27);
	gStyle->SetStatW(0.3);
	gStyle->SetStatX(1);
	gStyle->SetStatY(1);
	TH1I *PHHistNsegMin2 = MakePHHist(2);
	TH1I *PHHistNsegMin4 = MakePHHist(4);
	TH1D *PHHistMean = MakePHMeanHist(4);
	auto PHHistMaxY = PHHistNsegMin2->GetMaximum();
	auto PHHistminX = PHHistNsegMin2->GetBinCenter(PHHistNsegMin2->FindFirstBinAbove());
	auto PHHistmaxX = PHHistNsegMin2->GetBinCenter(PHHistNsegMin2->FindLastBinAbove());
	gPad->DrawFrame(PHHistminX - 0.5, 0, PHHistmaxX + 0.5, PHHistMaxY * 1.1, "PH;PH;N segments or N tracks");
	
	PHHistNsegMin2->Draw("sames");
	PHHistNsegMin4->Draw("same");
	PHHistMean->Draw("same");
	gPad->UseCurrentStyle();
	gPad->Update();
	PHHistNsegMin2->SetLineColor(kBlack);
	PHHistMean->SetFillColor(kYellow);
	TLegend *legPH = new TLegend(0.3, 0.87, 1.0, 0.93);
	legPH->AddEntry(PHHistNsegMin2, "nseg#geq2", "l");
	legPH->AddEntry(PHHistNsegMin4, "nseg#geq4", "l");
	legPH->AddEntry(PHHistMean, "mean, nseg#geq4", "lf");
	legPH->SetNColumns(3);
	legPH->Draw();

	// efficiency
	FnuEfficiency eff(pvr, title);
	eff.CalcEfficiency();
	// efficiency for each angle
	c.cd(5);
	TEfficiency *efficiencyForEachAngle = eff.GetEfficiencyForEachAngle();
	TString efficiencyForEachAngleOriginalTitle = efficiencyForEachAngle->GetTitle();
	efficiencyForEachAngle->SetTitle("efficiency for each angle");
	gPad->SetRightMargin(0.0);
	efficiencyForEachAngle->Draw();

	// efficiency for each plate
	c.cd(6);
	TEfficiency *efficiencyForEachPlate = eff.GetEfficiencyForEachPlate();
	TString efficiencyForEachPlateOriginalTitle = efficiencyForEachPlate->GetTitle();
	efficiencyForEachPlate->SetTitle("efficiency for each plate");
	gPad->SetRightMargin(0.0);
	efficiencyForEachPlate->Draw();

	// position resolution in x and y for each plate
	c.cd(7);
	c.SetGridx(1);
	gPad->SetRightMargin(0.0);
	gPad->SetTopMargin(0.13);
	TMultiGraph *mg2 = new TMultiGraph("mg2", "position resolution for each plate;plate;position resolution (#mum)");
	FnuResolution res(pvr, title);
	deltaXY = res.CalcDeltaXY(Xcenter, Ycenter, bin_width);
	// htree = res.MakeHistDeltaXY(deltaXY);
	posResPar = res.FitHistDeltaXY(deltaXY);
	// res.WriteHistDeltaXYWithFit("hist_deltaXY_with_fit.root");
	res.MakeGraphHistPositionResolution(posResPar);
	sigmaXGraph = res.GetSigmaXGraph();
	sigmaYGraph = res.GetSigmaYGraph();
	sigmaXGraph->SetMarkerSize(0.5);
	sigmaYGraph->SetMarkerSize(0.5);
	mg2->Add(sigmaXGraph);
	mg2->Add(sigmaYGraph);
	mg2->Draw("ap");
	TLegend *leg2 = new TLegend(0.3, 0.87, 1.0, 0.93);
	leg2->AddEntry(sigmaXGraph, "x direction", "p");
	leg2->AddEntry(sigmaYGraph, "y direction", "p");
	leg2->SetNColumns(2);
	leg2->Draw();

	// angular resolution in x and y for each plate
	c.cd(8);
	c.SetGridx(1);
	gPad->SetRightMargin(0.0);
	gPad->SetTopMargin(0.13);
	TMultiGraph *mgAngleResolution = new TMultiGraph("mgAngleResolution", "angular resolution for each plate;plate;angular resolution (rad)");

	angResPar = res.FitHistDeltaTXY(deltaXY);
	// res.WriteHistDeltaTXYWithFit("hist_deltaTXY_with_fit.root");
	res.MakeGraphHistAngleResolution(angResPar);
	sigmaTXGraph = res.GetSigmaTXGraph();
	sigmaTYGraph = res.GetSigmaTYGraph();
	sigmaTXGraph->SetMarkerSize(0.5);
	sigmaTYGraph->SetMarkerSize(0.5);
	mgAngleResolution->Add(sigmaTXGraph);
	mgAngleResolution->Add(sigmaTYGraph);
	mgAngleResolution->Draw("ap");
	TLegend *legAngleResolution = new TLegend(0.3, 0.87, 1.0, 0.93);
	legAngleResolution->AddEntry(sigmaTXGraph, "x direction", "p");
	legAngleResolution->AddEntry(sigmaTYGraph, "y direction", "p");
	legAngleResolution->SetNColumns(2);
	legAngleResolution->Draw();

	// npl and nseg
	c.cd(9);
	gStyle->SetPadRightMargin(0.0);
	gStyle->SetPadTopMargin(0.13);
	gStyle->SetOptStat("e");
	gStyle->SetStatH(0.09);
	gStyle->SetStatW(0.3);
	gStyle->SetStatX(1);
	gStyle->SetStatY(1);
	auto nplHistMaxY = nplHist->GetMaximum();
	auto nsegHistMaxY = nsegHist->GetMaximum();
	auto maxY = std::max(nplHistMaxY, nsegHistMaxY);
	auto minX = nplHist->GetBinCenter(nplHist->FindFirstBinAbove());
	auto maxX = nplHist->GetBinCenter(nplHist->FindLastBinAbove());
	gPad->DrawFrame(minX - (maxX - minX) * 0.05, 0.7, maxX + (maxX - minX) * 0.05, maxY * 2, "nseg and npl;nseg or npl;N tracks");
	nplHist->Draw("sames");
	nsegHist->Draw("same");
	gPad->UseCurrentStyle();
	gPad->Update();
	nsegHist->SetLineColor(kRed);
	TLegend *legNsegNpl = new TLegend(0.3, 0.87, 1.0, 0.93);
	legNsegNpl->AddEntry(nsegHist, "nseg", "l");
	legNsegNpl->AddEntry(nplHist, "npl", "l");
	legNsegNpl->SetNColumns(2);
	legNsegNpl->Draw();
	gPad->SetLogy();

	// first and last plate
	c.cd(10);
	gPad->SetTopMargin(0.13);
	gPad->SetRightMargin(0);
	// gStyle->SetOptStat("e");
	// gStyle->SetStatY(0.87);
	auto firstPlateHistMaxY = firstPlateHist->GetMaximum();
	auto lastPlateHistMaxY = lastPlateHist->GetMaximum();
	maxY = std::max(firstPlateHistMaxY, lastPlateHistMaxY);
	gPad->DrawFrame(plMin - (plMax - plMin) * 0.05, 0.7, plMax + (plMax - plMin) * 0.05, maxY * 2, "start and end plate;plate;N tracks");
	firstPlateHist->Draw("sames");
	lastPlateHist->Draw("same");
	TLegend *legFirstLast = new TLegend(0.3, 0.87, 1.0, 0.93);
	legFirstLast->AddEntry(firstPlateHist, "first plate", "l");
	legFirstLast->AddEntry(lastPlateHist, "last plate", "l");
	legFirstLast->SetNColumns(2);
	legFirstLast->Draw();
	gPad->UseCurrentStyle();
	lastPlateHist->SetLineColor(kRed);
	gPad->SetLogy();

	// second difference
	c.cd(11);
	int cellLengthMaxAvailable = (nPID - 1) / 2;
	if (cellLengthMaxAvailable >= 32)
		cellLengthMaxAvailable = 32;
	if (cellLengthMaxAvailable >= 16 && cellLengthMaxAvailable < 32)
		cellLengthMaxAvailable = 16;
	if (cellLengthMaxAvailable >= 8 && cellLengthMaxAvailable < 16)
		cellLengthMaxAvailable = 8;
	if (cellLengthMaxAvailable >= 4 && cellLengthMaxAvailable < 8)
		cellLengthMaxAvailable = 4;
	if (cellLengthMaxAvailable >= 2 && cellLengthMaxAvailable < 4)
		cellLengthMaxAvailable = 2;
	if (cellLengthMaxAvailable >= 1 && cellLengthMaxAvailable < 2)
		cellLengthMaxAvailable = 1;
	TTree *secondDifferenceTree = CalcSecondDifference(cellLengthMaxAvailable);
	MakeSecondDifferenceHist(secondDifferenceTree, cellLengthMaxAvailable);
	gPad->SetTopMargin(0.13);
	gPad->SetRightMargin(0);
	// gStyle->SetOptStat("e");
	// gStyle->SetStatY(0.87);
	TMultiGraph *mgSecondDifference = new TMultiGraph("mgSecondDifference", Form("mean of second difference (cell length %d);plate;mean of second difference (#mum)", cellLengthMaxAvailable));
	grErrorsSecondDifferenceMeanX->SetLineColor(kRed);
	grErrorsSecondDifferenceMeanY->SetLineColor(kBlue);
	grErrorsSecondDifferenceMeanX->SetMarkerColor(kRed);
	grErrorsSecondDifferenceMeanY->SetMarkerColor(kBlue);
	grErrorsSecondDifferenceMeanX->SetMarkerStyle(20);
	grErrorsSecondDifferenceMeanY->SetMarkerStyle(20);
	grErrorsSecondDifferenceMeanX->SetMarkerSize(0.5);
	grErrorsSecondDifferenceMeanY->SetMarkerSize(0.5);
	mgSecondDifference->Add(grErrorsSecondDifferenceMeanX);
	mgSecondDifference->Add(grErrorsSecondDifferenceMeanY);
	mgSecondDifference->Draw("ap");
	TLegend *legSecondDifference = new TLegend(0.3, 0.87, 1.0, 0.93);
	legSecondDifference->AddEntry(grErrorsSecondDifferenceMeanX, "x direction", "p");
	legSecondDifference->AddEntry(grErrorsSecondDifferenceMeanY, "y direction", "p");
	legSecondDifference->SetNColumns(2);
	legSecondDifference->Draw();
	gPad->UseCurrentStyle();

	c.cd(12);
	TText tx;
	tx.DrawTextNDC(0.1, 0.9, title);

	c.Print("summary_plot_" + title + ".pdf");

	// set original title
	positionHist->SetTitle(positionHistOriginalTitle);
	angleHistNarrow->SetTitle(angleHistNarrowOriginalTitle);
	efficiencyForEachAngle->SetTitle(efficiencyForEachAngleOriginalTitle);
	efficiencyForEachPlate->SetTitle(efficiencyForEachPlateOriginalTitle);

	// set default style
	gROOT->SetStyle("Modern");
	gStyle->SetOptStat();
	gStyle->SetStatY(0.935);
	positionHist->UseCurrentStyle();
	angleHistNarrow->UseCurrentStyle();
	efficiencyForEachPlate->UseCurrentStyle();
	efficiencyForEachAngle->UseCurrentStyle();
	sigmaXGraph->UseCurrentStyle();
	sigmaYGraph->UseCurrentStyle();
	nsegHist->UseCurrentStyle();
	nplHist->UseCurrentStyle();
	firstPlateHist->UseCurrentStyle();
	lastPlateHist->UseCurrentStyle();

	// res.PrintHistDeltaXY("hist_deltaXY_" + title + ".pdf");
	// TFile foutTreeHist("tree_hist_"+title+".root","recreate");
	// htree->Write();
}