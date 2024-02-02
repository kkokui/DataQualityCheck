#include "FnuQualityCheck.h"

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
	double bins_arr_angle[] = {0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5};
	SetBinsAngle(26, bins_arr_angle);
	double bins_arr_TXTY[] = {-0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.19, -0.18, -0.17, -0.16, -0.15, -0.14, -0.13, -0.12, -0.11, -0.10, -0.09, -0.08, -0.07, -0.06, -0.05, -0.04, -0.03, -0.02, -0.01, 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5};
	SetBinsTXTY(52, bins_arr_TXTY);
}

void FnuQualityCheck::SetBinsAngle(int nbins, double bins[])
{
	bins_vec_angle.assign(&bins[0], &bins[nbins + 1]);
}

void FnuQualityCheck::SetBinsTXTY(int nbins, double bins[])
{
	bins_vec_TXTY.assign(&bins[0], &bins[nbins + 1]);
}

FnuQualityCheck::~FnuQualityCheck()
{
}

void FnuQualityCheck::CalcDeltaXY(double Xcenter, double Ycenter, double bin_width)
{
	deltaXY = new TTree("deltaXY", "deltaXY");
	deltaXY->Branch("deltaX", &deltaXV);
	deltaXY->Branch("deltaY", &deltaYV);
	deltaXY->Branch("deltaTX", &deltaTXV);
	deltaXY->Branch("deltaTY", &deltaTYV);
	deltaXY->Branch("x", &xV);
	deltaXY->Branch("y", &yV);
	deltaXY->Branch("slopeX", &slopeXV);
	deltaXY->Branch("slopeY", &slopeYV);
	deltaXY->Branch("cross_the_line", &crossTheLineV);
	deltaXY->Branch("trid", &tridV);
	deltaXY->Branch("nseg", &nsegV);
	deltaXY->Branch("plate", &plate);
	// Loop over the tracks
	// for (int iPID = 2; iPID < nPID - 2; iPID++)
	for (int iPID = 0; iPID < nPID; iPID++)
	{
		EdbPattern *pattern = pvr->GetPattern(iPID);
		if (pattern == NULL)
			continue;
		plate = pattern->Plate();
		for (int itrk = 0; itrk < ntrk; itrk++)
		{
			EdbTrackP *t = pvr->GetTrack(itrk);
			int count = 0;
			int countForLSM = 0;
			double arrayX[5], arrayY[5], arrayZ[5];
			double arrayXForLSM[4], arrayYForLSM[4], arrayZForLSM[4];
			double probedSegmentX, probedSegmentY, probedSegmentZ;
			double probedSegmentTX, probedSegmentTY;
			nseg = t->N();
			for (int iseg = 0; iseg < nseg; iseg++)
			{
				EdbSegP *s = t->GetSegment(iseg);
				int sPID = s->PID();
				if(iPID==0) //first plate
				{
					if(sPID>iPID+4)
						break;
				}
				if(iPID==1) //second plate
				{
					if(sPID>iPID+3)
						break;
				}
				if(iPID>=2&&iPID<nPID-2) //mid plate
				{
					if (sPID > iPID + 2)
						break;
					if (sPID < iPID - 2)
						continue;
				}
				if(iPID==nPID-2) //second plate from last plate
				{
					if(sPID<iPID-3)
						continue;
				}
				if(iPID==nPID-1) //last plate
				{
					if(sPID<iPID-4)
						continue;
				}
				arrayX[count] = s->X();
				arrayY[count] = s->Y();
				arrayZ[count] = s->Z();
				count++;
				if (sPID == iPID)
				{
					probedSegmentX= s->X();
					probedSegmentY= s->Y();
					probedSegmentZ= s->Z();
					probedSegmentTX = s->TX();
					probedSegmentTY = s->TY();
				}else{
					arrayXForLSM[countForLSM] = s->X();
					arrayYForLSM[countForLSM] = s->Y();
					arrayZForLSM[countForLSM] = s->Z();
					countForLSM++;
				}
				if (count == 5)
					break;
			}
			if (count != 5)
				continue;
			int areaX[5];
			int areaY[5];
			for (int ipoint = 0; ipoint < 5; ipoint++)
			{
				areaX[ipoint] = (arrayX[ipoint] - (Xcenter - XYrange)) / bin_width;
				areaY[ipoint] = (arrayY[ipoint] - (Ycenter - XYrange)) / bin_width;
			}
			int cross_the_line = 0;
			for (int ipoint = 0; ipoint < 5 - 1; ipoint++)
			{
				if (areaX[ipoint] != areaX[ipoint + 1] || areaY[ipoint] != areaY[ipoint + 1])
				{
					cross_the_line = 1;
					break;
				}
			}
			double a0, slopeX, slopeY;
			CalcLSM(arrayZForLSM, arrayXForLSM, 4, a0, slopeX);
			double x3fit = a0 + slopeX * probedSegmentZ;
			CalcLSM(arrayZForLSM, arrayYForLSM, 4, a0, slopeY);
			double y3fit = a0 + slopeY * probedSegmentZ;

			// Calculate delta X and delta Y.
			deltaXV->push_back(probedSegmentX - x3fit);
			deltaYV->push_back(probedSegmentY - y3fit);
			deltaTXV->push_back(probedSegmentTX - slopeX);
			deltaTYV->push_back(probedSegmentTY - slopeY);
			xV->push_back(t->X());
			yV->push_back(t->Y());
			slopeXV->push_back(slopeX);
			slopeYV->push_back(slopeY);
			crossTheLineV->push_back(cross_the_line);
			tridV->push_back(t->ID());
			nsegV->push_back(nseg);
		}
		deltaXY->Fill();
		deltaXV->clear();
		deltaYV->clear();
		deltaTXV->clear();
		deltaTYV->clear();
		xV->clear();
		yV->clear();
		slopeXV->clear();
		slopeYV->clear();
		crossTheLineV->clear();
		tridV->clear();
		nsegV->clear();
	}
}

void FnuQualityCheck::FitDeltaXY()
{
	// Create histograms of delta x and y and fit them.
	// Before using this method, do CalcDeltaXY() to make deltaXY data.
	posResPar = new TTree("posResPar", "posResPar");

	posResPar->Branch("sigmaX", &sigmaX);
	posResPar->Branch("sigmaY", &sigmaY);
	posResPar->Branch("meanX", &meanX);
	posResPar->Branch("meanY", &meanY);
	posResPar->Branch("entries", &entries);
	posResPar->Branch("plate", &plate);

	hdeltaX = new TH1D("hdeltaX", "hdeltaX", 100, -2, 2);
	hdeltaY = new TH1D("hdeltaY", "hdeltaY", 100, -2, 2);
	htree = new TTree("htree", "htree");
	htree->Branch("hdeltaX", &hdeltaX);
	htree->Branch("hdeltaY", &hdeltaY);
	htree->Branch("plate", &plate);

	TF1 *f = new TF1("gaus", "gaus", -2, 2);
	f->SetParLimits(5, 0, 0.4);
	gStyle->SetOptFit();

	// int plMin = deltaXY->GetMinimum("pl");
	// int plMax = deltaXY->GetMaximum("pl");
	TBranch *deltaXB = deltaXY->GetBranch("deltaX");
	TBranch *deltaYB = deltaXY->GetBranch("deltaY");
	TBranch *slopeXB = deltaXY->GetBranch("slopeX");
	TBranch *slopeYB = deltaXY->GetBranch("slopeY");
	TBranch *plateB = deltaXY->GetBranch("plate");
	for (int ient = 0; ient < deltaXY->GetEntriesFast(); ient++)
	{
		deltaXB->GetEntry(ient);
		deltaYB->GetEntry(ient);
		slopeXB->GetEntry(ient);
		slopeYB->GetEntry(ient);
		plateB->GetEntry(ient);

		const auto slopeXMean = std::accumulate(slopeXV->begin(), slopeXV->end(), 0.0) / slopeXV->size();
		const auto slopeYMean = std::accumulate(slopeYV->begin(), slopeYV->end(), 0.0) / slopeYV->size();
		for (int i = 0; i < deltaXV->size(); i++)
		{
			if (fabs(deltaYV->at(i)) <= 2 && fabs(deltaXV->at(i)) <= 2 && fabs(slopeXV->at(i) - slopeXMean) < angleCut && fabs(slopeYV->at(i) - slopeYMean) < angleCut)
			{
				hdeltaX->Fill(deltaXV->at(i));
				hdeltaY->Fill(deltaYV->at(i));
			}
		}
		hdeltaX->SetTitle(Form("pl%d %s;deltaX (#mum);", plate, title.Data()));
		hdeltaY->SetTitle(Form("pl%d %s;deltaY (#mum);", plate, title.Data()));
		hdeltaX->Draw();
		f->SetParameters(1000, 0, 0.2);
		meanX = hdeltaX->GetMean();
		double RMSX = hdeltaX->GetRMS();
		hdeltaX->Fit(f, "Q", "", meanX - RMSX, meanX + RMSX);
		entries = hdeltaX->GetEntries();
		sigmaX = f->GetParameter(2);
		hdeltaY->Draw();
		f->SetParameters(1000, 0, 0.2);
		meanY = hdeltaY->GetMean();
		double RMSY = hdeltaY->GetRMS();
		hdeltaY->Fit(f, "Q", "", meanY - RMSY, meanY + RMSY);
		sigmaY = f->GetParameter(2);
		htree->Fill();
		posResPar->Fill();
		hdeltaX->Reset();
		hdeltaY->Reset();
	}
}

void FnuQualityCheck::FitDeltaXYAllPlatesTogether()
{
	// Create histograms of delta x and y and fit them.
	// Before using this method, do CalcDeltaXY() to make deltaXY data.
	posResPar = new TTree("posResPar", "posResPar");

	posResPar->Branch("sigmaX", &sigmaX);
	posResPar->Branch("sigmaY", &sigmaY);
	posResPar->Branch("meanX", &meanX);
	posResPar->Branch("meanY", &meanY);
	posResPar->Branch("entries", &entries);
	posResPar->Branch("plate", &plate);

	hdeltaX = new TH1D("hdeltaX", "hdeltaX", 100, -2, 2);
	hdeltaY = new TH1D("hdeltaY", "hdeltaY", 100, -2, 2);
	htree = new TTree("htree", "htree");
	htree->Branch("hdeltaX", &hdeltaX);
	htree->Branch("hdeltaY", &hdeltaY);
	htree->Branch("plate", &plate);

	TF1 *f = new TF1("gaus", "gaus", -2, 2);
	f->SetParLimits(5, 0, 0.4);
	gStyle->SetOptFit();

	// int plMin = deltaXY->GetMinimum("pl");
	// int plMax = deltaXY->GetMaximum("pl");
	TBranch *deltaXB = deltaXY->GetBranch("deltaX");
	TBranch *deltaYB = deltaXY->GetBranch("deltaY");
	TBranch *slopeXB = deltaXY->GetBranch("slopeX");
	TBranch *slopeYB = deltaXY->GetBranch("slopeY");
	TBranch *plateB = deltaXY->GetBranch("plate");
	for (int ient = 0; ient < deltaXY->GetEntriesFast(); ient++)
	{
		deltaXB->GetEntry(ient);
		deltaYB->GetEntry(ient);
		slopeXB->GetEntry(ient);
		slopeYB->GetEntry(ient);
		plateB->GetEntry(ient);

		const auto slopeXMean = std::accumulate(slopeXV->begin(), slopeXV->end(), 0.0) / slopeXV->size();
		const auto slopeYMean = std::accumulate(slopeYV->begin(), slopeYV->end(), 0.0) / slopeYV->size();
		for (int i = 0; i < deltaXV->size(); i++)
		{
			if (fabs(deltaYV->at(i)) <= 2 && fabs(deltaXV->at(i)) <= 2 && fabs(slopeXV->at(i) - slopeXMean) < angleCut && fabs(slopeYV->at(i) - slopeYMean) < angleCut)
			{
				hdeltaX->Fill(deltaXV->at(i));
				hdeltaY->Fill(deltaYV->at(i));
			}
		}
	}
	hdeltaX->SetTitle(Form("pl%d %s;deltaX (#mum);", plate, title.Data()));
	hdeltaY->SetTitle(Form("pl%d %s;deltaY (#mum);", plate, title.Data()));
	hdeltaX->Draw();
	f->SetParameters(1000, 0, 0.2);
	meanX = hdeltaX->GetMean();
	double RMSX = hdeltaX->GetRMS();
	hdeltaX->Fit(f, "Q", "", meanX - RMSX, meanX + RMSX);
	entries = hdeltaX->GetEntries();
	sigmaX = f->GetParameter(2);
	hdeltaY->Draw();
	f->SetParameters(1000, 0, 0.2);
	meanY = hdeltaY->GetMean();
	double RMSY = hdeltaY->GetRMS();
	hdeltaY->Fit(f, "Q", "", meanY - RMSY, meanY + RMSY);
	sigmaY = f->GetParameter(2);
	htree->Fill();
	posResPar->Fill();
}

void FnuQualityCheck::CalcLSM(double x[], double y[], int N, double &a0, double &a1)
{
	// y = a0 + a1*x
	int i;
	double A00 = 0, A01 = 0, A02 = 0, A11 = 0, A12 = 0;

	for (i = 0; i < N; i++)
	{
		A00 += 1.0;
		A01 += x[i];
		A02 += y[i];
		A11 += x[i] * x[i];
		A12 += x[i] * y[i];
	}
	a0 = (A02 * A11 - A01 * A12) / (A00 * A11 - A01 * A01);
	a1 = (A00 * A12 - A01 * A02) / (A00 * A11 - A01 * A01);
}
void FnuQualityCheck::MakePosResGraphHist()
{
	// Make graphs and histograms about position resolution.
	// This makes graphs of mean:plate, histograms of position resolution and graphs of position resolution : plate.

	// meanX and meanY : plate
	posResPar->Draw("meanX:plate");
	int N = posResPar->GetSelectedRows();
	meanXGraph = new TGraph(N, posResPar->GetV2(), posResPar->GetV1());
	meanXGraph->SetMarkerStyle(20);
	meanXGraph->SetMarkerColor(kRed);
	meanXGraph->SetNameTitle("meanYGraph", "mean Y (" + title + ");plate;mean (#mum)");
	posResPar->Draw("meanY:plate");
	N = posResPar->GetSelectedRows();
	meanYGraph = new TGraph(N, posResPar->GetV2(), posResPar->GetV1());
	meanYGraph->SetMarkerStyle(20);
	meanYGraph->SetMarkerColor(kBlue);
	meanYGraph->SetNameTitle("meanXGraph", "mean X (" + title + ");plate;mean (#mum)");

	// sigmaX and sigmaY
	sigmaXHist = new TH1D("sigmaXHist", "position resolution X (" + title + ");position resolution (#mum)", 100, 0, 1);
	posResPar->Draw("sigmaX>>sigmaXHist", "abs(meanX)<100&&entries>0");
	sigmaYHist = new TH1D("sigmaYHist", "position resolution Y (" + title + ");position resolution (#mum)", 100, 0, 1);
	posResPar->Draw("sigmaY>>sigmaYHist", "abs(meanY)<100&&entries>0");

	// sigmaX and sigmaY:pl
	posResPar->Draw("sigmaX:plate");
	N = posResPar->GetSelectedRows();
	sigmaXGraph = new TGraph(N, posResPar->GetV2(), posResPar->GetV1());
	sigmaXGraph->SetMarkerStyle(20);
	sigmaXGraph->SetMarkerColor(kRed);
	sigmaXGraph->SetNameTitle("sigmaXGraph", "position resolution X (" + title + ");plate;position resolution (#mum)");
	posResPar->Draw("sigmaY:plate");
	N = posResPar->GetSelectedRows();
	sigmaYGraph = new TGraph(N, posResPar->GetV2(), posResPar->GetV1());
	sigmaYGraph->SetMarkerStyle(20);
	sigmaYGraph->SetMarkerColor(kBlue);
	sigmaYGraph->SetNameTitle("sigmaYGraph", "position resolution Y (" + title + ");plate;position resolution (#mum)");
}

void FnuQualityCheck::WritePosResGraphHist(TString filename)
{
	// write graphs and histograms about position resolution to file.
	TFile fout(filename, "recreate");
	meanXGraph->Write();
	meanYGraph->Write();
	sigmaXHist->Write();
	sigmaYHist->Write();
	sigmaXGraph->Write();
	sigmaYGraph->Write();
	fout.Close();
}
void FnuQualityCheck::PrintPosResGraphHist(TString filename)
{
	// int plMax = posResPar->GetMaximum("pl");
	// int plMin = posResPar->GetMinimum("pl");
	TCanvas *c1 = new TCanvas();
	c1->Print(filename + "[");

	// meanX and meanY : plate
	c1->SetGridx(1);
	TMultiGraph *mg = new TMultiGraph("mg", "mean:plate " + title + ";plate;mean (#mum)");
	mg->Add(meanXGraph);
	mg->Add(meanYGraph);
	mg->Draw("ap");
	TLegend *leg = new TLegend(0.9, 0.8, 1.0, 0.9);
	leg->AddEntry(meanXGraph, "meanX", "p");
	leg->AddEntry(meanYGraph, "meanY", "p");
	leg->Draw();
	c1->Print(filename);
	c1->SetGridx(0);

	// sigmaX and sigmaY
	TList *l = new TList;
	l->Add(sigmaXHist);
	l->Add(sigmaYHist);
	TH1F *resolution = new TH1F("resolution", "position resolution (" + title + ");position resolution (#mum)", 100, 0, 1);
	resolution->Merge(l);
	resolution->Draw();
	c1->Print(filename);

	// sigmaX and sigmaY:pl
	c1->SetGridx(1);
	TMultiGraph *mg2 = new TMultiGraph("mg2", "sigma:plate (" + title + ");plate;sigma (#mum)");
	mg2->Add(sigmaXGraph);
	mg2->Add(sigmaYGraph);
	mg2->Draw("ap");
	TLegend *leg2 = new TLegend(0.9, 0.8, 1.0, 0.9);
	leg2->AddEntry(sigmaXGraph, "resolutionX", "p");
	leg2->AddEntry(sigmaYGraph, "resolutionY", "p");
	leg2->Draw();
	c1->Print(filename);
	c1->Clear();
	c1->Print(filename);
	c1->Print(filename + "]");
}
void FnuQualityCheck::PrintDeltaXYHist(TString filename)
{
	TCanvas c;
	// c.Print("pos_res/deltaxy_" + title + ".pdf[");
	c.Print(filename + "[");
	for (int ient = 0; ient < htree->GetEntriesFast(); ient++)
	{
		htree->GetEntry(ient);
		hdeltaX->Draw();
		// c.Print("pos_res/deltaxy_" + title + ".pdf");
		c.Print(filename);
		hdeltaY->Draw();
		// c.Print("pos_res/deltaxy_" + title + ".pdf");
		c.Print(filename);
		printf("Histograms for plate %d have been printed\n", plate);
	}
	// c.Print("pos_res/deltaxy_" + title + ".pdf]");
	c.Print(filename + "]");
}
void FnuQualityCheck::WritePosResPar(TString filename)
{
	// TFile fout("pos_res/sigmaPar_" + title + ".root", "recreate");
	TFile fout(filename, "recreate");
	posResPar->Write();
	fout.Close();
}
void FnuQualityCheck::WriteDeltaXY(TString filename)
{
	// TFile fout("deltaXY/tree_" + title + ".root", "recreate");
	TFile fout(filename, "recreate");
	deltaXY->Write();
	fout.Close();
}
void FnuQualityCheck::CalcMeanDeltaXY(double Xcenter, double Ycenter, double cellSize)
{
	// calculate the mean of the deltaX or deltaY in each area (cellSize*cellSize) for drawing the arrow of positional displacement
	// Xcenter and Ycenter are the center positions in the volume
	int nbinsx, nbinsy;
	nbinsx = nbinsy = ceil(XYrange * 2 / cellSize);
	deltaX2DHist = new TH2D("deltaX2DHist", "deltaX mean (" + title + ");x (#mum);y (#mum);deltaX mean (#mum)", nbinsx, Xcenter - XYrange, Xcenter - XYrange + cellSize * nbinsx, nbinsy, Ycenter - XYrange, Ycenter - XYrange + cellSize * nbinsy);
	deltaY2DHist = new TH2D("deltaY2DHist", "deltaY mean (" + title + ");x (#mum);y (#mum);deltaY mean (#mum)", nbinsx, Xcenter - XYrange, Xcenter - XYrange + cellSize * nbinsx, nbinsy, Ycenter - XYrange, Ycenter - XYrange + cellSize * nbinsy);
	entries2DHist = new TH2D("entries2DHist", "entries (" + title + ");x (#mum);y (#mum);entries", nbinsx, Xcenter - XYrange, Xcenter - XYrange + cellSize * nbinsx, nbinsy, Ycenter - XYrange, Ycenter - XYrange + cellSize * nbinsy);
	meanDeltaXY = new TTree("meanDeltaXY", "2D histograms of mean deltaXY for each plate");
	meanDeltaXY->Branch("deltaX2DHist", &deltaX2DHist);
	meanDeltaXY->Branch("deltaY2DHist", &deltaY2DHist);
	meanDeltaXY->Branch("entries2DHist", &entries2DHist);
	meanDeltaXY->Branch("plate", &plate);

	TBranch *deltaXB = deltaXY->GetBranch("deltaX");
	TBranch *deltaYB = deltaXY->GetBranch("deltaY");
	TBranch *xB = deltaXY->GetBranch("x");
	TBranch *yB = deltaXY->GetBranch("y");
	TBranch *slopeXB = deltaXY->GetBranch("slopeX");
	TBranch *slopeYB = deltaXY->GetBranch("slopeY");
	TBranch *plateB = deltaXY->GetBranch("plate");
	for (int ient = 0; ient < deltaXY->GetEntriesFast(); ient++)
	{
		deltaXB->GetEntry(ient);
		deltaYB->GetEntry(ient);
		xB->GetEntry(ient);
		yB->GetEntry(ient);
		slopeXB->GetEntry(ient);
		slopeYB->GetEntry(ient);
		plateB->GetEntry(ient);

		const auto slopeXMean = std::accumulate(slopeXV->begin(), slopeXV->end(), 0.0) / slopeXV->size();
		const auto slopeYMean = std::accumulate(slopeYV->begin(), slopeYV->end(), 0.0) / slopeYV->size();
		for (int i = 0; i < deltaXV->size(); i++)
		{
			if (fabs(deltaYV->at(i)) <= 2 && fabs(deltaXV->at(i)) <= 2 && fabs(slopeXV->at(i) - slopeXMean) < angleCut && fabs(slopeYV->at(i) - slopeYMean) < angleCut)
			{
				deltaX2DHist->Fill(xV->at(i), yV->at(i), deltaXV->at(i));
				deltaY2DHist->Fill(xV->at(i), yV->at(i), deltaYV->at(i));
				entries2DHist->Fill(xV->at(i), yV->at(i));
			}
		}
		deltaX2DHist->Divide(entries2DHist);
		deltaY2DHist->Divide(entries2DHist);
		meanDeltaXY->Fill();

		deltaX2DHist->Reset();
		deltaY2DHist->Reset();
		entries2DHist->Reset();
	}
}

void FnuQualityCheck::WriteMeanDeltaXY(TString filename)
{
	TFile fout(filename, "recreate");
	meanDeltaXY->Write();
	fout.Close();
}

int FnuQualityCheck::PrintMeanDeltaXYArrowPlot(TString filename)
{
	int nent = meanDeltaXY->GetEntriesFast();
	if (0 == nent)
	{
		printf("Do CalcMeanDeltaXY(double Xcenter, double Ycenter, double cellSize) before PrintMeanDeltaXYArrowPlot(TString filename)\n");
		return nent;
	}
	TArrow *arrow = new TArrow();
	double minXaxis = deltaX2DHist->GetXaxis()->GetXmin();
	double minYaxis = deltaX2DHist->GetYaxis()->GetXmin();
	double maxXaxis = deltaX2DHist->GetXaxis()->GetXmax();
	double maxYaxis = deltaX2DHist->GetYaxis()->GetXmax();
	double rangeXaxis = maxXaxis - minXaxis;
	double rangeYaxis = maxYaxis - minYaxis;
	int nbinsX = deltaX2DHist->GetXaxis()->GetNbins();
	int nbinsY = deltaX2DHist->GetYaxis()->GetNbins();
	TLatex *l = new TLatex;
	l->SetTextAlign(33);
	l->SetTextSize(0.03);
	l->SetTextFont(42);
	int scale = 2000; // Should not be hard coded...
	TCanvas c;
	c.Print(filename + "[");
	for (int ient = 0; ient < nent; ient++)
	{
		meanDeltaXY->GetEntry(ient);
		TH1F *fr = gPad->DrawFrame(minXaxis, minYaxis, maxXaxis, maxYaxis, Form("deltaXY pl%d ;x(#mum);y(#mum)", plate));
		arrow->DrawArrow(maxXaxis + rangeXaxis * 0.05, maxYaxis + rangeYaxis * 0.1, maxXaxis + rangeXaxis * 0.05 - scale * 0.5, maxYaxis + rangeYaxis * 0.1, 0.008, ">"); // equivalent to 0.5 μm.
		l->DrawLatex(maxXaxis + rangeXaxis * 0.05, maxYaxis + rangeYaxis * 0.07, "0.5 #mum");
		for (int ibinX = 1; ibinX <= nbinsX; ibinX++)
		{
			for (int ibinY = 1; ibinY <= nbinsY; ibinY++)
			{
				if (entries2DHist->GetBinContent(ibinX, ibinY) == 0)
				{
					continue;
				}
				double deltaXMean = deltaX2DHist->GetBinContent(ibinX, ibinY);
				double deltaYMean = deltaY2DHist->GetBinContent(ibinX, ibinY);
				double binCenterX = deltaX2DHist->GetXaxis()->GetBinCenter(ibinX);
				double binCenterY = deltaX2DHist->GetYaxis()->GetBinCenter(ibinY);
				arrow->DrawArrow(binCenterX, binCenterY, binCenterX + scale * deltaXMean, binCenterY + scale * deltaYMean, 0.008, ">");
			}
		}
		c.Print(filename);
	}
	c.Print(filename + "]");
	return nent;
}

TTree *FnuQualityCheck::MakeHistDeltaTXY(TTree *deltaXY)
{
	hdeltaTX = new TH1D("hdeltaTX", "hdeltaTX", 100, -0.02, 0.02);
	hdeltaTY = new TH1D("hdeltaTY", "hdeltaTY", 100, -0.02, 0.02);
	TTree *treeHistDeltaTXY = new TTree("treeHistDeltaTXY","tree for histograms of delta TX and TY");
	treeHistDeltaTXY->Branch("hdeltaTX", &hdeltaTX);
	treeHistDeltaTXY->Branch("hdeltaTY", &hdeltaTY);
	treeHistDeltaTXY->Branch("plate", &plate);

	deltaXY->SetBranchAddress("deltaTX",&deltaTXV);
	deltaXY->SetBranchAddress("deltaTY",&deltaTYV);
	deltaXY->SetBranchAddress("slopeX",&slopeXV);
	deltaXY->SetBranchAddress("slopeY",&slopeYV);
	deltaXY->SetBranchAddress("plate",&plate);
	TBranch *deltaTXB = deltaXY->GetBranch("deltaTX");
	TBranch *deltaTYB = deltaXY->GetBranch("deltaTY");
	TBranch *slopeXB = deltaXY->GetBranch("slopeX");
	TBranch *slopeYB = deltaXY->GetBranch("slopeY");
	TBranch *plateB = deltaXY->GetBranch("plate");

	for (int ient = 0; ient < deltaXY->GetEntriesFast(); ient++)
	{
		deltaTXB->GetEntry(ient);
		deltaTYB->GetEntry(ient);
		slopeXB->GetEntry(ient);
		slopeYB->GetEntry(ient);
		plateB->GetEntry(ient);
		const auto slopeXMean = std::accumulate(slopeXV->begin(), slopeXV->end(), 0.0) / slopeXV->size();
		const auto slopeYMean = std::accumulate(slopeYV->begin(), slopeYV->end(), 0.0) / slopeYV->size();
		for (int i = 0; i < deltaTXV->size(); i++)
		{
			if (fabs(deltaTYV->at(i)) <= 0.02 && fabs(deltaTXV->at(i)) <= 0.02 && fabs(slopeXV->at(i) - slopeXMean) < angleCut && fabs(slopeYV->at(i) - slopeYMean) < angleCut)
			{
				hdeltaTX->Fill(deltaTXV->at(i));
				hdeltaTY->Fill(deltaTYV->at(i));
			}
		}
		hdeltaTX->SetTitle(Form("pl%d %s;deltaX (#mum);", plate, title.Data()));
		hdeltaTY->SetTitle(Form("pl%d %s;deltaY (#mum);", plate, title.Data()));
		treeHistDeltaTXY->Fill();
		hdeltaTX->Reset();
		hdeltaTY->Reset();
	}
	return treeHistDeltaTXY;
}
void FnuQualityCheck::PrintDeltaTXYHist(TTree *treeHistDeltaTXY,TString filename)
{
	TCanvas c;
	// c.Print("pos_res/deltaxy_" + title + ".pdf[");
	c.Print(filename + "[");
	for (int ient = 0; ient < treeHistDeltaTXY->GetEntriesFast(); ient++)
	{
		treeHistDeltaTXY->GetEntry(ient);
		hdeltaTX->Draw();
		// c.Print("pos_res/deltaxy_" + title + ".pdf");
		c.Print(filename);
		hdeltaTY->Draw();
		// c.Print("pos_res/deltaxy_" + title + ".pdf");
		c.Print(filename);
		// printf("Histograms for plate %d have been printed\n", plate);
	}
	// c.Print("pos_res/deltaxy_" + title + ".pdf]");
	c.Print(filename + "]");
}
TTree *FnuQualityCheck::FitHistDeltaTXY(TTree *treeHistDeltaTXY)
{
	// Create histograms of delta tx and ty and fit them.
	TTree *angleResolutionPar = new TTree("angleResolutionPar","parameters for angle resolution");
	angleResolutionPar->Branch("sigmaTX", &sigmaTX);
	angleResolutionPar->Branch("sigmaTY", &sigmaTY);
	angleResolutionPar->Branch("meanTX", &meanTX);
	angleResolutionPar->Branch("meanTY", &meanTY);
	angleResolutionPar->Branch("entries", &entries);
	angleResolutionPar->Branch("plate", &plate);
	treeHistDeltaTXY->SetBranchAddress("hdeltaTX",&hdeltaTX);
	treeHistDeltaTXY->SetBranchAddress("hdeltaTY",&hdeltaTY);
	TF1 *f = new TF1("gaus", "gaus", -0.02, 0.02);
	// f->SetParLimits(5, 0, 0.4);
	gStyle->SetOptFit();
	for (int ient = 0; ient < treeHistDeltaTXY->GetEntriesFast(); ient++)
	{
		treeHistDeltaTXY->GetEntry(ient);
		hdeltaTX->Draw();
		f->SetParameters(1000, 0, 0.2);
		meanTX = hdeltaTX->GetMean();
		double RMSX = hdeltaTX->GetRMS();
		hdeltaTX->Fit(f, "Q", "", meanTX - RMSX, meanTX + RMSX);
		entries = hdeltaTX->GetEntries();
		sigmaTX = f->GetParameter(2);
		hdeltaTY->Draw();
		f->SetParameters(1000, 0, 0.2);
		meanTY = hdeltaTY->GetMean();
		double RMSY = hdeltaTY->GetRMS();
		hdeltaTY->Fit(f, "Q", "", meanTY - RMSY, meanTY + RMSY);
		entries = hdeltaTY->GetEntries();
		sigmaTY = f->GetParameter(2);
		angleResolutionPar->Fill();
	}
	return angleResolutionPar;
}
void FnuQualityCheck::MakeGraphHistAngleResolution(TTree *angleResolutionPar)
{
	// Make graphs and histograms about angle resolution.
	// This makes graphs of mean:plate, histograms of angle resolution and graphs of angle resolution : plate.

	// meanTX and meanTY : plate
	angleResolutionPar->Draw("meanTX:plate");
	int N = angleResolutionPar->GetSelectedRows();
	meanTXGraph = new TGraph(N, angleResolutionPar->GetV2(), angleResolutionPar->GetV1());
	meanTXGraph->SetMarkerStyle(20);
	meanTXGraph->SetMarkerColor(kRed);
	meanTXGraph->SetNameTitle("meanTYGraph", "mean TY (" + title + ");plate;mean (tan#theta)");
	angleResolutionPar->Draw("meanTY:plate");
	N = angleResolutionPar->GetSelectedRows();
	meanTYGraph = new TGraph(N, angleResolutionPar->GetV2(), angleResolutionPar->GetV1());
	meanTYGraph->SetMarkerStyle(20);
	meanTYGraph->SetMarkerColor(kBlue);
	meanTYGraph->SetNameTitle("meanTXGraph", "mean TX (" + title + ");plate;mean (tan#theta)");

	// sigmaTX and sigmaTY
	sigmaTXHist = new TH1D("sigmaTXHist", "angular resolution X (" + title + ");angular resolution (tan#theta)", 100, 0, 0.005);
	angleResolutionPar->Draw("sigmaTX>>sigmaTXHist", "abs(meanTX)<0.005&&entries>0");
	sigmaTYHist = new TH1D("sigmaTYHist", "angular resolution Y (" + title + ");angular resolution (tan#theta)", 100, 0, 0.005);
	angleResolutionPar->Draw("sigmaTY>>sigmaTYHist", "abs(meanTY)<0.005&&entries>0");

	// sigmaTX and sigmaTY:pl
	angleResolutionPar->Draw("sigmaTX:plate");
	N = angleResolutionPar->GetSelectedRows();
	sigmaTXGraph = new TGraph(N, angleResolutionPar->GetV2(), angleResolutionPar->GetV1());
	sigmaTXGraph->SetMarkerStyle(20);
	sigmaTXGraph->SetMarkerColor(kRed);
	sigmaTXGraph->SetNameTitle("sigmaTXGraph", "angular resolution X (" + title + ");plate;angular resolution (tan#theta)");
	angleResolutionPar->Draw("sigmaTY:plate");
	N = angleResolutionPar->GetSelectedRows();
	sigmaTYGraph = new TGraph(N, angleResolutionPar->GetV2(), angleResolutionPar->GetV1());
	sigmaTYGraph->SetMarkerStyle(20);
	sigmaTYGraph->SetMarkerColor(kBlue);
	sigmaTYGraph->SetNameTitle("sigmaTYGraph", "angular resolution Y (" + title + ");plate;angular resolution (tan#theta)");
}
void FnuQualityCheck::PrintGraphHistAngleResolution(TString filename)
{
	// int plMax = posResPar->GetMaximum("pl");
	// int plMin = posResPar->GetMinimum("pl");
	TCanvas *c1 = new TCanvas();
	c1->Print(filename + "[");

	// meanTX and meanTY : plate
	c1->SetGridx(1);
	TMultiGraph *mg = new TMultiGraph("mg", "mean:plate " + title + ";plate;mean (tan#theta)");
	mg->Add(meanTXGraph);
	mg->Add(meanTYGraph);
	mg->Draw("ap");
	TLegend *leg = new TLegend(0.9, 0.8, 1.0, 0.9);
	leg->AddEntry(meanTXGraph, "meanTX", "p");
	leg->AddEntry(meanTYGraph, "meanTY", "p");
	leg->Draw();
	c1->Print(filename);
	c1->SetGridx(0);

	// sigmaTX and sigmaTY
	TList *l = new TList;
	l->Add(sigmaTXHist);
	l->Add(sigmaTYHist);
	TH1F *resolution = new TH1F("resolution", "angular resolution (" + title + ");angular resolution (tan#theta)", 100, 0, 0.005);
	resolution->Merge(l);
	resolution->Draw();
	c1->Print(filename);

	// sigmaTX and sigmaTY:pl
	c1->SetGridx(1);
	TMultiGraph *mg2 = new TMultiGraph("mg2", "sigma:plate (" + title + ");plate;sigma (tan#theta)");
	mg2->Add(sigmaTXGraph);
	mg2->Add(sigmaTYGraph);
	mg2->Draw("ap");
	TLegend *leg2 = new TLegend(0.9, 0.8, 1.0, 0.9);
	leg2->AddEntry(sigmaTXGraph, "resolutionX", "p");
	leg2->AddEntry(sigmaTYGraph, "resolutionY", "p");
	leg2->Draw();
	c1->Print(filename);
	c1->Clear();
	c1->Print(filename);
	c1->Print(filename + "]");
}

void FnuQualityCheck::CalcEfficiency()
{
	double bins_angle[bins_vec_angle.size()];
	std::copy(bins_vec_angle.begin(), bins_vec_angle.end(), bins_angle);
	int nbins_angle = bins_vec_angle.size() - 1;
	double bins_TXTY[bins_vec_TXTY.size()];
	std::copy(bins_vec_TXTY.begin(), bins_vec_TXTY.end(), bins_TXTY);
	int nbins_TXTY = bins_vec_TXTY.size() - 1;

	eachAngleEfficiency = new TEfficiency("Eff_angle", Form("Efficiency for each angle (%s);tan#theta;efficiency", title.Data()), nbins_angle, bins_angle);
	eachPlateEfficiency = new TEfficiency("Eff_plate", Form("Efficiency for each plate (%s);plate;efficiency", title.Data()), plMax - plMin + 1, plMin - 0.5, plMax + 0.5);
	eachTXEfficiency = new TEfficiency("Eff_TX", Form("Efficiency for each TX (%s);tan#theta;efficiency", title.Data()), nbins_TXTY, bins_TXTY);
	eachTYEfficiency = new TEfficiency("Eff_TY", Form("Efficiency for each TY (%s);tan#theta;efficiency", title.Data()), nbins_TXTY, bins_TXTY);

	effInfo = new TTree("effInfo", "efficiency infomation");
	effInfo->Branch("trackID", &trackID);
	effInfo->Branch("x", &x);
	effInfo->Branch("y", &y);
	effInfo->Branch("angle", &angle);
	effInfo->Branch("TX", &TX);
	effInfo->Branch("TY", &TY);
	effInfo->Branch("plate", &plate);
	effInfo->Branch("nseg", &nseg);
	effInfo->Branch("W", &W);
	effInfo->Branch("hitsOnThePlate", &hitsOnThePlate);

	double x1, y1, z1, x2, y2, z2;
	for (int itrk = 0; itrk < ntrk; itrk++)
	{
		EdbTrackP *t = pvr->GetTrack(itrk);
		int nseg = t->N();

		for (int iPID = 0; iPID < nPID; iPID++)
		{
			EdbPattern *pattern = pvr->GetPattern(iPID);
			if (pattern == NULL)
				continue;
			int iplate = pattern->Plate();
			int counts = 0;
			hitsOnThePlate = 0;
			W = 0;
			for (int iseg = 0; iseg < nseg; iseg++)
			{
				EdbSegP *s = t->GetSegment(iseg);
				int sPID = s->PID();
				if (iPID == 0) // first plate
				{
					if (sPID > iPID + 4)
						break;
					if (sPID == iPID + 4 || sPID == iPID + 3)
						counts++;
					if (sPID == iPID + 1)
					{
						x1 = s->X();
						y1 = s->Y();
						z1 = s->Z();
						counts++;
					}
					if (sPID == iPID + 2)
					{
						x2 = s->X();
						y2 = s->Y();
						z2 = s->Z();
						counts++;
					}
				}
				if (iPID == 1) // second plate
				{
					if (sPID > iPID + 3)
						break;
					if (sPID == iPID + 3 || sPID == iPID + 2)
						counts++;
					if (sPID == iPID - 1)
					{
						x1 = s->X();
						y1 = s->Y();
						z1 = s->Z();
						counts++;
					}
					if (sPID == iPID + 1)
					{
						x2 = s->X();
						y2 = s->Y();
						z2 = s->Z();
						counts++;
					}
				}
				if (iPID >= 2 && iPID < nPID - 2) // mid plate
				{
					if (sPID > iPID + 2)
						break;
					if (sPID < iPID - 2)
						continue;
					if (sPID == iPID - 2 || sPID == iPID + 2)
						counts++;
					if (sPID == iPID - 1)
					{
						x1 = s->X();
						y1 = s->Y();
						z1 = s->Z();
						counts++;
					}
					if (sPID == iPID + 1)
					{
						x2 = s->X();
						y2 = s->Y();
						z2 = s->Z();
						counts++;
					}
				}
				if (iPID == nPID - 2) // second from last plate
				{
					if (sPID < iPID - 3)
						continue;
					if (sPID == iPID - 3 || sPID == iPID - 2)
						counts++;
					if (sPID == iPID - 1)
					{
						x1 = s->X();
						y1 = s->Y();
						z1 = s->Z();
						counts++;
					}
					if (sPID == iPID + 1)
					{
						x2 = s->X();
						y2 = s->Y();
						z2 = s->Z();
						counts++;
					}
				}
				if (iPID == nPID - 1) // last plate
				{
					if (sPID < iPID - 4)
						continue;
					if (sPID == iPID - 4 || sPID == iPID - 3)
						counts++;
					if (sPID == iPID - 2)
					{
						x1 = s->X();
						y1 = s->Y();
						z1 = s->Z();
						counts++;
					}
					if (sPID == iPID - 1)
					{
						x2 = s->X();
						y2 = s->Y();
						z2 = s->Z();
						counts++;
					}
				}
				if (sPID == iPID)
				{
					hitsOnThePlate = 1;
					W = s->W();
				}
			}
			if (counts == 4)
			{
				TX = (x2 - x1) / (z2 - z1);
				TY = (y2 - y1) / (z2 - z1);
				angle = sqrt(TX * TX + TY * TY);
				eachAngleEfficiency->Fill(hitsOnThePlate, angle);
				eachPlateEfficiency->Fill(hitsOnThePlate, iplate);
				eachTXEfficiency->Fill(hitsOnThePlate, TX);
				eachTYEfficiency->Fill(hitsOnThePlate, TY);
				trackID = t->ID();
				plate = iplate;
				if (iPID == 0)
				{
					x = 2 * x1 - x2;
					y = 2 * y1 - y2;
				}
				if (iPID >= 1 && iPID < nPID - 1)
				{
					x = (x1 + x2) / 2;
					y = (y1 + y2) / 2;
				}
				if (iPID == nPID - 1)
				{
					x = 2 * x2 - x1;
					y = 2 * y2 - y1;
				}
				effInfo->Fill();
			}
		}
	}
}

void FnuQualityCheck::PrintEfficiency(TString filename)
{
	// Plot efficiencies and print them.
	TCanvas *c = new TCanvas();
	c->Print(filename + "[");
	c->SetLogy(1);
	eachAngleEfficiency->GetCopyTotalHisto()->Draw();
	c->Print(filename);
	c->SetLogy(0);
	eachAngleEfficiency->Draw();
	c->Print(filename);
	eachPlateEfficiency->GetCopyTotalHisto()->Draw();
	c->Print(filename);
	eachPlateEfficiency->Draw();
	c->Print(filename);

	c->SetLogy(1);
	eachTXEfficiency->GetCopyTotalHisto()->Draw();
	c->Print(filename);
	c->SetLogy(0);
	eachTXEfficiency->Draw();
	c->Print(filename);

	c->SetLogy(1);
	eachTYEfficiency->GetCopyTotalHisto()->Draw();
	c->Print(filename);
	c->SetLogy(0);
	eachTYEfficiency->Draw();
	c->Print(filename);

	c->Print(filename + "]");
}

void FnuQualityCheck::WriteEfficiencyTree(TString filename)
{
	TFile fout(filename, "recreate");
	effInfo->Write();
	fout.Close();
}

void FnuQualityCheck::WriteEfficiency(TString filename)
{
	TFile fout(filename, "recreate");
	eachAngleEfficiency->Write();
	eachPlateEfficiency->Write();
	eachTXEfficiency->Write();
	eachTYEfficiency->Write();
	fout.Close();
}

void FnuQualityCheck::MakePositionHist()
{
	std::vector<double> positionXVec;
	std::vector<double> positionYVec;
	positionXVec.reserve(ntrk);
	positionYVec.reserve(ntrk);
	for (int itrk = 0; itrk < ntrk; itrk++)
	{
		EdbTrackP *t = pvr->GetTrack(itrk);
		if (t->N() < 5)
			continue;
		double X = t->X();
		double Y = t->Y();
		positionXVec.push_back(X);
		positionYVec.push_back(Y);
	}
	auto maxminXIterator = std::minmax_element(positionXVec.begin(), positionXVec.end());
	auto maxminYIterator = std::minmax_element(positionYVec.begin(), positionYVec.end());
	double minX = *maxminXIterator.first;
	double maxX = *maxminXIterator.second;
	double minY = *maxminYIterator.first;
	double maxY = *maxminYIterator.second;
	double marginX = (maxX - minX) / 10;
	double marginY = (maxY - minY) / 10;
	double minXAxis = minX - marginX;
	double maxXAxis = maxX + marginX;
	double minYAxis = minY - marginY;
	double maxYAxis = maxY + marginY;

	positionHist = new TH2D("positionHist", "position distribution (" + title + ");x (#mum);y (#mum);Ntracks / cm^{2}", 100, minXAxis, maxXAxis, 100, minYAxis, maxYAxis);
	double area = (maxXAxis - minXAxis) / 100 / 10000 * (maxYAxis - minYAxis) / 100 / 10000; // in cm^2
	for (int itrk = 0; itrk < positionXVec.size(); itrk++)
	{
		positionHist->Fill(positionXVec.at(itrk), positionYVec.at(itrk), 1 / area);
	}
}
void FnuQualityCheck::PrintPositionHist(TString filename)
{
	TCanvas ctemp;
	ctemp.SetRightMargin(0.15);
	positionHist->Draw("colz");
	positionHist->SetStats(0);
	positionHist->GetYaxis()->SetTitleOffset(1.5);
	ctemp.Print(filename);
}
void FnuQualityCheck::WritePositionHist(TString filename)
{
	TFile fout(filename, "recreate");
	positionHist->Write();
	fout.Close();
}
void FnuQualityCheck::MakeAngleHist()
{
	std::vector<double> angleXVec;
	std::vector<double> angleYVec;
	angleXVec.reserve(ntrk);
	angleYVec.reserve(ntrk);
	TH2D h2ForAngleCenterCalculation("h2ForAngleCenterCalculation", "angle;tan#theta_{x};tan#theta_{y}", 100, -0.02, 0.02, 100, -0.02, 0.02);
	// loop over the tracks
	for (int itrk = 0; itrk < ntrk; itrk++)
	{
		EdbTrackP *t = pvr->GetTrack(itrk);
		int nseg = t->N();
		// cut is nseg>=5
		if (nseg < 5)
			continue;
		// loop over the segments
		double X[nseg];
		double Y[nseg];
		double Z[nseg];
		for (int iseg = 0; iseg < nseg; iseg++)
		{
			X[iseg] = t->GetSegment(iseg)->X();
			Y[iseg] = t->GetSegment(iseg)->Y();
			Z[iseg] = t->GetSegment(iseg)->Z();
		}
		double TX, TY, a0;
		CalcLSM(Z, X, nseg, a0, TX);
		CalcLSM(Z, Y, nseg, a0, TY);
		h2ForAngleCenterCalculation.Fill(TX, TY);
		angleXVec.push_back(TX);
		angleYVec.push_back(TY);
	}
	int locmax, locmay, locmaz;
	h2ForAngleCenterCalculation.GetMaximumBin(locmax, locmay, locmaz);
	const double angleCenterX = h2ForAngleCenterCalculation.GetXaxis()->GetBinCenter(locmax);
	const double angleCenterY = h2ForAngleCenterCalculation.GetYaxis()->GetBinCenter(locmay);
	double halfRangeNarrow = 0.01;
	double minXAxisNarrow = angleCenterX - halfRangeNarrow;
	double maxXAxisNarrow = angleCenterX + halfRangeNarrow;
	double minYAxisNarrow = angleCenterY - halfRangeNarrow;
	double maxYAxisNarrow = angleCenterY + halfRangeNarrow;
	double halfRangeWide = 0.3;
	double minXAxisWide = angleCenterX - halfRangeWide;
	double maxXAxisWide = angleCenterX + halfRangeWide;
	double minYAxisWide = angleCenterY - halfRangeWide;
	double maxYAxisWide = angleCenterY + halfRangeWide;
	angleHistNarrow = new TH2D("angleHistNarrow", "angle distribution narrow (" + title + ");tan#theta_{x};tan#theta_{y};Ntracks", 200, minXAxisNarrow, maxXAxisNarrow, 200, minYAxisNarrow, maxYAxisNarrow);
	angleHistWide = new TH2D("angleHistWide", "angle distribution wide (" + title + ");tan#theta_{x};tan#theta_{y};Ntracks", 200, minXAxisWide, maxXAxisWide, 200, minYAxisWide, maxYAxisWide);
	for (int itrk = 0; itrk < angleXVec.size(); itrk++)
	{
		angleHistNarrow->Fill(angleXVec.at(itrk), angleYVec.at(itrk));
		angleHistWide->Fill(angleXVec.at(itrk), angleYVec.at(itrk));
	}
}
void FnuQualityCheck::PrintAngleHist(TString filename)
{
	TCanvas ctemp;
	ctemp.Print(filename + "[");
	ctemp.SetRightMargin(0.15);
	angleHistNarrow->Draw("colz");
	angleHistNarrow->SetStats(0);
	angleHistNarrow->GetYaxis()->SetTitleOffset(1.5);
	ctemp.Print(filename);
	ctemp.SetLogz();
	angleHistWide->Draw("colz");
	angleHistWide->SetStats(0);
	angleHistWide->GetYaxis()->SetTitleOffset(1.5);
	ctemp.Print(filename);
	ctemp.Print(filename + "]");
}
void FnuQualityCheck::WriteAngleHist(TString filename)
{
	TFile fout(filename, "recreate");
	angleHistNarrow->Write();
	angleHistWide->Write();
	fout.Close();
}

void FnuQualityCheck::MakeNsegHist()
{
	nsegHist = new TH1I("nsegHist", "nseg (" + title + ");nseg;Ntracks", nPID, 0.5, nPID + 0.5);
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
	nplHist = new TH1I("nplHist", "npl (" + title + ");npl;Ntracks", plMax - plMin + 1, 0.5, plMax - plMin + 1.5);
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
	firstPlateHist = new TH1I("firstPlateHist", "first plate (" + title + ");plate;Ntracks", plMax - plMin + 1, plMin - 0.5, plMax + 0.5);
	lastPlateHist = new TH1I("lastPlateHist", "last plate (" + title + ");plate;Ntracks", plMax - plMin + 1, plMin - 0.5, plMax + 0.5);
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
	secondDifferenceHist = new TH1D("secondDifferenceHist", Form("second difference (cell length %d);second differecne (#mum);entries", cellLength), 100, -30, 30);
	secondDifferenceTree->Draw("secondDiffX>>+secondDifferenceHist");
	secondDifferenceTree->Draw("secondDiffY>>+secondDifferenceHist");
}
void FnuQualityCheck::Summarize()
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

	// angle distribution wide
	c.cd(2);
	gStyle->SetPadRightMargin(0.25);
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

	// efficiency for each angle
	c.cd(4);
	TString eachAngleEfficiencyOriginalTitle = eachAngleEfficiency->GetTitle();
	eachAngleEfficiency->SetTitle("efficiency for each angle");
	gPad->SetRightMargin(0.0);
	eachAngleEfficiency->Draw();

	// efficiency for each plate
	c.cd(5);
	TString eachPlateEfficiencyOriginalTitle = eachPlateEfficiency->GetTitle();
	eachPlateEfficiency->SetTitle("efficiency for each plate");
	gPad->SetRightMargin(0.0);
	eachPlateEfficiency->Draw();

	// position resolution in x and y for each plate
	c.cd(6);
	c.SetGridx(1);
	gPad->SetRightMargin(0.0);
	gPad->SetTopMargin(0.13);
	TMultiGraph *mg2 = new TMultiGraph("mg2", "position resolution for each plate;plate;position resolution (#mum)");
	sigmaXGraph->SetMarkerSize(0.5);
	sigmaYGraph->SetMarkerSize(0.5);
	mg2->Add(sigmaXGraph);
	mg2->Add(sigmaYGraph);
	mg2->Draw("ap");
	TLegend *leg2 = new TLegend(0.3, 0.87, 1.0, 0.93);
	leg2->AddEntry(sigmaXGraph, "resolutionX", "p");
	leg2->AddEntry(sigmaYGraph, "resolutionY", "p");
	leg2->SetNColumns(2);
	leg2->Draw();

	// npl and nseg
	c.cd(7);
	gStyle->SetPadRightMargin(0.0);
	gStyle->SetPadTopMargin(0.13);
	gStyle->SetOptStat("e");
	gStyle->SetStatH(0.09);
	gStyle->SetStatW(0.2);
	gStyle->SetStatX(1);
	gStyle->SetStatY(1);
	auto nplHistMaxY = nplHist->GetMaximum();
	auto nsegHistMaxY = nsegHist->GetMaximum();
	auto maxY = std::max(nplHistMaxY, nsegHistMaxY);
	auto minX = nplHist->GetBinCenter(nplHist->FindFirstBinAbove());
	auto maxX = nplHist->GetBinCenter(nplHist->FindLastBinAbove());
	gPad->DrawFrame(minX - (maxX - minX) * 0.05, 0.7, maxX + (maxX - minX) * 0.05, maxY * 2, "nseg and npl;nseg or npl;# tracks");
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
	c.cd(8);
	gPad->SetTopMargin(0.13);
	gPad->SetRightMargin(0);
	// gStyle->SetOptStat("e");
	// gStyle->SetStatY(0.87);
	auto firstPlateHistMaxY = firstPlateHist->GetMaximum();
	auto lastPlateHistMaxY = lastPlateHist->GetMaximum();
	maxY = std::max(firstPlateHistMaxY, lastPlateHistMaxY);
	gPad->DrawFrame(plMin - (plMax - plMin) * 0.05, 0.7, plMax + (plMax - plMin) * 0.05, maxY * 2, "start and end plate;plate;# tracks");
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
	c.cd(9);
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
	secondDifferenceHist->Draw();
	gPad->UseCurrentStyle();

	c.cd(10);
	TText tx;
	tx.DrawTextNDC(0.1, 0.9, title);

	c.Print("summary_plot_"+title+".pdf");

	// set original title
	positionHist->SetTitle(positionHistOriginalTitle);
	angleHistNarrow->SetTitle(angleHistNarrowOriginalTitle);
	eachAngleEfficiency->SetTitle(eachAngleEfficiencyOriginalTitle);
	eachPlateEfficiency->SetTitle(eachPlateEfficiencyOriginalTitle);

	// set default style
	gROOT->SetStyle("Modern");
	gStyle->SetOptStat();
	gStyle->SetStatY(0.935);
	positionHist->UseCurrentStyle();
	angleHistNarrow->UseCurrentStyle();
	eachPlateEfficiency->UseCurrentStyle();
	eachAngleEfficiency->UseCurrentStyle();
	sigmaXGraph->UseCurrentStyle();
	sigmaYGraph->UseCurrentStyle();
	nsegHist->UseCurrentStyle();
	nplHist->UseCurrentStyle();
	firstPlateHist->UseCurrentStyle();
	lastPlateHist->UseCurrentStyle();
}