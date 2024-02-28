#include "FnuResolution.h"

#include <stdio.h>
#include <numeric>

#include <TStyle.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TArrow.h>
#include <TLatex.h>

FnuResolution::FnuResolution(EdbPVRec *pvr, TString title)
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

FnuResolution::~FnuResolution()
{
}

TTree *FnuResolution::CalcDeltaXY(double Xcenter, double Ycenter, double bin_width)
{
	TTree *deltaXY = new TTree("deltaXY", "deltaXY");
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
    return deltaXY;
}
TTree *FnuResolution::MakeHistDeltaXY(TTree *deltaXY)
{

	hdeltaX = new TH1D("hdeltaX", "hdeltaX", 100, -2, 2);
	hdeltaY = new TH1D("hdeltaY", "hdeltaY", 100, -2, 2);
	TTree *treeHistDeltaXY = new TTree("treeHistDeltaXY", "treeHistDeltaXY");
	treeHistDeltaXY->Branch("hdeltaX", &hdeltaX);
	treeHistDeltaXY->Branch("hdeltaY", &hdeltaY);
	treeHistDeltaXY->Branch("plate", &plate);

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
		
		treeHistDeltaXY->Fill();
		hdeltaX->Reset();
		hdeltaY->Reset();
	}
	return treeHistDeltaXY;

}

TTree *FnuResolution::FitHistDeltaXY(TTree *treeHistDeltaXY)
{
	// Fit histograms of delta x and y.
	// Before using this method, do CalcDeltaXY() to make deltaXY data.
	TTree *positionResolutionPar = new TTree("positionResolutionPar","parameters for position resolution");
    
	positionResolutionPar->Branch("sigmaX", &sigmaX);
	positionResolutionPar->Branch("sigmaY", &sigmaY);
	positionResolutionPar->Branch("meanX", &meanX);
	positionResolutionPar->Branch("meanY", &meanY);
	positionResolutionPar->Branch("entries", &entries);
	positionResolutionPar->Branch("plate", &plate);

	treeHistDeltaXY->SetBranchAddress("hdeltaX",&hdeltaX);
	treeHistDeltaXY->SetBranchAddress("hdeltaY",&hdeltaY);
	
	TF1 *f = new TF1("gaus", "gaus", -2, 2);
	f->SetParLimits(5, 0, 0.4);
	gStyle->SetOptFit();

	for (int ient = 0; ient < treeHistDeltaXY->GetEntriesFast(); ient++)
	{
		treeHistDeltaXY->GetEntry(ient);
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
		positionResolutionPar->Fill();
	}
    return positionResolutionPar;
}

void FnuResolution::CalcLSM(double x[], double y[], int N, double &a0, double &a1)
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
void FnuResolution::MakeGraphHistPositionResolution(TTree *positionResolutionPar)
{
	// Make graphs and histograms about position resolution.
	// This makes graphs of mean:plate, histograms of position resolution and graphs of position resolution : plate.

	// meanX and meanY : plate
	positionResolutionPar->Draw("meanX:plate");
	int N = positionResolutionPar->GetSelectedRows();
	meanXGraph = new TGraph(N, positionResolutionPar->GetV2(), positionResolutionPar->GetV1());
	meanXGraph->SetMarkerStyle(20);
	meanXGraph->SetMarkerColor(kRed);
	meanXGraph->SetNameTitle("meanYGraph", "mean Y (" + title + ");plate;mean (#mum)");
	positionResolutionPar->Draw("meanY:plate");
	N = positionResolutionPar->GetSelectedRows();
	meanYGraph = new TGraph(N, positionResolutionPar->GetV2(), positionResolutionPar->GetV1());
	meanYGraph->SetMarkerStyle(20);
	meanYGraph->SetMarkerColor(kBlue);
	meanYGraph->SetNameTitle("meanXGraph", "mean X (" + title + ");plate;mean (#mum)");

	// sigmaX and sigmaY
	sigmaXHist = new TH1D("sigmaXHist", "position resolution X (" + title + ");position resolution (#mum)", 100, 0, 1);
	positionResolutionPar->Draw("sigmaX>>sigmaXHist", "abs(meanX)<100&&entries>0");
	sigmaYHist = new TH1D("sigmaYHist", "position resolution Y (" + title + ");position resolution (#mum)", 100, 0, 1);
	positionResolutionPar->Draw("sigmaY>>sigmaYHist", "abs(meanY)<100&&entries>0");

	// sigmaX and sigmaY:pl
	positionResolutionPar->Draw("sigmaX:plate");
	N = positionResolutionPar->GetSelectedRows();
	sigmaXGraph = new TGraph(N, positionResolutionPar->GetV2(), positionResolutionPar->GetV1());
	sigmaXGraph->SetMarkerStyle(20);
	sigmaXGraph->SetMarkerColor(kRed);
	sigmaXGraph->SetNameTitle("sigmaXGraph", "position resolution X (" + title + ");plate;position resolution (#mum)");
	positionResolutionPar->Draw("sigmaY:plate");
	N = positionResolutionPar->GetSelectedRows();
	sigmaYGraph = new TGraph(N, positionResolutionPar->GetV2(), positionResolutionPar->GetV1());
	sigmaYGraph->SetMarkerStyle(20);
	sigmaYGraph->SetMarkerColor(kBlue);
	sigmaYGraph->SetNameTitle("sigmaYGraph", "position resolution Y (" + title + ");plate;position resolution (#mum)");
}

void FnuResolution::WriteGraphHistPositionResolution(TString filename)
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
void FnuResolution::PrintGraphHistPositionResolution(TString filename)
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
void FnuResolution::PrintHistDeltaXY(TTree *treeHistDeltaXY, TString filename)
{
	TCanvas c;
	// c.Print("pos_res/deltaxy_" + title + ".pdf[");
	c.Print(filename + "[");
	for (int ient = 0; ient < treeHistDeltaXY->GetEntriesFast(); ient++)
	{
		treeHistDeltaXY->GetEntry(ient);
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
void FnuResolution::CalcMeanDeltaXY(TTree *deltaXY, double Xcenter, double Ycenter, double cellSize)
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

void FnuResolution::WriteMeanDeltaXY(TString filename)
{
	TFile fout(filename, "recreate");
	meanDeltaXY->Write();
	fout.Close();
}

int FnuResolution::PrintMeanDeltaXYArrowPlot(TString filename)
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

TGraph *FnuResolution::GetSigmaXGraph() const
{
    return sigmaXGraph;
}

TGraph *FnuResolution::GetSigmaYGraph() const
{
    return sigmaYGraph;
}

TTree *FnuResolution::MakeHistDeltaTXY(TTree *deltaXY)
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
void FnuResolution::PrintHistDeltaTXY(TTree *treeHistDeltaTXY,TString filename)
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
TTree *FnuResolution::FitHistDeltaTXY(TTree *treeHistDeltaTXY)
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
void FnuResolution::MakeGraphHistAngleResolution(TTree *angleResolutionPar)
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
void FnuResolution::PrintGraphHistAngleResolution(TString filename)
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