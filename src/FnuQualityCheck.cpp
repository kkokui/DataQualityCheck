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

FnuQualityCheck::FnuQualityCheck(EdbPVRec *pvr, TString title)
	: pvr(pvr),
	  title(title),
	  nPID(pvr->Npatterns()),
	  plMin(pvr->GetPattern(0)->Plate()),
	  plMax(pvr->GetPattern(nPID - 1)->Plate()),
	  ntrk(pvr->Ntracks()),
	  deltaTXV(), deltaXV(), deltaYV(), deltaTYV(), xV(), yV(), slopeXV(), slopeYV(),
	  crossTheLineV(), tridV(), nsegV()
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
	double tx3;
	double ty3;
	// Loop over the tracks
	for (int iPID = 2; iPID < nPID - 2; iPID++)
	{
		for (int itrk = 0; itrk < ntrk; itrk++)
		{
			EdbTrackP *t = pvr->GetTrack(itrk);
			int count = 0;
			double x[5];
			double y[5];
			double z[5];
			nseg = t->N();
			for (int iseg = 0; iseg < nseg; iseg++)
			{
				EdbSegP *s = t->GetSegment(iseg);
				int sPID = s->PID();
				if (sPID > iPID + 2)
					break;
				if (sPID < iPID - 2)
					continue;
				x[sPID - iPID + 2] = s->X();
				y[sPID - iPID + 2] = s->Y();
				z[sPID - iPID + 2] = s->Z();
				count++;
				if (sPID == iPID)
				{
					tx3 = s->TX();
					ty3 = s->TY();
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
				areaX[ipoint] = (x[ipoint] - (Xcenter - XYrange)) / bin_width;
				areaY[ipoint] = (y[ipoint] - (Ycenter - XYrange)) / bin_width;
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
			double x_updown[4], y_updown[4], z_updown[4];
			for (int i = 0; i < 2; i++)
			{
				x_updown[i] = x[i];
				y_updown[i] = y[i];
				z_updown[i] = z[i];
			}
			for (int i = 2; i < 4; i++)
			{
				x_updown[i] = x[i + 1];
				y_updown[i] = y[i + 1];
				z_updown[i] = z[i + 1];
			}
			double a0, slopeX, slopeY;
			CalcLSM(z_updown, x_updown, 4, a0, slopeX);
			double x3fit = a0 + slopeX * z[2];
			CalcLSM(z_updown, y_updown, 4, a0, slopeY);
			double y3fit = a0 + slopeY * z[2];

			// Calculate delta X and delta Y.
			deltaXV->push_back(x[2] - x3fit);
			deltaYV->push_back(y[2] - y3fit);
			deltaTXV->push_back(tx3 - slopeX);
			deltaTYV->push_back(ty3 - slopeY);
			xV->push_back(t->X());
			yV->push_back(t->Y());
			slopeXV->push_back(slopeX);
			slopeYV->push_back(slopeY);
			crossTheLineV->push_back(cross_the_line);
			tridV->push_back(t->ID());
			nsegV->push_back(nseg);
		}
		plate = pvr->GetPattern(iPID)->Plate();
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
	double angcut = 0.01;
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
			if (fabs(deltaYV->at(i)) <= 2 && fabs(deltaXV->at(i)) <= 2 && fabs(slopeXV->at(i) - slopeXMean) < angcut && fabs(slopeYV->at(i) - slopeYMean) < angcut)
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
		for (int iPID = 2; iPID < nPID - 2; iPID++)
		{
			int iplate = pvr->GetPattern(iPID)->Plate();
			int counts = 0;
			hitsOnThePlate = 0;
			W = 0;
			for (int iseg = 0; iseg < nseg; iseg++)
			{
				EdbSegP *s = t->GetSegment(iseg);
				int sPID = s->PID();
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
				x = (x1 + x2) / 2;
				y = (y1 + y2) / 2;
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
		angleXVec.push_back(TX);
		angleYVec.push_back(TY);
	}
	const auto angleXMean = std::accumulate(angleXVec.begin(), angleXVec.end(), 0.0) / angleXVec.size();
	const auto angleYMean = std::accumulate(angleYVec.begin(), angleYVec.end(), 0.0) / angleYVec.size();
	double halfRangeNarrow = 0.008;
	double minXAxisNarrow = angleXMean - halfRangeNarrow;
	double maxXAxisNarrow = angleXMean + halfRangeNarrow;
	double minYAxisNarrow = angleYMean - halfRangeNarrow;
	double maxYAxisNarrow = angleYMean + halfRangeNarrow;
	double halfRangeWide = 0.3;
	double minXAxisWide = angleXMean - halfRangeWide;
	double maxXAxisWide = angleXMean + halfRangeWide;
	double minYAxisWide = angleYMean - halfRangeWide;
	double maxYAxisWide = angleYMean + halfRangeWide;
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
void FnuQualityCheck::PrintSummaryPlot()
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

	// angle distribution
	c.cd(2);
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

	// efficiency for each angle
	c.cd(3);
	TString eachAngleEfficiencyOriginalTitle = eachAngleEfficiency->GetTitle();
	eachAngleEfficiency->SetTitle("efficiency for each angle");
	gPad->SetRightMargin(0.0);
	eachAngleEfficiency->Draw();

	// efficiency for each plate
	c.cd(4);
	TString eachPlateEfficiencyOriginalTitle = eachPlateEfficiency->GetTitle();
	eachPlateEfficiency->SetTitle("efficiency for each plate");
	gPad->SetRightMargin(0.0);
	eachPlateEfficiency->Draw();

	// position resolution in x and y for each plate
	c.cd(5);
	c.SetGridx(1);
	gPad->SetRightMargin(0.0);
	gPad->SetTopMargin(0.13);
	TMultiGraph *mg2 = new TMultiGraph("mg2", "position resolution for each plate");
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
	c.cd(6);
	gStyle->SetPadRightMargin(0.0);
	gStyle->SetPadTopMargin(0.13);
	gStyle->SetOptStat("e");
	gStyle->SetStatH(0.05);
	gStyle->SetStatW(0.2);
	gStyle->SetStatX(1);
	gStyle->SetStatY(0.87);
	TString nplHistOriginalTitle = nplHist->GetTitle();
	nplHist->SetTitle("nseg and npl");
	nplHist->Draw();
	nsegHist->Draw("sames");
	gPad->UseCurrentStyle();
	gPad->Update();
	nsegHist->SetLineColor(kRed);
	TLegend *legNsegNpl = new TLegend(0.3, 0.87, 1.0, 0.93);
	legNsegNpl->AddEntry(nsegHist, "nseg", "l");
	legNsegNpl->AddEntry(nplHist, "npl", "l");
	legNsegNpl->SetNColumns(2);
	legNsegNpl->Draw();
	gPad->SetLogy();

	c.cd(7);
	TString firstPlateHistOriginalTitle = firstPlateHist->GetTitle();
	firstPlateHist->SetTitle("start and end plate");
	gPad->SetTopMargin(0.13);
	gPad->SetRightMargin(0);
	gStyle->SetOptStat("e");
	gStyle->SetStatY(0.87);
	firstPlateHist->Draw();
	lastPlateHist->Draw("same");
	TLegend *legFirstLast = new TLegend(0.3, 0.87, 1.0, 0.93);
	legFirstLast->AddEntry(firstPlateHist, "first plate", "l");
	legFirstLast->AddEntry(lastPlateHist, "last plate", "l");
	legFirstLast->SetNColumns(2);
	legFirstLast->Draw();
	gPad->UseCurrentStyle();
	lastPlateHist->SetLineColor(kRed);

	c.Print("summary_plot_test.pdf");

	// set original title
	positionHist->SetTitle(positionHistOriginalTitle);
	angleHistNarrow->SetTitle(angleHistNarrowOriginalTitle);
	eachAngleEfficiency->SetTitle(eachAngleEfficiencyOriginalTitle);
	eachPlateEfficiency->SetTitle(eachPlateEfficiencyOriginalTitle);
	nplHist->SetTitle(nplHistOriginalTitle);
	firstPlateHist->SetTitle(firstPlateHistOriginalTitle);

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