#include "FnuQualityCheck.h"

#include <stdio.h>

#include <EdbDataSet.h>
#include <TGraph.h>
#include <TH2.h>
#include <TObjArray.h>
#include <EdbPattern.h>
#include <TMath.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TLegend.h>
#include <TMultiGraph.h>

FnuQualityCheck::FnuQualityCheck(TString ititle)
{
	title = ititle;
	// file = new TFile("file.root","recreate");
}

FnuQualityCheck::~FnuQualityCheck()
{
}

void FnuQualityCheck::CalcDeltaXY(EdbPVRec *pvr, int ntrk, double Xcenter, double Ycenter, double bin_width)
{
	nPID = pvr->Npatterns();
	plMin = pvr->GetPattern(0)->Plate();
	plMax = pvr->GetPattern(nPID - 1)->Plate();
	deltaXY = new TTree("deltaXY", "deltaXY");
	// gDirectory->ls();
	double deltaX, deltaY, deltaTX, deltaTY, x_t, y_t, slopeX, slopeY;
	int plate, cross_the_line, trid, nseg;
	deltaXY->Branch("deltaX", &deltaX);
	deltaXY->Branch("deltaY", &deltaY);
	deltaXY->Branch("deltaTX", &deltaTX);
	deltaXY->Branch("deltaTY", &deltaTY);
	deltaXY->Branch("x", &x_t);
	deltaXY->Branch("y", &y_t);
	deltaXY->Branch("slopeX", &slopeX);
	deltaXY->Branch("slopeY", &slopeY);
	deltaXY->Branch("pl", &plate);
	deltaXY->Branch("cross_the_line", &cross_the_line);
	deltaXY->Branch("trid", &trid);
	deltaXY->Branch("nseg", &nseg);
	double tx3;
	double ty3;
	// Loop over the tracks
	for (int itrk = 0; itrk < ntrk; itrk++)
	{
		EdbTrackP *t = pvr->GetTrack(itrk);
		if (abs(t->TX() + 0.01) >= 0.01 || abs(t->TY() - 0.004) >= 0.01 || t->N() < 5)
			continue;

		trid = t->ID();
		nseg = t->N();
		TGraph grX;
		TGraph grY;

		for (int iPID = 2; iPID < nPID - 2; iPID++)
		{
			int count = 0;
			double x[5];
			double y[5];
			double z[5];
			for (int iseg = 0; iseg < nseg; iseg++)
			{
				EdbSegP *s = t->GetSegment(iseg);
				for (int ipoint = 0; ipoint < 5; ipoint++)
				{
					if (s->PID() == iPID - 2 + ipoint)
					{
						x[ipoint] = s->X();
						y[ipoint] = s->Y();
						z[ipoint] = s->Z();
						count++;
					}
				}
				if (s->PID() == iPID)
				{
					tx3 = s->TX();
					ty3 = s->TY();
				}
			}
			if (count != 5)
				continue;
			// printf("%d\n",grX.GetN());
			int areaX[5];
			int areaY[5];
			for (int ipoint = 0; ipoint < 5; ipoint++)
			{
				areaX[ipoint] = (x[ipoint] - (Xcenter - XYrange)) / bin_width;
				areaY[ipoint] = (y[ipoint] - (Ycenter - XYrange)) / bin_width;
			}
			cross_the_line = 0;
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
			double a0;
			lsm(z_updown, x_updown, 4, a0, slopeX);
			double x3fit = a0 + slopeX * z[2];
			lsm(z_updown, y_updown, 4, a0, slopeY);
			double y3fit = a0 + slopeY * z[2];

			// Calculate delta X and delta Y.
			deltaX = x[2] - x3fit;
			deltaY = y[2] - y3fit;
			deltaTX = tx3 - slopeX;
			deltaTY = ty3 - slopeY;
			x_t = t->X();
			y_t = t->Y();
			plate = pvr->GetPattern(iPID)->Plate();
			deltaXY->Fill();
		}
	}
}

void FnuQualityCheck::FitDeltaXY()
{
	// deltaXY->Show(1);
	double angcut = 0.01;
	// file->Delete("deltaXY");
	posResPar = new TTree("posResPar", "posResPar");
	// file->Add(posResPar);
	// gDirectory->ls();
	double sigmaX, sigmaY, meanX, meanY;
	int entries, pl;
	posResPar->Branch("sigmaX", &sigmaX);
	posResPar->Branch("sigmaY", &sigmaY);
	posResPar->Branch("meanX", &meanX);
	posResPar->Branch("meanY", &meanY);
	posResPar->Branch("entries", &entries);
	posResPar->Branch("pl", &pl);
	TF1 *f = new TF1("gaus", "gaus", -2, 2);
	f->SetParLimits(5, 0, 0.4);
	gStyle->SetOptFit();
	TCanvas *c1 = new TCanvas();
	c1->Print("pos_res/deltaxy_" + title + ".pdf[");

	TH1D *deltax = new TH1D("deltax", "deltax", 100, -2, 2);
	TH1D *deltay = new TH1D("deltay", "deltay", 100, -2, 2);

	// int plMin = deltaXY->GetMinimum("pl");
	// int plMax = deltaXY->GetMaximum("pl");

	for (int ipl = plMin; ipl <= plMax; ipl++)
	{
		if (deltaXY->GetEntries(Form("pl==%d", ipl)) == 0)
			continue;
		deltaXY->Draw("deltaX>>htempX", Form("abs(slopeX+0.01)<%f&&abs(slopeY)<%f&&pl==%d", angcut, angcut, ipl));
		TH1D *htempX = (TH1D *)gDirectory->Get("htempX");
		meanX = htempX->GetMean();
		delete htempX;
		deltaXY->Draw("deltaX>>deltax", Form("abs(deltaY)<=2&&abs(deltaX)<=2&&abs(slopeX+0.01)<%f&&abs(slopeY)<%f&&pl==%d", angcut, angcut, ipl));
		deltax->SetTitle(Form("pl%d %s;deltaX (#mum);", ipl, title.Data()));

		f->SetParameters(1000, 0, 0.2);

		deltax->Fit(f, "Q", "", -0.5, 0.5);
		sigmaX = f->GetParameter(2);

		c1->Print("pos_res/deltaxy_" + title + ".pdf");
		deltaXY->Draw("deltaY>>htempY", Form("abs(slopeX+0.01)<%f&&abs(slopeY)<%f&&pl==%d", angcut, angcut, ipl));
		TH1D *htempY = (TH1D *)gDirectory->Get("htempY");
		meanY = htempY->GetMean();
		entries = htempY->GetEntries();
		delete htempY;
		deltaXY->Draw("deltaY>>deltay", Form("abs(deltaY)<=2&&abs(deltaX)<=2&&abs(slopeX+0.01)<%f&&abs(slopeY)<%f&&pl==%d", angcut, angcut, ipl));
		deltay->SetTitle(Form("pl%d %s;deltaY (#mum);", ipl, title.Data()));
		f->SetParameters(1000, 0, 0.2);
		deltay->Fit(f, "Q", "", -0.5, 0.5);
		sigmaY = f->GetParameter(2);
		c1->Print("pos_res/deltaxy_" + title + ".pdf");
		pl = ipl;
		posResPar->Fill();
		printf("Histograms for plate %d have been printed\n", ipl);
	}
	c1->Print("pos_res/deltaxy_" + title + ".pdf]");
	// posResPar->Write();

	// posResPar->Draw("meanY:pl");
	// c1->Print("pos_res/test.pdf");
	// posResPar->Show(1);
	delete deltax;
	delete deltay;
}

void FnuQualityCheck::lsm(double x[], double y[], int N, double &a0, double &a1)
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

void FnuQualityCheck::PlotPosRes()
{
	// gSystem->Load("libTree");
	// int plMax = posResPar->GetMaximum("pl");
	// int plMin = posResPar->GetMinimum("pl");
	TCanvas *c1 = new TCanvas("test");
	c1->Print("pos_res/sigmaPar_" + title + ".pdf[");
	// gDirectory->ls();
	// posResPar->BranchRef();
	// TFile *fout1 = new TFile("pos_res/sigmaPar_" + title + ".root");
	// gDirectory->ls();
	// posResPar = (TTree *)gDirectory->Get("posResPar");

	// delete deltaXY;
	// meanY and meanX : plate
	// posResPar->Show(1);
	// posResPar->Print();
	posResPar->Draw("meanY:pl");
	// posResPar->GetV2();
	double *testv = posResPar->GetV1();

	// std::cout << testx[0] << std::endl;
	// std::cout << posResPar->GetV2()[0] << std::endl;
	TH1D *htest = new TH1D("name", "test", 100, 0, 100);
	// htest->Fill(50);
	// htest->Draw();
	// c1->Print("pos_res/sigmaPar_" + title + ".pdf");
	int N = posResPar->GetSelectedRows();
	TGraph *grY = new TGraph(N, posResPar->GetV2(), posResPar->GetV1());
	c1->SetGridx(1);
	// printf("debug\n");
	grY->SetMarkerStyle(20);
	grY->SetMarkerColor(kBlue);
	grY->SetTitle("meanY");
	posResPar->Draw("meanX:pl");
	N = posResPar->GetSelectedRows();
	TGraph *grX = new TGraph(N, posResPar->GetV2(), posResPar->GetV1());
	grX->SetMarkerStyle(20);
	grX->SetMarkerColor(kRed);
	grX->SetTitle("meanX");
	TMultiGraph *mg = new TMultiGraph("mg", "mean:plate " + title + ";plate;mean (#mum)");
	mg->Add(grX);
	mg->Add(grY);
	TH2D *fr = new TH2D("fr", "mean:plate (" + title + ");plate;mean (#mum)", 10, plMin, plMax, 10, -0.26, 0.26);
	// fr->Draw();
	fr->SetStats(0);
	mg->Draw("p");
	TLegend *leg = new TLegend(0.9, 0.8, 1.0, 0.9);
	leg->AddEntry(grX, "", "p");
	leg->AddEntry(grY, "", "p");
	leg->Draw();
	// htemp->SetTitle("meanY:meanX;meanX (#mum);meanY (#mum)");
	c1->Print("pos_res/sigmaPar_" + title + ".pdf");
	// sigmaX
	c1->SetGridx(0);
	posResPar->Draw("sigmaX>>h1(100,0,1)", "abs(meanX)<100&&entries>0");
	TH1F *h1 = (TH1F *)gDirectory->Get("h1");
	// c1->Print(Form("pos_res/sigmaPar_" + title + ".pdf"));
	// sigmaY
	posResPar->Draw("sigmaY>>h2(100,0,1)", "abs(meanY)<100&&entries>0");
	TH1F *h2 = (TH1F *)gDirectory->Get("h2");
	// c1->Print(Form("pos_res/sigmaPar_" + title + ".pdf"));
	// merge sigmaX and sigmaY
	TList *l = new TList;
	l->Add(h1);
	l->Add(h2);
	TH1F *h = new TH1F("Sigma", "position resolution (" + title + ");position resolution (#mum);number of areas", 100, 0, 1);
	h->Merge(l);
	h->Draw();
	c1->Print("pos_res/sigmaPar_" + title + ".pdf");
	// sigmaX and sigmaY:pl
	c1->SetGridx(1);
	posResPar->Draw("sigmaX:pl");
	N = posResPar->GetSelectedRows();
	TGraph *grsigX = new TGraph(N, posResPar->GetV2(), posResPar->GetV1());
	grsigX->SetMarkerStyle(20);
	grsigX->SetMarkerColor(kRed);
	grsigX->SetTitle("sigmaX");
	posResPar->Draw("sigmaY:pl");
	N = posResPar->GetSelectedRows();
	TGraph *grsigY = new TGraph(N, posResPar->GetV2(), posResPar->GetV1());
	grsigY->SetMarkerStyle(20);
	grsigY->SetMarkerColor(kBlue);
	grsigY->SetTitle("sigmaY");
	TMultiGraph *mg2 = new TMultiGraph("mg2", "sigma:plate (" + title + ");plate;sigma (#mum)");
	mg2->Add(grsigX);
	mg2->Add(grsigY);
	// TH2D *fr = new TH2D("fr",";plate;",10,0,52,10,-1.5,1.5);
	// fr->Draw();
	// fr->SetStats(0);
	mg2->Draw("p");
	TLegend *leg2 = new TLegend(0.9, 0.8, 1.0, 0.9);
	leg2->AddEntry(grsigX, "", "p");
	leg2->AddEntry(grsigY, "", "p");
	leg2->Draw();
	c1->Print("pos_res/sigmaPar_" + title + ".pdf");
	c1->Clear();
	c1->Print("pos_res/sigmaPar_" + title + ".pdf");
	c1->Print("pos_res/sigmaPar_" + title + ".pdf]");
}

void FnuQualityCheck::WritePosResPar()
{
	// TFile *fout1 = new TFile("pos_res/sigmaPar_" + title + ".root", "recreate");
	// posResPar->Write();
	// fout1->Close();
	// posResPar->Show(1);
}
