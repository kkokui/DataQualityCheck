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
#include <TCanvas.h>
#include <TSystem.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TEventList.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>
double sigmaX, sigmaY, meanX, meanY;
int entries, pl;

FnuQualityCheck::FnuQualityCheck(EdbPVRec *pvr, TString title)
{
	this->pvr = pvr;
	this->title = title;
	nPID = pvr->Npatterns();
	plMin = pvr->GetPattern(0)->Plate();
	plMax = pvr->GetPattern(nPID - 1)->Plate();
}

FnuQualityCheck::~FnuQualityCheck()
{
}
// void FnuQualityCheck::CalcDeltaXYFromRootFile(TString fname, double Xcenter, double Ycenter, TCut cut,double bin_width)
// {
// 	TFile fin(fname);
// 	TTree *tracks = (TTree*)fin.Get("tracks");
// 	tracks->Draw(">>elist",cut);
// 	TEventList *list = (TEventList*)gDirectory->Get("elist");
// 	nPID = tracks->GetMaximum("s.ePID")+1;
// 	// int plMin = deltaXY->GetMinimum("pl");
// 	// int plMax = deltaXY->GetMaximum("pl");
// 	deltaXY = new TTree("deltaXY", "deltaXY");
// 	double deltaX, deltaY, deltaTX, deltaTY, x_t, y_t, slopeX, slopeY;
// 	int plate, cross_the_line, trid, nseg;
// 	deltaXY->Branch("deltaX", &deltaX);
// 	deltaXY->Branch("deltaY", &deltaY);
// 	deltaXY->Branch("deltaTX", &deltaTX);
// 	deltaXY->Branch("deltaTY", &deltaTY);
// 	deltaXY->Branch("x", &x_t);
// 	deltaXY->Branch("y", &y_t);
// 	deltaXY->Branch("slopeX", &slopeX);
// 	deltaXY->Branch("slopeY", &slopeY);
// 	deltaXY->Branch("pl", &plate);
// 	deltaXY->Branch("cross_the_line", &cross_the_line);
// 	deltaXY->Branch("trid", &trid);
// 	deltaXY->Branch("nseg", &nseg);
// 	double tx3;
// 	double ty3;
// 	TBranch *bX = tracks->GetBranch("s.eX");
// 	TBranch *bY = tracks->GetBranch("s.eY");
// 	TBranch *bZ = tracks->GetBranch("s.eZ");
// 	TBranch *bTX = tracks->GetBranch("s.eTX");
// 	TBranch *bTY = tracks->GetBranch("s.eTY");
// 	TBranch *bnseg = tracks->GetBranch("nseg");
// 	for(int itrk = 0;itrk<list->GetN();itrk++)
// 	{
// 		int entr = list->GetEntry(itrk);

// 	}
// }
void FnuQualityCheck::CalcDeltaXY(int ntrk, double Xcenter, double Ycenter, double bin_width)
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

void FnuQualityCheck::PlotPosRes()
{
	// gSystem->Load("libTree");
	// int plMax = posResPar->GetMaximum("pl");
	// int plMin = posResPar->GetMinimum("pl");
	TCanvas *c1 = new TCanvas();
	c1->Print("pos_res/sigmaPar_" + title + ".pdf[");
	posResPar->Draw("meanY:plate");
	int N = posResPar->GetSelectedRows();
	TGraph *grY = new TGraph(N, posResPar->GetV2(), posResPar->GetV1());
	c1->SetGridx(1);
	grY->SetMarkerStyle(20);
	grY->SetMarkerColor(kBlue);
	grY->SetTitle("meanY");
	posResPar->Draw("meanX:plate");
	N = posResPar->GetSelectedRows();
	TGraph *grX = new TGraph(N, posResPar->GetV2(), posResPar->GetV1());
	grX->SetMarkerStyle(20);
	grX->SetMarkerColor(kRed);
	grX->SetTitle("meanX");
	TMultiGraph *mg = new TMultiGraph("mg", "mean:plate " + title + ";plate;mean (#mum)");
	mg->Add(grX);
	mg->Add(grY);
	// TH2D *fr = new TH2D("fr", "mean:plate (" + title + ");plate;mean (#mum)", 10, plMin, plMax, 10, -0.26, 0.26);
	// fr->Draw();
	// fr->SetStats(0);
	mg->Draw("ap");
	TLegend *leg = new TLegend(0.9, 0.8, 1.0, 0.9);
	leg->AddEntry(grX, "", "p");
	leg->AddEntry(grY, "", "p");
	leg->Draw();
	c1->Print("pos_res/sigmaPar_" + title + ".pdf");
	// c1->Print("pos_res/sigmaPar_" + title + ".C");
	// sigmaX
	c1->SetGridx(0);
	posResPar->Draw("sigmaX>>h1(100,0,1)", "abs(meanX)<100&&entries>0");
	TH1F *h1 = (TH1F *)gDirectory->Get("h1");
	// sigmaY
	posResPar->Draw("sigmaY>>h2(100,0,1)", "abs(meanY)<100&&entries>0");
	TH1F *h2 = (TH1F *)gDirectory->Get("h2");
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
	posResPar->Draw("sigmaX:plate");
	N = posResPar->GetSelectedRows();
	TGraph *grsigX = new TGraph(N, posResPar->GetV2(), posResPar->GetV1());
	grsigX->SetMarkerStyle(20);
	grsigX->SetMarkerColor(kRed);
	grsigX->SetTitle("sigmaX");
	posResPar->Draw("sigmaY:plate");
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
	mg2->Draw("ap");
	TLegend *leg2 = new TLegend(0.9, 0.8, 1.0, 0.9);
	leg2->AddEntry(grsigX, "", "p");
	leg2->AddEntry(grsigY, "", "p");
	leg2->Draw();
	c1->Print("pos_res/sigmaPar_" + title + ".pdf");
	c1->Clear();
	c1->Print("pos_res/sigmaPar_" + title + ".pdf");
	c1->Print("pos_res/sigmaPar_" + title + ".pdf]");
}
void FnuQualityCheck::PrintDeltaXYHist()
{
	TCanvas c;
	c.Print("pos_res/deltaxy_" + title + ".pdf[");
	for (int ient = 0; ient < htree->GetEntriesFast(); ient++)
	{
		htree->GetEntry(ient);
		hdeltaX->Draw();
		c.Print("pos_res/deltaxy_" + title + ".pdf");
		hdeltaY->Draw();
		c.Print("pos_res/deltaxy_" + title + ".pdf");
		printf("Histograms for plate %d have been printed\n", plate);
	}
	c.Print("pos_res/deltaxy_" + title + ".pdf]");
}
void FnuQualityCheck::WritePosResPar()
{
	TFile fout("pos_res/sigmaPar_" + title + ".root", "recreate");
	posResPar->Write();
	fout.Close();
}
void FnuQualityCheck::WriteDeltaXY()
{
	TFile fout("deltaXY/tree_" + title + ".root", "recreate");
	deltaXY->Write();
	fout.Close();
}

void FnuQualityCheck::CalcEfficiency()
{
	TEfficiency *pEff_angle = 0;
	TEfficiency *pEff_plate = 0;

	int ntrk = pvr->Ntracks();
	double bins[] = {0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5};
	int nbins = 26;
	double bins_pm[] = {-0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.19, -0.18, -0.17, -0.16, -0.15, -0.14, -0.13, -0.12, -0.11, -0.10, -0.09, -0.08, -0.07, -0.06, -0.05, -0.04, -0.03, -0.02, -0.01, 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5};
	int nbins_pm = 52;

	TH1D *h_angle_total = new TH1D("hist_angle_total", title + ";tan#theta;", nbins, bins);
	TH1D *h_angle_passed = new TH1D("hist_angle_passed", title + ";tan#theta;", nbins, bins);
	TH1D *h_plate_total = new TH1D("hist_plate_total", title + ";plate;", plMax - plMin, plMin, plMax);
	TH1D *h_plate_passed = new TH1D("hist_plate_passed", title + ";plate;", plMax - plMin, plMin, plMax);
	TH1D *h_TX_total = new TH1D("hist_TX_total", title + ";tan#theta;", nbins_pm, bins_pm);
	TH1D *h_TX_passed = new TH1D("hist_TX_passed", title + ";tan#theta;", nbins_pm, bins_pm);
	TH1D *h_TY_total = new TH1D("hist_TY_total", title + ";tan#theta;", nbins_pm, bins_pm);
	TH1D *h_TY_passed = new TH1D("hist_TY_passed", title + ";tan#theta;", nbins_pm, bins_pm);

	TFile treefout(Form("efficiency_output/tree_OtherInfo_%s.root", title.Data()), "recreate"); // to check the distribution of passed segments
	TTree *tree = new TTree("tree", "efficiency infomation");
	int trackID, pl, nseg, W, hitsOnThePlate;
	double x, y, angle, TX, TY;
	bool bpassed;
	tree->Branch("trackID", &trackID);
	tree->Branch("x", &x);
	tree->Branch("y", &y);
	tree->Branch("angle", &angle);
	tree->Branch("TX", &TX);
	tree->Branch("TY", &TY);
	tree->Branch("pl", &pl);
	tree->Branch("nseg", &nseg);
	tree->Branch("W", &W);
	tree->Branch("hitsOnThePlate", &hitsOnThePlate);

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

				if (s->PID() == iPID - 2 || s->PID() == iPID + 2)
					counts++;
				if (s->PID() == iPID - 1)
				{
					x1 = s->X();
					y1 = s->Y();
					z1 = s->Z();
					counts++;
				}
				if (s->PID() == iPID + 1)
				{
					x2 = s->X();
					y2 = s->Y();
					z2 = s->Z();
					counts++;
				}
				if (s->PID() == iPID)
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
				h_angle_total->Fill(angle);
				h_plate_total->Fill(iplate);
				h_TX_total->Fill(TX);
				h_TY_total->Fill(TY);
				if (hitsOnThePlate == 1)
				{
					h_angle_passed->Fill(angle);
					h_plate_passed->Fill(iplate);
					h_TX_passed->Fill(TX);
					h_TY_passed->Fill(TY);
				}
				trackID = t->ID();
				x = (x1 + x2) / 2;
				y = (y1 + y2) / 2;
				tree->Fill();
			}
		}
	}
	tree->Write();
	TCanvas *c = new TCanvas();
	c->Print(Form("efficiency_output/hist_efficiency_%s.pdf[", title.Data()));
	c->SetLogy(1);
	h_angle_total->Draw();
	c->Print(Form("efficiency_output/hist_efficiency_%s.pdf", title.Data()));
	c->SetLogy(0);
	pEff_angle = new TEfficiency(*h_angle_passed, *h_angle_total);
	pEff_angle->SetTitle(Form("Efficiency for each angle (%s);tan#theta;efficiency", title.Data()));
	pEff_angle->Draw();
	c->Print(Form("efficiency_output/hist_efficiency_%s.pdf", title.Data()));
	h_plate_total->Draw();
	c->Print(Form("efficiency_output/hist_efficiency_%s.pdf", title.Data()));
	pEff_plate = new TEfficiency(*h_plate_passed, *h_plate_total);
	pEff_plate->SetTitle(Form("Efficiency for each plate (%s);plate;efficiency", title.Data()));
	pEff_plate->Draw();
	c->Print(Form("efficiency_output/hist_efficiency_%s.pdf", title.Data()));

	c->SetLogy(1);
	h_TX_total->Draw();
	c->Print(Form("efficiency_output/hist_efficiency_%s.pdf", title.Data()));
	c->SetLogy(0);
	TGraphAsymmErrors *grEff_TX = new TGraphAsymmErrors(h_TX_passed, h_TX_total);
	grEff_TX->SetTitle(Form("Efficiency for each TX (%s);tan#theta;efficiency", title.Data()));
	grEff_TX->Draw("ap");
	c->Print(Form("efficiency_output/hist_efficiency_%s.pdf", title.Data()));

	c->SetLogy(1);
	h_TY_total->Draw();
	c->Print(Form("efficiency_output/hist_efficiency_%s.pdf", title.Data()));
	c->SetLogy(0);
	TGraphAsymmErrors *grEff_TY = new TGraphAsymmErrors(h_TY_passed, h_TY_total);
	grEff_TY->SetTitle(Form("Efficiency for each TY (%s);tan#theta;efficiency", title.Data()));
	grEff_TY->Draw("ap");
	c->Print(Form("efficiency_output/hist_efficiency_%s.pdf", title.Data()));

	c->Print(Form("efficiency_output/hist_efficiency_%s.pdf]", title.Data()));
	TFile fout(Form("efficiency_output/efficiency_%s.root", title.Data()), "recreate");
	pEff_angle->Write();
	pEff_plate->Write();
	grEff_TX->Write();
	grEff_TY->Write();
	fout.Close();

	FILE *ftxt = fopen("efficiency_output/efficiency.txt", "w"); // Make txt file for smearing parameters.
	for (int i = 1; i <= nbins; i++)
	{
		fprintf(ftxt, "%.3f, ", (bins[i] + bins[i - 1]) / 2);
	}
	fprintf(ftxt, "\n");
	for (int i = 1; i <= nbins; i++)
	{
		fprintf(ftxt, "%f ", pEff_angle->GetEfficiency(i));
	}
	fprintf(ftxt, "\n");
	fclose(ftxt);
}