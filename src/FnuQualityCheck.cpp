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

FnuQualityCheck::FnuQualityCheck(EdbPVRec *pvr, TString title)
	: pvr(pvr),
	  title(title),
	  nPID(pvr->Npatterns()),
	  plMin(pvr->GetPattern(0)->Plate()),
	  plMax(pvr->GetPattern(nPID - 1)->Plate()),
	  ntrk(pvr->Ntracks())
{
	double bins_arr_angle[] = {0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5};
	SetBinsAngle(26,bins_arr_angle);
	double bins_arr_TXTY[] = {-0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.19, -0.18, -0.17, -0.16, -0.15, -0.14, -0.13, -0.12, -0.11, -0.10, -0.09, -0.08, -0.07, -0.06, -0.05, -0.04, -0.03, -0.02, -0.01, 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5};
	SetBinsTXTY(52,bins_arr_TXTY);
}

void FnuQualityCheck::SetBinsAngle(int nbins, double bins[])
{
	bins_vec_angle.assign(&bins[0], &bins[nbins+1]);
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
				if(sPID>iPID+2)
					break;
				if(sPID<iPID-2)
					continue;
				x[sPID-iPID+2] = s->X();
				y[sPID-iPID+2] = s->Y();
				z[sPID-iPID+2] = s->Z();
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

void FnuQualityCheck::PlotPosRes(TString filename)
{
	// int plMax = posResPar->GetMaximum("pl");
	// int plMin = posResPar->GetMinimum("pl");
	TCanvas *c1 = new TCanvas();
	// c1->Print("pos_res/sigmaPar_" + title + ".pdf[");
	c1->Print(filename + "[");
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
	c1->Print(filename);
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
	c1->Print(filename);
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
	std::copy(bins_vec_angle.begin(),bins_vec_angle.end(),bins_angle);
	int nbins_angle = bins_vec_angle.size()-1;
	double bins_TXTY[bins_vec_TXTY.size()];
	std::copy(bins_vec_TXTY.begin(), bins_vec_TXTY.end(), bins_TXTY);
	int nbins_TXTY = bins_vec_TXTY.size()-1;

	pEff_angle = new TEfficiency("Eff_angle", Form("Efficiency for each angle (%s);tan#theta", title.Data()), nbins_angle, bins_angle);
	pEff_plate = new TEfficiency("Eff_plate", Form("Efficiency for each plate (%s);plate", title.Data()), plMax - plMin, plMin, plMax);
	pEff_TX = new TEfficiency("Eff_TX", Form("Efficiency for each TX (%s);tan#theta", title.Data()), nbins_TXTY, bins_TXTY);
	pEff_TY = new TEfficiency("Eff_TY", Form("Efficiency for each TY (%s);tan#theta", title.Data()), nbins_TXTY, bins_TXTY);

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
				if (sPID > iPID+2)
					break;
				if (sPID < iPID-2)
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
				pEff_angle->Fill(hitsOnThePlate, angle);
				pEff_plate->Fill(hitsOnThePlate, iplate);
				pEff_TX->Fill(hitsOnThePlate, TX);
				pEff_TY->Fill(hitsOnThePlate, TY);
				trackID = t->ID();
				x = (x1 + x2) / 2;
				y = (y1 + y2) / 2;
				effInfo->Fill();
			}
		}
	}
}

void FnuQualityCheck::PlotEfficiency(TString filename)
{
	// Plot efficiencies and print them.
	TCanvas *c = new TCanvas();
	c->Print(filename+"[");
	c->SetLogy(1);
	pEff_angle->GetCopyTotalHisto()->Draw();
	c->Print(filename);
	c->SetLogy(0);
	pEff_angle->Draw();
	c->Print(filename);
	pEff_plate->GetCopyTotalHisto()->Draw();
	c->Print(filename);
	pEff_plate->Draw();
	c->Print(filename);

	c->SetLogy(1);
	pEff_TX->GetCopyTotalHisto()->Draw();
	c->Print(filename);
	c->SetLogy(0);
	pEff_TX->Draw();
	c->Print(filename);

	c->SetLogy(1);
	pEff_TY->GetCopyTotalHisto()->Draw();
	c->Print(filename);
	c->SetLogy(0);
	pEff_TY->Draw();
	c->Print(filename);

	c->Print(filename+"]");
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
	pEff_angle->Write();
	pEff_plate->Write();
	pEff_TX->Write();
	pEff_TY->Write();
	fout.Close();
}