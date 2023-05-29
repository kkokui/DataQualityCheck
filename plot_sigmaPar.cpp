#include <stdio.h>
#include <iostream>
#include <TFile.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TTree.h>
#include <TSystem.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TMultiGraph.h>

int make_plot(TString filename, TString title)
{
	gSystem->Load("libTree");
	if (0 == TFile::Open(filename))
	{
		return 0;
	}
	TTree *par = (TTree *)gDirectory->Get("par");
	int plMax = par->GetMaximum("pl");
	int plMin = par->GetMinimum("pl");

	TString filename_short = filename;
	filename_short.ReplaceAll("pos_res/sigmaPar_", "");
	filename_short.ReplaceAll(".root", "");

	TCanvas *c1 = new TCanvas();
	c1->Print(Form("pos_res/sigmaPar_" + filename_short + ".pdf["));

	// meanY and meanX : plate
	par->Draw("meanY:pl");
	int N = par->GetSelectedRows();
	TGraph *grY = new TGraph(N, par->GetV2(), par->GetV1());
	grY->SetMarkerStyle(20);
	grY->SetMarkerColor(kBlue);
	grY->SetTitle("meanY");
	par->Draw("meanX:pl");
	N = par->GetSelectedRows();
	TGraph *grX = new TGraph(N, par->GetV2(), par->GetV1());
	grX->SetMarkerStyle(20);
	grX->SetMarkerColor(kRed);
	grX->SetTitle("meanX");
	TMultiGraph *mg = new TMultiGraph("mg", "mean:plate " + title + ";plate;mean (#mum)");
	mg->Add(grX);
	mg->Add(grY);
	TH2D *fr = new TH2D("fr", "mean:plate (" + title + ");plate;mean (#mum)", 10, plMin, plMax, 10, -0.26, 0.26);
	fr->Draw();
	fr->SetStats(0);
	mg->Draw("p");
	c1->SetGridx(1);
	TLegend *leg = new TLegend(0.7, 0.9, 0.9, 1.0);
	leg->AddEntry(grX, "", "p");
	leg->AddEntry(grY, "", "p");
	leg->Draw();
	// htemp->SetTitle("meanY:meanX;meanX (#mum);meanY (#mum)");
	c1->Print(Form("pos_res/sigmaPar_" + filename_short + ".pdf"));
	// sigmaX
	c1->SetGridx(0);
	par->Draw("sigmaX>>h1(100,0,1)", "abs(meanX)<100&&entries>0");
	TH1F *h1 = (TH1F *)gDirectory->Get("h1");
	// c1->Print(Form("pos_res/sigmaPar_" + filename_short + ".pdf"));
	// sigmaY
	par->Draw("sigmaY>>h2(100,0,1)", "abs(meanY)<100&&entries>0");
	TH1F *h2 = (TH1F *)gDirectory->Get("h2");
	// c1->Print(Form("pos_res/sigmaPar_" + filename_short + ".pdf"));
	// merge sigmaX and sigmaY
	TList *l = new TList;
	l->Add(h1);
	l->Add(h2);
	TH1F *h = new TH1F("Sigma", "position resolution (" + title + ");position resolution (#mum);number of areas", 100, 0, 1);
	h->Merge(l);
	h->Draw();
	c1->Print(Form("pos_res/sigmaPar_" + filename_short + ".pdf"));
	// sigmaX and sigmaY:pl
	c1->SetGridx(1);
	par->Draw("sigmaX:pl");
	N = par->GetSelectedRows();
	TGraph *grsigX = new TGraph(N, par->GetV2(), par->GetV1());
	grsigX->SetMarkerStyle(20);
	grsigX->SetMarkerColor(kRed);
	grsigX->SetTitle("sigmaX");
	par->Draw("sigmaY:pl");
	N = par->GetSelectedRows();
	TGraph *grsigY = new TGraph(N, par->GetV2(), par->GetV1());
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
	TLegend *leg2 = new TLegend(0.7, 0.9, 0.9, 1.0);
	leg2->AddEntry(grsigX, "", "p");
	leg2->AddEntry(grsigY, "", "p");
	leg2->Draw();
	c1->Print(Form("pos_res/sigmaPar_" + filename_short + ".pdf"));
	c1->Clear();
	c1->Print(Form("pos_res/sigmaPar_" + filename_short + ".pdf"));
	c1->Print(Form("pos_res/sigmaPar_" + filename_short + ".pdf]"));
}

int main(int argc, char *argv[])
{
	if (argc < 3)
	{
		printf("usage: %s pos_res/sigmaPar_BeforeAlign.root title\n", argv[0]);
		return 0;
	}
	TString filename = argv[1];
	TString title = argv[2];
	make_plot(filename, title);
}