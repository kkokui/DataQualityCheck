#include <stdio.h>
#include <iostream>
#include <TFile.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TTree.h>
#include <TSystem.h>
#include <TF1.h>
#include <TStyle.h>
struct TreeEntry
{
	double sigmaX, sigmaY, meanX, meanY;
	int entries, pl;
};
int main(int argc,char* argv[])
{
	if(argc<3)
	{
		printf("usage: %s deltaXY/tree_BeforeAlign.root title\n",argv[0]);
		return 0;
	}
	TString filename = argv[1];
	TString filename_short = argv[1];
	filename_short.ReplaceAll("deltaXY/tree_","");
	filename_short.ReplaceAll(".root","");
	TString title = argv[2];

	const float angcut = 0.01;
	
	if(0==TFile::Open(filename))
	{
		return 1;
	}
	TTree *tree = (TTree* )gDirectory->Get("tree");
	
	TFile::Open("pos_res/sigmaPar_"+filename_short+".root", "recreate");
	TTree *par = new TTree("par","parameter of position displacement");
	TreeEntry te;
	par->Branch("TreeEntry",&te,"sigmaX/D:sigmaY:meanX:meanY:entries/I:pl");
	// TNtupleD *par = new TNtupleD("par","par","sigmaX:sigmaY:entries:meanX:meanY:x:y:p2X:p2Y:pl");

	TF1 *f = new TF1("gaus","gaus",-2,2);
	f->SetParLimits(5,0,0.4);
	gStyle->SetOptFit();
	TCanvas *c1 = new TCanvas();
	c1->Print("pos_res/deltaxy_"+filename_short+".pdf[");
	
	TH1D *deltax = new TH1D("deltax","deltax",100,-2,2);
	TH1D *deltay = new TH1D("deltay","deltay",100,-2,2);
	
	int plMin = tree->GetMinimum("pl");
	int plMax = tree->GetMaximum("pl");
	
	for(int ipl=plMin;ipl<=plMax;ipl++){
		if(tree->GetEntries(Form("pl==%d",ipl))==0) continue;
		tree->Draw("deltaX>>htempX",Form("abs(tx+0.01)<%f&&abs(ty)<%f&&pl==%d",angcut,angcut,ipl));
		TH1F *htempX = (TH1F*)gDirectory->Get("htempX");
		te.meanX = htempX->GetMean();
		delete htempX;
		tree->Draw("deltaX>>deltax",Form("abs(deltaY)<=2&&abs(deltaX)<=2&&abs(tx+0.01)<%f&&abs(ty)<%f&&pl==%d",angcut,angcut,ipl) );
		deltax->SetTitle(Form("pl%d %s;deltaX (#mum);",ipl,title.Data()));
		
		f->SetParameters(1000,0,0.2);
		
		deltax->Fit(f,"Q","",-0.5,0.5);
		te.sigmaX = f->GetParameter(2);

		c1->Print("pos_res/deltaxy_" + filename_short + ".pdf");
		tree->Draw("deltaY>>htempY",Form("abs(tx+0.01)<%f&&abs(ty)<%f&&pl==%d",angcut,angcut,ipl));
		TH1F *htempY = (TH1F*)gDirectory->Get("htempY");
		te.meanY = htempY->GetMean();
		delete htempY;
		tree->Draw("deltaY>>deltay",Form("abs(deltaY)<=2&&abs(deltaX)<=2&&abs(tx+0.01)<%f&&abs(ty)<%f&&pl==%d",angcut,angcut,ipl));
		deltay->SetTitle(Form("pl%d %s;deltaY (#mum);",ipl,title.Data()));
		f->SetParameters(1000,0,0.2);
		deltay->Fit(f,"Q","",-0.5,0.5);
		te.sigmaY = f->GetParameter(2);
		c1->Print("pos_res/deltaxy_" + filename_short + ".pdf");
		te.entries = deltay->GetEntries();
		te.pl = ipl;
		par->Fill();
		printf("Histograms for plate %d have been printed\n", ipl);
	}
	c1->Print("pos_res/deltaxy_" + filename_short + ".pdf]");
	par->Write();
	delete par;
	delete tree;
}