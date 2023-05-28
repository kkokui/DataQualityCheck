#include <stdio.h>
#include <iostream>
#include <TFile.h>
#include <TCanvas.h>
#include <TArrow.h>
#include <TH2.h>
#include <TLatex.h>
#include <TTree.h>
#include <TSystem.h>

int main(int argc,char* argv[])
{
    if(argc<6)
    {
        printf("usage: %s deltaXY/tree_NAlign.root Xcenter Ycenter plMin plMax\n",argv[0]);
        return 0;
    }
    TString filename_deltaXY = argv[1];
    TString filename_short = argv[1];
    filename_short.ReplaceAll(".root","");
    filename_short.ReplaceAll("deltaXY/tree_","");
    double Xcenter, Ycenter;
    int plMin, plMax;
    sscanf(argv[2],"%lf",&Xcenter);
    sscanf(argv[3],"%lf",&Ycenter);
    sscanf(argv[4],"%d",&plMin);
    sscanf(argv[5],"%d",&plMax);

    const float angcut = 0.01;
    gSystem->Load("libTree");
    TFile::Open(filename_deltaXY);

    TTree *tree = (TTree *)gDirectory->Get("tree");
    
    TCanvas *c1 = new TCanvas();
    
	const double XYrange = 8500;
    
    c1->Print("deltaXY_XYdis/x_deltaX_"+filename_short+".pdf[");
    c1->Print("deltaXY_XYdis/y_deltaY_"+filename_short+".pdf[");
    c1->Print("deltaXY_XYdis/x_deltaY_"+filename_short+".pdf[");
    c1->Print("deltaXY_XYdis/y_deltaX_"+filename_short+".pdf[");
    
	TH2D *x_deltaXY = new TH2D("x_deltaXY","title",50,-5,5,100,Xcenter-XYrange,Xcenter+XYrange);
    TH2D *y_deltaXY = new TH2D("y_deltaXY","title",50,-5,5,100,Ycenter-XYrange,Ycenter+XYrange);
    //x_deltaXY->SetStats(0);
    x_deltaXY->GetYaxis()->SetTitleOffset(1.5);
    y_deltaXY->GetYaxis()->SetTitleOffset(1.5);
    //y_deltaXY->SetStats(0);
    for (int ipl = plMin; ipl <= plMax; ipl++)
    {
        if(tree->GetEntries(Form("pl==%d",ipl)))
        {
            //x:deltax
            tree->Draw("x:deltaX>>x_deltaXY",Form("pl==%d&&abs(tx + 0.01) < %f && abs(ty-0.004) < %f",ipl,angcut,angcut),"colz");
            x_deltaXY->SetTitle(Form("x:deltaX pl%d ;#deltaX (#mum);x (#mum)",ipl));
            c1->Print("deltaXY_XYdis/x_deltaX_"+filename_short+".pdf");
            //y:deltay
            tree->Draw("y:deltaY>>y_deltaXY",Form("pl==%d&&abs(tx + 0.01) < %f && abs(ty-0.004) < %f",ipl,angcut,angcut),"colz");
            y_deltaXY->SetTitle(Form("y:deltaY pl%d ;#deltaY (#mum);y (#mum)",ipl));
            c1->Print("deltaXY_XYdis/y_deltaY_"+filename_short+".pdf");
            //x:deltaY
            tree->Draw("x:deltaY>>x_deltaXY",Form("pl==%d&&abs(tx + 0.01) < %f && abs(ty-0.004) < %f",ipl,angcut,angcut),"colz");
            x_deltaXY->SetTitle(Form("x:deltaY pl%d ;#deltaY (#mum);x (#mum)",ipl));
            c1->Print("deltaXY_XYdis/x_deltaY_"+filename_short+".pdf");
            //y:deltaX
            tree->Draw("y:deltaX>>y_deltaXY",Form("pl==%d&&abs(tx + 0.01) < %f && abs(ty-0.004) < %f",ipl,angcut,angcut),"colz");
            y_deltaXY->SetTitle(Form("y:deltaX pl%d ;#deltaX (#mum);y (#mum)",ipl));
            c1->Print("deltaXY_XYdis/y_deltaX_"+filename_short+".pdf");
        }
    }
    c1->Print("deltaXY_XYdis/x_deltaX_"+filename_short+".pdf]");
    c1->Print("deltaXY_XYdis/y_deltaY_"+filename_short+".pdf]");
    c1->Print("deltaXY_XYdis/x_deltaY_"+filename_short+".pdf]");
    c1->Print("deltaXY_XYdis/y_deltaX_"+filename_short+".pdf]");

    c1->Print("deltaXY_XYdis/Arrow_deltaXY_"+filename_short+".pdf[");
    double x,y,deltaX,deltaY;
    int pl;
    tree->SetBranchAddress("x",&x);
    tree->SetBranchAddress("y",&y);
    tree->SetBranchAddress("deltaX",&deltaX);
    tree->SetBranchAddress("deltaY",&deltaY);
    tree->SetBranchAddress("pl",&pl);
    
    TArrow *arr = new TArrow();
    TH2F *frame = new TH2F("frame","title",10,Xcenter-XYrange-1000,Xcenter+XYrange+1000,10,Ycenter-XYrange-1000,Ycenter+XYrange+1000);
    frame->SetStats(0);
    frame->GetYaxis()->SetTitleOffset(1.5);
    TLatex *l = new TLatex;
    l->SetTextAlign(33);
    l->SetTextSize(0.03);
    l->SetTextFont(42);
    double bin_width=2000;
    int scale = 2000;
    for(int ipl = plMin;ipl<=plMax;ipl++)
    {
        if(tree->GetEntries(Form("pl==%d",ipl)))
        {
            frame->Draw();
            frame->SetTitle(Form("deltaXY pl%d ;x(#mum);y(#mum)",ipl));
            arr->DrawArrow(Xcenter+XYrange+1000,Ycenter+XYrange+2000,Xcenter+XYrange+1000-scale*0.5,Ycenter+XYrange+2000,0.008,">"); //equivalent to 0.5 μm.
            l->DrawLatex(Xcenter+XYrange+1000,Ycenter+XYrange+1500,"0.5 #mum");
            for (double iY = Ycenter-XYrange+bin_width/2; iY <= Ycenter+XYrange; iY += bin_width)
            {
                for (double iX = Xcenter-XYrange+bin_width/2; iX <= Xcenter+XYrange; iX += bin_width)
                {
                    float deltaXMean=0;
                    float deltaYMean=0;
                    int i=0;
                    for(int iEntry = 0;iEntry<tree->GetEntries();iEntry++)
                    {
                        tree->GetEntry(iEntry);
                        if(pl==ipl&&fabs(x-iX)<bin_width/2&&fabs(y-iY)<bin_width/2)
                        {
                            deltaXMean+=deltaX;
                            deltaYMean+=deltaY;
                            i++;
                        }
                    }
                    if(i==0) continue;
                    deltaXMean/=i;
                    deltaYMean/=i;
                    
                    arr->DrawArrow(iX,iY,iX+scale*deltaXMean,iY+scale*deltaYMean,0.008,">");
                }
            }
            c1->Print("deltaXY_XYdis/Arrow_deltaXY_"+filename_short+".pdf");
        }
    }
    c1->Print("deltaXY_XYdis/Arrow_deltaXY_"+filename_short+".pdf]");
    
}