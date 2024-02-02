#include <TStyle.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TText.h>
#include <TFile.h>

int print_all_position_resolution()
{
    gStyle->SetStatW(0.6);
    gStyle->SetStatH(0.6);
    gStyle->SetOptStat("mr");
    gStyle->SetLabelSize(0.06, "xy");
    gStyle->SetLabelOffset(0.01, "XYZ");
    gStyle->SetTitleSize(0.1, "xy");
    gStyle->SetTitleSize(0.08, "");
    gStyle->SetPadBottomMargin(0.22);
    gStyle->SetPadRightMargin(0.02);
    gStyle->SetPadLeftMargin(0.07);
    gStyle->SetPadTopMargin(0.02);
    gStyle->SetStatX(0.98);
    gStyle->SetStatY(0.98);
    TCanvas *c1 = new TCanvas();
    c1->Divide(6, 5);
    int isubpad = 0;
    for (int iRobust = 10; iRobust >= 6; iRobust--)
    {
        isubpad++;
        c1->cd(isubpad);
        TText *t = new TText(0.5, 0, Form("robust factor %.1f", iRobust * 0.1));
        t->SetNDC(1);
        t->SetTextAlign(21);
        t->SetTextSize(0.15);
        t->Draw();
    }
    isubpad++;
    double binWidthArr[] = {5000, 2000, 1000, 500};
    for (int i = 0; i < 4; i++)
    {
        for (int iRobust = 10; iRobust >= 6; iRobust--)
        {
            isubpad++;
            c1->cd(isubpad);
            TFile *fin = new TFile(Form("/data/Users/kokui/FASERnu/F222/zone4/temp/TFD/vert32063_pl053_167_new/reco32_065000_050000/v15/pos_res/graph_hist_after_align_binWidth%.0f_robustFactor%.1f.root", binWidthArr[i], iRobust * 0.1));
            TH1D *sigmaXHist = (TH1D *)gDirectory->Get("sigmaXHist");
            TH1D *sigmaYHist = (TH1D *)gDirectory->Get("sigmaYHist");
            TList *l = new TList;
            l->Add(sigmaXHist);
            l->Add(sigmaYHist);
            // TH1F *resolution = new TH1F("resolution", Form("binWidth%.0f_robustFactor%.1f;position resolution (#mum)", binWidthArr[i], iRobust * 0.1), 100, 0, 1);
            TH1F *resolution = new TH1F("resolution", ";position resolution (#mum)", 100, 0, 1);
            resolution->Merge(l);
            resolution->Draw();
        }
        isubpad++;
        c1->cd(isubpad);
        TText *t = new TText(0.07, 0.07, Form("%.1f mm division", binWidthArr[i] / 1000));
        t->SetNDC(1);
        // t->SetTextAlign(12);
        t->SetTextSize(0.13);
        t->SetTextAngle(90);
        t->Draw();
    }
    c1->Print("pos_res/all_position_resolution.pdf");
    return isubpad;
}