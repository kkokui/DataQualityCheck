#include <TStyle.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TText.h>
#include <TFile.h>
// #include <TGraph.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TLegend.h>

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
    double meanResolution[20];
    double meanResolutionError[20];
    int iCondition = 0;
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
            meanResolution[iCondition] = resolution->GetMean();
            meanResolutionError[iCondition] = resolution->GetMeanError();
            iCondition++;
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

    gStyle->SetStatW(0.6);
    gStyle->SetStatH(0.6);
    gStyle->SetStatX(0.98);
    gStyle->SetStatY(0.91);
    gStyle->SetTitleX(0.5);
    gStyle->SetTitleY(1.0);
    gStyle->SetTitleSize(0.08, "");  // ヒストグラムタイトルの文字サイズ(""を書かないとX軸タイトルのサイズが変わる。)
    gStyle->SetTitleSize(0.06);      // X軸タイトルの文字サイズ("X"って書いてもいい。)
    gStyle->SetTitleSize(0.06, "Y"); // Y軸タイトルの文字サイズ
    gStyle->SetLabelSize(0.04);      // X軸目盛の文字サイズ
    gStyle->SetLabelSize(0.04, "Y"); // Y軸目盛の文字サイズ
    gStyle->SetPadLeftMargin(0.07);
    gStyle->SetPadRightMargin(0.02);
    gStyle->SetPadBottomMargin(0.12); // 下側の余白サイズ
    gStyle->SetPadTopMargin(0.09);

    TCanvas *c2 = new TCanvas("c2", "multiple resolutions with different binWidth", 1500, 1000);
    c2->Divide(3, 2);
    isubpad = 1;

    c2->cd(isubpad);

    TFile *fin = new TFile("/data/Users/kokui/FASERnu/F222/zone4/temp/TFD/vert32063_pl053_167_new/reco32_065000_050000/v15/pos_res/graph_hist_before_align.root");
    TH1D *sigmaXHist = (TH1D *)gDirectory->Get("sigmaXHist");
    TH1D *sigmaYHist = (TH1D *)gDirectory->Get("sigmaYHist");
    TList *l = new TList;
    l->Add(sigmaXHist);
    l->Add(sigmaYHist);
    TH1F *resolution = new TH1F("resolution", "before align;position resolution (#mum)", 100, 0, 1);
    resolution->Merge(l);
    float xmin = 0;
    float ymin = 0;
    float xmax = 1;
    float ymax = 82;
    gPad->DrawFrame(xmin, ymin, xmax, ymax, "before align;position resolution (#mum)");
    resolution->Draw("sames");
    // for (int irobust = 10; irobust >= 6; irobust--)
    int iRobust = 7;
    for (int ibinwidth = 0; ibinwidth < 4; ibinwidth++)
    {
        isubpad++;
        c2->cd(isubpad);
        // float robustFactor = irobust * 0.1;
        // fin = new TFile(Form("deltaXY_XYdis/mean_deltaXY_after_align_binWidth500_robustFactor%.1f.root", robustFactor));
        fin = new TFile(Form("/data/Users/kokui/FASERnu/F222/zone4/temp/TFD/vert32063_pl053_167_new/reco32_065000_050000/v15/pos_res/graph_hist_after_align_binWidth%.0f_robustFactor%.1f.root", binWidthArr[ibinwidth], iRobust * 0.1));
        TH1D *sigmaXHist = (TH1D *)gDirectory->Get("sigmaXHist");
        TH1D *sigmaYHist = (TH1D *)gDirectory->Get("sigmaYHist");
        TList *l = new TList;
        l->Add(sigmaXHist);
        l->Add(sigmaYHist);
        TH1F *resolution = new TH1F("resolution", Form("%.1f mm division;position resolution (#mum)", binWidthArr[ibinwidth] / 1000), 100, 0, 1);
        resolution->Merge(l);
        gPad->DrawFrame(xmin, ymin, xmax, ymax, Form("%.1f mm division;position resolution (#mum)", binWidthArr[ibinwidth] / 1000));
        resolution->Draw("sames");
    }
    // c2->Print("multiple_arrow_plots_different_robust_factor.pdf");
    c2->Print("pos_res/multiple_resolutions_different_bin_width.pdf");

    gStyle->SetPadLeftMargin(0.12);
    gStyle->SetTitleOffset(1, "XY");
    TCanvas *c3 = new TCanvas();
    TH1F *frResolutionVsRobust = gPad->DrawFrame(0.55, 0.205, 1.05, 0.3, ";robust factor;position resolution (#mum)");
    frResolutionVsRobust->GetXaxis()->SetNdivisions(505);
    frResolutionVsRobust->GetYaxis()->SetNdivisions(505);
    TGraphErrors *graphResolutionVsRobust[4];
    // TLegend *legendResolutionVsRobust = new TLegend(0.3, 0.21);
    TLegend *legendResolutionVsRobust = new TLegend(0.16, 0.7,0.4,0.9);
    for (int iBinWidth = 0; iBinWidth < 4; iBinWidth++)
    {
        graphResolutionVsRobust[iBinWidth] = new TGraphErrors();
        graphResolutionVsRobust[iBinWidth]->SetName(Form("%.1f mm division", binWidthArr[iBinWidth] / 1000));
        int iPoint = 0;
        for (int iRobust = 10; iRobust >= 6; iRobust--)
        {
            graphResolutionVsRobust[iBinWidth]->SetPoint(iPoint, iRobust * 0.1, meanResolution[iBinWidth * 5 + iPoint]);
            graphResolutionVsRobust[iBinWidth]->SetPointError(iPoint, 0, meanResolutionError[iBinWidth * 5 + iPoint]);
                iPoint++;
        }
        graphResolutionVsRobust[iBinWidth]->SetLineColor(iBinWidth + 1);
        graphResolutionVsRobust[iBinWidth]->Draw("l");
        legendResolutionVsRobust->AddEntry(graphResolutionVsRobust[iBinWidth], Form("%.1f mm division", binWidthArr[iBinWidth] / 1000), "l");
    }
    legendResolutionVsRobust->Draw();
    c3->Print("pos_res/position_resolution_vs_robust.pdf");

    return isubpad;
}