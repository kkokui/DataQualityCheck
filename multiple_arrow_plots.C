#include <TFile.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TDirectory.h>
#include <TH2.h>
#include <TArrow.h>
#include <TLatex.h>
#include <TStyle.h>

int arrow_plot(TTree *meanDeltaXY, int ient, const char *title)
{
    int nent = meanDeltaXY->GetEntriesFast();
    if (0 == nent)
    {
        printf("nent = 0\n");
        return nent;
    }
    TH2D *deltaX2DHist = new TH2D();
    TH2D *deltaY2DHist = new TH2D();
    TH2D *entries2DHist = new TH2D();

    meanDeltaXY->SetBranchAddress("deltaX2DHist", &deltaX2DHist);
    meanDeltaXY->SetBranchAddress("deltaY2DHist", &deltaY2DHist);
    meanDeltaXY->SetBranchAddress("entries2DHist", &entries2DHist);
    meanDeltaXY->GetEntry(ient);
    // double minXaxis = deltaX2DHist->GetXaxis()->GetXmin();
    // double minYaxis = deltaX2DHist->GetYaxis()->GetXmin();
    // double maxXaxis = deltaX2DHist->GetXaxis()->GetXmax();
    // double maxYaxis = deltaX2DHist->GetYaxis()->GetXmax();
    double minXaxis = 60000;
    double minYaxis = deltaX2DHist->GetYaxis()->GetXmin();
    double maxXaxis = deltaX2DHist->GetXaxis()->GetXmax();
    double maxYaxis = 56000;
    double rangeXaxis = maxXaxis - minXaxis;
    double rangeYaxis = maxYaxis - minYaxis;
    int nbinsX = deltaX2DHist->GetXaxis()->GetNbins();
    int nbinsY = deltaX2DHist->GetYaxis()->GetNbins();
    TLatex *l = new TLatex;
    l->SetTextAlign(33);
    l->SetTextSize(0.05);
    l->SetTextFont(42);
    int scale = 2000; // Should not be hard coded...

    TH1F *fr = gPad->DrawFrame(minXaxis, minYaxis, maxXaxis, maxYaxis, title);
    TArrow *arrow = new TArrow();
    // arrow->SetAngle(80);
    arrow->SetLineWidth(1);
    arrow->DrawArrow(maxXaxis + rangeXaxis * 0.00, maxYaxis + rangeYaxis * 0.1, maxXaxis + rangeXaxis * 0.00 - scale * 0.5, maxYaxis + rangeYaxis * 0.1, 0.003, ">"); // equivalent to 0.5 μm.
    l->DrawLatex(maxXaxis + rangeXaxis * 0.00, maxYaxis + rangeYaxis * 0.07, "0.5 #mum");
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
            arrow->DrawArrow(binCenterX, binCenterY, binCenterX + scale * deltaXMean, binCenterY + scale * deltaYMean, 0.003, ">");
        }
    }
    return nent;
}

void multiple_arrow_plots()
{
    gStyle->SetTitleX(0.5);
    gStyle->SetTitleY(1.0);
    gStyle->SetTitleSize(0.08, "");  // ヒストグラムタイトルの文字サイズ(""を書かないとX軸タイトルのサイズが変わる。)
    gStyle->SetTitleSize(0.06);      // X軸タイトルの文字サイズ("X"って書いてもいい。)
    gStyle->SetTitleSize(0.06, "Y"); // Y軸タイトルの文字サイズ
    gStyle->SetLabelSize(0.04);      // X軸目盛の文字サイズ
    gStyle->SetLabelSize(0.04, "Y"); // Y軸目盛の文字サイズ
    gStyle->SetPadLeftMargin(0.2);
    gStyle->SetPadRightMargin(0);
    gStyle->SetPadBottomMargin(0.12); // 下側の余白サイズ
    gStyle->SetPadTopMargin(0.09);

    TCanvas *c0 = new TCanvas("c0", "multiple arrow plot", 1500, 1000);
    c0->Divide(3, 2);
    int isubpad = 1;

    c0->cd(isubpad);

    TFile *fin = new TFile("deltaXY_XYdis/mean_deltaXY_before_align.root");
    TTree *meanDeltaXY = (TTree *)gDirectory->Get("meanDeltaXY");
    arrow_plot(meanDeltaXY, 2, "before align;x (#mum);y (#mum)");
    delete fin;
    // for (int irobust = 10; irobust >= 6; irobust--)
    float binWidth[] = {5000,2000,1000,500};
    for (int ibinwidth = 0; ibinwidth <4; ibinwidth++)
    {
        isubpad++;
        c0->cd(isubpad);
        // float robustFactor = irobust * 0.1;
        // fin = new TFile(Form("deltaXY_XYdis/mean_deltaXY_after_align_binWidth500_robustFactor%.1f.root", robustFactor));
        fin = new TFile(Form("deltaXY_XYdis/mean_deltaXY_after_align_binWidth%.0f_robustFactor0.7.root", binWidth[ibinwidth]));
        meanDeltaXY = (TTree *)gDirectory->Get("meanDeltaXY");
        // arrow_plot(meanDeltaXY, 2, Form("robust factor %.1f;x (#mum);y (#mum)", robustFactor));
        arrow_plot(meanDeltaXY, 2, Form("%.1f mm division;x (#mum);y (#mum)", binWidth[ibinwidth] / 1000));
        delete fin;
    }
    // c0->Print("multiple_arrow_plots_different_robust_factor.pdf");
    c0->Print("multiple_arrow_plots_different_bin_width.pdf");
}