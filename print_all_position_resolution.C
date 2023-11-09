int print_all_position_resolution()
{
    gStyle->SetStatW(0.5);
    gStyle->SetStatH(0.4);
    gStyle->SetOptStat("mr");
    gStyle->SetLabelSize(0.05, "xy");
    gStyle->SetTitleSize(0.06, "xy");
    gStyle->SetTitleSize(0.08, "");
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadRightMargin(0);
    gStyle->SetPadLeftMargin(0.07);
    gStyle->SetStatX(1);
    gStyle->SetStatY(0.9);
    TCanvas *c1 = new TCanvas();
    c1->Divide(5, 4);
    int iDividedCanvas = 0;
    double binWidthArr[] = {5000, 2000, 1000, 500};
    for (int i = 0; i < 4; i++)
    {
        for (int iRobust = 10; iRobust >= 6; iRobust--)
        {
            iDividedCanvas++;
            c1->cd(iDividedCanvas);
            TFile *fin = new TFile(Form("/data/Users/kokui/FASERnu/F222/zone4/temp/TFD/vert32063_pl053_167_new/reco32_065000_050000/v15/pos_res/graph_hist_after_align_binWidth%.0f_robustFactor%.1f.root", binWidthArr[i], iRobust * 0.1));
            TH1D *sigmaXHist = (TH1D *)gDirectory->Get("sigmaXHist");
            TH1D *sigmaYHist = (TH1D *)gDirectory->Get("sigmaYHist");
            TList *l = new TList;
            l->Add(sigmaXHist);
            l->Add(sigmaYHist);
            TH1F *resolution = new TH1F("resolution", Form("binWidth%.0f_robustFactor%.1f;position resolution (#mum)", binWidthArr[i], iRobust * 0.1), 100, 0, 1);
            resolution->Merge(l);
            resolution->Draw();
        }
    }
    c1->Print("pos_res/all_position_resolution.pdf");
    return iDividedCanvas;
}