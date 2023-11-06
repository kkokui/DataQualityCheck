int read_momentum_from_nt()
{
    gStyle->SetPadRightMargin(0.15);
    gStyle->SetPadLeftMargin(0.15);
    TChain *ntBefore = new TChain("nt");
    ntBefore->Add("momentum_output/nt_before_align_*");
    TCanvas *c = new TCanvas();
    ntBefore->Draw("Prec_Coord:Prec_Coord>>hist(200,0,7500,200,0,7500)", "16==nicell&&angle_diff_max<1", "colz");
    TH2F *hist = (TH2F *)gDirectory->Get("hist");
    hist->SetStats(0);
    hist->SetTitle("Prec_Coord before alignment : Prec_Coord before alignment;P_{before} (GeV);P_{before} (GeV);entries");
    c->SetLogz();
    c->SetRealAspectRatio();
    c->Print("momentum_output/momentum_before_after.pdf[");
    double binWidthArr[] = {5000,2000,1000,500};
    int nPlot = 0;
    for (int i=0;i<4;i++)
    {
        for (int iRobust = 10; iRobust >= 6; iRobust--)
        {
            TChain *ntAfter = new TChain("nt");
            ntAfter->Add(Form("momentum_output/nt_after_align_binWidth%.0f_robustFactor%.1f_*",binWidthArr[i],iRobust*0.1));

            ntAfter->AddFriend(ntBefore,"before");
            ntAfter->Draw("Prec_Coord:before.Prec_Coord>>hist(200,0,7500,200,0,7500)", "16==nicell&&before.angle_diff_max<1", "colz");
            TH2F *hist = (TH2F*)gDirectory->Get("hist");
            hist->SetStats(0);
            hist->SetTitle(Form("Prec_Coord after alignment : Prec_Coord before alignment (binWidth%.0f_robustFactor%.1f);P_{before} (GeV);P_{after} (GeV);entries", binWidthArr[i] ,iRobust * 0.1));
            // ntAfter->Draw("trid:before.trid","","colz"); // to check if the track ID is the same
            c->Print("momentum_output/momentum_before_after.pdf");
            nPlot++;
        }
    }
    c->Print("momentum_output/momentum_before_after.pdf]");
    return nPlot;
}
