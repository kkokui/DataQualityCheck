#include <TStyle.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TH2F.h>

int read_momentum_from_nt()
{
    gStyle->SetPadRightMargin(0.15);
    gStyle->SetPadLeftMargin(0.15);
    TChain *ntBefore = new TChain("nt");
    ntBefore->Add("momentum_output/nt_before_align_icellmax32*");
    TCanvas *c = new TCanvas();
    ntBefore->Draw("Prec_Coord:Prec_Coord>>hist(200,0,7500,200,0,7500)", "32==nicell&&angle_diff_max<1", "colz");
    TH2F *hist = (TH2F *)gDirectory->Get("hist");
    hist->SetStats(0);
    hist->SetTitle("Prec_Coord before alignment : Prec_Coord before alignment;P_{before} (GeV);P_{before} (GeV);entries");
    c->SetLogz();
    c->SetRealAspectRatio();
    c->Print("momentum_output/momentum_before_after_icellmax32.pdf[");
    double binWidthArr[] = {5000, 2000, 1000, 500};
    int nPlot = 0;
    for (int i = 0; i < 4; i++)
    {
        for (int iRobust = 10; iRobust >= 6; iRobust--)
        {
            TChain *ntAfter = new TChain("nt");
            int nFilesAdded = ntAfter->Add(Form("momentum_output/nt_after_align_binWidth%.0f_robustFactor%.1f_icellmax32*", binWidthArr[i], iRobust * 0.1));
            if (0 == nFilesAdded)
            {
                printf("momentum_output/nt_after_align_binWidth%.0f_robustFactor%.1f_icellmax32* does not exist", binWidthArr[i], iRobust * 0.1);
                continue;
            }
            ntAfter->AddFriend(ntBefore, "before");
            ntAfter->Draw("Prec_Coord:before.Prec_Coord>>hist(200,0,7500,200,0,7500)", "32==nicell&&before.angle_diff_max<1", "colz");
            TH2F *hist = (TH2F *)gDirectory->Get("hist");
            hist->SetStats(0);
            hist->SetTitle(Form("momentum before and after alignment (binWidth%.0f_robustFactor%.1f);P_{before} (GeV);P_{after} (GeV);entries", binWidthArr[i], iRobust * 0.1));
            // ntAfter->Draw("trid:before.trid","","colz"); // to check if the track ID is the same
            c->Print("momentum_output/momentum_before_after_icellmax32.pdf");
            int n7000Before = ntAfter->GetEntries("32==nicell&&before.angle_diff_max<1&&before.Prec_Coord>=6999");
            int n7000After = ntAfter->GetEntries("32==nicell&&before.angle_diff_max<1&&Prec_Coord>=6999");
            printf("Number of 7000 GeV BEFORE alignment is %d, Number of 7000 GeV AFTER alignment is %d\n", n7000Before, n7000After);
            nPlot++;
        }
    }
    c->Print("momentum_output/momentum_before_after_icellmax32.pdf]");
    return nPlot;
}
