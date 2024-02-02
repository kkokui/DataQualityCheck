#include <TStyle.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TText.h>

int read_momentum_from_nt()
{
    gStyle->SetPadRightMargin(0.21);
    gStyle->SetPadLeftMargin(0.2);
    gStyle->SetPadBottomMargin(0.22);
    gStyle->SetPadTopMargin(0.02);
    gStyle->SetLabelOffset(0.01,"XYZ");
    gStyle->SetTitleOffset(1.3, "Y");
    gStyle->SetTitleOffset(1.0, "Z");
    gStyle->SetTitleSize(0.08, "XYZ");
    gStyle->SetLabelSize(0.06, "XYZ");
    gStyle->SetNdivisions(505, "xyz");

    TChain *ntBefore = new TChain("nt");
    ntBefore->Add("momentum_output/nt_before_align_icellmax32*");
    ntBefore->Draw("Prec_Coord:Prec_Coord>>hist(200,0,7500,200,0,7500)", "32==nicell&&angle_diff_max<1", "colz");
    TH2F *hist = (TH2F *)gDirectory->Get("hist");
    hist->SetStats(0);
    hist->SetTitle("Prec_Coord before alignment : Prec_Coord before alignment;P_{before} (GeV);P_{before} (GeV);entries");
    // c->SetLogz();
    // c->SetRealAspectRatio();
    // c->Print("momentum_output/momentum_before_after_icellmax32.pdf[");
    TCanvas *c = new TCanvas("c", "multiple plots", 750, 500);
    c->Divide(6, 5,0.005,0.005);
    int isubpad = 0;
    for (int iRobust = 10; iRobust >= 6; iRobust--)
    {
        isubpad++;
        c->cd(isubpad);
        TText *t = new TText(0.5,0,Form("robust factor %.1f",iRobust*0.1));
        t->SetNDC(1);
        t->SetTextAlign(21);
        t->SetTextSize(0.15);
        t->Draw();
    }
    isubpad++;

    double binWidthArr[] = {5000, 2000, 1000, 500};
    int nPlot = 0;
    for (int i = 0; i < 4; i++)
    {
        for (int iRobust = 10; iRobust >= 6; iRobust--)
        {
            isubpad++;
            c->cd(isubpad);

            gPad->SetLogz();

            TChain *ntAfter = new TChain("nt");
            int nFilesAdded = ntAfter->Add(Form("momentum_output/nt_after_align_binWidth%.0f_robustFactor%.1f_icellmax32_*", binWidthArr[i], iRobust * 0.1));
            if (0 == nFilesAdded)
            {
                printf("momentum_output/nt_after_align_binWidth%.0f_robustFactor%.1f_icellmax32_* does not exist\n", binWidthArr[i], iRobust * 0.1);
                continue;
            }
            ntAfter->AddFriend(ntBefore, "before");
            ntAfter->Draw(Form("Prec_Coord:before.Prec_Coord>>hist%d(200,0,7500,200,0,7500)", isubpad), "32==nicell&&before.angle_diff_max<1", "colz");
            TH2F *hist = (TH2F *)gDirectory->Get(Form("hist%d", isubpad));
            hist->SetStats(0);
            // hist->SetTitle(Form("momentum before and after alignment (binWidth%.0f_robustFactor%.1f);P_{before} (GeV);P_{after} (GeV);entries", binWidthArr[i], iRobust * 0.1));
            hist->SetTitle(";P before (GeV);P after (GeV);entries");
            // ntAfter->Draw("trid:before.trid","","colz"); // to check if the track ID is the same
            // c->Print("momentum_output/momentum_before_after_icellmax32.pdf");
            int n7000Before = ntAfter->GetEntries("32==nicell&&before.angle_diff_max<1&&before.Prec_Coord>=6999");
            int n7000After = ntAfter->GetEntries("32==nicell&&before.angle_diff_max<1&&Prec_Coord>=6999");
            printf("Number of 7000 GeV BEFORE alignment is %d, Number of 7000 GeV AFTER alignment is %d\n", n7000Before, n7000After);
            TText *tN = new TText(500,6000,Form("%d",n7000After));
            tN->SetTextSize(0.15);
            tN->Draw();
            nPlot++;
            // break;
        }
        isubpad++;
        c->cd(isubpad);
        TText *t = new TText(0.07, 0.05, Form("%.1f mm division", binWidthArr[i] / 1000));
        t->SetNDC(1);
        // t->SetTextAlign(12);
        t->SetTextSize(0.13);
        t->SetTextAngle(90);
        t->Draw();
        // break;
    }
    c->Print("momentum_output/momentum_before_after_icellmax32.pdf");
    return nPlot;
}
