#include <stdio.h>
#include <map>
#include <TH3.h>
#include <TProfile2D.h>
#include <TCanvas.h>
#include <TStyle.h>
void fill_hist(TH3D *h3Efficiency, int zone, double bottomLeftX, double bottomLeftY)
{
    FILE *fp = fopen(Form("efficiency_zone%d.txt", zone), "rt");
    char buf[256];
    for (; fgets(buf, sizeof(buf), fp);)
    {
        if (ferror(fp) || feof(fp))
            break;

        int plate, reco;
        double x, y, efficiency;
        sscanf(buf, "%d %d %lf %lf %lf", &plate, &reco, &x, &y, &efficiency);
        double areaCenterX = (reco - 1) % 9 * 15000 + 5000;
        double areaCenterY = (reco - 1) / 9 * 15000 + 5000;
        if (x == areaCenterX - 8500)
            continue;
        if (x == areaCenterX + 8500)
            continue;
        if (y == areaCenterY - 8500)
            continue;
        if (y == areaCenterY + 8500)
            continue;

        int bin = h3Efficiency->FindBin(x + bottomLeftX, y + bottomLeftY, plate);
        if (efficiency != 0)
        {
            h3Efficiency->SetBinContent(bin, efficiency);
        }
    }
}

void read_all_efficiency_text()
{
    const int plMin = 10; // for zones 3 or 4
    const int plMax = 348; // for zones 3 or 4
    int npl = plMax - plMin + 1;
    const double binWidth = 1000;

    // const double minX = -4000; //for multiple zones (3 and 4)
    // const double maxX = 254000; //for multiple zones (3 and 4)
    // const double minY = 87500; //for multiple zones (3 and 4)
    // const double maxY = 192500; //for multiple zones (3 and 4)

    const double minX = -4000; // for one zone
    const double maxX = 134000; // for one zone
    const double minY = -4000; // for one zone
    const double maxY = 104000; // for one zone

    int nBinsX = (maxX - minX) / binWidth;
    int nBinsY = (maxY - minY) / binWidth;
    TH3D *h3Efficiency = new TH3D("h3Efficiency", "efficiency;x;y;plate", nBinsX, minX, maxX, nBinsY, minY, maxY, npl, plMin - 0.5, plMax + 0.5);

    // fill_hist(h3Efficiency, 3, 0, 90000); // for multiple zones
    // fill_hist(h3Efficiency, 4, 122000, 90000); // for multiple zones

    fill_hist(h3Efficiency, 4, 0, 0);

    TProfile2D *prof2EfficiencyMap = new TProfile2D("prof2EfficiencyMap", "efficiency;x (#mum);y (#mum);efficiency", nBinsX, minX, maxX, nBinsY, minY, maxY, "");
    TProfile2D *prof2EfficiencyPlateVsX = new TProfile2D("prof2EfficiencyPlateVsX", "efficiency for each plate and x;plate;x (#mum);efficiency", npl, plMin - 0.5, plMax + 0.5, nBinsX, minX, maxX, "");
    TProfile2D *prof2EfficiencyPlateVsY = new TProfile2D("prof2EfficiencyPlateVsY", "efficiency for each plate and y;plate;y (#mum);efficiency", npl, plMin - 0.5, plMax + 0.5, nBinsY, minY, maxY, "");
    for (int iBinPl = 1; iBinPl <= npl; iBinPl++)
    {
        for (int iBinX = 1; iBinX <= nBinsX; iBinX++)
        {
            for (int iBinY = 1; iBinY <= nBinsY; iBinY++)
            {
                double efficiency = h3Efficiency->GetBinContent(iBinX, iBinY, iBinPl);
                double x = h3Efficiency->GetXaxis()->GetBinCenter(iBinX);
                double y = h3Efficiency->GetYaxis()->GetBinCenter(iBinY);
                int plate = h3Efficiency->GetZaxis()->GetBinCenter(iBinPl);
                prof2EfficiencyMap->Fill(x, y, efficiency);
                prof2EfficiencyPlateVsX->Fill(plate,x,efficiency);
                prof2EfficiencyPlateVsY->Fill(plate,y,efficiency);
            }
        }
    }
    gStyle->SetPadRightMargin(0.15);
    TCanvas *c0 = new TCanvas("c0","efficiency map");
    prof2EfficiencyMap->Draw("colz");
    prof2EfficiencyMap->SetStats(0);
    prof2EfficiencyMap->SetMinimum(0.0001);
    prof2EfficiencyMap->SetMaximum(1);

    TCanvas *c1 = new TCanvas("c1", "efficiency for each plate and x (or y)", 1200, 500);
    c1->Divide(2);
    c1->cd(1);
    prof2EfficiencyPlateVsX->Draw("colz");
    prof2EfficiencyPlateVsX->SetStats(0);
    prof2EfficiencyPlateVsX->SetMinimum(0.0001);
    prof2EfficiencyPlateVsX->SetMaximum(1);
    c1->cd(2);
    prof2EfficiencyPlateVsY->Draw("colz");
    prof2EfficiencyPlateVsY->SetStats(0);
    prof2EfficiencyPlateVsY->SetMinimum(0.0001);
    prof2EfficiencyPlateVsY->SetMaximum(1);
}