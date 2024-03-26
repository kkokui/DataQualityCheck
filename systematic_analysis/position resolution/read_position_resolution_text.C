int read_position_resolution_text(const char *filename)
{
    gStyle->SetNdivisions(506, "Y");
    TGaxis::SetMaxDigits(4);
    gStyle->SetOptStat("emr");
    // gStyle->SetTitleOffset(0.8, "Y");
    const int plMin = 10;
    const int plMax = 348;
    const int nDivisionX = 1;
    const int nDivisionY = 1;
    TH2D *histPositionResolutionX = new TH2D("histPositionResolutionX", "position resolution x;plate;reco", plMax - plMin + 1, plMin - 0.5, plMax + 0.5, 63, 1 - 0.5, 63 + 0.5);
    TH2D *histPositionResolutionY = new TH2D("histPositionResolutionY", "position resolution y;plate;reco", plMax - plMin + 1, plMin - 0.5, plMax + 0.5, 63, 1 - 0.5, 63 + 0.5);
    TH2D *histPositionResolutionXY = new TH2D("histPositionResolutionXY", "position resolution;x (#mum);y (#mum)", 9 * nDivisionX, -2500, 132500, 7 * nDivisionY, -2500, 102500);
    TH2D *histPositionResolutionPID = new TH2D("histPositionResolutionPID", "position_resolution;PID;subVolume*63+reco", 30, -0.5, 29.5, 13 * 63, 20 * 63 + 1 - 0.5, 32 * 63 + 63 + 0.5);
    TH1D *histPositionResolutionXDistribution = new TH1D("histPositionResolutionXDistribution", "position resolution x;position resolution (#mum);", 100, 0, 1.6);
    TH1D *histPositionResolutionYDistribution = new TH1D("histPositionResolutionYDistribution", "position resolution y;position resolution (#mum);", 100, 0, 1.6);
    TH1D *histPositionResolutionDistribution = new TH1D("histPositionResolutionDistribution", "position resolution x and y;position resolution (#mum);", 100, 0, 1.6);
    int plateOld = -1;
    int subVolume = 0;
    FILE *fp = fopen(filename, "rt");
    char buf[256];
    std::map<std::tuple<int, int, int, double, double>, std::pair<double, double>> mapResolution;
    for (; fgets(buf, sizeof(buf), fp);)
    {
        if (ferror(fp) || feof(fp))
            break;
        int subVolume, reco, plate, entries;
        double x, y, meanX, rmsX, sigmaX, meanY, rmsY, sigmaY;
        sscanf(buf, "%d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf", &subVolume, &reco, &plate, &entries, &x, &y, &meanX, &rmsX, &sigmaX, &meanY, &rmsY, &sigmaY);
        if (entries < 2)
            continue;
        mapResolution[std::make_tuple(reco, subVolume, plate, x, y)] = std::make_pair(sigmaX, sigmaY);
    }
    for (auto itrMap = mapResolution.begin(); itrMap != mapResolution.end(); itrMap++)
    {
        double reco = std::get<0>(itrMap->first);
        double subVolume = std::get<1>(itrMap->first);
        double plate = std::get<2>(itrMap->first);
        double x = std::get<3>(itrMap->first);
        double y = std::get<4>(itrMap->first);
        double sigmaX = itrMap->second.first;
        double sigmaY = itrMap->second.second;
        int bin = histPositionResolutionX->FindBin(plate, reco);
        histPositionResolutionX->SetBinContent(bin, sigmaX);

        bin = histPositionResolutionY->FindBin(plate, reco);
        histPositionResolutionY->SetBinContent(bin, sigmaY);
        // printf("bin = %d",bin);
        // if (0 == histPositionResolution->GetBinContent(bin))
        // {
        //     // printf("fill entry");
        //     histPositionResolution->Fill(plate, reco, PositionResolution);
        // }
        if (290 == plate)
        {
            int binXY = histPositionResolutionXY->FindBin(x, y);
            if (0 == histPositionResolutionXY->GetBinContent(binXY))
            {
                histPositionResolutionXY->Fill(x, y, sigmaX);
            }
        }
        // if(plateOld+1!=plate)
        // {
        //     subVolume=plate/10;
        //     // printf("%d\n",subVolume);
        // }
        // plateOld = plate;
        // histPositionResolutionPID->Fill(plate-(subVolume*10+1),subVolume*63+reco,PositionResolution);
    }

    double ratioLessThan0point6 = 0;

    int nBinsX = histPositionResolutionX->GetXaxis()->GetNbins();
    int nBinsY = histPositionResolutionX->GetYaxis()->GetNbins();
    for (int iBinX = 1; iBinX <= nBinsX; iBinX++)
    {
        for (int iBinY = 1; iBinY <= nBinsY; iBinY++)
        {
            double resolution = histPositionResolutionX->GetBinContent(iBinX, iBinY);
            if (resolution != 0)
            {
                histPositionResolutionXDistribution->Fill(resolution);
                histPositionResolutionDistribution->Fill(resolution);
                if (resolution < 0.6)
                    ratioLessThan0point6 += 1;
            }
        }
    }

    nBinsX = histPositionResolutionY->GetXaxis()->GetNbins();
    nBinsY = histPositionResolutionY->GetYaxis()->GetNbins();
    for (int iBinX = 1; iBinX <= nBinsX; iBinX++)
    {
        for (int iBinY = 1; iBinY <= nBinsY; iBinY++)
        {
            double resolution = histPositionResolutionY->GetBinContent(iBinX, iBinY);
            if (resolution != 0)
            {
                histPositionResolutionYDistribution->Fill(resolution);
                histPositionResolutionDistribution->Fill(resolution);
                if (resolution < 0.6)
                    ratioLessThan0point6 += 1;
            }
        }
    }

    ratioLessThan0point6 /= histPositionResolutionDistribution->GetEntries();
    printf("ratio of the resolution less than 0.6 is %f\n", ratioLessThan0point6);

    histPositionResolutionX->SetMaximum(1.5);
    histPositionResolutionX->SetMinimum(0.0);
    histPositionResolutionX->Draw("colz");
    histPositionResolutionX->SetStats(0);
    TCanvas *c1 = new TCanvas();
    histPositionResolutionY->SetMaximum(1.5);
    histPositionResolutionY->SetMinimum(0.0);
    histPositionResolutionY->Draw("colz");
    histPositionResolutionY->SetStats(0);

    TCanvas *c2 = new TCanvas();
    gPad->DrawFrame(5000, 5100, 125000, 94900, "position_resolution;x (#mum);y (#mum)");
    histPositionResolutionXY->Draw("samecolz");
    histPositionResolutionXY->SetStats(0);
    gPad->RedrawAxis();

    TCanvas *c3 = new TCanvas();
    histPositionResolutionPID->SetMaximum(1.0);
    histPositionResolutionPID->SetMinimum(0.0);
    histPositionResolutionPID->Draw("colz");
    histPositionResolutionPID->SetStats(0);

    TCanvas *c4 = new TCanvas();
    histPositionResolutionXDistribution->Draw();

    TCanvas *c5 = new TCanvas();
    histPositionResolutionYDistribution->Draw();

    TCanvas *c6 = new TCanvas();
    histPositionResolutionDistribution->Draw();
    fclose(fp);
    return histPositionResolutionX->GetEntries();
}
