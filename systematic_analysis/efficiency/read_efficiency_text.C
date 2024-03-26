int read_efficiency_text(const char *filename)
{
    gStyle->SetNdivisions(506, "Y");
    TGaxis::SetMaxDigits(4);
    // gStyle->SetTitleOffset(0.8, "Y");
    const int plMin = 10;
    const int plMax = 348;
    const int nDivisionX = 1;
    const int nDivisionY = 1;
    TH2D *histEfficiency = new TH2D("histEfficiency", "Hit Efficiency;plate;reco", plMax - plMin + 1, plMin - 0.5, plMax + 0.5, 63, 1 - 0.5, 63 + 0.5);
    TH2D *histEfficiencyPlateX = new TH2D("histEfficiencyPlateX", "Hit efficiency for each plate and x;plate;x (#mum)", plMax - plMin + 1, plMin - 0.5, plMax + 0.5, 120, 5000, 125000);
    TH2D *histNPlateX = new TH2D("histNPlateX", "N;plate;x (#mum)", plMax - plMin + 1, plMin - 0.5, plMax + 0.5, 120, 5000, 125000);
    TH2D *histEfficiencyPlateY = new TH2D("histEfficiencyPlateY", "Hit efficiency for each plate and y;plate;y (#mum)", plMax - plMin + 1, plMin - 0.5, plMax + 0.5, 90, 5000, 95000);
    TH2D *histNPlateY = new TH2D("histNPlateY", "N;plate;y (#mum)", plMax - plMin + 1, plMin - 0.5, plMax + 0.5, 90, 5000, 95000);
    TH2D *histEfficiencyXYEachArea = new TH2D("histEfficiencyXYEachArea", "Hit Efficiency;x (#mum);y (#mum)", 9 * nDivisionX, -2500, 132500, 7 * nDivisionY, -2500, 102500);
    TH2D *histNXYEachArea = new TH2D("histNXYEachArea", "N;x (#mum);y (#mum)", 9 * nDivisionX, -2500, 132500, 7 * nDivisionY, -2500, 102500);

    double binWidth = 1000;
    double minX = -4000;
    double maxX = 134000;
    double minY = -4000;
    double maxY = 104000;
    int nBinsX = (maxX - minX) / binWidth;
    int nBinsY = (maxY - minY) / binWidth;
    TH2D *histEfficiencyXY = new TH2D("histEfficiencyXY", "Hit Efficiency;x (#mum);y (#mum)", nBinsX, minX, maxX, nBinsY, minY, maxY);
    TH2D *histNXY = new TH2D("histNXY", "N;x (#mum);y (#mum)", nBinsX, minX, maxX, nBinsY, minY, maxY);
    TH2D *histEfficiencyPID = new TH2D("histEfficiencyPID", "Hit efficiency;PID;subVolume*63+reco", 30, -0.5, 29.5, 13 * 63, 20 * 63 + 1 - 0.5, 32 * 63 + 63 + 0.5);
    TH1D *histEfficiencyDistribution = new TH1D("histEfficiencyDistribution", "hit efficiency;efficiency", 100, 0, 1);
    std::map<std::tuple<int, double, double>, std::pair<int, double>> mapEfficiency;
    int plateOld = -1;
    int subVolume = 0;
    FILE *fp = fopen(filename, "rt");
    char buf[256];
    for (; fgets(buf, sizeof(buf), fp);)
    {
        if (ferror(fp) || feof(fp))
            break;
        int plate, reco;
        double x, y, Efficiency;
        sscanf(buf, "%d %d %lf %lf %lf", &plate, &reco, &x, &y, &Efficiency);
        // if (x == 5000)
        //     continue;
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
        mapEfficiency[std::make_tuple(plate, x, y)] = std::make_pair(reco, Efficiency);
        int bin = histEfficiency->FindBin(plate, reco);
        histEfficiency->SetBinContent(bin, Efficiency);
        // printf("bin = %d",bin);
        // if (0 == histEfficiency->GetBinContent(bin))
        // {
        //     // printf("fill entry");
        //     histEfficiency->Fill(plate, reco, Efficiency);
        // }
        // int binXY = histEfficiencyXY->FindBin(x, y);

        if (plateOld + 1 != plate)
        {
            subVolume = plate / 10;
            // printf("%d\n",subVolume);
        }
        plateOld = plate;
        histEfficiencyPID->Fill(plate - (subVolume * 10 + 1), subVolume * 63 + reco, Efficiency);
    }
    // TFile *file = new TFile("treeEfficiency_30plates_effAreaEachPlate_division10_edge_area.root", "recreate");
    TTree *treeEfficiency = new TTree("treeEfficiency", "eff");
    int plate, reco;
    double x, y, efficiency;
    treeEfficiency->Branch("plate", &plate);
    treeEfficiency->Branch("reco", &reco);
    treeEfficiency->Branch("x", &x);
    treeEfficiency->Branch("y", &y);
    treeEfficiency->Branch("efficiency", &efficiency);
    for (auto itrMap = mapEfficiency.begin(); itrMap != mapEfficiency.end(); itrMap++)
    {
        plate = std::get<0>(itrMap->first);
        x = std::get<1>(itrMap->first);
        y = std::get<2>(itrMap->first);
        reco = itrMap->second.first;
        efficiency = itrMap->second.second;
        histEfficiencyPlateX->Fill(plate, x, efficiency);
        histNPlateX->Fill(plate, x);
        histEfficiencyPlateY->Fill(plate, y, efficiency);
        histNPlateY->Fill(plate, y);
        histEfficiencyXYEachArea->Fill(x, y, efficiency);
        histNXYEachArea->Fill(x, y);
        histEfficiencyXY->Fill(x, y, efficiency);
        histNXY->Fill(x, y);
        treeEfficiency->Fill();
    }
    // treeEfficiency->Write();

    int nBinsPlate = histEfficiency->GetXaxis()->GetNbins();
    int nBinsReco = histEfficiency->GetYaxis()->GetNbins();
    for (int iBinX = 1; iBinX <= nBinsPlate; iBinX++)
    {
        for (int iBinY = 1; iBinY <= nBinsReco; iBinY++)
        {
            double efficiency = histEfficiency->GetBinContent(iBinX, iBinY);
            if (efficiency != 0)
                histEfficiencyDistribution->Fill(efficiency);
        }
    }

    histEfficiency->SetMaximum(1.0);
    histEfficiency->SetMinimum(0.0);
    histEfficiency->Draw("colz");
    histEfficiency->SetStats(0);

    TCanvas *c0 = new TCanvas("c0","efficiency for each plate and x (or y)",1200,500);
    c0->Divide(2);
    c0->cd(1);
    histEfficiencyPlateX->Divide(histNPlateX);
    histEfficiencyPlateX->SetMaximum(1.0);
    histEfficiencyPlateX->SetMinimum(0.0);
    histEfficiencyPlateX->Draw("colz");
    histEfficiencyPlateX->SetStats(0);
    c0->cd(2);
    histEfficiencyPlateY->Divide(histNPlateY);
    histEfficiencyPlateY->SetMaximum(1.0);
    histEfficiencyPlateY->SetMinimum(0.0);
    histEfficiencyPlateY->Draw("colz");
    histEfficiencyPlateY->SetStats(0);

    TCanvas *c1 = new TCanvas();
    gPad->DrawFrame(5000, 5100, 125000, 94900, "Hit efficiency;x (#mum);y (#mum)");
    histEfficiencyXYEachArea->Divide(histNXYEachArea);
    histEfficiencyXYEachArea->SetMaximum(1.0);
    histEfficiencyXYEachArea->SetMinimum(0.0);
    histEfficiencyXYEachArea->Draw("samecolz");
    histEfficiencyXYEachArea->SetStats(0);
    gPad->RedrawAxis();

    TCanvas *c2 = new TCanvas();
    // gPad->DrawFrame(5000, 5100, 125000, 94900, "Hit efficiency;x (#mum);y (#mum)");
    histEfficiencyXY->Divide(histNXY);
    histEfficiencyXY->SetMaximum(1.0);
    histEfficiencyXY->SetMinimum(0.0);
    histEfficiencyXY->Draw("colz");
    histEfficiencyXY->SetStats(0);
    // gPad->RedrawAxis();

    TCanvas *c3 = new TCanvas();
    histEfficiencyPID->SetMaximum(1.0);
    histEfficiencyPID->SetMinimum(0.0);
    histEfficiencyPID->Draw("colz");
    histEfficiencyPID->SetStats(0);

    TCanvas *c4 = new TCanvas();
    histEfficiencyDistribution->Draw();

    fclose(fp);
    return histEfficiency->GetEntries();
}
