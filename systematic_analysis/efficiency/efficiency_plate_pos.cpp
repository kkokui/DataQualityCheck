#include <stdio.h>
#include <fstream>
#include <sstream>
#include <EdbDataSet.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TH1.h>
#include <TGraphAsymmErrors.h>

int calc_eff_area(EdbPVRec *pvr, double &x_max, double &x_min, double &y_max, double &y_min)
{
    int ntrk = pvr->Ntracks();
    int ntrk_cut = 0;
    int nseg_max = 30;
    std::vector<double> vec_x;
    std::vector<double> vec_y;
    while (1000 > ntrk_cut && nseg_max >= 5)
    {
        for (int itrk = 0; itrk < ntrk; itrk++)
        {
            EdbTrackP *t = pvr->GetTrack(itrk);
            if (nseg_max == t->N())
            {
                EdbSegP *s0 = t->GetSegmentFirst();
                double x = s0->X();
                double y = s0->Y();
                vec_x.push_back(x);
                vec_y.push_back(y);
            }
        }
        ntrk_cut = vec_x.size();
        nseg_max--;
    }
    std::sort(vec_x.begin(), vec_x.end());
    std::sort(vec_y.begin(), vec_y.end());
    x_min = vec_x.at(ntrk_cut * 5 / 1000);
    x_max = vec_x.at(ntrk_cut - 1 - ntrk_cut * 5 / 1000);
    y_min = vec_y.at(ntrk_cut * 5 / 1000);
    y_max = vec_y.at(ntrk_cut - 1 - ntrk_cut * 5 / 1000);
    if (x_min < 5000)
    {
        x_min = 5000;
    }
    if (y_min < 5000)
    {
        y_min = 5000;
    }
    if (x_max > 125000)
    {
        x_max = 125000;
    }
    if (y_max > 95000)
    {
        y_max = 95000;
    }
    return 0;
}

int main(int argc, char *argv[])
{
    if (argc < 5)
    {
        printf("Usage : ./efficiency_plate_pos linked_tracks.root title reco binWidth\n");
        return 1;
    }
    TString filename_linked_tracks = argv[1];
    FILE *ftest = fopen(Form("%s",filename_linked_tracks.Data()),"r");
    if(ftest==NULL)
    {
        printf("%s does not exist\n",filename_linked_tracks.Data());
        return 1;
    }
    fclose(ftest);

    TString title = argv[2]; // used for title of histograms and name of output file
    int reco = atoi(argv[3]);
    int binWidth = atoi(argv[4]);
    double Xcenter = (reco - 1) % 9 * 15000 + 5000;
    double Ycenter = (reco - 1) / 9 * 15000 + 5000;
    double range = 9000;
    // TString cut = "Entry$<5000";
    TString cut = "nseg>=5";

    EdbDataProc *dproc = new EdbDataProc;
    EdbPVRec *pvr = new EdbPVRec;
    dproc->ReadTracksTree(*pvr, filename_linked_tracks, cut);

    int plMin = pvr->GetPatternByPID(0)->Plate();
    int plMax = plMin + pvr->Npatterns() - 1;
    int npl = plMax - plMin + 1;
    int ntrk = pvr->Ntracks();

    int pl, hitsOnThePlate;

    double x_min;
    double x_max;
    double y_min;
    double y_max;
    // calc_eff_area(pvr, x_max, x_min, y_max, y_min);

    double x1, y1, x2, y2;

    TCanvas *c = new TCanvas();
    // c->Print(Form("test.pdf[", title.Data()));
    for (int iplate = plMin + 2; iplate <= plMax - 2; iplate++)
    {
        // x_min = 5000;
        // x_max = 125000;
        // y_min = 5000;
        // y_max = 95000;

        // TEfficiency *eff2D = new TEfficiency("eff2D", Form("Efficiency for each plate (%s);plate;efficiency", title.Data()), Ndivision1D, Xcenter - range, Xcenter + range, Ndivision1D, Ycenter - range, Ycenter + range);
        int nBins = range*2/binWidth;
        TEfficiency *eff2D = new TEfficiency("eff2D", Form("Efficiency for each plate (%s);plate;efficiency", title.Data()), nBins, Xcenter - range, Xcenter + range, nBins, Ycenter - range, Ycenter + range);
        // if (iplate != 10)
        //    continue;

        // std::vector<double> vecX3, vecY3;
        // for (int itrk = 0; itrk < ntrk; itrk++)
        // {
        //     EdbTrackP *t = pvr->GetTrack(itrk);
        //     int nseg = t->N();
        //     EdbSegP *s0 = t->GetSegmentFirst();
        //     double x0 = s0->X();
        //     double y0 = s0->Y();
        //     if (x0 > x_max)
        //         continue;
        //     if (x0 < x_min)
        //         continue;
        //     if (y0 > y_max)
        //         continue;
        //     if (y0 < y_min)
        //         continue;
        //     for (int iseg = 0; iseg < nseg; iseg++)
        //     {
        //         EdbSegP *s = t->GetSegment(iseg);
        //         if (iplate == s->Plate())
        //         {
        //             vecX3.push_back(s->X());
        //             vecY3.push_back(s->Y());
        //             break;
        //         }
        //     }
        // }
        // int ntrkParPlate = vecX3.size();
        // if (0 != ntrkParPlate)
        // {
        //     std::sort(vecX3.begin(), vecX3.end());
        //     std::sort(vecY3.begin(), vecY3.end());
        //     x_min = vecX3.at(ntrkParPlate * 5 / 1000);
        //     x_max = vecX3.at(ntrkParPlate - 1 - ntrkParPlate * 5 / 1000);
        //     y_min = vecY3.at(ntrkParPlate * 5 / 1000);
        //     y_max = vecY3.at(ntrkParPlate - 1 - ntrkParPlate * 5 / 1000);
        // }

        for (int itrk = 0; itrk < ntrk; itrk++)
        {
            EdbTrackP *t = pvr->GetTrack(itrk);
            int nseg = t->N();
            EdbSegP *s0 = t->GetSegmentFirst();
            double x0 = s0->X();
            double y0 = s0->Y();
            
            // if (x0 > x_max)
            //     continue;
            // if (x0 < x_min)
            //     continue;
            // if (y0 > y_max)
            //     continue;
            // if (y0 < y_min)
            //     continue;

            int counts = 0;
            hitsOnThePlate = 0;
            for (int iseg = 0; iseg < nseg; iseg++)
            {
                EdbSegP *s = t->GetSegment(iseg);

                if (s->Plate() == iplate - 2 || s->Plate() == iplate + 2)
                    counts++;
                if (s->Plate() == iplate - 1)
                {
                    x1 = s->X();
                    y1 = s->Y();
                    counts++;
                }
                if (s->Plate() == iplate + 1)
                {
                    x2 = s->X();
                    y2 = s->Y();
                    counts++;
                }
                if (s->Plate() == iplate)
                {
                    hitsOnThePlate = 1;
                }
            }
            if (counts == 4)
            {
                double x = (x1 + x2) / 2;
                double y = (y1 + y2) / 2;
                eff2D->Fill(hitsOnThePlate, x, y);
            }
        }
        eff2D->Draw("colz");
        c->Update();
        // c->Print(Form("test.pdf", title.Data()));
        FILE *ftxt = fopen(Form("efficiency_txt_output/%s_reco%02d.txt", title.Data(), reco), "a"); // Make txt file for efficiency
        if(ftxt==NULL)
        {
            printf("failed to open output file\n");
            return 1;
        }
        TH2D *h2 = (TH2D *)eff2D->GetPaintedHistogram();
        for (int ix = 1; ix <= nBins; ix++)
        {
            for (int iy = 1; iy <= nBins; iy++)
            {
                int gbin = h2->GetBin(ix, iy);
                double x = h2->GetXaxis()->GetBinCenter(ix);
                double y = h2->GetYaxis()->GetBinCenter(iy);
                double eff = h2->GetBinContent(gbin);
                fprintf(ftxt, "%d %d %f %f %f\n", iplate, reco, x, y, eff);
            }
        }
        fclose(ftxt);
        delete eff2D;
    }
    // c->Print(Form("test.pdf]", title.Data()));

    return 0;
}
