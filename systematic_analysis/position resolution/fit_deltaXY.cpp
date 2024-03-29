#include <stdio.h>
#include <algorithm>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TF1.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TCut.h>

int sub_volume;
int reco;

double CalcMaxValue(int N, double *arr)
{
    double MaxValue = arr[0];
    for (int i = 0; i < N; i++)
    {
        if (MaxValue < arr[i])
        {
            MaxValue = arr[i];
        }
    }
    return MaxValue;
}
double CalcMinValue(int N, double *arr)
{
    double MinValue = arr[0];
    for (int i = 0; i < N; i++)
    {
        if (MinValue > arr[i])
        {
            MinValue = arr[i];
        }
    }
    return MinValue;
}
void calc_eff_area(TTree *tree, int reco, double &range_x_max, double &range_x_min, double &range_y_max, double &range_y_min)
{
    tree->Draw("y:x", "nseg>=15");
    range_x_max = CalcMaxValue(tree->GetSelectedRows(), tree->GetV2());
    range_x_min = CalcMinValue(tree->GetSelectedRows(), tree->GetV2());
    range_y_max = CalcMaxValue(tree->GetSelectedRows(), tree->GetV1());
    range_y_min = CalcMinValue(tree->GetSelectedRows(), tree->GetV1());
}
int calc_eff_area(TTree *tree, double &x_max, double &x_min, double &y_max, double &y_min)
{
    int nent = tree->GetEntriesFast();
    double x, y;
    int nseg;
    tree->SetBranchAddress("x", &x);
    tree->SetBranchAddress("y", &y);
    tree->SetBranchAddress("nseg", &nseg);
    int ntrk_cut = 0;
    int nseg_max = 30;
    std::vector<double> vec_x;
    std::vector<double> vec_y;
    while (1000 > ntrk_cut && nseg_max >= 5)
    {
        for (int ient = 0; ient < nent; ient++)
        {
            tree->GetEntry(ient);
            if (nseg_max == nseg)
            {
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
    x_max = vec_x.at(ntrk_cut - ntrk_cut * 5 / 1000);
    y_min = vec_y.at(ntrk_cut * 5 / 1000);
    y_max = vec_y.at(ntrk_cut - ntrk_cut * 5 / 1000);
    return 0;
}
void fit_deltaXY(TTree *tree, double division_center_x, double division_center_y, TCut cut, TCanvas *c1, FILE *ftxt, TString title)
{
    for (int i = 3; i <= 28; i++)
    {
        int pl = sub_volume * 10 + i;
        TCut plcut = Form("pl==%d", pl);
        double fit_range;
        double fit_center;
        TF1 *fit_func_x = new TF1("fit_func_x", "gaus", -2, 2);
        TF1 *fit_func_y = new TF1("fit_func_y", "gaus", -2, 2);
        fit_func_x->SetParameter(2, 0.3);
        fit_func_y->SetParameter(2, 0.3);
        fit_func_x->SetParLimits(2, 0., 1.5);
        fit_func_y->SetParLimits(2, 0., 1.5);
        tree->Draw("deltaX>>hx", cut + plcut);
        if(tree->GetSelectedRows()==0) continue;
        TH1D *hx = (TH1D *)gDirectory->Get("hx");
        // tree->Draw("deltaX>>hx_fit(100,-2,2)","abs(tx+0.01)<0.01&&abs(ty)<0.01");
        tree->Draw("deltaX>>hx_fit(100,-3,3)", cut + plcut);
        TH1D *hx_fit = (TH1D *)gDirectory->Get("hx_fit");
        fit_range = hx_fit->GetRMS();
        fit_center = hx_fit->GetMean();
        hx_fit->Fit(fit_func_x, "QB", "", fit_center-fit_range,fit_center+ fit_range);
        gStyle->SetOptFit();
        // c1->Print("0.9_peak_test/"+title+".pdf");

        tree->Draw("deltaY>>hy", cut + plcut);
        TH1D *hy = (TH1D *)gDirectory->Get("hy");
        // tree->Draw("deltaY>>hy_fit(100,-2,2)","abs(tx+0.01)<0.01&&abs(ty)<0.01");
        tree->Draw("deltaY>>hy_fit(100,-3,3)", cut + plcut);
        TH1D *hy_fit = (TH1D *)gDirectory->Get("hy_fit");
        fit_range = hy_fit->GetRMS();
        fit_center = hy_fit->GetMean();
        // fit_center = 0;
        hy_fit->Fit(fit_func_y, "QB", "", fit_center-fit_range, fit_center+fit_range);

        fprintf(ftxt, "%d %d %d %.0f %f %f %f %f %f %f %f %f\n", sub_volume, reco, pl, hx_fit->GetEntries(), division_center_x, division_center_y, hx->GetMean(), hx->GetRMS(), fit_func_x->GetParameter(2), hy->GetMean(), hy->GetRMS(), fit_func_y->GetParameter(2));
    }
}
int main(int argc, char *argv[])
{
    if (argc < 5)
    {
        printf("usage: ./fit_deltaXY deltaXY/tree_BeforeAlign.root title sub_volume reco\n");
        return 1;
    }
    gSystem->Load("libTree");
    TFile *fin = new TFile(argv[1]);
    if(fin->IsZombie())
    {
        return 0;
    }
    TString title = argv[2];
    sub_volume = atoi(argv[3]);
    reco = atoi(argv[4]);
    double Xcenter = (reco - 1) % 9 * 15000 + 5000;
    double Ycenter = (reco - 1) / 9 * 15000 + 5000;
    TCanvas *c1 = new TCanvas();
    TTree *tree = (TTree *)gDirectory->Get("tree");
    double range_x_max = 125000;
    double range_x_min = 11000;
    double range_y_max = 95000;
    double range_y_min = 5000;
    // calc_eff_area(tree, reco, range_x_max,range_x_min,range_y_max,range_y_min);
    // calc_eff_area(tree, range_x_max, range_x_min, range_y_max, range_y_min);
    TCut rangeXY = Form("x>%.0f&&x<%.0f&&y>%.0f&&y<%.0f", range_x_min, range_x_max, range_y_min, range_y_max);
    const int ndivision1D = 1;
    TCut divisionX[ndivision1D];
    TCut divisionY[ndivision1D];
    // c1->Print("0.9_peak_test/"+title+".pdf[");
    FILE *ftxt = fopen(Form("position_resolution_txt_output/%s_reco%02d.txt", title.Data(), reco), "a"); // Make txt file for efficiency
    // FILE *ftxt = fopen(Form("0.9_peak_test/%s.txt", title.Data(), reco), "a"); // Make txt file for efficiency
    for (int i = 0; i < ndivision1D; i++)
    {
        for (int j = 0; j < ndivision1D; j++)
        {
            double xmin = Xcenter - 8500 + 17000. / ndivision1D * i;
            double xmax = Xcenter - 8500 + 17000. / ndivision1D * (i + 1);
            divisionX[i] = Form("%.0f<x&&x<%.0f", xmin, xmax);
            double ymin = Ycenter - 8500 + 17000. / ndivision1D * j;
            double ymax = Ycenter - 8500 + 17000. / ndivision1D * (j + 1);
            divisionY[j] = Form("%.0f<y&&y<%.0f", ymin, ymax);
            double division_center_x = (xmin + xmax) / 2;
            double division_center_y = (ymin + ymax) / 2;
            fit_deltaXY(tree, division_center_x, division_center_y, rangeXY + divisionX[i] + divisionY[j], c1, ftxt,title);
        }
    }
    fclose(ftxt);

    // c1->Print("0.9_peak_test/"+title+".pdf]");
}
