#include <stdio.h>
#include <algorithm>

#include <TTree.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TCut.h>
#include <TH2.h>

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        printf("usage: ./getEntriesWithCut filename\n");
        return 1;
    }
    const char *filename = argv[1];
    int sub_v, reco, sub_v2;
    float Xcenter, Ycenter;
    char *p;
    char q[100];
    strcpy(q, filename);
    p = strtok(q, "pl_reco/ToReprocess");
    sscanf(p, "%d", &sub_v);
    p = strtok(NULL, "pl_reco/ToReprocess");
    sscanf(p, "%d", &sub_v2);
    p = strtok(NULL, "pl_reco/ToReprocess");
    sscanf(p, "%d", &reco);
    p = strtok(NULL, "pl_reco/ToReprocess");
    sscanf(p, "%f", &Xcenter);
    p = strtok(NULL, "pl_reco/ToReprocess");
    sscanf(p, "%f", &Ycenter);
    // TFile *fin = new TFile(Form("/mnt/Disk4/F222/zone3/%s", filename));
    TFile *fin = new TFile(Form("/mnt/Disk5/F222/zone3/%s", filename));
    if (fin->IsZombie())
    {
        return 0;
    }
    TTree *tracks = (TTree *)gDirectory->Get("tracks");
    TCanvas *c1 = new TCanvas();
    // tracks->Draw("s[0].eY:s[0].eX", "s[0].ePID<=5&&nseg>=20");
    int Nselected = 0;
    int iseg = 14;
    double edge_x_max = 125000;
    double edge_x_min = 11000;
    double edge_y_max = 95000;
    double edge_y_min = 5000;
    TCut edgeCut = Form("s[0].eX>%.0f&&s[0].eX<%.0f&&s[0].eY>%.0f&&s[0].eY<%.0f", edge_x_min, edge_x_max, edge_y_min, edge_y_max);
    double no_overlap_x_max = Xcenter + 7500;
    double no_overlap_x_min = Xcenter - 7500;
    double no_overlap_y_max = Ycenter + 7500;
    double no_overlap_y_min = Ycenter - 7500;
    TCut noOverlapXY = Form("s[0].eX>%.0f&&s[0].eX<%.0f&&s[0].eY>%.0f&&s[0].eY<%.0f", no_overlap_x_min, no_overlap_x_max, no_overlap_y_min, no_overlap_y_max);
    while (Nselected < 100)
    {
        tracks->Draw("s[0].eY:s[0].eX", Form("s[%d].ePID<=14", iseg) + edgeCut + noOverlapXY);
        Nselected = tracks->GetSelectedRows();
        iseg--;
        if (iseg < 5)
            break;
    }
    // c1->Print("test.pdf[");
    // c1->Print("test.pdf");
    double range_x_max = no_overlap_x_max;
    double range_x_min = no_overlap_x_min;
    double range_y_max = no_overlap_y_max;
    double range_y_min = no_overlap_y_min;
    int selectedRows = tracks->GetSelectedRows();
    // printf("selectedRows is %d\n",selectedRows);
    if (0 != selectedRows)
    {
        double *arrayX = tracks->GetV2();
        std::vector<double> vectorX(arrayX, arrayX + selectedRows);
        double *arrayY = tracks->GetV1();
        std::vector<double> vectorY(arrayY, arrayY + selectedRows);
        std::sort(vectorX.begin(), vectorX.end());
        std::sort(vectorY.begin(), vectorY.end());
        range_x_max = vectorX.at(selectedRows - 1 - selectedRows * 50 / 1000);
        range_x_min = vectorX.at(selectedRows * 50 / 1000);
        range_y_max = vectorY.at(selectedRows - 1 - selectedRows * 50 / 1000);
        range_y_min = vectorY.at(selectedRows * 50 / 1000);
    }
    TCut rangeXY = Form("s[0].eX>%f&&s[0].eX<%f&&s[0].eY>%f&&s[0].eY<%f", range_x_min, range_x_max, range_y_min, range_y_max);
    TCut muonCut = "s[0].ePID<=4&&nseg>=4";

    // tracks->Draw("s[0].eY:s[0].eX",muonCut,"colz");
    // TBox *box = new TBox(range_x_min,range_y_min,range_x_max,range_y_max);
    // box->SetLineColor(1);
    // box->SetFillStyle(0);
    // box->Draw("same");
    // c1->Print("test.pdf");

    tracks->SetAlias("tx", "(s[nseg-1].eX-s[0].eX)/(s[nseg-1].eZ-s[0].eZ)");
    tracks->SetAlias("ty", "(s[nseg-1].eY-s[0].eY)/(s[nseg-1].eZ-s[0].eZ)");
    tracks->Draw("ty:tx>>h(200,-0.02,0.02,200,-0.02,0.02)", muonCut + rangeXY, "colz");
    TH2D *h = (TH2D *)gDirectory->Get("h");
    // c1->Print("test.pdf");
    // c1->Print("test.pdf]");
    int locmax, locmay, locmaz;
    h->GetMaximumBin(locmax, locmay, locmaz);
    double beamCenterX = h->GetXaxis()->GetBinCenter(locmax);
    double beamCenterY = h->GetYaxis()->GetBinCenter(locmay);
    // double anglecut = 0.01;
    double anglecut = 0.1;
    TCut angleCenterCut = Form("abs(tx-%f)<%f&&abs(ty-%f)<%f", beamCenterX, anglecut, beamCenterY, anglecut);

    int n = tracks->GetEntries(muonCut + rangeXY + angleCenterCut);
    // int n = tracks->GetEntries(muonCut + rangeXY);
    double area = (range_x_max - range_x_min) * (range_y_max - range_y_min) / 10000 / 10000;
    printf("%d, %f\n", n, area);
    FILE *fp = fopen("30plates_zone3_xMin11000_AngleCut0.1.txt", "at");
    fprintf(fp, "%d %d %.0f %.0f %f\n", sub_v, reco, Xcenter, Ycenter, n / area);
    fclose(fp);
    printf("%d %d %.0f %.0f %f\n", sub_v, reco, Xcenter, Ycenter, n / area);
    return 0;
}
