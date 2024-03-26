#include <TH2.h>
#include <TTree.h>
#include <TStyle.h>
#include <TGaxis.h>
#include <TCanvas.h>
#include <TPaletteAxis.h>
#include <TGraph.h>
#include <TText.h>

double lumi = 9.523; // fb^-1
void fill_with_density(TH2D *h2, int zone, double bottomLeftX, double bottomLeftY)
{
    TTree *t = new TTree("t", "track_density.txt");
    // t->ReadFile(Form("15plates_and_30plates_zone%d.txt", zone), "sub_v/I:reco/I:Xcenter/F:Ycenter/F:track_density/D");
    t->ReadFile(Form("track_density_zone%d.txt", zone), "sub_v/I:reco/I:Xcenter/F:Ycenter/F:track_density/D");
    int reco;
    int sub_v;
    float Xcenter;
    float Ycenter;
    double track_density;
    t->SetBranchAddress("reco", &reco);
    t->SetBranchAddress("sub_v", &sub_v);
    t->SetBranchAddress("Xcenter", &Xcenter);
    t->SetBranchAddress("Ycenter", &Ycenter);
    t->SetBranchAddress("track_density", &track_density);
    for (int ient = 0; ient < t->GetEntriesFast(); ient++)
    {
        t->GetEntry(ient);
        if (sub_v / 10 != 18)
            continue;
        int bin = h2->FindBin(Xcenter + bottomLeftX, Ycenter + bottomLeftY);
        h2->SetBinContent(bin, track_density / lumi);
    }
}

void plot_track_density()
{
    gStyle->SetPadRightMargin(0.2);
    // gStyle->SetNdivisions(508, "Y");
    gStyle->SetNdivisions(508, "Z");
    gStyle->SetTitleOffset(1.5, "z");
    TGaxis::SetMaxDigits(4);
    TH2D *h2;
    int zone = 1;
    // h2 = new TH2D("h2", Form("track density (zone%d);x (#mum);y (#mum);tracks/cm^{2}/fb^{-1}", zone), 9, -2500, 132500, 7, -2500, 102500);
    h2 = new TH2D("h2", "track density;x (#mum);y (#mum);tracks/cm^{2}/fb^{-1}", 17, -2500, 252500, 13, -2500, 192500);
    fill_with_density(h2, 3, 0, 90000);
    fill_with_density(h2, 1, 0, 0);
    fill_with_density(h2, 4, 120000, 90000);
    h2->SetMinimum(h2->GetMinimum(0));
    // h2->SetMinimum(14422.177);
    // h2->SetMaximum(16428.558);
    TCanvas *c1 = new TCanvas();
    // gPad->DrawFrame(0, 0, 130000, 100000, Form("track density (zone%d);x (#mum);y (#mum);tracks/cm^{2}", zone));
    // gPad->DrawFrame(5000, 5000, 125000, 95000, Form("track density (zone%d);x (#mum);y (#mum);tracks/cm^{2}", zone));
    TH1F *frame = gPad->DrawFrame(5000, 5000, 245000, 185000, "track density;x (#mum);y (#mum);tracks/cm^{2}");
    frame->GetXaxis()->SetNdivisions(509);
    frame->GetYaxis()->SetNdivisions(507);
    h2->Draw("samecolz");
    // h2->Draw("colz");
    h2->SetStats(0);
    gPad->Update();
    TPaletteAxis *palette = (TPaletteAxis *)h2->GetListOfFunctions()->FindObject("palette");
    palette->SetX1NDC(0.85);
    palette->SetX2NDC(0.90);
    gPad->RedrawAxis();
    TGraph *pointLOS = new TGraph;
    pointLOS->SetPoint(0, 115000, 109000);
    pointLOS->SetMarkerStyle(20);
    pointLOS->Draw("P");
    TText *textLOS = new TText(115000, 109000, "LOS ");
    textLOS->SetTextAlign(33);
    textLOS->Draw();
    c1->SetRealAspectRatio();
    // gStyle->SetPadRightMargin(0.15);
    TTree *t = new TTree("t_trackDensityVsPlate", "track_density.txt");
    t->ReadFile("track_density_zone4.txt", "sub_v/I:reco/I:Xcenter/F:Ycenter/F:track_density/D");
    TCanvas *c2 = new TCanvas();
    gPad->DrawFrame(0,6000,350,16500);
    t->SetMarkerStyle(6);
    t->Draw(Form("track_density/%.3f:sub_v",lumi),"","same");
    // gPad->DrawFrame(0,)
    // c1->Print("track_density_2D.pdf");
}
