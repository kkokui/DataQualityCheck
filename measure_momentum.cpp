#include <stdio.h>

#include <FnuMomCoord.hpp>
#include <FnuDivideAlign.h>
#include <FnuQualityCheck.h>

// #include <EdbDataSet.h>

// #include <TString.h>

int main(int argc, char *argv[])
{
    if (argc < 8)
    {
        printf("usage: ./align_and_measure_momentum linked_tracks.root title cut Xcenter Ycenter binwidth robustfactor\n");
        return 1;
    }
    TString filename_linked_tracks = argv[1];
    TString title = argv[2];
    TString cut = argv[3];
    double Xcenter, Ycenter, binWidth, robustFactor;
    sscanf(argv[4],"%lf",&Xcenter);
    sscanf(argv[5],"%lf",&Ycenter);
    sscanf(argv[6],"%lf",&binWidth);
    sscanf(argv[7],"%lf",&robustFactor);
    EdbDataProc *dproc = new EdbDataProc;
    EdbPVRec *pvr = new EdbPVRec;
    dproc->ReadTracksTree(*pvr, filename_linked_tracks, cut);
    TObjArray *tracks = pvr->GetTracks();
    int ntrk = tracks->GetEntriesFast();
    if (0 == ntrk)
    {
        return 1;
    }
    FnuDivideAlign align;
    align.SetRobustFactor(robustFactor);
    align.SetBinWidth(binWidth);
    // align.Align(tracks, Xcenter, Ycenter, pvr->Npatterns());

    // FnuMomCoord mom;
    // mom.ReadParFile("~/LEPP/FASERnu/Tools/FnuMomCoord/par/Data_up_to_200plates.txt");
    // TH1D hmom("hmom", "momentum", 100, 0, 7000);
    // for (int itrk = 0; itrk < ntrk; itrk++)
    // {
    //     EdbTrackP *t = (EdbTrackP *)tracks->At(itrk);
    //     double p = mom.CalcMomentum(t);
    //     hmom.Fill(p);
    // }
    // TCanvas c;
    // hmom.Draw();
    // c.Print("momentum_no_align.pdf");
    FnuQualityCheck qc(pvr,"akirec");
    qc.CalcEfficiency();
    // qc.CalcDeltaXY(Xcenter,Ycenter,binWidth);
    // qc.FitDeltaXY();
    // qc.PlotPosRes("posres_no_align");
}