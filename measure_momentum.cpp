#include <FnuMomCoord.hpp>

void CalcAllTrackMomentum(TObjArray *tracks,TString title)
{
    FnuMomCoord mom;
    // printf("Momentum calculation started\n");
    mom.ReadParFile("/home/kokui/LEPP/FASERnu/Tools/FnuMomCoord/par/Data_up_to_100plates_mod1.txt");
    // TCanvas *c = new TCanvas();
    // c->Print("test.pdf[");
    TH1D hmom("hmom", "momentum", 100, 0, 7000);
    int ntrk = tracks->GetEntriesFast();
    for (int itrk = 0; itrk < ntrk; itrk++)
    {
        EdbTrackP *t = (EdbTrackP *)tracks->At(itrk);
        // tracks with angle diff max larger than 1 mrad is ignored.
        // FnuMomCoord::CalcTrackAngleDiffMax() returns mrad.
        // if(mom.CalcTrackAngleDiffMax(t)>=1)
        // {
        //     continue;
        // }
        double p = mom.CalcMomentum(t);
        // mom.DrawMomGraphCoord(t,c,"test");
        hmom.Fill(p);
        if(0==itrk%1000)
        {
            printf("%3d%%\r",itrk*100/ntrk);
            fflush(stdout);
        }
    }
    // c->Print("test.pdf]");
    printf("100%% done\n", ntrk, ntrk);
    mom.WriteRootFile("momentum_output/nt_"+title);
}

int main(int argc, char *argv[])
{
    if (argc < 4)
    {
        printf("usage: ./align_and_measure_momentum linked_tracks.root title cut\n");
        return 1;
    }
    TString filename_linked_tracks = argv[1];
    TString title = argv[2];
    TString cut = argv[3];
    EdbDataProc *dproc = new EdbDataProc;
    EdbPVRec *pvr = new EdbPVRec;
    dproc->ReadTracksTree(*pvr, filename_linked_tracks, cut);
    TObjArray *tracks = pvr->GetTracks();
    int ntrk = tracks->GetEntriesFast();
    if (0 == ntrk)
    {
        return 1;
    }

    CalcAllTrackMomentum(tracks,title);
}