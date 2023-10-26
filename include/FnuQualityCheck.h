#include <EdbDataSet.h>

class FnuQualityCheck
{
    private:
        EdbPVRec *pvr;
        TFile *file;
        TTree *deltaXY;
        TTree *posResPar;
        TTree *htree;
        TTree *effInfo;
        TString title;
        int ntrk;
        int nPID;
        double XYrange;
        int plMin;
        int plMax;
        // variables for TTree
        int plate;
        std::vector<double> *deltaXV=0, *deltaYV=0, *deltaTXV=0, *deltaTYV=0, *xV=0, *yV=0, *slopeXV=0, *slopeYV=0;
        std::vector<int> *crossTheLineV=0, *tridV=0, *nsegV=0;
        TH1D *hdeltaX;
        TH1D *hdeltaY;
        int trackID, nseg, W, hitsOnThePlate;
        double x, y, angle, TX, TY;

    public:
        FnuQualityCheck(EdbPVRec *pvr, TString title);
        ~FnuQualityCheck();
        // methods for position resolution
        // void CalcDeltaXYFromRootFile(TString fname = "linked_tracks.root", double Xcenter, double Ycenter, TCut cut = "nseg>=5", double bin_width);
        void CalcDeltaXY(double Xcenter, double Ycenter, double bin_width);
        void FitDeltaXY();
        void CalcLSM(double x[],double y[], int N, double &a0, double &a1);
        void PlotPosRes(TString filename);
        void WritePosResPar(TString filename);
        void WriteDeltaXY(TString filename);
        void PrintDeltaXYHist(TString filename);
        // methods for efficiency
        void CalcEfficiency();
        void PrintEfficiency();
        void WriteEfficiencyTree(TString filename);
};
