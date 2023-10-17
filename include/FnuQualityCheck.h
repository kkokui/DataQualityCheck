#include <EdbDataSet.h>

class FnuQualityCheck
{
    private:
        TFile *file;
        TTree *deltaXY;
        TTree *posResPar;
        TTree *htree;
        TString title;
        int nPID;
        double XYrange;
        int plMin;
        int plMax;
        double deltaX, deltaY, deltaTX, deltaTY, x_t, y_t, slopeX, slopeY;
        int plate, cross_the_line, trid, nseg;
        std::vector<double> deltaXV, deltaYV, deltaTXV, deltaTYV, xV, yV, slopeXV, slopeYV;
        std::vector<int> crossTheLineV, tridV, nsegV;
        TH1D *hdeltaX;
        TH1D *hdeltaY;

    public:
        FnuQualityCheck(TString ititle);
        ~FnuQualityCheck();
        // void CalcDeltaXYFromRootFile(TString fname = "linked_tracks.root", double Xcenter, double Ycenter, TCut cut = "nseg>=5", double bin_width);
        void CalcDeltaXY(EdbPVRec *pvr, int ntrk, double Xcenter, double Ycenter, double bin_width);
        void FitDeltaXY();
        void lsm(double x[],double y[], int N, double &a0, double &a1);
        void PlotPosRes();
        void WritePosResPar();
        void WriteDeltaXY();
        void PrintDeltaXYHist();
};
