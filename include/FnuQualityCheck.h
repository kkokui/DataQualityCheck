#include <EdbDataSet.h>

class FnuQualityCheck
{
    private:
        TFile *file;
        TTree *deltaXY;
        TTree *posResPar;
        TString title;
        int nPID;
        double XYrange;
        int plMin;
        int plMax;
    public:
        FnuQualityCheck(TString ititle);
        ~FnuQualityCheck();
        void CalcDeltaXY(EdbPVRec *pvr, int ntrk, double Xcenter, double Ycenter, double bin_width);
        void FitDeltaXY();
        void lsm(double x[],double y[], int N, double &a0, double &a1);
        void PlotPosRes();
        void WritePosResPar();
};
