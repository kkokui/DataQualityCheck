#include <TVirtualFitter.h>
#include <TObjArray.h>
#include <EdbPattern.h>

const int NPIDMAX = 800;

class FnuDivideAlign {
    private:
        double binWidth;
        TTree *alignPar;
        TVirtualFitter *minuit;
        double p[NPIDMAX*2];
        int nPID;
        double rangeXY;
        // values for TTree
        double iXBranchValue, iYBranchValue, shiftXBranchValue, shiftYBranchValue;
        int pidBranchValue;

    public:
        FnuDivideAlign();
        ~FnuDivideAlign();
        void SetBinWidth(double bwidth);
        void SetRobustFactor(float rfactor);
        double GetBinWidth();
        float GetRobustFactor();
        void CalcAlignPar(TObjArray *tracks,double iX, double iY, int fixflag);
        int CountPassedSeg(EdbTrackP *t, double iX, double iY);
        void ApplyAlign(EdbTrackP *t, double iX, double iY);
        int Align(TObjArray *tracks,double Xcenter, double Ycenter,int nPatterns);
        void WriteAlignPar(TString filename = "alignPar.root");
};