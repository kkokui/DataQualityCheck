#include <EdbDataSet.h>

class FnuResolution
{
    private:
    EdbPVRec *pvr;
    TTree *meanDeltaXY;
    TString title;
    int ntrk;
    int nPID;
    double XYrange;
    double angleCut;
    int plMin;
    int plMax;
    TGraph *meanXGraph, *meanYGraph, *sigmaXGraph, *sigmaYGraph;
    TGraph *meanTXGraph, *meanTYGraph, *sigmaTXGraph, *sigmaTYGraph;
    TH1D *sigmaXHist, *sigmaYHist;
    TH1D *sigmaTXHist, *sigmaTYHist;

    int plate, nseg;
    std::vector<double> *deltaXV, *deltaYV, *deltaTXV, *deltaTYV, *xV, *yV, *slopeXV, *slopeYV;
    std::vector<int> *crossTheLineV, *tridV, *nsegV;
    double sigmaX, sigmaY, meanX, meanY;
    double sigmaTX, sigmaTY, meanTX, meanTY;
    int entries;
    TH1D *hdeltaX,*hdeltaY;
    TH1D *hdeltaTX,*hdeltaTY;
    TH2D *deltaX2DHist;
    TH2D *deltaY2DHist;
    TH2D *entries2DHist;
    public:
    FnuResolution(EdbPVRec *pvr, TString title);
    ~FnuResolution();
    TTree *CalcDeltaXY(double Xcenter, double Ycenter, double bin_width);
    TTree *MakeHistDeltaXY(TTree *deltaXY);
    TTree *FitHistDeltaXY(TTree *treeHistDeltaXY);
    void CalcLSM(double x[], double y[], int N, double &a0, double &a1);
    void MakeGraphHistPositionResolution(TTree *positionResolutionPar);
    void WriteGraphHistPositionResolution(TString filename);
    void PrintGraphHistPositionResolution(TString filename);
    void PrintHistDeltaXY(TTree *treeHistDeltaXY, TString filename);
    void CalcMeanDeltaXY(TTree *deltaXY, double Xcenter, double Ycenter, double cellSize);
    int PrintMeanDeltaXYArrowPlot(TString filename);
    void WriteMeanDeltaXY(TString filename);
    TGraph *GetSigmaXGraph() const;
    TGraph *GetSigmaYGraph() const;

    // angular resolution
    TTree *MakeHistDeltaTXY(TTree *deltaXY);
    TTree *FitHistDeltaTXY(TTree *treeHistDetlaTXY);
    void PrintHistDeltaTXY(TTree *treeHistDetlaTXY, TString filename);
    void MakeGraphHistAngleResolution(TTree *angleResolutionPar);
    void PrintGraphHistAngleResolution(TString filename);
};