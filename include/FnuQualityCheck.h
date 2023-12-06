#include <EdbDataSet.h>
#include <TEfficiency.h>

class FnuQualityCheck
{
private:
    EdbPVRec *pvr;
    TFile *file;
    TTree *deltaXY;
    TTree *posResPar;
    TTree *htree;
    TTree *meanDeltaXY;
    TTree *effInfo;
    TString title;
    int ntrk;
    int nPID;
    double XYrange;
    double angleCut;
    int plMin;
    int plMax;
    TGraph *meanXGraph, *meanYGraph, *sigmaXGraph, *sigmaYGraph;
    TH1D *sigmaXHist, *sigmaYHist;
    std::vector<double> bins_vec_angle;
    std::vector<double> bins_vec_TXTY;
    TEfficiency *eachAngleEfficiency, *eachPlateEfficiency, *eachTXEfficiency, *eachTYEfficiency;
    TH2D *positionHist;
    TH2D *angleHistWide;
    TH2D *angleHistNarrow;
    TH1I *nsegHist;
    TH1I *nplHist;
    TH1I *firstPlateHist;
    TH1I *lastPlateHist;

    // variables for TTree
    int plate;
    std::vector<double> *deltaXV, *deltaYV, *deltaTXV, *deltaTYV, *xV, *yV, *slopeXV, *slopeYV;
    std::vector<int> *crossTheLineV, *tridV, *nsegV;
    double sigmaX, sigmaY, meanX, meanY;
    int entries;
    TH1D *hdeltaX;
    TH1D *hdeltaY;
    int trackID, nseg, W, hitsOnThePlate;
    double x, y, angle, TX, TY;
    TH2D *deltaX2DHist;
    TH2D *deltaY2DHist;
    TH2D *entries2DHist;
    double *secondDifferenceX, *secondDifferenceY, *slopeX, *slopeY;
    int entriesParPlate;
    TH1D *secondDifferenceHist;

public:
    FnuQualityCheck(EdbPVRec *pvr, TString title);
    ~FnuQualityCheck();
    // methods for position resolution
    void CalcDeltaXY(double Xcenter, double Ycenter, double bin_width);
    void FitDeltaXY();
    void FitDeltaXYAllPlatesTogether();
    void CalcLSM(double x[], double y[], int N, double &a0, double &a1);
    void MakePosResGraphHist();
    void WritePosResGraphHist(TString filename);
    void PrintPosResGraphHist(TString filename);
    void WritePosResPar(TString filename);
    void WriteDeltaXY(TString filename);
    void PrintDeltaXYHist(TString filename);
    void CalcMeanDeltaXY(double Xcenter, double Ycenter, double cellSize);
    int PrintMeanDeltaXYArrowPlot(TString filename);
    void WriteMeanDeltaXY(TString filename);
    // methods for efficiency
    void CalcEfficiency();
    void SetBinsAngle(int nbins, double bins[]);
    void SetBinsTXTY(int nbins, double bins[]);
    void PrintEfficiency(TString filename);
    void WriteEfficiencyTree(TString filename);
    void WriteEfficiency(TString filename);
    // methods for position distribution
    void MakePositionHist();
    void PrintPositionHist(TString filename);
    void WritePositionHist(TString filename);
    // methods for angle distribution
    void MakeAngleHist();
    void PrintAngleHist(TString filename);
    void WriteAngleHist(TString filename);
    // methods for nseg
    void MakeNsegHist();
    void PrintNsegHist(TString filename);
    void WriteNsegHist(TString filename);
    // methods for npl
    void MakeNplHist();
    void PrintNplHist(TString filename);
    void WriteNplHist(TString filename);
    // methods for start and end plate
    void MakeFirstLastPlateHist();
    void PrintFirstLastPlateHist(TString filename);
    void WriteFirstLastPlateHist(TString filename);
    // methods for second difference
    TTree *CalcSecondDifference(int cellLength);
    void MakeSecondDifferenceHist(TTree *secondDifferenceTree, int cellLength);
    // methods for summary plot
    void Summarize();
};
