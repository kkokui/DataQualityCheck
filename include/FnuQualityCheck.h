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
    TTree *effInfo;
    TString title;
    int ntrk;
    int nPID;
    double XYrange;
    int plMin;
    int plMax;
    TGraph *meanXGraph, *meanYGraph, *sigmaXGraph, *sigmaYGraph;
    TH1D *sigmaXHist, *sigmaYHist;
    std::vector<double> bins_vec_angle;
    std::vector<double> bins_vec_TXTY;
    TEfficiency *pEff_angle, *pEff_plate, *pEff_TX, *pEff_TY;
    TH2D *positionHist;
    TH2D *angleHistWide;
    TH2D *angleHistNarrow;
    TH1D *nsegHist;
    TH1D *nplHist;

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

public:
    FnuQualityCheck(EdbPVRec *pvr, TString title);
    ~FnuQualityCheck();
    // methods for position resolution
    void CalcDeltaXY(double Xcenter, double Ycenter, double bin_width);
    void FitDeltaXY();
    void CalcLSM(double x[], double y[], int N, double &a0, double &a1);
    void MakePosResGraphHist();
    void WritePosResGraphHist(TString filename);
    void PrintPosResGraphHist(TString filename);
    void WritePosResPar(TString filename);
    void WriteDeltaXY(TString filename);
    void PrintDeltaXYHist(TString filename);
    // methods for efficiency
    void CalcEfficiency();
    void SetBinsAngle(int nbins, double bins[]);
    void SetBinsTXTY(int nbins, double bins[]);
    void PlotEfficiency(TString filename);
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
    // methods for nseg and npl
    void MakeNsegHist();
    void PrintNsegHist(TString filename);
    void WriteNsegHist(TString filename);
    void MakeNplHist();
    void PrintNplHist(TString filename);
    void WriteNplHist(TString filename);

    // methods for summary plot
    void PrintSummaryPlot();
};
