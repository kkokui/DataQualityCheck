#include <EdbDataSet.h>
#include <TEfficiency.h>

class FnuQualityCheck
{
private:
    EdbPVRec *pvr;
    TFile *file;
    TTree *deltaXY;
    TTree *posResPar;
    TTree *angResPar;
    TTree *htree;
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
    TH2D *positionHist;
    TH2D *angleHistWide;
    TH2D *angleHistNarrow;
    TH1I *PHHistNsegMin2;
    TH1I *PHHistNsegMin4;
    TH1D *PHHistMean;
    TEfficiency *efficiencyForEachAngle;
    TEfficiency *efficiencyForEachPlate;
    TH1I *nsegHist;
    TH1I *nplHist;
    TH1I *firstPlateHist;
    TH1I *lastPlateHist;

    // variables for TTree
    int plate;
    std::vector<double> *deltaXV, *deltaYV, *deltaTXV, *deltaTYV, *xV, *yV, *slopeXV, *slopeYV;
    std::vector<int> *crossTheLineV, *tridV, *nsegV;
    double sigmaX, sigmaY, meanX, meanY;
    double sigmaTX, sigmaTY, meanTX, meanTY;
    int entries;
    TH1D *hdeltaX, *hdeltaY;
    TH1D *hdeltaTX, *hdeltaTY;
    TH2D *deltaX2DHist;
    TH2D *deltaY2DHist;
    TH2D *entries2DHist;
    double *secondDifferenceX, *secondDifferenceY, *slopeX, *slopeY;
    int entriesParPlate;
    TH1D *secondDifferenceHist;
    TGraphErrors *grErrorsSecondDifferenceMeanX;
    TGraphErrors *grErrorsSecondDifferenceMeanY;

public:
    FnuQualityCheck(EdbPVRec *pvr, TString title);
    ~FnuQualityCheck();

    // methods for PH
    TH1I *MakePHHist(int nsegMin);
    TH1D *MakePHMeanHist(int nsegMin);

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
    void Summarize(double Xcenter, double Ycenter, double bin_width);
    void WriteSummaryPlot(TString filename);
};
