#include <EdbDataSet.h>

class FnuPositionDistribution
{
private:
    EdbPVRec *pvr;
TString title;
    int ntrk;
    int nPID;
    double XYrange;
    double angleCut;
    int plMin;
    int plMax;
    TH2D *positionHist;

public:
    FnuPositionDistribution(EdbPVRec *pvr, TString title);
    ~FnuPositionDistribution();

    TH2D *MakePositionHist();
    void PrintPositionHist(TString filename);
    void WritePositionHist(TString filename);
};