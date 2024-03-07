#include <EdbDataSet.h>

class FnuAngleDistribution
{
private:
    EdbPVRec *pvr;
    TString title;
    int ntrk;
    std::vector<double> angleXVec;
	std::vector<double> angleYVec;
    double angleCenterX;
    double angleCenterY;
    TH2D *angleHistWide;
    TH2D *angleHistNarrow;
public:
    FnuAngleDistribution(EdbPVRec *pvr, TString title);
    ~FnuAngleDistribution();
    void CalcLSM(double x[], double y[], int N, double &a0, double &a1);
    void CalcAngle();
    TH2D *MakeAngleHistWide();
    TH2D *MakeAngleHistWide(std::vector<double> angleXVec,std::vector<double> angleYVec);
    TH2D *MakeAngleHistNarrow();
    TH2D *MakeAngleHistNarrow(std::vector<double> angleXVec,std::vector<double> angleYVec);
    void PrintAngleHist(TString filename);
    void WriteAngleHist(TString filename);

};