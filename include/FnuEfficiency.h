#include <EdbDataSet.h>
#include <TEfficiency.h>

class FnuEfficiency
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

    TTree *effInfo;

    std::vector<double> bins_vec_angle;
    std::vector<double> bins_vec_TXTY;
    TEfficiency *efficiencyForEachAngle, *efficiencyForEachPlate, *efficiencyForEachTX, *efficiencyForEachTY;

    int plate, trackID, nseg, W, hitsOnThePlate;
    double x, y, angle, TX, TY;
public:
    FnuEfficiency(EdbPVRec *pvr, TString title);
    ~FnuEfficiency();

    void CalcEfficiency();
    void SetBinsAngle(int nbins, double bins[]);
    void SetBinsTXTY(int nbins, double bins[]);
    void PrintEfficiency(TString filename);
    void WriteEfficiencyTree(TString filename);
    void WriteEfficiency(TString filename);
    TEfficiency *GetEfficiencyForEachAngle() const;
    TEfficiency *GetEfficiencyForEachPlate() const;
};