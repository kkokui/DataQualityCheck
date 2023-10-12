#include <TVirtualFitter.h>
#include <TObjArray.h>
#include <EdbPattern.h>
// #include <EdbDataSet.h>
struct cudaSegment{
	int flag, pid;
	float x,y,z;
};

const int NPIDMAX = 800;

struct cudaTrack{
	float x,y,z,tx,ty,tx_first8,ty_first8;
	int nseg;
	cudaSegment segments[NPIDMAX];
};

class FnuDivideAlign {
    private:
        double binWidth;
        TTree *alignPar;

        TVirtualFitter *gMinuit;
        double p[NPIDMAX*2] ={};
        int nPID;
        double XYrange;

        // Data buffer for the GPU process
        cudaTrack* h_tracks;
        int* indexArray;
    public:
        FnuDivideAlign();
        ~FnuDivideAlign();
        void SetBinWidth(double bwidth);
        void SetRobustFactor(float rfactor);
        double GetBinWidth();
        float GetRobustFactor();
        void calc_align_par(TObjArray *tracks,double iX, double iY, int fixflag);
        int count_passed_seg(EdbTrackP *t, double iX, double iY);
        void apply_align(EdbTrackP *t, double iX, double iY);
        int dedicated_align(TObjArray *tracks,double Xcenter, double Ycenter,int nPatterns);
        void WriteAlignPar(TString filename);
};