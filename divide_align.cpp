#include "FnuDivideAlign.h"
#include <EdbDataSet.h>
double rangeXY = 8500;
void ApplyAlignBicubic(EdbSegP *s, double Xcenter, double Ycenter, double binWidth)
{
	// Apply alignment to segments which passed a divided area

	int pid = s->PID();

	// get shift value calculated by divide align
	// double

	// get x and y values of the segment
	double segmentPositionX = s->X();
	double segmentPositionY = s->Y();
	// detect nearest middle position of division
	int nearestMiddlePointNumberX = (segmentPositionX - (Xcenter - rangeXY) - binWidth / 2) / binWidth;
	int nearestMiddlePointNumberY = (segmentPositionY - (Ycenter - rangeXY) - binWidth / 2) / binWidth;
	double nearestMiddlePointPositionX = Xcenter - rangeXY + nearestMiddlePointNumberX * binWidth;
	double nearestMiddlePointPositionY = Ycenter - rangeXY + nearestMiddlePointNumberY * binWidth;
	// detect 16 reference values
	double referencePositionX[4][4];
	double referencePositionY[4][4];

	double referenceShiftX[4][4];
	double referenceShiftY[4][4];
	// if(segmentPositionX<nearestBinPositionX)
	// {

	// 	referenceShiftX[0][0] =
	// }
	// calculate alignment parameter
	// apply alignment
}
double GetWeight(double distance, double sharpFactor)
{
	double distanceAbs = fabs(distance);
	if (distanceAbs <= 1)
	{
		return (sharpFactor + 2) * distanceAbs * distanceAbs * distanceAbs - (sharpFactor + 3) * distanceAbs * distanceAbs + 1;
	}
	else if (distanceAbs <= 2)
	{
		return sharpFactor * distanceAbs * distanceAbs * distanceAbs - 5 * sharpFactor * distanceAbs * distanceAbs + 8 * sharpFactor * distanceAbs - 4 * sharpFactor;
	}
	else
	{
		return 0;
	}
	
}
int main(int argc, char *argv[])
{
	if (argc < 6)
	{
		printf("Usage: ./calc_dxy linked_tracks.root title reco binWidth robustFactor\n");
		return 1;
	}

	TString filename_linked_tracks = argv[1];
	TString title = argv[2];
	double binWidth;
	float robustFactor;
	int reco;
	sscanf(argv[3], "%d", &reco);
	double Xcenter = (reco - 1) % 9 * 15000 + 5000;
	double Ycenter = (reco - 1) / 9 * 15000 + 5000;
	sscanf(argv[4], "%lf", &binWidth);
	sscanf(argv[5], "%f", &robustFactor);

	EdbDataProc *dproc = new EdbDataProc;
	EdbPVRec *pvr = new EdbPVRec;
	dproc->ReadTracksTree(*pvr, filename_linked_tracks, "1");
	// dproc->ReadTracksTree(*pvr, filename_linked_tracks, "Entry$<5000");
	// dproc->ReadTracksTree(*pvr, filename_linked_tracks, "Entry$%100==0");

	TObjArray *tracks = pvr->GetTracks();
	int nPatterns = pvr->Npatterns();
	int ntrk = tracks->GetEntriesFast();

	if (ntrk == 0)
	{
		printf("ntrk==0\n");
		return 0;
	}

	FnuDivideAlign align;
	align.SetRobustFactor(robustFactor);
	align.SetBinWidth(binWidth);
	// align.Align(tracks, Xcenter, Ycenter, nPatterns);
	// align.WriteAlignPar("align_output/alignPar_" + title + ".root");

	TFile *alignParFile = new TFile(Form("align_output/alignPar_binWidth%.0f_robustFactor%.1f.root", binWidth, robustFactor));
	if (alignParFile->IsZombie())
	{
		printf("Failed to open alignment parameter file\n");
		return 1;
	}
	TTree *alignPar = (TTree *)gDirectory->Get("alignPar");
	std::vector<std::vector<std::vector<double>>> shiftXEachDivision;
	std::vector<std::vector<std::vector<double>>> shiftYEachDivision;
	shiftXEachDivision.resize(nPatterns);
	shiftYEachDivision.resize(nPatterns);
	int nXID = (alignPar->GetMaximum("iX") - alignPar->GetMinimum("iX"))/binWidth +1;
	int nYID = (alignPar->GetMaximum("iY") - alignPar->GetMinimum("iY"))/binWidth +1;
	double miniX = alignPar->GetMinimum("iX");
	double miniY = alignPar->GetMinimum("iY");
	int nDivisionX = rangeXY * 2 / binWidth + 1;
	int nDivisionY = rangeXY * 2 / binWidth + 1;
	for (int pid = 0; pid < nPatterns; pid++)
	{
		shiftXEachDivision[pid].resize(nXID);
		shiftYEachDivision[pid].resize(nXID);
		for (int iXID = 0; iXID < nXID; iXID++)
		{
			shiftXEachDivision[pid][iXID].resize(nYID);
			shiftYEachDivision[pid][iXID].resize(nYID);
		}
	}
	double iX, iY, shiftX, shiftY;
	int pid;
	alignPar->SetBranchAddress("iX", &iX);
	alignPar->SetBranchAddress("iY", &iY);
	alignPar->SetBranchAddress("shiftX", &shiftX);
	alignPar->SetBranchAddress("shiftY", &shiftY);
	alignPar->SetBranchAddress("pid", &pid);
	for (int ient = 0; ient < alignPar->GetEntries(); ient++)
	{
		alignPar->GetEntry(ient);
		int iXID = (iX - miniX) / binWidth;
		int iYID = (iY - miniY) / binWidth;
		shiftXEachDivision[pid][iXID][iYID] = shiftX;
		shiftYEachDivision[pid][iXID][iYID] = shiftY;
	}
	TFile *testShiftAfterFile = new TFile("testShiftAfter_" + title + ".root", "recreate");
	if (testShiftAfterFile->IsZombie())
	{
		printf("Failed to open the file for the test of bicubic\n");
		return 1;
	}
	TTree *testShiftAfter = new TTree("testShiftAfter", "test of shiftX after bicubic interpolation");
	double segmentPositionXForTree, segmentPositionYForTree, shiftXForThisSegForTree, shiftYForThisSegForTree, shiftXNearest, shiftYNearest;
	int pidForTree;
	testShiftAfter->Branch("X", &segmentPositionXForTree);
	testShiftAfter->Branch("Y", &segmentPositionYForTree);
	testShiftAfter->Branch("shiftX", &shiftXForThisSegForTree);
	testShiftAfter->Branch("shiftY", &shiftYForThisSegForTree);
	testShiftAfter->Branch("shiftXNearest", &shiftXNearest);
	testShiftAfter->Branch("shiftYNearest", &shiftYNearest);
	testShiftAfter->Branch("pid", &pidForTree);


	TTree *testWeight = new TTree("testWeight","test");
	double weightXForTree,weightYForTree,distanceXForTree,distanceYForTree;
	testWeight->Branch("weightX",&weightXForTree);
	testWeight->Branch("weightY",&weightYForTree);
	testWeight->Branch("distanceX",&distanceXForTree);
	testWeight->Branch("distanceY",&distanceYForTree);
	for (int itrk = 0; itrk < ntrk; itrk++)
	{
		EdbTrackP *t = (EdbTrackP *)tracks->At(itrk);
		int nseg = t->N();
		for (int iseg = 0; iseg < nseg; iseg++)
		{
			EdbSegP *s = t->GetSegment(iseg);
			// Apply alignment to segments which passed a divided area

			int pid = s->PID();

			// get x and y values of the segment
			double segmentPositionX = s->X();
			double segmentPositionY = s->Y();
			// get the segment position in a grid of width 1
			double segmentPositionGridX = (segmentPositionX - miniX) / binWidth;
			double segmentPositionGridY = (segmentPositionY - miniY) / binWidth;
			int referencePositionGridX[4];
			for (int i = 0; i < 4; i++)
			{
				if (segmentPositionGridX <= 1)
				{
					referencePositionGridX[i] = i;
				}
				else if (segmentPositionGridX < nXID - 2)
				{
					referencePositionGridX[i] = (int)segmentPositionGridX - 1 + i;
				}
				else
				{
					referencePositionGridX[i] = nXID - 4 + i;
				}
			}
			int referencePositionGridY[4];
			for (int i = 0; i < 4; i++)
			{
				if (segmentPositionGridY <= 1)
				{
					referencePositionGridY[i] = i;
				}
				else if (segmentPositionGridY < nYID - 2)
				{
					referencePositionGridY[i] = (int)segmentPositionGridY - 1 + i;
				}
				else
				{
					referencePositionGridY[i] = nYID - 4 + i;
				}
			}
			// detect nearest middle position of division
			// int nearestMiddlePointNumberX = (segmentPositionX - (Xcenter - rangeXY) - binWidth / 2) / binWidth;
			// int nearestMiddlePointNumberY = (segmentPositionY - (Ycenter - rangeXY) - binWidth / 2) / binWidth;
			// double nearestMiddlePointPositionX = Xcenter - rangeXY + nearestMiddlePointNumberX * binWidth;
			// double nearestMiddlePointPositionY = Ycenter - rangeXY + nearestMiddlePointNumberY * binWidth;

			// calculate 16 distances of references
			double referenceDistanceX[4];
			double referenceDistanceY[4];
			for (int i = 0; i < 4; i++)
			{
				referenceDistanceX[i] = fabs(segmentPositionGridX - referencePositionGridX[i]);
				referenceDistanceY[i] = fabs(segmentPositionGridY - referencePositionGridY[i]);
			}

			double weightX[4];
			double weightY[4];
			double sharpFactor = -0.5;
			for (int i = 0; i < 4; i++)
			{
				weightX[i] = GetWeight(referenceDistanceX[i], sharpFactor);
				weightY[i] = GetWeight(referenceDistanceY[i], sharpFactor);
				weightXForTree=weightX[i];
				weightYForTree=weightY[i];
				distanceXForTree = referenceDistanceX[i];
				distanceYForTree = referenceDistanceY[i];
				testWeight->Fill();
			}
			double shiftXForThisSeg = 0;
			double shiftYForThisSeg = 0;
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					shiftXForThisSeg += shiftXEachDivision[pid][referencePositionGridX[i]][referencePositionGridY[j]] * weightX[i] * weightY[j];
					shiftYForThisSeg += shiftYEachDivision[pid][referencePositionGridX[i]][referencePositionGridY[j]] * weightX[i] * weightY[j];
				}
			}
			segmentPositionXForTree = segmentPositionX;
			segmentPositionYForTree = segmentPositionY;
			shiftXForThisSegForTree = shiftXForThisSeg;
			shiftYForThisSegForTree = shiftYForThisSeg;
			int nearestXID = (segmentPositionX - (miniX-binWidth/2)) / binWidth;
			int nearestYID = (segmentPositionY - (miniY-binWidth/2)) / binWidth;
			shiftXNearest = shiftXEachDivision[pid][nearestXID][nearestYID];
			shiftYNearest = shiftYEachDivision[pid][nearestXID][nearestYID];
			pidForTree = pid;
			testShiftAfter->Fill();

			// apply alignment
			s->SetX(s->X() + shiftXForThisSeg);
			s->SetY(s->Y() + shiftYForThisSeg);
		}
	}
	testShiftAfter->Write();
	testWeight->Write();

	TObjArray *tracks_t = new TObjArray;
	for (int itrk = 0; itrk < ntrk; itrk++)
	{
		EdbTrackP *t = (EdbTrackP *)tracks->At(itrk);
		tracks_t->Add(t);
	}
	dproc->MakeTracksTree(*tracks_t, 0, 0, Form("/data/Users/kokui/FASERnu/F222/zone4/temp/TFD/vert32063_pl053_167_new/reco32_065000_050000/v15/linked_tracks_after_align_binWidth%.0f_robustFactor%.1f_bicubic.root", binWidth, robustFactor));

	return 0;
}
