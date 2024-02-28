#include "FnuPositionDistribution.h"

#include <TCanvas.h>

FnuPositionDistribution::FnuPositionDistribution(EdbPVRec *pvr, TString title)
	: pvr(pvr),
	  title(title),
	  nPID(pvr->Npatterns()),
	  plMin(pvr->GetPattern(0)->Plate()),
	  plMax(pvr->GetPattern(nPID - 1)->Plate()),
	  ntrk(pvr->Ntracks()),
	  angleCut(0.01),
	  XYrange(8500)
{
}

FnuPositionDistribution::~FnuPositionDistribution()
{
}

TH2D *FnuPositionDistribution::MakePositionHist()
{
	std::vector<double> positionXVec;
	std::vector<double> positionYVec;
	positionXVec.reserve(ntrk);
	positionYVec.reserve(ntrk);
	for (int itrk = 0; itrk < ntrk; itrk++)
	{
		EdbTrackP *t = pvr->GetTrack(itrk);
		if (t->N() < 5)
			continue;
		double X = t->X();
		double Y = t->Y();
		positionXVec.push_back(X);
		positionYVec.push_back(Y);
	}
	auto maxminXIterator = std::minmax_element(positionXVec.begin(), positionXVec.end());
	auto maxminYIterator = std::minmax_element(positionYVec.begin(), positionYVec.end());
	double minX = *maxminXIterator.first;
	double maxX = *maxminXIterator.second;
	double minY = *maxminYIterator.first;
	double maxY = *maxminYIterator.second;
	double marginX = (maxX - minX) / 10;
	double marginY = (maxY - minY) / 10;
	double minXAxis = minX - marginX;
	double maxXAxis = maxX + marginX;
	double minYAxis = minY - marginY;
	double maxYAxis = maxY + marginY;

	positionHist = new TH2D("positionHist", "position distribution (" + title + ");x (#mum);y (#mum);Ntracks / cm^{2}", 100, minXAxis, maxXAxis, 100, minYAxis, maxYAxis);
	double area = (maxXAxis - minXAxis) / 100 / 10000 * (maxYAxis - minYAxis) / 100 / 10000; // in cm^2
	for (int itrk = 0; itrk < positionXVec.size(); itrk++)
	{
		positionHist->Fill(positionXVec.at(itrk), positionYVec.at(itrk), 1 / area);
	}
    return positionHist;
}
void FnuPositionDistribution::PrintPositionHist(TString filename)
{
	TCanvas ctemp;
	ctemp.SetRightMargin(0.15);
	positionHist->Draw("colz");
	positionHist->SetStats(0);
	positionHist->GetYaxis()->SetTitleOffset(1.5);
	ctemp.Print(filename);
}
void FnuPositionDistribution::WritePositionHist(TString filename)
{
	TFile fout(filename, "recreate");
	positionHist->Write();
	fout.Close();
}