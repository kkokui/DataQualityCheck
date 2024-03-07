#include "FnuAngleDistribution.h"

#include <TCanvas.h>

FnuAngleDistribution::FnuAngleDistribution(EdbPVRec *pvr, TString title)
	: pvr(pvr),
	  title(title),
	  ntrk(pvr->Ntracks())
{
}

FnuAngleDistribution::~FnuAngleDistribution()
{
}

void FnuAngleDistribution::CalcLSM(double x[], double y[], int N, double &a0, double &a1)
{
	// y = a0 + a1*x
	int i;
	double A00 = 0, A01 = 0, A02 = 0, A11 = 0, A12 = 0;

	for (i = 0; i < N; i++)
	{
		A00 += 1.0;
		A01 += x[i];
		A02 += y[i];
		A11 += x[i] * x[i];
		A12 += x[i] * y[i];
	}
	a0 = (A02 * A11 - A01 * A12) / (A00 * A11 - A01 * A01);
	a1 = (A00 * A12 - A01 * A02) / (A00 * A11 - A01 * A01);
}

void FnuAngleDistribution::CalcAngle()
{
    // calculate track angles and fill them to vector
	angleXVec.reserve(ntrk);
	angleYVec.reserve(ntrk);
    TH2D h2ForAngleCenterCalculation("h2ForAngleCenterCalculation", "angle;tan#theta_{x};tan#theta_{y}", 100, -0.02, 0.02, 100, -0.02, 0.02);
	// loop over the tracks
	for (int itrk = 0; itrk < ntrk; itrk++)
	{
		EdbTrackP *t = pvr->GetTrack(itrk);
		int nseg = t->N();
		// cut is nseg>=5
		if (nseg < 5)
			continue;
		// loop over the segments
		double X[nseg];
		double Y[nseg];
		double Z[nseg];
		for (int iseg = 0; iseg < nseg; iseg++)
		{
			X[iseg] = t->GetSegment(iseg)->X();
			Y[iseg] = t->GetSegment(iseg)->Y();
			Z[iseg] = t->GetSegment(iseg)->Z();
		}
		double TX, TY, a0;
		CalcLSM(Z, X, nseg, a0, TX);
		CalcLSM(Z, Y, nseg, a0, TY);
		h2ForAngleCenterCalculation.Fill(TX, TY);
		angleXVec.push_back(TX);
		angleYVec.push_back(TY);
	}
	int locmax, locmay, locmaz;
	h2ForAngleCenterCalculation.GetMaximumBin(locmax, locmay, locmaz);
	angleCenterX = h2ForAngleCenterCalculation.GetXaxis()->GetBinCenter(locmax);
	angleCenterY = h2ForAngleCenterCalculation.GetYaxis()->GetBinCenter(locmay);
}

TH2D *FnuAngleDistribution::MakeAngleHistWide()
{
    // Make 2d histogram of angle in wide range
	CalcAngle();
	double halfRangeWide = 0.3;
	double minXAxisWide = angleCenterX - halfRangeWide;
	double maxXAxisWide = angleCenterX + halfRangeWide;
	double minYAxisWide = angleCenterY - halfRangeWide;
	double maxYAxisWide = angleCenterY + halfRangeWide;
	angleHistWide = new TH2D("angleHistWide", "angle distribution wide (" + title + ");tan#theta_{x};tan#theta_{y};Ntracks", 200, minXAxisWide, maxXAxisWide, 200, minYAxisWide, maxYAxisWide);
	for (int itrk = 0; itrk < angleXVec.size(); itrk++)
	{
		angleHistWide->Fill(angleXVec.at(itrk), angleYVec.at(itrk));
	}
    return angleHistWide;
}

TH2D *FnuAngleDistribution::MakeAngleHistWide(std::vector<double> angleXVec,std::vector<double> angleYVec)
{
    // Make 2d histogram of angle in wide range from argment vector
	const double halfRangeWide = 0.3;
	double minXAxisWide = angleCenterX - halfRangeWide;
	double maxXAxisWide = angleCenterX + halfRangeWide;
	double minYAxisWide = angleCenterY - halfRangeWide;
	double maxYAxisWide = angleCenterY + halfRangeWide;
	angleHistWide = new TH2D("angleHistWide", "angle distribution wide (" + title + ");tan#theta_{x};tan#theta_{y};Ntracks", 200, minXAxisWide, maxXAxisWide, 200, minYAxisWide, maxYAxisWide);
	for (int itrk = 0; itrk < angleXVec.size(); itrk++)
	{
		angleHistWide->Fill(angleXVec.at(itrk), angleYVec.at(itrk));
	}
    return angleHistWide;
}

TH2D *FnuAngleDistribution::MakeAngleHistNarrow()
{
    // Make 2d histogram of angle in narrow range
	CalcAngle();
	const double halfRangeNarrow = 0.01;
	double minXAxisNarrow = angleCenterX - halfRangeNarrow;
	double maxXAxisNarrow = angleCenterX + halfRangeNarrow;
	double minYAxisNarrow = angleCenterY - halfRangeNarrow;
	double maxYAxisNarrow = angleCenterY + halfRangeNarrow;
	angleHistNarrow = new TH2D("angleHistNarrow", "angle distribution narrow (" + title + ");tan#theta_{x};tan#theta_{y};Ntracks", 200, minXAxisNarrow, maxXAxisNarrow, 200, minYAxisNarrow, maxYAxisNarrow);
	for (int itrk = 0; itrk < angleXVec.size(); itrk++)
	{
		angleHistNarrow->Fill(angleXVec.at(itrk), angleYVec.at(itrk));
	}
    return angleHistNarrow;
}

TH2D *FnuAngleDistribution::MakeAngleHistNarrow(std::vector<double> angleXVec,std::vector<double> angleYVec)
{
    // Make 2d histogram of angle in narrow range from argment vector
	const double halfRangeNarrow = 0.01;
	double minXAxisNarrow = angleCenterX - halfRangeNarrow;
	double maxXAxisNarrow = angleCenterX + halfRangeNarrow;
	double minYAxisNarrow = angleCenterY - halfRangeNarrow;
	double maxYAxisNarrow = angleCenterY + halfRangeNarrow;
	angleHistNarrow = new TH2D("angleHistNarrow", "angle distribution narrow (" + title + ");tan#theta_{x};tan#theta_{y};Ntracks", 200, minXAxisNarrow, maxXAxisNarrow, 200, minYAxisNarrow, maxYAxisNarrow);
	for (int itrk = 0; itrk < angleXVec.size(); itrk++)
	{
		angleHistNarrow->Fill(angleXVec.at(itrk), angleYVec.at(itrk));
	}
    return angleHistNarrow;
}

void FnuAngleDistribution::PrintAngleHist(TString filename)
{
	TCanvas ctemp;
	ctemp.Print(filename + "[");
	ctemp.SetRightMargin(0.15);
	angleHistNarrow->Draw("colz");
	angleHistNarrow->SetStats(0);
	angleHistNarrow->GetYaxis()->SetTitleOffset(1.5);
	ctemp.Print(filename);
	ctemp.SetLogz();
	angleHistWide->Draw("colz");
	angleHistWide->SetStats(0);
	angleHistWide->GetYaxis()->SetTitleOffset(1.5);
	ctemp.Print(filename);
	ctemp.Print(filename + "]");
}
void FnuAngleDistribution::WriteAngleHist(TString filename)
{
	TFile fout(filename, "recreate");
	angleHistNarrow->Write();
	angleHistWide->Write();
	fout.Close();
}
