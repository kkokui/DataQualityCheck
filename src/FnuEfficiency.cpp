#include "FnuEfficiency.h"

#include <TGraphAsymmErrors.h>
#include <TCanvas.h>

FnuEfficiency::FnuEfficiency(EdbPVRec *pvr, TString title)
	: pvr(pvr),
	  title(title),
	  nPID(pvr->Npatterns()),
	  plMin(pvr->GetPattern(0)->Plate()),
	  plMax(pvr->GetPattern(nPID - 1)->Plate()),
	  ntrk(pvr->Ntracks()),
	  angleCut(0.01),
	  XYrange(8500)
{
	double bins_arr_angle[] = {0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5};
	SetBinsAngle(26, bins_arr_angle);
	double bins_arr_TXTY[] = {-0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.19, -0.18, -0.17, -0.16, -0.15, -0.14, -0.13, -0.12, -0.11, -0.10, -0.09, -0.08, -0.07, -0.06, -0.05, -0.04, -0.03, -0.02, -0.01, 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5};
	SetBinsTXTY(52, bins_arr_TXTY);
}

FnuEfficiency::~FnuEfficiency()
{
}

void FnuEfficiency::SetBinsAngle(int nbins, double bins[])
{
	bins_vec_angle.assign(&bins[0], &bins[nbins + 1]);
}

void FnuEfficiency::SetBinsTXTY(int nbins, double bins[])
{
	bins_vec_TXTY.assign(&bins[0], &bins[nbins + 1]);
}

void FnuEfficiency::CalcEfficiency()
{
	double bins_angle[bins_vec_angle.size()];
	std::copy(bins_vec_angle.begin(), bins_vec_angle.end(), bins_angle);
	int nbins_angle = bins_vec_angle.size() - 1;
	double bins_TXTY[bins_vec_TXTY.size()];
	std::copy(bins_vec_TXTY.begin(), bins_vec_TXTY.end(), bins_TXTY);
	int nbins_TXTY = bins_vec_TXTY.size() - 1;

	efficiencyForEachAngle = new TEfficiency("Eff_angle", Form("Efficiency for each angle (%s);tan#theta;efficiency", title.Data()), nbins_angle, bins_angle);
	efficiencyForEachPlate = new TEfficiency("Eff_plate", Form("Efficiency for each plate (%s);plate;efficiency", title.Data()), plMax - plMin + 1, plMin - 0.5, plMax + 0.5);
	efficiencyForEachTX = new TEfficiency("Eff_TX", Form("Efficiency for each TX (%s);tan#theta;efficiency", title.Data()), nbins_TXTY, bins_TXTY);
	efficiencyForEachTY = new TEfficiency("Eff_TY", Form("Efficiency for each TY (%s);tan#theta;efficiency", title.Data()), nbins_TXTY, bins_TXTY);

	effInfo = new TTree("effInfo", "efficiency infomation");
	effInfo->Branch("trackID", &trackID);
	effInfo->Branch("x", &x);
	effInfo->Branch("y", &y);
	effInfo->Branch("angle", &angle);
	effInfo->Branch("TX", &TX);
	effInfo->Branch("TY", &TY);
	effInfo->Branch("plate", &plate);
	effInfo->Branch("nseg", &nseg);
	effInfo->Branch("W", &W);
	effInfo->Branch("hitsOnThePlate", &hitsOnThePlate);

	double x1, y1, z1, x2, y2, z2;
	for (int itrk = 0; itrk < ntrk; itrk++)
	{
		EdbTrackP *t = pvr->GetTrack(itrk);
		int nseg = t->N();

		for (int iPID = 0; iPID < nPID; iPID++)
		{
			EdbPattern *pattern = pvr->GetPattern(iPID);
			if (pattern == NULL)
				continue;
			int iplate = pattern->Plate();
			int counts = 0;
			hitsOnThePlate = 0;
			W = 0;
			for (int iseg = 0; iseg < nseg; iseg++)
			{
				EdbSegP *s = t->GetSegment(iseg);
				int sPID = s->PID();
				if (iPID == 0) // first plate
				{
					if (sPID > iPID + 4)
						break;
					if (sPID == iPID + 4 || sPID == iPID + 3)
						counts++;
					if (sPID == iPID + 1)
					{
						x1 = s->X();
						y1 = s->Y();
						z1 = s->Z();
						counts++;
					}
					if (sPID == iPID + 2)
					{
						x2 = s->X();
						y2 = s->Y();
						z2 = s->Z();
						counts++;
					}
				}
				if (iPID == 1) // second plate
				{
					if (sPID > iPID + 3)
						break;
					if (sPID == iPID + 3 || sPID == iPID + 2)
						counts++;
					if (sPID == iPID - 1)
					{
						x1 = s->X();
						y1 = s->Y();
						z1 = s->Z();
						counts++;
					}
					if (sPID == iPID + 1)
					{
						x2 = s->X();
						y2 = s->Y();
						z2 = s->Z();
						counts++;
					}
				}
				if (iPID >= 2 && iPID < nPID - 2) // mid plate
				{
					if (sPID > iPID + 2)
						break;
					if (sPID < iPID - 2)
						continue;
					if (sPID == iPID - 2 || sPID == iPID + 2)
						counts++;
					if (sPID == iPID - 1)
					{
						x1 = s->X();
						y1 = s->Y();
						z1 = s->Z();
						counts++;
					}
					if (sPID == iPID + 1)
					{
						x2 = s->X();
						y2 = s->Y();
						z2 = s->Z();
						counts++;
					}
				}
				if (iPID == nPID - 2) // second from last plate
				{
					if (sPID < iPID - 3)
						continue;
					if (sPID == iPID - 3 || sPID == iPID - 2)
						counts++;
					if (sPID == iPID - 1)
					{
						x1 = s->X();
						y1 = s->Y();
						z1 = s->Z();
						counts++;
					}
					if (sPID == iPID + 1)
					{
						x2 = s->X();
						y2 = s->Y();
						z2 = s->Z();
						counts++;
					}
				}
				if (iPID == nPID - 1) // last plate
				{
					if (sPID < iPID - 4)
						continue;
					if (sPID == iPID - 4 || sPID == iPID - 3)
						counts++;
					if (sPID == iPID - 2)
					{
						x1 = s->X();
						y1 = s->Y();
						z1 = s->Z();
						counts++;
					}
					if (sPID == iPID - 1)
					{
						x2 = s->X();
						y2 = s->Y();
						z2 = s->Z();
						counts++;
					}
				}
				if (sPID == iPID)
				{
					hitsOnThePlate = 1;
					W = s->W();
				}
			}
			if (counts == 4)
			{
				TX = (x2 - x1) / (z2 - z1);
				TY = (y2 - y1) / (z2 - z1);
				angle = sqrt(TX * TX + TY * TY);
				efficiencyForEachAngle->Fill(hitsOnThePlate, angle);
				efficiencyForEachPlate->Fill(hitsOnThePlate, iplate);
				efficiencyForEachTX->Fill(hitsOnThePlate, TX);
				efficiencyForEachTY->Fill(hitsOnThePlate, TY);
				trackID = t->ID();
				plate = iplate;
				if (iPID == 0)
				{
					x = 2 * x1 - x2;
					y = 2 * y1 - y2;
				}
				if (iPID >= 1 && iPID < nPID - 1)
				{
					x = (x1 + x2) / 2;
					y = (y1 + y2) / 2;
				}
				if (iPID == nPID - 1)
				{
					x = 2 * x2 - x1;
					y = 2 * y2 - y1;
				}
				effInfo->Fill();
			}
		}
	}
}

void FnuEfficiency::PrintEfficiency(TString filename)
{
	// Plot efficiencies and print them.
	TCanvas *c = new TCanvas();
	c->Print(filename + "[");
	c->SetLogy(1);
	efficiencyForEachAngle->GetCopyTotalHisto()->Draw();
	c->Print(filename);
	c->SetLogy(0);
	efficiencyForEachAngle->Draw();
	c->Print(filename);
	efficiencyForEachPlate->GetCopyTotalHisto()->Draw();
	c->Print(filename);
	efficiencyForEachPlate->Draw();
	c->Print(filename);

	c->SetLogy(1);
	efficiencyForEachTX->GetCopyTotalHisto()->Draw();
	c->Print(filename);
	c->SetLogy(0);
	efficiencyForEachTX->Draw();
	c->Print(filename);

	c->SetLogy(1);
	efficiencyForEachTY->GetCopyTotalHisto()->Draw();
	c->Print(filename);
	c->SetLogy(0);
	efficiencyForEachTY->Draw();
	c->Print(filename);

	c->Print(filename + "]");
}

void FnuEfficiency::WriteEfficiencyTree(TString filename)
{
	TFile fout(filename, "recreate");
	effInfo->Write();
	fout.Close();
}

void FnuEfficiency::WriteEfficiency(TString filename)
{
	TFile fout(filename, "recreate");
	efficiencyForEachAngle->Write();
	efficiencyForEachPlate->Write();
	efficiencyForEachTX->Write();
	efficiencyForEachTY->Write();
	fout.Close();
}

TEfficiency *FnuEfficiency::GetEfficiencyForEachAngle() const
{
    return efficiencyForEachAngle;
}

TEfficiency *FnuEfficiency::GetEfficiencyForEachPlate() const
{
    return efficiencyForEachPlate;
}