#include <stdio.h>
#include <fstream>
#include <sstream>
#include <EdbDataSet.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TH1.h>

int main(int argc , char *argv[]){
	if(argc<3){
		printf("Usage : ./efficiency linked_tracks.root title\n");
		return 1;
	}
	TString filename_linked_tracks = argv[1];
	TString title = argv[2]; // used for title of histograms and name of output file
	TString cut = "nseg>=5";

	EdbDataProc *dproc = new EdbDataProc;
	EdbPVRec *pvr = new EdbPVRec;
	dproc->ReadTracksTree(*pvr, filename_linked_tracks, cut);

	int plMin = pvr->GetPatternByPID(0)->Plate();
	int plMax = plMin + pvr->Npatterns() - 1;
	
	TEfficiency *pEff_angle =0;
	TEfficiency *pEff_plate =0;
	
	int ntrk = pvr->Ntracks();
	double bins[] = {0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.25,0.3,0.35,0.4,0.45,0.5};
	int nbins = 26;
	
	TH1D *h_angle_total = new TH1D("hist_angle_total", title+";tan#theta;", nbins, bins);
	TH1D *h_angle_passed = new TH1D("hist_angle_passed", title+";tan#theta;", nbins, bins);
	TH1D *h_plate_total = new TH1D("hist_plate_total", title+";plate;", plMax-plMin, plMin,plMax);
	TH1D *h_plate_passed = new TH1D("hist_plate_passed", title + ";plate;", plMax - plMin, plMin, plMax);

	TFile treefout(Form("efficiency_output/tree_not_passed_%s.root", title.Data()), "recreate"); // to check the distribution of passed segments
	TNtuple *tree = new TNtuple("tree","title","trackID:x:y:angle:pl:nseg");
	int trackID,pl,nseg;
	double x, y, angle;
	tree->Branch("trackID",&trackID);
	tree->Branch("x",&x);
	tree->Branch("y",&y);
	tree->Branch("angle",&angle);
	tree->Branch("pl",&pl);
	tree->Branch("nseg",&nseg);

	double x1, y1, z1, x2, y2, z2;
	for(int itrk=0; itrk<ntrk; itrk++){
		EdbTrackP *t = pvr->GetTrack(itrk);
		int nseg = t->N();
		for(int iplate=plMin;iplate<=plMax;iplate++){
			int counts = 0;
			int hitsOnThePlate = 0;
			for(int iseg=0;iseg<nseg;iseg++){
				EdbSegP *s = t->GetSegment(iseg);
				
				if(s->Plate()==iplate-2||s->Plate()==iplate+2) counts++;
				if(s->Plate()==iplate-1){
					x1=s->X();
					y1=s->Y();
					z1=s->Z();
					counts++;
				}
				if(s->Plate()==iplate+1){
					x2=s->X();
					y2=s->Y();
					z2=s->Z();
					counts++;
				}
				if(s->Plate()==iplate)
				{
					hitsOnThePlate=1;
				}
			}
			if(counts==4) {
				angle = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))/(z2-z1);
				h_angle_total->Fill(angle);
				h_plate_total->Fill(iplate);
				if(hitsOnThePlate==1) 
				{
					h_angle_passed->Fill(angle);
					h_plate_passed->Fill(iplate);
				}
				else
				{
					tree->Fill(t->ID(),(x1+x2)/2,(y1+y2)/2,angle,iplate,nseg);
				}
			}
			
		}
		
	}
	tree->Write();
	TCanvas *c = new TCanvas();
	c->Print(Form("efficiency_output/hist_efficiency_%s.pdf[", title.Data()));
	c->SetLogy(1);
	h_angle_total->Draw();
	c->Print(Form("efficiency_output/hist_efficiency_%s.pdf",title.Data()));
	c->SetLogy(0);
	pEff_angle = new TEfficiency(*h_angle_passed, *h_angle_total);
	pEff_angle->SetTitle(Form("Efficiency for each angle (%s);tan#theta;efficiency",title.Data()));
	pEff_angle->Draw();
	c->Print(Form("efficiency_output/hist_efficiency_%s.pdf",title.Data()));
	h_plate_total->Draw();
	c->Print(Form("efficiency_output/hist_efficiency_%s.pdf",title.Data()));
	pEff_plate = new TEfficiency(*h_plate_passed, *h_plate_total);
	pEff_plate->SetTitle(Form("Efficiency for each plate (%s);plate;efficiency", title.Data()));
	pEff_plate->Draw();
	c->Print(Form("efficiency_output/hist_efficiency_%s.pdf",title.Data()));
	c->Print(Form("efficiency_output/hist_efficiency_%s.pdf]", title.Data()));
	TFile fout(Form("efficiency_output/efficiency_%s.root", title.Data()), "recreate");
	pEff_angle->Write();
	pEff_plate->Write();
	fout.Close();

	FILE *ftxt = fopen("efficiency_output/efficiency.txt","w"); //Make txt file for smearing parameters.
	for(int i=1;i<=nbins;i++)
	{
		fprintf(ftxt, "%.3f, ",(bins[i]+bins[i-1])/2);
	}
	fprintf(ftxt,"\n");
	for(int i=1;i<=nbins;i++)
	{
		fprintf(ftxt,"%f ",pEff_angle->GetEfficiency(i));
	}
	fprintf(ftxt,"\n");
	fclose(ftxt);
	return 0;
}