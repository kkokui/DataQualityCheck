#include <FnuQualityCheck.h>

#include <stdio.h>
#include <chrono>

#include <EdbDataSet.h>

int main(int argc, char *argv[])
{
	auto start = std::chrono::system_clock::now();
	if (argc < 5)
	{
		printf("Usage: ./quality_check linked_tracks.root title Xcenter Ycenter\n");
		return 1;
	}

	TString filename_linked_tracks = argv[1];
	TString title = argv[2];
	double Xcenter, Ycenter, bin_width;
	bin_width = 20000;
	sscanf(argv[3], "%lf", &Xcenter);
	sscanf(argv[4], "%lf", &Ycenter);
	// sscanf(argv[5], "%lf", &bin_width);

	EdbDataProc *dproc = new EdbDataProc;
	EdbPVRec *pvr = new EdbPVRec;

	dproc->ReadTracksTree(*pvr, filename_linked_tracks, "nseg>=2");
	// dproc->ReadTracksTree(*pvr, filename_linked_tracks, "npl>=64");
	// dproc->ReadTracksTree(*pvr, filename_linked_tracks, "Entry$%20000==0");
	// dproc->ReadTracksTree(*pvr, filename_linked_tracks, "Entry$<5000");

	TObjArray *tracks = pvr->GetTracks();
	int ntrk = tracks->GetEntriesFast();

	if (ntrk == 0)
	{
		printf("ntrk==0\n");
		return 0;
	}
	FnuQualityCheck qc(pvr, title);
	qc.Summarize(Xcenter, Ycenter, bin_width);

	auto end = std::chrono::system_clock::now();
	auto dur = end - start;
	auto secondsPassed = std::chrono::duration_cast<std::chrono::seconds>(dur).count();
	std::cout << secondsPassed << "seconds\n";
	return 0;
}
