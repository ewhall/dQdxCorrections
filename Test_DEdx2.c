#include "DEdx_Test.h"
#include "Cal_Test.h"

int Test_DEdx2() {
	std::string track_file = "/Users/mwhall/Tools/DLLEEdata/MCC9/5E19/tracker_reco.root";
	std::string vertex_variables_file = "/Users/mwhall/Tools/DLLEEdata/MCC9/5E19/dllee_vertex.root";
	std::string truth_file = "/Users/mwhall/Tools/DLLEEdata/MCC9/Overlay/MC_Info.root";
	std::string correction_file = "/Users/mwhall/dllee_unified/LArCV/app/CalibMaps/CalibrationMaps_MCC9.root";
	std::string spline_file = "/Users/mwhall/dllee_unified/LArCV/app/Reco3D/Proton_Muon_Range_dEdx_LAr_TSplines.root";

	using namespace larlite;
	larlite::storage_manager mgr;	

	mgr.set_io_mode(mgr.kREAD);
	mgr.add_in_filename(track_file);
	if(!mgr.open()) {
		std::cerr << "Failed to open ROOT file. Aborting." << std::endl;
	    return 1;
	}

	std::string dir_name = "test map";
	DEdx_Test DEDX;

	//TFile f("CalibrationFile.root","UPDATE");
	// f.mkdir(dir_name);
	// f.cd(dir_name);

	DEDX.Initialize();
	DEDX.LoadSelectionFiles(vertex_variables_file);
	DEDX.LoadCalibrationFile(correction_file);
	//DEDX.LoadTruthFile(truth_file);
	int x = 0;
	while(mgr.next_event() && x < 1000000) {
		DEDX.AllVertexFill(&mgr);
		x++;
	}

	std::cout << "Loading Splines" << std::endl;

	DEDX.LoadDEdxSplines(spline_file);




	DEDX.Tracks_All();
	DEDX.RecoFit_new();
	//DEDX.All_DEdxPlots();
	//DEDX.LoadCalibrationFile(correction_file);



	return 0;
}