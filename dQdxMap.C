#include "Cal_Test.h"
#include <iostream>
#include <fstream>

int dQdxMap() { 
	std::string track_file = "/Users/mwhall/Tools/DLLEEdata/MCC9/BNB_EXT/tracker_reco.root";
	std::string vertex_variables_file = "/Users/mwhall/Tools/DLLEEdata/MCC9/BNB_EXT/dllee_vertex.root";
	std::string spline_file = "/Users/mwhall/dllee_unified/LArCV/app/Reco3D/Proton_Muon_Range_dEdx_LAr_TSplines.root";


	
	using namespace larlite;
	larlite::storage_manager mgr;

	mgr.set_io_mode(mgr.kREAD);
	mgr.add_in_filename(track_file);
	if(!mgr.open()) {
		std::cerr << "Failed to open ROOT file. Aborting." << std::endl;
	    return 1;
	}

	Cal_Test Cal;
	Cal.Initialize();
	Cal.LoadSelectionFiles(vertex_variables_file);
	int ievt {0};
	int x = 0;
	while(mgr.next_event() && x < 1000000) {
		Cal.AllVertexFill(&mgr);
		if (ievt % 1000 == 0) std::cout << ievt << std::endl;
		ievt++;
		x++;
	}

	Cal.ImageCorrectionMap();





	return 0;
	}
