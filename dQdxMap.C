#include "Cal_Test.h"
#include <iostream>
#include <fstream>

// void LoadSelectionFile();

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

	// TFile *fSelection = TFile::Open(Form("%s",argv[2]),"READ");
	// if(!(fSelection->IsOpen())){
	// 	std::cout << "ERROR : could not open selection file!" << std::endl;
 //    }	
 //    std::cout << "selection file opened" << std::endl;

 //    TTree *NuMuVertexVariable = (TTree*)fSelection->Get("NuMuVertexVariables");
 //    if(NuMuVertexVariable==nullptr)std::cout << "ERROR, no ttree loaded " << std::endl;
 //    else std::cout << "ttree loaded with " <<NuMuVertexVariable->GetEntries() << " entries" << std::endl;
	Cal_Test Cal;
	//Cal.Initialize();
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

	//Cal.RunInfo();
	//Cal.LookAtX();
	//Cal.LookAtZ();
	Cal.LoadDEdxSplines(spline_file);
	Cal.MakedQdxPlots();
	Cal.Make2DCorrectionMap();
	Cal.Finalize();
	// Cal.Make2DCorrectionMap();
	//Cal.ImageCorrectionMap();





	// while (mgr.next_event()) {
	// 	const larlite::event_track *track = (larlite::event_track*)mgr.get_data<larlite::event_track>("trackReco");
	// 	const larlite::vertex *vertex = (larlite::event_vertex*)mgr.get_data<larlite::event_vertex>("trackReco");
	// 	int run = mgr.run_id();
	// 	int run = mgr.subrun_id();
	// 	int event = mgr.event_id();

	// 	for(size_t vtx_idx {0}; vtx_idx < vertex->size(); vtx_idx++) {

	// 	}

	// 	int size = track->size();

	// 	}
	return 0;
	}
