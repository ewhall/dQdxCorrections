#ifndef DEDX_TEST_CXX
#define DEDX_TEST_CXX

#include "DEdx_Test.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TLine.h"
#include "TLatex.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TVirtualFitter.h"
#include "TGraph.h"
#include "TGraphErrors.h"	
#include "TStyle.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "DataFormat/hit.h"
#include "DataFormat/track.h"
#include "DataFormat/vertex.h"
#include "LArUtil/Geometry.h"
#include "LArUtil/GeometryHelper.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>

double RecombinationFunction(double *x, double *par){

	double alpha = par[0];
	double beta = par[1];
	double C = par[2]; //
	double rho = par[4]; // density of liquid argon
	double e_field = par[5]; // electric field of microboone
	double W = par[3]; // ionization energy of liquid argon
	double dedx = x[0];

	double dqdx = (1./C)*((rho*e_field)/(beta*W))*log(dedx*(beta/(rho*e_field))+alpha);

	return dqdx;
}

double dEdxFunction(double *x, double *par){

	double alpha = par[0];
	double beta = par[1];
	double C = par[2]; //
	double rho = par[4]; // density of liquid argon
	double e_field = par[5]; // electric field of microboone
	double W = par[3]; // ionization energy of liquid argon
	double dqdx = x[0];

	double dedx = (exp((dqdx/C)*((beta*W)/(rho*e_field))) - alpha)/(beta/(rho*e_field));

	return dedx;
}

double deltaR(double vx, double vy, double vz, std::vector<double> vec) {
		double dx = vec[0] - vx;
		double dy = vec[1] - vy;
		double dz = vec[2] - vz;
		double dr = sqrt(pow(dx,2.0) + pow(dy,2.0) + pow(dz,2.0));
		return dr;
}

double findE(std::vector<std::vector<double>> truthID) {
	std::vector<double> temp_pdg = truthID[0];
	std::vector<double> temp_energy = truthID[1];
	std::vector<double> temp_energy_p {0};	
	int n = temp_pdg.size();
	for (size_t idx {0}; idx <= n; idx++) {
		if (temp_pdg[idx] == 2212) {temp_energy_p.push_back(temp_energy[idx]);}
	}
	return *std::max_element(temp_energy_p.begin(), temp_energy_p.end());
}

bool InFid(double x, double y, double z) {
	if(x < 0+15 || x > 256.35-15)return 0;
	if(y < -116.5+15 || y > 116.5-25)return 0;
	if(z < 0+15 || z > 1036.8-15)return 0;
	if(y > tan(30*TMath::Pi()/180)*z-115)return 0;
	if(y > 605-tan(30*TMath::Pi()/180)*z && y < 615-tan(30*TMath::Pi()/180)*z)return 0;
	if(y > 415-tan(30*TMath::Pi()/180)*z && y < 430-tan(30*TMath::Pi()/180)*z)return 0;
	if(y > 155-tan(30*TMath::Pi()/180)*z && y < 160-tan(30*TMath::Pi()/180)*z)return 0;
	if(y > 231-tan(30*TMath::Pi()/180)*z && y < 235-tan(30*TMath::Pi()/180)*z)return 0;
	if(y > tan(30*TMath::Pi()/180)*z-200 && y < tan(30*TMath::Pi()/180)*z-180)return 0;
	if(y > tan(30*TMath::Pi()/180)*z-380 && y < tan(30*TMath::Pi()/180)*z-370)return 0;
	if(y > tan(30*TMath::Pi()/180)*z-350 && y < tan(30*TMath::Pi()/180)*z-345)return 0;
	if(y > tan(30*TMath::Pi()/180)*z-555 && y < tan(30*TMath::Pi()/180)*z-549)return 0;
	if(y > tan(30*TMath::Pi()/180)*z-600 && y < tan(30*TMath::Pi()/180)*z-605)return 0;
	if(y > tan(30*TMath::Pi()/180)*z-625 && y < tan(30*TMath::Pi()/180)*z-630)return 0;
	if(z > 53  && z < 57 )return 0;
	if(z > 91  && z < 96 )return 0;
	if(z > 120 && z < 125)return 0;
	if(z > 245 && z < 249)return 0;
	if(z > 287 && z < 293)return 0;
	if(z > 398 && z < 404)return 0;
	if(z > 412 && z < 418)return 0;
	if(z > 700 && z < 740)return 0;
	if(z > 806 && z < 811)return 0;
	if(z > 820 && z < 826)return 0;
	if(z > 873 && z < 879)return 0;
	else return 1;
}



namespace larlite {

	bool DEdx_Test::Initialize() {

		gStyle->SetOptStat(0);

		dEdx_Info = true;

		xmin = 0.0;	xmax = 256.25;
		ymin = -116.5; ymax = 116.5;
		zmin = 0.0; zmax = 1036.8;
		Nz = 300; Ny = 75, Nx = 75; // number of bins to be made along that axis
		plane = larlite::geo::kZ; // change this to expand to 3 planes

		alpha_ex = 0.93;
		beta_ex = 0.212;
		Wion = 23.6e-6;
		E_Field = 0.273;
		LAr_Density = 1.38;
		first_run = 0;
		last_run = 1233567;

		fRecombinationExpected = new TF1("fRecombinationExpected",RecombinationFunction,0,1000,6);
        fRecombinationExpected->SetParameters(alpha_ex,beta_ex,1,Wion,LAr_Density,E_Field);
        fRecombinationExpected->SetParNames("alpha","beta","C","W_{ion}","#rho_{LAr}","E_{field}");
        fRecombinationExpected->FixParameter(0, alpha_ex);
        fRecombinationExpected->FixParameter(1, beta_ex);
        fRecombinationExpected->FixParameter(3, Wion);
        fRecombinationExpected->FixParameter(4,LAr_Density);
        fRecombinationExpected->FixParameter(5,E_Field);

        avgVal =  fRecombinationExpected->Eval(2.104) / 10.584908; //* (11.389267/10.352675);
        avgVal_cor = fRecombinationExpected->Eval(2.104) / 8.45; //* (11.389267/10.352675); 
        std::cout << avgVal << std::endl;

        avgVal_U = fRecombinationExpected->Eval(2.104) / 10.875;
        avgVal_V = fRecombinationExpected->Eval(2.104) / 10.125;
        avgVal_Y = fRecombinationExpected->Eval(2.104) / 10; //9.875;
        avgVal_all = fRecombinationExpected->Eval(2.104) / (32.375/3);
        avgVal_U_Corr = 5975.21; //5183.17;//fRecombinationExpected->Eval(2.104) / 10.875;
        avgVal_V_Corr = 6003.31; //5567.1; ////fRecombinationExpected->Eval(2.104) / 10.125;
        avgVal_Y_Corr = 5559.66; //6351.2; //fRecombinationExpected->Eval(2.104) / 8.875;
        avgVal_all_Corr = 5743.59; //6234.13; //fRecombinationExpected->Eval(2.104) / 9.04167;

        std::cout << "Scaling Factor (Plane 0): " << avgVal_U_Corr << std::endl;
        std::cout << "Scaling Factor (Plane 1): " << avgVal_V_Corr << std::endl;
        std::cout << "Scaling Factor (Plane 2): " << avgVal_Y_Corr << std::endl;
        std::cout << "Scaling Factor (3 Planes Averaged): " << avgVal_all_Corr << std::endl;
		return true;
	}

	bool DEdx_Test::Finalize() {

		UncorrectedTracks();
		CorrectedTracks();
		RecombinationFit();
		MakeDEDXPlots();
		MakeRecoFill();

		TFile f("CalibrationFile.root","RECREATE");

		TCanvas *c1 = new TCanvas("c1","c1",200,10,750,1000);
		c1->Divide(1,3);
		c1->cd(1); hUncorrectedDQdx->Draw("colz");
		c1->cd(2); hUncorrectedDQdx_m->Draw("colz");
		c1->cd(3); hUncorrectedDQdx_p->Draw("colz");

		TCanvas *C = new TCanvas("C","C",200,10,300,500);
		C->cd(); h5e19_2->Draw();

		TCanvas *c2 = new TCanvas("c2","c2",200,10,750,1000);
		c2->Divide(1,3);
		c2->cd(1); hCorrectedDQdx->Draw("colz");
		c2->cd(2); hCorrectedDQdx_m->Draw("colz");
		c2->cd(3); hCorrectedDQdx_p->Draw("colz");

		// // TCanvas *c3 = new TCanvas("c3","c3",200,10,750,300);
		// // c3->cd();
		// // CorrectionMap3D->Draw("BOX2 Z");

		TCanvas *c3 = new TCanvas("c3","c3",200,10,750,1000);
		c3->Divide(1,3);
		c3->cd(1); hCorrectedDEdx->Draw("colz");
		TProfile *prot = hCorrectedDEdx_m->ProfileX("prot");
		TProfile *muo = hCorrectedDEdx_p->ProfileX("muo");
		sMuonRange2dEdx->SetLineColor(3);
		sMuonRange2dEdx->SetLineWidth(3);
		sMuonRange2dEdx->Draw("same");
		sProtonRange2dEdx->SetLineColor(2);
		sProtonRange2dEdx->SetLineWidth(3);
		sProtonRange2dEdx->Draw("same");
		c3->cd(2); hCorrectedDEdx_m->Draw("colz");
		sMuonRange2dEdx->SetLineColor(3);
		sMuonRange2dEdx->SetLineWidth(3);
		prot->SetMarkerColor(6);
		prot->Draw("same");
		sMuonRange2dEdx->Draw("same");
		c3->cd(3); hCorrectedDEdx_p->Draw("colz");
		sProtonRange2dEdx->SetLineColor(2);
		sProtonRange2dEdx->SetLineWidth(3);
		muo->SetMarkerColor(6);
		muo->Draw("same");
		sProtonRange2dEdx->Draw("same");

		TCanvas *c4 = new TCanvas("c4","c4",200,10,750,600);
		c4->Divide(1,2);
		// c4->cd(1);
		// hDQdx_p_cor->SetLineColor(2);
		// hDQdx_p_cor->SetLineWidth(5);
		// hDQdx_p_cor->Draw();
		// hDQdx_m_cor->SetLineColor(3);
		// hDQdx_m_cor->SetLineWidth(5);
		// hDQdx_m_cor->Draw("same");
		c4->cd(2);
		hDQdx_p_un->SetLineColor(2);
		hDQdx_p_un->SetLineWidth(5);
		hDQdx_p_un->Draw();
		hDQdx_m_un->SetLineColor(3);
		hDQdx_m_un->SetLineWidth(5);
		hDQdx_m_un->Draw("same");

		TCanvas *c5 = new TCanvas("c5","c5",200,10,750,600);
		c5->Divide(1,2);
		c5->cd(1);
		hDEdx_m->SetLineColor(3);
		hDEdx_m->SetLineWidth(5);
		hDEdx_m->Draw();
		hDEdx_p->SetLineColor(2);
		hDEdx_p->SetLineWidth(5);
		hDEdx_p->Draw("same");



		// CorrectionMap3D->Write();
		// RecombinationParameters->Write();



		return true;
	}

	void DEdx_Test::LoadCalibrationFile(std::string file_loc) {

		TFile *fCalibration = TFile::Open(Form("%s",file_loc.c_str()),"READ");
		if (!(fCalibration->IsOpen())){
			std::cout << "Error: could not open calibration file" << std::endl;
		}
		std::cout << "calibration file opened" << std::endl;

		CorrectionMap3D = (TH3D*)fCalibration->Get("hImageCalibrationMap_02");
		CorrectionMap3D_U = (TH3D*)fCalibration->Get("hImageCalibrationMap_00");
		CorrectionMap3D_V = (TH3D*)fCalibration->Get("hImageCalibrationMap_01");
		// RecoFit = (TH2D*)fCalibration->Get("RecoFit");
		// hRecoFit = (TH2D*) RecoFit->Clone("hRecoFit");
		// hRecoFit->SetDirectory(0);
		// hRecoFit->SetName("hRecoFit");

		if (CorrectionMap3D==nullptr) std::cout << "Error, no Correction Map loaded" << std::endl;

		int emptyBin {0};
		int fullBin {0};
		for (size_t idx {0}; idx < CorrectionMap3D->GetSize(); idx++) {
			double calFac = CorrectionMap3D->GetBinContent(idx);
			if (calFac == 0.0) emptyBin++;
			if (calFac != 0.0) fullBin++;
			//std::cout << idx << ": " << calFac << std::endl;
		}

		std::cout << emptyBin << " " << fullBin << std::endl;



	}

	void DEdx_Test::LoadRecombinationFile(std::string file_loc) {

		TFile *fRecombination = TFile::Open(Form("%s",file_loc.c_str()),"READ");
		if (!(fRecombination->IsOpen())) {
			std::cout << "Error: could not open recombination file" << std::endl;
		}
		std::cout << "recombination file opened" << std::endl;

		CorrectionMap3D = (TH3D*)fRecombination->Get("CorrectionMap3D");
		if (CorrectionMap3D==nullptr) std::cout << "Error, no correction map loaded" << std::endl;

		TTree *RecombinationParameters = (TTree*)fRecombination->Get("RecombinationParameters");
		if (RecombinationParameters==nullptr) std::cout << "Error, no ttree loaded" << std::endl;


	}

	void DEdx_Test::LoadSelectionFiles(std::string file_loc) {

		TFile *fSelection = TFile::Open(Form("%s",file_loc.c_str()),"READ");
		if(!(fSelection->IsOpen())){
			std::cout << "ERROR : could not open selection file!" << std::endl;
	    }	
	    std::cout << "selection file opened" << std::endl;

	    TTree *NuMuVertexVariable = (TTree*)fSelection->Get("NuMuVertexVariables");
	    if(NuMuVertexVariable==nullptr)std::cout << "ERROR, no ttree loaded " << std::endl;
	    else std::cout << "ttree loaded with " <<NuMuVertexVariable->GetEntries() << " entries" << std::endl;

	    int run_sel, subrun_sel, event_sel, vtxid_sel, N5cmTracks_sel, InFiducial_sel, PassCuts_sel, goodReco_sel;
	    float cosmicLL_sel;
	    std::vector<double> *TrackLength_v_sel = 0;

	    NuMuVertexVariable->SetBranchAddress("run",&run_sel);
		NuMuVertexVariable->SetBranchAddress("subrun",&subrun_sel);
		NuMuVertexVariable->SetBranchAddress("event",&event_sel);
		NuMuVertexVariable->SetBranchAddress("vtxid",&vtxid_sel);

	    NuMuVertexVariable->SetBranchAddress("InFiducial",&InFiducial_sel);
	    NuMuVertexVariable->SetBranchAddress("Good3DReco",&goodReco_sel);
	    NuMuVertexVariable->SetBranchAddress("N5cmTracks",&N5cmTracks_sel);
		NuMuVertexVariable->SetBranchAddress("TrackLength_v",&TrackLength_v_sel);
		NuMuVertexVariable->SetBranchAddress("CosmicLL",&cosmicLL_sel);
		NuMuVertexVariable->SetBranchAddress("PassCuts",&PassCuts_sel);


	    for (Long64_t i = 0; i < NuMuVertexVariable->GetEntries();i++) { 
	    	NuMuVertexVariable->GetEntry(i);
	    	//std::cout << "pass" << std::endl;
	    	if (InFiducial_sel != 1 || goodReco_sel == 0 || N5cmTracks_sel != 2 || TrackLength_v_sel->size() != 2) continue;
	    	//std::cout << InFiducial_sel << std::endl;
	    	if (TrackLength_v_sel->at(0) > 35 && TrackLength_v_sel->at(1) > 35 
	    		&& TrackLength_v_sel->at(0)+TrackLength_v_sel->at(1) > 200
	    		&& cosmicLL_sel < 3) {

	    		std::vector<int> info;
	    		info.push_back(run_sel);
	    		info.push_back(subrun_sel);
	    		info.push_back(event_sel);
	    		info.push_back(vtxid_sel);
	    		SelectEvtID_cosmic.push_back(info);
	    	} 
	    	else if (cosmicLL_sel > 4 && PassCuts_sel == 1) {
	    		std::vector<int> info;
	    		info.push_back(run_sel);
	    		info.push_back(subrun_sel);
	    		info.push_back(event_sel);
	    		info.push_back(vtxid_sel);
	    		SelectEvtID_proton.push_back(info);

	    	}
	    }
    	std::cout << SelectEvtID_cosmic.size() << " & " << SelectEvtID_proton.size() << std::endl;
        std::cout << "tree looped through" << std::endl;
        fSelection->Close();
	}

	void DEdx_Test::AllVertexFill(storage_manager * mgr) {

		//std::cout << "Filling Vertex Data" << std::endl;

		const larlite::event_vertex *vertex_v = (larlite::event_vertex *)mgr->get_data<larlite::event_vertex>("trackReco_sceadded");
		run = mgr->run_id();
		subrun = mgr->subrun_id();
		event = mgr->event_id();
		vtx_id = -1;

		if(vertex_v->size() == 0){return;}

        larlite::event_track* ev_trk=nullptr;
    	auto const& vtx_to_trk = mgr->find_one_ass(vertex_v->id(), ev_trk, vertex_v->id().second);
    	if(!ev_trk) throw larlite::DataFormatException("Could not find associated track data product!");

    	for(int vertex_index=0;vertex_index<vertex_v->size();vertex_index++){
            vtx_id = vertex_index;
            selectedProton = false;
            selectedCosmic = false;
            std::string key_info = std::to_string(run) + "_" + std::to_string(subrun) + "_" + std::to_string(event) + "_" + std::to_string(vtx_id);
            std::string truth_info_key = std::to_string(run) + "_" + std::to_string(subrun) + "_" + std::to_string(event); 
            if (dEdx_Info) dEdxInfo = dEdx_Map[key_info];
            if (dEdx_Info) TruthInfo = Truth_Map[truth_info_key];

        	SelectedVertex(run, subrun,event,vtx_id);


            //if(!selectedCosmic)continue;
            if(!selectedProton && !selectedCosmic)continue;

            if(thisVertex.size()!=0)thisVertex.clear();
            double sumLength = 0;
            for(auto const& trk_index : vtx_to_trk[vertex_index]) {
                thisVertex.push_back( (*ev_trk)[trk_index]);
                sumLength+=(*ev_trk)[trk_index].Length();
            }
            if(selectedCosmic){AllVertex_m.push_back(thisVertex);}
            if(selectedProton){
            	AllVertex_p.push_back(thisVertex);
            	if (dEdx_Info) dEdxInfo_p.push_back(dEdxInfo);
            	if (dEdx_Info) TruthInfo_p.push_back(TruthInfo);
            	std::vector<double> info;
            	if (SelectTruthID_proton.size() != 0) {
	            	for (size_t i {0}; i < SelectTruthID_proton.size(); i++) {
	            		if (run != SelectTruthID_proton[i][0]) continue;
	            		if (subrun != SelectTruthID_proton[i][1]) continue;
	            		if (event != SelectTruthID_proton[i][2]) continue;
	            		double vtx_x = SelectTruthID_proton[i][3];
	            		double vtx_y = SelectTruthID_proton[i][4];
	            		double vtx_z = SelectTruthID_proton[i][5];
	            		info.push_back(vtx_x);
	            		info.push_back(vtx_y);
	            		info.push_back(vtx_z);
	            		//std::cout << "Finding proton truth energy" << std::endl;
	            		if (SelectTruthID_energy[i][1].size() != 0) {
	            			double nrg = findE(SelectTruthID_energy[i]);
	            			info.push_back(nrg);
	            		}
	            		else {info.push_back(0.0);}
					}
					AllTruth_p.push_back(info);
            	}
            }

    	}// vtx_id

    	//std::cout << AllVertex_p.size() << " and " << AllTruth_p.size() << std::endl;
    	return;
	}

	void DEdx_Test::SelectedVertex(int run, int subrun, int event, int vtx_id) {

		selectedProton = false;
        selectedCosmic = false;
        for(size_t i = 0;i<SelectEvtID_cosmic.size();i++){
            if(run!= SelectEvtID_cosmic[i][0])continue;
            if(subrun!=SelectEvtID_cosmic[i][1])continue;
            if(event!=SelectEvtID_cosmic[i][2])continue;
            if(vtx_id!=SelectEvtID_cosmic[i][3])continue;
            if((run >= 5119 && run <=5600) || (run >= 5750 && run <= 5955)) {
            selectedCosmic=true;
        }
        }

        for(size_t i = 0;i<SelectEvtID_proton.size();i++){
            if(run!= SelectEvtID_proton[i][0])continue;
            if(subrun!=SelectEvtID_proton[i][1])continue;
            if(event!=SelectEvtID_proton[i][2])continue;
            if(vtx_id!=SelectEvtID_proton[i][3])continue;
            selectedProton=true;
        }

        return;
	}

	void DEdx_Test::LoadDEdxSplines(std::string file_loc) {

		std::cout << "Load dEdx Splines." << std::endl;
		TFile *fSplines = TFile::Open(Form("%s",file_loc.c_str()),"READ");
		sProtonRange2dEdx = (TSpline3*)fSplines->Get("sProtonRange2dEdx");
		sMuonRange2dEdx = (TSpline3*)fSplines->Get("sMuonRange2dEdx");
		sProtonRange2T = (TSpline3*)fSplines->Get("sProtonRange2T");

	}

	void DEdx_Test::UncorrectedTracks() {
		TH2D* hUnDQdx[3];
		TH2D* hUnDQdx_p[3];
		TH2D* hUnDQdx_m[3];
		for (size_t iPlane {0}; iPlane < 3; iPlane++) {
			hUnDQdx[iPlane] = new TH2D(Form("hUnDQdx_%zu",iPlane),Form("hUnDQdx_%zu",iPlane),400,0,400,300,500000);
			hUnDQdx_p[iPlane] = new TH2D(Form("hUnDQdx_p_%zu",iPlane),Form("hUnDQdx_p_%zu",iPlane),400,0,400,300,500000);
			hUnDQdx_m[iPlane] = new TH2D(Form("hUnDQdx_m_%zu",iPlane),Form("hUnDQdx_m_%zu",iPlane),400,0,400,300,500000);

		}
		hUncorrectedDQdx = new TH2D("UncorrectedDQdx","UncorrectedDQdx", 400, 0, 400, 300, 0, 500000);
		hUncorrectedDQdx_p = new TH2D("UncorrectedDQdx_p","UncorrectedDQdx_p", 400, 0, 400, 300, 0, 500000);
		hUncorrectedDQdx_m = new TH2D("UncorrectedDQdx_m","UncorrectedDQdx_m", 400, 0, 400, 300, 0, 500000);
		hDQdx_p_un = new TH1D("DQdx_un_p","DQdx_un_p",20,0,10);
		hDQdx_m_un = new TH1D("DQdx_un_m","DQdx_un_m",20,0,10);
		hdr = new TH1D("dr","dr",100, 0, 20);
		hEnergy = new TH1D("proton_energy","proton_energy",100,0,1000);
		hEnergy_p = new TH2D("Energy_p","Reconstructed_v_Truth_Energy",50,0,1000,50,0,1000);
		hWeird = new TH3D("Weird","",Nz, zmin, zmax, Nx, xmin, xmax, Ny, ymin, ymax);
		hUncorrectedADC = new TH1D("ADC","ADC",200,0,50);
		h5e19_2 = new TH1D("h5e19_2","h5e19_2",Nx,xmin,xmax);
		h5e19_2_raw = new TH1D("h5e19_2_raw","h5e19_2_raw",Nx,xmin,xmax);

		hUncorrectedRecombinationPlot = new TH2D("UnCorrectedRecombinationPlot","UncorrectedRecombinationPlot", 60, 0, 15, 50, 0, 500000);
		
		std::vector<std::vector<double>> binFill {};
		std::vector<double> totalCount {};
		binFill.resize(h5e19_2_raw->GetSize());

		for ( size_t ivert = 0; ivert < AllVertex_p.size(); ivert++) {
			double avg_dQdx[2] {0,0};
			double avg_dQdx_10cm[2] {0,0};
			double v_x = AllVertex_p[ivert][0].Vertex().X();
			double v_y = AllVertex_p[ivert][0].Vertex().Y();
			double v_z = AllVertex_p[ivert][0].Vertex().Z();
			double dr = deltaR(v_x, v_y, v_z, AllTruth_p[ivert]);
			bool InFiducial = InFid(v_x, v_y, v_z);
			double dr = 1;
			hdr->Fill(dr);
			// hEnergy->Fill(AllTruth_p[ivert][3]);




			if (dr > 5 || InFiducial == 0) continue;
				for ( size_t itrack = 0; itrack < AllVertex_p[ivert].size(); itrack++) {
					int points {0};
					int points_10cm {0};
					auto curr_vert = AllVertex_p[ivert][itrack];
					for (size_t iloc = 0; iloc < curr_vert.NumberTrajectoryPoints(); iloc++) {
						double raw_dQdx = curr_vert.DQdxAtPoint(iloc,plane);
						TVector3 point = curr_vert.LocationAtPoint(iloc);
						double reslength = curr_vert.Length(iloc);
						if (reslength < 1) continue;
						if (raw_dQdx <= 0 || (raw_dQdx < 37.01 && raw_dQdx > 36.99)) continue;
						points += 1;
						avg_dQdx[itrack] += raw_dQdx;
						if (InFid(point.X(), point.Y(), point.Z()) ==  false) continue;
						//hUncorrectedDQdx->Fill(reslength, raw_dQdx*avgVal);
						if (reslength < 10) { 
							avg_dQdx_10cm[itrack] += raw_dQdx; 
							points_10cm += 1;
						}
					}
					avg_dQdx[itrack] /= points;
				}
				int kProton = 0;
				int kMuon = 1;
				if (AllVertex_p[ivert][kProton].Length() < 5 || AllVertex_p[ivert][kMuon].Length() < 5) continue;
				if (AllVertex_p[ivert][kProton].Length() < AllVertex_p[ivert][kMuon].Length()) {
					if (avg_dQdx[kProton] < avg_dQdx[kMuon]) continue;
					hDQdx_p_un->Fill(avg_dQdx_10cm[kProton]*avgVal);
					hDQdx_m_un->Fill(avg_dQdx_10cm[kMuon]*avgVal);
					double length_t = AllVertex_p[ivert][kProton].Length();
					double reco_e = sProtonRange2T->Eval(length_t);
					//hEnergy_p->Fill(AllTruth_p[ivert][3], reco_e);
					//double per_e = (sProtonRange2T->Eval(length_t))/(AllTruth_p[ivert][3]);
					//if (per_e < 0.8 || per_e > 1.2) continue;
					for (size_t iloc = 0; iloc < AllVertex_p[ivert][kProton].NumberTrajectoryPoints(); iloc++) {
						auto curr_vert = AllVertex_p[ivert][kProton];
						double raw_dQdx = curr_vert.DQdxAtPoint(iloc,plane);
						TVector3 point = curr_vert.LocationAtPoint(iloc);
						//hUncorrectedADC->Fill(raw_dQdx);
						if (raw_dQdx < 37.01 && raw_dQdx > 36.99) {
							//std::cout << raw_dQdx << std::endl;
							std::vector<double> info {0};
							info.push_back(point.X());
							info.push_back(point.Y());
							info.push_back(point.Z());
							weird_wires.push_back(info);
						}
						double reslength = curr_vert.Length(iloc);
						if (reslength < 0.25) continue;
						if (raw_dQdx <= 0 || (raw_dQdx == 37 )) continue;
						if (InFid(point.X(), point.Y(), point.Z()) ==  false) continue;
						hUncorrectedDQdx_p->Fill(reslength, raw_dQdx*avgVal);
						hUncorrectedRecombinationPlot->Fill(sProtonRange2dEdx->Eval(reslength),raw_dQdx*avgVal);
					}
					for (size_t iloc = 0; iloc < AllVertex_p[ivert][kMuon].NumberTrajectoryPoints(); iloc++) {
						auto curr_vert = AllVertex_p[ivert][kMuon];
						double raw_dQdx = curr_vert.DQdxAtPoint(iloc,plane);
						TVector3 point = curr_vert.LocationAtPoint(iloc);
						if (raw_dQdx < 37.01 && raw_dQdx > 36.99) {
							//std::cout << raw_dQdx << std::endl;
							std::vector<double> info {0};
							info.push_back(point.X());
							info.push_back(point.Y());
							info.push_back(point.Z());
							weird_wires.push_back(info);
						}
						double reslength = curr_vert.Length(iloc);
						if (reslength < 0.25) continue;
						//if (curr_vert.End().Z() > v_z) continue;
						int binX = h5e19_2_raw->FindBin(point.X());
						binFill[binX].push_back(curr_vert.DQdxAtPoint(iloc,plane));
						totalCount.push_back(curr_vert.DQdxAtPoint(iloc,plane));
						h5e19_2_raw->Fill(point.X(), curr_vert.DQdxAtPoint(iloc,plane));
						
						if (raw_dQdx <= 0 || (raw_dQdx == 37)) continue;
						if (InFid(point.X(), point.Y(), point.Z()) ==  false) continue;
						hUncorrectedDQdx_m->Fill(reslength, raw_dQdx*avgVal);

						//hUncorrectedRecombinationPlot->Fill(sMuonRange2dEdx->Eval(reslength),raw_dQdx*avgVal);
						if (reslength > 20) {
							//hUncorrectedADC->Fill(raw_dQdx);
						}

					}

				}
				kProton = 1;
				kMuon = 0;
				if (AllVertex_p[ivert][kProton].Length() < AllVertex_p[ivert][kMuon].Length()) {
					if (avg_dQdx[kProton] < avg_dQdx[kMuon]) continue;
					hDQdx_p_un->Fill(avg_dQdx_10cm[kProton]*avgVal);
					hDQdx_m_un->Fill(avg_dQdx_10cm[kMuon]*avgVal);
					double length_t = AllVertex_p[ivert][kProton].Length();
					double reco_e = sProtonRange2T->Eval(length_t);
					//hEnergy_p->Fill(AllTruth_p[ivert][3], reco_e);
					//double per_e = (sProtonRange2T->Eval(length_t))/(AllTruth_p[ivert][3]);
					//if (per_e < 0.8 || per_e > 1.2) continue;
					for (size_t iloc = 0; iloc < AllVertex_p[ivert][kProton].NumberTrajectoryPoints(); iloc++) {
						auto curr_vert = AllVertex_p[ivert][kProton];
						double raw_dQdx = curr_vert.DQdxAtPoint(iloc,plane);
						TVector3 point = curr_vert.LocationAtPoint(iloc);
						double reslength = curr_vert.Length(iloc);
						//hUncorrectedADC->Fill(raw_dQdx);
						if (raw_dQdx < 37.01 && raw_dQdx > 36.99) {
							//std::cout << raw_dQdx << std::endl;
							std::vector<double> info {0};
							info.push_back(point.X());
							info.push_back(point.Y());
							info.push_back(point.Z());
							weird_wires.push_back(info);
						}
						if (reslength < 0.25) continue;
						if (raw_dQdx <= 0 || (raw_dQdx == 37)) continue;
						if (InFid(point.X(), point.Y(), point.Z()) ==  false) continue;
						hUncorrectedDQdx_p->Fill(reslength, raw_dQdx*avgVal);
						hUncorrectedRecombinationPlot->Fill(sProtonRange2dEdx->Eval(reslength),raw_dQdx*avgVal);
					}

					for (size_t iloc = 0; iloc < AllVertex_p[ivert][kMuon].NumberTrajectoryPoints(); iloc++) {
						auto curr_vert = AllVertex_p[ivert][kMuon];
						double raw_dQdx = curr_vert.DQdxAtPoint(iloc,plane);
						TVector3 point = curr_vert.LocationAtPoint(iloc);
						double reslength = curr_vert.Length(iloc);
						if (raw_dQdx < 37.01 && raw_dQdx > 36.99) {
							//std::cout << raw_dQdx << std::endl;
							std::vector<double> info {0};
							info.push_back(point.X());
							info.push_back(point.Y());
							info.push_back(point.Z());
							weird_wires.push_back(info);
						}
						if (reslength < 0.25) continue;
						//if (curr_vert.End().Z() > v_z) continue;
						int binX = h5e19_2_raw->FindBin(point.X());
						binFill[binX].push_back(curr_vert.DQdxAtPoint(iloc,plane));
						totalCount.push_back(curr_vert.DQdxAtPoint(iloc,plane));
						h5e19_2_raw->Fill(point.X(), curr_vert.DQdxAtPoint(iloc,plane));
						
						if (raw_dQdx <= 0 || (raw_dQdx == 37)) continue;
						if (InFid(point.X(), point.Y(), point.Z()) ==  false) continue;
						hUncorrectedDQdx_m->Fill(reslength, raw_dQdx*avgVal);
						//hUncorrectedRecombinationPlot->Fill(sMuonRange2dEdx->Eval(reslength),raw_dQdx*avgVal);
						if (reslength > 20) {
							//hUncorrectedADC->Fill(raw_dQdx);
						}
					
				}

			
		}

	}
		for (size_t idx {0}; idx < binFill.size() ; idx++) {
			if (binFill[idx].size() > 5) {
				double avg_dQdx = h5e19_2_raw->GetBinContent(idx) / binFill[idx].size();
				h5e19_2->SetBinContent(idx, avg_dQdx);
			}
		}



	TLine *line = new TLine(0,0,1000,1000);
	TCanvas *dR = new TCanvas("dR","dR", 200, 10, 300, 750);
	dR->Divide(1,2);
	dR->cd(1); hdr->Draw();
	dR->cd(2); hUncorrectedADC->Draw();
	TCanvas *dEng = new TCanvas("dEng","dEng", 200, 10, 300, 750);
	dEng->Divide(1,2);
	dEng->cd(1); hEnergy->Draw();
	dEng->cd(2); hEnergy_p->Draw("colz");
	line->Draw("same"); 

	for (size_t idx {0}; idx < weird_wires.size(); idx++) {
		//std::cout << weird_wires[idx][0] << ", " << weird_wires[idx][1] << ", " << weird_wires[idx][2] << std::endl;
		hWeird->Fill(weird_wires[idx][2], weird_wires[idx][0], weird_wires[idx][1]);
	}
	std::cout << weird_wires.size() << std::endl;
	TCanvas *dWires = new TCanvas("Wires","  ", 200, 10, 300, 750);
	dWires->cd();
	hWeird->Draw("BOX2");;


}


	void DEdx_Test::CorrectedTracks() { 

		std::cout << "Made Reco Plot" << std::endl;
		const int ix {41}; 
		double reco_start {2.75}, reco_end {13};
		double add_len = (reco_end - reco_start) / ix;
		double bin_start = reco_start;
		double bin_end = reco_start + add_len;
		std::vector<double> RecoBins {};
		TH1D* hRecoPlotBin[ix];
		for (int idx {0}; idx < ix; idx++) {
			RecoBins.push_back(bin_end);
			hRecoPlotBin[idx] = new TH1D(Form("hRecoPlotBin_%02d",idx),Form("hRecoPlotBin_%02d",idx),60,0,500000);
			bin_start+=add_len;
			bin_end+=add_len;
		}

		hCorrectedDQdx = new TH2D("CorrectedDQdx","CorrectedDQdx", 400, 0, 400, 300, 0, 500000);
		hCorrectedDQdx_p = new TH2D("CorrectedDQdx_p","CorrectedDQdx_p", 400, 0, 400, 300, 0, 500000);
		hCorrectedDQdx_m = new TH2D("CorrectedDQdx_m","CorrectedDQdx_m", 400, 0, 400, 300, 0, 500000);
		hDQdx_p_cor = new TH1D("DQdx_cor_p","DQdx_cor_p",100,0,500000);
		hDQdx_m_cor = new TH1D("DQdx_cor_m","DQdx_cor_m",100,0,500000);

		hCorrectedRecombinationPlot = new TH2D("CorrectedRecombinationPlot","CorrectedRecombinationPlot", 60, 0, 15, 60, 0, 500000);
		
		for ( size_t ivert = 0; ivert < AllVertex_p.size(); ivert++) {
			double avg_dQdx[2] {0,0};
			double v_x = AllVertex_p[ivert][0].Vertex().X();
			double v_y = AllVertex_p[ivert][0].Vertex().Y();
			double v_z = AllVertex_p[ivert][0].Vertex().Z();
			bool InFiducial = InFid(v_x, v_y, v_z);
			//double dr = deltaR(v_x, v_y, v_z, AllTruth_p[ivert]);
			double dr = 1;

			if (dr > 5 || InFiducial == 0) continue;
			for ( size_t itrack = 0; itrack < AllVertex_p[ivert].size(); itrack++) {
				int points {0};
				auto curr_vert = AllVertex_p[ivert][itrack];
				for (size_t iloc = 0; iloc < curr_vert.NumberTrajectoryPoints(); iloc++) {
					double raw_dQdx = curr_vert.DQdxAtPoint(iloc,plane);
					TVector3 point = curr_vert.LocationAtPoint(iloc);
					double reslength = curr_vert.Length(iloc);
					double bin = CorrectionMap3D->FindBin(point.Z(), point.X(), point.Y());
					double calFac = CorrectionMap3D->GetBinContent(bin);
					double cor_dQdx = raw_dQdx*avgVal_cor*calFac;
					// std::cout << "(" << point.X() << ", " << point.Y() << ", " << point.Z() << ") : " << raw_dQdx << " : " << calFac << std::endl;
					if (reslength < 0.25) continue;
					if (raw_dQdx <= 0 || (raw_dQdx < 37.01 && raw_dQdx > 36.99) || calFac == 0) continue;
					//std::cout << "Passed." << std::endl;
					points += 1;
					avg_dQdx[itrack] += cor_dQdx;
					if (InFid(point.X(), point.Y(), point.Z()) ==  false) continue;
					hCorrectedDQdx->Fill(reslength, cor_dQdx);
				}
				avg_dQdx[itrack] /= points;
			}
			int kProton = 0;
			int kMuon = 1;
			if (AllVertex_p[ivert][kProton].Length() < 5 || AllVertex_p[ivert][kMuon].Length() < 5) continue;
			if (AllVertex_p[ivert][kProton].Length() < AllVertex_p[ivert][kMuon].Length()) {
				if (avg_dQdx[kProton] < avg_dQdx[kMuon]) continue;
				//std::cout << "Proton: " << avg_dQdx[kProton] << " and Muon: " << avg_dQdx[kMuon] << std::endl;
				hDQdx_p_cor->Fill(avg_dQdx[kProton]);
				hDQdx_m_cor->Fill(avg_dQdx[kMuon]);
				for (size_t iloc = 0; iloc < AllVertex_p[ivert][kProton].NumberTrajectoryPoints(); iloc++) {
					auto curr_vert = AllVertex_p[ivert][kProton];
					double raw_dQdx = curr_vert.DQdxAtPoint(iloc,plane);
					TVector3 point = curr_vert.LocationAtPoint(iloc);
					double reslength = curr_vert.Length(iloc);
					int bin = CorrectionMap3D->FindBin(point.Z(), point.X(), point.Y());
					double calFac_p = CorrectionMap3D->GetBinContent(bin);
					double cor_dQdx = raw_dQdx*calFac_p*avgVal_cor;
					if (reslength < 0.25) continue;
					if (cor_dQdx <= 0 || (raw_dQdx < 37.01 && raw_dQdx > 36.99)) continue;
					if (InFid(point.X(), point.Y(), point.Z()) ==  false) continue;
					hCorrectedDQdx_p->Fill(reslength, cor_dQdx);
					double dEdx_theory = sProtonRange2dEdx->Eval(reslength);
					hCorrectedRecombinationPlot->Fill(dEdx_theory,cor_dQdx);
					int index {1};
					int bin_num {1};
					//std::cout << "check which bin" << std::endl;
					while (dEdx_theory > RecoBins[index] && index < 40) {
						bin_num = index;
						index++;

					}
					//std::cout << "put in bin" << std::endl;
					//std::cout << "Bin number = " << bin_num << "and dEdx = " << dEdx_theory << "and RecoBins = "<< RecoBins[index] << std::endl;
					hRecoPlotBin[bin_num]->Fill(cor_dQdx);
					//std::cout << "we put in bin" << std::endl;

				}
				for (size_t iloc = 0; iloc < AllVertex_p[ivert][kMuon].NumberTrajectoryPoints(); iloc++) {
					auto curr_vert = AllVertex_p[ivert][kMuon];
					double raw_dQdx = curr_vert.DQdxAtPoint(iloc,plane);
					TVector3 point = curr_vert.LocationAtPoint(iloc);
					double reslength = curr_vert.Length(iloc);
					int bin = CorrectionMap3D->FindBin(point.Z(), point.X(), point.Y());
					double calFac_m = CorrectionMap3D->GetBinContent(bin);
					double cor_dQdx = raw_dQdx*calFac_m*avgVal_cor;
					if (reslength < 0.25) continue;
					if (cor_dQdx <= 0 || (raw_dQdx < 37.01 && raw_dQdx > 36.99)) continue;
					if (InFid(point.X(), point.Y(), point.Z()) ==  false) continue;
					hCorrectedDQdx_m->Fill(reslength, cor_dQdx);
					// if (reslength < 5 || reslength > 75) continue;
					//hCorrectedRecombinationPlot->Fill(sMuonRange2dEdx->Eval(reslength),cor_dQdx);
					if (reslength > 50) {
						hUncorrectedADC->Fill(raw_dQdx*calFac_m);
						}
				
			
	}

		kProton = 1;
		kMuon = 0;
		if (AllVertex_p[ivert][kProton].Length() < AllVertex_p[ivert][kMuon].Length()) {
			if (avg_dQdx[kProton] < avg_dQdx[kMuon]) continue;
			//std::cout << "Proton: " << avg_dQdx[kProton] << " and Muon: " << avg_dQdx[kMuon] << std::endl;
			hDQdx_p_cor->Fill(avg_dQdx[kProton]);
			hDQdx_m_cor->Fill(avg_dQdx[kMuon]);
			for (size_t iloc = 0; iloc < AllVertex_p[ivert][kProton].NumberTrajectoryPoints(); iloc++) {
				auto curr_vert = AllVertex_p[ivert][kProton];
				double raw_dQdx = curr_vert.DQdxAtPoint(iloc,plane);
				TVector3 point = curr_vert.LocationAtPoint(iloc);
				double reslength = curr_vert.Length(iloc);
				int bin = CorrectionMap3D->FindBin(point.Z(), point.X(), point.Y());
				double calFac_p = CorrectionMap3D->GetBinContent(bin);
				double cor_dQdx = raw_dQdx*calFac_p*avgVal_cor;
				if (reslength < 0.25) continue;
				if (cor_dQdx <= 0 || (raw_dQdx < 37.01 && raw_dQdx > 36.99)) continue;
				if (InFid(point.X(), point.Y(), point.Z()) ==  false) continue;
				hCorrectedDQdx_p->Fill(reslength, cor_dQdx);
				double dEdx_theory = sProtonRange2dEdx->Eval(reslength);
				hCorrectedRecombinationPlot->Fill(sProtonRange2dEdx->Eval(reslength),cor_dQdx);
				int index {1};
				int bin_num {1};
				//std::cout << "check which bin" << std::endl;
				while (dEdx_theory < RecoBins[index] && index < 40) {
					bin_num = index;
					index++;
				}
				//std::cout << "put in bin" << std::endl;
				hRecoPlotBin[bin_num]->Fill(cor_dQdx);
				//std::cout << "Bin number = " << bin_num << "and dEdx = " << dEdx_theory << std::endl;
				//std::cout << "we put in bin" << std::endl;

			}
			for (size_t iloc = 0; iloc < AllVertex_p[ivert][kMuon].NumberTrajectoryPoints(); iloc++) {
				auto curr_vert = AllVertex_p[ivert][kMuon];
				double raw_dQdx = curr_vert.DQdxAtPoint(iloc,plane);
				TVector3 point = curr_vert.LocationAtPoint(iloc);
				double reslength = curr_vert.Length(iloc);
				int bin = CorrectionMap3D->FindBin(point.Z(), point.X(), point.Y());
				double calFac_m = CorrectionMap3D->GetBinContent(bin);
				double cor_dQdx = raw_dQdx*calFac_m*avgVal_cor;
				if (reslength < 0.25) continue;
				if (cor_dQdx <= 0 || (raw_dQdx < 37.01 && raw_dQdx > 36.99)) continue;
				if (InFid(point.X(), point.Y(), point.Z()) ==  false) continue;
				hCorrectedDQdx_m->Fill(reslength, cor_dQdx);
				// if (reslength < 5 || reslength > 75) continue;
				//hCorrectedRecombinationPlot->Fill(sMuonRange2dEdx->Eval(reslength),cor_dQdx);
				if (reslength > 50) {
						hUncorrectedADC->Fill(raw_dQdx*calFac_m);
					}
				}
			}
		}


		}
		//int ix = 16;
		TCanvas *cRecoBins = new TCanvas("cRecoBins","cRecoBins",800,800);
		cRecoBins->Divide(5,5);
		int in = 1;
		for (int idx{5}; idx < ix-7; idx++) {
			double mean = hRecoPlotBin[idx]->GetMean();
			double rms = hRecoPlotBin[idx]->GetRMS();
			TF1 *f1 = new TF1("f1","gaus",0,500000);
			hRecoPlotBin[idx]->Fit("f1","","",mean-rms, mean+rms);
			if (idx < 23) {
				cRecoBins->cd(in);
				hRecoPlotBin[idx]->Draw();
				f1->Draw("same");
				in++;
			}
				mu.push_back(f1->GetParameter(1));
				mu_error.push_back(f1->GetParError(1));
				sigma.push_back(f1->GetParameter(2));
				sigma_error.push_back(f1->GetParError(2));
			double val = RecoBins[idx] - 0.25/2;
			ebin.push_back(val);
			ey.push_back(0.25);

			}
		TCanvas *ADC = new TCanvas("ADC","ADC",400,400);
		ADC->cd();
		//double mean = hUncorrectedADC->GetMean();
		//double rms = hUncorrectedADC->GetRMS();
		//TF1 *f2 = new TF1("f2","gaus",4,mean+rms);
		//hUncorrectedADC->Fit("f2","","",4, mean+rms);
		hUncorrectedADC->Draw();
		//f2->Draw("same");
		TFile *ADCFile = new TFile("ADCFile.root","NEW");
		hUncorrectedADC->Write("hCorrectedADC");
		ADCFile->Close();
 
	}

	void DEdx_Test::RecombinationFit() {

		fRecombination = new TF1("Recombination",RecombinationFunction,0,1000,6);
		fRecombination->SetParameters(alpha_ex, beta_ex, 1, Wion, LAr_Density, E_Field);
		fRecombination->SetParNames("alpha","beta","C","W_{ion}","#rho_{LAr}","E_{field}");
		fRecombination->FixParameter(2,1);
		fRecombination->FixParameter(3,Wion);
		fRecombination->FixParameter(4,LAr_Density);
		fRecombination->FixParameter(5,E_Field);

		TGraphErrors *g2 = new TGraphErrors(ebin.size(),&ebin[0],&mu[0],&ey[0],&mu_error[0]);

		TCanvas *cRecombination = new TCanvas("Recombination","Recombination",800,600);
		TF1 *fgaus = new TF1("fgaus","gaus(0)",0,1000);
		fgaus->SetParameters(100,200,20,10,200,100);

		TProfile *P_cor = hCorrectedRecombinationPlot->ProfileX("P_cor");
		g2->Fit(fRecombination,"re","e",2.75,15);
		hCorrectedRecombinationPlot->Draw("colz");
		fRecombination->SetLineWidth(3);
		fRecombination->Draw("same");
		fRecombinationExpected->SetLineColor(3);
		fRecombinationExpected->SetLineWidth(3);
		fRecombinationExpected->Draw("same");
		//P_cor->Draw("sames");
		g2->SetLineColor(6);
		g2->Draw("same");

		

		RecombinationParameters = new TTree("RecombinationParameters","RecombinationParameters");
		RecombinationParameters->Branch("alpha",&alpha_fit);
		RecombinationParameters->Branch("beta",&beta_fit);
		RecombinationParameters->Branch("alpha_sigma",&alpha_sigma);
		RecombinationParameters->Branch("beta_sigma",&beta_sigma);

		alpha_fit = fRecombination->GetParameter(0);
		beta_fit = fRecombination->GetParameter(1);
		alpha_sigma = fRecombination->GetParError(0);
		beta_sigma = fRecombination->GetParError(1);

		RecombinationParameters->Fill();


		//RooDouble *fuck = new RooDouble(alpha_fit);

		// fdEdxConvert = new TF1("fDEdxConvert",dEdxFunction,0,1000000,6);
		// fdEdxConvert->SetParameters(alpha_fit,beta_fit,1,Wion,LAr_Density,E_Field);

		TCanvas *cRecombination_un = new TCanvas("Recombination_Un","Recombination_Un",800,600);

		fRecombination_un = new TF1("Recombination_un",RecombinationFunction,0,1000,6);
		fRecombination_un->SetParameters(alpha_ex, beta_ex, 1, Wion, LAr_Density, E_Field);
		fRecombination_un->SetParNames("alpha","beta","C","W_{ion}","#rho_{LAr}","E_{field}");
		fRecombination_un->SetParLimits(0,0,1);
		fRecombination_un->SetParLimits(1,0,1);

		fRecombination_un->FixParameter(2,1);
		fRecombination_un->FixParameter(3,Wion);
		fRecombination_un->FixParameter(4,LAr_Density);
		fRecombination_un->FixParameter(5,E_Field);

		TProfile *P = hUncorrectedRecombinationPlot->ProfileX("P");
		P->Fit(fRecombination_un,"re","e",2.75,15);
		hUncorrectedRecombinationPlot->Draw("colz");
		fRecombination_un->SetLineWidth(4);
		fRecombination_un->Draw("same");		
		fRecombinationExpected->SetLineColor(3);
		fRecombinationExpected->SetLineWidth(4);
		fRecombinationExpected->Draw("same");
		P->SetMarkerSize(5);
		P->Draw("sames");

		alpha_fit_un = fRecombination_un->GetParameter(0);
		beta_fit_un = fRecombination_un->GetParameter(1);

		fdEdxConvert = new TF1("fDEdxConvert",dEdxFunction,0,1000000,6);
		fdEdxConvert->SetParameters(alpha_fit,beta_fit,1,Wion,LAr_Density,E_Field);


	}

	void DEdx_Test::MakeDEDXPlots() {
		hCorrectedDEdx = new TH2D("CorrectedDEdx","CorrectedDEdx",400,0,400,300,0,30);
		hCorrectedDEdx_p = new TH2D("CorrectedDEdx_p","CorrectedDEdx_p",400,0,400,300,0,30);
		hCorrectedDEdx_m = new TH2D("CorrectedDEdx_m","CorrectedDEdx_m",400,0,400,300,0,30);
		hDEdx_p = new TH1D("DEdx_p","DQdx_p",40,0,30);
		hDEdx_m = new TH1D("DEdx_m","DQdx_m",40,0,30);


		for ( size_t ivert = 0; ivert < AllVertex_p.size(); ivert++) {
			double v_x = AllVertex_p[ivert][0].Vertex().X();
			double v_y = AllVertex_p[ivert][0].Vertex().Y();
			double v_z = AllVertex_p[ivert][0].Vertex().Z();
			bool InFiducial = InFid(v_x, v_y, v_z);
			// double dr = deltaR(v_x, v_y, v_z, AllTruth_p[ivert]);
			double dr = 1;

			if (dr > 5 || InFiducial == 0) continue;
			double avg_dEdx[2] {0,0};
			double dEdx_5cm[2] {0,0};
			for ( size_t itrack = 0; itrack < AllVertex_p[ivert].size(); itrack++) {
				int points {0};
				int points_5cm {0};
				auto curr_vert = AllVertex_p[ivert][itrack];
				for (size_t iloc = 0; iloc < curr_vert.NumberTrajectoryPoints(); iloc++) {
					double raw_dQdx = curr_vert.DQdxAtPoint(iloc,plane);
					TVector3 point = curr_vert.LocationAtPoint(iloc);
					double reslength = curr_vert.Length(iloc);
					double bin = CorrectionMap3D->FindBin(point.Z(), point.X(), point.Y());
					double calFac = CorrectionMap3D->GetBinContent(bin);
					double cor_dQdx = raw_dQdx*avgVal_cor*calFac;
					double cor_dEdx = fdEdxConvert->Eval(cor_dQdx);
					// std::cout << "(" << point.X() << ", " << point.Y() << ", " << point.Z() << ") : " << raw_dQdx << " : " << calFac << std::endl;
					if (reslength < 0.25) continue;
					if (cor_dQdx <= 0 || (raw_dQdx < 37.01 && raw_dQdx > 36.99)) continue;
					if (InFid(point.X(), point.Y(), point.Z()) ==  false) continue;
					//std::cout << "Passed." << std::endl;
					points += 1;
					avg_dEdx[itrack] += cor_dEdx;
					hCorrectedDEdx->Fill(reslength, cor_dEdx);
					if (reslength < 10) { 
						dEdx_5cm[itrack] +=cor_dEdx; 
						points_5cm += 1;
					}
				}
				avg_dEdx[itrack] /= points;
				dEdx_5cm[itrack] /= points_5cm;
			}
			int kProton = 0;
			int kMuon = 1;
			if (AllVertex_p[ivert][kProton].Length() < 5 || AllVertex_p[ivert][kMuon].Length() < 5) continue;
			if (AllVertex_p[ivert][kProton].Length() < AllVertex_p[ivert][kMuon].Length()) {
				if (avg_dEdx[kProton] < avg_dEdx[kMuon]) continue;
				hDEdx_p->Fill(dEdx_5cm[kProton]);
				hDEdx_m->Fill(dEdx_5cm[kMuon]);

				for (size_t iloc = 0; iloc < AllVertex_p[ivert][kProton].NumberTrajectoryPoints(); iloc++) {
					auto curr_vert = AllVertex_p[ivert][kProton];
					double raw_dQdx = curr_vert.DQdxAtPoint(iloc,plane);
					TVector3 point = curr_vert.LocationAtPoint(iloc);
					double reslength = curr_vert.Length(iloc);
					int bin = CorrectionMap3D->FindBin(point.Z(), point.X(), point.Y());
					double calFac_p = CorrectionMap3D->GetBinContent(bin);
					double cor_dQdx = raw_dQdx*avgVal_cor*calFac_p;
					double cor_dEdx = fdEdxConvert->Eval(cor_dQdx);
					if (reslength < 2) continue;
					if (cor_dQdx <= 0 || (raw_dQdx < 37.01 && raw_dQdx > 36.99)) continue;
					if (InFid(point.X(), point.Y(), point.Z()) ==  false) continue;
					hCorrectedDEdx_p->Fill(reslength, cor_dEdx);
					//hCorrectedRecombinationPlot->Fill(sProtonRange2dEdx->Eval(reslength),cor_dEdx);

				}
				for (size_t iloc = 0; iloc < AllVertex_p[ivert][kMuon].NumberTrajectoryPoints(); iloc++) {
					auto curr_vert = AllVertex_p[ivert][kMuon];
					double raw_dQdx = curr_vert.DQdxAtPoint(iloc,plane);
					TVector3 point = curr_vert.LocationAtPoint(iloc);
					double reslength = curr_vert.Length(iloc);
					int bin = CorrectionMap3D->FindBin(point.Z(), point.X(), point.Y());
					double calFac_m = CorrectionMap3D->GetBinContent(bin);
					double cor_dQdx = raw_dQdx*avgVal_cor*calFac_m;
					double cor_dEdx = fdEdxConvert->Eval(cor_dQdx);
					if (cor_dQdx <= 0 || (raw_dQdx < 37.01 && raw_dQdx > 36.99)) continue;
					if (InFid(point.X(), point.Y(), point.Z()) ==  false) continue;
					hCorrectedDEdx_m->Fill(reslength, cor_dEdx);
					// if (reslength < 5 || reslength > 75) continue;
					// hCorrectedRecombinationPlot->Fill(sMuonRange2dEdx->Eval(reslength),cor_dQdx);
					}
			}

			kProton = 1;
			kMuon = 0;
			if (AllVertex_p[ivert][kProton].Length() < AllVertex_p[ivert][kMuon].Length()) {
				if (avg_dEdx[kProton] < avg_dEdx[kMuon]) continue;
				//std::cout << "Proton: " << avg_dQdx[kProton] << " and Muon: " << avg_dQdx[kMuon] << std::endl;
				hDEdx_p->Fill(dEdx_5cm[kProton]);
				hDEdx_m->Fill(dEdx_5cm[kMuon]);
				for (size_t iloc = 0; iloc < AllVertex_p[ivert][kProton].NumberTrajectoryPoints(); iloc++) {
					auto curr_vert = AllVertex_p[ivert][kProton];
					double raw_dQdx = curr_vert.DQdxAtPoint(iloc,plane);
					TVector3 point = curr_vert.LocationAtPoint(iloc);
					double reslength = curr_vert.Length(iloc);
					int bin = CorrectionMap3D->FindBin(point.Z(), point.X(), point.Y());
					double calFac_p = CorrectionMap3D->GetBinContent(bin);
					double cor_dQdx = raw_dQdx*avgVal_cor*calFac_p;
					double cor_dEdx = fdEdxConvert->Eval(cor_dQdx);
					if (reslength < 1) continue;
					if (cor_dQdx <= 0 || (raw_dQdx < 37.01 && raw_dQdx > 36.99)) continue;
					if (InFid(point.X(), point.Y(), point.Z()) ==  false) continue;
					hCorrectedDEdx_p->Fill(reslength, cor_dEdx);

				}
				for (size_t iloc = 0; iloc < AllVertex_p[ivert][kMuon].NumberTrajectoryPoints(); iloc++) {
					auto curr_vert = AllVertex_p[ivert][kMuon];
					double raw_dQdx = curr_vert.DQdxAtPoint(iloc,plane);
					TVector3 point = curr_vert.LocationAtPoint(iloc);
					double reslength = curr_vert.Length(iloc);
					int bin = CorrectionMap3D->FindBin(point.Z(), point.X(), point.Y());
					double calFac_m = CorrectionMap3D->GetBinContent(bin);
					double cor_dQdx = raw_dQdx*avgVal_cor*calFac_m;
					double cor_dEdx = fdEdxConvert->Eval(cor_dQdx);
					if (reslength < 1) continue;
					if (cor_dQdx <= 0 || (raw_dQdx < 37.01 && raw_dQdx > 36.99)) continue;
					if (InFid(point.X(), point.Y(), point.Z()) ==  false) continue;
					hCorrectedDEdx_m->Fill(reslength, cor_dEdx);
					// if (reslength < 5 || reslength > 75) continue;
					}
				}


		}
	}

	void DEdx_Test::LoadTruthFile(std::string file_loc) {

		// TFile *fMCTruth = TFile::Open(Form("%s",file_loc.c_str()),"READ");
		// if(!(fMCTruth->IsOpen())){
		// 	std::cout << "ERROR : could not open selection file!" << std::endl;
	 //    }	
	 //    std::cout << "selection file opened" << std::endl;

	 //    TTree *MCTruth = (TTree*)fMCTruth->Get("EventMCINFO_DL");
	 //    if(MCTruth==nullptr)std::cout << "ERROR, no ttree loaded " << std::endl;
	 //    else std::cout << "ttree loaded with " <<MCTruth->GetEntries() << " entries" << std::endl;

	 //    int run_sel, subrun_sel, event_sel, entry_sel; 
	 //    double parentSCEX_sel, parentSCEY_sel, parentSCEZ_sel;

	 //    std::vector<double> *energyk_sel {nullptr};
	 //    std::vector<double> *pdg_sel {nullptr};

	 //    std::vector<int> 	*daughterPdg_v {nullptr};
	 //    std::vector<double> *daughterX_v {nullptr};
	 //    std::vector<double> *daughterY_v {nullptr};
	 //    std::vector<double> *daughterZ_v {nullptr};
	 //    std::vector<double> *daughterX_end_v {nullptr};
	 //    std::vector<double> *daghterY_end_v {nullptr};
	 //    std::vector<double> *daughterZ_end_v {nullptr};


	 //    MCTruth->SetBranchAddress("parentSCEX",&parentSCEX_sel);
	 //    MCTruth->SetBranchAddress("parentSCEY",&parentSCEY_sel);
	 //    MCTruth->SetBranchAddress("parentSCEZ",&parentSCEZ_sel);
	 //    MCTruth->SetBranchAddress("run",&run_sel);
	 //    MCTruth->SetBranchAddress("subrun",&subrun_sel);
	 //    MCTruth->SetBranchAddress("event",&event_sel); 
	 //    MCTruth->SetBranchAddress("daughter_energyk_v",&energyk_sel);
	 //    MCTruth->SetBranchAddress("daughterPdg_v",&pdg_sel);

	 //    MCTruth->SetBranchAddress("daughterX_v",&daughterX_v);
	 //    MCTruth->SetBranchAddress("daughterY_v",&daughterY_v);
	 //    MCTruth->SetBranchAddress("daughterZ_v",&daughterZ_v);
	 //    MCTruth->SetBranchAddress("daughterX_end_v",&daughterX_end_v);
	 //    MCTruth->SetBranchAddress("daughterY_end_v",&daughterY_end_v);
	 //    MCTruth->SetBranchAddress("daughterZ_end_v",&daughterZ_end_v);

	 //    // need to add in some vertex_id counter so I can match with the other variables properly.  Look at vertex_id above


	 //    for (Long64_t i {0}; i < MCTruth->GetEntries(); i++) {
	 //    	MCTruth->GetEntry(i);
	 //    	std::vector<std::vector<double>> vertex;
	 //    	std::string key = std::to_string(run_sel) + "_" + std::to_string(subrun_sel) + "_" + std::to_string(event_sel); // + "_" + std::to_string(vtx_id_sel)
	 //    	vertex.push_back(*pdg_sel);
	 //    	vertex.push_back(*daughterX_v);
	 //    	vertex.push_back(*daughterY_v);
	 //    	vertex.push_back(*daughterZ_v);
	 //    	vertex.push_back(*daughterX_end_v);
	 //    	vertex.push_back(*daughterY_end_v);
	 //    	vertex.push_back(*daughterZ_end_v);
	 //    	Truth_Map[key] = vertex;
	 //    }
	 //    for (Long64_t i {0}; i < MCTruth->GetEntries(); i++) {
	 //    	MCTruth->GetEntry(i);
	 //    	std::vector<double> info;
	 //    	std::vector<std::vector<double>> info2;
	 //    	info.push_back(run_sel);
	 //    	info.push_back(subrun_sel);
	 //    	info.push_back(event_sel);
	 //    	info.push_back(parentSCEX_sel);
	 //    	info.push_back(parentSCEY_sel);
	 //    	info.push_back(parentSCEZ_sel);
	 //    	info2.push_back(*pdg_sel);
	 //    	info2.push_back(*energyk_sel);

	 //    	for (size_t idx {0}; idx < SelectEvtID_proton.size(); idx++) {
	 //    		if (info[0] == SelectEvtID_proton[idx][0] && info[1] == SelectEvtID_proton[idx][1] && info[2] == SelectEvtID_proton[idx][2]) {
	 //    			SelectTruthID_proton.push_back(info);
	 //    			SelectTruthID_energy.push_back(info2);
	 //    			//std::cout << "Found Entry" << std::endl;
	 //    		}
	 //    	}
	 //    }

	    //for (Long64_t i = 0; i < NuMuVertexVariable->GetEntries();i++) { 
	    	//NuMuVertexVariable->GetEntry(i);
	    	//std::cout << "pass" << std::endl;
	    	//if (InFiducial_sel != 1 || goodReco_sel == 0 || N5cmTracks_sel != 2 || TrackLength_v_sel->size() != 2) continue;
	    	//std::cout << InFiducial_sel << std::endl;
	    	//if (TrackLength_v_sel->at(0) > 35 && TrackLength_v_sel->at(1) > 35 
	    		//&& TrackLength_v_sel->at(0)+TrackLength_v_sel->at(1) > 200
	    		//&& cosmicLL_sel < 3) {

	    		//std::vector<int> info;
	    		//info.push_back(run_sel);
	    		//info.push_back(subrun_sel);
	    		//info.push_back(event_sel);
	    		//info.push_back(vtxid_sel);
	    		//SelectEvtID_cosmic.push_back(info);
	    	//} 
	    	//else if (cosmicLL_sel > 4 && PassCuts_sel == 1) {
	    		//std::vector<int> info;
	    		//info.push_back(run_sel);
	    		//info.push_back(subrun_sel);
	    		//info.push_back(event_sel);
	    		//info.push_back(vtxid_sel);
	    		//SelectEvtID_proton.push_back(info);

	    	//}
	    //}
    	//std::cout << SelectEvtID_cosmic.size() << " & " << SelectEvtID_proton.size() << std::endl;
        //std::cout << "tree looped through" << std::endl;
        // fMCTruth->Close();
	
	}

	void DEdx_Test::MakeRecoFill(){
		// Open into read/write mode
		TFile *f = new TFile("RecoFile.root","UPDATE");
		// TTree *T = new TTree("Run1","Run1");
		std::cout << "Opening TTree" << std::endl;
		TTree *T = (TTree*) f->Get("Run1");

		std::cout << "Setting Variables" << std::endl;
		double alpha_cor, alpha_un;
		double beta_cor, beta_un;
		double scale_cor, scale_un;
		double run_start, run_final;
		TH3D *CalMap {nullptr};

		std::cout << "Setting Branches" << std::endl;
		T->Branch("alpha_un",&alpha_un);
		std::cout << "Setting 2nd Branches" << std::endl;

		T->SetBranchAddress("alpha_cor",&alpha_cor);
		T->SetBranchAddress("beta_un",&beta_un);
		T->SetBranchAddress("beta_cor",&beta_cor);
		T->SetBranchAddress("scale_cor",&scale_cor);
		T->SetBranchAddress("scale_un",&scale_un);
		T->SetBranchAddress("run_start",&run_start);
		T->SetBranchAddress("run_final",&run_final);
		std::cout << "Setting CalMapBranch" << std::endl;
		T->SetBranchAddress("CalMap",&CalMap);

		std::cout << "Filling Branches" << std::endl;
		alpha_un = alpha_fit_un;
		alpha_cor = alpha_fit;
		beta_un = beta_fit_un;
		beta_cor = beta_fit;
		scale_cor = avgVal_cor;
		scale_un = avgVal;
		run_start = first_run;
		run_final = last_run;
		std::cout << "Filling Histogram" << std::endl;
		CalMap = (TH3D*) CorrectionMap3D->Clone();

		std::cout << "Filling Tree" << std::endl;
		T->Fill();
		std::cout << "Writing Tree" << std::endl;
		T->Write();



	}

	void DEdx_Test::CreateDEDXMap(std::string file_loc) {

		TFile *fMCTruth = TFile::Open(Form("%s",file_loc.c_str()),"READ");
		if(!(fMCTruth->IsOpen())){ std::cout << "Error: could not open selection file!" << std::endl; }
		else { std::cout << "selection file opened" << std::endl; }

		TTree *MCTruth = (TTree*) fMCTruth->Get("dEdxTree");
		if (MCTruth==nullptr)std::cout << "Error, no TTree loaded" << std::endl;
		else std::cout << "TTree loaded with " << MCTruth->GetEntries() << std::endl;

		int run_sel, subrun_sel, event_sel, vtx_id_sel;

		std::vector<double> *Avg_dEdx {nullptr};
		std::vector<double> *Avg_dEdx_5cm {nullptr};
		std::vector<double> *Avg_dEdx_10cm {nullptr};
		std::vector<double> *Avg_dEdx_un {nullptr};
		std::vector<double> *Avg_dEdx_5cm_un {nullptr};
		std::vector<double> *Avg_dEdx_10cm_un {nullptr};

		std::vector<double> *Avg_corrdQdx {nullptr};
		std::vector<double> *Avg_corrdQdx_5cm {nullptr};
		std::vector<double> *Avg_corrdQdx_10cm {nullptr};
		std::vector<double> *Avg_corrdQdx_un {nullptr};
		std::vector<double> *Avg_corrdQdx_5cm_un {nullptr};
		std::vector<double> *Avg_corrdQdx_10cm_un {nullptr};

		std::vector<double> *chi2_m_hypothesis {nullptr};
		std::vector<double> *chi2_p_hypothesis {nullptr};
		std::vector<double> *chi2_invert_p_hypothesis {nullptr};
		std::vector<double> *chi2_invert_m_hypothesis {nullptr};
		std::vector<double> *chi2_m_hypothesis_un {nullptr};
		std::vector<double> *chi2_p_hypothesis_un {nullptr};
		std::vector<double> *chi2_invert_p_hypothesis_un {nullptr};
		std::vector<double> *chi2_invert_m_hypothesis_un {nullptr};

		MCTruth->SetBranchAddress("Avg_dEdx",&Avg_dEdx);
		MCTruth->SetBranchAddress("Avg_dEdx_5cm",&Avg_dEdx_5cm);
		MCTruth->SetBranchAddress("Avg_dEdx_10cm",&Avg_dEdx_10cm);
		MCTruth->SetBranchAddress("Avg_dEdx_un",&Avg_dEdx_un);
		MCTruth->SetBranchAddress("Avg_dEdx_5cm_un",&Avg_dEdx_5cm_un);
		MCTruth->SetBranchAddress("Avg_dEdx_10cm_un",&Avg_dEdx_10cm_un);

		MCTruth->SetBranchAddress("Avg_corrdQdx",&Avg_corrdQdx);
		MCTruth->SetBranchAddress("Avg_corrdQdx_5cm",&Avg_corrdQdx_5cm);
		MCTruth->SetBranchAddress("Avg_corrdQdx_10cm",&Avg_corrdQdx_10cm);
		MCTruth->SetBranchAddress("Avg_corrdQdx_un",&Avg_corrdQdx_un);
		MCTruth->SetBranchAddress("Avg_corrdQdx_5cm_un",&Avg_corrdQdx_5cm_un);
		MCTruth->SetBranchAddress("Avg_corrdQdx_10cm_un",&Avg_corrdQdx_10cm_un);

		MCTruth->SetBranchAddress("chi2_m_hypothesis",&chi2_m_hypothesis);
		MCTruth->SetBranchAddress("chi2_p_hypothesis",&chi2_p_hypothesis);
		MCTruth->SetBranchAddress("chi2_invert_p_hypothesis",&chi2_invert_p_hypothesis);
		MCTruth->SetBranchAddress("chi2_invert_m_hypothesis",&chi2_invert_m_hypothesis);
		MCTruth->SetBranchAddress("chi2_m_hypothesis_un",&chi2_m_hypothesis_un);
		MCTruth->SetBranchAddress("chi2_p_hypothesis_un",&chi2_p_hypothesis_un);
		MCTruth->SetBranchAddress("chi2_invert_p_hypothesis_un",&chi2_invert_p_hypothesis_un);
		MCTruth->SetBranchAddress("chi2_invert_m_hypothesis_un",&chi2_invert_m_hypothesis_un);

		MCTruth->SetBranchAddress("run",&run_sel);
		MCTruth->SetBranchAddress("subrun",&subrun_sel);
		MCTruth->SetBranchAddress("event",&event_sel);
		MCTruth->SetBranchAddress("vtx_id",&vtx_id_sel);

		for (Long64_t i {0}; i < MCTruth->GetEntries(); i++) {
			MCTruth->GetEntry(i);
			std::vector<std::vector<double>> vertex;
			std::string key = std::to_string(run_sel) + "_" + std::to_string(subrun_sel) + "_" + std::to_string(event_sel) + "_" + std::to_string(vtx_id_sel);
			vertex.push_back(*Avg_dEdx);
			vertex.push_back(*Avg_dEdx_5cm);
			vertex.push_back(*Avg_dEdx_10cm);
			vertex.push_back(*Avg_corrdQdx);
			vertex.push_back(*Avg_corrdQdx_5cm);
			vertex.push_back(*Avg_corrdQdx_10cm);
			vertex.push_back(*chi2_m_hypothesis);
			vertex.push_back(*chi2_p_hypothesis);
			vertex.push_back(*chi2_invert_m_hypothesis);
			vertex.push_back(*chi2_invert_p_hypothesis);
			vertex.push_back(*Avg_dEdx_un);
			vertex.push_back(*Avg_dEdx_5cm_un);
			vertex.push_back(*Avg_dEdx_10cm_un);
			vertex.push_back(*Avg_corrdQdx_un);
			vertex.push_back(*Avg_corrdQdx_5cm_un);
			vertex.push_back(*Avg_corrdQdx_10cm_un);
			vertex.push_back(*chi2_m_hypothesis_un);
			vertex.push_back(*chi2_p_hypothesis_un);
			vertex.push_back(*chi2_invert_m_hypothesis_un);
			vertex.push_back(*chi2_invert_p_hypothesis_un);
			dEdx_Map[key] = vertex;
		}
		std::cout << "dEdx Map Made" << std::endl;
	}

	void DEdx_Test::DEdxPlots() {
		hAvg_dEdx_p = new TH1D("hAvg_dEdx_p","hAvg_dEdx_p", 20, 0, 11);
		hAvg_dEdx_p_5cm = new TH1D("hAvg_dEdx_p_5cm","hAvg_dEdx_p_5cm", 20, 0, 11);
		hAvg_dEdx_p_10cm = new TH1D("hAvg_dEdx_p_10cm","hAvg_dEdx_p_10cm", 20, 0, 11);
		hAvg_dEdx_p_un = new TH1D("hAvg_dEdx_p_un","hAvg_dEdx_p_un", 20, 0, 11);
		hAvg_dEdx_p_5cm_un = new TH1D("hAvg_dEdx_p_5cm_un","hAvg_dEdx_p_5cm_un", 20, 0, 11);
		hAvg_dEdx_p_10cm_un = new TH1D("hAvg_dEdx_p_10cm_un","hAvg_dEdx_p_10cm_un", 20, 0, 11);

		hAvg_dEdx_m = new TH1D("hAvg_dEdx_m","hAvg_dEdx_m", 20, 0, 11);
		hAvg_dEdx_m_5cm = new TH1D("hAvg_dEdx_m_5cm","hAvg_dEdx_m_5cm", 20, 0, 11);
		hAvg_dEdx_m_10cm = new TH1D("hAvg_dEdx_m_10cm","hAvg_dEdx_m_10cm", 20, 0, 11);
		hAvg_dEdx_m_un = new TH1D("hAvg_dEdx_m_un","hAvg_dEdx_m_un", 20, 0, 11);
		hAvg_dEdx_m_5cm_un = new TH1D("hAvg_dEdx_m_5cm_un","hAvg_dEdx_m_5cm_un", 20, 0, 11);
		hAvg_dEdx_m_10cm_un = new TH1D("hAvg_dEdx_m_10cm_un","hAvg_dEdx_m_10cm_un", 20, 0, 11);

		hAvg_corrdQdx_p = new TH1D("hAvg_corrdQdx_p","hAvg_corrdQdx_p", 40, 0, 25);
		hAvg_corrdQdx_p_5cm = new TH1D("hAvg_corrdQdx_p_5cm","hAvg_corrdQdx_p_5cm", 40, 0, 25);
		hAvg_corrdQdx_p_10cm = new TH1D("hAvg_corrdQdx_p_10cm","hAvg_corrdQdx_p_10cm", 40, 0, 25);
		hAvg_corrdQdx_p_un = new TH1D("hAvg_corrdQdx_p_un","hAvg_corrdQdx_p_un", 40, 0, 25);
		hAvg_corrdQdx_p_5cm_un = new TH1D("hAvg_corrdQdx_p_5cm_un","hAvg_corrdQdx_p_5cm_un", 40, 0, 25);
		hAvg_corrdQdx_p_10cm_un = new TH1D("hAvg_corrdQdx_p_10cm_un","hAvg_corrdQdx_p_10cm_un", 40, 0, 25);

		hAvg_corrdQdx_m = new TH1D("hAvg_corrdQdx_m","hAvg_corrdQdx_m", 40, 0, 25);
		hAvg_corrdQdx_m_5cm = new TH1D("hAvg_corrdQdx_m_5cm","hAvg_corrdQdx_m_5cm", 40, 0, 25);
		hAvg_corrdQdx_m_10cm = new TH1D("hAvg_corrdQdx_m_10cm","hAvg_corrdQdx_m_10cm", 40, 0, 25);
		hAvg_corrdQdx_m_un = new TH1D("hAvg_corrdQdx_m_un","hAvg_corrdQdx_m_un", 40, 0, 25);
		hAvg_corrdQdx_m_5cm_un = new TH1D("hAvg_corrdQdx_m_5cm_un","hAvg_corrdQdx_m_5cm_un", 40, 0, 25);
		hAvg_corrdQdx_m_10cm_un = new TH1D("hAvg_corrdQdx_m_10cm_un","hAvg_corrdQdx_m_10cm_un", 40, 0, 25);

		hChi2_m_m = new TH1D("hChi2_m","hChi2_m",50, -1, 30);
		hChi2_m_p = new TH1D("hChi2_m","hChi2_m",50, -1, 30);
		hChi2_p_p = new TH1D("hChi2_p","hChi2_p",50,-1,30);
		hChi2_p_m = new TH1D("hChi2_p","hChi2_p",50,-1,30);

		hChi2_Inv_m_m = new TH1D("hChi2_Inv_m_m","hChi2_Inv_m_m",50,-1,30);
		hChi2_Inv_m_p = new TH1D("hChi2_Inv_m_p","hChi2_Inv_m_p",50,-1,30);
		hChi2_Inv_p_m = new TH1D("hChi2_Inv_p_m","hChi2_Inv_p_m",50,-1,30);
		hChi2_Inv_p_p = new TH1D("hChi2_Inv_p_p","hChi2_Inv_p_p",50,-1,30);

		hChi2_2D_m = new TH2D("hChi2_2D_m","hChi2_2D_m",50,-1,30,50,-1,30);
		hChi2_2D_p = new TH2D("hChi2_2D_p","hChi2_2D_p",50,-1,30,50,-1,30);
		

		for ( size_t ivert = 0; ivert < AllVertex_p.size(); ivert++) {
			double avg_dQdx[2] {0,0};
			double avg_dQdx_10cm[2] {0,0};
			double v_x = AllVertex_p[ivert][0].Vertex().X();
			double v_y = AllVertex_p[ivert][0].Vertex().Y();
			double v_z = AllVertex_p[ivert][0].Vertex().Z();
			//double dr = deltaR(v_x, v_y, v_z, AllTruth_p[ivert]);
			bool InFiducial = InFid(v_x, v_y, v_z);
			double dr = 1;
			// hdr->Fill(dr);
			// hEnergy->Fill(AllTruth_p[ivert][3]);





			if (dr > 5 || InFiducial == 0) continue;
				for ( size_t itrack = 0; itrack < AllVertex_p[ivert].size(); itrack++) {
					int points {0};
					int points_10cm {0};
					auto curr_vert = AllVertex_p[ivert][itrack];
					for (size_t iloc = 0; iloc < curr_vert.NumberTrajectoryPoints(); iloc++) {
						double raw_dQdx = curr_vert.DQdxAtPoint(iloc,plane);
						TVector3 point = curr_vert.LocationAtPoint(iloc);
						double reslength = curr_vert.Length(iloc);
						if (reslength < 1) continue;
						if (raw_dQdx <= 0 || (raw_dQdx < 37.01 && raw_dQdx > 36.99)) continue;
						points += 1;
						avg_dQdx[itrack] += raw_dQdx;
						if (InFid(point.X(), point.Y(), point.Z()) ==  false) continue;
						if (reslength < 10) { 
							avg_dQdx_10cm[itrack] += raw_dQdx; 
							points_10cm += 1;
						}
					}
					avg_dQdx[itrack] /= points;
				}
				int kProton = 0;
				int kMuon = 1;
				if (AllVertex_p[ivert][kProton].Length() < 5 || AllVertex_p[ivert][kMuon].Length() < 5) continue;
				if (AllVertex_p[ivert][kProton].Length() < AllVertex_p[ivert][kMuon].Length()) {
					if (avg_dQdx[kProton] < avg_dQdx[kMuon]) continue;
					//std::cout << "filling dedx corr" << std::endl;
					if (dEdxInfo_p[ivert][0][kProton] < 7) {
						hAvg_dEdx_p->Fill(dEdxInfo_p[ivert][0][kProton]*1.05);
					}
					else {hAvg_dEdx_p->Fill(dEdxInfo_p[ivert][0][kProton]*0.95);}
					//hAvg_dEdx_p->Fill(dEdxInfo_p[ivert][0][kProton]);
					if (dEdxInfo_p[ivert][0][kMuon] < 2.5) {
						hAvg_dEdx_m->Fill(dEdxInfo_p[ivert][0][kMuon]*1.05);
					}
					else {hAvg_dEdx_m->Fill(dEdxInfo_p[ivert][0][kMuon]*0.95);}
					//hAvg_dEdx_m->Fill(dEdxInfo_p[ivert][0][kMuon]);
					if (dEdxInfo_p[ivert][1][kProton] < 9) {
						hAvg_dEdx_p_5cm->Fill(dEdxInfo_p[ivert][1][kProton]*1.05);
					}
					else {hAvg_dEdx_p_5cm->Fill(dEdxInfo_p[ivert][1][kProton]*0.95);}
					//hAvg_dEdx_p_5cm->Fill(dEdxInfo_p[ivert][1][kProton]);
					if (dEdxInfo_p[ivert][1][kMuon] < 4.5) {
						hAvg_dEdx_m_5cm->Fill(dEdxInfo_p[ivert][1][kMuon]*1.05);
					}
					else { hAvg_dEdx_m_5cm->Fill(dEdxInfo_p[ivert][1][kMuon]*0.95);}
					//hAvg_dEdx_m_5cm->Fill(dEdxInfo_p[ivert][1][kMuon]);
					if (dEdxInfo_p[ivert][2][kProton] < 8.5) {
						hAvg_dEdx_p_10cm->Fill(dEdxInfo_p[ivert][2][kProton]*1.05);
					}
					else { hAvg_dEdx_p_10cm->Fill(dEdxInfo_p[ivert][2][kProton]*0.95); }
					//hAvg_dEdx_p_10cm->Fill(dEdxInfo_p[ivert][2][kProton]);
					if ( dEdxInfo_p[ivert][2][kMuon] < 4.25) {
						hAvg_dEdx_m_10cm->Fill(dEdxInfo_p[ivert][2][kMuon]*1.05);
					}
					else {hAvg_dEdx_m_10cm->Fill(dEdxInfo_p[ivert][2][kMuon]*0.95);}
					//hAvg_dEdx_m_10cm->Fill(dEdxInfo_p[ivert][2][kMuon]);
					//std::cout << "filling dqdx corr" << std::endl;
					hAvg_corrdQdx_p->Fill(dEdxInfo_p[ivert][3][kProton]);
					hAvg_corrdQdx_m->Fill(dEdxInfo_p[ivert][3][kMuon]);
					hAvg_corrdQdx_p_5cm->Fill(dEdxInfo_p[ivert][4][kProton]);
					hAvg_corrdQdx_m_5cm->Fill(dEdxInfo_p[ivert][4][kMuon]);
					hAvg_corrdQdx_p_10cm->Fill(dEdxInfo_p[ivert][5][kProton]);
					hAvg_corrdQdx_m_10cm->Fill(dEdxInfo_p[ivert][5][kMuon]);
					//std::cout << "filling chi2 corr" << std::endl;
					hChi2_m_m->Fill(dEdxInfo_p[ivert][6][kMuon]);
					hChi2_m_p->Fill(dEdxInfo_p[ivert][6][kProton]);
					hChi2_p_p->Fill(dEdxInfo_p[ivert][7][kProton]);
					hChi2_p_m->Fill(dEdxInfo_p[ivert][7][kMuon]);
					hChi2_Inv_m_m->Fill(dEdxInfo_p[ivert][9][kMuon]);
					hChi2_Inv_m_p->Fill(dEdxInfo_p[ivert][9][kProton]);
					hChi2_Inv_p_m->Fill(dEdxInfo_p[ivert][8][kMuon]);
					hChi2_Inv_p_p->Fill(dEdxInfo_p[ivert][8][kProton]);
					//std::cout << "filling dedx un" << std::endl;
					hAvg_dEdx_p_un->Fill(dEdxInfo_p[ivert][10][kProton]);
					hAvg_dEdx_m_un->Fill(dEdxInfo_p[ivert][10][kMuon]);
					hAvg_dEdx_p_5cm_un->Fill(dEdxInfo_p[ivert][11][kProton]);
					hAvg_dEdx_m_5cm_un->Fill(dEdxInfo_p[ivert][11][kMuon]);
					hAvg_dEdx_p_10cm_un->Fill(dEdxInfo_p[ivert][12][kProton]);
					hAvg_dEdx_m_10cm_un->Fill(dEdxInfo_p[ivert][12][kMuon]);
					// hAvg_corrdQdx_p_un->Fill(dEdxInfo_p[ivert][13][kProton]);
					// hAvg_corrdQdx_m_un->Fill(dEdxInfo_p[ivert][13][kMuon]);
					// hAvg_corrdQdx_p_5cm_un->Fill(dEdxInfo_p[ivert][14][kProton]);
					// hAvg_corrdQdx_m_5cm_un->Fill(dEdxInfo_p[ivert][14][kMuon]);
					// hAvg_corrdQdx_p_10cm_un->Fill(dEdxInfo_p[ivert][15][kProton]);
					// hAvg_corrdQdx_m_10cm_un->Fill(dEdxInfo_p[ivert][15][kMuon]);
					// hChi2_2D_p->Fill(dEdxInfo_p[ivert][6][kProton],dEdxInfo_p[ivert][7][kProton]);
					// hChi2_2D_m->Fill(dEdxInfo_p[ivert][6][kMuon],dEdxInfo_p[ivert][7][kMuon]);


					}

				kProton = 1;
				kMuon = 0;
				if (AllVertex_p[ivert][kProton].Length() < AllVertex_p[ivert][kMuon].Length()) {
					if (avg_dQdx[kProton] < avg_dQdx[kMuon]) continue;
					hAvg_dEdx_p->Fill(dEdxInfo_p[ivert][0][kProton]);
					hAvg_dEdx_m->Fill(dEdxInfo_p[ivert][0][kMuon]);
					hAvg_dEdx_p_5cm->Fill(dEdxInfo_p[ivert][1][kProton]);
					hAvg_dEdx_m_5cm->Fill(dEdxInfo_p[ivert][1][kMuon]);
					hAvg_dEdx_p_10cm->Fill(dEdxInfo_p[ivert][2][kProton]);
					hAvg_dEdx_m_10cm->Fill(dEdxInfo_p[ivert][2][kMuon]);
					hAvg_corrdQdx_p->Fill(dEdxInfo_p[ivert][3][kProton]);
					hAvg_corrdQdx_m->Fill(dEdxInfo_p[ivert][3][kMuon]);
					hAvg_corrdQdx_p_5cm->Fill(dEdxInfo_p[ivert][4][kProton]);
					hAvg_corrdQdx_m_5cm->Fill(dEdxInfo_p[ivert][4][kMuon]);
					hAvg_corrdQdx_p_10cm->Fill(dEdxInfo_p[ivert][5][kProton]);
					hAvg_corrdQdx_m_10cm->Fill(dEdxInfo_p[ivert][5][kMuon]);
					hChi2_m_m->Fill(dEdxInfo_p[ivert][6][kMuon]);
					hChi2_m_p->Fill(dEdxInfo_p[ivert][6][kProton]);
					hChi2_p_p->Fill(dEdxInfo_p[ivert][7][kProton]);
					hChi2_p_m->Fill(dEdxInfo_p[ivert][7][kMuon]);
					hChi2_Inv_m_m->Fill(dEdxInfo_p[ivert][9][kMuon]);
					hChi2_Inv_m_p->Fill(dEdxInfo_p[ivert][9][kProton]);
					hChi2_Inv_p_m->Fill(dEdxInfo_p[ivert][8][kMuon]);
					hChi2_Inv_p_p->Fill(dEdxInfo_p[ivert][8][kProton]);	
					hAvg_dEdx_p_un->Fill(dEdxInfo_p[ivert][10][kProton]);
					hAvg_dEdx_m_un->Fill(dEdxInfo_p[ivert][10][kMuon]);
					hAvg_dEdx_p_5cm_un->Fill(dEdxInfo_p[ivert][11][kProton]);
					hAvg_dEdx_m_5cm_un->Fill(dEdxInfo_p[ivert][11][kMuon]);
					hAvg_dEdx_p_10cm_un->Fill(dEdxInfo_p[ivert][12][kProton]);
					hAvg_dEdx_m_10cm_un->Fill(dEdxInfo_p[ivert][12][kMuon]);
					// hAvg_corrdQdx_p_un->Fill(dEdxInfo_p[ivert][13][kProton]);
					// hAvg_corrdQdx_m_un->Fill(dEdxInfo_p[ivert][13][kMuon]);
					// hAvg_corrdQdx_p_5cm_un->Fill(dEdxInfo_p[ivert][14][kProton]);
					// hAvg_corrdQdx_m_5cm_un->Fill(dEdxInfo_p[ivert][14][kMuon]);
					// hAvg_corrdQdx_p_10cm_un->Fill(dEdxInfo_p[ivert][15][kProton]);
					// hAvg_corrdQdx_m_10cm_un->Fill(dEdxInfo_p[ivert][15][kMuon]);	
					// hChi2_2D_p->Fill(dEdxInfo_p[ivert][6][kProton],dEdxInfo_p[ivert][7][kProton]);
					// hChi2_2D_m->Fill(dEdxInfo_p[ivert][6][kMuon],dEdxInfo_p[ivert][7][kMuon]);			
				}

			
		}

		TCanvas *c1 = new TCanvas("c1","c1",200,10,600,1000);
		c1->Divide(1,3);
		c1->cd(1);
		hAvg_dEdx_m->SetLineColor(4);
		hAvg_dEdx_m->SetLineWidth(3);
		hAvg_dEdx_m->Draw();
		hAvg_dEdx_p->SetLineColor(2);
		hAvg_dEdx_p->SetLineWidth(3);
		hAvg_dEdx_p->Draw("same");
		hAvg_dEdx_m_un->SetLineColor(4);
		hAvg_dEdx_m_un->SetLineWidth(1);
		hAvg_dEdx_m_un->SetLineStyle(9);
		hAvg_dEdx_m_un->Draw("same");
		hAvg_dEdx_p_un->SetLineColor(2);
		hAvg_dEdx_p_un->SetLineWidth(1);
		hAvg_dEdx_p_un->SetLineStyle(9);
		hAvg_dEdx_p_un->Draw("same");

		c1->cd(2); 
		hAvg_dEdx_m_5cm->SetLineColor(4);
		hAvg_dEdx_m_5cm->SetLineWidth(3);
		hAvg_dEdx_m_5cm->Draw();
		hAvg_dEdx_p_5cm->SetLineColor(2);
		hAvg_dEdx_p_5cm->SetLineWidth(3);
		hAvg_dEdx_p_5cm->Draw("same");
		hAvg_dEdx_m_5cm_un->SetLineColor(4);
		hAvg_dEdx_m_5cm_un->SetLineWidth(1);
		hAvg_dEdx_m_5cm_un->SetLineStyle(9);
		hAvg_dEdx_m_5cm_un->Draw("same");
		hAvg_dEdx_p_5cm_un->SetLineColor(2);
		hAvg_dEdx_p_5cm_un->SetLineWidth(1);
		hAvg_dEdx_p_5cm_un->SetLineStyle(9);
		hAvg_dEdx_p_5cm_un->Draw("same");

		c1->cd(3);
		hAvg_dEdx_p_10cm->Draw();
		hAvg_dEdx_m_10cm->SetLineColor(4);
		hAvg_dEdx_m_10cm->SetLineWidth(3);
		hAvg_dEdx_m_10cm->Draw();
		hAvg_dEdx_p_10cm->SetLineWidth(3);
		hAvg_dEdx_p_10cm->SetLineColor(2);
		hAvg_dEdx_p_10cm->Draw("same");
		hAvg_dEdx_m_10cm_un->SetLineColor(4);
		hAvg_dEdx_m_10cm_un->SetLineWidth(1);
		hAvg_dEdx_m_10cm_un->SetLineStyle(9);
		hAvg_dEdx_m_10cm_un->Draw("same");
		hAvg_dEdx_p_10cm_un->SetLineColor(2);
		hAvg_dEdx_p_10cm_un->SetLineWidth(1);
		hAvg_dEdx_p_10cm_un->SetLineStyle(9);
		hAvg_dEdx_p_10cm_un->Draw("same");

		TCanvas *c2 = new TCanvas("c2","c2",200,10,600,1000);
		c2->Divide(1,3);
		c2->cd(1);
		hAvg_corrdQdx_m->SetLineColor(4);
		hAvg_corrdQdx_p->SetLineColor(2);
		hAvg_corrdQdx_m->SetLineWidth(3);
		hAvg_corrdQdx_p->SetLineWidth(3);
		hAvg_corrdQdx_m->Draw();
		hAvg_corrdQdx_p->Draw("same");
		// hAvg_corrdQdx_m_un->SetLineColor(4);
		// hAvg_corrdQdx_m_un->SetLineWidth(1);
		// hAvg_corrdQdx_m_un->SetLineStyle(9);
		// hAvg_corrdQdx_m_un->Draw("same");
		// hAvg_corrdQdx_p_un->SetLineColor(2);
		// hAvg_corrdQdx_p_un->SetLineWidth(1);
		// hAvg_corrdQdx_p_un->SetLineStyle(9);
		// hAvg_corrdQdx_p_un->Draw("same");

		c2->cd(2);
		hAvg_corrdQdx_p_5cm->SetLineColor(2);
		hAvg_corrdQdx_m_5cm->SetLineColor(4);
		hAvg_corrdQdx_p_5cm->SetLineWidth(3);
		hAvg_corrdQdx_m_5cm->SetLineWidth(3);
		hAvg_corrdQdx_m_5cm->Draw();
		hAvg_corrdQdx_p_5cm->Draw("same");
		// hAvg_corrdQdx_m_5cm_un->SetLineColor(4);
		// hAvg_corrdQdx_m_5cm_un->SetLineWidth(1);
		// hAvg_corrdQdx_m_5cm_un->SetLineStyle(9);
		// hAvg_corrdQdx_m_5cm_un->Draw("same");
		// hAvg_corrdQdx_p_5cm_un->SetLineColor(2);
		// hAvg_corrdQdx_p_5cm_un->SetLineWidth(1);
		// hAvg_corrdQdx_p_5cm_un->SetLineStyle(9);
		// hAvg_corrdQdx_p_5cm_un->Draw("same");

		c2->cd(3);
		hAvg_corrdQdx_p_10cm->SetLineColor(2);
		hAvg_corrdQdx_m_10cm->SetLineColor(4);
		hAvg_corrdQdx_p_10cm->SetLineWidth(3);
		hAvg_corrdQdx_m_10cm->SetLineWidth(3);
		hAvg_corrdQdx_m_10cm->Draw();
		hAvg_corrdQdx_p_10cm->Draw("same");
		// hAvg_corrdQdx_m_10cm_un->SetLineColor(4);
		// hAvg_corrdQdx_m_10cm_un->SetLineWidth(1);
		// hAvg_corrdQdx_m_10cm_un->SetLineStyle(9);
		// hAvg_corrdQdx_m_10cm_un->Draw("same");
		// hAvg_corrdQdx_p_10cm_un->SetLineColor(2);
		// hAvg_corrdQdx_p_10cm_un->SetLineWidth(1);
		// hAvg_corrdQdx_p_10cm_un->SetLineStyle(9);
		// hAvg_corrdQdx_p_10cm_un->Draw("same");

		TCanvas *c3 = new TCanvas("c3","c3",200,10, 750, 750);
		c3->Divide(1,2);
		c3->cd(1);
		hChi2_m_m->SetLineWidth(3); hChi2_m_p->SetLineWidth(3);
		hChi2_m_m->SetLineColor(4); hChi2_m_m->Draw();
		hChi2_m_p->SetLineColor(2); hChi2_m_p->Draw("same");
		c3->cd(2);
		hChi2_p_m->SetLineWidth(3); hChi2_p_p->SetLineWidth(3);
		hChi2_p_m->SetLineColor(4); hChi2_p_m->Draw();
		hChi2_p_p->SetLineColor(2); hChi2_p_p->Draw("same");

		TCanvas *c4 = new TCanvas("c4","c4",200,10, 750, 750);
		c4->Divide(1,2);
		c4->cd(1);
		hChi2_Inv_m_m->SetLineWidth(3); hChi2_Inv_m_p->SetLineWidth(3);
		hChi2_Inv_m_m->SetLineColor(4); hChi2_Inv_m_m->Draw();
		hChi2_Inv_m_p->SetLineColor(2); hChi2_Inv_m_p->Draw("same");
		c4->cd(2);
		hChi2_Inv_p_m->SetLineWidth(3); hChi2_Inv_p_p->SetLineWidth(3);
		hChi2_Inv_p_m->SetLineColor(4); hChi2_Inv_p_m->Draw();
		hChi2_Inv_p_p->SetLineColor(2); hChi2_Inv_p_p->Draw("same");

		TCanvas *c5 = new TCanvas("c5","c5",200,10,400,400);
		c5->cd();
		hChi2_2D_m->SetMarkerStyle(3);
		hChi2_2D_m->SetMarkerColor(4);
		hChi2_2D_p->SetMarkerStyle(3);
		hChi2_2D_p->SetMarkerColor(2);
		hChi2_2D_m->Draw();
		hChi2_2D_p->Draw("same");
		
	}

	void DEdx_Test::Tracks_All() {
		const int ix {41};
		double reco_start {2.75}, reco_end{13};
		double add_len = (reco_end - reco_start) / ix;
		double bin_start = reco_start;
		double bin_end = reco_start + add_len;
		std::vector<double> RecoBins {};
		TH1D* hCorrRecoPlotBin_Y[ix];
		TH1D* hUnRecoPlotBin_Y[ix];
		TH1D* hCorrRecoPlotBin_All[ix];
		TH1D* hUnRecoPlotBin_All[ix];
		TH1D* hADC_Y_Un = new TH1D("hADC_Y_Un","hADC_Y_Un",100,0,250000);
		TH1D* hADC_Y_Corr = new TH1D("hADC_Y_Corr","hADC_Y_Corr",100,0,250000);
		TH1D* hADC_All_Un = new TH1D("hADC_All_Un","hADC_All_Un",100,0,250000);
		TH1D* hADC_All_Cor = new TH1D("hADC_All_Cor","hADC_Y_Cor",100,0,250000);

		for (int idx {0}; idx < ix; idx++) {
			RecoBins.push_back(bin_end);
			hCorrRecoPlotBin_Y[idx] = new TH1D(Form("hCorrRecoPlotBin_Y%02d",idx),Form("hCorrRecoPlotBin_Y%02d",idx),60,0,500000);
			hUnRecoPlotBin_Y[idx] = new TH1D(Form("hUnRecoPlotBin_Y%02d",idx),Form("hUnRecoPlotBin_Y%02d",idx),60,0,500000);
			hCorrRecoPlotBin_All[idx] = new TH1D(Form("hCorrRecoPlotBin_All%02d",idx),Form("hCorrRecoPlotBin_All%02d",idx),60,0,500000);
			hUnRecoPlotBin_All[idx] = new TH1D(Form("hUnRecoPlotBin_All%02d",idx),Form("hUnRecoPlotBin_All%02d",idx),60,0,500000);
			bin_start+=add_len;
			bin_end+=add_len;
		}

		TH2D* hUnDQdx[3];
		TH2D* hUnDQdx_p[3];
		TH2D* hUnDQdx_m[3];
		TH1D* hADC[3];
		TH1D* hADC_Corr[3];
		for (size_t iPlane {1}; iPlane < 4; iPlane++) {
			hUnDQdx[iPlane] = new TH2D(Form("hUnDQdx_%zu",iPlane),Form("hUnDQdx_%zu",iPlane),400,0,400,300, 0, 500000);
			hUnDQdx_p[iPlane] = new TH2D(Form("hUnDQdx_p_%zu",iPlane),Form("hUnDQdx_p_%zu",iPlane),400,0,400,300, 0, 500000);
			hUnDQdx_m[iPlane] = new TH2D(Form("hUnDQdx_m_%zu",iPlane),Form("hUnDQdx_m_%zu",iPlane),400,0,400,300, 0, 500000);
			hADC[iPlane] = new TH1D(Form("hADC_%zu",iPlane),Form("hADC_%zu",iPlane),200,0,30);
			hADC_Corr[iPlane] = new TH1D(Form("hADC_Corr_%zu",iPlane),Form("hADC_Corr_%zu",iPlane),200,0,30);
		}
		hUnDQdx_All = new TH2D("hUnDQdx_All","hUnDQdx_All", 400, 0, 400, 300, 0, 500000);
		hUnDQdx_All_p = new TH2D("hUnDQdx_All_p","hUnDQdx_All_p", 400, 0, 400, 300, 0, 500000);
		hUnDQdx_All_m = new TH2D("hUnDQdx_All_m","hUnDQdx_All_m", 400, 0, 400, 300, 0, 500000);
		hADC_All = new TH1D("hADC_All","hADC_All",200,0,30);
		hADC_Corr_All = new TH1D("hADC_Corr_All","hADC_Corr_All",200,0,30);

		TH2D* hCorrDQdx[3];
		TH2D* hCorrDQdx_p[3];
		TH2D* hCorrDQdx_m[3];
		for (size_t iPlane {1}; iPlane < 4; iPlane++) {
			hCorrDQdx[iPlane] = new TH2D(Form("hCorrDQdx_%zu",iPlane),Form("hCorrDQdx_%zu",iPlane),400,0,400,300, 0, 500000);
			hCorrDQdx_p[iPlane] = new TH2D(Form("hCorrDQdx_p_%zu",iPlane),Form("hCorrDQdx_p_%zu",iPlane),400,0,400,300, 0, 500000);
			hCorrDQdx_m[iPlane] = new TH2D(Form("hCorrDQdx_m_%zu",iPlane),Form("hCorrDQdx_m_%zu",iPlane),400,0,400,300, 0, 500000);
		}
		hCorrDQdx_All = new TH2D("hCorrDQdx_All","hCorrDQdx_All", 400, 0, 400, 300, 0, 500000);
		hCorrDQdx_All_p = new TH2D("hCorrDQdx_All_p","hCorrDQdx_All_p", 400, 0, 400, 300, 0, 500000);
		hCorrDQdx_All_m = new TH2D("hCorrDQdx_All_m","hCorrDQdx_All_m", 400, 0, 400, 300, 0, 500000);

		hUncorrectedRecombinationPlot_Y = new TH2D("hUncorrectedRecombinationPlot_Y","hUncorrectedRecombinationPlot_Y",60,0,15,50,0,500000);
		hCorrectedRecombinationPlot_Y = new TH2D("hCorrectedRecombinationPlot_Y","hCorrectedRecombinationPlot_Y",60,0,15,50,0,500000);
		hUncorrectedRecombinationPlot_All = new TH2D("hUncorrectedRecombinationPlot_All","hUncorrectedRecombinationPlot_All", 60, 0, 15, 50, 0, 500000);
		hCorrectedRecombinationPlot_All = new TH2D("hCorrectedRecombinationPlot_All","hCorrectedRecombinationPlot_All", 60, 0, 15, 50, 0, 500000);




		for ( size_t ivert = 0; ivert < AllVertex_p.size(); ivert++) {
			double avg_dQdx[2] {0,0};
			double avg_dQdx_10cm[2] {0,0};
			double v_x = AllVertex_p[ivert][0].Vertex().X();
			double v_y = AllVertex_p[ivert][0].Vertex().Y();
			double v_z = AllVertex_p[ivert][0].Vertex().Z();
			//double dr = deltaR(v_x, v_y, v_z, AllTruth_p[ivert]);
			bool InFiducial = InFid(v_x, v_y, v_z);
			double dr = 1;

			if (dr > 5 || InFiducial == 0) continue;
				for ( size_t itrack = 0; itrack < AllVertex_p[ivert].size(); itrack++) {
					int points {0};
					int points_10cm {0};
					auto curr_vert = AllVertex_p[ivert][itrack];
					for (size_t iloc = 0; iloc < curr_vert.NumberTrajectoryPoints(); iloc++) {
						double raw_dQdx_U = curr_vert.DQdxAtPoint(iloc,larlite::geo::kU);
						double raw_dQdx_V = curr_vert.DQdxAtPoint(iloc,larlite::geo::kV);
						double raw_dQdx_Y = curr_vert.DQdxAtPoint(iloc,larlite::geo::kZ);
						TVector3 point = curr_vert.LocationAtPoint(iloc);
						double reslength = curr_vert.Length(iloc);
						if (reslength < 1) continue;
						if (raw_dQdx_U <= 0 || raw_dQdx_V <= 0 || raw_dQdx_Y <= 0 || raw_dQdx_U == 37. || raw_dQdx_V == 37. || raw_dQdx_Y == 37.) continue;
						points += 1;
						avg_dQdx[itrack] += raw_dQdx_U + raw_dQdx_V + raw_dQdx_Y;
						if (InFid(point.X(), point.Y(), point.Z()) ==  false) continue;
						//hUncorrectedDQdx->Fill(reslength, raw_dQdx*avgVal);
						}
					avg_dQdx[itrack] /= points;
					}
				int kProton = 0;
				int kMuon = 1;
				if (AllVertex_p[ivert][kProton].Length() < 5 || AllVertex_p[ivert][kMuon].Length() < 5) continue;
				if (AllVertex_p[ivert][kProton].Length() < AllVertex_p[ivert][kMuon].Length()) {
					if (avg_dQdx[kProton] < avg_dQdx[kMuon]) continue;
					for (size_t iloc = 0; iloc < AllVertex_p[ivert][kProton].NumberTrajectoryPoints(); iloc++) {
						auto curr_vert = AllVertex_p[ivert][kProton];
						double raw_dQdx_U = curr_vert.DQdxAtPoint(iloc,larlite::geo::kU);
						double raw_dQdx_V = curr_vert.DQdxAtPoint(iloc,larlite::geo::kV);
						double raw_dQdx_Y = curr_vert.DQdxAtPoint(iloc,larlite::geo::kZ);
						double raw_dQdx = raw_dQdx_U + raw_dQdx_V + raw_dQdx_Y;
						TVector3 point = curr_vert.LocationAtPoint(iloc);
						int bin = CorrectionMap3D->FindBin(point.X(),point.Y(),point.Z());
						double corr_dQdx_U = raw_dQdx_U*CorrectionMap3D_U->GetBinContent(bin);
						double corr_dQdx_V = raw_dQdx_V*CorrectionMap3D_V->GetBinContent(bin);
						double corr_dQdx_Y = raw_dQdx_Y*CorrectionMap3D->GetBinContent(bin);
						double corr_dQdx = corr_dQdx_U + corr_dQdx_V + corr_dQdx_Y;
						double reslength = curr_vert.Length(iloc);
						if (reslength < 0.25) continue;
						if (raw_dQdx_U > 0 && (raw_dQdx_U > 37.01 || raw_dQdx_U < 36.99)) { 
							hUnDQdx_p[1]->Fill(reslength, raw_dQdx_U*avgVal_U);
							hCorrDQdx_p[1]->Fill(reslength, corr_dQdx_U*avgVal_U); 
						}
						if (raw_dQdx_V > 0 && (raw_dQdx_V > 37.01 || raw_dQdx_V < 36.99)) { 
							hUnDQdx_p[2]->Fill(reslength, raw_dQdx_V*avgVal_V); 
							hCorrDQdx_p[2]->Fill(reslength, corr_dQdx_V*avgVal_V);
						}
						if (raw_dQdx_Y > 0 && (raw_dQdx_Y > 37.01 || raw_dQdx_Y < 36.99)) { 
							hUnDQdx_p[3]->Fill(reslength, raw_dQdx_Y*avgVal_Y); 
							hCorrDQdx_p[3]->Fill(reslength, corr_dQdx_Y*avgVal_Y_Corr);
							hUncorrectedRecombinationPlot_Y->Fill(sProtonRange2dEdx->Eval(reslength),raw_dQdx_Y*avgVal_Y);
							hCorrectedRecombinationPlot_Y->Fill(sProtonRange2dEdx->Eval(reslength),corr_dQdx_Y*avgVal_Y_Corr);
							double dEdx_theory_Y = sProtonRange2dEdx->Eval(reslength);
							int index_Y {1};
							int bin_num_Y {1};
							while (dEdx_theory_Y > RecoBins[index_Y] && index_Y < 40) {
								bin_num_Y = index_Y; 
								index_Y++;
							}
							hCorrRecoPlotBin_Y[bin_num_Y]->Fill(corr_dQdx_Y*avgVal_Y_Corr);
							hUnRecoPlotBin_Y[bin_num_Y]->Fill(raw_dQdx_Y*avgVal_Y);
						}
						if (raw_dQdx_U <= 0 || raw_dQdx_V <= 0 || raw_dQdx_Y <= 0 || (raw_dQdx_U < 37.01 && raw_dQdx_U > 36.99) || (raw_dQdx_V < 37.01 && raw_dQdx_V > 36.99) || (raw_dQdx_Y < 37.01 && raw_dQdx_Y > 36.99)) continue;
						if (InFid(point.X(), point.Y(), point.Z()) ==  false) continue;
						hUnDQdx_All_p->Fill(reslength, raw_dQdx*avgVal_all/3);
						hCorrDQdx_All_p->Fill(reslength, corr_dQdx*avgVal_all_Corr/3);
						hUncorrectedRecombinationPlot_All->Fill(sProtonRange2dEdx->Eval(reslength),raw_dQdx*avgVal_all/3);
						hCorrectedRecombinationPlot_All->Fill(sProtonRange2dEdx->Eval(reslength),corr_dQdx*avgVal_all_Corr/3);
						double dEdx_theory = sProtonRange2dEdx->Eval(reslength);
						int index {1};
						int bin_num {1};
						while (dEdx_theory > RecoBins[index] && index < 40) {
							bin_num = index; 
							index++;
						}
						hCorrRecoPlotBin_All[bin_num]->Fill(corr_dQdx*avgVal_all_Corr/3);
						hUnRecoPlotBin_All[bin_num]->Fill(raw_dQdx*avgVal_all/3);				
					}
					for (size_t iloc = 0; iloc < AllVertex_p[ivert][kMuon].NumberTrajectoryPoints(); iloc++) {
						auto curr_vert = AllVertex_p[ivert][kMuon];
						double raw_dQdx_U = curr_vert.DQdxAtPoint(iloc,larlite::geo::kU);
						double raw_dQdx_V = curr_vert.DQdxAtPoint(iloc,larlite::geo::kV);
						double raw_dQdx_Y = curr_vert.DQdxAtPoint(iloc,larlite::geo::kZ);
						double raw_dQdx = raw_dQdx_U + raw_dQdx_V + raw_dQdx_Y;
						TVector3 point = curr_vert.LocationAtPoint(iloc);
						int bin = CorrectionMap3D->FindBin(point.X(),point.Y(),point.Z());
						double corr_dQdx_U = raw_dQdx_U*CorrectionMap3D_U->GetBinContent(bin);
						double corr_dQdx_V = raw_dQdx_V*CorrectionMap3D_V->GetBinContent(bin);
						double corr_dQdx_Y = raw_dQdx_Y*CorrectionMap3D->GetBinContent(bin);
						double corr_dQdx = corr_dQdx_U + corr_dQdx_V + corr_dQdx_Y;
						double reslength = curr_vert.Length(iloc);
						if (reslength < 0.25) continue;
						if (raw_dQdx_U > 0 && (raw_dQdx_U > 37.01 || raw_dQdx_U < 36.99)) { 
							hUnDQdx_m[1]->Fill(reslength, raw_dQdx_U*avgVal_U);
							hCorrDQdx_m[1]->Fill(reslength, corr_dQdx_U*avgVal_U);
							if (reslength > 20) {hADC[1]->Fill(raw_dQdx_U); hADC_Corr[1]->Fill(corr_dQdx_U); }
						}
						if (raw_dQdx_V > 0 && (raw_dQdx_V > 37.01 || raw_dQdx_V < 36.99)) { 
							hUnDQdx_m[2]->Fill(reslength, raw_dQdx_V*avgVal_V);
							hCorrDQdx_m[2]->Fill(reslength, corr_dQdx_V*avgVal_V);
							if (reslength > 20) {hADC[2]->Fill(raw_dQdx_V);  hADC_Corr[2]->Fill(corr_dQdx_V); }
						}
						if (raw_dQdx_Y > 0 && (raw_dQdx_Y > 37.01 || raw_dQdx_Y < 36.99)) { 
							hUnDQdx_m[3]->Fill(reslength, raw_dQdx_Y*avgVal_Y); 
							hCorrDQdx_m[3]->Fill(reslength, corr_dQdx_Y*avgVal_Y_Corr);
							if (reslength > 20) {hADC[3]->Fill(raw_dQdx_Y);  hADC_Corr[3]->Fill(corr_dQdx_Y); hADC_Y_Un->Fill(raw_dQdx_Y*avgVal_Y); hADC_Y_Corr->Fill(corr_dQdx_Y*avgVal_Y_Corr);}
						}
						if (raw_dQdx_U <= 0 || raw_dQdx_V <= 0 || raw_dQdx_Y <= 0 || raw_dQdx_U == 37. || raw_dQdx_V == 37. || raw_dQdx_Y == 37.) continue;
						if (InFid(point.X(), point.Y(), point.Z()) ==  false) continue;
						hUnDQdx_All_m->Fill(reslength, raw_dQdx*avgVal_all/3);
						hCorrDQdx_All_m->Fill(reslength, corr_dQdx*avgVal_all/3);
						if (reslength > 20) { hADC_All->Fill(raw_dQdx/3);  hADC_Corr_All->Fill(corr_dQdx/3); hADC_All_Un->Fill(raw_dQdx*avgVal_all/3); hADC_All_Cor->Fill(corr_dQdx*avgVal_all_Corr/3);}
					}
				}
				kProton = 1;
				kMuon = 0;
				if (AllVertex_p[ivert][kProton].Length() < AllVertex_p[ivert][kMuon].Length()) {
					if (avg_dQdx[kProton] < avg_dQdx[kMuon]) continue;
					for (size_t iloc = 0; iloc < AllVertex_p[ivert][kProton].NumberTrajectoryPoints(); iloc++) {
						auto curr_vert = AllVertex_p[ivert][kProton];
						double raw_dQdx_U = curr_vert.DQdxAtPoint(iloc,larlite::geo::kU);
						double raw_dQdx_V = curr_vert.DQdxAtPoint(iloc,larlite::geo::kV);
						double raw_dQdx_Y = curr_vert.DQdxAtPoint(iloc,larlite::geo::kZ);
						double raw_dQdx = raw_dQdx_U + raw_dQdx_V + raw_dQdx_Y;
						TVector3 point = curr_vert.LocationAtPoint(iloc);
						int bin = CorrectionMap3D->FindBin(point.X(),point.Y(),point.Z());
						double corr_dQdx_U = raw_dQdx_U*CorrectionMap3D_U->GetBinContent(bin);
						double corr_dQdx_V = raw_dQdx_V*CorrectionMap3D_V->GetBinContent(bin);
						double corr_dQdx_Y = raw_dQdx_Y*CorrectionMap3D->GetBinContent(bin);
						double corr_dQdx = corr_dQdx_U + corr_dQdx_V + corr_dQdx_Y;
						double reslength = curr_vert.Length(iloc);
						if (reslength < 0.25) continue;
						if (raw_dQdx_U > 0 && (raw_dQdx_U > 37.01 || raw_dQdx_U < 36.99)) { 
							hUnDQdx_p[1]->Fill(reslength, raw_dQdx_U*avgVal_U);
							hCorrDQdx_p[1]->Fill(reslength, corr_dQdx_U*avgVal_U); 
						}
						if (raw_dQdx_V > 0 && (raw_dQdx_V > 37.01 || raw_dQdx_V < 36.99)) { 
							hUnDQdx_p[2]->Fill(reslength, raw_dQdx_V*avgVal_V); 
							hCorrDQdx_p[2]->Fill(reslength, corr_dQdx_V*avgVal_V);
						}
						if (raw_dQdx_Y > 0 && (raw_dQdx_Y > 37.01 || raw_dQdx_Y < 36.99)) { 
							hUnDQdx_p[3]->Fill(reslength, raw_dQdx_Y*avgVal_Y); 
							hCorrDQdx_p[3]->Fill(reslength, corr_dQdx_Y*avgVal_Y_Corr);
							hCorrectedRecombinationPlot_Y->Fill(sProtonRange2dEdx->Eval(reslength),corr_dQdx_Y*avgVal_Y_Corr);
							hUncorrectedRecombinationPlot_Y->Fill(sProtonRange2dEdx->Eval(reslength),raw_dQdx_Y*avgVal_Y);
							double dEdx_theory_Y = sProtonRange2dEdx->Eval(reslength);
							int index_Y {1};
							int bin_num_Y {1};
							while (dEdx_theory_Y > RecoBins[index_Y] && index_Y < 40) {
								bin_num_Y = index_Y; 
								index_Y++;
							}
							hCorrRecoPlotBin_Y[bin_num_Y]->Fill(corr_dQdx_Y*avgVal_Y_Corr);
							hUnRecoPlotBin_Y[bin_num_Y]->Fill(raw_dQdx_Y*avgVal_Y);
						}
						if (raw_dQdx_U <= 0 || raw_dQdx_V <= 0 || raw_dQdx_Y <= 0 || (raw_dQdx_U < 37.01 && raw_dQdx_U > 36.99) || (raw_dQdx_V < 37.01 && raw_dQdx_V > 36.99) || (raw_dQdx_Y < 37.01 && raw_dQdx_Y > 36.99)) continue;
						if (InFid(point.X(), point.Y(), point.Z()) ==  false) continue;
						hUnDQdx_All_p->Fill(reslength, raw_dQdx*avgVal_all/3);
						hCorrDQdx_All_p->Fill(reslength, corr_dQdx*avgVal_all_Corr/3);
						hUncorrectedRecombinationPlot_All->Fill(sProtonRange2dEdx->Eval(reslength),raw_dQdx*avgVal_all/3);
						hCorrectedRecombinationPlot_All->Fill(sProtonRange2dEdx->Eval(reslength),corr_dQdx*avgVal_all_Corr/3);
						double dEdx_theory = sProtonRange2dEdx->Eval(reslength);
						int index {1};
						int bin_num {1};
						while (dEdx_theory > RecoBins[index] && index < 40) {
							bin_num = index; 
							index++;
						}
						hCorrRecoPlotBin_All[bin_num]->Fill(corr_dQdx*avgVal_all_Corr/3);
						hUnRecoPlotBin_All[bin_num]->Fill(raw_dQdx*avgVal_all/3);

					}
					for (size_t iloc = 0; iloc < AllVertex_p[ivert][kMuon].NumberTrajectoryPoints(); iloc++) {
						auto curr_vert = AllVertex_p[ivert][kMuon];
						double raw_dQdx_U = curr_vert.DQdxAtPoint(iloc,larlite::geo::kU);
						double raw_dQdx_V = curr_vert.DQdxAtPoint(iloc,larlite::geo::kV);
						double raw_dQdx_Y = curr_vert.DQdxAtPoint(iloc,larlite::geo::kZ);
						double raw_dQdx = raw_dQdx_U + raw_dQdx_V + raw_dQdx_Y;
						TVector3 point = curr_vert.LocationAtPoint(iloc);
						int bin = CorrectionMap3D->FindBin(point.X(),point.Y(),point.Z());
						double corr_dQdx_U = raw_dQdx_U*CorrectionMap3D_U->GetBinContent(bin);
						double corr_dQdx_V = raw_dQdx_V*CorrectionMap3D_V->GetBinContent(bin);
						double corr_dQdx_Y = raw_dQdx_Y*CorrectionMap3D->GetBinContent(bin);
						double corr_dQdx = corr_dQdx_U + corr_dQdx_V + corr_dQdx_Y;
						double reslength = curr_vert.Length(iloc);
						if (reslength < 0.25) continue;
						if (raw_dQdx_U > 0 && (raw_dQdx_U > 37.01 || raw_dQdx_U < 36.99)) { 
							hUnDQdx_m[1]->Fill(reslength, raw_dQdx_U*avgVal_U);
							hCorrDQdx_m[1]->Fill(reslength, corr_dQdx_U*avgVal_U);
							if (reslength > 20) { hADC[1]->Fill(raw_dQdx_U);  hADC_Corr[1]->Fill(corr_dQdx_U); }
						}
						if (raw_dQdx_V > 0 && (raw_dQdx_V > 37.01 || raw_dQdx_V < 36.99)) { 
							hUnDQdx_m[2]->Fill(reslength, raw_dQdx_V*avgVal_V);
							hCorrDQdx_m[2]->Fill(reslength, corr_dQdx_V*avgVal_V);
							if (reslength > 20) { hADC[2]->Fill(raw_dQdx_V);  hADC_Corr[2]->Fill(corr_dQdx_V); }
						}
						if (raw_dQdx_Y > 0 && (raw_dQdx_Y > 37.01 || raw_dQdx_Y < 36.99)) { 
							hUnDQdx_m[3]->Fill(reslength, raw_dQdx_Y*avgVal_Y); 
							hCorrDQdx_m[3]->Fill(reslength, corr_dQdx_Y*avgVal_Y_Corr);
							if (reslength > 20) {hADC[3]->Fill(raw_dQdx_Y);  hADC_Corr[3]->Fill(corr_dQdx_Y); hADC_Y_Un->Fill(raw_dQdx_Y*avgVal_Y); hADC_Y_Corr->Fill(corr_dQdx_Y*avgVal_Y_Corr);}
						}
						if (raw_dQdx_U <= 0 || raw_dQdx_V <= 0 || raw_dQdx_Y <= 0 || raw_dQdx_U == 37. || raw_dQdx_V == 37. || raw_dQdx_Y == 37.) continue;
						if (InFid(point.X(), point.Y(), point.Z()) ==  false) continue;
						hUnDQdx_All_m->Fill(reslength, raw_dQdx*avgVal_all/3);
						hCorrDQdx_All_m->Fill(reslength, corr_dQdx*avgVal_all/3);
						if (reslength > 20) { hADC_All->Fill(raw_dQdx/3);  hADC_Corr_All->Fill(corr_dQdx/3); hADC_All_Un->Fill(raw_dQdx*avgVal_all/3); hADC_All_Cor->Fill(corr_dQdx*avgVal_all_Corr/3);}
					}
				}
		}	
		fRecombinationExpected = new TF1("fRecombinationExpected",RecombinationFunction,0,1000,6);
        fRecombinationExpected->SetParameters(alpha_ex,beta_ex,1,Wion,LAr_Density,E_Field);
        fRecombinationExpected->SetParNames("alpha","beta","C","W_{ion}","#rho_{LAr}","E_{field}");
        fRecombinationExpected->FixParameter(0, alpha_ex);
        fRecombinationExpected->FixParameter(1, beta_ex);
        fRecombinationExpected->FixParameter(3, Wion);
        fRecombinationExpected->FixParameter(4,LAr_Density);
        fRecombinationExpected->FixParameter(5,E_Field);

		TCanvas *cRecoBins = new TCanvas("cRecoBins","cRecoBins",200,10,500,1000);
		cRecoBins->Divide(1,5);
		TF1 *f2[42];
		for (int idx{1}; idx < ix - 7; idx++) {
			double mean_Un_Y = hUnRecoPlotBin_Y[idx]->GetMean();
			double rms_Un_Y = hUnRecoPlotBin_Y[idx]->GetRMS();
			double mean_Corr_Y = hCorrRecoPlotBin_Y[idx]->GetMean();
			double rms_Corr_Y = hCorrRecoPlotBin_Y[idx]->GetRMS();
			double mean_Un_All = hUnRecoPlotBin_All[idx]->GetMean();
			double rms_Un_All = hUnRecoPlotBin_All[idx]->GetRMS();
			double mean_Corr_All = hCorrRecoPlotBin_All[idx]->GetMean();
			double rms_Corr_All = hCorrRecoPlotBin_All[idx]->GetRMS();
			TF1 *f1 = new TF1("f1","gaus",0,500000);
			f2[idx] = new TF1(Form("f2%02d",idx),"gaus",0,500000);
			TF1 *f3 = new TF1("f3","gaus",0,500000);
			TF1 *f4 = new TF1("f4","gaus",0,500000);
			hUnRecoPlotBin_Y[idx]->Fit("f1","","",mean_Un_Y-rms_Un_Y, mean_Un_Y+rms_Un_Y);
			hCorrRecoPlotBin_Y[idx]->Fit(Form("f2%02d",idx),"","",mean_Corr_Y-rms_Corr_Y, mean_Corr_Y+rms_Corr_Y);
			hUnRecoPlotBin_All[idx]->Fit("f3","","",mean_Un_All-rms_Un_All, mean_Un_All+rms_Un_All);
			hCorrRecoPlotBin_All[idx]->Fit("f4","","",mean_Corr_All-rms_Corr_All, mean_Corr_All+rms_Corr_All);

			double val = RecoBins[idx] - 0.25/2;

			double ArgoVal = fRecombinationExpected->Eval(val);
			double ArgoVal_min = ArgoVal - (ArgoVal*.1);
			double ArgoVal_max = ArgoVal + (ArgoVal*.1);

			if (f2[idx]->GetParameter(1) < ArgoVal_min || f2[idx]->GetParameter(1) > ArgoVal_max ) {continue;}
			ebin.push_back(val);
			ey.push_back(0.25);

			mu_Un_Y.push_back(f1->GetParameter(1));
			mu_error_Un_Y.push_back(f1->GetParError(1));
			sigma_Un_Y.push_back(f1->GetParameter(2));
			sigma_error_Un_Y.push_back(f1->GetParError(2));

			mu_Corr_Y.push_back(f2[idx]->GetParameter(1));
			mu_error_Corr_Y.push_back(f2[idx]->GetParError(1));
			sigma_Corr_Y.push_back(f2[idx]->GetParameter(2));
			sigma_error_Corr_Y.push_back(f2[idx]->GetParError(2));

			mu_Un_All.push_back(f3->GetParameter(1));
			mu_error_Un_All.push_back(f3->GetParError(1));
			sigma_Un_All.push_back(f3->GetParameter(2));
			sigma_error_Un_All.push_back(f3->GetParError(2));

			mu_Corr_All.push_back(f4->GetParameter(1));	
			mu_error_Corr_All.push_back(f4->GetParError(1));
			sigma_Corr_All.push_back(f4->GetParameter(2));
			sigma_error_Corr_All.push_back(f4->GetParError(2));

			
		}

		cRecoBins->cd(1); hCorrRecoPlotBin_Y[4]->Draw();
		cRecoBins->cd(2); hCorrRecoPlotBin_Y[8]->Draw();
		cRecoBins->cd(3); hCorrRecoPlotBin_Y[12]->Draw();
		cRecoBins->cd(4); hCorrRecoPlotBin_Y[16]->Draw();
		cRecoBins->cd(5); hCorrRecoPlotBin_Y[20]->Draw();

		double RawU = hADC[1]->GetBinCenter(hADC[1]->GetMaximumBin());
		double meanU = hADC[1]->GetMean();
		double CorrU = hADC_Corr[1]->GetBinCenter(hADC_Corr[1]->GetMaximumBin());
		double meanCorrU = hADC_Corr[1]->GetMean();
		double RawV = hADC[2]->GetBinCenter(hADC[2]->GetMaximumBin());
		double meanV = hADC[2]->GetMean();
		double CorrV = hADC_Corr[2]->GetBinCenter(hADC_Corr[2]->GetMaximumBin());
		double meanCorrV = hADC_Corr[2]->GetMean();
		double RawY = hADC[3]->GetBinCenter(hADC[3]->GetMaximumBin());
		double meanY = hADC[3]->GetMean();
		double CorrY = hADC_Corr[3]->GetBinCenter(hADC_Corr[3]->GetMaximumBin());
		double meanCorrY = hADC_Corr[3]->GetMean();
		double RawAll = (hADC_All->GetBinCenter(hADC_All->GetMaximumBin()));
		double meanAll = hADC_All->GetMean();
		double CorrAll = hADC_Corr_All->GetBinCenter(hADC_Corr_All->GetMaximumBin());
		double meanCorrAll = hADC_Corr_All->GetMean();

		for (size_t iPlane{0}; iPlane < 3; iPlane++) {
			double mean = hADC_Corr[iPlane+1]->GetMean();
			double rms = hADC_Corr[iPlane+1]->GetRMS();
			TF1 *f1 = new TF1("f1","gaus",0,30);
			hADC_Corr[iPlane+1]->Fit("f1","","",mean-rms, mean+rms);
			std::cout << "mean: " << f1->GetParameter(0) << " and error: " << f1->GetParError(0) << std::endl;
		}

		double mean = hADC_Corr_All->GetMean();
		double rms = hADC_Corr_All->GetRMS();
		TF1 *f1 = new TF1("f1","gaus",0,30);
		hADC_Corr_All->Fit("f1","","",mean-rms, mean+rms);
		std::cout << "mean: " << f1->GetParameter(0) << " and error: " << f1->GetParError(0) << std::endl;



		std::cout << "RawU ADC: " << RawU << " and MeanU ADC: " << meanU << std::endl;
		std::cout << "CorrU ADC: " << CorrU << " and MeanCorrU ADC: " << meanCorrU << std::endl;
		std::cout << "RawV ADC: " << RawV << " and MeanV ADC: " << meanV << std::endl;
		std::cout << "CorrV ADC: " << CorrV << " and MeanCorrV ADC: " << meanCorrV << std::endl;
		std::cout << "RawY ADC: " << RawY << " and MeanY ADC: " << meanY << std::endl;
		std::cout << "CorrY ADC: " << CorrY << " and MeanCorrY ADC: " << meanCorrY << std::endl;
		std::cout << "RawAll ADC: " << RawAll << " and MeanAll ADC: " << meanAll << std::endl;
		std::cout << "CorrAll ADC: " << CorrAll << " and MeanCorrAll ADC: " << meanCorrAll << std::endl;

		TCanvas *ADC = new TCanvas("ADC","ADC",200,10,500,1000);
		ADC->Divide(1,4);
		for (size_t iPlane{1}; iPlane < 4; iPlane++) {
			ADC->cd(iPlane);  hADC_Corr[iPlane]->Draw(); hADC[iPlane]->SetLineColor(2); hADC[iPlane]->Draw("same");
		}
		ADC->cd(4); hADC_Corr_All->Draw(); hADC_All->Draw("same");

		TCanvas *cADC = new TCanvas("cADC","",500,400);
		cADC->cd();
		cADC->SetLeftMargin(0.13);
		cADC->SetBottomMargin(0.13);
		hADC[1]->SetLineColor(4);
		hADC[1]->SetLineWidth(4);
		hADC[2]->SetLineWidth(4);
		hADC[3]->SetLineWidth(4);
		hADC[2]->SetLineColor(3);
		hADC[3]->SetLineColor(2);
		hADC_All->SetLineColor(6);
		hADC_All->SetLineWidth(4);
		hADC[1]->GetXaxis()->SetTitle("dQ/dx (ADC/cm)");
		hADC[1]->GetYaxis()->SetTitle("Entries");
		hADC[1]->GetXaxis()->SetTitleOffset(1);
		hADC[1]->GetXaxis()->SetTitleSize(0.06);
		hADC[1]->GetXaxis()->SetLabelSize(0.05);
		hADC[1]->GetYaxis()->SetTitleOffset(1);
		hADC[1]->GetYaxis()->SetTitleSize(0.06);
		hADC[1]->GetYaxis()->SetLabelSize(0.05);
		hADC_All->SetTitle("");
		hADC_All->Draw();
		hADC[1]->Draw("same");
		hADC[3]->Draw("same");
		hADC[2]->Draw("same");
		TLegend *legADC = new TLegend(0.65,0.65,0.87,0.87);
    	legADC->AddEntry(hADC[1],"Plane 0","L");
		legADC->AddEntry(hADC[2],"Plane 1","L");
		legADC->AddEntry(hADC[3],"Plane 2","L");
		legADC->AddEntry(hADC_All,"All Planes","L");
		legADC->Draw("same");

		TCanvas *cADC_Corr = new TCanvas("cADC_Corr","",500,400);
		cADC_Corr->SetLeftMargin(0.13);
		cADC_Corr->SetBottomMargin(0.13);
		cADC_Corr->cd();
		cADC_Corr->SetLeftMargin(0.13);
		cADC_Corr->SetBottomMargin(0.13);
		hADC_Corr[1]->SetLineColor(4);
		hADC_Corr[1]->SetLineWidth(4);
		hADC_Corr[2]->SetLineWidth(4);
		hADC_Corr[3]->SetLineWidth(4);
		hADC_Corr[2]->SetLineColor(3);
		hADC_Corr[3]->SetLineColor(2);
		hADC_Corr_All->SetLineColor(1);
		hADC_Corr_All->SetLineWidth(4);
		hADC_Corr_All->GetXaxis()->SetTitle("dQ/dx (ADC/cm)");
		hADC_Corr_All->GetYaxis()->SetTitle("Entries");
		hADC_Corr_All->GetXaxis()->SetTitleOffset(1);
		hADC_Corr_All->GetXaxis()->SetTitleSize(0.06);
		hADC_Corr_All->GetXaxis()->SetLabelSize(0.05);
		hADC_Corr_All->GetYaxis()->SetTitleOffset(1);
		hADC_Corr_All->GetYaxis()->SetTitleSize(0.06);
		hADC_Corr_All->GetYaxis()->SetLabelSize(0.05);
		hADC_Corr_All->Draw();
		//hADC_Corr[1]->Draw("same");
		hADC_Corr[3]->Draw("same");
		//hADC_Corr[2]->Draw("same");
		cADC_Corr->Update();
		TLine *All_Corr = new TLine(9.04167,0,9.04167,cADC_Corr->GetUymax());
		TLine *Y_Corr = new TLine(8.875,0,8.875,cADC_Corr->GetUymax());
		All_Corr->SetLineColor(1);
		All_Corr->SetLineStyle(6);
		All_Corr->SetLineWidth(3);
		Y_Corr->SetLineColor(2);
		Y_Corr->SetLineStyle(6);
		Y_Corr->SetLineWidth(3);
		Y_Corr->Draw("same");
		All_Corr->Draw("same");

		TLegend *legADC_Corr = new TLegend(0.55,0.55,0.87,0.87);
    	//legADC_Corr->AddEntry(hADC_Corr[1],"Plane 0","L");
		//legADC_Corr->AddEntry(hADC_Corr[2],"Plane 1","L");
		legADC_Corr->AddEntry(hADC_Corr[3],"Plane 2","L");
		legADC_Corr->AddEntry(hADC_Corr_All,"All Planes","L");
		legADC_Corr->AddEntry(Y_Corr,"Most Probable Value (Plane 0)","L");
		legADC_Corr->AddEntry(All_Corr,"Most Probable Value (All Planes)","L");
		legADC_Corr->Draw("same");

		TCanvas *cADC_After_Y = new TCanvas("cADC_After_Y","",500,400);
		cADC_After_Y->cd();
		cADC_After_Y->SetLeftMargin(0.13);
		cADC_After_Y->SetBottomMargin(0.13);
		hADC_Y_Corr->SetLineColor(2);
		hADC_Y_Corr->SetLineWidth(4);
		hADC_Y_Un->SetLineColor(4);
		hADC_Y_Un->SetLineWidth(4);
		hADC_Y_Corr->Draw();
		hADC_Y_Un->Draw("same");

		TCanvas *cADC_After_All = new TCanvas("cADC_After_All","",500,400);
		cADC_After_All->cd();
		cADC_After_All->SetLeftMargin(0.13);
		cADC_After_All->SetBottomMargin(0.13);
		hADC_All_Cor->SetLineColor(2);
		hADC_All_Cor->SetLineWidth(4);
		hADC_All_Un->SetLineColor(4);
		hADC_All_Un->SetLineWidth(4);
		hADC_All_Cor->Draw();
		hADC_All_Un->Draw("same");



		TCanvas *Trackz = new TCanvas("Trackz","Trackz",200,10,750,1000);
		Trackz->Divide(2,4);
		int canvas_count {1};
		for (size_t iPlane {1}; iPlane < 4; iPlane++) {
			Trackz->cd(canvas_count);
			hUnDQdx_p[iPlane]->GetZaxis()->SetRangeUser(0,40);
			hUnDQdx_p[iPlane]->SetTitle("");
			hUnDQdx_p[iPlane]->GetXaxis()->SetTitle("Residual Range (cm)");
			hUnDQdx_p[iPlane]->GetYaxis()->SetTitle("Uncorrected dQ/dx (e/cm)");
			hUnDQdx_p[iPlane]->Draw("colz");
			canvas_count++;
			Trackz->cd(canvas_count);
			hUnDQdx_m[iPlane]->GetZaxis()->SetRangeUser(0, 80);
			hUnDQdx_m[iPlane]->SetTitle("");
			hUnDQdx_m[iPlane]->GetXaxis()->SetTitle("Residual Range (cm)");
			hUnDQdx_m[iPlane]->GetYaxis()->SetTitle("Uncorrected dQ/dx (e/cm)");
			hUnDQdx_m[iPlane]->Draw("colz");
			canvas_count++;
		}
		Trackz->cd(7);
		hUnDQdx_All_p->GetZaxis()->SetRangeUser(0,40);
		hUnDQdx_All_p->GetXaxis()->SetTitle("Residual Range (cm)");
		hUnDQdx_All_p->GetYaxis()->SetTitle("Uncorrected dQ/dx (e/cm)");
		hUnDQdx_All_p->Draw("colz");
		Trackz->cd(8);
		hUnDQdx_All_m->GetZaxis()->SetRangeUser(0,80);
		hUnDQdx_All_m->GetXaxis()->SetTitle("Residual Range (cm)");
		hUnDQdx_All_m->GetYaxis()->SetTitle("Uncorrected dQ/dx (e/cm)");
		hUnDQdx_All_m->Draw("colz");


		TCanvas *Trackz2 = new TCanvas("Trackz2","Trackz2",200,10,750,1000);
		Trackz2->Divide(2,4);
		canvas_count = 1;
		for (size_t iPlane {1}; iPlane < 4; iPlane++) {
			Trackz2->cd(canvas_count);
			hCorrDQdx_p[iPlane]->GetZaxis()->SetRangeUser(0,40);
			hCorrDQdx_p[iPlane]->GetXaxis()->SetTitle("Residual Range (cm)");
			hCorrDQdx_p[iPlane]->GetXaxis()->SetTitle("Residual Range (cm)");
			hCorrDQdx_p[iPlane]->GetYaxis()->SetTitle("Corrected dQ/dx (e/cm)");
			hCorrDQdx_p[iPlane]->Draw("colz");
			canvas_count++;
			Trackz2->cd(canvas_count);
			hCorrDQdx_m[iPlane]->GetZaxis()->SetRangeUser(0,80);
			hCorrDQdx_m[iPlane]->GetXaxis()->SetTitle("Residual Range (cm)");
			hCorrDQdx_m[iPlane]->GetYaxis()->SetTitle("Corrected dQ/dx (e/cm)");
			hCorrDQdx_m[iPlane]->Draw("colz");
			canvas_count++;
		}
		Trackz2->cd(7);
		hCorrDQdx_All_p->GetZaxis()->SetRangeUser(0,40);
		hCorrDQdx_All_p->GetXaxis()->SetTitle("Residual Range (cm)");
		hCorrDQdx_All_p->GetYaxis()->SetTitle("Corrected dQ/dx (e/cm)");
		hCorrDQdx_All_p->Draw("colz");
		Trackz2->cd(8);
		hCorrDQdx_All_m->GetZaxis()->SetRangeUser(0,80);
		hCorrDQdx_All_m->GetXaxis()->SetTitle("Residual Range (cm)");
		hCorrDQdx_All_m->GetYaxis()->SetTitle("Corrected dQ/dx (e/cm)");
		hCorrDQdx_All_m->Draw("colz");
		}

	void DEdx_Test::RecoFit_new() {
		fRecombination = new TF1("Recombination",RecombinationFunction,0,1000,6);
		fRecombination->SetParameters(alpha_ex, beta_ex, 1, Wion, LAr_Density, E_Field);
		fRecombination->SetParNames("alpha","beta","C","W_{ion}","#rho_{LAr}","E_{field}");
		fRecombination->FixParameter(2,1);
		fRecombination->FixParameter(3,Wion);
		fRecombination->FixParameter(4,LAr_Density);
		fRecombination->FixParameter(5,E_Field);

		TGraphErrors *g2 = new TGraphErrors(ebin.size(),&ebin[0],&mu_Corr_Y[0],&ey[0],&mu_error_Corr_Y[0]);

		TCanvas *cRecombination = new TCanvas("CorrectedRecombinationPlot_Y","CorrectedRecombinationPlot_Y",800,600);
		TF1 *fgaus = new TF1("fgaus","gaus(0)",0,1000);
		fgaus->SetParameters(100,200,20,10,200,100);

		TProfile *P_cor = hCorrectedRecombinationPlot_Y->ProfileX("P_cor");
		g2->Fit(fRecombination,"re","e",2.75,15);
		hCorrectedRecombinationPlot_Y->GetYaxis()->SetTitle("Correctected dQ/dx (e/cm)");
		hCorrectedRecombinationPlot_Y->GetXaxis()->SetTitle("Expected dE/dx (MeV/cm)");
		hCorrectedRecombinationPlot_Y->GetXaxis()->SetTitleOffset(1);
		hCorrectedRecombinationPlot_Y->GetXaxis()->SetTitleSize(0.06);
		hCorrectedRecombinationPlot_Y->GetXaxis()->SetLabelSize(0.05);
		hCorrectedRecombinationPlot_Y->GetYaxis()->SetTitleOffset(1);
		hCorrectedRecombinationPlot_Y->GetYaxis()->SetTitleSize(0.06);
		hCorrectedRecombinationPlot_Y->GetYaxis()->SetLabelSize(0.05);
		//hCorrectedRecombinationPlot_Y->GetZaxis()->SetRangeUser(0, 220);
		hCorrectedRecombinationPlot_Y->Draw("colz");
		fRecombination->SetLineWidth(3);
		fRecombination->Draw("same");
		fRecombinationExpected->SetLineColor(3);
		fRecombinationExpected->SetLineWidth(3);
		fRecombinationExpected->Draw("same");
		//P_cor->Draw("sames");
		g2->SetLineColor(6);
		g2->Draw("same");

		fRecombination2 = new TF1("Recombination2",RecombinationFunction,0,1000,6);
		fRecombination2->SetParameters(alpha_ex, beta_ex, 1, Wion, LAr_Density, E_Field);
		fRecombination2->SetParNames("alpha","beta","C","W_{ion}","#rho_{LAr}","E_{field}");
		fRecombination2->FixParameter(2,1);
		fRecombination2->FixParameter(3,Wion);
		fRecombination2->FixParameter(4,LAr_Density);
		fRecombination2->FixParameter(5,E_Field);

		TGraphErrors *g3 = new TGraphErrors(ebin.size(),&ebin[0],&mu_Un_Y[0],&ey[0],&mu_error_Un_Y[0]);

		TCanvas *cRecombination2 = new TCanvas("UnCorrectedRecombinationPlot_Y","UnCorrectedRecombinationPlot_Y",800,600);
		TF1 *fgaus2 = new TF1("fgaus2","gaus(0)",0,1000);
		fgaus2->SetParameters(100,200,20,10,200,100);

		TProfile *P_cor2 = hUncorrectedRecombinationPlot_Y->ProfileX("P_cor2");
		g3->Fit(fRecombination2,"re","e",2.75,15);
		hUncorrectedRecombinationPlot_Y->GetXaxis()->SetTitle("Uncorrectected dQ/dx (e/cm)");
		hUncorrectedRecombinationPlot_Y->GetYaxis()->SetTitle("Expected dE/dx (MeV/cm)");

		hUncorrectedRecombinationPlot_Y->Draw("colz");
		fRecombination2->SetLineWidth(3);
		fRecombination2->Draw("same");
		fRecombinationExpected->SetLineColor(3);
		fRecombinationExpected->SetLineWidth(3);
		fRecombinationExpected->Draw("same");
		//P_cor->Draw("sames");
		g3->SetLineColor(6);
		g3->Draw("same");

		fRecombination3 = new TF1("Recombination3",RecombinationFunction,0,1000,6);
		fRecombination3->SetParameters(alpha_ex, beta_ex, 1, Wion, LAr_Density, E_Field);
		fRecombination3->SetParNames("alpha","beta","C","W_{ion}","#rho_{LAr}","E_{field}");
		fRecombination3->FixParameter(2,1);
		fRecombination3->FixParameter(3,Wion);
		fRecombination3->FixParameter(4,LAr_Density);
		fRecombination3->FixParameter(5,E_Field);

		TGraphErrors *g4 = new TGraphErrors(ebin.size(),&ebin[0],&mu_Un_All[0],&ey[0],&mu_error_Un_All[0]);

		TCanvas *cRecombination3 = new TCanvas("UnCorrectedRecombinationPlot_All","UnCorrectedRecombinationPlot_All",800,600);
		TF1 *fgaus3 = new TF1("fgaus3","gaus(0)",0,1000);
		fgaus3->SetParameters(100,200,20,10,200,100);

		TProfile *P_cor3 = hUncorrectedRecombinationPlot_Y->ProfileX("P_cor3");
		g4->Fit(fRecombination3,"re","e",2.75,15);
		hUncorrectedRecombinationPlot_All->GetYaxis()->SetTitle("Uncorrectected dQ/dx (e/cm)");
		hUncorrectedRecombinationPlot_All->GetXaxis()->SetTitle("Expected dE/dx (MeV/cm)");
		//hUncorrectedRecombinationPlot_All->GetZaxis()->SetRangeUser(0, 220);
		hUncorrectedRecombinationPlot_All->GetXaxis()->SetTitleOffset(1);
		hUncorrectedRecombinationPlot_All->GetXaxis()->SetTitleSize(0.06);
		hUncorrectedRecombinationPlot_All->GetXaxis()->SetLabelSize(0.05);
		hUncorrectedRecombinationPlot_All->GetYaxis()->SetTitleOffset(1);
		hUncorrectedRecombinationPlot_All->GetYaxis()->SetTitleSize(0.06);
		hUncorrectedRecombinationPlot_All->GetYaxis()->SetLabelSize(0.05);
		hUncorrectedRecombinationPlot_All->Draw("colz");
		fRecombination3->SetLineWidth(3);
		fRecombination3->Draw("same");
		fRecombinationExpected->SetLineColor(3);
		fRecombinationExpected->SetLineWidth(3);
		fRecombinationExpected->Draw("same");
		//P_cor->Draw("sames");
		g4->SetLineColor(6);
		g4->Draw("same");
		TLegend *legdEdx = new TLegend(0.55,0.55,0.87,0.87);
    	//legADC_Corr->AddEntry(hADC_Corr[1],"Plane 0","L");
		//legADC_Corr->AddEntry(hADC_Corr[2],"Plane 1","L");
		legdEdx->AddEntry(fRecombination3,"Best Fit","L");
		legdEdx->AddEntry(fRecombinationExpected,"Expected Fit","L");
		legdEdx->AddEntry(g4,"Most Probable Values","L");
		legdEdx->Draw("same");

		fRecombination4 = new TF1("Recombination4",RecombinationFunction,0,1000,6);
		fRecombination4->SetParameters(alpha_ex, beta_ex, 1, Wion, LAr_Density, E_Field);
		fRecombination4->SetParNames("alpha","beta","C","W_{ion}","#rho_{LAr}","E_{field}");
		fRecombination4->FixParameter(2,1);
		fRecombination4->FixParameter(3,Wion);
		fRecombination4->FixParameter(4,LAr_Density);
		fRecombination4->FixParameter(5,E_Field);

		TGraphErrors *g5 = new TGraphErrors(ebin.size(),&ebin[0],&mu_Corr_All[0],&ey[0],&mu_error_Corr_All[0]);

		TCanvas *cRecombination4 = new TCanvas("CorrectedRecombinationPlot_All","CorrectedRecombinationPlot_Y",800,600);
		TF1 *fgaus4 = new TF1("fgaus4","gaus(0)",0,1000);
		fgaus2->SetParameters(100,200,20,10,200,100);

		TProfile *P_cor4 = hCorrectedRecombinationPlot_Y->ProfileX("P_cor4");
		g5->Fit(fRecombination4,"re","e",2.75,15);
		hCorrectedRecombinationPlot_All->GetYaxis()->SetTitle("Correctected dQ/dx (e/cm)");
		hCorrectedRecombinationPlot_All->GetXaxis()->SetTitle("Expected dE/dx (MeV/cm)");
		//hCorrectedRecombinationPlot_All->GetZaxis()->SetRangeUser(0, 220);
		hCorrectedRecombinationPlot_All->GetXaxis()->SetTitleOffset(1);
		hCorrectedRecombinationPlot_All->GetXaxis()->SetTitleSize(0.06);
		hCorrectedRecombinationPlot_All->GetXaxis()->SetLabelSize(0.05);
		hCorrectedRecombinationPlot_All->GetYaxis()->SetTitleOffset(1);
		hCorrectedRecombinationPlot_All->GetYaxis()->SetTitleSize(0.06);
		hCorrectedRecombinationPlot_All->GetYaxis()->SetLabelSize(0.05);
		hCorrectedRecombinationPlot_All->Draw("colz");
		fRecombination4->SetLineWidth(3);
		fRecombination4->Draw("same");
		fRecombinationExpected->SetLineColor(3);
		fRecombinationExpected->SetLineWidth(3);
		fRecombinationExpected->Draw("same");
		//P_cor->Draw("sames");
		g5->SetLineColor(6);
		g5->Draw("same");
		TLegend *legdEdx_Corr = new TLegend(0.55,0.55,0.87,0.87);
    	//legADC_Corr->AddEntry(hADC_Corr[1],"Plane 0","L");
		//legADC_Corr->AddEntry(hADC_Corr[2],"Plane 1","L");
		legdEdx_Corr->AddEntry(fRecombination4,"Best Fit","L");
		legdEdx_Corr->AddEntry(fRecombinationExpected,"Expected Fit","L");
		legdEdx_Corr->AddEntry(g5,"Most Probable Values","L");
		legdEdx_Corr->Draw("same");

		alpha_fit_Y = fRecombination->GetParameter(0);
		beta_fit_Y = fRecombination->GetParameter(1);
		alpha_fit_All = fRecombination4->GetParameter(0);
		beta_fit_All = fRecombination4->GetParameter(1);

		fdEdxConvert_Y = new TF1("fDEdxConvert_Y",dEdxFunction,0,1000000,6);
		fdEdxConvert_Y->SetParameters(alpha_fit_Y,beta_fit_Y,1,Wion,LAr_Density,E_Field);
		fdEdxConvert_All = new TF1("fDEdxConvert_All",dEdxFunction,0,1000000,6);
		fdEdxConvert_All->SetParameters(alpha_fit_All,beta_fit_All,1,Wion,LAr_Density,E_Field);


	}

	void DEdx_Test::All_DEdxPlots() {
		std::cout << "DEDX PLOT TIME" << std::endl;
		const int ix {30};
		double rr_start {0}, rr_end {100};
		double add_len = (rr_end - rr_start) / ix;
		double bin_start = rr_start;
		double bin_end = rr_start + add_len;
		std::vector<double> ResBins {};
		TH1D* hResPlotBin_All[ix];
		TH1D* hResPlotBin_Y[ix];
		TH1D* hResPlotBin_All_m[ix];
		TH1D* hResPlotBin_Y_m[ix];
		for (int idx {0}; idx < ix; idx++){
			ResBins.push_back(bin_end);
			hResPlotBin_All[idx] = new TH1D(Form("hResPlotBin_All_%02d",idx),Form("hResPlotBin_All_%02d",idx),70,0,35);
			hResPlotBin_Y[idx] = new TH1D(Form("hResPlotBin_Y_%02d",idx),Form("hResPlotBin_Y_%02d",idx),70,0,35);
			hResPlotBin_All_m[idx] = new TH1D(Form("hResPlotBin_All_m_%02d",idx),Form("hResPlotBin_All_m_%02d",idx),70,0,35);
			hResPlotBin_Y_m[idx] = new TH1D(Form("hResPlotBin_Y_m_%02d",idx),Form("hResPlotBin_Y_m_%02d",idx),70,0,35);
			bin_start+=add_len;
			bin_end+=add_len;
		}

		std::cout << "resplot stuff made" << std::endl;

		hCorrectedDEdx_Y = new TH2D("hCorrectedDEdx_Y","hCorrectedDEdx_Y",400,0,400,200,0,20);
		hCorrectedDEdx_Y_p = new TH2D("hCorrectedDEdx_Y_p","hCorrectedDEdx_Y_p",25,0,50,100,0,20);
		hCorrectedDEdx_Y_m = new TH2D("hCorrectedDEdx_Y_m","hCorrectedDEdx_Y_m",50,0,100,50,0,20);
		hCorrectedDEdx_All = new TH2D("hCorrectedDEdx_All","hCorrectedDEdx_All",400,0,400,300,0,30);
		hCorrectedDEdx_All_p = new TH2D("hCorrectedDEdx_All_p","hCorrectedDEdx_All_p",25,0,50,100,0,20);
		hCorrectedDEdx_All_m = new TH2D("hCorrectedDEdx_All_m","hCorrectedDEdx_All_m",50,0,100,50,0,20);

		for ( size_t ivert = 0; ivert < AllVertex_p.size(); ivert++) {
			double avg_dQdx[2] {0,0};
			double avg_dQdx_10cm[2] {0,0};
			double v_x = AllVertex_p[ivert][0].Vertex().X();
			double v_y = AllVertex_p[ivert][0].Vertex().Y();
			double v_z = AllVertex_p[ivert][0].Vertex().Z();
			//double dr = deltaR(v_x, v_y, v_z, AllTruth_p[ivert]);
			bool InFiducial = InFid(v_x, v_y, v_z);
			double dr = 1;

			if (dr > 5 || InFiducial == 0) continue;
				for ( size_t itrack = 0; itrack < AllVertex_p[ivert].size(); itrack++) {
					int points {0};
					int points_10cm {0};
					auto curr_vert = AllVertex_p[ivert][itrack];
					for (size_t iloc = 0; iloc < curr_vert.NumberTrajectoryPoints(); iloc++) {
						double raw_dQdx_U = curr_vert.DQdxAtPoint(iloc,larlite::geo::kU);
						double raw_dQdx_V = curr_vert.DQdxAtPoint(iloc,larlite::geo::kV);
						double raw_dQdx_Y = curr_vert.DQdxAtPoint(iloc,larlite::geo::kZ);
						TVector3 point = curr_vert.LocationAtPoint(iloc);
						double reslength = curr_vert.Length(iloc);
						if (reslength < 1) continue;
						if (raw_dQdx_U <= 0 || raw_dQdx_V <= 0 || raw_dQdx_Y <= 0 || raw_dQdx_U == 37. || raw_dQdx_V == 37. || raw_dQdx_Y == 37.) continue;
						points += 1;
						avg_dQdx[itrack] += raw_dQdx_U + raw_dQdx_V + raw_dQdx_Y;
						if (InFid(point.X(), point.Y(), point.Z()) ==  false) continue;
						//hUncorrectedDQdx->Fill(reslength, raw_dQdx*avgVal);
						}
					avg_dQdx[itrack] /= points;
					}
				int kProton = 0;
				int kMuon = 1;
				if (AllVertex_p[ivert][kProton].Length() < 5 || AllVertex_p[ivert][kMuon].Length() < 5) continue;
				if (AllVertex_p[ivert][kProton].Length() < AllVertex_p[ivert][kMuon].Length()) {
					if (avg_dQdx[kProton] < avg_dQdx[kMuon]) continue;
					for (size_t iloc = 0; iloc < AllVertex_p[ivert][kProton].NumberTrajectoryPoints(); iloc++) {
						auto curr_vert = AllVertex_p[ivert][kProton];
						double raw_dQdx_U = curr_vert.DQdxAtPoint(iloc,larlite::geo::kU);
						double raw_dQdx_V = curr_vert.DQdxAtPoint(iloc,larlite::geo::kV);
						double raw_dQdx_Y = curr_vert.DQdxAtPoint(iloc,larlite::geo::kZ);
						double raw_dQdx = raw_dQdx_U + raw_dQdx_V + raw_dQdx_Y;
						TVector3 point = curr_vert.LocationAtPoint(iloc);
						int bin = CorrectionMap3D->FindBin(point.X(),point.Y(),point.Z());
						double corr_dQdx_U = raw_dQdx_U*CorrectionMap3D_U->GetBinContent(bin);
						double corr_dQdx_V = raw_dQdx_V*CorrectionMap3D_V->GetBinContent(bin);
						double corr_dQdx_Y = raw_dQdx_Y*CorrectionMap3D->GetBinContent(bin);
						double corr_dQdx = corr_dQdx_U + corr_dQdx_V + corr_dQdx_Y;
						double reslength = curr_vert.Length(iloc);
						if (reslength < 0.25) continue;
						if (raw_dQdx_Y > 0 && (raw_dQdx_Y > 37.01 || raw_dQdx_Y < 36.99)) {  
							double corr_dEdx_Y = fdEdxConvert_Y->Eval(corr_dQdx_Y * avgVal_Y_Corr);
							hCorrectedDEdx_Y_p->Fill(reslength,corr_dEdx_Y);
							hCorrectedDEdx_Y->Fill(reslength,corr_dEdx_Y);
							int index {1};
							while (reslength > ResBins[index] && index < 25 ) {
								index++;
							}
							hResPlotBin_Y[index]->Fill(corr_dEdx_Y);
						}
						if (raw_dQdx_U <= 0 || raw_dQdx_V <= 0 || raw_dQdx_Y <= 0 || (raw_dQdx_U < 37.01 && raw_dQdx_U > 36.99) || (raw_dQdx_V < 37.01 && raw_dQdx_V > 36.99) || (raw_dQdx_Y < 37.01 && raw_dQdx_Y > 36.99)) continue;
						if (InFid(point.X(), point.Y(), point.Z()) ==  false) continue;
							double corr_dEdx = fdEdxConvert_All->Eval(corr_dQdx * avgVal_all_Corr/3);
							hCorrectedDEdx_All_p->Fill(reslength,corr_dEdx);
							hCorrectedDEdx_All->Fill(reslength,corr_dEdx);
							int index {1};
							while (reslength > ResBins[index] && index < 25) {
								index++;
							}
							hResPlotBin_All[index]->Fill(corr_dEdx);
						}

					}
					for (size_t iloc = 0; iloc < AllVertex_p[ivert][kMuon].NumberTrajectoryPoints(); iloc++) {
						auto curr_vert = AllVertex_p[ivert][kMuon];
						double raw_dQdx_U = curr_vert.DQdxAtPoint(iloc,larlite::geo::kU);
						double raw_dQdx_V = curr_vert.DQdxAtPoint(iloc,larlite::geo::kV);
						double raw_dQdx_Y = curr_vert.DQdxAtPoint(iloc,larlite::geo::kZ);
						double raw_dQdx = raw_dQdx_U + raw_dQdx_V + raw_dQdx_Y;
						TVector3 point = curr_vert.LocationAtPoint(iloc);
						int bin = CorrectionMap3D->FindBin(point.X(),point.Y(),point.Z());
						double corr_dQdx_U = raw_dQdx_U*CorrectionMap3D_U->GetBinContent(bin);
						double corr_dQdx_V = raw_dQdx_V*CorrectionMap3D_V->GetBinContent(bin);
						double corr_dQdx_Y = raw_dQdx_Y*CorrectionMap3D->GetBinContent(bin);
						double corr_dQdx = corr_dQdx_U + corr_dQdx_V + corr_dQdx_Y;
						double reslength = curr_vert.Length(iloc);
						if (reslength < 0.25) continue;
						if (raw_dQdx_Y > 0 && (raw_dQdx_Y > 37.01 || raw_dQdx_Y < 36.99)) {  
							double corr_dEdx_Y = fdEdxConvert_Y->Eval(corr_dQdx_Y* avgVal_Y_Corr);
							hCorrectedDEdx_Y_m->Fill(reslength,corr_dEdx_Y);
							hCorrectedDEdx_Y->Fill(reslength,corr_dEdx_Y);
							int index {1};
							while (reslength > ResBins[index] && index < 25 ) {
								index++;
							}
							hResPlotBin_Y_m[index]->Fill(corr_dEdx_Y);
						}
						if (raw_dQdx_U <= 0 || raw_dQdx_V <= 0 || raw_dQdx_Y <= 0 || (raw_dQdx_U < 37.01 && raw_dQdx_U > 36.99) || (raw_dQdx_V < 37.01 && raw_dQdx_V > 36.99) || (raw_dQdx_Y < 37.01 && raw_dQdx_Y > 36.99)) continue;
						if (InFid(point.X(), point.Y(), point.Z()) ==  false) continue;
							double corr_dEdx = fdEdxConvert_All->Eval(corr_dQdx* avgVal_all_Corr/3);
							hCorrectedDEdx_All_m->Fill(reslength,corr_dEdx);
							hCorrectedDEdx_All->Fill(reslength,corr_dEdx);
							int index {1};
							while (reslength > ResBins[index] && index < 25 ) {
								index++;
							}
							hResPlotBin_All_m[index]->Fill(corr_dEdx);
						}	
				
				kProton = 0;
				kMuon = 1;
				if (AllVertex_p[ivert][kProton].Length() < 5 || AllVertex_p[ivert][kMuon].Length() < 5) continue;
				if (AllVertex_p[ivert][kProton].Length() < AllVertex_p[ivert][kMuon].Length()) {
					if (avg_dQdx[kProton] < avg_dQdx[kMuon]) continue;
					for (size_t iloc = 0; iloc < AllVertex_p[ivert][kProton].NumberTrajectoryPoints(); iloc++) {
						auto curr_vert = AllVertex_p[ivert][kProton];
						double raw_dQdx_U = curr_vert.DQdxAtPoint(iloc,larlite::geo::kU);
						double raw_dQdx_V = curr_vert.DQdxAtPoint(iloc,larlite::geo::kV);
						double raw_dQdx_Y = curr_vert.DQdxAtPoint(iloc,larlite::geo::kZ);
						double raw_dQdx = raw_dQdx_U + raw_dQdx_V + raw_dQdx_Y;
						TVector3 point = curr_vert.LocationAtPoint(iloc);
						int bin = CorrectionMap3D->FindBin(point.X(),point.Y(),point.Z());
						double corr_dQdx_U = raw_dQdx_U*CorrectionMap3D_U->GetBinContent(bin);
						double corr_dQdx_V = raw_dQdx_V*CorrectionMap3D_V->GetBinContent(bin);
						double corr_dQdx_Y = raw_dQdx_Y*CorrectionMap3D->GetBinContent(bin);
						double corr_dQdx = corr_dQdx_U + corr_dQdx_V + corr_dQdx_Y;
						double reslength = curr_vert.Length(iloc);
						if (reslength < 0.25) continue;
						if (raw_dQdx_Y > 0 && (raw_dQdx_Y > 37.01 || raw_dQdx_Y < 36.99)) {  
							double corr_dEdx_Y = fdEdxConvert_Y->Eval(corr_dQdx_Y* avgVal_Y_Corr);
							//std::cout << corr_dEdx_Y << std::endl;
							hCorrectedDEdx_Y_p->Fill(reslength,corr_dEdx_Y);
							hCorrectedDEdx_Y->Fill(reslength,corr_dEdx_Y);
							int index {1};
							while (reslength > ResBins[index] && index < 25) {
								index++;
							}
							hResPlotBin_Y[index]->Fill(corr_dEdx_Y);
						}
						if (raw_dQdx_U <= 0 || raw_dQdx_V <= 0 || raw_dQdx_Y <= 0 || (raw_dQdx_U < 37.01 && raw_dQdx_U > 36.99) || (raw_dQdx_V < 37.01 && raw_dQdx_V > 36.99) || (raw_dQdx_Y < 37.01 && raw_dQdx_Y > 36.99)) continue;
						if (InFid(point.X(), point.Y(), point.Z()) ==  false) continue;
							double corr_dEdx = fdEdxConvert_All->Eval(corr_dQdx* avgVal_all_Corr/3);
							hCorrectedDEdx_All_p->Fill(reslength,corr_dEdx);
							hCorrectedDEdx_All->Fill(reslength,corr_dEdx);
							int index {1};
							while (reslength > ResBins[index] && index < 25) {
								index++;
							}
							hResPlotBin_All[index]->Fill(corr_dEdx);
						}			
					}
					for (size_t iloc = 0; iloc < AllVertex_p[ivert][kMuon].NumberTrajectoryPoints(); iloc++) {
						auto curr_vert = AllVertex_p[ivert][kMuon];
						double raw_dQdx_U = curr_vert.DQdxAtPoint(iloc,larlite::geo::kU);
						double raw_dQdx_V = curr_vert.DQdxAtPoint(iloc,larlite::geo::kV);
						double raw_dQdx_Y = curr_vert.DQdxAtPoint(iloc,larlite::geo::kZ);
						double raw_dQdx = raw_dQdx_U + raw_dQdx_V + raw_dQdx_Y;
						TVector3 point = curr_vert.LocationAtPoint(iloc);
						int bin = CorrectionMap3D->FindBin(point.X(),point.Y(),point.Z());
						double corr_dQdx_U = raw_dQdx_U*CorrectionMap3D_U->GetBinContent(bin);
						double corr_dQdx_V = raw_dQdx_V*CorrectionMap3D_V->GetBinContent(bin);
						double corr_dQdx_Y = raw_dQdx_Y*CorrectionMap3D->GetBinContent(bin);
						double corr_dQdx = corr_dQdx_U + corr_dQdx_V + corr_dQdx_Y;
						double reslength = curr_vert.Length(iloc);
						if (reslength < 0.25) continue;
						if (raw_dQdx_Y > 0 && (raw_dQdx_Y > 37.01 || raw_dQdx_Y < 36.99)) {  
							double corr_dEdx_Y = fdEdxConvert_Y->Eval(corr_dQdx_Y* avgVal_Y_Corr);
							hCorrectedDEdx_Y_m->Fill(reslength,corr_dEdx_Y);
							hCorrectedDEdx_Y->Fill(reslength,corr_dEdx_Y);
							int index {1};
							while (reslength > ResBins[index] && index < 25 ) {
								index++;
							}
							hResPlotBin_Y_m[index]->Fill(corr_dEdx_Y);
						}
						if (raw_dQdx_U <= 0 || raw_dQdx_V <= 0 || raw_dQdx_Y <= 0 || (raw_dQdx_U < 37.01 && raw_dQdx_U > 36.99) || (raw_dQdx_V < 37.01 && raw_dQdx_V > 36.99) || (raw_dQdx_Y < 37.01 && raw_dQdx_Y > 36.99)) continue;
						if (InFid(point.X(), point.Y(), point.Z()) ==  false) continue;
							double corr_dEdx = fdEdxConvert_All->Eval(corr_dQdx* avgVal_all_Corr/3);
							hCorrectedDEdx_All_m->Fill(reslength,corr_dEdx);
							hCorrectedDEdx_All->Fill(reslength,corr_dEdx);
							int index {1};
							while (reslength > ResBins[index] && index < 25 ) {
								index++;
							}
							hResPlotBin_All_m[index]->Fill(corr_dEdx);
					}	
				}

			for (int idx {1}; idx < ix; idx++) {
				double mean_Y = hResPlotBin_Y[idx]->GetMean();
				double rms_Y = hResPlotBin_Y[idx]->GetRMS();
				double mean_All = hResPlotBin_All[idx]->GetMean();
				double rms_All = hResPlotBin_All[idx]->GetRMS();
				double mean_Y_m = hResPlotBin_Y_m[idx]->GetMean();
				double rms_Y_m = hResPlotBin_Y_m[idx]->GetRMS();
				double mean_All_m = hResPlotBin_All_m[idx]->GetMean();
				double rms_All_m = hResPlotBin_All_m[idx]->GetRMS();
				TF1 *f1 = new TF1("f1","gaus",0,50);
				TF1 *f2 = new TF1("f2","gaus",0,50);
				TF1 *f3 = new TF1("f3","gaus",0,50);
				TF1 *f4 = new TF1("f4","gaus",0,50);
				hResPlotBin_Y[idx]->Fit("f1","","",mean_Y-rms_Y, mean_Y+rms_Y);
				std::cout << "f1_p1: " << f1->GetParameter(1) << std::endl;
				hResPlotBin_All[idx]->Fit("f2","","",mean_All-rms_All,mean_All-rms_All);
				std::cout << "f2_p1: " << f2->GetParameter(1) << std::endl;
				hResPlotBin_Y_m[idx]->Fit("f3","","",mean_Y_m-rms_Y_m, mean_Y_m+rms_Y_m);
				std::cout << "f3_p1: " << f3->GetParameter(1) << std::endl;
				hResPlotBin_All_m[idx]->Fit("f4","","",mean_All_m-rms_All_m,mean_All_m-rms_All_m);


				mu_Y.push_back(f1->GetParameter(1));
				mu_error_Y.push_back(f1->GetParameter(2));
				mu_All.push_back(f2->GetParameter(1));
				mu_error_All.push_back(f2->GetParameter(2));
				mu_Y_m.push_back(f3->GetParameter(1));
				mu_error_Y_m.push_back(f3->GetParameter(2));
				mu_All_m.push_back(f4->GetParameter(1));
				mu_error_All_m.push_back(f4->GetParameter(2));

				std::cout << "muon: " <<  mu_Y_m[idx-1] << std::endl;

				double val = ResBins[idx] - 25/ix;
				rbin.push_back(val);
				ry.push_back(50/ix);
				std::cout << rbin[idx-1] << std::endl;

			}
			TProfile *prot_All = hCorrectedDEdx_All_p->ProfileX("prot_All");
			TProfile *prot_Y = hCorrectedDEdx_Y_p->ProfileX("prot_Y");
			TProfile *muo_All = hCorrectedDEdx_All_m->ProfileX("muo_All");
			TProfile *muo_Y = hCorrectedDEdx_Y_m->ProfileX("muo_Y");

			TGraphErrors *r1 = new TGraphErrors(rbin.size(),&rbin[0],&mu_Y[0],&ry[0],&mu_error_Y[0]);
			TGraphErrors *r2 = new TGraphErrors(rbin.size(),&rbin[0],&mu_All[0],&ry[0],&mu_error_All[0]);
			TGraphErrors *r3 = new TGraphErrors(rbin.size(),&rbin[0],&mu_Y_m[0],&ry[0],&mu_error_Y_m[0]);
			TGraphErrors *r4 = new TGraphErrors(rbin.size(),&rbin[0],&mu_All_m[0],&ry[0],&mu_error_All_m[0]);

			hCorrectedDEdx_All->GetXaxis()->SetTitle("Residual Range (cm)");
			hCorrectedDEdx_All->GetYaxis()->SetTitle("Corrected dE/dx (MeV/cm)");
			hCorrectedDEdx_All_p->GetXaxis()->SetTitle("Residual Range (cm)");
			hCorrectedDEdx_All_p->GetYaxis()->SetTitle("Corrected dE/dx (MeV/cm)");
			hCorrectedDEdx_All_m->GetXaxis()->SetTitle("Residual Range (cm)");
			hCorrectedDEdx_All_m->GetYaxis()->SetTitle("Corrected dE/dx (MeV/cm)");
			hCorrectedDEdx_Y->GetXaxis()->SetTitle("Residual Range (cm)");
			hCorrectedDEdx_Y->GetYaxis()->SetTitle("Corrected dE/dx (MeV/cm)");
			hCorrectedDEdx_Y_p->GetXaxis()->SetTitle("Residual Range (cm)");
			hCorrectedDEdx_Y_p->GetYaxis()->SetTitle("Corrected dE/dx (MeV/cm)");
			hCorrectedDEdx_Y_m->GetXaxis()->SetTitle("Residual Range (cm)");
			hCorrectedDEdx_Y_m->GetYaxis()->SetTitle("Corrected dE/dx (MeV/cm)");


			prot_All->SetMarkerColor(6);
			prot_Y->SetMarkerColor(6);
			muo_All->SetMarkerColor(6);
			muo_Y->SetMarkerColor(6);
			sMuonRange2dEdx->SetLineColor(3);
			sMuonRange2dEdx->SetLineWidth(3);
			sProtonRange2dEdx->SetLineColor(2);
			sProtonRange2dEdx->SetLineWidth(3);
			TCanvas *cDEdxPlots = new TCanvas("DEdxPlots","DEdxPlots",750,750);
			cDEdxPlots->Divide(2,3);
			cDEdxPlots->cd(1); hCorrectedDEdx_All->Draw("colz"); sMuonRange2dEdx->Draw("same"); sProtonRange2dEdx->Draw("same");
			cDEdxPlots->cd(2); hCorrectedDEdx_Y->Draw("colz"); sMuonRange2dEdx->Draw("same"); sProtonRange2dEdx->Draw("same");
			cDEdxPlots->cd(3); hCorrectedDEdx_All_p->Draw("colz"); sProtonRange2dEdx->Draw("same"); r2->Draw("same");//prot_All->Draw("same");
			cDEdxPlots->cd(4); hCorrectedDEdx_Y_p->Draw("colz"); sProtonRange2dEdx->Draw("same"); r1->Draw("same");//prot_Y->Draw("same");
			cDEdxPlots->cd(5); hCorrectedDEdx_All_m->Draw("colz"); sMuonRange2dEdx->Draw("same"); r4->Draw("same");//muo_All->Draw("same");
			cDEdxPlots->cd(6); hCorrectedDEdx_Y_m->Draw("colz"); sMuonRange2dEdx->Draw("same"); r3->Draw("same");//muo_Y->Draw("same");

			TCanvas *stfu = new TCanvas("stfu","stfu",750,400);
			stfu->Divide(2,1);
			stfu->cd(1); hCorrectedDEdx_All_p->Draw("colz"); sProtonRange2dEdx->Draw("same");
			stfu->cd(2); hCorrectedDEdx_All_m->Draw("colz"); sMuonRange2dEdx->Draw("same");

		}



};










#endif
