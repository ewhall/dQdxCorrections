#ifndef CAL_TEST_CXX
#define CAL_TEST_CXX

#include "Cal_Test.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TLine.h"
#include "TLatex.h"
#include "TH2D.h"
#include "TVirtualFitter.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
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

bool InFiducial(double x, double y, double z) {
	/* 
		To remove points along the edge of the detector.  Does not remove points 
		in the deadwire regions, like in the nu_mu and nu_e selections.
	*/
	if(x < 0+15 || x > 256.35-15)return 0;
	if(y < -116.5+15 || y > 116.5-25)return 0;
	if(z < 0+15 || z > 1036.8-15)return 0;
	else return 1;
}

double RecombinationFun(double *x, double *par){
	/*
		Recombination Function, note not really used in this code -i think-
	*/
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


namespace larlite {

	bool Cal_Test::Initialize() {
		/* 
			Set detector size, number of bins, recombination parameters and expected function.
		*/ 
		xmin = 0.0;	xmax = 256.25;
		ymin = -116.5; ymax = 116.5;
		zmin = 0.0; zmax = 1036.8;
		Nz = 200; Ny = 50, Nx = 50; // number of bins to be made along that axis
		plane = larlite::geo::kZ; // change this to expand to 3 planes

		alpha_ex = 0.93;
		beta_ex = 0.212;
		Wion = 23.6e-6;
		E_Field = 0.273;
		LAr_Density = 1.38;

		gStyle->SetOptStat(0);

		fRecombinationExpected = new TF1("fRecombinationExpected",RecombinationFun,0,1000,6);
        fRecombinationExpected->SetParameters(alpha_ex,beta_ex,1,Wion,LAr_Density,E_Field);
        fRecombinationExpected->SetParNames("alpha","beta","C","W_{ion}","#rho_{LAr}","E_{field}");
        fRecombinationExpected->FixParameter(0, alpha_ex);
        fRecombinationExpected->FixParameter(1, beta_ex);
        fRecombinationExpected->FixParameter(3, Wion);
        fRecombinationExpected->FixParameter(4,LAr_Density);
        fRecombinationExpected->FixParameter(5,E_Field);

		return true;
	}

	bool Cal_Test::Finalize() {
		MakedQdxPlots();
		std::cout << "Made dQdxPlots" << std::endl;
		Make2DCorrectionMap();
		std::cout << "Made 2D Correction Map" << std::endl;
		MakeXCorrectionMap();
		std::cout << "Make XCorrection Map" << std::endl;
		Make3DCorretionMap();


		std::cout << "Step 1" << std::endl;
		MakeVoxelCorrectionMap();
		std::cout << "Step 2" << std::endl;

		MakeCorrecteddQdxPlots();
		CosmicTracks();

		std::cout << "Step 3" << std::endl;
		TFile f("CalibrationFile.root","RECREATE");
		TCanvas *c1 = new TCanvas("c1", "c1", 200, 10, 750, 1000);
		c1->Divide(1,3);
		c1->cd(1); hUncorrecteddQdx->Draw("colz");	
		c1->cd(2); hUncorrecteddQdx_m->Draw("colz"); 
		c1->cd(3); hUncorrecteddQdx_p->Draw("colz"); 
		c1->Update();

		std::cout << "Step 4" << std::endl;
		// Make2DCorrectionMap();
		TCanvas *c2 = new TCanvas("c2", "c2", 200, 10, 750, 1000);
		c2->Divide(1,4);
		c2->cd(1); hHitMap2D->Draw("colz"); //hHitMap2D->Write();
		c2->cd(2); hRawDQdxMap2D->Draw("colz"); //hRawDQdxMap2D->Write();
		c2->cd(3); hAvgDQdxMap2D->GetZaxis()->SetRangeUser(0,25); //hAvgDQdxMap2D->Draw("colz"); hAvgDQdxMap2D->Write();
		c2->cd(4); hMedianDQdxMap2D->Draw("colz"); //hMedianDQdxMap2D->Write();
		c2->Update();
		c2->Print("hitmap.pdf");

		std::cout << "Step 5" << std::endl;
		TCanvas *c3 = new TCanvas("c3", "c3", 200, 10, 750, 1000);
		c3->cd();
		hCalibrationMap2D->GetZaxis()->SetRangeUser(0,2);
		hCalibrationMap2D->Draw("colz");
		c3->Print("correction_map.pdf");

		std::cout << "Step 6" << std::endl;
		TCanvas *c4 = new TCanvas("c4", "c4", 200, 10, 750, 1000);
		c4->Divide(1,4);
		c4->cd(1); hHitMap1D->Draw("colz"); //hHitMap1D->Write();
		c4->cd(2); hRawDQdxMap1D->Draw("colz"); //hRawDQdxMap1D->Write();
		c4->cd(4); hAvgDQdxMap1D->Draw("colz"); //hMedianDQdxMap1D->Write();
		c4->Update();
		c4->Print("hitmap2.pdf");

		std::cout << "Step 7" << std::endl;
		TCanvas *c5 = new TCanvas("c5","c5", 200, 10, 750, 300);
		c5->cd();
		hCalibrationMap1D->Draw();
		TCanvas *c6 = new TCanvas("c6","c6", 200, 10, 750, 300);
		c6->cd();
		hCalibrationMap3D->SetLineWidth(0);
		hCalibrationMap3D->Draw("BOX2 Z"); //hCalibrationMap3D->Write();
		hCalibrationMap3D->Write();

		std::cout << "Step 8" << std::endl;
		TCanvas *c7 = new TCanvas("c7", "c7", 200, 10, 750, 1000);
		c7->Divide(1,3);
		hCorrecteddQdx->GetZaxis()->SetRangeUser(0,60);
		hCorrecteddQdx_m->GetZaxis()->SetRangeUser(0,30);
		hCorrecteddQdx_p->GetZaxis()->SetRangeUser(0,30);
		c7->cd(1); hCorrecteddQdx->Draw("colz"); 	 //hCorrecteddQdx->Write();
		c7->cd(2); hCorrecteddQdx_m->Draw("colz"); //hCorrecteddQdx_m->Write();
		c7->cd(3); hCorrecteddQdx_p->Draw("colz"); //hCorrecteddQdx_p->Write();
		c7->Print("Corrected_dQdx.pdf");

		std::cout << "Step 9" << std::endl;
		TCanvas *c8 = new TCanvas("c8","c8", 200, 10, 750, 300);
		c8->cd();
		hCorrectionMap3D->SetLineWidth(0);
		hAvgDQdxMap2D->GetZaxis()->SetRangeUser(0,20); //hAvgDQdxMap2D->Draw("colz");

		hCorrectionMap3D->Draw("BOX2 Z"); hCorrectionMap3D->Write();


		std::cout << "Hello." << std::endl;

		return true;
	}

	void Cal_Test::LoadSelectionFiles(std::string file_loc) {
		/*
			Takes in NuMuVertexVariables tree (FinalVertexVariables, DLLEE_Vertex_Variables) 
			and selects all events that pass associated cuts into cosmic or proton like tracks
			by (run, subrun, event, vtxid) to be used in AllVertexFill.
			--- needs to be updates for DL Merged files ---
		*/
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

	    //NuMuVertexVariable->SetBranchAddress("InFiducial",&InFiducial_sel);
	    NuMuVertexVariable->SetBranchAddress("Reco_goodness_v",&goodReco_sel);
	    NuMuVertexVariable->SetBranchAddress("NtracksReco",&N5cmTracks_sel);
		NuMuVertexVariable->SetBranchAddress("Length_v",&TrackLength_v_sel);
		NuMuVertexVariable->SetBranchAddress("possibleCosmic",&cosmicLL_sel);
		//NuMuVertexVariable->SetBranchAddress("PassCuts",&PassCuts_sel);


	    for (Long64_t i = 0; i < NuMuVertexVariable->GetEntries();i++) { 
	    	NuMuVertexVariable->GetEntry(i);
	    	//std::cout << "pass" << std::endl;
	    	if (goodReco_sel == 0 || N5cmTracks_sel != 2 || TrackLength_v_sel->size() != 2) continue;
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
	    	if (cosmicLL_sel < 3 && TrackLength_v_sel->at(0) > 35 && TrackLength_v_sel->at(1) > 35) { //} && PassCuts_sel = 1) {
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

	void Cal_Test::SelectedVertex(int run, int subrun, int event, int vtx_id) {
		/*
			Called by AllVertexFill
			Just determines if the current event is in SelectEvtID_cosmic or 
			SelectEvtID_proton and changes selectedProton/selectedCosmic value to true.
			Also contains selectedCosmic_during, _after, _before to look at differences
			in the EXT from before, during, and after the beam taking window in run1
			(this functionality is really only used in "add_here")
		*/

        selectedProton = false;
        selectedCosmic = false;
        for(size_t i = 0;i<SelectEvtID_cosmic.size();i++){
            if(run!= SelectEvtID_cosmic[i][0])continue;
            if(subrun!=SelectEvtID_cosmic[i][1])continue;
            if(event!=SelectEvtID_cosmic[i][2])continue;
            if(vtx_id!=SelectEvtID_cosmic[i][3])continue;
            if((run >= 5119 && run <= 5600) || (run >= 5750 && run <= 5955) ) {selectedCosmic_during=true;}
            if( run > 5955) {selectedCosmic_after=true;}
            if(run < 5119) {selectedCosmic_before=true;}
            selectedCosmic=true;
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

	void Cal_Test::AllVertexFill(storage_manager * mgr) {

		/*
			Saves selected larlite tracks from "trackReco_sceadded" tree to AllVertex_{whatever} 
		*/

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
            selectedCosmic_during = false;
            selectedCosmic_after = false;
            selectedCosmic_before = false;
            selectedCosmic = false;

        	SelectedVertex(run, subrun,event,vtx_id);

            //if(!selectedCosmic)continue;
            if(!selectedProton && !selectedCosmic)continue;

            if(thisVertex.size()!=0)thisVertex.clear();
            double sumLength = 0;
            for(auto const& trk_index : vtx_to_trk[vertex_index]) {
                thisVertex.push_back( (*ev_trk)[trk_index]);
                sumLength+=(*ev_trk)[trk_index].Length();
            }
            if(selectedCosmic_during){AllVertex_m_during.push_back(thisVertex);}
            if(selectedCosmic_after){AllVertex_m_after.push_back(thisVertex);}
            if(selectedProton){AllVertex_p.push_back(thisVertex);}
            if(selectedCosmic){
            	AllVertex_m.push_back(thisVertex);
            	std::vector<int> info {};
            	info.push_back(run);
            	info.push_back(subrun);
            	info.push_back(event);
            	AllVertex_run.push_back(info); // Used in RunInfo
            }
            if(selectedCosmic_before){AllVertex_m_before.push_back(thisVertex);}

    	}// vtx_id
    	return;
	}

	void Cal_Test::LoadDEdxSplines(std::string file_loc) {
		/*
		Loads in dE/dx splines from LArCV folders.
		*/

		std::cout << "Load dEdx Splines." << std::endl;
		TFile *fSplines = TFile::Open(Form("%s",file_loc.c_str()),"READ");
		sProtonRange2dEdx = (TSpline3*)fSplines->Get("sProtonRange2dEdx");
		sMuonRange2dEdx = (TSpline3*)fSplines->Get("sMuonRange2dEdx");
		sProtonRange2T = (TSpline3*)fSplines->Get("sProtonRange2T");
	}

	void Cal_Test::LoadCalibrationFile(std::string file_loc) {
		/*
		Loads in 3D calibration maps for use in DataGainFactor()
		*/
		TFile *fCalibration = TFile::Open(Form("%s",file_loc.c_str()),"READ");
		if (!(fCalibration->IsOpen())){
			std::cout << "Error: could not open calibration file" << std::endl;
		}
		else std::cout << "Calibration File Opened" << std::endl;

		CorrectionMap3D_Y = (TH3D*)fCalibration->Get("hImageCalibrationMap_02");
		CorrectionMap3D_U = (TH3D*)fCalibration->Get("hImageCalibrationMap_00");
		CorrectionMap3D_V = (TH3D*)fCalibration->Get("hImageCalibrationMap_01");

		if (CorrectionMap3D_Y==nullptr) std::cout << "Error, no Correction Map loaded" << std::endl;
	}



	void Cal_Test::MakeDQdxMaps() {

		//TFile f("dQdx.root","RECREATE");
		//TCanvas *c1 = new TCanvas("c1", "c1", 200, 10, 4000, 10000);
		//c1->Divide(1,3);

		hUncorrecteddQdx = new TH2D("hUncorrecteddQdx","Uncorrected, Uncalibrated (All Tracks);Residual Length (cm);Uncalibrated dQdx",400,0,400,300,0,300);
		hUncorrecteddQdx_p = new TH2D("hUncorrecteddQdx_p","Uncorrected, Uncalibrated (Protons);Residual Length (cm);Uncalibrated dQdx",400,0,400,300,0,300);
		hUncorrecteddQdx_m = new TH2D("hUncorrecteddQdx_m","Uncorrected, Uncalibrated (Muons);Residual Length (cm);Uncalibrated dQdx",400,0,400,300,0,300);

		double avg_dQdx[2] {0,0};
		for (size_t ivert = 0; ivert < AllVertex_p.size(); ivert++) {
			for (size_t itrack = 0; itrack < AllVertex_p[ivert].size(); itrack++) {
				auto curr_vert = AllVertex_p[ivert][itrack];
				for (size_t iloc = 0; iloc < curr_vert.NumberTrajectoryPoints(); iloc++) {
					double raw_dQdx = curr_vert.DQdxAtPoint(iloc,plane);
					TVector3 point = curr_vert.LocationAtPoint(iloc);
					double reslength = curr_vert.Length(iloc);
					if (reslength < 5) continue;
					if (raw_dQdx <= 0 || raw_dQdx == 111) continue;
					avg_dQdx[itrack] += raw_dQdx;
					hUncorrecteddQdx->Fill(reslength, raw_dQdx);
				}
				avg_dQdx[itrack] /= AllVertex_p[ivert][itrack].NumberTrajectoryPoints();
			}
			int kProton = 0;
			int kMuon = 1;
			if (AllVertex_p[ivert][kProton].Length() < 5 || AllVertex_p[ivert][kMuon].Length() < 5) continue;
			if (AllVertex_p[ivert][kProton].Length() > AllVertex_p[ivert][kMuon].Length()) continue;
			if (avg_dQdx[kProton] < avg_dQdx[kMuon]) continue;
			for (size_t iloc = 0; iloc < AllVertex_p[ivert][kProton].NumberTrajectoryPoints(); iloc++){
				auto curr_vert = AllVertex_p[ivert][kProton];
				double raw_dQdx = curr_vert.DQdxAtPoint(iloc,plane);
				TVector3 point = curr_vert.LocationAtPoint(iloc);
				double reslength = curr_vert.Length(iloc);
				if (reslength < 2) continue;
				if (raw_dQdx <= 0 || raw_dQdx == 111) continue;
				hUncorrecteddQdx_p->Fill(reslength, raw_dQdx);
			}
			for (size_t iloc = 0; iloc < AllVertex_p[ivert][kMuon].NumberTrajectoryPoints(); iloc++){
				auto curr_vert = AllVertex_p[ivert][kMuon];
				double raw_dQdx = curr_vert.DQdxAtPoint(iloc,plane);
				TVector3 point = curr_vert.LocationAtPoint(iloc);
				double reslength = curr_vert.Length(iloc);
				if (reslength < 2) continue;
				if (raw_dQdx <= 0 || raw_dQdx == 111) continue;
				hUncorrecteddQdx_m->Fill(reslength, raw_dQdx);
			}
		}
		std::cout << hUncorrecteddQdx_p->GetEntries() << std::endl;
		std::cout << hUncorrecteddQdx_m->GetEntries() << std::endl;

		// c1->cd(1); hUncorrecteddQdx->Draw("colz");
		// c1->cd(2); hUncorrecteddQdx_m->Draw("colz");
		// c1->cd(3); hUncorrecteddQdx_p->Draw("colz"); //hUncorrecteddQdx_m->Write();
		// c1->Update();
		// c1->Print("graphs.pdf");
		// c2->Print("graphs2.pdf");
		// hUncorrecteddQdx->Write();
		// hUncorrecteddQdx_p->Write();
		// hUncorrecteddQdx_m->Write();	


			// c1->Print("graphs.pdf");
			// c2->Print("graphs2.pdf");
			// c1->cd(1); hUncorrecteddQdx->DrawClone("colz");
			// c1->cd(2); hUncorrecteddQdx_p->DrawClone("colz");
			// c1->cd(3); hUncorrecteddQdx_m->DrawClone("colz");
			// c1->Update();
			// c2->cd(1); ProjX_p->DrawClone();
			// c2->cd(2); ProjY_p->DrawClone();
			// c2->cd(3); ProjX_m->Draw();
			// c2->cd(4); ProjY_m->Draw();
			// c2->Update();
			// c1->Print("graphs.pdf");
			// c2->Print("graphs2.pdf");
			// hUncorrecteddQdx->Write();
			// hUncorrecteddQdx_p->Write();	
			// hUncorrecteddQdx_m->Write();
	}

	void Cal_Test::DataGainFactor() {
		// TH1D* hADC[4];
		// TH1D* hADC_Corr[4];
		//hStart = new TH1D("hStart","hStart",200, -116.5, 116.6);
		//hEnd = new TH2D("hEnd","hEnd",200,-116.5,116.5);

		int Num_Stop {0};
		TF1 *fRecombinationExpected = new TF1("fRecombinationExpected",RecombinationFun,0,1000,6);
		fRecombinationExpected->SetParameters(alpha_ex,beta_ex,1,Wion,LAr_Density,E_Field);
		double MIP_dQdx = fRecombinationExpected->Eval(2.104);
		for (size_t iPlane{1}; iPlane < 5; iPlane++) { 
			hADC[iPlane] = new TH1D(Form("hADC_%zu",iPlane-1),Form("hADC_%zu",iPlane-1),200,0,30);
			hADC_Corr[iPlane] = new TH1D(Form("hADC_Corr_%zu",iPlane-1),Form("hADC_Corr_%zu",iPlane-1),200,0,30);
		}

		for (size_t ivert {0}; ivert < AllVertex_p.size(); ivert++) {
			bool FidIn[2];
			std::vector<std::vector<double>> end_points;
			for (size_t itrack {0}; itrack < AllVertex_p[ivert].size(); itrack++) {
				auto curr_vert = AllVertex_p[ivert][itrack];
				double end_X = curr_vert.End().X();
				double end_Y = curr_vert.End().Y();
				double end_Z = curr_vert.End().Z();
				std::vector<double> end_point = {end_X, end_Y, end_Z};
				end_points.push_back(end_point);
			} 
			FidIn[0] = InFiducial(end_points[0][0],end_points[0][1],end_points[0][2]);
			FidIn[1] = InFiducial(end_points[1][0],end_points[1][1],end_points[1][2]);
			if (FidIn[0] == FidIn[1]) continue;
			int top_track, bottom_track;
			if (FidIn[0] == 1) {top_track = 1; bottom_track = 0;}
			else {top_track = 0; bottom_track = 1;}
			if (end_points[top_track][1] < 106.) continue;
			Num_Stop++;
			for (size_t iloc {0}; iloc < AllVertex_p[ivert][top_track].NumberTrajectoryPoints(); iloc++) {
				auto curr_vert = AllVertex_p[ivert][top_track];
				double raw_dQdx_U = curr_vert.DQdxAtPoint(iloc,larlite::geo::kU);
				double raw_dQdx_V = curr_vert.DQdxAtPoint(iloc,larlite::geo::kV);
				double raw_dQdx_Y = curr_vert.DQdxAtPoint(iloc,larlite::geo::kZ);
				double raw_dQdx = (raw_dQdx_U + raw_dQdx_V + raw_dQdx_Y)/3;
				TVector3 point = curr_vert.LocationAtPoint(iloc); 
				int bin = CorrectionMap3D_Y->FindBin(point.X(),point.Y(),point.Z());
				double corr_dQdx_U = raw_dQdx_U * CorrectionMap3D_U->GetBinContent(bin);
				double corr_dQdx_V = raw_dQdx_V * CorrectionMap3D_V->GetBinContent(bin);
				double corr_dQdx_Y = raw_dQdx_Y * CorrectionMap3D_Y->GetBinContent(bin);
				double corr_dQdx = (corr_dQdx_U + corr_dQdx_V + corr_dQdx_Y)/3;
				if (raw_dQdx_U > 0 && raw_dQdx_U != 37. && corr_dQdx_U > 0) {
					hADC[1]->Fill(raw_dQdx_U); hADC_Corr[1]->Fill(corr_dQdx_U);
				}
				if (raw_dQdx_V > 0 && raw_dQdx_V != 37. && corr_dQdx_V > 0) {
					hADC[2]->Fill(raw_dQdx_V); hADC_Corr[2]->Fill(corr_dQdx_V);
				}
				if (raw_dQdx_Y > 0 && raw_dQdx_Y != 37. && corr_dQdx_Y > 0) {
					hADC[3]->Fill(raw_dQdx_Y); hADC_Corr[3]->Fill(corr_dQdx_Y);
				}
				if (raw_dQdx_U <= 0 || raw_dQdx_V <= 0 || raw_dQdx_Y <= 0 || raw_dQdx_U == 37. || raw_dQdx_V == 37. || raw_dQdx_Y == 37. || corr_dQdx_V <= 0 || corr_dQdx_U <= 0 || corr_dQdx_Y <= 0) continue;
				hADC[4]->Fill(raw_dQdx); hADC_Corr[4]->Fill(corr_dQdx);
			}
			for (size_t iloc {0}; iloc < AllVertex_p[ivert][bottom_track].NumberTrajectoryPoints(); iloc++) {
				auto curr_vert = AllVertex_p[ivert][bottom_track];
				double reslength = curr_vert.Length(iloc);
				if (reslength < 40) continue;
				double raw_dQdx_U = curr_vert.DQdxAtPoint(iloc,larlite::geo::kU);
				double raw_dQdx_V = curr_vert.DQdxAtPoint(iloc,larlite::geo::kV);
				double raw_dQdx_Y = curr_vert.DQdxAtPoint(iloc,larlite::geo::kZ);
				double raw_dQdx = (raw_dQdx_U + raw_dQdx_V + raw_dQdx_Y)/3;
				TVector3 point = curr_vert.LocationAtPoint(iloc);
				int bin = CorrectionMap3D_Y->FindBin(point.X(),point.Y(),point.Z());
				double corr_dQdx_U = raw_dQdx_U * CorrectionMap3D_U->GetBinContent(bin);
				double corr_dQdx_V = raw_dQdx_V * CorrectionMap3D_V->GetBinContent(bin);
				double corr_dQdx_Y = raw_dQdx_Y * CorrectionMap3D_Y->GetBinContent(bin);
				double corr_dQdx = (corr_dQdx_U + corr_dQdx_V + corr_dQdx_Y)/3;
				if (raw_dQdx_U > 0 && raw_dQdx_U != 37. && corr_dQdx_U > 0) {
					hADC[1]->Fill(raw_dQdx_U); hADC_Corr[1]->Fill(corr_dQdx_U);
				}
				if (raw_dQdx_V > 0 && raw_dQdx_V != 37. && corr_dQdx_V > 0) {
					hADC[2]->Fill(raw_dQdx_V); hADC_Corr[2]->Fill(corr_dQdx_V);
				}
				if (raw_dQdx_Y > 0 && raw_dQdx_Y != 37. && corr_dQdx_Y > 0) {
					hADC[3]->Fill(raw_dQdx_Y); hADC_Corr[3]->Fill(corr_dQdx_Y);
				}
				if (raw_dQdx_U <= 0 || raw_dQdx_V <= 0 || raw_dQdx_Y <= 0 || raw_dQdx_U == 37. || raw_dQdx_V == 37. || raw_dQdx_Y == 37. || corr_dQdx_V <= 0 || corr_dQdx_U <= 0 || corr_dQdx_Y <= 0) continue;
				hADC[4]->Fill(raw_dQdx); hADC_Corr[4]->Fill(corr_dQdx);
			}
		}

		// for (size_t iPlane{1}; iPlane < 5; iPlane++) {
		// 	double max = hADC[iPlane]->GetBinCenter(hADC[iPlane]->GetMaximumBin());
		// 	double rms = 1 * hADC[iPlane]->GetRMS();
		// 	double max_Corr = hADC_Corr[iPlane]->GetBinCenter(hADC_Corr[iPlane]->GetMaximumBin());
		// 	double rms_Corr = hADC_Corr[iPlane]->GetRMS();
		// 	f1[iPlane] = new TF1(Form("fgaus_%zu",iPlane),"gaus(0)",max - rms,max + rms);
		// 	hADC[iPlane]->Fit(Form("fgaus_%zu",iPlane),"q","", max - rms, max + rms);
		// 	f1_corr[iPlane] = new TF1(Form("fgaus_Corr_%zu",iPlane),"gaus(0)",max_Corr - rms_Corr,max_Corr + rms_Corr);
		// 	hADC_Corr[iPlane]->Fit(Form("fgaus_Corr_%zu",iPlane),"q","", max_Corr - rms_Corr, max_Corr + rms_Corr);
		// 	std::cout << "Uncorrected " << iPlane << " Plane: " << MIP_dQdx/f1[iPlane]->GetParameter(1) << " +/- " << (MIP_dQdx/f1[iPlane]->GetParameter(1))/f1[iPlane]->GetParameter(1) * f1[iPlane]->GetParError(1) << std::endl;
		// 	std::cout << "Corrected " << iPlane << " Plane: " << MIP_dQdx/f1_corr[iPlane]->GetParameter(1) << " +/- " << (MIP_dQdx/f1[iPlane]->GetParameter(1))/f1_corr[iPlane]->GetParameter(1) * f1_corr[iPlane]->GetParError(1) << std::endl;
		// }

		TCanvas *cADC = new TCanvas("cADC","cADC",350,800);
		TCanvas *cADC_Corr = new TCanvas("cADC_Corr","cADC_Corr",350,800);
		cADC->Divide(1,4);
		cADC_Corr->Divide(1,4);
		for (size_t iPlane {1}; iPlane < 5; iPlane++) {
			cADC->cd(iPlane);
			hADC[iPlane]->Draw();
			cADC_Corr->cd(iPlane);
			hADC_Corr[iPlane]->Draw();

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
		hADC_Corr[4]->SetLineColor(1);
		hADC_Corr[4]->SetLineWidth(4);
		hADC_Corr[4]->GetXaxis()->SetTitle("dQ/dx (ADC/cm)");
		hADC_Corr[4]->GetYaxis()->SetTitle("Entries");
		hADC_Corr[4]->GetXaxis()->SetTitleOffset(1);
		hADC_Corr[4]->GetXaxis()->SetTitleSize(0.06);
		hADC_Corr[4]->GetXaxis()->SetLabelSize(0.05);
		hADC_Corr[4]->GetYaxis()->SetTitleOffset(1);
		hADC_Corr[4]->GetYaxis()->SetTitleSize(0.06);
		hADC_Corr[4]->GetYaxis()->SetLabelSize(0.05);
		hADC_Corr[4]->Draw();
		hADC_Corr[1]->Draw("same");
		hADC_Corr[3]->Draw("same");
		hADC_Corr[2]->Draw("same");
		cADC_Corr->Update();
		// TLine *All_Corr = new TLine(9.04167,0,9.04167,cADC_Corr->GetUymax());
		// TLine *Y_Corr = new TLine(8.875,0,8.875,cADC_Corr->GetUymax());
		// All_Corr->SetLineColor(1);
		// All_Corr->SetLineStyle(6);
		// All_Corr->SetLineWidth(3);
		// Y_Corr->SetLineColor(2);
		// Y_Corr->SetLineStyle(6);
		// Y_Corr->SetLineWidth(3);
		// Y_Corr->Draw("same");
		// All_Corr->Draw("same");

		TLegend *legADC_Corr = new TLegend(0.55,0.55,0.87,0.87);
    	legADC_Corr->AddEntry(hADC_Corr[1],"Plane 0","L");
		legADC_Corr->AddEntry(hADC_Corr[2],"Plane 1","L");
		legADC_Corr->AddEntry(hADC_Corr[3],"Plane 2","L");
		legADC_Corr->AddEntry(hADC_Corr[4],"All Planes","L");
		//legADC_Corr->AddEntry(Y_Corr,"Most Probable Value (Plane 0)","L");
		//legADC_Corr->AddEntry(All_Corr,"Most Probable Value (All Planes)","L");
		legADC_Corr->Draw("same");
		}



		//cADC->SaveAs("ADC_new.png");

		//std::cout << "Number of Stopping Tracks: " << Num_Stop << std::endl;


		return;
	}

	void Cal_Test::ImageCorrectionMap() {

		/*
		Takes selected cosmics and creates calibration maps in x,y,z
		*/	

		//creating some variables, edit Xbin, Ybin, Zbin to change map binning size
		int Xbin{50}, Ybin{50}, Zbin{250};
		larlite::geo::View_t views[3];
		views[0] = larlite::geo::kU;
		views[1] = larlite::geo::kV;
		views[2] = larlite::geo::kZ;

		std::vector<std::vector<std::vector<double>>> fillBin;
		std::vector<std::vector<std::vector<double>>> fillBin_ZY;
		std::vector<std::vector<double>> totalCount;

		fillBin.resize(3);
		fillBin_ZY.resize(3);
		totalCount.resize(3);


		
		hXBins = new TH1D("hXBins","hXBins",Xbin,xmin,xmax);

		for (int planes {0}; planes < 3; planes++) {
			h2DImageHitMap[planes+1] = new TH2D(Form("h2DImageHitMap_%02d",planes),Form("h2DImageHitMap_%02d",planes),Ybin,ymin,ymax,Zbin,zmin,zmax);
			h2DImageHitMap_ZY[planes+1] = new TH2D(Form("h2DImageHitMap_ZY_%02d",planes),Form("h2DImageHitMap_ZY_%02d",planes),Zbin,zmin,zmax,Ybin,ymin,ymax);
			hImageCalibrationMap[planes+1] = new TH3D(Form("hImageCalibrationMap_%02d",planes),Form("hImageCalibrationMap_%02d",planes),Xbin,xmin,xmax,Ybin,ymin,ymax,Zbin,zmin,zmax);
			h2DImageMap_ZY[planes+1] = new TH2D(Form("h2DImageMap_ZY_%02d",planes),Form("plane_%02d",planes),Zbin,zmin,zmax,Ybin,ymin,ymax);
			h2DImageMap[planes+1] = new TH2D(Form("h2DImageMap_%02d",planes),Form("h2DImageMap_%02d",planes),Ybin,ymin,ymax,Zbin,zmin,zmax);
			h2DRawImageMap_ZY[planes+1] = new TH2D(Form("h2DRawImageMap_ZY_%02d",planes),Form("h2DRawImageMap_ZY_%02d",planes),Zbin,zmin,zmax,Ybin,ymin,ymax);
			h2DRawImageMap[planes+1] = new TH2D(Form("h2DRawImageMap_%02d",planes),Form("h2DRawImageMap_%02d",planes),Ybin,ymin,ymax,Zbin,zmin,zmax);
			h2DImageCalMap[planes+1] = new TH2D(Form("h2DImageCalMap_%02d",planes),Form("h2DImageCalMap_%02d",planes),Ybin,ymin,ymax,Zbin,zmin,zmax);
			h2DImageCalMap_ZY[planes+1] = new TH2D(Form("h2DImageCalMap_ZY_%02d",planes),Form("plane_%02d",planes),Zbin,zmin,zmax,Ybin,ymin,ymax);
			h2DTotalCount[planes+1] = new TH1D(Form("h2DTotalCount_%02d",planes),"",50,0,35);
			fillBin[planes].resize(h2DImageHitMap[planes+1]->GetSize());
			fillBin_ZY[planes].resize(h2DImageHitMap_ZY[planes+1]->GetSize());
		}	
		
		//Looping through all the selected tracks and adding there dqdx value to 3d location
		//note: simply uses two histograms, a total charge and a hit number, average dqdx for a bin
		// is the total charge divided by the hit number

		for (size_t ivert{0}; ivert < AllVertex_m.size(); ivert++) {
			for (size_t itrack{0}; itrack < AllVertex_m[ivert].size(); itrack++) {
				auto curr_vert = AllVertex_m[ivert][itrack];
				for (size_t iloc {0}; iloc < curr_vert.NumberTrajectoryPoints(); iloc++) {
					//std::cout << "Debug beta" << std::endl;
					TVector3 point = curr_vert.LocationAtPoint(iloc);
					//std::cout << "Debug alpha" << std::endl;
					int bin = h2DImageHitMap[1]->FindBin(point.Y(), point.Z());
					int bin_ZY = h2DImageHitMap_ZY[1]->FindBin(point.Z(), point.Y());
					for (int planes{0}; planes < 3; planes++) {
						if (curr_vert.DQdxAtPoint(iloc,views[planes]) > 1000
						|| curr_vert.DQdxAtPoint(iloc, views[planes]) < 3
						|| curr_vert.Length(iloc) < 20) continue;
						int curr = h2DImageHitMap[planes+1]->GetBinContent(bin);
						int curr_ZY = h2DImageHitMap_ZY[planes+1]->GetBinContent(bin_ZY);
						h2DImageHitMap[planes+1]->SetBinContent(bin, curr+1);
						h2DImageHitMap_ZY[planes+1]->SetBinContent(bin_ZY, curr+1);
						//std::cout << "Debug A" << std::endl;
						fillBin[planes][bin].push_back(curr_vert.DQdxAtPoint(iloc,views[planes]));
						fillBin_ZY[planes][bin_ZY].push_back(curr_vert.DQdxAtPoint(iloc,views[planes]));
						h2DRawImageMap[planes+1]->Fill(point.Y(),point.Z(),curr_vert.DQdxAtPoint(iloc,views[planes]));
						h2DRawImageMap_ZY[planes+1]->Fill(point.Z(),point.Y(),curr_vert.DQdxAtPoint(iloc,views[planes]));
						//std::cout << "Debug B" << std::endl;
						totalCount[planes].push_back(curr_vert.DQdxAtPoint(iloc,views[planes]));
						//std::cout << "Debug C" << std::endl;
						h2DTotalCount[planes+1]->Fill(curr_vert.DQdxAtPoint(iloc,views[planes]));
					}
				}
			} 
		}

		std::vector<double> f1Mean;
		f1Mean.resize(3);
		TCanvas *cADCfit = new TCanvas("cADCfit","cADCfit",250,800);
		cADCfit->Divide(1,3);

		// plotting stuff
		for (int planes{0}; planes < 3; planes++) {
			double rms = h2DTotalCount[planes+1]->GetRMS();
			double mean = h2DTotalCount[planes+1]->GetMean();
			double max = h2DTotalCount[planes+1]->GetMaximumBin();
			double maxVal = h2DTotalCount[planes+1]->GetBinCenter(max);
			std::cout << "Plane " << planes << "Median ADC Value: " << maxVal << "." << std::endl;
			//TF1 *f1 = new TF1("f1","gaus");
			//h2DTotalCount[planes+1]->Fit("f1","","",mean-rms,mean+rms);
			//std::cout << "Debug A" << std::endl;
			f1Mean[planes] = maxVal;
			//std::cout << "Debug B" << std::endl;
			cADCfit->cd(planes+1);
			//std::cout << "Debug C" << std::endl;
			h2DTotalCount[planes+1]->Draw();
			//f1->Draw("same");
		}


		TCanvas *cADC = new TCanvas("cADC","",500,400);
		cADC->cd();
		cADC->SetLeftMargin(0.13);
		cADC->SetBottomMargin(0.13);
		h2DTotalCount[1]->SetLineColor(4);
		h2DTotalCount[1]->SetLineWidth(4);
		h2DTotalCount[2]->SetLineWidth(4);
		h2DTotalCount[3]->SetLineWidth(4);
		h2DTotalCount[2]->SetLineColor(3);
		h2DTotalCount[3]->SetLineColor(2);
		h2DTotalCount[1]->GetXaxis()->SetTitle("dQ/dx (ADC/cm)");
		h2DTotalCount[1]->GetYaxis()->SetTitle("Entries");
		h2DTotalCount[1]->GetXaxis()->SetTitleOffset(1);
		h2DTotalCount[1]->GetXaxis()->SetTitleSize(0.06);
		h2DTotalCount[1]->GetXaxis()->SetLabelSize(0.05);
		h2DTotalCount[1]->GetYaxis()->SetTitleOffset(1);
		h2DTotalCount[1]->GetYaxis()->SetTitleSize(0.06);
		h2DTotalCount[1]->GetYaxis()->SetLabelSize(0.05);
		h2DTotalCount[1]->Draw();
		h2DTotalCount[2]->Draw("same");
		h2DTotalCount[3]->Draw("same");
		cADC->Update();
		TLine *Line0 = new TLine(f1Mean[0],0,f1Mean[0],cADC->GetUymax());
		TLine *Line1 = new TLine(f1Mean[1],0,f1Mean[1],cADC->GetUymax());
		TLine *Line2 = new TLine(f1Mean[2],0,f1Mean[2],cADC->GetUymax());
		Line0->SetLineStyle(9);
		Line0->SetLineWidth(3);
		Line0->SetLineColor(4);
		Line1->SetLineStyle(5);
		Line1->SetLineWidth(3);
		Line1->SetLineColor(3);
		Line2->SetLineStyle(6);
		Line2->SetLineWidth(3);
		Line2->SetLineColor(2);

		Line0->Draw("same");
		Line1->Draw("same");
		Line2->Draw("same");


		TLegend *legADC = new TLegend(0.55,0.55,0.87,0.87);
    	legADC->AddEntry(h2DTotalCount[1],"Plane 0","L");
		legADC->AddEntry(h2DTotalCount[2],"Plane 1","L");
		legADC->AddEntry(h2DTotalCount[3],"Plane 2","L");
		legADC->AddEntry(Line0,"Average MIP dQ/dx (Plane 0)","L");
		legADC->AddEntry(Line1,"Average MIP dQ/dx (Plane 1)","L");
		legADC->AddEntry(Line2,"Average MIP dQ/dx (Plane 2)","L");
		legADC->Draw("same");

		//cADC->SaveAs("TotalADC_3Planes.pdf");

		TCanvas *cAvgCalMap = new TCanvas("cAvgCalMap","cAvgCalMap",400,800);
		TCanvas *cCalMap = new TCanvas("cCalMap","cCalMap",400,800);

		cCalMap->Divide(1,3);
		cAvgCalMap->Divide(1,3);

		// making the calibration map 

		for (size_t idx {0}; idx < h2DImageHitMap[1]->GetSize(); idx++) {
			for (int planes{0}; planes < 3; planes++) {
				if (h2DImageHitMap[planes+1]->GetBinContent(idx) < 5) {
					h2DImageMap[planes+1]->SetBinContent(idx,0);
					h2DImageCalMap[planes+1]->SetBinContent(idx,1);
				}
				else {
					double N = h2DImageHitMap[planes+1]->GetBinContent(idx);
					double tot = h2DRawImageMap[planes+1]->GetBinContent(idx);
					double val = tot/N;
					h2DImageMap[planes+1]->SetBinContent(idx,val);
					h2DImageCalMap[planes+1]->SetBinContent(idx,f1Mean[planes] / val);

				}
			}
		}

		// more plotting 

		TCanvas *plane0 = new TCanvas("plane0","plane0",700,450);
		TCanvas *plane1 = new TCanvas("plane1","plane1",700,450);
		TCanvas *plane2 = new TCanvas("plane2","plane2",700,450);
		TCanvas *plane00 = new TCanvas("plane00","plane00",700,450);
		TCanvas *plane11 = new TCanvas("plane11","plane11",700,450);
		TCanvas *plane22 = new TCanvas("plane22","plane22",700,450);

		for (size_t idx {0}; idx < h2DImageHitMap_ZY[1]->GetSize(); idx++) {
			for (int planes{0}; planes < 3; planes++) {
				if (h2DImageHitMap_ZY[planes+1]->GetBinContent(idx) < 5) {
					h2DImageMap_ZY[planes+1]->SetBinContent(idx,0);
					h2DImageCalMap_ZY[planes+1]->SetBinContent(idx,0);
				}
				else {
					double N = h2DImageHitMap_ZY[planes+1]->GetBinContent(idx);
					double tot = h2DRawImageMap_ZY[planes+1]->GetBinContent(idx);
					double val = tot/N;
					h2DImageMap_ZY[planes+1]->SetBinContent(idx,val);
					h2DImageCalMap_ZY[planes+1]->SetBinContent(idx,f1Mean[planes] / val);

				}
			}
		}
		for (size_t iPlane{1}; iPlane < 4; iPlane++) {
			h2DImageMap_ZY[iPlane]->GetXaxis()->SetTitle("Z Coordinate (cm)");
			h2DImageCalMap_ZY[iPlane]->GetXaxis()->SetTitle("Z Coordinate (cm)");
			h2DImageMap_ZY[iPlane]->GetYaxis()->SetTitle("Y Coordinate (cm)");
			h2DImageCalMap_ZY[iPlane]->GetYaxis()->SetTitle("Y Coordinate (cm)");
			h2DImageCalMap_ZY[iPlane]->GetZaxis()->SetRangeUser(0,3);
			h2DImageMap_ZY[iPlane]->GetXaxis()->SetTitleOffset(0.75);
			h2DImageMap_ZY[iPlane]->GetXaxis()->SetTitleSize(0.06);
			h2DImageMap_ZY[iPlane]->GetXaxis()->SetLabelSize(0.05);
			h2DImageMap_ZY[iPlane]->GetYaxis()->SetTitleOffset(0.75);
			h2DImageMap_ZY[iPlane]->GetYaxis()->SetTitleSize(0.06);
			h2DImageMap_ZY[iPlane]->GetYaxis()->SetLabelSize(0.05);
			h2DImageCalMap_ZY[iPlane]->GetXaxis()->SetTitleOffset(0.75);
			h2DImageCalMap_ZY[iPlane]->GetXaxis()->SetTitleSize(0.06);
			h2DImageCalMap_ZY[iPlane]->GetXaxis()->SetLabelSize(0.05);
			h2DImageCalMap_ZY[iPlane]->GetYaxis()->SetTitleOffset(0.75);
			h2DImageCalMap_ZY[iPlane]->GetYaxis()->SetTitleSize(0.06);
			h2DImageCalMap_ZY[iPlane]->GetYaxis()->SetLabelSize(0.05);

		}

		plane0->cd(); plane0->SetLeftMargin(0.13); plane0->SetBottomMargin(0.13); h2DImageMap_ZY[1]->Draw("colz");
		plane00->cd(); plane00->SetLeftMargin(0.13); plane00->SetBottomMargin(0.13); h2DImageCalMap_ZY[1]->Draw("colz");
		plane1->cd(); plane1->SetLeftMargin(0.13); plane1->SetBottomMargin(0.13); h2DImageMap_ZY[2]->Draw("colz");
		plane11->cd(); plane11->SetLeftMargin(0.13); plane11->SetBottomMargin(0.13); h2DImageCalMap_ZY[2]->Draw("colz");
		plane2->cd(); plane2->SetLeftMargin(0.13); plane2->SetBottomMargin(0.13); h2DImageMap_ZY[3]->Draw("colz");
		plane22->cd(); plane22->SetLeftMargin(0.13); plane22->SetBottomMargin(0.13); h2DImageCalMap_ZY[3]->Draw("colz");

		// plane0->SaveAs(Form("%s.pdf",plane0->GetName()));
		// plane00->SaveAs(Form("%s.pdf",plane00->GetName()));
		// plane1->SaveAs(Form("%s.pdf",plane1->GetName()));
		// plane11->SaveAs(Form("%s.pdf",plane11->GetName()));
		// plane2->SaveAs(Form("%s.pdf",plane2->GetName()));
		// plane22->SaveAs(Form("%s.pdf",plane22->GetName()));


		for (int planes{0}; planes < 3; planes++) {
			cAvgCalMap->cd(planes+1);
			h2DImageMap[planes+1]->Draw("colz");
		}
		for (int planes{0}; planes < 3; planes++) {
			cCalMap->cd(planes+1);
			h2DImageCalMap[planes+1]->GetZaxis()->SetRangeUser(0.,2.);
			h2DImageCalMap[planes+1]->Draw("colz");
		}

		for ( size_t ybin {0}; ybin < Ybin; ybin++ ) {
			for ( size_t zbin {0}; zbin < Zbin; zbin++ ) {
				for ( size_t xbin {0}; xbin < Xbin; xbin++ ) {
					int bin = hImageCalibrationMap[1]->GetBin(xbin,ybin,zbin);
					for (int planes {0}; planes < 3; planes++) {
						double cal = h2DImageCalMap[planes+1]->GetBinContent(ybin,zbin);
						hImageCalibrationMap[planes+1]->SetBinContent(bin,cal);
					}
				}
			}
		}

		// saving calibration maps
		TFile f("CalibrationMaps.root","RECREATE");
		hImageCalibrationMap[1]->Write();
		hImageCalibrationMap[2]->Write();
		hImageCalibrationMap[3]->Write();


	}


	void Cal_Test::Make2DCorrectionMap() {
		larlite::geo::View_t plane = larlite::geo::kZ;

		// Use a Vector of Vectors of doubles to store bin values, so each vector of doubles cooresponds to a bin, then pushback additonal values to that vector
		// can then find the mean or avg of that bin as needed
		hHitMap2D = new TH2D("hHitMap2D","hHitMap2D", Nz, zmin, zmax, Ny, ymin, ymax);
		hRawDQdxMap2D = new TH2D("RawDQdxMap2D","RawDQdxMap2D", Nz, zmin, zmax, Ny, ymin, ymax);
		hAvgDQdxMap2D = new TH2D("AvgDQdxMap2D","AvgDQdxMap2D", Nz, zmin, zmax, Ny, ymin, ymax);
		hMedianDQdxMap2D = new TH2D("MedianDQdxMap2D","MedianDQdxMap2D", Nz, zmin, zmax, Ny, ymin, ymax);
		hCalibrationMap2D = new TH2D("CalibrationMap2D","CalibrationMap2D", Nz, zmin, zmax, Ny, ymin, ymax);
		hLowHit2D = new TH2D("hLowHit2D","hLowHit2D", Nz, zmin, zmax, Ny, ymin, ymax);
		hTotalCount = new TH1D("TotalADC","TotalADC", 50, 0, 50);
		hRecoFit = new TH2D("RecoFit","RecoFit",60,0,15,50,0,500000);



		binFill.resize(hHitMap2D->GetSize());
		std::vector<double> totalCount {};
	

		for (size_t ivert{0}; ivert < AllVertex_m.size(); ivert++) {
			for (size_t itrack{0}; itrack < AllVertex_m[ivert].size(); itrack++) {
				auto curr_vert = AllVertex_m[ivert][itrack];
				for (size_t iloc {0}; iloc < curr_vert.NumberTrajectoryPoints(); iloc++) {
					if (curr_vert.DQdxAtPoint(iloc,plane) > 111
						|| curr_vert.DQdxAtPoint(iloc, plane) == 0
						|| curr_vert.Length(iloc) < 20) continue;
					TVector3 point = curr_vert.LocationAtPoint(iloc);
					hTotalCount->Fill(curr_vert.DQdxAtPoint(iloc,plane));

				}
			}
		}
		double rms = hTotalCount->GetRMS();
		double mean = hTotalCount->GetMean();
		TF1 *f1 = new TF1("f1","gaus", 6, 18);
		TF1 *f2 = new TF1("f2","gaus");
		hTotalCount->Fit("f1","","", 6, 18);
		//hTotalCount->Fit("f1","","");
		double f1Mean = f1->GetParameter(1);
		double f1Sigma = f1->GetParameter(2);

		TCanvas *fit = new TCanvas("fit","fit",200,10,300,300);
		fit->cd();
		hTotalCount->Draw();
		f1->SetLineWidth(4);
		f1->Draw("same");



		std::cout << f1Mean << " and " << f1Sigma << std::endl;
		minADC = f1Mean - 1.5*f1Sigma;
		maxADC = f1Mean + 1.5*f1Sigma;

		medVal = f1Mean;

		std::cout << "ADC range:" << minADC << " - " << maxADC << std::endl;

		for (size_t ivert{0}; ivert < AllVertex_m.size(); ivert++) {
			for (size_t itrack{0}; itrack < AllVertex_m[ivert].size(); itrack++) {
				auto curr_vert = AllVertex_m[ivert][itrack];
				//std::cout << curr_vert.Run() << std::endl;
				for (size_t iloc {0}; iloc < curr_vert.NumberTrajectoryPoints(); iloc++) {
					if (curr_vert.DQdxAtPoint(iloc,plane) > 1000
						|| curr_vert.DQdxAtPoint(iloc, plane) < 3
						|| curr_vert.Length(iloc) < 20) continue;
					TVector3 point = curr_vert.LocationAtPoint(iloc);
					double res_length = curr_vert.Length(iloc);
					auto bin = hHitMap2D->FindBin(point.Z(), point.Y());
					int curr = hHitMap2D->GetBinContent(bin);
					//hRecoFit->Fill(sMuonRange2dEdx->Eval(res_length),(curr_vert.DQdxAtPoint(iloc,plane)*avgCor));
					hHitMap2D->SetBinContent(bin, curr+1);
					hRawDQdxMap2D->Fill(point.Z(), point.Y(), curr_vert.DQdxAtPoint(iloc, plane));
					binFill[bin].push_back(curr_vert.DQdxAtPoint(iloc,plane));
					totalCount.push_back(curr_vert.DQdxAtPoint(iloc,plane));
					//if (curr_vert.DQdxAtPoint(iloc,plane) < 1) hLowHit2D->Fill(point.Z(),point.Y());
				}
			}
		}

		std::cout << "Hit Map 2D Filled" << std::endl;

		if (hHitMap2D->GetSize() != hRawDQdxMap2D->GetSize()) {
			std::cout << "Histogram's don't have same number of bins!";
		}
		else {
			for (size_t idx {0}; idx < hHitMap2D->GetSize() ; idx++) {
				if (hHitMap2D->GetBinContent(idx) > 5) {
					double med_dQdx {0.0};
					double avg_dQdx = hRawDQdxMap2D->GetBinContent(idx) / hHitMap2D->GetBinContent(idx);
					hAvgDQdxMap2D->SetBinContent(idx, avg_dQdx);
					auto temp = binFill[idx];
					std::sort(temp.begin(),temp.end());
					if (temp.size() % 2 == 1) {
						med_dQdx = temp[temp.size() / 2];
					}
					else {
						med_dQdx = 0.5*(temp[temp.size() / 2 - 1] + temp[temp.size()/2]);
					}
					hMedianDQdxMap2D->SetBinContent(idx, med_dQdx);
				}
				// else {
				// 	hAvgDQdxMap2D->SetBinContent(idx, 0);
				// 	hMedianDQdxMap2D->SetBinContent(idx, 0);

				// }
			}
		}

		std::cout << "Avg DQdx Map 2D Filled" << std::endl;	

		// std::sort(totalCount.begin(), totalCount.end());
		//medVal = 0.0;
		//if (totalCount.size() % 2 == 1) medVal = totalCount[totalCount.size() / 2];
		//else medVal = 0.5*(totalCount[totalCount.size() / 2 - 1] + totalCount[totalCount.size()/2]); 

		std::cout << "Avg Value Found" << std::endl;

		for (size_t idx {0}; idx < hAvgDQdxMap2D->GetSize(); idx++) {
			if (hAvgDQdxMap2D->GetBinContent(idx) != 0) {
				hCalibrationMap2D->SetBinContent(idx, medVal / hAvgDQdxMap2D->GetBinContent(idx));
			}
			//else hCalibrationMap2D->SetBinContent(idx, 1);
		}

		std::cout << "2D Calibration Map Made!" << std::endl;
	}

	void Cal_Test::LookAtX() { 
		h5e19_avg = new TH1D("h5e19","h5e19",Nx,xmin,xmax);
		h5e19_raw = new TH1D("h5e19_raw","h5e19_raw",Nx,xmin,xmax);
		hAfter_avg = new TH1D("hAfter","hAfter",Nx,xmin,xmax);
		hAfter_raw = new TH1D("hAfter_raw","hAfter_raw",Nx,xmin,xmax);
		hBefore_avg = new TH1D("hBefore_avg","hBefore_avg",Nx,xmin,xmax);
		hBefore_raw = new TH1D("hBefore_raw","hBefore_raw",Nx,xmin,xmax);

		std::vector<std::vector<double>> binFill_during {}, binFill_after {}, binFill_before;
		std::vector<double> totalCount_during {}, totalCount_after {}, totalCount_before {};
		binFill_during.resize(h5e19_avg->GetSize());
		binFill_after.resize(h5e19_avg->GetSize());
		binFill_before.resize(h5e19_avg->GetSize());


		for (size_t ivert{0}; ivert < AllVertex_m_during.size(); ivert++) {
			for (size_t itrack{0}; itrack < AllVertex_m_during[ivert].size(); itrack++) {
				auto curr_vert = AllVertex_m_during[ivert][itrack];
				for (size_t iloc {0}; iloc < curr_vert.NumberTrajectoryPoints(); iloc++) {
					if (curr_vert.DQdxAtPoint(iloc,plane) > 1000
						|| curr_vert.DQdxAtPoint(iloc, plane) < 3
						|| curr_vert.Length(iloc) < 20) continue;
					TVector3 point = curr_vert.LocationAtPoint(iloc);
					int binX = h5e19_raw->FindBin(point.X());
					h5e19_raw->Fill(point.X(), curr_vert.DQdxAtPoint(iloc, plane));
					binFill_during[binX].push_back(curr_vert.DQdxAtPoint(iloc, plane));
					totalCount_during.push_back(curr_vert.DQdxAtPoint(iloc,plane));
					//if (curr_vert.DQdxAtPoint(iloc,plane) < 1) hLowHit1D->Fill(point.X());

				}
			}
		}
	
		for (size_t idx {0}; idx < binFill_during.size() ; idx++) {
			if (binFill_during[idx].size() > 5) {
				double avg_dQdx = h5e19_raw->GetBinContent(idx) / binFill_during[idx].size();
				h5e19_avg->SetBinContent(idx, avg_dQdx);
			}
		}

		for (size_t ivert{0}; ivert < AllVertex_m_after.size(); ivert++) {
			for (size_t itrack{0}; itrack < AllVertex_m_after[ivert].size(); itrack++) {
				auto curr_vert = AllVertex_m_after[ivert][itrack];
				for (size_t iloc {0}; iloc < curr_vert.NumberTrajectoryPoints(); iloc++) {
					if (curr_vert.DQdxAtPoint(iloc,plane) > 1000
						|| curr_vert.DQdxAtPoint(iloc, plane) < 3
						|| curr_vert.Length(iloc) < 20) continue;
					TVector3 point = curr_vert.LocationAtPoint(iloc);
					int binX = h5e19_raw->FindBin(point.X());
					hAfter_raw->Fill(point.X(), curr_vert.DQdxAtPoint(iloc, plane));
					binFill_after[binX].push_back(curr_vert.DQdxAtPoint(iloc, plane));
					totalCount_after.push_back(curr_vert.DQdxAtPoint(iloc,plane));
					//if (curr_vert.DQdxAtPoint(iloc,plane) < 1) hLowHit1D->Fill(point.X());

				}
			}
		}
	}

	void Cal_Test::LookAtZ() { 
		h5e19_avg = new TH1D("h5e19","h5e19",Nz,zmin,zmax);
		h5e19_raw = new TH1D("h5e19_raw","h5e19_raw",Nz,zmin,zmax);
		hAfter_avg = new TH1D("hAfter","hAfter",Nz,zmin,zmax);
		hAfter_raw = new TH1D("hAfter_raw","hAfter_raw",Nz,zmin,zmax);
		hBefore_avg = new TH1D("hBefore_avg","hBefore_avg",Nx,xmin,xmax);
		hBefore_raw = new TH1D("hBefore_raw","hBefore_raw",Nx,xmin,xmax);

		std::vector<std::vector<double>> binFill_during {}, binFill_after {}, binFill_before;
		std::vector<double> totalCount_during {}, totalCount_after {}, totalCount_before {};
		binFill_during.resize(h5e19_avg->GetSize());
		binFill_after.resize(h5e19_avg->GetSize());
		binFill_before.resize(h5e19_avg->GetSize());


		for (size_t ivert{0}; ivert < AllVertex_m_during.size(); ivert++) {
			for (size_t itrack{0}; itrack < AllVertex_m_during[ivert].size(); itrack++) {
				auto curr_vert = AllVertex_m_during[ivert][itrack];
				for (size_t iloc {0}; iloc < curr_vert.NumberTrajectoryPoints(); iloc++) {
					if (curr_vert.DQdxAtPoint(iloc,plane) > 1000
						|| curr_vert.DQdxAtPoint(iloc, plane) < 3
						|| curr_vert.Length(iloc) < 20) continue;
					TVector3 point = curr_vert.LocationAtPoint(iloc);
					int binX = h5e19_raw->FindBin(point.Z());
					h5e19_raw->Fill(point.Z(), curr_vert.DQdxAtPoint(iloc, plane));
					binFill_during[binX].push_back(curr_vert.DQdxAtPoint(iloc, plane));
					totalCount_during.push_back(curr_vert.DQdxAtPoint(iloc,plane));
					//if (curr_vert.DQdxAtPoint(iloc,plane) < 1) hLowHit1D->Fill(point.X());

				}
			}
		}
	
		for (size_t idx {0}; idx < binFill_during.size() ; idx++) {
			if (binFill_during[idx].size() > 5) {
				double avg_dQdx = h5e19_raw->GetBinContent(idx) / binFill_during[idx].size();
				h5e19_avg->SetBinContent(idx, avg_dQdx);
			}
		}

		for (size_t ivert{0}; ivert < AllVertex_m_after.size(); ivert++) {
			for (size_t itrack{0}; itrack < AllVertex_m_after[ivert].size(); itrack++) {
				auto curr_vert = AllVertex_m_after[ivert][itrack];
				for (size_t iloc {0}; iloc < curr_vert.NumberTrajectoryPoints(); iloc++) {
					if (curr_vert.DQdxAtPoint(iloc,plane) > 1000
						|| curr_vert.DQdxAtPoint(iloc, plane) < 3
						|| curr_vert.Length(iloc) < 20) continue;
					TVector3 point = curr_vert.LocationAtPoint(iloc);
					int binX = h5e19_raw->FindBin(point.Z());
					hAfter_raw->Fill(point.Z(), curr_vert.DQdxAtPoint(iloc, plane));
					binFill_after[binX].push_back(curr_vert.DQdxAtPoint(iloc, plane));
					totalCount_after.push_back(curr_vert.DQdxAtPoint(iloc,plane));
					//if (curr_vert.DQdxAtPoint(iloc,plane) < 1) hLowHit1D->Fill(point.X());

				}
			}
		}
	
		for (size_t idx {0}; idx < binFill_after.size() ; idx++) {
			if (binFill_after[idx].size() > 5) {
				double avg_dQdx = hAfter_raw->GetBinContent(idx) / binFill_after[idx].size();
				hAfter_avg->SetBinContent(idx, avg_dQdx);
			}
		}

		for (size_t ivert{0}; ivert < AllVertex_m_before.size(); ivert++) {
			for (size_t itrack{0}; itrack < AllVertex_m_before[ivert].size(); itrack++) {
				auto curr_vert = AllVertex_m_before[ivert][itrack];
				for (size_t iloc {0}; iloc < curr_vert.NumberTrajectoryPoints(); iloc++) {
					if (curr_vert.DQdxAtPoint(iloc,plane) > 1000
						|| curr_vert.DQdxAtPoint(iloc, plane) < 3
						|| curr_vert.Length(iloc) < 20) continue;
					TVector3 point = curr_vert.LocationAtPoint(iloc);
					int binX = h5e19_raw->FindBin(point.X());
					hBefore_raw->Fill(point.X(), curr_vert.DQdxAtPoint(iloc, plane));
					binFill_before[binX].push_back(curr_vert.DQdxAtPoint(iloc, plane));
					totalCount_before.push_back(curr_vert.DQdxAtPoint(iloc,plane));
					//if (curr_vert.DQdxAtPoint(iloc,plane) < 1) hLowHit1D->Fill(point.X());

				}
			}
		}
	
		for (size_t idx {0}; idx < binFill_before.size() ; idx++) {
			if (binFill_before[idx].size() > 5) {
				double avg_dQdx = hBefore_raw->GetBinContent(idx) / binFill_before[idx].size();
				hBefore_avg->SetBinContent(idx, avg_dQdx);
			}
		}


	}

	void Cal_Test::MakeXCorrectionMap() {
		hHitMap1D = new TH1D("HitMap1D","HitMap1D", Nx, xmin, xmax);
		hRawDQdxMap1D = new TH1D("RawDQdxMap1D","RawDQdxMap1D", Nx, xmin, xmax);
		hAvgDQdxMap1D = new TH1D("AvgDQdxMap1D","AvgDQdxMap1D", Nx, xmin, xmax);
		hMedianDQdxMap1D = new TH1D("MedianDQdxMap1D","MedianDQdxMap1D", Nx, xmin, xmax);
		hCalibrationMap1D = new TH1D("CalibrationMap1D","CalibrationMap1D", Nx, xmin, xmax);
		hLowHit1D = new TH1D("LowHit1D","LowHit1D", Nx, xmin, xmax);



		std::vector<std::vector<double>> binFill {};
		std::vector<double> totalCount {};
		binFill.resize(hHitMap1D->GetSize());


		for (size_t ivert{0}; ivert < AllVertex_m.size(); ivert++) {
			for (size_t itrack{0}; itrack < AllVertex_m[ivert].size(); itrack++) {
				auto curr_vert = AllVertex_m[ivert][itrack];
				for (size_t iloc {0}; iloc < curr_vert.NumberTrajectoryPoints(); iloc++) {
					if (curr_vert.DQdxAtPoint(iloc,plane) > 20
						|| curr_vert.DQdxAtPoint(iloc, plane) < 5
						|| curr_vert.Length(iloc) < 20) continue;
					TVector3 point = curr_vert.LocationAtPoint(iloc);
					int bin = hHitMap2D->FindBin(point.Z(), point.Y());
					int binX = hHitMap1D->FindBin(point.X());
					int curr = hHitMap1D->GetBinContent(bin);
					double calFac = hCalibrationMap2D->GetBinContent(bin);
					hHitMap1D->Fill(point.X());
					hRawDQdxMap1D->Fill(point.X(), curr_vert.DQdxAtPoint(iloc, plane));
					binFill[binX].push_back(curr_vert.DQdxAtPoint(iloc, plane));
					totalCount.push_back(curr_vert.DQdxAtPoint(iloc,plane));
					//if (curr_vert.DQdxAtPoint(iloc,plane) < 1) hLowHit1D->Fill(point.X());




				}
			}
		}
	
		for (size_t idx {0}; idx < binFill.size() ; idx++) {
			if (binFill[idx].size() > 5) {
				double med_dQdx {0.0};
				double avg_dQdx = hRawDQdxMap1D->GetBinContent(idx) / binFill[idx].size();
				hAvgDQdxMap1D->SetBinContent(idx, avg_dQdx);
				auto temp = binFill[idx];
				std::sort(temp.begin(),temp.end());
				if (temp.size() % 2 == 1) {
					med_dQdx = temp[temp.size() / 2];
				}
				else {
					med_dQdx = 0.5*(temp[temp.size() / 2 - 1] + temp[temp.size()/2]);
				}
				hMedianDQdxMap1D->SetBinContent(idx, med_dQdx);
			}
			// else {
			// 	hAvgDQdxMap1D->SetBinContent(idx, 0);
			// 	hMedianDQdxMap1D->SetBinContent(idx, 0);

			// }
		}
		double avgVal_X = std::accumulate(totalCount.begin(), totalCount.end(), 0.0) / totalCount.size();
		// std::sort(totalCount.begin(), totalCount.end());
		// double medVal {0.0};
		// if (totalCount.size() % 2 == 1) medVal = totalCount[totalCount.size() / 2];
		// else medVal = 0.5*(totalCount[totalCount.size() / 2 - 1] + totalCount[totalCount.size()/2]); 


		for (size_t idx {0}; idx < hAvgDQdxMap1D->GetSize(); idx++) {
			if (hAvgDQdxMap1D->GetBinContent(idx) != 0) {
				hCalibrationMap1D->SetBinContent(idx, medVal / hAvgDQdxMap1D->GetBinContent(idx));
			}
			else hCalibrationMap1D->SetBinContent(idx, 1);	
		}
	}

	void Cal_Test::Make3DCorretionMap() {

		hCalibrationMap3D = new TH3D("CalibrationMap3D","CalibrationMap3D",Nz,zmin,zmax,Nx,xmin,xmax,Ny,ymin,ymax);
		for (size_t idx {0}; idx < hCalibrationMap3D->GetSize(); idx++) {
			hCalibrationMap3D->SetBinContent(idx, 0.0);
		}

		for ( size_t xbin {0}; xbin < Nx; xbin++) {
			for ( size_t ybin {0}; ybin < Ny; ybin++) {
				for ( size_t zbin {0}; zbin < Nz; zbin++) {
					int bin = hCalibrationMap3D->GetBin(zbin,xbin,ybin);
					double calYZ = hCalibrationMap2D->GetBinContent(zbin,ybin);
					double calX = hCalibrationMap1D->GetBinContent(xbin);
					if (hCalibrationMap3D->GetBinContent(bin) == 0 && calYZ*calX < 1.8) {
						hCalibrationMap3D->SetBinContent(bin, calYZ*calX);
					}
					else {
						hCalibrationMap3D->SetBinContent(bin, 1.0);
					}

				}
			}
		}
	}

	void Cal_Test::MakeCorrecteddQdxPlots() {

		hCorrecteddQdx = new TH2D("CorrecteddQdx","Corrected, Calibrated (All Tracks);Residual Length (cm);Calibrated dQdx",400,0,400,300,0,300);
		hCorrecteddQdx_p = new TH2D("CorrecteddQdx_p","Corrected, Calibrated (Protons);Residual Length (cm);Calibrated dQdx",400,0,400,300,0,300);
		hCorrecteddQdx_m = new TH2D("CorrecteddQdx_m","Corrected, Calibrated (Muons);Residual Length (cm);Calibrated dQdx",400,0,400,300,0,300);

		double avg_dQdx[2] {0,0};

		for (size_t ivert = 0; ivert < AllVertex_p.size(); ivert++) {
			for (size_t itrack = 0; itrack < AllVertex_p[ivert].size(); itrack++) {
				auto curr_vert = AllVertex_p[ivert][itrack];
				for (size_t iloc = 0; iloc < curr_vert.NumberTrajectoryPoints(); iloc++) {
					double raw_dQdx = curr_vert.DQdxAtPoint(iloc,plane);
					TVector3 point = curr_vert.LocationAtPoint(iloc);
					double reslength = curr_vert.Length(iloc);
					if (reslength < 5) continue;
					if (raw_dQdx <= 0 || raw_dQdx == 111) continue;
					int bin = hCalibrationMap3D->FindBin(point.Z(), point.X(), point.Y());
					double calFac = hCalibrationMap3D->GetBinContent(bin);
					double corr_dQdx = raw_dQdx*calFac;
					avg_dQdx[itrack] += corr_dQdx;
					hCorrecteddQdx->Fill(reslength, corr_dQdx);
				}
				avg_dQdx[itrack] /= AllVertex_p[ivert][itrack].NumberTrajectoryPoints();
			}
			int kProton = 0;
			int kMuon = 1;
			if (AllVertex_p[ivert][kProton].Length() < 5 || AllVertex_p[ivert][kMuon].Length() < 5) continue;
			if (AllVertex_p[ivert][kProton].Length() > AllVertex_p[ivert][kMuon].Length()) continue;
			if (avg_dQdx[kProton] < avg_dQdx[kMuon]) continue;
			for (size_t iloc = 0; iloc < AllVertex_p[ivert][kProton].NumberTrajectoryPoints(); iloc++){
				auto curr_vert = AllVertex_p[ivert][kProton];
				double raw_dQdx = curr_vert.DQdxAtPoint(iloc,plane);
				TVector3 point = curr_vert.LocationAtPoint(iloc);
				double reslength = curr_vert.Length(iloc);
				if (reslength < 2) continue;
				if (raw_dQdx <= 0 || raw_dQdx == 111) continue;
				int bin = hCalibrationMap3D->FindBin(point.Z(), point.X(), point.Y());
				double calFac = hCalibrationMap3D->GetBinContent(bin);
				double corr_dQdx = raw_dQdx*calFac;
				hCorrecteddQdx_p->Fill(reslength, corr_dQdx);
			}
			for (size_t iloc = 0; iloc < AllVertex_p[ivert][kMuon].NumberTrajectoryPoints(); iloc++){
				auto curr_vert = AllVertex_p[ivert][kMuon];
				double raw_dQdx = curr_vert.DQdxAtPoint(iloc,plane);
				TVector3 point = curr_vert.LocationAtPoint(iloc);
				double reslength = curr_vert.Length(iloc);
				if (reslength < 2) continue;
				if (raw_dQdx <= 0 || raw_dQdx == 111) continue;
				int bin = hCalibrationMap3D->FindBin(point.Z(), point.X(), point.Y());
				double calFac = hCalibrationMap3D->GetBinContent(bin);
				double corr_dQdx = raw_dQdx*calFac;
				hCorrecteddQdx_m->Fill(reslength, corr_dQdx);
			}
		}
		std::cout << hCorrecteddQdx_p->GetEntries() << std::endl;
		std::cout << hCorrecteddQdx_m->GetEntries() << std::endl;
	}

	void Cal_Test::MakeVoxelCorrectionMap() {
		hHitMap3D = new TH3D("HitMap3D","HitMap3D", Nz, zmin, zmax, Nx, xmin, xmax, Ny, ymin, ymax);
		hAvgDQdxMap3D = new TH3D("AverageDQdx3D","AverageDQdx3D", Nz, zmin, zmax, Nx, xmin, xmax, Ny, ymin, ymax);
		hCorrectionMap3D = new TH3D("CorrectionMap3D","CorrectionMap3D", Nz, zmin, zmax, Nx, xmin, xmax, Ny, ymin, ymax);
		hRawDQdxMap3D = new TH3D("RawDQdx3D","RawDQdx3D", Nz, zmin, zmax, Nx, xmin, xmax, Ny, ymin, ymax);
		hLowHit3D = new TH3D("LowHit3D","LowHit3D", Nz, zmin, zmax, Nx, xmin, xmax, Ny, ymin, ymax);


		std::vector<std::vector<double>> binFill;
		std::vector<double> totalCount;
		binFill.resize(hHitMap3D->GetSize());

		for (size_t ivert{0}; ivert < AllVertex_m.size(); ivert++) {
			for (size_t itrack{0}; itrack < AllVertex_m[ivert].size(); itrack++) {
				auto curr_vert = AllVertex_m[ivert][itrack];
				for (size_t iloc {0}; iloc < curr_vert.NumberTrajectoryPoints(); iloc++) {
					if (curr_vert.DQdxAtPoint(iloc,plane) > 100 // maxADC
						|| curr_vert.DQdxAtPoint(iloc, plane) < 5
						|| curr_vert.Length(iloc) < 5) continue;
					TVector3 point = curr_vert.LocationAtPoint(iloc);
					auto bin = hHitMap3D->FindBin(point.Z(), point.X(), point.Y());
					int curr = hHitMap3D->GetBinContent(bin);
					hHitMap3D->Fill(point.Z(), point.X(), point.Y());
					hRawDQdxMap3D->Fill(point.Z(), point.X(), point.Y(), curr_vert.DQdxAtPoint(iloc, plane));
					binFill[bin].push_back(curr_vert.DQdxAtPoint(iloc,plane));
					totalCount.push_back(curr_vert.DQdxAtPoint(iloc,plane));
					//hTotalCount->Fill(curr_vert.DQdxAtPoint(iloc,plane));
					//if (curr_vert.DQdxAtPoint(iloc,plane) < 1) hLowHit3D->Fill(point.Z(),point.X(),point.Y());
				}
			}
		}

		if (hHitMap3D->GetSize() != hRawDQdxMap3D->GetSize()) {
			std::cout << "Histogram's don't have same number of bins!";
		std::cout << binFill.size() << " " << totalCount.size() << std::endl;
		}
		else {
			for (size_t idx {0}; idx < binFill.size() ; idx++) {
				if (binFill[idx].size() > 1) {
					//std::cout << "pass \n";
					double med_dQdx {0.0};
					double avg_dQdx = std::accumulate(binFill[idx].begin(), binFill[idx].end(), 0.0) / binFill[idx].size();
					hAvgDQdxMap3D->SetBinContent(idx, avg_dQdx);
					auto temp = binFill[idx];

				}
			}
		}

		// medVal = std::accumulate(totalCount.begin(), totalCount.end(), 0.0) / totalCount.size();

		for (size_t idx {0}; idx < hAvgDQdxMap3D->GetSize(); idx++) {
			if (hAvgDQdxMap3D->GetBinContent(idx) != 0.0) {
				double calFac = medVal/ hAvgDQdxMap3D->GetBinContent(idx);
				//std::cout << calFac << std::endl;
				if (calFac < 2.0) 
					hCorrectionMap3D->SetBinContent(idx, calFac);
				else {
					hCorrectionMap3D->SetBinContent(idx, 1);
				}
			}

			//else hCorrectionMap3D->SetBinContent(idx, 1);
		}
	}

	void Cal_Test::CosmicTracks() {

		hCosmicdQdx = new TH1D("hCosmicdQdx","hCosmicdQdx",100,0,50);
		hCosmicdQdx_cor = new TH1D("hCosmicdQdx_cor","hCosmicdQdx_cor",100,0,50);
		hWeirdWires = new TH1D("hWeirdWires","hWeirdWires",50,36.9999999,37.0000001);
		hWeird = new TH3D("Weird","??",Nz, zmin, zmax, Nx, xmin, xmax, Ny, ymin, ymax);

		for (size_t ivert = 0; ivert < AllVertex_m.size(); ivert++) {
			for (size_t itrack = 0; itrack < AllVertex_m[ivert].size(); itrack++) {
				auto curr_vert = AllVertex_m[ivert][itrack];
				for (size_t iloc = 0; iloc < curr_vert.NumberTrajectoryPoints(); iloc++) {
					double raw_dQdx = curr_vert.DQdxAtPoint(iloc,plane);
					TVector3 point = curr_vert.LocationAtPoint(iloc);
					double reslength = curr_vert.Length(iloc);
					int bin = hCorrectionMap3D->FindBin(point.Z(), point.X(), point.Y());
					double calFact = hCalibrationMap3D->GetBinContent(bin);
					if (calFact == 1) continue;
					//if (reslength < 5) continue;
					if (raw_dQdx <= 0 || raw_dQdx == 111) continue;
					hCosmicdQdx->Fill(raw_dQdx);
					hCosmicdQdx_cor->Fill(raw_dQdx*calFact);
					hWeirdWires->Fill(raw_dQdx);
					if (raw_dQdx > 36.9 && raw_dQdx < 37.1) {
						hWeird->Fill(point.Z(),point.X(),point.Y());
						//std::cout << raw_dQdx << std::endl;
					}

				}
			}
			
		}
		TF1 *fgaus = new TF1("fgaus","gaus(0)",5,20);
		hCosmicdQdx_cor->Fit("fgaus","","", 5, 20);
		TCanvas *WeirdWires = new TCanvas("WW","WW",200, 10, 750, 300);
		WeirdWires->cd();
		hWeirdWires->Draw();
	}

	void RecombinationFit() {

	}

	void Cal_Test::RunInfo() { 
		hADC_Total = new TH1D("hADC_Total","hADC_Total",100,4500,8500);
		hADC_Hits = new TH1D("hADC_Hits","hADC_Hits",100,4500,8500);
		hADC_Average = new TH1D("hADC_Average","hADC_Average",100,4500,8500);
	



		for (size_t ivert{0}; ivert < AllVertex_m.size(); ivert++) {
			int binX = hADC_Total->FindBin(AllVertex_run[ivert][0]);
			//std::cout << binX << std::endl;
			for (size_t itrack{0}; itrack < AllVertex_m[ivert].size(); itrack++) {
				auto curr_vert = AllVertex_m[ivert][itrack];
				for (size_t iloc {0}; iloc < curr_vert.NumberTrajectoryPoints(); iloc++) {
					if (curr_vert.DQdxAtPoint(iloc,plane) > 1000
						|| curr_vert.DQdxAtPoint(iloc, plane) < 0
						|| curr_vert.Length(iloc) < 30) continue;
					hADC_Total->Fill(AllVertex_run[ivert][0], curr_vert.DQdxAtPoint(iloc, plane));
					hADC_Hits->Fill(AllVertex_run[ivert][0]);

				}
			}
		}
		int size = hADC_Hits->GetSize();
		for (size_t idx {0}; idx <= size; idx++) {
			double Avg_ADC = hADC_Total->GetBinContent(idx) / hADC_Hits->GetBinContent(idx);
			std::cout << Avg_ADC << std::endl;
			if (hADC_Hits->GetBinContent(idx) == 0) continue;
			hADC_Average->SetBinContent(idx,Avg_ADC);
		}
		hADC_Average->Draw();
	}
}


#endif