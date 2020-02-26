#ifndef DEDX_TEST_H
#define DEDX_TEST_H

#include "Analysis/ana_base.h"
#include "DataFormat/track.h"
#include "LArUtil/Geometry.h"
#include "LArUtil/GeometryHelper.h"
#include "TGraph2D.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TImage.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TSpline.h"
#include <sys/stat.h>
#include <unistd.h>
#include <vector>
#include <string>
#include "RooDouble.h"

namespace larlite {

	class DEdx_Test : public ana_base {

	public: 
		DEdx_Test(){ _name="DEdx_Test"; _fout=0;}
		virtual ~DEdx_Test(){}
		virtual bool Initialize();
		virtual bool Finalize();
		void LoadCalibrationFile(std::string file_loc);
		void LoadSelectionFiles(std::string file_loc);
		void LoadDEdxSplines(std::string file_loc);
		void AllVertexFill(storage_manager * mgr);
		void SelectedVertex(int run, int subrun, int event, int vtx_id);
		void CorrectedTracks();
		void UncorrectedTracks();
		void RecombinationFit();
		void MakeDEDXPlots();
		void LoadRecombinationFile(std::string file_loc);
		void LoadTruthFile(std::string file_loc);
		void MakeRecoFill();
		void CreateDEDXMap(std::string file_loc);
		void DEdxPlots();
                void Tracks_All();
                void RecoFit_new();
                void All_DEdxPlots();





	protected:
		int run, subrun, event, vtx_id;
		int Nz, Ny, Nx;
		double xmin, xmax, ymin, ymax, zmin, zmax;
		bool selectedProton, selectedCosmic;
		double Wion, E_Field, LAr_Density, alpha_ex, beta_ex;
		double avgVal, avgVal_cor;
                double avgVal_U, avgVal_V, avgVal_Y, avgVal_all;
                double avgVal_U_Corr, avgVal_V_Corr, avgVal_Y_Corr, avgVal_all_Corr;
		double alpha_fit, beta_fit, alpha_sigma, beta_sigma, alpha_fit_un, beta_fit_un, first_run, last_run;
                double alpha_fit_Y, beta_fit_Y, alpha_fit_All, beta_fit_All;

		larlite::geo::View_t plane;	


		std::vector<std::vector<int>> SelectEvtID_proton;
		std::vector<std::vector<std::vector<double>>> SelectTruthID_energy;
		std::vector<std::vector<double>> SelectTruthID_proton;
		std::vector<std::vector<int>> SelectEvtID_cosmic;
		std::vector<std::vector<double>> AllTruth_p;
		std::vector<larlite::track> thisVertex;
		std::vector<std::vector<larlite::track>> AllVertex_m;
        std::vector<std::vector<larlite::track>> AllVertex_p;
        std::vector<std::vector<std::vector<double>>> dEdxInfo_p;
        std::vector<std::vector<std::vector<double>>> TruthInfo_p;
        std::vector<std::vector<double>> weird_wires;
        std::vector<double> mu;
        std::vector<double> sigma;
        std::vector<double> mu_error;
        std::vector<double> sigma_error;
        std::vector<double> ebin;
        std::vector<double> ey;
        std::vector<double> rbin;
        std::vector<double> ry;
        std::vector<double> mu_Y;
        std::vector<double> mu_error_Y;
        std::vector<double> mu_All;
        std::vector<double> mu_error_All;
        std::vector<double> mu_Y_m;
        std::vector<double> mu_error_Y_m;
        std::vector<double> mu_All_m;
        std::vector<double> mu_error_All_m;

        TSpline3 *sProtonRange2dEdx;
        TSpline3 *sMuonRange2dEdx;
        TSpline3 *sProtonRange2T;

        TH3D *CorrectionMap3D;
        TH3D *CorrectionMap3D_U;
        TH3D *CorrectionMap3D_V;
        TH3D *hWeird;

        TH2D *hRecombinationPlot;
        TH2D *hUncorrectedRecombinationPlot;
        TH2D *hUncorrectedDQdx;
        TH2D *hUncorrectedDQdx_p;
        TH2D *hUncorrectedDQdx_m;
        TH2D *RecoFit;
        TH2D *hRecoFit;

        TH2D *hCorrectedRecombinationPlot;
        TH2D *hCorrectedDQdx;
        TH2D *hCorrectedDQdx_p;
        TH2D *hCorrectedDQdx_m;

        TH2D *hCorrectedDEdx;
        TH2D *hCorrectedDEdx_p;
        TH2D *hCorrectedDEdx_m;
        TH2D *hEnergy_p;

        TF1 *fRecombination;
        TF1 *fRecombination2;
        TF1 *fRecombination3;
        TF1 *fRecombination4;
        TF1 *fRecombinationExpected;
        TF1 *fRecombination_un;
        TF1 *fdEdxConvert;
        TF1 *fdEdxConvert_Y;
        TF1 *fdEdxConvert_All;

        TH1D *hDQdx_p_cor;
        TH1D *hDQdx_m_cor;
        TH1D *hDQdx_p_un;
        TH1D *hDQdx_m_un;
        TH1D *hDEdx_p;
        TH1D *hDEdx_m;
        TH1D *hdr;
        TH1D *hEnergy;
        TH1D *hUncorrectedADC;
        TH1D *hADC_Corr_All;

        TH1D *hAvg_dEdx_p;
        TH1D *hAvg_dEdx_p_5cm;
        TH1D *hAvg_dEdx_p_10cm;
        TH1D *hAvg_dEdx_m;
        TH1D *hAvg_dEdx_m_5cm;
        TH1D *hAvg_dEdx_m_10cm;


        TH1D *hAvg_corrdQdx_p;
        TH1D *hAvg_corrdQdx_p_5cm;
        TH1D *hAvg_corrdQdx_p_10cm;
        TH1D *hAvg_corrdQdx_m;
        TH1D *hAvg_corrdQdx_m_5cm;
        TH1D *hAvg_corrdQdx_m_10cm;
        TH1D *hChi2_m_m;
        TH1D *hChi2_m_p;
        TH1D *hChi2_p_m;
        TH1D *hChi2_p_p;
        TH1D *hChi2_Inv_m_m;
        TH1D *hChi2_Inv_m_p;
        TH1D *hChi2_Inv_p_m;
        TH1D *hChi2_Inv_p_p;


        TH1D *hAvg_dEdx_p_un;
        TH1D *hAvg_dEdx_p_5cm_un;
        TH1D *hAvg_dEdx_p_10cm_un;
        TH1D *hAvg_dEdx_m_un;
        TH1D *hAvg_dEdx_m_5cm_un;
        TH1D *hAvg_dEdx_m_10cm_un;
        TH1D *hAvg_corrdQdx_p_un;
        TH1D *hAvg_corrdQdx_p_5cm_un;
        TH1D *hAvg_corrdQdx_p_10cm_un;
        TH1D *hAvg_corrdQdx_m_un;
        TH1D *hAvg_corrdQdx_m_5cm_un;
        TH1D *hAvg_corrdQdx_m_10cm_un;
        TH1D *hChi2_m_m_un;
        TH1D *hChi2_m_p_un;
        TH1D *hChi2_p_m_un;
        TH1D *hChi2_p_p_un;
        TH1D *hChi2_Inv_m_m_un;
        TH1D *hChi2_Inv_m_p_un;
        TH1D *hChi2_Inv_p_m_un;
        TH1D *hChi2_Inv_p_p_un;
        TH2D *hChi2_2D_p;
        TH2D *hChi2_2D_m;

        TH2D *hUnDQdx_All;
        TH2D *hUnDQdx_All_p;
        TH2D *hUnDQdx_All_m;
        TH1D *hADC_All;

        TH2D* hCorrDQdx_All;
        TH2D* hCorrDQdx_All_p;
        TH2D* hCorrDQdx_All_m;

        TH2D *hUncorrectedRecombinationPlot_All;
        TH2D *hCorrectedRecombinationPlot_All;
        TH2D *hCorrectedRecombinationPlot_Y;
        TH2D *hUncorrectedRecombinationPlot_Y;

        TH2D *hCorrectedDEdx_Y;
        TH2D *hCorrectedDEdx_Y_p;
        TH2D *hCorrectedDEdx_Y_m;
        TH2D *hCorrectedDEdx_All;
        TH2D *hCorrectedDEdx_All_p;
        TH2D *hCorrectedDEdx_All_m;

        std::vector<double> mu_Un_Y, mu_error_Un_Y, sigma_Un_Y, sigma_error_Un_Y;
        std::vector<double> mu_Corr_Y, mu_error_Corr_Y, sigma_Corr_Y, sigma_error_Corr_Y;
        std::vector<double> mu_Un_All, mu_error_Un_All, sigma_Un_All, sigma_error_Un_All;
        std::vector<double> mu_Corr_All, mu_error_Corr_All, sigma_Corr_All, sigma_error_Corr_All;





        TH1D *h5e19_2;
        TH1D *h5e19_2_raw;


        TTree *RecombinationParameters;

        std::map<std::string, std::vector<std::vector<double>>> dEdx_Map;
        std::map<std::string, std::vector<std::vector<double>>> Truth_Map;
        bool dEdx_Info;
        std::vector<std::vector<double>> dEdxInfo;
        std::vector<std::vector<double>> TruthInfo;



	};



}






#endif