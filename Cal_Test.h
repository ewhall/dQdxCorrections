#ifndef CAL_TEST_H
#define CAL_TEST_H

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
#include "TFile.h"
#include "TTree.h"
#include "TSpline.h"
#include <sys/stat.h>
#include <unistd.h>
#include <vector>
#include <string>

namespace larlite {

	class Cal_Test : public ana_base {

	public:

		Cal_Test(){ _name="Cal_Test"; _fout=0;}
		virtual ~Cal_Test(){}
		virtual bool Initialize();
		virtual bool Finalize();
		void LoadSelectionFiles(std::string file_loc);
		void AllVertexFill(storage_manager * mgr);
		void SelectedVertex(int run, int subrun, int event, int vtx_id);
		void MakeDQdxPlots();
		void MakeCorrecteddQdxPlots();
		void Make2DCorrectionMap();
		void MakeXCorrectionMap();
		void Make3DCorretionMap();
		void MakeVoxelCorrectionMap();
		void CosmicTracks();
		void LoadDEdxSplines(std::string file_loc);
		void LoadCalibrationFile(std::string file_loc);
		void LookAtX();
		void LookAtZ();
		void RunInfo();
		void ImageCorrectionMap();
		void DataGainFactor();

	protected:
		int run, subrun, event, vtx_id;
		int Nz, Ny, Nx;
		double xmin, xmax, ymin, ymax, zmin, zmax;
		double avgVal, medVal;
		double minADC, maxADC;
		bool selectedProton, selectedCosmic, selectedCosmic_during, selectedCosmic_after, selectedCosmic_before;
		double Wion, E_Field, LAr_Density, alpha_ex, beta_ex;
		double avgCor;
		double alpha_fit, beta_fit;

		larlite::geo::View_t plane; // change this to expand to 3 planes


		std::vector<std::vector<int>> SelectEvtID_proton;
		std::vector<std::vector<int>> SelectEvtID_cosmic;
		std::vector<larlite::track> thisVertex;
		std::vector<std::vector<larlite::track>> AllVertex_m;
		std::vector<std::vector<int>> AllVertex_run;
		std::vector<std::vector<larlite::track>> AllVertex_m_during;
		std::vector<std::vector<larlite::track>> AllVertex_m_after;
		std::vector<std::vector<larlite::track>> AllVertex_m_before;
        std::vector<std::vector<larlite::track>> AllVertex_p;
        std::vector<std::vector<double>> binFill;

        TH2D *hUncorrecteddQdx;
        TH2D *hUncorrecteddQdx_p;
        TH2D *hUncorrecteddQdx_m;
        TH2D *hCorrecteddQdx;
        TH2D *hCorrecteddQdx_p;
        TH2D *hCorrecteddQdx_m;
        TH2D *hHitMap2D;
        TH2D *hRawDQdxMap2D;
        TH2D *hAvgDQdxMap2D;
        TH2D *hMedianDQdxMap2D;
        TH2D *hCalibrationMap2D;
		TH1D *hHitMap1D;
		TH1D *hRawDQdxMap1D;
		TH1D *hAvgDQdxMap1D;
		TH1D *hMedianDQdxMap1D;
		TH1D *hCalibrationMap1D;
		TH3D *hCalibrationMap3D;

        TSpline3 *sProtonRange2dEdx;
        TSpline3 *sMuonRange2dEdx;
        TSpline3 *sProtonRange2T;

		TH3D *hHitMap3D;
		TH3D *hLowHits3D;
		TH3D *hRawDQdxMap3D;
		TH3D *hAvgDQdxMap3D;
		TH3D *hCorrectionMap3D;
		TH3D *hWeird;
		TH1D *hADC[4];
		TH1D *hADC_Corr[4];

		TH1D *hStart;
		TH1D *hEnd;

		TH2D *hLowHit2D;
		TH3D *hLowHit3D;
		TH1D *hLowHit1D;
		TH1D *hTotalCount;
		TH1D *hCosmicdQdx;
		TH1D *hCosmicdQdx_cor;
		TH1D *hWeirdWires;
		TH2D *hRecoFit;

        TF1 *fRecombinationExpected;
        TF1 *f1[4];
        TF1 *f1_corr[4];

        TH1D *h5e19_avg;
        TH1D *h5e19_raw;
        TH1D *hAfter_avg;
        TH1D *hAfter_raw;
        TH1D *hBefore_avg;
        TH1D *hBefore_raw;

        TH1D *hADC_Total;
        TH1D *hADC_Hits;
        TH1D *hADC_Average;

        TH3D *hImageCalibrationMap[3];
        TH2D *h2DImageMap[3];
        TH2D *h2DImageMap_ZY[3];
        TH2D *h2DRawImageMap_ZY[3];
        TH1D *h2DTotalCount[3];
        TH2D *h2DRawImageMap[3];
        TH2D *h2DImageCalMap[3];
        TH2D *h2DImageCalMap_ZY[3];
        TH1D *hXBins;
        TH2D *h2DImageHitMap[3];
        TH2D *h2DImageHitMap_ZY[3];

        TH3D *CorrectionMap3D_Y;
        TH3D *CorrectionMap3D_U;
        TH3D *CorrectionMap3D_V;

        //larlite::geo::View_t views[3];



	};
}


#endif