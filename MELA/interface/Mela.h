#ifndef MELA_Mela_h
#define MELA_Mela_h

#include "TLorentzVector.h"
#include <vector>
#include <TRandom3.h>


class TFile;
class TH1F;
class TH2F;
class TH3F;
class RooRealVar;
class RooAbsPdf;
class RooArgSet;
class AngularPdfFactory;
class VectorPdfFactory;
class TensorPdfFactory;
class RooqqZZ_JHU_ZgammaZZ_fast;
class newZZMatrixElement;
class TGraph;
class SuperMELA;

#include <ZZMatrixElement/MELA/interface/TVar.hh>
#include <ZZMatrixElement/MELA/interface/TEvtProb.hh>
#include <ZZMatrixElement/MELA/src/AngularPdfFactory.h>
#include <ZZMatrixElement/MELA/src/VectorPdfFactory.h>
#include <ZZMatrixElement/MELA/src/TensorPdfFactory.h>
#include <ZZMatrixElement/MELA/src/RooqqZZ_JHU_ZgammaZZ_fast.h>


class Mela{

public:
  
  // Mela(){};
  Mela(int LHCsqrts=8, float mh=125); // higgs mass for supermela
  ~Mela();
  
  void setProcess(TVar::Process myModel, TVar::MatrixElement myME, TVar::Production myProduction);
  void setMelaHiggsWidth(float myHiggsWidth=-1);
  void setMelaLeptonInterference(TVar::LeptonInterference myLepInterf = TVar::DefaultLeptonInterf);
  void setRemoveLeptonMasses(bool MasslessLeptonSwitch = false);
  void resetMCFM_EWKParameters(double ext_Gf, double ext_aemmz, double ext_mW, double ext_mZ, double ext_xW, int ext_ewscheme=3);
  void setRenFacScaleMode(TVar::EventScaleScheme renormalizationSch, TVar::EventScaleScheme factorizationSch, double ren_sf, double fac_sf);

  void computeP(float mZZ, float mZ1, float mZ2, // input kinematics
		float costhetastar,
		float costheta1,
		float costheta2,
		float phi,
		float phi1,
		int flavor,
		float& prob,            // output probability
		bool useConstant=true
		);

  void computeP(float mZZ, float mZ1, float mZ2, // input kinematics
		float costhetastar,
		float costheta1, 
		float costheta2,
		float phi,
		float phi1,
		int flavor,
    double selfDHvvcoupl[SIZE_HVV][2],
		float& prob                   // output probability
		);

  void computeD_CP(float mZZ, float mZ1, float mZ2, // input kinematics
		float costhetastar,
		float costheta1,
		float costheta2,
		float phi,
		float phi1,
		int flavor,
		TVar::MatrixElement myME,
		TVar::Process myType,
		float& prob
		);

  void computeP_selfDspin2(float mZZ, float mZ1, float mZ2, // spin 2 arbitray coupling 
        float costhetastar,
        float costheta1,
        float costheta2,
        float phi,
        float phi1,
        int flavor,
        double selfDGggcoupl[SIZE_GGG][2],
        double selfDGvvcoupl[SIZE_GVV][2],
        float& prob
		);

  void computeP_selfDspin1(float mZZ, float mZ1, float mZ2, // arbitrary spin 1 coupling 
        float costhetastar,
        float costheta1,
        float costheta2,
        float phi,
        float phi1,
        int flavor,
        double selfDZvvcoupl[SIZE_ZVV][2],
        float& prob
		);

  void computeP(float mZZ, float mZ1, float mZ2, // input kinematics
		float costhetastar,
		float costheta1, 
		float costheta2,
		float phi,
		float phi1,
		int flavor,
		double couplingvals[SIZE_HVV_FREENORM],
		float& prob                   // output probability
		);
    
  void computeP(TLorentzVector Z1_lept1, int Z1_lept1Id,  // input 4-vectors
		TLorentzVector Z1_lept2, int Z1_lept2Id,  // 
		TLorentzVector Z2_lept1, int Z2_lept1Id,
		TLorentzVector Z2_lept2, int Z2_lept2Id,
		float& prob                             // output probability
		);

  void computeP(TLorentzVector Z1_lept1, int Z1_lept1Id,  // input 4-vectors
		TLorentzVector Z1_lept2, int Z1_lept2Id,  // 
		TLorentzVector Z2_lept1, int Z2_lept1Id,
		TLorentzVector Z2_lept2, int Z2_lept2Id,
		double couplingvals[SIZE_HVV_FREENORM],
		float& prob                             // output probability
		);

  void computeProdDecP(
    TLorentzVector jet[2],
    TLorentzVector Higgs_daughter[4],
    int jet_pdgid[2],
    int Higgs_daughter_pdgid[4],

    double selfDHvvcoupl[SIZE_HVV][2],
    double selfDHwwcoupl[SIZE_HVV][2],
    float& prob
    );

  void computeProdP(TLorentzVector Jet1, int Jet1_Id,
		TLorentzVector Jet2, int Jet2_Id,
		TLorentzVector Decay1, int Decay1_Id,
		TLorentzVector Decay2, int Decay2_Id,
		double selfDHggcoupl[SIZE_HGG][2],
		double selfDHvvcoupl[SIZE_HVV_VBF][2],
		double selfDHwwcoupl[SIZE_HWW_VBF][2],
		float& prob
		);

  void computeProdP(TLorentzVector Jet1, int Jet1_Id,
		TLorentzVector Jet2, int Jet2_Id,
		TLorentzVector Decay1, int Decay1_Id,
		TLorentzVector Decay2, int Decay2_Id,
		float& prob
		);

  void computeProdP(
		TLorentzVector V_daughter[2],
		TLorentzVector Higgs_daughter[4],
		int V_daughter_pdgid[2],
		int Higgs_daughter_pdgid[4],
		bool includeHiggsDecay,
		double selfDHvvcoupl[SIZE_HVV_VBF][2],
		float& prob
		);

  void computeProdP(
    TLorentzVector vTTH[6],
    TLorentzVector Higgs,
    int ttbar_daughters_pdgid[6],
    double selfDHvvcoupl[SIZE_TTH][2],
    float& prob,
    int topDecay=0,
    int topProcess=2
    );

  void computeProdP(
    TLorentzVector vTTH[6],
    TLorentzVector Higgs,
    int ttbar_daughters_pdgid[6],
    float& prob,
    int topDecay=0,
    int topProcess=2
    );

  void computeProdP(
    TLorentzVector p_first, int id_first,
    TLorentzVector p_second, int id_second,
    TLorentzVector Higgs,
    float& prob
    );

  void get_PAux(float& prob){ prob = auxiliaryProb; }; // SuperProb
  MelaOutputRecord getFullMERecord(); // Full parton-by-parton ME record (11x11 matrix of MEs + 2 parton luminosity arrays of size 11)

  void computePM4l(float mZZ,
		TVar::LeptonFlavor flavor,
		TVar::SuperMelaSyst syst, 
		float& prob
		); //SuperMela

  void computePM4l(TLorentzVector Z1_lept1, int Z1_lept1Id,  // input 4-vectors
		TLorentzVector Z1_lept2, int Z1_lept2Id,  // 
		TLorentzVector Z2_lept1, int Z2_lept1Id,
		TLorentzVector Z2_lept2, int Z2_lept2Id,
		TVar::SuperMelaSyst syst, 
		float& prob
		); //SuperMela
  
  
  // Ordering of Z1/Z2 according to internal convention
  void checkZorder(float& z1mass, float& z2mass, float& costhetastar, float& costheta1, float& costheta2, float& phi, float& phistar1);

  // Access ZZMEs Calculate4Momentum
  std::vector<TLorentzVector> calculate4Momentum(double Mx,double M1,double M2,double theta,double theta1,double theta2,double Phi1,double Phi);

  // Calculation weight to correct for lepton interference
  void computeWeight(float mZZ, float mZ1, float mZ2, 
		float costhetastar,
		float costheta1, 
		float costheta2,
		float phi,
		float phi1,
		// return variables:
		float& w
		);

  void computeWeight(float mZZ, float mZ1, float mZ2, 
		float costhetastar,
		float costheta1, 
		float costheta2,
		float phi,
		float phi1,
		double couplingvals[SIZE_HVV_FREENORM],
		// return variables:
		float& w
		);

  RooAbsPdf* pdf;
  AngularPdfFactory* spin0Model;
  VectorPdfFactory* spin1Model;
  TensorPdfFactory* spin2Model;
  RooqqZZ_JHU_ZgammaZZ_fast* qqZZmodel;
  SuperMELA* super;
  TRandom3 *myR; // random number for resolution systematics


  RooRealVar* mzz_rrv;
  RooRealVar* z1mass_rrv;
  RooRealVar* z2mass_rrv;
  RooRealVar* costhetastar_rrv; 
  RooRealVar* costheta1_rrv;
  RooRealVar* costheta2_rrv;
  RooRealVar* phi_rrv;
  RooRealVar* phi1_rrv;
  RooRealVar* pt_rrv;  
  RooRealVar* upFrac_rrv;

  TGraph* vaScale_4e;
  TGraph* vaScale_4mu;
  TGraph* vaScale_2e2mu;
  TGraph* DggZZ_scalefactor;
  void setCTotalBkgGraphs(TFile* fcontainer, TGraph* tgC[]);
  void constructDggr(float mzz, int flavor, float bkg_VAMCFM_noscale, float ggzz_VAMCFM_noscale, float ggHZZ_prob_pure_noscale, float ggHZZ_prob_int_noscale, float& myDggr);
  float getConstant(int flavor, float mZZ, bool useOldggZZConstants=false);

  TGraph* tgtotalbkg[3];
  void computeD_gg(float mZZ, float mZ1, float mZ2, // input kinematics
        float costhetastar,
        float costheta1,
        float costheta2,
        float phi,
        float phi1,
        int flavor,
        TVar::MatrixElement myME,
        TVar::Process myType ,
        float& prob
		);

private:

  // 
  // data memmbers 
  // 
  bool usePowhegTemplate_;
  int LHCsqrts;
  TVar::Process myModel_;
  TVar::MatrixElement myME_;
  TVar::Production myProduction_;
  newZZMatrixElement* ZZME;
  
  float auxiliaryProb;
  
  // 
  // functions 
  // 

  void reset_PAux(){ auxiliaryProb=1.; }; // SuperProb reset
};

#endif

