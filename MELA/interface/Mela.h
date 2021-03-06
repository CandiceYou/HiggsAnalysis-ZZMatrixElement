/*
************* HEADER: CMS MELA interface to MCFM/JHUGen-MELA *************
Please see the ../src/Mela.cc file for the instructions.
*/

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
class ScalarPdfFactory_ggH;
class VectorPdfFactory;
class TensorPdfFactory;
class RooqqZZ_JHU_ZgammaZZ_fast;
class newZZMatrixElement;
class TGraph;
class SuperMELA;

#include <ZZMatrixElement/MELA/interface/TVar.hh>
#include <ZZMatrixElement/MELA/interface/TEvtProb.hh>
#include <ZZMatrixElement/MELA/interface/ScalarPdfFactory_ggH.h>
#include <ZZMatrixElement/MELA/interface/VectorPdfFactory.h>
#include <ZZMatrixElement/MELA/interface/TensorPdfFactory_HVV.h>
#include <ZZMatrixElement/MELA/interface/RooqqZZ_JHU_ZgammaZZ_fast.h>

class Mela{

public:

  // Mela(){};
  Mela(double LHCsqrts_=13., double mh_=125., TVar::VerbosityLevel verbosity_=TVar::ERROR); // higgs mass for supermela
  ~Mela();

  void setProcess(TVar::Process myModel, TVar::MatrixElement myME, TVar::Production myProduction);
  void setVerbosity(TVar::VerbosityLevel verbosity_=TVar::ERROR);
  void setMelaLeptonInterference(TVar::LeptonInterference myLepInterf=TVar::DefaultLeptonInterf);
  void setRemoveLeptonMasses(bool MasslessLeptonSwitch=true);
  void setRemoveJetMasses(bool MasslessLeptonSwitch=true);
  void setMelaPrimaryHiggsMass(double myHiggsMass);
  void setMelaHiggsMass(double myHiggsMass, int index=0);
  void setMelaHiggsWidth(double myHiggsWidth=-1, int index=0);
  void setMelaHiggsMassWidth(double myHiggsMass, double myHiggsWidth, int index);
  void setRenFacScaleMode(TVar::EventScaleScheme renormalizationSch, TVar::EventScaleScheme factorizationSch, double ren_sf, double fac_sf);
  void setCandidateDecayMode(TVar::CandidateDecayMode mode);
  void setCurrentCandidateFromIndex(unsigned int icand); // Switches to another candidate
  void setCurrentCandidate(MELACandidate* cand); // Switches to another candidate
  void setInputEvent(
    SimpleParticleCollection_t* pDaughters,
    SimpleParticleCollection_t* pAssociated=0,
    SimpleParticleCollection_t* pMothers=0,
    bool isGen=false
    ); // Adds another candidate
  void resetInputEvent(); // Reset the input candidates. Important to call in order to clean up TEvtProb!
  void setTempCandidate(
    SimpleParticleCollection_t* pDaughters,
    SimpleParticleCollection_t* pAssociated=0,
    SimpleParticleCollection_t* pMothers=0,
    bool isGen=false
    ); // Adds a temp. candidate
  void appendTopCandidate(SimpleParticleCollection_t* TopDaughters); // Adds a top

  // Function to set EW parameters in MCFM/JHUGen
  void resetMCFM_EWKParameters(double ext_Gf, double ext_aemmz, double ext_mW, double ext_mZ, double ext_xW, int ext_ewscheme=3);


  MelaIO* getIORecord(); // Full parton-by-parton ME record
  MELACandidate* getCurrentCandidate();
  int getCurrentCandidateIndex();
  int getNCandidates();
  std::vector<MELATopCandidate*>* getTopCandidateCollection();


  float getConstant(bool useOldggZZConstants=false);
  float getConstant_JHUGenUndecayed();
  float getConstant_m4l(bool useOldggZZConstants=false);
  void get_PAux(float& prob); // SuperProb

  void computeDecayAngles(
    float& qH,
    float& m1,
    float& m2,
    float& costheta1,
    float& costheta2,
    float& Phi,
    float& costhetastar,
    float& Phi1
    );

  void computeP_selfDspin0(
    double selfDHvvcoupl_input[nSupportedHiggses][SIZE_HVV][2],
    float& prob,
    bool useConstant=true
    );
  void computeP_selfDspin1(
    double selfDZqqcoupl_input[SIZE_ZQQ][2],
    double selfDZvvcoupl_input[SIZE_ZVV][2],
    float& prob,
    bool useConstant=true
    );
  void computeP_selfDspin1(
    double selfDZvvcoupl_input[SIZE_ZVV][2],
    float& prob,
    bool useConstant=true
    );
  void computeP_selfDspin2(
    double selfDGggcoupl_input[SIZE_GGG][2],
    double selfDGqqcoupl_input[SIZE_GQQ][2],
    double selfDGvvcoupl_input[SIZE_GVV][2],
    float& prob,
    bool useConstant=true
    );
  void computeP_selfDspin2(
    double selfDGggcoupl_input[SIZE_GGG][2],
    double selfDGvvcoupl_input[SIZE_GVV][2],
    float& prob,
    bool useConstant=true
    );
  void computeP(
    double selfDHvvcoupl_freenorm_input[SIZE_HVV_FREENORM],
    float& prob,
    bool useConstant=true
    );
  void computeP(
    float& prob,
    bool useConstant=true
    );

  void computeD_CP(
    TVar::MatrixElement myME,
    TVar::Process myType,
    float& prob
    );

  //****VVH Spin-0****//
  void computeProdDecP(
    double selfDHvvcoupl_input[nSupportedHiggses][SIZE_HVV][2],
    double selfDHwwcoupl_input[nSupportedHiggses][SIZE_HVV][2],
    float& prob,
    bool useConstant=true
    );
  void computeProdDecP(
    float& prob,
    bool useConstant=true
    );

  //****HJ/HJJ/VBF Spin-0****//
  void computeProdP(
    double selfDHggcoupl_input[SIZE_HGG][2],
    double selfDHvvcoupl_input[nSupportedHiggses][SIZE_HVV][2],
    double selfDHwwcoupl_input[nSupportedHiggses][SIZE_HVV][2],
    float& prob,
    bool useConstant=true
    );
  void computeProdP(
    float& prob,
    bool useConstant=true
    );

  //****VH Spin-0****//
  void computeProdP_VH(
    double selfDHvvcoupl_input[nSupportedHiggses][SIZE_HVV][2],
    float& prob,
    bool includeHiggsDecay=false,
    bool useConstant=true
    );
  void computeProdP_VH(
    float& prob,
    bool includeHiggsDecay=false,
    bool useConstant=true
    );

  //***ttH Spin-0****//
  void computeProdP_ttH(
    float& prob,
    int topProcess=2,
    int topDecay=0,
    bool useConstant=true
    );

  // Calculation weight to correct for fermion interference
  void compute4FermionWeight(float& w);

  // Calculation of X propagator
  void getXPropagator(TVar::ResonancePropagatorScheme scheme, float& prop);

  //***SuperMela***//
  void computePM4l(
    TVar::SuperMelaSyst syst,
    float& prob
    );


  // Access ZZMEs Calculate4Momentum
  std::vector<TLorentzVector> calculate4Momentum(double Mx, double M1, double M2, double theta, double theta1, double theta2, double Phi1, double Phi);


  RooAbsPdf* pdf;
  ScalarPdfFactory_ggH* ggSpin0Model;
  VectorPdfFactory* spin1Model;
  TensorPdfFactory_HVV* spin2Model;
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
  RooRealVar* Y_rrv;
  RooRealVar* upFrac_rrv;

  TGraph* vaScale_4e;
  TGraph* vaScale_4mu;
  TGraph* vaScale_2e2mu;
  TGraph* DggZZ_scalefactor;
  void setCTotalBkgGraphs(TFile* fcontainer, TGraph* tgC[]);
  void constructDggr(
    float mzz,
    int flavor, // (1,2,3) for (4e,4mu,2e2mu)
    float bkg_VAMCFM_noscale,
    float ggzz_VAMCFM_noscale,
    float ggHZZ_prob_pure_noscale,
    float ggHZZ_prob_int_noscale,
    float& myDggr
    );

  TGraph* tgtotalbkg[3];
  void computeD_gg(
    TVar::MatrixElement myME,
    TVar::Process myType,
    float& prob
    );

  // Self-define arrays are now members of MELA.
  // There are a lot of them!
  //****Spin-0****//
  double selfDHvvcoupl_freenorm[SIZE_HVV_FREENORM];
  double selfDHqqcoupl[SIZE_HQQ][2];
  double selfDHggcoupl[SIZE_HGG][2];
  // The first dimension (of size [nSupportedHiggses=2]) supports a second resonance present in MCFM
  double selfDHzzcoupl[nSupportedHiggses][SIZE_HVV][2];
  double selfDHwwcoupl[nSupportedHiggses][SIZE_HVV][2];
  double selfDHzzLambda_qsq[nSupportedHiggses][4][3];
  double selfDHwwLambda_qsq[nSupportedHiggses][4][3];
  int selfDHzzCLambda_qsq[nSupportedHiggses][3];
  int selfDHwwCLambda_qsq[nSupportedHiggses][3];
  bool differentiate_HWW_HZZ;
  //****Spin-1****//
  double selfDZqqcoupl[SIZE_ZQQ][2];
  double selfDZvvcoupl[SIZE_ZVV][2];
  //****Spin-2****//
  double selfDGqqcoupl[SIZE_GQQ][2];
  double selfDGggcoupl[SIZE_GGG][2];
  double selfDGvvcoupl[SIZE_GVV][2];
  // That is a lot of them!

protected:

  //
  // Data memmbers
  //
  double LHCsqrts;
  TVar::Process myModel_;
  TVar::MatrixElement myME_;
  TVar::Production myProduction_;
  TVar::LeptonInterference myLepInterf_;
  TVar::VerbosityLevel myVerbosity_;

  newZZMatrixElement* ZZME;

  float auxiliaryProb;

  MELACandidate* melaCand; // Pointer to persistent TEvtProb object

  //
  // Functions
  //
  bool configureAnalyticalPDFs();
  void reset_SelfDCouplings();
  void reset_PAux(); // SuperProb reset
  void reset_CandRef();
};

#endif

