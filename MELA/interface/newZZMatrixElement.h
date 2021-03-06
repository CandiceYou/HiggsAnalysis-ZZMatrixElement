#ifndef newZZMatrixElement_newZZMatrixElement_h
#define newZZMatrixElement_newZZMatrixElement_h

#include <vector>
#include <ZZMatrixElement/MELA/interface/TVar.hh>
#include <ZZMatrixElement/MELA/interface/TEvtProb.hh>


class  newZZMatrixElement{
public:

  newZZMatrixElement(
    const char* pathtoPDFSet,
    int PDFMember,
    const char* pathtoHiggsCSandWidth, // path to the textfiles of the HiggsCSandWidth package
    double ebeam,
    TVar::VerbosityLevel verbosity
    );

  ~newZZMatrixElement();

  void computeXS(
    TVar::Process process_,
    TVar::MatrixElement me_,
    TVar::Production production_,
    float &mevalue
    );

  void computeProdXS_VVHVV(
    TVar::Process process_,
    TVar::MatrixElement me_,
    TVar::Production production_,
    float& mevalue
    );

  void computeProdXS_JJH(
    TVar::Process process_,
    TVar::MatrixElement me_,
    TVar::Production production_,
    float &mevalue
    );

  void computeProdXS_JH(
    TVar::Process process_,
    TVar::MatrixElement me_,
    TVar::Production production_,
    float &mevalue
    );

  void computeProdXS_VH(
    TVar::Process process_,
    TVar::MatrixElement me_,
    TVar::Production production_,
    float &mevalue,
    bool includeHiggsDecay=false
    );

  void computeProdXS_ttH(
    TVar::Process process_,
    TVar::MatrixElement me_,
    TVar::Production production_,
    float &mevalue,
    int topProcess,
    int topDecay=0
    );

  void get_XPropagator(TVar::ResonancePropagatorScheme scheme, float& prop);

  // Set-functions
  void set_Process(TVar::Process process_, TVar::MatrixElement me_, TVar::Production production_); // Sets variables in Xcal2 as well
  void set_Verbosity(TVar::VerbosityLevel verbosity_); // Sets variables in Xcal2 as well
  void set_LeptonInterference(TVar::LeptonInterference myLepInterf); // Sets variables in Xcal2 as well
  //
  void set_TempCandidate(
    SimpleParticleCollection_t* pDaughters,
    SimpleParticleCollection_t* pAssociated=0,
    SimpleParticleCollection_t* pMothers=0,
    bool isGen=false
    ); // Sets melaCand in Xcal2 to a temporary candidate, without pushing this candidate to candList of Xcal2 for storage and deletion at a later stage
  //
  void set_RenFacScaleMode(TVar::EventScaleScheme renormalizationSch, TVar::EventScaleScheme factorizationSch, double ren_sf, double fac_sf); // Sets variables exclusive to Xcal2
  void set_LHAgrid(const char* path, int pdfmember=0); // Sets variable exclusive to Xcal2
  void set_PrimaryHiggsMass(double mh);
  void set_CandidateDecayMode(TVar::CandidateDecayMode mode); // Sets variables exclusive to Xcal2
  void set_CurrentCandidateFromIndex(unsigned int icand); // Sets variables exclusive to Xcal2
  void set_CurrentCandidate(MELACandidate* cand); // Sets variables exclusive to Xcal2
  void set_InputEvent(
    SimpleParticleCollection_t* pDaughters,
    SimpleParticleCollection_t* pAssociated=0,
    SimpleParticleCollection_t* pMothers=0,
    bool isGen=false
    ); // Sets variables exclusive to Xcal2
  void append_TopCandidate(SimpleParticleCollection_t* TopDaughters); // Sets variable exclusive to Xcal2
  //
  void set_mHiggs(double mh_, int index); // Does not set any variables in Xcal2!
  void set_wHiggs(double gah_, int index); // Does not set any variables in Xcal2!
  void set_mHiggs_wHiggs(double mh_, double gah_, int index); // Does not set any variables in Xcal2!
  //

  // Reset-functions
  void reset_MCFM_EWKParameters(double ext_Gf, double ext_aemmz, double ext_mW, double ext_mZ, double ext_xW, int ext_ewscheme=3);
  void resetPerEvent(); // Resets variables and owned objects that are per-event
  void reset_InputEvent(); // Resets all candidates in Xcal2, to be called at the end of each event after all computations are done

  void set_SpinZeroCouplings(
    double selfDHvvcoupl_freenorm[SIZE_HVV_FREENORM],
    double selfDHqqcoupl[SIZE_HQQ][2],
    double selfDHggcoupl[SIZE_HGG][2],
    double selfDHzzcoupl[nSupportedHiggses][SIZE_HVV][2],
    double selfDHwwcoupl[nSupportedHiggses][SIZE_HVV][2],
    double selfDHzzLambda_qsq[nSupportedHiggses][4][3],
    double selfDHwwLambda_qsq[nSupportedHiggses][4][3],
    int selfDHzzCLambda_qsq[nSupportedHiggses][3],
    int selfDHwwCLambda_qsq[nSupportedHiggses][3],
    bool diffHWW = false
    );
  void set_SpinOneCouplings(
    double selfDZqqcoupl[SIZE_ZQQ][2],
    double selfDZvvcoupl[SIZE_ZVV][2]
    );
  void set_SpinTwoCouplings(
    double selfDGqqcoupl[SIZE_GQQ][2],
    double selfDGggcoupl[SIZE_GGG][2],
    double selfDGvvcoupl[SIZE_GVV][2]
    );

  // Compute four-momenta from angles - not cos(theta) - only 
  std::vector<TLorentzVector> Calculate4Momentum(double Mx, double M1, double M2, double theta, double theta1, double theta2, double Phi1, double Phi);

  // Get-functions
  MelaIO* get_IORecord();
  MELACandidate* get_CurrentCandidate();
  int get_CurrentCandidateIndex();
  int get_NCandidates();
  std::vector<MELATopCandidate*>* get_TopCandidateCollection(); // Just so that the user can set MELATopCandidate::passSelection=true or false to omit some tops, in case tere are more than two

protected:
  
  TVar::VerbosityLevel processVerbosity;
  TVar::LeptonInterference processLeptonInterference;
  TVar::Process processModel;
  TVar::MatrixElement processME;
  TVar::Production processProduction;

  double EBEAM;
  double mHiggs[nSupportedHiggses];
  double wHiggs[nSupportedHiggses];
  TEvtProb Xcal2;

  SpinZeroCouplings* selfD_SpinZeroCouplings;
  SpinOneCouplings* selfD_SpinOneCouplings;
  SpinTwoCouplings* selfD_SpinTwoCouplings;

  MELACandidate* melaCand; // Pointer to current candidate object of Xcal2
  std::vector<MELAParticle*> tmpPartList; // Vector of pointers to the owned, temporary MELAParticles
  // Having a temporary top candidate list does not make much sense at the moment
  //std::vector<MELATopCandidate*> tmpTopCandList; // Vector of pointers to the owned, temporary MELATopCandidates
  std::vector<MELACandidate*> tmpCandList; // Vector of pointers to the owned, temporary MELACandidates

};
#endif
