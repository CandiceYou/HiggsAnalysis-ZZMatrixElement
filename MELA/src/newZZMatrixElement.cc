#include <ZZMatrixElement/MELA/interface/newZZMatrixElement.h>
#include "TLorentzRotation.h"


using namespace std;
using namespace TUtil;


newZZMatrixElement::newZZMatrixElement(
  const char* pathtoPDFSet,
  int PDFMember,
  const char* pathtoHiggsCSandWidth,
  double ebeam,
  TVar::VerbosityLevel verbosity
  ) :
  processVerbosity(verbosity),
  processLeptonInterference(TVar::DefaultLeptonInterf),
  EBEAM(ebeam),
  Xcal2(pathtoHiggsCSandWidth, ebeam, pathtoPDFSet, PDFMember, verbosity),
  melaCand(0)
{
  if (processVerbosity>=TVar::DEBUG) cout << "Begin newZZMatrixElement constructor" << endl;

  // Set default parameters explicitly
  set_PrimaryHiggsMass(125.);
  set_mHiggs(125., 0); set_wHiggs(-1., 0);
  set_mHiggs(-1., 1); set_wHiggs(0., 1);

  selfD_SpinZeroCouplings = Xcal2.GetSelfDSpinZeroCouplings();
  selfD_SpinOneCouplings = Xcal2.GetSelfDSpinOneCouplings();
  selfD_SpinTwoCouplings = Xcal2.GetSelfDSpinTwoCouplings();

  if (processVerbosity>=TVar::DEBUG) cout << "End newZZMatrixElement constructor" << endl;
}
newZZMatrixElement::~newZZMatrixElement(){
  if (processVerbosity>=TVar::DEBUG) cout << "Begin newZZMatrixElement destructor" << endl;
  /*std::cout << "End of newZZME" << std::endl;*/
  resetPerEvent();
  reset_InputEvent();
  if (processVerbosity>=TVar::DEBUG) cout << "End newZZMatrixElement destructor" << endl;
}


vector<TLorentzVector> newZZMatrixElement::Calculate4Momentum(double Mx, double M1, double M2, double theta, double theta1, double theta2, double Phi1, double Phi){
  double phi1, phi2;
  phi1=TMath::Pi()-Phi1;
  phi2=Phi1+Phi;

  double gamma1=1, gamma2=1, beta1=0, beta2=0;

  if (M1>0. && Mx>0.){
    gamma1=(Mx*Mx+M1*M1-M2*M2)/(2*Mx*M1);
    beta1=sqrt(1.-1./(gamma1*gamma1));
  }
  if (M2>0. && Mx>0.){
    gamma2=(Mx*Mx-M1*M1+M2*M2)/(2*Mx*M2);
    beta2=sqrt(1.-1./(gamma2*gamma2));
  }

  //gluon 4 vectors
  TLorentzVector p1CM(0, 0, Mx/2, Mx/2);
  TLorentzVector p2CM(0, 0, -Mx/2, Mx/2);

  //vector boson 4 vectors
  TLorentzVector kZ1(gamma1*M1*sin(theta)*beta1, 0, gamma1*M1*cos(theta)*beta1, gamma1*M1);
  TLorentzVector kZ2(-gamma2*M2*sin(theta)*beta2, 0, -gamma2*M2*cos(theta)*beta2, gamma2*M2);

  //Rotation and Boost matrices. Note gamma1*beta1*M1=gamma2*beta2*M2.

  TLorentzRotation Z1ToZ, Z2ToZ;

  Z1ToZ.Boost(0, 0, beta1);
  Z2ToZ.Boost(0, 0, beta2);
  Z1ToZ.RotateY(theta);
  Z2ToZ.RotateY(TMath::Pi()+theta);


  //fermion 4 vectors in vector boson rest frame

  TLorentzVector p3Z1((M1/2)*sin(theta1)*cos(phi1), (M1/2)*sin(theta1)*sin(phi1), (M1/2)*cos(theta1), (M1/2)*1);
  TLorentzVector p4Z1(-(M1/2)*sin(theta1)*cos(phi1), -(M1/2)*sin(theta1)*sin(phi1), -(M1/2)*cos(theta1), (M1/2)*1);
  TLorentzVector p5Z2((M2/2)*sin(theta2)*cos(phi2), (M2/2)*sin(theta2)*sin(phi2), (M2/2)*cos(theta2), (M2/2)*1);
  TLorentzVector p6Z2(-(M2/2)*sin(theta2)*cos(phi2), -(M2/2)*sin(theta2)*sin(phi2), -(M2/2)*cos(theta2), (M2/2)*1);

  // fermions 4 vectors in CM frame

  TLorentzVector p3CM, p4CM, p5CM, p6CM;

  p3CM=Z1ToZ*p3Z1;
  p4CM=Z1ToZ*p4Z1;
  p5CM=Z2ToZ*p5Z2;
  p6CM=Z2ToZ*p6Z2;

  vector<TLorentzVector> p;

  p.push_back(p3CM);
  p.push_back(p4CM);
  p.push_back(p5CM);
  p.push_back(p6CM);

  return p;
}

// Set-functions that set variables that belong to Xcal2 in addition to setting variables that belong to newZZMatrixElement
void newZZMatrixElement::set_Process(TVar::Process process_, TVar::MatrixElement me_, TVar::Production production_){
  processModel = process_; Xcal2.SetProcess(processModel);
  processProduction = production_; Xcal2.SetProduction(processProduction);
  processME = me_; Xcal2.SetMatrixElement(processME);
}
void newZZMatrixElement::set_Verbosity(TVar::VerbosityLevel verbosity_){ processVerbosity = verbosity_; Xcal2.SetVerbosity(verbosity_); }
void newZZMatrixElement::set_LeptonInterference(TVar::LeptonInterference lepInterf_){ processLeptonInterference = lepInterf_; Xcal2.SetLeptonInterf(processLeptonInterference); }
// Set-functions that set variables that exclusively belong to Xcal2
void newZZMatrixElement::set_LHAgrid(const char* path, int pdfmember){ Xcal2.Set_LHAgrid(path, pdfmember); }
void newZZMatrixElement::set_RenFacScaleMode(TVar::EventScaleScheme renormalizationSch, TVar::EventScaleScheme factorizationSch, double ren_sf, double fac_sf){
  Xcal2.SetRenFacScaleMode(renormalizationSch, factorizationSch, ren_sf, fac_sf);
}
void newZZMatrixElement::set_CandidateDecayMode(TVar::CandidateDecayMode mode){ Xcal2.SetCandidateDecayMode(mode); }
void newZZMatrixElement::set_PrimaryHiggsMass(double mh){ Xcal2.SetPrimaryHiggsMass(mh); }
void newZZMatrixElement::set_CurrentCandidateFromIndex(unsigned int icand){ Xcal2.SetCurrentCandidateFromIndex(icand); }
void newZZMatrixElement::set_CurrentCandidate(MELACandidate* cand){ Xcal2.SetCurrentCandidate(cand); }
void newZZMatrixElement::set_InputEvent(
  SimpleParticleCollection_t* pDaughters,
  SimpleParticleCollection_t* pAssociated,
  SimpleParticleCollection_t* pMothers,
  bool isGen
  ){ // Adds a new candidate to Xcal2
  Xcal2.SetInputEvent(
    pDaughters,
    pAssociated,
    pMothers,
    isGen
    );
}
// Sets melaCand in Xcal2 to a temporary candidate, without pushing this candidate to candList of Xcal2
void newZZMatrixElement::set_TempCandidate(
  SimpleParticleCollection_t* pDaughters,
  SimpleParticleCollection_t* pAssociated,
  SimpleParticleCollection_t* pMothers,
  bool isGen
  ){
  MELACandidate* cand = ConvertVectorFormat(
    pDaughters,
    pAssociated,
    pMothers,
    isGen,
    &tmpPartList, &tmpCandList // push_back is done automatically
    );
  if (cand!=0){
    melaCand=cand;
    set_CurrentCandidate(melaCand);
  }
}
// Adds a top candidate
void newZZMatrixElement::append_TopCandidate(SimpleParticleCollection_t* TopDaughters){ Xcal2.AppendTopCandidate(TopDaughters); }
// Set-functions that do not set anything that belongs to Xcal2
void newZZMatrixElement::set_mHiggs(double mh_, int index){
  if (index<nSupportedHiggses && index>=0) mHiggs[index] = mh_;
  else cerr << "newZZMatrixElement::set_mHiggs: Only resonances 0 (regular) and 1 (additional, possibly high-mass) are supported" << endl;
}
void newZZMatrixElement::set_wHiggs(double gah_, int index){
  if (index<nSupportedHiggses && index>=0) wHiggs[index] = (double)gah_;
  else cerr << "newZZMatrixElement::set_wHiggs: Only resonances 0 (regular) and 1 (additional, possibly high-mass) are supported" << endl;
}
void newZZMatrixElement::set_mHiggs_wHiggs(double mh_, double gah_, int index){
  if (index<nSupportedHiggses && index>=0){
    mHiggs[index] = mh_;
    wHiggs[index] = gah_;
  }
  else cerr << "newZZMatrixElement::set_mHiggs_wHiggs: Only resonances 0 (regular) and 1 (additional, possibly high-mass) are supported" << endl;
}
// reset_MCFM_EWKParameters resets the MCFM EW parameters to those specified. This function is a wrapper around the TEvtProb version.
void newZZMatrixElement::reset_MCFM_EWKParameters(double ext_Gf, double ext_aemmz, double ext_mW, double ext_mZ, double ext_xW, int ext_ewscheme){
  Xcal2.ResetMCFM_EWKParameters(ext_Gf, ext_aemmz, ext_mW, ext_mZ, ext_xW, ext_ewscheme);
}
//
// resetPerEvent resets the mass, width and lepton interference settings and deletes the temporary input objects newZZMatrixElement owns.
void newZZMatrixElement::resetPerEvent(){
  // Protection against forgetfulness; custom width has to be set per-computation
  set_mHiggs(Xcal2.GetPrimaryHiggsMass(), 0); // Sets mHiggs[0]
  if (wHiggs[0]>=0.) set_wHiggs(-1., 0);

  if (mHiggs[1]>=0.) set_mHiggs(-1., 1);
  if (wHiggs[1]>0.) set_wHiggs(0., 1);

  Xcal2.SetHiggsMass(mHiggs[0], -1, -1);

  // Return back to default lepton interference settings after each calculation
  if (processLeptonInterference!=TVar::DefaultLeptonInterf) set_LeptonInterference(TVar::DefaultLeptonInterf);

  // Delete the temporary input objects owned
  for (unsigned int ic=0; ic<tmpCandList.size(); ic++){ if (tmpCandList.at(ic)!=0) delete tmpCandList.at(ic); }
  //for (unsigned int itc=0; itc<tmpTopCandList.size(); itc++){ if (tmpTopCandList.at(itc)!=0) delete tmpTopCandList.at(itc); }
  for (unsigned int ip=0; ip<tmpPartList.size(); ip++){ if (tmpPartList.at(ip)!=0) delete tmpPartList.at(ip); }
  melaCand=0;
}
// Resets all candidates in Xcal2, to be called at the end of each event after all computations are done
void newZZMatrixElement::reset_InputEvent(){ Xcal2.ResetInputEvent(); }


MelaIO* newZZMatrixElement::get_IORecord(){ return Xcal2.GetIORecord(); }
MELACandidate* newZZMatrixElement::get_CurrentCandidate(){ return Xcal2.GetCurrentCandidate(); }
int newZZMatrixElement::get_CurrentCandidateIndex(){ return Xcal2.GetCurrentCandidateIndex(); }
int newZZMatrixElement::get_NCandidates(){ return Xcal2.GetNCandidates(); }
vector<MELATopCandidate*>* newZZMatrixElement::get_TopCandidateCollection(){ return Xcal2.GetTopCandidates(); }


void newZZMatrixElement::set_SpinZeroCouplings(
  double selfDHvvcoupl_freenorm[SIZE_HVV_FREENORM],
  double selfDHqqcoupl[SIZE_HQQ][2],
  double selfDHggcoupl[SIZE_HGG][2],
  double selfDHzzcoupl[nSupportedHiggses][SIZE_HVV][2],
  double selfDHwwcoupl[nSupportedHiggses][SIZE_HVV][2],
  double selfDHzzLambda_qsq[nSupportedHiggses][4][3],
  double selfDHwwLambda_qsq[nSupportedHiggses][4][3],
  int selfDHzzCLambda_qsq[nSupportedHiggses][3],
  int selfDHwwCLambda_qsq[nSupportedHiggses][3],
  bool diffHWW
  ){
  Xcal2.AllowSeparateWWCouplings(diffHWW);
  for (int ic=0; ic<SIZE_HVV_FREENORM; ic++) selfD_SpinZeroCouplings->SetHVVFreeNormCouplings(ic, selfDHvvcoupl_freenorm[ic]);
  for (int ic=0; ic<SIZE_HQQ; ic++) selfD_SpinZeroCouplings->SetHQQCouplings(ic, selfDHqqcoupl[ic][0], selfDHqqcoupl[ic][1]);
  for (int ic=0; ic<SIZE_HGG; ic++) selfD_SpinZeroCouplings->SetHGGCouplings(ic, selfDHggcoupl[ic][0], selfDHggcoupl[ic][1]);

  for (int jh=1; jh<=(int)nSupportedHiggses; jh++){
    for (int ic=0; ic<SIZE_HVV; ic++) selfD_SpinZeroCouplings->SetHVVCouplings(ic, selfDHzzcoupl[jh-1][ic][0], selfDHzzcoupl[jh-1][ic][1], false, jh);
    for (int ic=0; ic<SIZE_HVV; ic++) selfD_SpinZeroCouplings->SetHVVCouplings(ic, selfDHwwcoupl[jh-1][ic][0], selfDHwwcoupl[jh-1][ic][1], true, jh);
    for (int ik=0; ik<3; ik++){
      for (int ig=0; ig<4; ig++){
        selfD_SpinZeroCouplings->SetHVVLambdaQ2(ig, ik, selfDHzzLambda_qsq[jh-1][ig][ik], false, jh);
        selfD_SpinZeroCouplings->SetHVVLambdaQ2(ig, ik, selfDHwwLambda_qsq[jh-1][ig][ik], true, jh);
      }
      selfD_SpinZeroCouplings->SetHVVSignCQ2(ik, selfDHzzCLambda_qsq[jh-1][ik], false, jh);
      selfD_SpinZeroCouplings->SetHVVSignCQ2(ik, selfDHwwCLambda_qsq[jh-1][ik], true, jh);
    }
  }
}
void newZZMatrixElement::set_SpinOneCouplings(
  double selfDZqqcoupl[SIZE_ZQQ][2],
  double selfDZvvcoupl[SIZE_ZVV][2]
  ){
  for (int ic=0; ic<SIZE_ZQQ; ic++) selfD_SpinOneCouplings->SetZQQCouplings(ic, selfDZqqcoupl[ic][0], selfDZqqcoupl[ic][1]);
  for (int ic=0; ic<SIZE_ZVV; ic++) selfD_SpinOneCouplings->SetZVVCouplings(ic, selfDZvvcoupl[ic][0], selfDZvvcoupl[ic][1]);
}
void newZZMatrixElement::set_SpinTwoCouplings(
  double selfDGqqcoupl[SIZE_GQQ][2],
  double selfDGggcoupl[SIZE_GGG][2],
  double selfDGvvcoupl[SIZE_GVV][2]
  ){
  for (int ic=0; ic<SIZE_GQQ; ic++) selfD_SpinTwoCouplings->SetGQQCouplings(ic, selfDGqqcoupl[ic][0], selfDGqqcoupl[ic][1]);
  for (int ic=0; ic<SIZE_GGG; ic++) selfD_SpinTwoCouplings->SetGGGCouplings(ic, selfDGggcoupl[ic][0], selfDGggcoupl[ic][1]);
  for (int ic=0; ic<SIZE_GVV; ic++) selfD_SpinTwoCouplings->SetGVVCouplings(ic, selfDGvvcoupl[ic][0], selfDGvvcoupl[ic][1]);
}



// Higgs + 0 jets dedicated function (with no Higgs decay)
void newZZMatrixElement::computeXS(
  TVar::Process process_,
  TVar::MatrixElement me_,
  TVar::Production production_,
  float &mevalue
  ){
  set_Process(process_, me_, production_);
  melaCand = get_CurrentCandidate();

  if (melaCand!=0){
    double zzmass = melaCand->m();
    if (
      processME==TVar::MCFM
      ||
      (processModel==TVar::bkgZZ_SMHiggs && processProduction == TVar::ZZGG && processME==TVar::JHUGen)
      ){
      for (int jh=0; jh<(int)nSupportedHiggses; jh++) Xcal2.SetHiggsMass(mHiggs[jh], wHiggs[jh], jh+1);
    }
    else Xcal2.SetHiggsMass(zzmass, wHiggs[0], -1);

    mevalue = Xcal2.XsecCalc_XVV(
      processModel, processProduction,
      processVerbosity
      );
  }

  resetPerEvent();
  return;
}

// VBF+VH dedicated function (production(+)decay)
void newZZMatrixElement::computeProdXS_VVHVV(
  TVar::Process process_,
  TVar::MatrixElement me_,
  TVar::Production production_,
  float& mevalue
  ){
  set_Process(process_, me_, production_);
  melaCand = get_CurrentCandidate();

  if (melaCand!=0){
    double zzmass = melaCand->m();
    if (
      processME==TVar::MCFM
      ||
      (processModel==TVar::bkgZZ_SMHiggs && processProduction == TVar::ZZGG && processME==TVar::JHUGen)
      ){
      for (int jh=0; jh<(int)nSupportedHiggses; jh++) Xcal2.SetHiggsMass(mHiggs[jh], wHiggs[jh], jh+1);
    }
    else Xcal2.SetHiggsMass(zzmass, wHiggs[0], -1);

    mevalue = Xcal2.XsecCalc_VVXVV(
      processModel, processProduction,
      processVerbosity
      );
  }

  resetPerEvent();
  return;
}

// Higgs + 2 jet dedicated function (with no Higgs decay)
void newZZMatrixElement::computeProdXS_JJH(
  TVar::Process process_,
  TVar::MatrixElement me_,
  TVar::Production production_,
  float &mevalue
  ){
  set_Process(process_, me_, production_);
  melaCand = get_CurrentCandidate();

  if (melaCand!=0){
    mevalue  = Xcal2.XsecCalcXJJ(
      processModel, processProduction,
      processVerbosity
      );
  }

  resetPerEvent();
  return;
}

// Higgs + 1 jet: Only SM is supported for now.
void newZZMatrixElement::computeProdXS_JH(
  TVar::Process process_,
  TVar::MatrixElement me_,
  TVar::Production production_,
  float &mevalue
  ){
  set_Process(process_, me_, production_);
  melaCand = get_CurrentCandidate();

  if (melaCand!=0){
    mevalue  = Xcal2.XsecCalcXJ(
      processModel, processProduction,
      processVerbosity
      );
  }

  resetPerEvent();
  return;
}

// Dedicated VH function (with no Higgs decay)
void newZZMatrixElement::computeProdXS_VH(
  TVar::Process process_,
  TVar::MatrixElement me_,
  TVar::Production production_,
  float &mevalue,
  bool includeHiggsDecay
  ){
  set_Process(process_, me_, production_);
  melaCand = get_CurrentCandidate();

  if (melaCand!=0){
    double zzmass = melaCand->m();
    if (processME==TVar::MCFM){ for (int jh=0; jh<(int)nSupportedHiggses; jh++) Xcal2.SetHiggsMass(mHiggs[jh], wHiggs[jh], jh+1); }
    else Xcal2.SetHiggsMass(zzmass, wHiggs[0], -1);

    mevalue  = Xcal2.XsecCalc_VX(
      processModel, processProduction,
      processVerbosity,
      includeHiggsDecay
      );
  }

  resetPerEvent();
  return;
}

// Dedicated ttH/bbH function (with no Higgs decay)
// topProcess=0 for gg, =1 for qqb initial states
// topDecay=0 (default) for no top decay, =1 to include t/tb->b/bb + W+/W-(->ffb). =1 not relevant for the bbH process.
void newZZMatrixElement::computeProdXS_ttH(
  TVar::Process process_,
  TVar::MatrixElement me_,
  TVar::Production production_,
  float &mevalue,
  int topProcess,
  int topDecay
  ){
  set_Process(process_, me_, production_);
  melaCand = get_CurrentCandidate();

  if (melaCand!=0){
    double zzmass = melaCand->m();
    if (processME==TVar::MCFM){ for (int jh=0; jh<(int)nSupportedHiggses; jh++) Xcal2.SetHiggsMass(mHiggs[jh], wHiggs[jh], jh+1); }
    else Xcal2.SetHiggsMass(zzmass, wHiggs[0], -1);

    if (production_ == TVar::ttH || production_ ==TVar::bbH){
      mevalue  = Xcal2.XsecCalc_TTX(
        processModel, processProduction,
        processVerbosity,
        topProcess, topDecay
        );
    }
  }

  resetPerEvent();
  return;
}

// Higgs propagator
void newZZMatrixElement::get_XPropagator(TVar::ResonancePropagatorScheme scheme, float& prop){
  prop=0.;
  melaCand = get_CurrentCandidate();

  if (melaCand!=0) prop=Xcal2.GetXPropagator(scheme);

  resetPerEvent();
  return;
}


