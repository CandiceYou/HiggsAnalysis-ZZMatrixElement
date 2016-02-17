//-----------------------------------------------------------------------------
//
// Class EventProb Module
//
//   EventProb Module
//
// March 21 2011
// S. Jindariani (sergo@fnal.gov)
// Y. Gao (ygao@fnal.gov)
// K. Burkett (burkett@fnal.gov)
//-----------------------------------------------------------------------------

#include "ZZMatrixElement/MELA/interface/TVar.hh"
#include "ZZMatrixElement/MELA/interface/TEvtProb.hh"


ClassImp(TEvtProb)

using namespace std;

//-----------------------------------------------------------------------------
// Constructors and Destructor
//-----------------------------------------------------------------------------
TEvtProb::TEvtProb(
  const char* path, double ebeam, const char* pathtoPDFSet, int PDFMember
  ) :
  EBEAM(ebeam){
  mcfm_init_((char *)"input.DAT", (char *)"./");
  SetEwkCouplingParameters();
  energy_.sqrts = 2.*EBEAM;
  coupling_();
  string path_string = path;
  myCSW_ = new HiggsCSandWidth_MELA(path_string);
  //std::cout << path << std::endl;
  SetLeptonInterf(TVar::DefaultLeptonInterf);
  spinzerohiggs_anomcoupl_.LambdaBSM=1000;
  spinzerohiggs_anomcoupl_.Lambda_z1=10000;
  spinzerohiggs_anomcoupl_.Lambda_z2=10000;
  spinzerohiggs_anomcoupl_.Lambda_z3=10000;
  spinzerohiggs_anomcoupl_.Lambda_z4=10000;
  spinzerohiggs_anomcoupl_.Lambda_Q=10000;

  InitJHUGenMELA(pathtoPDFSet, PDFMember);

  ResetCouplings();
  ResetRenFacScaleMode();
}


TEvtProb::~TEvtProb(){
  delete myCSW_;
}

/*
void TEvtProb::ResetMCFM_EWKParameters(double ext_Gf, double ext_aemmz, double ext_mW, double ext_mZ){
	ewinput_.Gf_inp = ext_Gf;
	ewinput_.aemmz_inp = ext_aemmz;
	ewinput_.wmass_inp = ext_mW;
	ewinput_.zmass_inp = ext_mZ;
	ewinput_.xw_inp = 1.-pow(ext_mW/ext_mZ,2);
	coupling_();
}
*/

void TEvtProb::ResetMCFM_EWKParameters(double ext_Gf, double ext_aemmz, double ext_mW, double ext_mZ, double ext_xW, int ext_ewscheme){
  if (ext_ewscheme<-1 || ext_ewscheme>3) ext_ewscheme=3;
  ewinput_.Gf_inp = ext_Gf;
  ewinput_.aemmz_inp = ext_aemmz;
  ewinput_.wmass_inp = ext_mW;
  ewinput_.zmass_inp = ext_mZ;
  ewinput_.xw_inp = ext_xW;
  ewscheme_.ewscheme = ext_ewscheme;
  coupling_();
}

// Set NNPDF driver path
void TEvtProb::Set_LHAgrid(const char* path, int pdfmember){
  char path_nnpdf_c[200];
  sprintf(path_nnpdf_c, "%s", path);
  int pathLength = strlen(path_nnpdf_c);
  nnpdfdriver_(path_nnpdf_c, &pathLength);
  nninitpdf_(&pdfmember);
}

//
// Directly calculate the ZZ->4l differential cross-section 
// 
double TEvtProb::XsecCalc(
  TVar::Process proc, TVar::Production production, const hzz4l_event_type &hzz4l_event,
  TVar::VerbosityLevel verbosity
  ){
  ResetIORecord();

  //Initialize Process
  SetProcess(proc);
  SetProduction(production);
  bool forceUseMCFM = (_matrixElement == TVar::MCFM || _process == TVar::bkgZZ_SMHiggs);

  int flavor = abs(hzz4l_event.PdgCode[0]) == abs(hzz4l_event.PdgCode[2]) ? 1 : 3;
  bool needBSMHiggs=false;
  if (forceUseMCFM){
    for (int vv = 0; vv < SIZE_HVV; vv++){
      if ((selfDSpinZeroCoupl.Hzzcoupl)[vv][1] != 0 || (selfDSpinZeroCoupl.Hzzcoupl)[vv][0] != 0){ needBSMHiggs = true; break; } // Only possible if selfDefine is called.
    }
    if (needBSMHiggs) SetLeptonInterf(TVar::InterfOn); // All anomalous coupling computations have to use lepton interference
    SetMCFMSpinZeroVVCouplings(needBSMHiggs, (selfDSpinZeroCoupl.Hzzcoupl), (selfDSpinZeroCoupl.Hzzcoupl));
    My_choose(_process, _production, _leptonInterf, flavor);
  }

  //constants
  double sqrts = 2.*EBEAM;
  double W=sqrts*sqrts;

  //Weight calculation
  double flux=1.;
  double dXsec=0.;

  mcfm_event_type mcfm_event;
  // assign the right initial momentum
  // assumes the events are boosted to have 0 transverse momenta
  TLorentzVector totalMom = hzz4l_event.p[0] + hzz4l_event.p[1] + hzz4l_event.p[2] + hzz4l_event.p[3];
  double sysPz= totalMom.Z();
  double sysE = totalMom.T();
  mcfm_event.p[2].SetPxPyPzE(hzz4l_event.p[0].Px(), hzz4l_event.p[0].Py(), hzz4l_event.p[0].Pz(), hzz4l_event.p[0].Energy());
  mcfm_event.p[3].SetPxPyPzE(hzz4l_event.p[1].Px(), hzz4l_event.p[1].Py(), hzz4l_event.p[1].Pz(), hzz4l_event.p[1].Energy());
  mcfm_event.p[4].SetPxPyPzE(hzz4l_event.p[2].Px(), hzz4l_event.p[2].Py(), hzz4l_event.p[2].Pz(), hzz4l_event.p[2].Energy());
  mcfm_event.p[5].SetPxPyPzE(hzz4l_event.p[3].Px(), hzz4l_event.p[3].Py(), hzz4l_event.p[3].Pz(), hzz4l_event.p[3].Energy());

  mcfm_event.PdgCode[0] = 21;
  mcfm_event.PdgCode[1] = 21;
  mcfm_event.PdgCode[2] = hzz4l_event.PdgCode[0];
  mcfm_event.PdgCode[3] = hzz4l_event.PdgCode[1];
  mcfm_event.PdgCode[4] = hzz4l_event.PdgCode[2];
  mcfm_event.PdgCode[5] = hzz4l_event.PdgCode[3];

  //Matrix Element evaluation in qX=qY=0 frame
  //Evaluate f(x1)f(x2)|M(q)|/x1/x2 
  // 
  double qX = totalMom.X();
  double qY = totalMom.Y();
  double qE = sysE;

  if ((qX*qX+qY*qY)>0){
    TVector3 boostV(qX/qE, qY/qE, 0);
    for (int ipt=2; ipt<6; ipt++) mcfm_event.p[ipt].Boost(-boostV);
    totalMom.Boost(-boostV);
  }

  sysPz= totalMom.Z();
  sysE = totalMom.T();
  double pz0 = (sysE+sysPz)/2.;
  double pz1 = -(sysE-sysPz)/2.;
  mcfm_event.p[0].SetPxPyPzE(0., 0., pz0, TMath::Abs(pz0));
  mcfm_event.p[1].SetPxPyPzE(0., 0., pz1, TMath::Abs(pz1));

  //event selections in Lab Frame
  double msqjk=0;
  if (forceUseMCFM) msqjk = SumMatrixElementPDF(_process, _production, _matrixElement, &event_scales, &mcfm_event, &RcdME, &flux, EBEAM, (selfDSpinZeroCoupl.Hvvcoupl_freenorm));
  else if (_matrixElement == TVar::JHUGen){
    AllowSeparateWWCouplings(false); // HZZ couplings are used for both in spin-0
    // all the possible couplings
    double Hggcoupl[SIZE_HGG][2] ={ { 0 } };
    double Hvvcoupl[SIZE_HVV][2] ={ { 0 } };
    double HvvLambda_qsq[4][3] ={ { 0 } };
    int HvvCLambda_qsq[3] ={ 0 };

    double Zqqcoupl[SIZE_ZQQ][2] ={ { 0 } };
    double Zvvcoupl[SIZE_ZVV][2] ={ { 0 } };

    double Gqqcoupl[SIZE_GQQ][2] ={ { 0 } };
    double Gggcoupl[SIZE_GGG][2] ={ { 0 } };
    double Gvvcoupl[SIZE_GVV][2] ={ { 0 } };

    // 
    // set spin 0 default numbers
    // 
    // By default set the Spin 0 couplings for SM case (0+m)
    Hggcoupl[0][0]=1.0;  Hggcoupl[0][1]=0.0;
    Hvvcoupl[0][0]=1.0;  Hvvcoupl[0][1]=0.0;
    for (int ic=0; ic<4; ic++){
      for (int ik=0; ik<3; ik++) HvvLambda_qsq[ic][ik]=100.;
    }
    // 
    // set spin 2 default numbers (2+m)
    // 
    Gqqcoupl[0][0]=1.0;  Gqqcoupl[0][1]=0.0;
    Gqqcoupl[1][0]=1.0;  Gqqcoupl[1][1]=0.0;
    Gggcoupl[0][0]=1.0;  Gggcoupl[0][1]=0.0;
    Gvvcoupl[0][0]=1.0;  Gvvcoupl[0][1]=0.0;
    Gvvcoupl[4][0]=1.0;  Gvvcoupl[4][1]=0.0;
    //
    // set spin 1 default numbers (1-)
    //
    Zqqcoupl[0][0]=1.0;  Zqqcoupl[0][1]=0.0;
    Zqqcoupl[1][0]=1.0;  Zqqcoupl[1][1]=0.0;
    Zvvcoupl[0][0]=1.0;  Zvvcoupl[0][1]=0.0; // 1-

    bool isSpinZero = false;
    bool isSpinOne = false;
    bool isSpinTwo = false;

    if (_process == TVar::HSMHiggs) isSpinZero = true; // Already the default
    // 0-
    else if (_process == TVar::H0minus) {
      Hvvcoupl[0][0] = 0.0;
      Hvvcoupl[1][0] = 0.0;
      Hvvcoupl[2][0] = 0.0;
      Hvvcoupl[3][0] = 1.0;
      isSpinZero = true;
    }
    // 0h+
    else if (_process == TVar::H0hplus) {
      Hvvcoupl[0][0] = 0.0;
      Hvvcoupl[1][0] = 1.0;
      Hvvcoupl[2][0] = 0.0;
      Hvvcoupl[3][0] = 0.0;
      isSpinZero = true;
    }
    // 0+L1
    else if (_process == TVar::H0_g1prime2){
      Hvvcoupl[0][0] = 0.;
      Hvvcoupl[11][0] = -12046.01;
      isSpinZero = true;
    }
    else if (_process == TVar::H0_Zgs) {
      Hvvcoupl[4][0] = 0.0688;
      Hvvcoupl[0][0] = 0.;
      isSpinZero = true;
    }
    else if (_process == TVar::H0_gsgs) {
      Hvvcoupl[7][0] = -0.0898;
      Hvvcoupl[0][0] = 0.;
      isSpinZero = true;
    }
    else if (_process == TVar::H0_Zgs_PS) {
      Hvvcoupl[6][0] = 0.0855;
      Hvvcoupl[0][0] = 0.;
      isSpinZero = true;
    }
    else if (_process == TVar::H0_gsgs_PS) {
      Hvvcoupl[9][0] = -0.0907;
      Hvvcoupl[0][0] = 0.;
      isSpinZero = true;
    }
    else if (_process == TVar::H0_Zgsg1prime2) {
      Hvvcoupl[30][0] = -7591.914; // +- 6.613
      Hvvcoupl[0][0] = 0.;
      isSpinZero = true;
    }
    else if (_process == TVar::SelfDefine_spin0){
      for (int j=0; j<2; j++){
        for (int i=0; i<SIZE_HGG; i++) Hggcoupl[i][j] = (selfDSpinZeroCoupl.Hggcoupl)[i][j];
        for (int i=0; i<SIZE_HVV; i++) Hvvcoupl[i][j] = (selfDSpinZeroCoupl.Hzzcoupl)[i][j];
      }
      for (int j=0; j<3; j++){
        for (int i=0; i<4; i++) HvvLambda_qsq[i][j] = (selfDSpinZeroCoupl.HzzLambda_qsq)[i][j];
        HvvCLambda_qsq[j] = (selfDSpinZeroCoupl.HzzCLambda_qsq)[j];
      }
      isSpinZero = true;
    }


    else if (_process == TVar::H2_g1g5) isSpinTwo = true; // Already the default
    // 2h-
    else if (_process == TVar::H2_g8){
      // gg production coupling constants
      Gggcoupl[0][0]=0.0;  Gggcoupl[0][1]=0.0;
      Gggcoupl[1][0]=0.0;  Gggcoupl[1][1]=0.0;
      Gggcoupl[2][0]=0.0;  Gggcoupl[2][1]=0.0;
      Gggcoupl[3][0]=0.0;  Gggcoupl[3][1]=0.0;
      Gggcoupl[4][0]=1.0;  Gggcoupl[4][1]=0.0; // 2h-

      // Graviton->ZZ coupling constants 
      Gvvcoupl[0][0]=0.0;  Gvvcoupl[0][1]=0.0;
      Gvvcoupl[1][0]=0.0;  Gvvcoupl[1][1]=0.0;
      Gvvcoupl[2][0]=0.0;  Gvvcoupl[2][1]=0.0;
      Gvvcoupl[3][0]=0.0;  Gvvcoupl[3][1]=0.0;
      Gvvcoupl[4][0]=0.0;  Gvvcoupl[4][1]=0.0;
      Gvvcoupl[5][0]=0.0;  Gvvcoupl[5][1]=0.0;
      Gvvcoupl[6][0]=0.0;  Gvvcoupl[6][1]=0.0;
      Gvvcoupl[7][0]=1.0;  Gvvcoupl[7][1]=0.0; // 2h-
      Gvvcoupl[8][0]=0.0;  Gvvcoupl[8][1]=0.0;
      Gvvcoupl[9][0]=0.0;  Gvvcoupl[9][1]=0.0;

      isSpinTwo = true;
    }
    // 2h+
    else if (_process == TVar::H2_g4){
      // gg production coupling constants
      Gggcoupl[0][0]=0.0;  Gggcoupl[0][1]=0.0;
      Gggcoupl[1][0]=0.0;  Gggcoupl[1][1]=0.0;
      Gggcoupl[2][0]=0.0;  Gggcoupl[2][1]=0.0;
      Gggcoupl[3][0]=1.0;  Gggcoupl[3][1]=0.0; // 2h+
      Gggcoupl[4][0]=0.0;  Gggcoupl[4][1]=0.0;

      // Graviton->ZZ coupling constants 
      Gvvcoupl[0][0]=0.0;  Gvvcoupl[0][1]=0.0;
      Gvvcoupl[1][0]=0.0;  Gvvcoupl[1][1]=0.0;
      Gvvcoupl[2][0]=0.0;  Gvvcoupl[2][1]=0.0;
      Gvvcoupl[3][0]=1.0;  Gvvcoupl[3][1]=0.0; // 2h+
      Gvvcoupl[4][0]=0.0;  Gvvcoupl[4][1]=0.0;
      Gvvcoupl[5][0]=0.0;  Gvvcoupl[5][1]=0.0;
      Gvvcoupl[6][0]=0.0;  Gvvcoupl[6][1]=0.0;
      Gvvcoupl[7][0]=0.0;  Gvvcoupl[7][1]=0.0;
      Gvvcoupl[8][0]=0.0;  Gvvcoupl[8][1]=0.0;
      Gvvcoupl[9][0]=0.0;  Gvvcoupl[9][1]=0.0;

      isSpinTwo = true;
    }
    // 2b+
    else if (_process == TVar::H2_g5){
      // gg production coupling constants
      Gggcoupl[0][0]=1.0;  Gggcoupl[0][1]=0.0;  // 2b+
      Gggcoupl[1][0]=0.0;  Gggcoupl[1][1]=0.0;
      Gggcoupl[2][0]=0.0;  Gggcoupl[2][1]=0.0;
      Gggcoupl[3][0]=0.0;  Gggcoupl[3][1]=0.0;
      Gggcoupl[4][0]=0.0;  Gggcoupl[4][1]=0.0;

      // Graviton->ZZ coupling constants 
      Gvvcoupl[0][0]=0.0;  Gvvcoupl[0][1]=0.0;
      Gvvcoupl[1][0]=0.0;  Gvvcoupl[1][1]=0.0;
      Gvvcoupl[2][0]=0.0;  Gvvcoupl[2][1]=0.0;
      Gvvcoupl[3][0]=0.0;  Gvvcoupl[3][1]=0.0;
      Gvvcoupl[4][0]=1.0;  Gvvcoupl[4][1]=0.0; // 2b+
      Gvvcoupl[5][0]=0.0;  Gvvcoupl[5][1]=0.0;
      Gvvcoupl[6][0]=0.0;  Gvvcoupl[6][1]=0.0;
      Gvvcoupl[7][0]=0.0;  Gvvcoupl[7][1]=0.0;
      Gvvcoupl[8][0]=0.0;  Gvvcoupl[8][1]=0.0;
      Gvvcoupl[9][0]=0.0;  Gvvcoupl[9][1]=0.0;

      isSpinTwo = true;
    }
    else if (_process == TVar::H2_g2){
      // gg production coupling constants
      Gggcoupl[0][0]=0.0;  Gggcoupl[0][1]=0.0;
      Gggcoupl[1][0]=1.0;  Gggcoupl[1][1]=0.0;
      Gggcoupl[2][0]=0.0;  Gggcoupl[2][1]=0.0;
      Gggcoupl[3][0]=0.0;  Gggcoupl[3][1]=0.0;
      Gggcoupl[4][0]=0.0;  Gggcoupl[4][1]=0.0;

      // Graviton->ZZ coupling constants
      Gvvcoupl[0][0]=0.0;  Gvvcoupl[0][1]=0.0;
      Gvvcoupl[1][0]=1.0;  Gvvcoupl[1][1]=0.0;
      Gvvcoupl[2][0]=0.0;  Gvvcoupl[2][1]=0.0;
      Gvvcoupl[3][0]=0.0;  Gvvcoupl[3][1]=0.0;
      Gvvcoupl[4][0]=0.0;  Gvvcoupl[4][1]=0.0;
      Gvvcoupl[5][0]=0.0;  Gvvcoupl[5][1]=0.0;
      Gvvcoupl[6][0]=0.0;  Gvvcoupl[6][1]=0.0;
      Gvvcoupl[7][0]=0.0;  Gvvcoupl[7][1]=0.0;
      Gvvcoupl[8][0]=0.0;  Gvvcoupl[8][1]=0.0;
      Gvvcoupl[9][0]=0.0;  Gvvcoupl[9][1]=0.0;

      isSpinTwo = true;
    }
    // 2h3plus
    else if (_process == TVar::H2_g3){
      // gg production coupling constants
      Gggcoupl[0][0]=0.0;  Gggcoupl[0][1]=0.0;
      Gggcoupl[1][0]=0.0;  Gggcoupl[1][1]=0.0;
      Gggcoupl[2][0]=1.0;  Gggcoupl[2][1]=0.0;
      Gggcoupl[3][0]=0.0;  Gggcoupl[3][1]=0.0;
      Gggcoupl[4][0]=0.0;  Gggcoupl[4][1]=0.0;

      // Graviton->ZZ coupling constants
      Gvvcoupl[0][0]=0.0;  Gvvcoupl[0][1]=0.0;
      Gvvcoupl[1][0]=0.0;  Gvvcoupl[1][1]=0.0;
      Gvvcoupl[2][0]=1.0;  Gvvcoupl[2][1]=0.0;
      Gvvcoupl[3][0]=0.0;  Gvvcoupl[3][1]=0.0;
      Gvvcoupl[4][0]=0.0;  Gvvcoupl[4][1]=0.0;
      Gvvcoupl[5][0]=0.0;  Gvvcoupl[5][1]=0.0;
      Gvvcoupl[6][0]=0.0;  Gvvcoupl[6][1]=0.0;
      Gvvcoupl[7][0]=0.0;  Gvvcoupl[7][1]=0.0;
      Gvvcoupl[8][0]=0.0;  Gvvcoupl[8][1]=0.0;
      Gvvcoupl[9][0]=0.0;  Gvvcoupl[9][1]=0.0;

      isSpinTwo = true;
    }
    // 2h6+
    else if (_process == TVar::H2_g6){
      // gg production coupling constants
      Gggcoupl[0][0]=1.0;  Gggcoupl[0][1]=0.0;
      Gggcoupl[1][0]=0.0;  Gggcoupl[1][1]=0.0;
      Gggcoupl[2][0]=0.0;  Gggcoupl[2][1]=0.0;
      Gggcoupl[3][0]=0.0;  Gggcoupl[3][1]=0.0;
      Gggcoupl[4][0]=0.0;  Gggcoupl[4][1]=0.0;

      // Graviton->ZZ coupling constants
      Gvvcoupl[0][0]=0.0;  Gvvcoupl[0][1]=0.0;
      Gvvcoupl[1][0]=0.0;  Gvvcoupl[1][1]=0.0;
      Gvvcoupl[2][0]=0.0;  Gvvcoupl[2][1]=0.0;
      Gvvcoupl[3][0]=0.0;  Gvvcoupl[3][1]=0.0;
      Gvvcoupl[4][0]=0.0;  Gvvcoupl[4][1]=0.0;
      Gvvcoupl[5][0]=1.0;  Gvvcoupl[5][1]=0.0;
      Gvvcoupl[6][0]=0.0;  Gvvcoupl[6][1]=0.0;
      Gvvcoupl[7][0]=0.0;  Gvvcoupl[7][1]=0.0;
      Gvvcoupl[8][0]=0.0;  Gvvcoupl[8][1]=0.0;
      Gvvcoupl[9][0]=0.0;  Gvvcoupl[9][1]=0.0;

      isSpinTwo = true;
    }
    // 2h7plus
    else if (_process == TVar::H2_g7){
      // gg production coupling constants
      Gggcoupl[0][0]=1.0;  Gggcoupl[0][1]=0.0;
      Gggcoupl[1][0]=0.0;  Gggcoupl[1][1]=0.0;
      Gggcoupl[2][0]=0.0;  Gggcoupl[2][1]=0.0;
      Gggcoupl[3][0]=0.0;  Gggcoupl[3][1]=0.0;
      Gggcoupl[4][0]=0.0;  Gggcoupl[4][1]=0.0;

      // Graviton->ZZ coupling constants
      Gvvcoupl[0][0]=0.0;  Gvvcoupl[0][1]=0.0;
      Gvvcoupl[1][0]=0.0;  Gvvcoupl[1][1]=0.0;
      Gvvcoupl[2][0]=0.0;  Gvvcoupl[2][1]=0.0;
      Gvvcoupl[3][0]=0.0;  Gvvcoupl[3][1]=0.0;
      Gvvcoupl[4][0]=0.0;  Gvvcoupl[4][1]=0.0;
      Gvvcoupl[5][0]=0.0;  Gvvcoupl[5][1]=0.0;
      Gvvcoupl[6][0]=1.0;  Gvvcoupl[6][1]=0.0;
      Gvvcoupl[7][0]=0.0;  Gvvcoupl[7][1]=0.0;
      Gvvcoupl[8][0]=0.0;  Gvvcoupl[8][1]=0.0;
      Gvvcoupl[9][0]=0.0;  Gvvcoupl[9][1]=0.0;

      isSpinTwo = true;
    }
    // 2h9minus
    else if (_process == TVar::H2_g9){
      // gg production coupling constants
      Gggcoupl[0][0]=0.0;  Gggcoupl[0][1]=0.0;
      Gggcoupl[1][0]=0.0;  Gggcoupl[1][1]=0.0;
      Gggcoupl[2][0]=0.0;  Gggcoupl[2][1]=0.0;
      Gggcoupl[3][0]=0.0;  Gggcoupl[3][1]=0.0;
      Gggcoupl[4][0]=1.0;  Gggcoupl[4][1]=0.0;

      // Graviton->ZZ coupling constants
      Gvvcoupl[0][0]=0.0;  Gvvcoupl[0][1]=0.0;
      Gvvcoupl[1][0]=0.0;  Gvvcoupl[1][1]=0.0;
      Gvvcoupl[2][0]=0.0;  Gvvcoupl[2][1]=0.0;
      Gvvcoupl[3][0]=0.0;  Gvvcoupl[3][1]=0.0;
      Gvvcoupl[4][0]=0.0;  Gvvcoupl[4][1]=0.0;
      Gvvcoupl[5][0]=0.0;  Gvvcoupl[5][1]=0.0;
      Gvvcoupl[6][0]=0.0;  Gvvcoupl[6][1]=0.0;
      Gvvcoupl[7][0]=0.0;  Gvvcoupl[7][1]=0.0;
      Gvvcoupl[8][0]=1.0;  Gvvcoupl[8][1]=0.0;
      Gvvcoupl[9][0]=0.0;  Gvvcoupl[9][1]=0.0;

      isSpinTwo = true;
    }
    // 2h10minus
    else if (_process == TVar::H2_g10){
      // gg production coupling constants
      Gggcoupl[0][0]=0.0;  Gggcoupl[0][1]=0.0;
      Gggcoupl[1][0]=0.0;  Gggcoupl[1][1]=0.0;
      Gggcoupl[2][0]=0.0;  Gggcoupl[2][1]=0.0;
      Gggcoupl[3][0]=0.0;  Gggcoupl[3][1]=0.0;
      Gggcoupl[4][0]=1.0;  Gggcoupl[4][1]=0.0;

      // Graviton->ZZ coupling constants
      Gvvcoupl[0][0]=0.0;  Gvvcoupl[0][1]=0.0;
      Gvvcoupl[1][0]=0.0;  Gvvcoupl[1][1]=0.0;
      Gvvcoupl[2][0]=0.0;  Gvvcoupl[2][1]=0.0;
      Gvvcoupl[3][0]=0.0;  Gvvcoupl[3][1]=0.0;
      Gvvcoupl[4][0]=0.0;  Gvvcoupl[4][1]=0.0;
      Gvvcoupl[5][0]=0.0;  Gvvcoupl[5][1]=0.0;
      Gvvcoupl[6][0]=0.0;  Gvvcoupl[6][1]=0.0;
      Gvvcoupl[7][0]=0.0;  Gvvcoupl[7][1]=0.0;
      Gvvcoupl[8][0]=0.0;  Gvvcoupl[8][1]=0.0;
      Gvvcoupl[9][0]=1.0;  Gvvcoupl[9][1]=0.0;

      isSpinTwo = true;
    }
    else if (_process == TVar::SelfDefine_spin2){
      for (int j=0; j<2; j++){
        for (int i=0; i<SIZE_GGG; i++) Gggcoupl[i][j] = (selfDSpinTwoCoupl.Gggcoupl)[i][j];
        for (int i=0; i<SIZE_GQQ; i++) Gqqcoupl[i][j] = (selfDSpinTwoCoupl.Gqqcoupl)[i][j];
        for (int i=0; i<SIZE_GVV; i++) Gvvcoupl[i][j] = (selfDSpinTwoCoupl.Gvvcoupl)[i][j];
      }
      isSpinTwo = true;
    }

    // 1+
    else if (_process == TVar::H1plus) {
      // z->vv coupling constants
      Zvvcoupl[0][0]=0.0;  Zvvcoupl[0][1]=0.0;
      Zvvcoupl[1][0]=1.0;  Zvvcoupl[1][1]=0.0; // 1+
      isSpinOne = true;
    }
    // 1-
    else if (_process == TVar::H1minus) isSpinOne = true; // Already the default
    else if (_process == TVar::SelfDefine_spin1){
      for (int j=0; j<2; j++){
        for (int i=0; i<SIZE_ZQQ; i++) Zqqcoupl[i][j] = (selfDSpinOneCoupl.Zqqcoupl)[i][j];
        for (int i=0; i<SIZE_ZVV; i++) Zvvcoupl[i][j] = (selfDSpinOneCoupl.Zvvcoupl)[i][j];
      }
      isSpinOne = true;
    }

    if (isSpinZero){
      SetJHUGenSpinZeroGGCouplings(Hggcoupl);
      SetJHUGenSpinZeroVVCouplings(Hvvcoupl, HvvCLambda_qsq, HvvLambda_qsq, false);
    }
    else if (isSpinOne) SetJHUGenSpinOneCouplings(Zqqcoupl, Zvvcoupl);
    else if (isSpinTwo) SetJHUGenSpinTwoCouplings(Gggcoupl, Gvvcoupl, Gqqcoupl);
    else{
      cerr << "TEvtProb::XsecCalc: JHUGen ME is not spin zero, one or two! The process is described by Process: " << _process << ", Production: " << _production << ", and ME: " << _matrixElement << endl;
      msqjk = 0;
    }
    if (isSpinZero || isSpinOne || isSpinTwo) msqjk = JHUGenMatEl(_process, _production, &mcfm_event);
  } // end of JHUGen matrix element calculations

  if (msqjk<=0){ mcfm_event.pswt=0; }

  flux=fbGeV2/(mcfm_event.p[0].Energy()*mcfm_event.p[1].Energy())	/(4*W);
  //dXsec=msqjk*flux;
  dXsec=msqjk;


  if (verbosity >= TVar::DEBUG){
    cout <<"Process " << TVar::ProcessName(_process)
      << " TEvtProb::XsecCalc(): dXsec=" << dXsec
      << " Msq=" << msqjk
      << " flux=" << flux
      << endl;
  }

  ResetCouplings(); // Should come first
  if (forceUseMCFM){ // Set defaults. Should come next...
    if (needBSMHiggs) SetLeptonInterf(TVar::DefaultLeptonInterf);
    SetMCFMSpinZeroVVCouplings(false, (selfDSpinZeroCoupl.Hzzcoupl), (selfDSpinZeroCoupl.Hzzcoupl)); // ... because of this!
  }
  ResetRenFacScaleMode();
  return dXsec;
}

double TEvtProb::XsecCalc_VVXVV(
  TVar::Process proc, TVar::Production production,
  const hzz4l_event_type &hzz4l_event,
  TVar::VerbosityLevel verbosity
  ){
  ResetIORecord();

  //Initialize Process
  SetProcess(proc);
  SetProduction(production);
  bool forceUseMCFM = (_matrixElement == TVar::MCFM || _process == TVar::bkgZZ_SMHiggs);

  int flavor = abs(hzz4l_event.PdgCode[0]) == abs(hzz4l_event.PdgCode[2]) ? 1 : 3;
  bool needBSMHiggs=false;
  // _process == TVar::bkgZZ_SMHiggs && _matrixElement == TVar::JHUGen is still MCFM
  if (forceUseMCFM){ // Always uses MCFM
    for (int vv = 0; vv < SIZE_HVV; vv++){
      if (
        (selfDSpinZeroCoupl.Hzzcoupl)[vv][1] != 0 || (selfDSpinZeroCoupl.Hzzcoupl)[vv][0] != 0 ||
        (selfDSpinZeroCoupl.Hwwcoupl)[vv][1] != 0 || (selfDSpinZeroCoupl.Hwwcoupl)[vv][0] != 0
        ){ needBSMHiggs = true; break; } // Only possible if selfDefine is called.
    }
    if (needBSMHiggs) SetLeptonInterf(TVar::InterfOn); // All anomalous coupling computations have to use lepton interference
    SetMCFMSpinZeroVVCouplings(needBSMHiggs, (selfDSpinZeroCoupl.Hzzcoupl), (selfDSpinZeroCoupl.Hwwcoupl));
    My_choose(_process, _production, _leptonInterf, flavor);
  }

  //constants
  double sqrts = 2.*EBEAM;
  double W=sqrts*sqrts;

  //Weight calculation
  double flux=1.;
  double dXsec=0.;

  mcfm_event_type mcfm_event;
  // assign the right initial momentum
  // assumes the events are boosted to have 0 transverse momenta
  TLorentzVector totalMom = hzz4l_event.p[0] + hzz4l_event.p[1] + hzz4l_event.p[2] + hzz4l_event.p[3] + hzz4l_event.extraParticle_p[0] + hzz4l_event.extraParticle_p[1];
  double sysPz= totalMom.Z();
  double sysE = totalMom.T();
  mcfm_event.p[2].SetPxPyPzE(hzz4l_event.p[0].Px(), hzz4l_event.p[0].Py(), hzz4l_event.p[0].Pz(), hzz4l_event.p[0].Energy());
  mcfm_event.p[3].SetPxPyPzE(hzz4l_event.p[1].Px(), hzz4l_event.p[1].Py(), hzz4l_event.p[1].Pz(), hzz4l_event.p[1].Energy());
  mcfm_event.p[4].SetPxPyPzE(hzz4l_event.p[2].Px(), hzz4l_event.p[2].Py(), hzz4l_event.p[2].Pz(), hzz4l_event.p[2].Energy());
  mcfm_event.p[5].SetPxPyPzE(hzz4l_event.p[3].Px(), hzz4l_event.p[3].Py(), hzz4l_event.p[3].Pz(), hzz4l_event.p[3].Energy());
  mcfm_event.p[6].SetPxPyPzE(hzz4l_event.extraParticle_p[0].Px(), hzz4l_event.extraParticle_p[0].Py(), hzz4l_event.extraParticle_p[0].Pz(), hzz4l_event.extraParticle_p[0].Energy());
  mcfm_event.p[7].SetPxPyPzE(hzz4l_event.extraParticle_p[1].Px(), hzz4l_event.extraParticle_p[1].Py(), hzz4l_event.extraParticle_p[1].Pz(), hzz4l_event.extraParticle_p[1].Energy());

  mcfm_event.PdgCode[0] = 1;
  mcfm_event.PdgCode[1] = -1;
  mcfm_event.PdgCode[2] = hzz4l_event.PdgCode[0];
  mcfm_event.PdgCode[3] = hzz4l_event.PdgCode[1];
  mcfm_event.PdgCode[4] = hzz4l_event.PdgCode[2];
  mcfm_event.PdgCode[5] = hzz4l_event.PdgCode[3];
  mcfm_event.PdgCode[6] = hzz4l_event.extraParticle_PdgCode[0];
  mcfm_event.PdgCode[7] = hzz4l_event.extraParticle_PdgCode[1];

  //Matrix Element evaluation in qX=qY=0 frame
  //Evaluate f(x1)f(x2)|M(q)|/x1/x2 
  // 
  double qX = totalMom.X();
  double qY = totalMom.Y();
  double qE = sysE;

  if ((qX*qX+qY*qY)>0){
    TVector3 boostV(qX/qE, qY/qE, 0);
    for (int ipt=2; ipt<8; ipt++) mcfm_event.p[ipt].Boost(-boostV);
    totalMom.Boost(-boostV);
  }

  sysPz= totalMom.Z();
  sysE = totalMom.T();
  double pz0 = (sysE+sysPz)/2.;
  double pz1 = -(sysE-sysPz)/2.;
  mcfm_event.p[0].SetPxPyPzE(0., 0., pz0, TMath::Abs(pz0));
  mcfm_event.p[1].SetPxPyPzE(0., 0., pz1, TMath::Abs(pz1));

  //event selections in Lab Frame
  double msqjk=0;
  if (forceUseMCFM) msqjk = SumMatrixElementPDF(_process, _production, _matrixElement, &event_scales, &mcfm_event, &RcdME, &flux, EBEAM, (selfDSpinZeroCoupl.Hvvcoupl_freenorm));
  if (msqjk<=0){ mcfm_event.pswt=0; }
  flux=fbGeV2/(mcfm_event.p[0].Energy()*mcfm_event.p[1].Energy())	/(4*W);
  //dXsec=msqjk*flux;
  dXsec=msqjk;

  if (verbosity >= TVar::DEBUG){
    cout <<"Process " << TVar::ProcessName(_process)
      << " TEvtProb::XsecCalc(): dXsec=" << dXsec
      << " Msq=" << msqjk
      << " flux=" << flux
      << endl;
  }

  ResetCouplings(); // Should come first
  if (forceUseMCFM){ // Set defaults. Should come next...
    if (needBSMHiggs) SetLeptonInterf(TVar::DefaultLeptonInterf);
    SetMCFMSpinZeroVVCouplings(false, (selfDSpinZeroCoupl.Hzzcoupl), (selfDSpinZeroCoupl.Hzzcoupl)); // ... because of this!
  }
  ResetRenFacScaleMode();
  return dXsec;
}

// Cross-section calculations for H + 2 jets
double TEvtProb::XsecCalcXJJ(TVar::Process proc, TVar::Production production, TLorentzVector p4[3], TVar::VerbosityLevel verbosity){
  ResetIORecord();
  double dXsec = 0;

  // Initialize Process
  SetProcess(proc);
  SetProduction(production);
  //constants
  //double sqrts = 2.*EBEAM;
  //double W=sqrts*sqrts;

  // first/second number is the real/imaginary part  
  double Hggcoupl[SIZE_HGG][2] ={ { 0 } };
  double Hzzcoupl[SIZE_HVV_VBF][2] ={ { 0 } };
  double Hwwcoupl[SIZE_HVV_VBF][2] ={ { 0 } };
  double HzzLambda_qsq[4][3] ={ { 0 } };
  double HwwLambda_qsq[4][3] ={ { 0 } };
  int HzzCLambda_qsq[3] ={ 0 };
  int HwwCLambda_qsq[3] ={ 0 };

  Hggcoupl[0][0]=1.0;  Hggcoupl[0][1]=0.0; // g2 
  Hzzcoupl[0][0]=1.0;  Hzzcoupl[0][1]=0.0; // g1
  Hwwcoupl[0][0]=1.0;  Hwwcoupl[0][1]=0.0; // g1
  for (int ic=0; ic<4; ic++){
    for (int ik=0; ik<3; ik++){
      HzzLambda_qsq[ic][ik]=100.;
      HwwLambda_qsq[ic][ik]=100.;
    }
  }

  // 0-
  if (_process == TVar::H0minus){
    Hggcoupl[0][0] = 0.0;
    Hggcoupl[2][0] = 1.0;

    Hzzcoupl[0][0] = 0.0;
    Hzzcoupl[3][0] = 1.0;
    Hwwcoupl[0][0] = 0.0;
    Hwwcoupl[3][0] = 1.0;
  }
  // 0+h
  else if (_process == TVar::H0hplus) { // No need to re-set ggcoupl
    Hzzcoupl[0][0] = 0.0;
    Hzzcoupl[1][0] = 1.0;
    Hwwcoupl[0][0] = 0.0;
    Hwwcoupl[1][0] = 1.0;
  }
  // 0+L1
  else if (_process == TVar::H0_g1prime2){ // No need to re-set ggcoupl
    Hzzcoupl[0][0] = 0.;
    Hzzcoupl[5][0] = -12046.01;
    Hwwcoupl[0][0] = 0.;
    Hwwcoupl[5][0] = -12046.01;
  }
  else if (_process == TVar::SelfDefine_spin0){
    for (int j=0; j<2; j++){
      for (int i=0; i<SIZE_HGG; i++) Hggcoupl[i][j] = (selfDSpinZeroCoupl.Hggcoupl)[i][j];
      for (int i=0; i<SIZE_HVV_VBF; i++){
        Hzzcoupl[i][j] = (selfDSpinZeroCoupl.Hzzcoupl_NoGamma)[i][j];
        Hwwcoupl[i][j] = (selfDSpinZeroCoupl.Hwwcoupl_NoGamma)[i][j];
      }
    }
    for (int j=0; j<3; j++){
      for (int i=0; i<4; i++){
        HzzLambda_qsq[i][j] = (selfDSpinZeroCoupl.HzzLambda_qsq)[i][j];
        HwwLambda_qsq[i][j] = (selfDSpinZeroCoupl.HwwLambda_qsq)[i][j];
      }
      HzzCLambda_qsq[j] = (selfDSpinZeroCoupl.HzzCLambda_qsq)[j];
      HwwCLambda_qsq[j] = (selfDSpinZeroCoupl.HwwCLambda_qsq)[j];
    }
  }
  SetJHUGenSpinZeroGGCouplings(Hggcoupl);
  SetJHUGenSpinZeroVVCouplings_NoGamma(Hzzcoupl, HzzCLambda_qsq, HzzLambda_qsq, false);
  SetJHUGenSpinZeroVVCouplings_NoGamma(Hwwcoupl, HwwCLambda_qsq, HwwLambda_qsq, true);

  // input kinematics 
  //  !----- p1 and p2 used to get hadronic s
  //  !----- P(p1)+P(p2) -> j(p3) + j(p4) + H(p5)
  // p[0] -> p1
  // p[1] -> p2
  // p[2] -> p3
  // p[3] -> p4
  // p[4] -> p5
  TLorentzVector p[5];
  p[2].SetPxPyPzE(p4[0].Px(), p4[0].Py(), p4[0].Pz(), p4[0].E());
  p[3].SetPxPyPzE(p4[1].Px(), p4[1].Py(), p4[1].Pz(), p4[1].E());
  p[4].SetPxPyPzE(p4[2].Px(), p4[2].Py(), p4[2].Pz(), p4[2].E());

  TLorentzVector pCoM = p[2] + p[3] + p[4];
  double qX = pCoM.X();
  double qY = pCoM.Y();
  double qE = pCoM.E();
  double qPt = (qX*qX+qY*qY);
  if ((qPt)>0){
    TVector3 boostV(-qX/qE, -qY/qE, 0); // Unit boost vector
    for (int ipt=2; ipt<5; ipt++) p[ipt].Boost(boostV);
  }

  // assign the right initial momentum
  // assumes the events are boosted to have 0 transverse momenta
  double sysPz= (pCoM).Pz();
  double sysE = (pCoM).Energy();
  double pz0 = (sysE+sysPz)/2.;
  double pz1 = -(sysE-sysPz)/2.;
  p[0].SetPxPyPzE(0., 0., pz0, TMath::Abs(pz0));
  p[1].SetPxPyPzE(0., 0., pz1, TMath::Abs(pz1));


  // calculate the matrix element squared
  if (_matrixElement == TVar::JHUGen){
    dXsec = HJJMatEl(_process, _production, _matrixElement, &event_scales, &RcdME, p, verbosity, EBEAM);
    if (verbosity >= TVar::DEBUG) std::cout <<"Process " << TVar::ProcessName(_process) << " TEvtProb::XsecCalc(): dXsec=" << dXsec << endl;
  }
  else std::cout << "Non-JHUGen vbfMELA is not supported!" << std::endl;

  ResetCouplings();
  ResetRenFacScaleMode();
  return dXsec;
}

// Cross-section calculations for H (SM) + 1 jet
double TEvtProb::XsecCalcXJ(TVar::Process proc, TVar::Production production, TLorentzVector p4[2], TVar::VerbosityLevel verbosity){
  ResetIORecord();
  double dXsec = 0;

  // Initialize Process
  SetProcess(proc);
  SetProduction(production);
  //constants
  //double sqrts = 2.*EBEAM;
  //double W=sqrts*sqrts;

  // first/second number is the real/imaginary part  


  // input kinematics 
  //  !----- p1 and p2 used to get hadronic s
  //  !----- P(p1)+P(p2) -> H(p3) + j(p4)
  // p[0] -> p1
  // p[1] -> p2
  // p[2] -> p3
  // p[3] -> p4
  // p[4] -> 0
  TLorentzVector p[5];
  for (int mom = 0; mom < 2; mom++){
    p[mom+2].SetPxPyPzE(p4[mom].Px(), p4[mom].Py(), p4[mom].Pz(), p4[mom].E());
  }
  p[4].SetPxPyPzE(0, 0, 0, 0);

  TLorentzVector pCoM = p[2] + p[3];
  double qX = pCoM.X();
  double qY = pCoM.Y();
  double qE = pCoM.E();
  double qPt = (qX*qX+qY*qY);
  if ((qPt)>0){
    TVector3 boostV(-qX/qE, -qY/qE, 0); // Unit boost vector
    for (int ipt=2; ipt<4; ipt++) p[ipt].Boost(boostV);
  }

  // assign the right initial momentum
  // assumes the events are boosted to have 0 transverse momenta
  double sysPz = (pCoM).Pz();
  double sysE = (pCoM).Energy();
  double pz0 = (sysE+sysPz)/2.;
  double pz1 = -(sysE-sysPz)/2.;
  p[0].SetPxPyPzE(0., 0., pz0, TMath::Abs(pz0));
  p[1].SetPxPyPzE(0., 0., pz1, TMath::Abs(pz1));


  // calculate the matrix element squared
  if (_matrixElement == TVar::JHUGen){
    dXsec = HJJMatEl(_process, _production, _matrixElement, &event_scales, &RcdME, p, verbosity, EBEAM);
    if (verbosity >= TVar::DEBUG) std::cout << "Process " << TVar::ProcessName(_process) << " TEvtProb::XsecCalc(): dXsec=" << dXsec << endl;
  }
  else std::cout << "Non-JHUGen HJ is not supported!" << std::endl;

  ResetRenFacScaleMode();
  return dXsec;
}


double TEvtProb::XsecCalc_VX(
  TVar::Process proc, TVar::Production production, vh_event_type &vh_event,
  TVar::VerbosityLevel verbosity
  ){
  ResetIORecord();
  MelaIO RcdME_temp;

  //Initialize Process
  SetProcess(proc);
  SetProduction(production);
  AllowSeparateWWCouplings(false); // HZZ couplings are used for both in spin-0

  // 0: Higgs; 1,2: V-daughters
  // PDG ID for V daughters could be passed as 0. While mothers cannot really be gluons, specifying 0 will mean averaging over u,d,c,s,b.
  bool inclusiveHadronicJets = (vh_event.PdgCode[1]==0 || vh_event.PdgCode[2]==0);
  const double N_Q=5.;

  //Weight calculation
  double msqjk=0;

  if (
    (abs(vh_event.PdgCode[2]) == 12 ||
    abs(vh_event.PdgCode[2]) == 14 ||
    abs(vh_event.PdgCode[2]) == 16) && _production == TVar::WH
    ){ // First daughter of W has to be a neutrino

    vh_event.p[1] = vh_event.p[1] + vh_event.p[2];
    vh_event.p[2] = vh_event.p[1] - vh_event.p[2];
    vh_event.p[1] = vh_event.p[1] - vh_event.p[2];
    vh_event.PdgCode[1] = vh_event.PdgCode[1] + vh_event.PdgCode[2];
    vh_event.PdgCode[2] = vh_event.PdgCode[1] - vh_event.PdgCode[2];
    vh_event.PdgCode[1] = vh_event.PdgCode[1] - vh_event.PdgCode[2];
  }

  int Vdecay_id[6] ={
    vh_event.PdgCode[1], vh_event.PdgCode[2],
    vh_event.PdgCode_Hdecay[0],
    vh_event.PdgCode_Hdecay[1],
    vh_event.PdgCode_Hdecay[2],
    vh_event.PdgCode_Hdecay[3]
  };

  TLorentzVector pCoM = vh_event.p[0] + vh_event.p[1] + vh_event.p[2];
  double qX = pCoM.X();
  double qY = pCoM.Y();
  double qE = pCoM.E();
  double qPt = (qX*qX+qY*qY);
  if ((qPt)>0){
    TVector3 boostV(-qX/qE, -qY/qE, 0); // Unit boost vector
    for (int ipt=0; ipt<3; ipt++) vh_event.p[ipt].Boost(boostV);
    for (int ipt=0; ipt<4; ipt++) vh_event.pHdecay[ipt].Boost(boostV);
  }

  // assign the right initial momentum
  // assumes the events are boosted to have 0 transverse momenta
  double sysPz= (vh_event.p[0] + vh_event.p[1] + vh_event.p[2]).Pz();
  double sysE = (vh_event.p[0] + vh_event.p[1] + vh_event.p[2]).Energy();
  double pz0 = (sysE+sysPz)/2.;
  double pz1 = -(sysE-sysPz)/2.;
  TLorentzVector pVH[5];
  TLorentzVector pHdaughter[4];
  pVH[2].SetPxPyPzE(vh_event.p[0].Px(), vh_event.p[0].Py(), vh_event.p[0].Pz(), vh_event.p[0].Energy());
  pVH[3].SetPxPyPzE(vh_event.p[1].Px(), vh_event.p[1].Py(), vh_event.p[1].Pz(), vh_event.p[1].Energy());
  pVH[4].SetPxPyPzE(vh_event.p[2].Px(), vh_event.p[2].Py(), vh_event.p[2].Pz(), vh_event.p[2].Energy());
  pVH[0].SetPxPyPzE(0, 0, pz0, TMath::Abs(pz0));
  pVH[1].SetPxPyPzE(0, 0, pz1, TMath::Abs(pz1));
  for (int ipt=0; ipt<4; ipt++) pHdaughter[ipt]=vh_event.pHdecay[ipt];

  if (_matrixElement == TVar::JHUGen) {
    SetJHUGenDistinguishWWCouplings(false); // HZZ couplings are used for both in spin-0
    // Set Couplings at the HVV* vertex
    double Hvvcoupl[SIZE_HVV_VBF][2] ={ { 0 } };
    double HvvLambda_qsq[4][3] ={ { 0 } };
    int HvvCLambda_qsq[3] ={ 0 };

    // By default set the Spin 0 couplings for SM case
    Hvvcoupl[0][0]=1.0;  Hvvcoupl[0][1]=0.0;   // first/second number is the real/imaginary part
    for (int ic=0; ic<4; ic++){
      for (int ik=0; ik<3; ik++) HvvLambda_qsq[ic][ik]=100.;
    }

    // 0-
    if (_process == TVar::H0minus) {
      Hvvcoupl[0][0] = 0.0;
      Hvvcoupl[1][0] = 0.0;
      Hvvcoupl[2][0] = 0.0;
      Hvvcoupl[3][0] = 1.0;
    }
    // 0h+
    else if (_process == TVar::H0hplus) {
      Hvvcoupl[0][0] = 0.0;
      Hvvcoupl[1][0] = 1.0;
      Hvvcoupl[2][0] = 0.0;
      Hvvcoupl[3][0] = 0.0;
    }
    // 0+L1
    else if (_process == TVar::H0_g1prime2){
      Hvvcoupl[0][0] = 0.;
      Hvvcoupl[5][0] = -12046.01;
    }
    else if (_process == TVar::SelfDefine_spin0){
      for (int i=0; i<SIZE_HVV_VBF; i++){
        for (int j=0; j<2; j++) Hvvcoupl[i][j] = (selfDSpinZeroCoupl.Hzzcoupl_NoGamma)[i][j];
      }
      for (int j=0; j<3; j++){
        for (int i=0; i<4; i++) HvvLambda_qsq[i][j] = (selfDSpinZeroCoupl.HzzLambda_qsq)[i][j];
        HvvCLambda_qsq[j] = (selfDSpinZeroCoupl.HzzCLambda_qsq)[j];
      }
    }
    SetJHUGenSpinZeroVVCouplings_NoGamma(Hvvcoupl, HvvCLambda_qsq, HvvLambda_qsq, false);

    if (!inclusiveHadronicJets) msqjk = VHiggsMatEl(_process, _production, &event_scales, &RcdME, pVH, pHdaughter, Vdecay_id, verbosity, EBEAM);
    else{
      for (int outgoing1 = -nf; outgoing1 <= nf; outgoing1++){
        if (_production == TVar::ZH){
          if (outgoing1 <= 0) continue;
          Vdecay_id[0] = outgoing1;
          Vdecay_id[1] = -outgoing1;
          msqjk += (VHiggsMatEl(_process, _production, &event_scales, &RcdME_temp, pVH, pHdaughter, Vdecay_id, verbosity, EBEAM)) / N_Q; // Average over quark flavors
          RcdME.addMERecord(&RcdME_temp, 1./N_Q);
          RcdME_temp.reset();
        }
        else if (_production == TVar::WH){
          if (outgoing1 == 0) continue;
          if (outgoing1 == 2 || outgoing1 == 4){ // u or c to d-bar, b-bar or s-bar
            for (int outgoing2 = -nf; outgoing2 < 0; outgoing2++){
              if (abs(outgoing2) == abs(outgoing1)) continue;
              Vdecay_id[0] = outgoing1;
              Vdecay_id[1] = outgoing2;
              msqjk += (VHiggsMatEl(_process, _production, &event_scales, &RcdME_temp, pVH, pHdaughter, Vdecay_id, verbosity, EBEAM)) / 12.; // Average over quark flavors; CAUTION about 12: Depends on nf, (nf+1)*(nf-1)/2 or nf**2/2
              RcdME.addMERecord(&RcdME_temp, 1./N_Q);
              RcdME_temp.reset();
            }
          }
          if (outgoing1 == -2 || outgoing1 == -4){ // u-bar or c-bar to d, b or s
            for (int outgoing2 = 1; outgoing2 < nf + 1; outgoing2++){
              if (abs(outgoing2) == abs(outgoing1)) continue;
              Vdecay_id[0] = outgoing1;
              Vdecay_id[1] = outgoing2;
              msqjk += (VHiggsMatEl(_process, _production, &event_scales, &RcdME_temp, pVH, pHdaughter, Vdecay_id, verbosity, EBEAM)) / 12.; // Average over quark flavors; CAUTION about 12: Depends on nf
              RcdME.addMERecord(&RcdME_temp, 1./N_Q);
              RcdME_temp.reset();
            }
          }
        }
      }
    }
  } // end of JHUGen matrix element calculations
  else{
    std::cout << "Non-JHUGen VH is not supported!" << std::endl;
    msqjk = 0;
  }

  ResetCouplings();
  ResetRenFacScaleMode();
  return msqjk;
}

// Cross-section calculations for ttbar -> H
double TEvtProb::XsecCalc_TTX(
  TVar::Process proc, TVar::Production production,
  tth_event_type &tth_event,
  int topDecay, int topProcess,
  TVar::VerbosityLevel verbosity
  ){
  // Set Couplings at the TTH vertex
  double Hqqcoupl[SIZE_HQQ][2]={ { 0 } };

  //Initialize Process
  SetProcess(proc);
  SetProduction(production);

  double msq=0;

  TLorentzVector pCoM(0, 0, 0, 0);
  const int nV = 7;
  for (int vv=0; vv<nV; vv++) pCoM = pCoM + tth_event.p[vv];
  double qX = pCoM.X();
  double qY = pCoM.Y();
  double qE = pCoM.E();
  double qPt = (qX*qX+qY*qY);
  if ((qPt)>0){
    TVector3 boostV(-qX/qE, -qY/qE, 0); // Unit boost vector
    for (int vv=0; vv<nV; vv++) tth_event.p[vv].Boost(boostV);
  }

  TLorentzVector p_tbar(0, 0, 0, 0);
  TLorentzVector p_t(0, 0, 0, 0);

  bool unknownTopFlavor=false;
  int indexTTBAR=0;
  if (
    (tth_event.PdgCode_tdecay[0][0]==0 || tth_event.PdgCode_tdecay[1][0]==0) ||
    (tth_event.PdgCode_tdecay[0][0]>0 && tth_event.PdgCode_tdecay[1][0]>0) ||
    (tth_event.PdgCode_tdecay[0][0]<0 && tth_event.PdgCode_tdecay[1][0]<0)
    ) unknownTopFlavor=true;
  else if (tth_event.PdgCode_tdecay[0][0]>0) indexTTBAR=1;

  for (int vv=1; vv<4; vv++) p_tbar = p_tbar + tth_event.p[vv + 3*indexTTBAR];
  for (int vv=1; vv<4; vv++) p_t = p_t + tth_event.p[vv + 3*(1-indexTTBAR)];

  TLorentzVector pTTH[11];
  pTTH[0].SetXYZT(0, 0, EBEAM, EBEAM);
  pTTH[1].SetXYZT(0, 0, -EBEAM, EBEAM);
  pTTH[2].SetXYZT(tth_event.p[0].X(), tth_event.p[0].Y(), tth_event.p[0].Z(), tth_event.p[0].T());
  pTTH[3].SetXYZT(p_tbar.X(), p_tbar.Y(), p_tbar.Z(), p_tbar.T());
  pTTH[4].SetXYZT(p_t.X(), p_t.Y(), p_t.Z(), p_t.T());
  for (int vv=1; vv<4; vv++) pTTH[vv+4].SetXYZT(tth_event.p[vv + 3*indexTTBAR].X(), tth_event.p[vv + 3*indexTTBAR].Y(), tth_event.p[vv + 3*indexTTBAR].Z(), tth_event.p[vv + 3*indexTTBAR].T());
  for (int vv=1; vv<4; vv++) pTTH[vv+7].SetXYZT(tth_event.p[vv + 3*(1-indexTTBAR)].X(), tth_event.p[vv + 3*(1-indexTTBAR)].Y(), tth_event.p[vv + 3*(1-indexTTBAR)].Z(), tth_event.p[vv + 3*(1-indexTTBAR)].T());


  if (_matrixElement == TVar::JHUGen){
    // By default set the Spin 0 couplings for SM case
    Hqqcoupl[0][0]=1.0;  Hqqcoupl[0][1]=0.0;   // first/second number is the real/imaginary part
    for (int i = 1; i<SIZE_HQQ; i++){ for (int com=0; com<2; com++) Hqqcoupl[i][com] = 0; }

    // 0-
    if (_process == TVar::H0minus) {
      Hqqcoupl[0][0] = 0.0;
      Hqqcoupl[1][0] = 1.0;
    }
    else if (_process == TVar::SelfDefine_spin0){
      for (int i=0; i<SIZE_HQQ; i++){
        for (int j=0; j<2; j++) Hqqcoupl[i][j] = (selfDSpinZeroCoupl.Hqqcoupl)[i][j];
      }
    }
    SetJHUGenSpinZeroQQCouplings(Hqqcoupl);

    double topMass = 173.2;
    double topWidth = 2.;
    if (_production == TVar::bbH){
      topMass = 4.75;
      topWidth = 0;
    }
    msq = TTHiggsMatEl(_production, pTTH, _hmass, _hwidth, topMass, topWidth, topDecay, topProcess, verbosity);
    if (unknownTopFlavor){
      //cout << "Unknown top flavor" << endl;
      pTTH[4].SetXYZT(p_tbar.X(), p_tbar.Y(), p_tbar.Z(), p_tbar.T());
      pTTH[3].SetXYZT(p_t.X(), p_t.Y(), p_t.Z(), p_t.T());
      for (int vv=1; vv<4; vv++) pTTH[vv+7].SetXYZT(tth_event.p[vv + 3*indexTTBAR].X(), tth_event.p[vv + 3*indexTTBAR].Y(), tth_event.p[vv + 3*indexTTBAR].Z(), tth_event.p[vv + 3*indexTTBAR].T());
      for (int vv=1; vv<4; vv++) pTTH[vv+4].SetXYZT(tth_event.p[vv + 3*(1-indexTTBAR)].X(), tth_event.p[vv + 3*(1-indexTTBAR)].Y(), tth_event.p[vv + 3*(1-indexTTBAR)].Z(), tth_event.p[vv + 3*(1-indexTTBAR)].T());
      msq += TTHiggsMatEl(_production, pTTH, _hmass, _hwidth, topMass, topWidth, topDecay, topProcess, verbosity);
      msq *= 0.5;
    }
  }
  else std::cout << "Non-JHUGen ttH is not supported!" << std::endl;

  ResetCouplings();
  ResetRenFacScaleMode();
  return msq;
}


void TEvtProb::SetRenFacScaleMode(TVar::EventScaleScheme renormalizationSch, TVar::EventScaleScheme factorizationSch, double ren_sf, double fac_sf){
  event_scales.renomalizationScheme = renormalizationSch;
  event_scales.factorizationScheme = factorizationSch;
  event_scales.ren_scale_factor = ren_sf;
  event_scales.fac_scale_factor = fac_sf;
}


// this appears to be some kind of 
// way of setting MCFM parameters through
// an interface defined in TMCFM.hh
void TEvtProb::SetHiggsMass(double mass, float wHiggs){
  _hmass = mass;
  masses_mcfm_.hmass = _hmass;
  if (wHiggs < 0.){
    _hwidth = myCSW_->HiggsWidth(0, min(mass, 1000.));
    masses_mcfm_.hwidth = _hwidth;
  }
  else{
    _hwidth = wHiggs;
    masses_mcfm_.hwidth = _hwidth;
  }
  SetJHUGenHiggsMassWidth(_hmass, _hwidth);
}
