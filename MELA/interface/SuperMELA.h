#ifndef MELA_SuperMela_h
#define MELA_SuperMela_h


#include <Riostream.h>
#include <string>
#include <fstream>

#include "TLorentzVector.h"
#include "ZZMatrixElement/MELA/interface/Mela.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZZ4LRooPdfs.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZZ2L2QRooPdfs.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooCBShape.h"
#include "RooFFTConvPdf.h"



class SuperMELA {

 public:
  SuperMELA(double mH=125,string channel="4mu",int LHCsqrts=8);
  ~SuperMELA();
  void init();

  double GetSigShapeSystematic(string parName);
  double GetSigShapeParameter(string parName);

  void SetVerbosity(bool verb=true){verbose_=verb;}
  void SetDecayChannel(string myChan);
  void SetMH(double myMH){
    mHVal_=myMH;
    mH_rrv_->setVal(mHVal_);
    if(verbose_)std::cout<<"Setting MH to "<<mHVal_<<std::endl;
    init();
  }

  void SetPathToCards(string dirToCards){ pathToCards_=dirToCards;
    if(verbose_)std::cout<<"New path to cards is "<<pathToCards_.c_str()<<std::endl;}

  std::pair<double,double>   M4lProb(double m4l);
  std::pair<double,double>   M4lProb(std::pair<double, double>);

 private:

  void readSigParsFromFile(string &str_mean_CB,string &str_sigma_CB ,string &str_n_CB ,string &str_alpha_CB ,string &str_n2_CB ,string &str_alpha2_CB);
  void readBkgParsFromFile(std::vector<double> &apars );
  void readSigSystFromFile(string &str_mean_CB_err_e, string &str_mean_CB_err_m,
			   string &str_sigma_CB_err_e, string &str_sigma_CB_err_m);

  void calc_mZZ_range(const double mHVal,double &low_M,double &high_M );
  bool checkChannel();
  ///data members
  double mHVal_;
  double sqrts_;
  double lowMH_,highMH_;
  string strChan_; int ch_;
  bool verbose_;
  string pathToCards_;
  //signal m4l shape
  RooRealVar *m4l_rrv_;//this one is the variable!
  RooRealVar *mH_rrv_;//this one is a fixed param !
  RooRealVar *mean_dummy_,*sigma_dummy_,*alpha_dummy_,*n_dummy_;
  RooFormulaVar *n_CB_, *alpha_CB_, *n2_CB_, *alpha2_CB_,*mean_CB_,*sigma_CB_,*meanTOT_CB_;//,*gamma_BW_;
  RooFormulaVar *mean_CB_err_, *sigma_CB_err_;
  // RooCBShape *sig_CB_;
  RooDoubleCB *sig_CB_;
  RooRelBWUFParam *sig_BW_;//this is defined in HZZ4LRooPdfs.h
  RooFFTConvPdf *sig_FFT_;
  RooRealVar *mean_BW_,*width_BW_;
  double norm_sig_CB_, norm_sig_FFT_;

  //qqZZ background m4l shape
  RooRealVar *a0_qqZZ_, *a1_qqZZ_, *a2_qqZZ_, *a3_qqZZ_,
             *a4_qqZZ_, *a5_qqZZ_, *a6_qqZZ_, *a7_qqZZ_,
             *a8_qqZZ_, *a9_qqZZ_, *a10_qqZZ_,*a11_qqZZ_,
             *a12_qqZZ_,   *a13_qqZZ_  ;
  RooqqZZPdf_v2 *qqZZ_pdf_;       
  double norm_bkg_qqZZ_;


};


#endif
