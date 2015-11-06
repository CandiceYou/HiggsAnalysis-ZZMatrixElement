#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include "TMath.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include <ZZMatrixElement/MELA/interface/Mela.h>
#include <ZZMatrixElement/MELA/src/computeAngles.h>

using namespace RooFit;

using namespace std;

namespace {
  const bool doMuEffCorr    = true;
  const bool doEleEffCorr   = true;
  const bool doHighmassCorr = true; 
  const bool doHqTCorr      = false;
  const bool saveJets       = true;
  const bool applySystMu    = false;
  const bool applySystEle   = false;

  const int pwhg_flag = 0;    // 0 means standard high mass weights
                                // 1 means CPS+ + Interference
                                // 2 means CPS- + Interference
                                // 3 means CPS  + Interference+
                                // 4 means CPS  + Interference-

  const float Run2011AFraction = 0.465;
  const float Zmass = 91.1876;
  const float gamZ = 2.5;
  const float R1Val = 0.15;
  const float R2Val = 0.15;
  const int para = 1;
  const bool use_acc=false;
  const float M_muon = 0.105658389;
  const float M_electron = 0.00051099907;
  const float M_tau = 1.777;
  const int PDG_electron=11,PDG_muon=13,PDG_tau=15;
}



float getMCFMMELAWeight(Mela& myMela, int lepId[4], float angularOrdered[8], double ggcoupl[2],int useConstant=0){
	float myprob=1.0;
	int myflavor=-1;
	if(abs(lepId[0])==abs(lepId[1]) &&
		abs(lepId[0])==abs(lepId[2]) &&
		abs(lepId[0])==abs(lepId[3])){
			if(abs(lepId[0])==11) myflavor=1;
			else myflavor=2;
	}
	else myflavor=3;

	if(myflavor>=0) myMela.computeP(angularOrdered[0],angularOrdered[1],angularOrdered[2],angularOrdered[3],
		angularOrdered[4],angularOrdered[5],angularOrdered[6],angularOrdered[7],
	    myflavor,
	    myprob,
		useConstant
		);
	return myprob;
};
float getJHUGenMELAWeight(Mela& myMela, int lepId[4], float angularOrdered[8], double selfDHvvcoupl[SIZE_HVV][2]){
	float myprob=1.0;
	int myflavor=-1;
	if(abs(lepId[0])==abs(lepId[1]) &&
		abs(lepId[0])==abs(lepId[2]) &&
		abs(lepId[0])==abs(lepId[3])){
			if(abs(lepId[0])==11) myflavor=1;
			else myflavor=2;
	}
	else myflavor=3;

	if(myflavor>=0) myMela.computeP(angularOrdered[0],angularOrdered[1],angularOrdered[2],angularOrdered[3],
		angularOrdered[4],angularOrdered[5],angularOrdered[6],angularOrdered[7],
	    myflavor,
	    selfDHvvcoupl,
	    myprob
		);
	return myprob;
};
float getSuperMELA(Mela& myMela, int lepId[4], float mZZ, TVar::SuperMelaSyst syst){
	float myprob=1.0;
	TVar::LeptonFlavor myflavor=TVar::Flavor_Dummy;
	if(abs(lepId[0])==abs(lepId[1]) &&
		abs(lepId[0])==abs(lepId[2]) &&
		abs(lepId[0])==abs(lepId[3])){
			if(abs(lepId[0])==11) myflavor=TVar::Flavor_4e;
			else myflavor=TVar::Flavor_4mu;
	}
	else myflavor=TVar::Flavor_2e2mu;

	if(myflavor!=TVar::Flavor_Dummy) myMela.computePM4l(mZZ,
	    myflavor,
	    syst,
	    myprob
		);
	return myprob;
};


void testME_FullMELA(int erg_tev=8, float mPOLE=125.6){
	TVar::VerbosityLevel verbosity = TVar::INFO;

	Mela mela(erg_tev,mPOLE);
	TGraph* vaScale_4e = mela.vaScale_4e;
	TGraph* vaScale_4mu = mela.vaScale_4mu;
	TGraph* vaScale_2e2mu = mela.vaScale_2e2mu;
	TGraph* DggZZ_scalefactor = mela.DggZZ_scalefactor;


	double p0plus_VAJHU;
	double bkg_VAMCFM,bkg_STU_VAMCFM,bkg_S_VAMCFM,bkg_TU_VAMCFM;
	double ggzz_VAMCFM,ggZZ_prob_Total,ggHZZ_prob_Bare;
	double bkg_VAMCFM_noscale , ggzz_VAMCFM_noscale , ggHZZ_prob_pure , ggHZZ_prob_int , ggHZZ_prob_int_alt , ggHZZ_prob_pure_noscale , ggHZZ_prob_int_noscale;
	double qqScale , ggScale;
	double p0plus_m4l, p0plus_m4l_ScaleUp, p0plus_m4l_ScaleDown, p0plus_m4l_ResUp, p0plus_m4l_ResDown;
	double bkg_m4l, bkg_m4l_ScaleUp, bkg_m4l_ScaleDown, bkg_m4l_ResUp, bkg_m4l_ResDown;
/*
*    Row   *    ZZMass *    Z1Mass *    Z2Mass * helcosthe * helcosthe *    helphi * costhetas * phistarZ1 *
************************************************************************************************************
*   364271 * 285.97302 * 92.097068 * 93.300392 * -0.550634 * -0.428076 * -2.060295 * 0.3669659 * -0.000344 *
*/
/*
	float mzz = 285.97302; 
//	float mzz = 91.; 
	float m1 = 92.097068;
//	float m1 = 41.471450;
	float m2 = 93.300392;
	float h1 = -0.550634;
	float h2 = -0.428076;
	float phi = -2.060295;
	float hs = 0.3669659;
	float phi1 = -0.000344;
*/
	float mzz = 126.; 
//	float mzz = 91.; 
	float m1 = 91.471450;
//	float m1 = 41.471450;
	float m2 = 12.139782;
	float h1 = 0.2682896;
	float h2 = 0.1679779;
	float phi = 1.5969792;
	float hs = -0.727181;
	float phi1 = 1.8828257;
	
	int mflavor = 3;

//	int lepIdOrdered[4]={ 11,-11,13,-13 };
	int lepIdOrdered[4]={ 11,-11,11,-11 };
	float angularOrdered[8]={mzz,m1,m2,hs,h1,h2,phi,phi1};

	double selfDHvvcoupl[SIZE_HVV][2]={{0.}};
	double ggvvcoupl[2]={0,0};

	mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
	p0plus_m4l = getSuperMELA(mela,lepIdOrdered,mzz,TVar::SMSyst_None);

	mela.setProcess(TVar::bkgZZ, TVar::JHUGen, TVar::ZZGG);
	bkg_m4l = getSuperMELA(mela,lepIdOrdered,mzz,TVar::SMSyst_None);

	mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
	selfDHvvcoupl[0][0]=1;
	p0plus_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

	mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB);
	mela.setMelaLeptonInterference(TVar::InterfOn);
	bkg_VAMCFM = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |qqZZ|**2
	mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB_STU);
	mela.setMelaLeptonInterference(TVar::InterfOn);
	bkg_STU_VAMCFM = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |qqZZ|**2
	mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB_S);
	mela.setMelaLeptonInterference(TVar::InterfOn);
	bkg_S_VAMCFM = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |qqZZ|**2
	mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB_TU);
	mela.setMelaLeptonInterference(TVar::InterfOn);
	bkg_TU_VAMCFM = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |qqZZ|**2

	mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZGG);
	mela.setMelaLeptonInterference(TVar::InterfOn);
	ggzz_VAMCFM = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |ggZZ|**2

	mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
	mela.setMelaLeptonInterference(TVar::InterfOn);
	ggHZZ_prob_pure = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |ggHZZ|**2

	mela.setProcess(TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
	mela.setMelaLeptonInterference(TVar::InterfOn);
	ggZZ_prob_Total = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |ggZZ + ggHZZ|**2
	ggHZZ_prob_int = ggZZ_prob_Total - ggHZZ_prob_pure - ggzz_VAMCFM;

	ggScale=1.0;
	qqScale=1.0;
	if(abs(lepIdOrdered[0])==abs(lepIdOrdered[1]) &&
		abs(lepIdOrdered[0])==abs(lepIdOrdered[2]) &&
		abs(lepIdOrdered[0])==abs(lepIdOrdered[3])){
			if(abs(lepIdOrdered[0])==11){
				if(mzz > 900) qqScale = vaScale_4e->Eval(900.);
				else if (mzz <  100 ) qqScale = vaScale_4e->Eval(100.);
				else qqScale = vaScale_4e->Eval(mzz);

				if(mzz > 900) ggScale = vaScale_4e->Eval(900.);
				else if (mzz <  110 ) ggScale = vaScale_4e->Eval(110.);
				else ggScale = vaScale_4e->Eval(mzz);
			}
			else{
				if(mzz > 900) qqScale = vaScale_4mu->Eval(900.);
				else if (mzz <  100 ) qqScale = vaScale_4mu->Eval(100.);
				else qqScale = vaScale_4mu->Eval(mzz);

				if(mzz > 900) ggScale = vaScale_4mu->Eval(900.);
				else if (mzz <  110 ) ggScale = vaScale_4mu->Eval(110.);
				else ggScale = vaScale_4mu->Eval(mzz);
			};
	}
	else{
		if(mzz > 900) qqScale = vaScale_2e2mu->Eval(900.);
		else if (mzz <  100 ) qqScale = vaScale_2e2mu->Eval(100.);
		else qqScale = vaScale_2e2mu->Eval(mzz);

		if(mzz > 900) ggScale = vaScale_2e2mu->Eval(900.);
		else if (mzz <  110 ) ggScale = vaScale_2e2mu->Eval(110.);
		else ggScale = vaScale_2e2mu->Eval(mzz);
	};
	if(mzz > 900) ggScale /= DggZZ_scalefactor->Eval(900.);
	else if (mzz <  110 ) ggScale /= DggZZ_scalefactor->Eval(110.);
	else ggScale /= DggZZ_scalefactor->Eval(mzz);

//	ggZZ_prob_Total /= ggScale;
//	ggzz_VAMCFM /= ggScale;
//	ggHZZ_prob_pure /= ggScale;
//	ggHZZ_prob_int /= ggScale;
//	bkg_VAMCFM /= qqScale;

	ofstream tout("./testME_SingleEvent.txt");

	tout << "|ggZZ+ggHZZ|**2\t|ggZZ|**2\t|ggHZZ|**2\t|ggHZZ|**2 JHUGEN\tInterf.\t|qqZZ|**2\tPSig_m4l\tPBkg_m4l" << endl;
	tout << ggZZ_prob_Total << '\t'
		 << ggzz_VAMCFM << '\t'
		 << ggHZZ_prob_pure << '\t'
		 << p0plus_VAJHU << '\t'
		 << ggHZZ_prob_int << '\t'
		 << bkg_VAMCFM << '\t'
		 << p0plus_m4l << '\t'
		 << bkg_m4l << endl;

	tout << bkg_STU_VAMCFM << '\t'
		 << bkg_S_VAMCFM << '\t'
		 << bkg_TU_VAMCFM << endl;

	tout.close();
};

void testME_FullMELA_MultipleMELA(int erg_tev=8, float mPOLE=125.6){
  TVar::VerbosityLevel verbosity = TVar::INFO;

  ofstream tout("./testME_SingleEvent_MultipleMELA.txt");


  Mela mela(erg_tev, mPOLE);
  Mela mela2(erg_tev, mPOLE);


  double bkg_VAMCFM, ggzz_VAMCFM, ggZZ_prob_Total, ggHZZ_prob_pure;
  double bkg_VAMCFM2, ggzz_VAMCFM2, ggZZ_prob_Total2, ggHZZ_prob_pure2;
  double bkg_VAMCFM3, ggzz_VAMCFM3, ggZZ_prob_Total3, ggHZZ_prob_pure3;

  float mzz = 126.;
  float m1 = 91.471450;
  float m2 = 12.139782;
  float h1 = 0.2682896;
  float h2 = 0.1679779;
  float phi = 1.5969792;
  float hs = -0.727181;
  float phi1 = 1.8828257;

  int mflavor = 3;

  int lepIdOrdered[4]={ 11, -11, 11, -11 };
  float angularOrdered[8]={ mzz, m1, m2, hs, h1, h2, phi, phi1 };

  double selfDHvvcoupl[SIZE_HVV][2]={ { 0. } };
  double ggvvcoupl[2]={ 0, 0 };

  selfDHvvcoupl[0][0]=1;

  mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB);
  mela.setMelaLeptonInterference(TVar::InterfOn);
  bkg_VAMCFM = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |qqZZ|**2

  mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZGG);
  mela.setMelaLeptonInterference(TVar::InterfOn);
  ggzz_VAMCFM = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |ggZZ|**2

  mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
  mela.setMelaLeptonInterference(TVar::InterfOn);
  ggHZZ_prob_pure = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |ggHZZ|**2

  mela.setProcess(TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
  mela.setMelaLeptonInterference(TVar::InterfOn);
  ggZZ_prob_Total = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |ggZZ + ggHZZ|**2


  mela2.resetMCFM_EWKParameters(1.16639E-05, 1./128., 80.399, 91.1876, 0.23119);

  mela2.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB);
  mela2.setMelaLeptonInterference(TVar::InterfOn);
  bkg_VAMCFM2 = getMCFMMELAWeight(mela2, lepIdOrdered, angularOrdered, ggvvcoupl); // |qqZZ|**2

  mela2.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZGG);
  mela2.setMelaLeptonInterference(TVar::InterfOn);
  ggzz_VAMCFM2 = getMCFMMELAWeight(mela2, lepIdOrdered, angularOrdered, ggvvcoupl); // |ggZZ|**2

  mela2.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
  mela2.setMelaLeptonInterference(TVar::InterfOn);
  ggHZZ_prob_pure2 = getMCFMMELAWeight(mela2, lepIdOrdered, angularOrdered, ggvvcoupl); // |ggHZZ|**2

  mela2.setProcess(TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
  mela2.setMelaLeptonInterference(TVar::InterfOn);
  ggZZ_prob_Total2 = getMCFMMELAWeight(mela2, lepIdOrdered, angularOrdered, ggvvcoupl); // |ggZZ + ggHZZ|**2


  mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB);
  mela.setMelaLeptonInterference(TVar::InterfOn);
  bkg_VAMCFM3 = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |qqZZ|**2

  mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZGG);
  mela.setMelaLeptonInterference(TVar::InterfOn);
  ggzz_VAMCFM3 = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |ggZZ|**2

  mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
  mela.setMelaLeptonInterference(TVar::InterfOn);
  ggHZZ_prob_pure3 = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |ggHZZ|**2

  mela.setProcess(TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
  mela.setMelaLeptonInterference(TVar::InterfOn);
  ggZZ_prob_Total3 = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |ggZZ + ggHZZ|**2

  tout << "|ggZZ+ggHZZ|**2\t|ggZZ|**2\t|ggHZZ|**2\t|ggHZZ|**2 JHUGEN\tInterf.\t|qqZZ|**2\tPSig_m4l\tPBkg_m4l" << endl;
  tout << ggZZ_prob_Total << '\t'
    << ggzz_VAMCFM << '\t'
    << ggHZZ_prob_pure << '\t'
    << bkg_VAMCFM << '\t'
    << endl;
  tout << ggZZ_prob_Total2 << '\t'
    << ggzz_VAMCFM2 << '\t'
    << ggHZZ_prob_pure2 << '\t'
    << bkg_VAMCFM2 << '\t'
    << endl;
  tout << ggZZ_prob_Total3 << '\t'
    << ggzz_VAMCFM3 << '\t'
    << ggHZZ_prob_pure3 << '\t'
    << bkg_VAMCFM3 << '\t'
    << endl;

  tout.close();
}



void testME_FullMELA_MC(){
	int erg_tev=8;
	float mPOLE=125.6;
	char TREE_NAME[] = "GenTree";

	TVar::VerbosityLevel verbosity = TVar::INFO;

	Mela mela(erg_tev,mPOLE);
//	mela.setProcess(TVar::SelfDefine, TVar::JHUGen, TVar::GG);
//	mela.setProcess(TVar::HZZ_4l, TVar::MCFM, TVar::GG);
	TGraph* vaScale_4e = mela.vaScale_4e;
	TGraph* vaScale_4mu = mela.vaScale_4mu;
	TGraph* vaScale_2e2mu = mela.vaScale_2e2mu;
	TGraph* DggZZ_scalefactor = mela.DggZZ_scalefactor;

	char cinput[]="/afs/cern.ch/work/u/usarica/WidthAnalysis/LHC_8TeV/2mu2e/BSI/HZZ4l-125_6-8TeV-BSI.root";
	TFile* finput = new TFile(cinput);

	ofstream tout("./testME_MC_125p6_BSIEvents.txt");
	tout << "|ggZZ+ggHZZ|**2\t|ggZZ|**2\t|ggHZZ|**2\t|ggHZZ|**2 JHUGEN\tInterf.\t|qqZZ|**2\tPSig_m4l\tPBkg_m4l\t"
		 << "Lep1Id\tLep2Id\tLep3Id\tLep4Id\t"
		 << "mzz\tm1\tm2\ths\th1\th2\tphi\tphi1"
		 << endl;

	double p0plus_VAJHU;
	double bkg_VAMCFM;
	double ggzz_VAMCFM,ggZZ_prob_Total,ggHZZ_prob_Bare;
	double ggHZZ_prob_pure_NEW , ggHZZ_prob_int_NEW,ggzz_VAMCFM_NEW,ggZZ_prob_Total_NEW;
	double bkg_VAMCFM_noscale , ggzz_VAMCFM_noscale , ggHZZ_prob_pure , ggHZZ_prob_int , ggHZZ_prob_int_alt , ggHZZ_prob_pure_noscale , ggHZZ_prob_int_noscale;
	double qqScale , ggScale;
	double p0plus_m4l, p0plus_m4l_ScaleUp, p0plus_m4l_ScaleDown, p0plus_m4l_ResUp, p0plus_m4l_ResDown;
	double bkg_m4l, bkg_m4l_ScaleUp, bkg_m4l_ScaleDown, bkg_m4l_ResUp, bkg_m4l_ResDown;

	float mzz = 126.; 
	float m1 = 91.471450;
	float m2 = 12.139782;
	float h1 = 0.2682896;
	float h2 = 0.1679779;
	float phi = 1.5969792;
	float hs = -0.727181;
	float phi1 = 1.8828257;
	float GenLep1Id,GenLep2Id,GenLep3Id,GenLep4Id;

	TTree* tree = (TTree*) finput->Get(TREE_NAME);
	tree->SetBranchAddress("GenZZMass", &mzz);
	tree->SetBranchAddress("GenZ1Mass", &m1);
	tree->SetBranchAddress("GenZ2Mass", &m2);
	tree->SetBranchAddress("GenhelcosthetaZ1", &h1);
	tree->SetBranchAddress("GenhelcosthetaZ2", &h2);
	tree->SetBranchAddress("Genhelphi", &phi);
	tree->SetBranchAddress("Gencosthetastar", &hs);
	tree->SetBranchAddress("GenphistarZ1", &phi1);
	tree->SetBranchAddress("GenLep1Id", &GenLep1Id);
	tree->SetBranchAddress("GenLep2Id", &GenLep2Id);
	tree->SetBranchAddress("GenLep3Id", &GenLep3Id);
	tree->SetBranchAddress("GenLep4Id", &GenLep4Id);

	int nEntries = tree->GetEntries();
	for(int ev=0;ev<20;ev++){
		tree->GetEntry(ev);

//		int lepIdOrdered[4]={ GenLep1Id,GenLep2Id,GenLep3Id,GenLep4Id };
		int lepIdOrdered[4]={ 11,-11,11,-11 };
		float angularOrdered[8]={mzz,m1,m2,hs,h1,h2,phi,phi1};

		double selfDHvvcoupl[SIZE_HVV][2]={{0.}};
		double ggvvcoupl[2]={0,0};

		mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
		p0plus_m4l = getSuperMELA(mela,lepIdOrdered,mzz,TVar::SMSyst_None);

		mela.setProcess(TVar::bkgZZ, TVar::JHUGen, TVar::ZZGG);
		bkg_m4l = getSuperMELA(mela,lepIdOrdered,mzz,TVar::SMSyst_None);

		mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
		selfDHvvcoupl[0][0]=1;
		p0plus_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB);
		bkg_VAMCFM = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |qqZZ|**2

		mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZGG);
		ggzz_VAMCFM = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |ggZZ|**2

		mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
		ggHZZ_prob_pure = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |ggHZZ|**2

		mela.setProcess(TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
		ggZZ_prob_Total = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |ggZZ + ggHZZ|**2
		ggHZZ_prob_int = ggZZ_prob_Total - ggHZZ_prob_pure - ggzz_VAMCFM;

		mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZGG);
	mela.setMelaLeptonInterference(TVar::InterfOn);
		ggzz_VAMCFM_NEW = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |ggZZ|**2

		mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
	mela.setMelaLeptonInterference(TVar::InterfOn);
		ggHZZ_prob_pure_NEW = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |ggHZZ|**2

		mela.setProcess(TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
	mela.setMelaLeptonInterference(TVar::InterfOn);
		ggZZ_prob_Total_NEW = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |ggZZ + ggHZZ|**2
		ggHZZ_prob_int_NEW = ggZZ_prob_Total_NEW - ggHZZ_prob_pure_NEW - ggzz_VAMCFM_NEW;


		ggScale=1.0;
		qqScale=1.0;
		if(abs(lepIdOrdered[0])==abs(lepIdOrdered[1]) &&
			abs(lepIdOrdered[0])==abs(lepIdOrdered[2]) &&
			abs(lepIdOrdered[0])==abs(lepIdOrdered[3])){
				if(abs(lepIdOrdered[0])==11){
					if(mzz > 900) qqScale = vaScale_4e->Eval(900.);
					else if (mzz <  100 ) qqScale = vaScale_4e->Eval(100.);
					else qqScale = vaScale_4e->Eval(mzz);

					if(mzz > 900) ggScale = vaScale_4e->Eval(900.);
					else if (mzz <  110 ) ggScale = vaScale_4e->Eval(110.);
					else ggScale = vaScale_4e->Eval(mzz);
				}
				else{
					if(mzz > 900) qqScale = vaScale_4mu->Eval(900.);
					else if (mzz <  100 ) qqScale = vaScale_4mu->Eval(100.);
					else qqScale = vaScale_4mu->Eval(mzz);

					if(mzz > 900) ggScale = vaScale_4mu->Eval(900.);
					else if (mzz <  110 ) ggScale = vaScale_4mu->Eval(110.);
					else ggScale = vaScale_4mu->Eval(mzz);
				};
		}
		else{
			if(mzz > 900) qqScale = vaScale_2e2mu->Eval(900.);
			else if (mzz <  100 ) qqScale = vaScale_2e2mu->Eval(100.);
			else qqScale = vaScale_2e2mu->Eval(mzz);

			if(mzz > 900) ggScale = vaScale_2e2mu->Eval(900.);
			else if (mzz <  110 ) ggScale = vaScale_2e2mu->Eval(110.);
			else ggScale = vaScale_2e2mu->Eval(mzz);
		};
		if(mzz > 900) ggScale /= DggZZ_scalefactor->Eval(900.);
		else if (mzz <  110 ) ggScale /= DggZZ_scalefactor->Eval(110.);
		else ggScale /= DggZZ_scalefactor->Eval(mzz);

//		ggZZ_prob_Total /= ggScale;
//		ggzz_VAMCFM /= ggScale;
//		ggHZZ_prob_pure /= ggScale;
//		ggHZZ_prob_int /= ggScale;
//		bkg_VAMCFM /= qqScale;

		tout << ggZZ_prob_Total << '\t'
			 << ggzz_VAMCFM << '\t'
			 << ggHZZ_prob_pure << '\t'
			 << p0plus_VAJHU << '\t'
			 << ggHZZ_prob_int << '\t'
			 << bkg_VAMCFM << '\t'
			 << p0plus_m4l << '\t'
			 << bkg_m4l << '\t'
			 << GenLep1Id << '\t'
			 << GenLep2Id << '\t'
			 << GenLep3Id << '\t'
			 << GenLep4Id << '\t';
		for(int t=0;t<8;t++) tout << angularOrdered[t] << '\t';
		tout << endl;
	};

	
	tout.close();
	finput->Close();
};


void testME_FullMELA_MC_CustomWidth(){
	int erg_tev=8;
	float mPOLE=125.6;
	float wPOLE=4.15e-3;
	char TREE_NAME[] = "GenTree";

	TVar::VerbosityLevel verbosity = TVar::INFO;

	Mela mela(erg_tev,mPOLE);
//	mela.setProcess(TVar::SelfDefine, TVar::JHUGen, TVar::GG);
//	mela.setProcess(TVar::HZZ_4l, TVar::MCFM, TVar::GG);
	TGraph* vaScale_4e = mela.vaScale_4e;
	TGraph* vaScale_4mu = mela.vaScale_4mu;
	TGraph* vaScale_2e2mu = mela.vaScale_2e2mu;
	TGraph* DggZZ_scalefactor = mela.DggZZ_scalefactor;

	char cinput[]="/afs/cern.ch/work/u/usarica/WidthAnalysis/LHC_8TeV/2mu2e/BSI/HZZ4l-125_6-8TeV-BSI.root";
	TFile* finput = new TFile(cinput);

	ofstream tout("./testME_MC_125p6_BSIEvents_CustomWidth.txt");
	tout << "|ggZZ+ggHZZ|**2\t|ggZZ|**2\t|ggHZZ|**2\t|ggHZZ|**2 JHUGEN\tInterf.\t|qqZZ|**2\tPSig_m4l\tPBkg_m4l\t"
		 << "Lep1Id\tLep2Id\tLep3Id\tLep4Id\t"
		 << "mzz\tm1\tm2\ths\th1\th2\tphi\tphi1"
		 << endl;

	double p0plus_VAJHU;
	double bkg_VAMCFM;
	double ggzz_VAMCFM,ggZZ_prob_Total,ggHZZ_prob_Bare;
	double bkg_VAMCFM_noscale , ggzz_VAMCFM_noscale , ggHZZ_prob_pure , ggHZZ_prob_int , ggHZZ_prob_int_alt , ggHZZ_prob_pure_noscale , ggHZZ_prob_int_noscale;
	double qqScale , ggScale;
	double p0plus_m4l, p0plus_m4l_ScaleUp, p0plus_m4l_ScaleDown, p0plus_m4l_ResUp, p0plus_m4l_ResDown;
	double bkg_m4l, bkg_m4l_ScaleUp, bkg_m4l_ScaleDown, bkg_m4l_ResUp, bkg_m4l_ResDown;

	float mzz = 126.; 
	float m1 = 91.471450;
	float m2 = 12.139782;
	float h1 = 0.2682896;
	float h2 = 0.1679779;
	float phi = 1.5969792;
	float hs = -0.727181;
	float phi1 = 1.8828257;
	float GenLep1Id,GenLep2Id,GenLep3Id,GenLep4Id;

	TTree* tree = (TTree*) finput->Get(TREE_NAME);
	tree->SetBranchAddress("GenZZMass", &mzz);
	tree->SetBranchAddress("GenZ1Mass", &m1);
	tree->SetBranchAddress("GenZ2Mass", &m2);
	tree->SetBranchAddress("GenhelcosthetaZ1", &h1);
	tree->SetBranchAddress("GenhelcosthetaZ2", &h2);
	tree->SetBranchAddress("Genhelphi", &phi);
	tree->SetBranchAddress("Gencosthetastar", &hs);
	tree->SetBranchAddress("GenphistarZ1", &phi1);
	tree->SetBranchAddress("GenLep1Id", &GenLep1Id);
	tree->SetBranchAddress("GenLep2Id", &GenLep2Id);
	tree->SetBranchAddress("GenLep3Id", &GenLep3Id);
	tree->SetBranchAddress("GenLep4Id", &GenLep4Id);

	int nEntries = tree->GetEntries();
	for(int ev=0;ev<20;ev++){
		tree->GetEntry(ev);

//		int lepIdOrdered[4]={ GenLep1Id,GenLep2Id,GenLep3Id,GenLep4Id };
		int lepIdOrdered[4]={ 11,-11,11,-11 };
		float angularOrdered[8]={mzz,m1,m2,hs,h1,h2,phi,phi1};

		double selfDHvvcoupl[SIZE_HVV][2]={{0.}};
		double ggvvcoupl[2]={0,0};

		mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
		mela.setMelaHiggsWidth(wPOLE);
		p0plus_m4l = getSuperMELA(mela,lepIdOrdered,mzz,TVar::SMSyst_None);

		mela.setProcess(TVar::bkgZZ, TVar::JHUGen, TVar::ZZGG);
		mela.setMelaHiggsWidth(wPOLE);
		bkg_m4l = getSuperMELA(mela,lepIdOrdered,mzz,TVar::SMSyst_None);

		mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
		mela.setMelaHiggsWidth(wPOLE);
		selfDHvvcoupl[0][0]=1;
		p0plus_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB);
		mela.setMelaHiggsWidth(wPOLE);
		bkg_VAMCFM = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |qqZZ|**2

		mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZGG);
		mela.setMelaHiggsWidth(wPOLE);
		ggzz_VAMCFM = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |ggZZ|**2

		mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
		mela.setMelaHiggsWidth(wPOLE);
		ggHZZ_prob_pure = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |ggHZZ|**2

		mela.setProcess(TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
		mela.setMelaHiggsWidth(wPOLE);
		ggZZ_prob_Total = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |ggZZ + ggHZZ|**2
		ggHZZ_prob_int = ggZZ_prob_Total - ggHZZ_prob_pure - ggzz_VAMCFM;

		ggScale=1.0;
		qqScale=1.0;
		if(abs(lepIdOrdered[0])==abs(lepIdOrdered[1]) &&
			abs(lepIdOrdered[0])==abs(lepIdOrdered[2]) &&
			abs(lepIdOrdered[0])==abs(lepIdOrdered[3])){
				if(abs(lepIdOrdered[0])==11){
					if(mzz > 900) qqScale = vaScale_4e->Eval(900.);
					else if (mzz <  100 ) qqScale = vaScale_4e->Eval(100.);
					else qqScale = vaScale_4e->Eval(mzz);

					if(mzz > 900) ggScale = vaScale_4e->Eval(900.);
					else if (mzz <  110 ) ggScale = vaScale_4e->Eval(110.);
					else ggScale = vaScale_4e->Eval(mzz);
				}
				else{
					if(mzz > 900) qqScale = vaScale_4mu->Eval(900.);
					else if (mzz <  100 ) qqScale = vaScale_4mu->Eval(100.);
					else qqScale = vaScale_4mu->Eval(mzz);

					if(mzz > 900) ggScale = vaScale_4mu->Eval(900.);
					else if (mzz <  110 ) ggScale = vaScale_4mu->Eval(110.);
					else ggScale = vaScale_4mu->Eval(mzz);
				};
		}
		else{
			if(mzz > 900) qqScale = vaScale_2e2mu->Eval(900.);
			else if (mzz <  100 ) qqScale = vaScale_2e2mu->Eval(100.);
			else qqScale = vaScale_2e2mu->Eval(mzz);

			if(mzz > 900) ggScale = vaScale_2e2mu->Eval(900.);
			else if (mzz <  110 ) ggScale = vaScale_2e2mu->Eval(110.);
			else ggScale = vaScale_2e2mu->Eval(mzz);
		};
		if(mzz > 900) ggScale /= DggZZ_scalefactor->Eval(900.);
		else if (mzz <  110 ) ggScale /= DggZZ_scalefactor->Eval(110.);
		else ggScale /= DggZZ_scalefactor->Eval(mzz);

//		ggZZ_prob_Total /= ggScale;
//		ggzz_VAMCFM /= ggScale;
//		ggHZZ_prob_pure /= ggScale;
//		ggHZZ_prob_int /= ggScale;
//		bkg_VAMCFM /= qqScale;

		tout << ggZZ_prob_Total << '\t'
			 << ggzz_VAMCFM << '\t'
			 << ggHZZ_prob_pure << '\t'
			 << p0plus_VAJHU << '\t'
			 << ggHZZ_prob_int << '\t'
			 << bkg_VAMCFM << '\t'
			 << p0plus_m4l << '\t'
			 << bkg_m4l << '\t'
			 << GenLep1Id << '\t'
			 << GenLep2Id << '\t'
			 << GenLep3Id << '\t'
			 << GenLep4Id << '\t';
		for(int t=0;t<8;t++) tout << angularOrdered[t] << '\t';
		tout << endl;
	};

	
	tout.close();
	finput->Close();
};


void testME_FullMELA_FullMC_CustomWidth(){
	int erg_tev=8;
	float mPOLE=125.6;
	float wPOLE=4.15e-3;
	char TREE_NAME[] = "GenTree";

//	TVar::VerbosityLevel verbosity = TVar::INFO;

	Mela mela(erg_tev,mPOLE);
//	mela.setProcess(TVar::SelfDefine, TVar::JHUGen, TVar::GG);
//	mela.setProcess(TVar::HZZ_4l, TVar::MCFM, TVar::GG);
	TGraph* vaScale_4e = mela.vaScale_4e;
	TGraph* vaScale_4mu = mela.vaScale_4mu;
	TGraph* vaScale_2e2mu = mela.vaScale_2e2mu;
	TGraph* DggZZ_scalefactor = mela.DggZZ_scalefactor;

	char cinput[]="/afs/cern.ch/work/u/usarica/HZZ4l-125p6-FullAnalysis/LHC_8TeV/GenSignal/HZZ4lTree_powheg15jhuGenV3-0L1f05ph0H125.6_Generated.root";
	char coutput[]="HZZ4lTree_powheg15jhuGenV3-0L1f05ph0H125.6_Generated_TestOnly.root";
	TFile* finput = new TFile(cinput,"read");
	TFile* foutput = new TFile(coutput,"recreate");

	double p0plus_VAJHU;
	double p0minus_VAJHU;
	double p0_g1prime2_VAJHU;
	double pg1g4_VAJHU;
	double pg1g4_pi2_VAJHU;
	double pg1g1prime2_VAJHU;

	double p0plus_VAJHU_NEW;
	double p0minus_VAJHU_NEW;
	double p0_g1prime2_VAJHU_NEW;
	double pg1g4_VAJHU_NEW;
	double pg1g4_pi2_VAJHU_NEW;
	double pg1g1prime2_VAJHU_NEW;

	double bkg_VAMCFM;
	double ggzz_VAMCFM,ggZZ_prob_Total,ggHZZ_prob_Bare;
	double bkg_VAMCFM_noscale , ggzz_VAMCFM_noscale , ggHZZ_prob_pure , ggHZZ_prob_int , ggHZZ_prob_int_alt , ggHZZ_prob_pure_noscale , ggHZZ_prob_int_noscale;
	double qqScale , ggScale;
	double p0plus_m4l, p0plus_m4l_ScaleUp, p0plus_m4l_ScaleDown, p0plus_m4l_ResUp, p0plus_m4l_ResDown;
	double bkg_m4l, bkg_m4l_ScaleUp, bkg_m4l_ScaleDown, bkg_m4l_ResUp, bkg_m4l_ResDown;

	float mzz = 126.; 
	float m1 = 91.471450;
	float m2 = 12.139782;
	float h1 = 0.2682896;
	float h2 = 0.1679779;
	float phi = 1.5969792;
	float hs = -0.727181;
	float phi1 = 1.8828257;
	float GenLep1Id,GenLep2Id,GenLep3Id,GenLep4Id;

	TTree* tree = (TTree*) finput->Get(TREE_NAME);
	tree->SetBranchAddress("GenZZMass", &mzz);
	tree->SetBranchAddress("GenZ1Mass", &m1);
	tree->SetBranchAddress("GenZ2Mass", &m2);
	tree->SetBranchAddress("GenhelcosthetaZ1", &h1);
	tree->SetBranchAddress("GenhelcosthetaZ2", &h2);
	tree->SetBranchAddress("Genhelphi", &phi);
	tree->SetBranchAddress("Gencosthetastar", &hs);
	tree->SetBranchAddress("GenphistarZ1", &phi1);
	tree->SetBranchAddress("GenLep1Id", &GenLep1Id);
	tree->SetBranchAddress("GenLep2Id", &GenLep2Id);
	tree->SetBranchAddress("GenLep3Id", &GenLep3Id);
	tree->SetBranchAddress("GenLep4Id", &GenLep4Id);

	TTree* newtree = new TTree("TestTree","");
	newtree->Branch("p0plus_VAJHU",&p0plus_VAJHU);
	newtree->Branch("p0minus_VAJHU",&p0minus_VAJHU);
	newtree->Branch("p0_g1prime2_VAJHU",&p0_g1prime2_VAJHU);
	newtree->Branch("pg1g4_VAJHU",&pg1g4_VAJHU);
	newtree->Branch("pg1g4_pi2_VAJHU",&pg1g4_pi2_VAJHU);
	newtree->Branch("pg1g1prime2_VAJHU",&pg1g1prime2_VAJHU);
	newtree->Branch("p0plus_VAJHU_NEW",&p0plus_VAJHU_NEW);
	newtree->Branch("p0minus_VAJHU_NEW",&p0minus_VAJHU_NEW);
	newtree->Branch("p0_g1prime2_VAJHU_NEW",&p0_g1prime2_VAJHU_NEW);
	newtree->Branch("pg1g4_VAJHU_NEW",&pg1g4_VAJHU_NEW);
	newtree->Branch("pg1g4_pi2_VAJHU_NEW",&pg1g4_pi2_VAJHU_NEW);
	newtree->Branch("pg1g1prime2_VAJHU_NEW",&pg1g1prime2_VAJHU_NEW);

	int nEntries = tree->GetEntries();
	double selfDHvvcoupl[SIZE_HVV][2] = { { 0. } };
	mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
	for (int ev = 0; ev < nEntries; ev++){
		tree->GetEntry(ev);

		int lepIdOrdered[4] = { GenLep1Id, GenLep2Id, GenLep3Id, GenLep4Id };
//		int lepIdOrdered[4]={ 11,-11,11,-11 };
		float angularOrdered[8] = { mzz, m1, m2, hs, h1, h2, phi, phi1 };

		selfDHvvcoupl[0][0] = 1;
		selfDHvvcoupl[3][0] = 0;
		selfDHvvcoupl[3][1] = 0;
		selfDHvvcoupl[11][0] = 0;
		p0plus_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		selfDHvvcoupl[0][0] = 0;
		selfDHvvcoupl[3][0] = 1;
		selfDHvvcoupl[3][1] = 0;
		selfDHvvcoupl[11][0] = 0;
		p0minus_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		selfDHvvcoupl[0][0] = 0;
		selfDHvvcoupl[3][0] = 0;
		selfDHvvcoupl[3][1] = 0;
		selfDHvvcoupl[11][0] = 1;
		p0_g1prime2_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		selfDHvvcoupl[0][0] = 1;
		selfDHvvcoupl[3][0] = 2.521;
		selfDHvvcoupl[3][1] = 0;
		selfDHvvcoupl[11][0] = 0;
		pg1g4_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		selfDHvvcoupl[0][0] = 1;
		selfDHvvcoupl[3][0] = 0;
		selfDHvvcoupl[3][1] = 2.521;
		selfDHvvcoupl[11][0] = 0;
		pg1g4_pi2_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		selfDHvvcoupl[0][0] = 1;
		selfDHvvcoupl[3][0] = 0;
		selfDHvvcoupl[3][1] = 0;
		selfDHvvcoupl[11][0] = -12046.01;
		pg1g1prime2_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		selfDHvvcoupl[0][0] = 1;
		selfDHvvcoupl[3][0] = 0;
		selfDHvvcoupl[3][1] = 0;
		selfDHvvcoupl[11][0] = 0;
		mela.setMelaHiggsWidth(wPOLE);
		p0plus_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		selfDHvvcoupl[0][0] = 0;
		selfDHvvcoupl[3][0] = 1;
		selfDHvvcoupl[3][1] = 0;
		selfDHvvcoupl[11][0] = 0;
		mela.setMelaHiggsWidth(wPOLE);
		p0minus_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		selfDHvvcoupl[0][0] = 0;
		selfDHvvcoupl[3][0] = 0;
		selfDHvvcoupl[3][1] = 0;
		selfDHvvcoupl[11][0] = 1;
		mela.setMelaHiggsWidth(wPOLE);
		p0_g1prime2_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		selfDHvvcoupl[0][0] = 1;
		selfDHvvcoupl[3][0] = 2.521;
		selfDHvvcoupl[3][1] = 0;
		selfDHvvcoupl[11][0] = 0;
		mela.setMelaHiggsWidth(wPOLE);
		pg1g4_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		selfDHvvcoupl[0][0] = 1;
		selfDHvvcoupl[3][0] = 0;
		selfDHvvcoupl[3][1] = 2.521;
		selfDHvvcoupl[11][0] = 0;
		mela.setMelaHiggsWidth(wPOLE);
		pg1g4_pi2_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		selfDHvvcoupl[0][0] = 1;
		selfDHvvcoupl[3][0] = 0;
		selfDHvvcoupl[3][1] = 0;
		selfDHvvcoupl[11][0] = -12046.01;
		mela.setMelaHiggsWidth(wPOLE);
		pg1g1prime2_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		newtree->Fill();
	};


	foutput->WriteTObject(newtree);
	foutput->Close();
	finput->Close();
};


void testME_FullMELA_FullMC_LeptonInterf(){
	int erg_tev=8;
	float mPOLE=125.6;
	float wPOLE=4.15e-3;
	char TREE_NAME[] = "SelectedTree";

	Mela mela(erg_tev,mPOLE);

//	char cinput[]="/afs/cern.ch/work/u/usarica/WidthAnalysis/LHC_8TeV/2mu2e/QQZZ/HZZ4l-125_6-8TeV-QQZZ.root";
//	char cinput[]="/afs/cern.ch/work/u/usarica/WidthAnalysis/LHC_8TeV/2mu2e/Sig/HZZ4l-125_6-8TeV-Sig.root";
//  char cinput[]="/afs/cern.ch/work/u/usarica/WidthAnalysis/LHC_8TeV/4e/Sig/HZZ4l-125_6-8TeV-Sig.root";
  TString cinput_main = "/scratch0/hep/ianderso/CJLST/140519/PRODFSR";
  if (erg_tev==8) cinput_main.Append("_8TeV");
  TString cinput = Form("%s/4e/HZZ4lTree_ggTo4e_Contin-MCFM67.root", cinput_main.Data());
  char coutput[]="HZZ4l-125_6-8TeV-Sig_TestOnly.root";
	TFile* finput = new TFile(cinput,"read");
	TFile* foutput = new TFile(coutput,"recreate");

	double p0plus_VAJHU,ggHZZ_prob_pure,ggHZZ_prob_int,ggZZ_prob_Total;
	double bkg_VAMCFM,ggzz_VAMCFM;
	double p0plus_VAJHU_NEW,ggHZZ_prob_pure_NEW,ggHZZ_prob_int_NEW,ggZZ_prob_Total_NEW;
	double bkg_VAMCFM_NEW,ggzz_VAMCFM_NEW;
	double bkg_VAMCFM_STU,bkg_VAMCFM_TU,bkg_VAMCFM_S;
  float Dgg10;

	double p0plus_m4l, p0plus_m4l_ScaleUp, p0plus_m4l_ScaleDown, p0plus_m4l_ResUp, p0plus_m4l_ResDown;
	double bkg_m4l, bkg_m4l_ScaleUp, bkg_m4l_ScaleDown, bkg_m4l_ResUp, bkg_m4l_ResDown;

	float mzz = 126.; 
	float m1 = 91.471450;
	float m2 = 12.139782;
	float h1 = 0.2682896;
	float h2 = 0.1679779;
	float phi = 1.5969792;
	float hs = -0.727181;
	float phi1 = 1.8828257;
	float GenLep1Id,GenLep2Id,GenLep3Id,GenLep4Id;

	TTree* tree = (TTree*) finput->Get(TREE_NAME);
	tree->SetBranchAddress("ZZMass", &mzz);
	tree->SetBranchAddress("Z1Mass", &m1);
	tree->SetBranchAddress("Z2Mass", &m2);
	tree->SetBranchAddress("helcosthetaZ1", &h1);
	tree->SetBranchAddress("helcosthetaZ2", &h2);
	tree->SetBranchAddress("helphi", &phi);
	tree->SetBranchAddress("costhetastar", &hs);
	tree->SetBranchAddress("phistarZ1", &phi1);


	TTree* newtree = new TTree("TestTree","");
	newtree->Branch("ZZMass", &mzz);
	newtree->Branch("Z1Mass", &m1);
	newtree->Branch("Z2Mass", &m2);
	newtree->Branch("helcosthetaZ1", &h1);
	newtree->Branch("helcosthetaZ2", &h2);
	newtree->Branch("helphi", &phi);
	newtree->Branch("costhetastar", &hs);
	newtree->Branch("phistarZ1", &phi1);

	newtree->Branch("p0plus_VAJHU",&p0plus_VAJHU);
	newtree->Branch("bkg_VAMCFM",&bkg_VAMCFM);
	newtree->Branch("ggzz_VAMCFM",&ggzz_VAMCFM);
	newtree->Branch("ggHZZ_prob_pure",&ggHZZ_prob_pure);
	newtree->Branch("ggHZZ_prob_int",&ggHZZ_prob_int);

	newtree->Branch("p0plus_VAJHU_NEW",&p0plus_VAJHU_NEW);
	newtree->Branch("p0plus_VAJHU_NEW",&p0plus_VAJHU_NEW);
	newtree->Branch("bkg_VAMCFM_NEW",&bkg_VAMCFM_NEW);
	newtree->Branch("ggzz_VAMCFM_NEW",&ggzz_VAMCFM_NEW);
	newtree->Branch("ggHZZ_prob_pure_NEW",&ggHZZ_prob_pure_NEW);
	newtree->Branch("ggHZZ_prob_int_NEW",&ggHZZ_prob_int_NEW);

	newtree->Branch("bkg_VAMCFM_STU",&bkg_VAMCFM_STU);
	newtree->Branch("bkg_VAMCFM_S",&bkg_VAMCFM_S);
	newtree->Branch("bkg_VAMCFM_TU",&bkg_VAMCFM_TU);

  newtree->Branch("Dgg10", &Dgg10);

	int nEntries = tree->GetEntries();
	double selfDHvvcoupl[SIZE_HVV][2] = { { 0. } };
	double ggvvcoupl[2]={0,0};
	mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
//	for (int ev = 0; ev < nEntries; ev++){
  int ctr=0;
  int ev = 0;
  while (ctr<10000){
    if (ev==tree->GetEntries()) break;
    tree->GetEntry(ev);
    ev++;
    if (mzz>200) ctr++;
//		int lepIdOrdered[4] = { GenLep1Id, GenLep2Id, GenLep3Id, GenLep4Id };
		int lepIdOrdered[4]={ 11,-11,11,-11 };
		float angularOrdered[8] = { mzz, m1, m2, hs, h1, h2, phi, phi1 };
/*
		mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
		selfDHvvcoupl[0][0]=1;
		p0plus_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
*/

		mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB);
		bkg_VAMCFM = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |qqZZ|**2

		mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZGG);
		ggzz_VAMCFM = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |ggZZ|**2

		mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
		ggHZZ_prob_pure = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |ggHZZ|**2

		mela.setProcess(TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
		ggZZ_prob_Total = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |ggZZ + ggHZZ|**2
		ggHZZ_prob_int = ggZZ_prob_Total - ggHZZ_prob_pure - ggzz_VAMCFM;

		mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
		selfDHvvcoupl[0][0]=1;
		mela.setMelaLeptonInterference(TVar::InterfOn);
		p0plus_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB);
		mela.setMelaLeptonInterference(TVar::InterfOn);
		bkg_VAMCFM_NEW = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |qqZZ|**2

  	mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZGG);
		mela.setMelaLeptonInterference(TVar::InterfOn);
		ggzz_VAMCFM_NEW = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |ggZZ|**2

		mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
		mela.setMelaLeptonInterference(TVar::InterfOn);
		ggHZZ_prob_pure_NEW = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |ggHZZ|**2

		mela.setProcess(TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
		mela.setMelaLeptonInterference(TVar::InterfOn);
		ggZZ_prob_Total_NEW = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |ggZZ + ggHZZ|**2
		ggHZZ_prob_int_NEW = ggZZ_prob_Total_NEW - ggHZZ_prob_pure_NEW - ggzz_VAMCFM_NEW;

		mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB_STU);
		bkg_VAMCFM_STU = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |qqZZ|**2

		mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB_S);
		bkg_VAMCFM_S = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |qqZZ|**2

		mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB_TU);
		bkg_VAMCFM_TU = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |qqZZ|**2

		int lepIdOrderedp[4]={ 13,-13,11,-11 };
		mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
		selfDHvvcoupl[0][0]=1;
		p0plus_VAJHU = getJHUGenMELAWeight(mela, lepIdOrderedp, angularOrdered, selfDHvvcoupl)*0.5;

    mela.setProcess(TVar::D_gg10, TVar::MCFM, TVar::ZZGG);
    mela.computeD_gg(angularOrdered[0], angularOrdered[1], angularOrdered[2], angularOrdered[3],
      angularOrdered[4], angularOrdered[5], angularOrdered[6], angularOrdered[7], 1, TVar::MCFM, TVar::D_gg10, Dgg10);

		newtree->Fill();
	};


	foutput->WriteTObject(newtree);
	foutput->Close();
	finput->Close();
};


void testME_FullMELA_FullMC_ZZQQBSTU(){
	int erg_tev=8;
	float mPOLE=125.6;
	float wPOLE=4.15e-3;
	char TREE_NAME[] = "GenTree";

//	TVar::VerbosityLevel verbosity = TVar::INFO;

	Mela mela(erg_tev,mPOLE);
//	mela.setProcess(TVar::SelfDefine, TVar::JHUGen, TVar::GG);
//	mela.setProcess(TVar::HZZ_4l, TVar::MCFM, TVar::GG);
	TGraph* vaScale_4e = mela.vaScale_4e;
	TGraph* vaScale_4mu = mela.vaScale_4mu;
	TGraph* vaScale_2e2mu = mela.vaScale_2e2mu;
	TGraph* DggZZ_scalefactor = mela.DggZZ_scalefactor;

//	char cinput[]="/afs/cern.ch/work/u/usarica/WidthAnalysis/LHC_8TeV/2mu2e/QQZZ/HZZ4l-125_6-8TeV-QQZZ.root";
	char cinput[]="/afs/cern.ch/work/u/usarica/WidthAnalysis/LHC_8TeV/4e/QQZZ/HZZ4l-125_6-8TeV-QQZZ.root";
	char coutput[]="HZZ4l-125_6-8TeV-QQZZ_TestOnly.root";
	TFile* finput = new TFile(cinput,"read");
	TFile* foutput = new TFile(coutput,"recreate");

	double p0plus_VAJHU,ggHZZ_prob_pure,ggHZZ_prob_int,ggZZ_prob_Total;
	double bkg_VAMCFM,ggzz_VAMCFM;
	double p0plus_VAJHU_NEW,ggHZZ_prob_pure_NEW,ggHZZ_prob_int_NEW,ggZZ_prob_Total_NEW;
	double bkg_VAMCFM_NEW,ggzz_VAMCFM_NEW;
	double bkg_VAMCFM_STU,bkg_VAMCFM_TU,bkg_VAMCFM_S;

	double p0plus_m4l, p0plus_m4l_ScaleUp, p0plus_m4l_ScaleDown, p0plus_m4l_ResUp, p0plus_m4l_ResDown;
	double bkg_m4l, bkg_m4l_ScaleUp, bkg_m4l_ScaleDown, bkg_m4l_ResUp, bkg_m4l_ResDown;

	float mzz = 126.; 
	float m1 = 91.471450;
	float m2 = 12.139782;
	float h1 = 0.2682896;
	float h2 = 0.1679779;
	float phi = 1.5969792;
	float hs = -0.727181;
	float phi1 = 1.8828257;
	float GenLep1Id,GenLep2Id,GenLep3Id,GenLep4Id;

	TTree* tree = (TTree*) finput->Get(TREE_NAME);
	tree->SetBranchAddress("GenZZMass", &mzz);
	tree->SetBranchAddress("GenZ1Mass", &m1);
	tree->SetBranchAddress("GenZ2Mass", &m2);
	tree->SetBranchAddress("GenhelcosthetaZ1", &h1);
	tree->SetBranchAddress("GenhelcosthetaZ2", &h2);
	tree->SetBranchAddress("Genhelphi", &phi);
	tree->SetBranchAddress("Gencosthetastar", &hs);
	tree->SetBranchAddress("GenphistarZ1", &phi1);
	tree->SetBranchAddress("GenLep1Id", &GenLep1Id);
	tree->SetBranchAddress("GenLep2Id", &GenLep2Id);
	tree->SetBranchAddress("GenLep3Id", &GenLep3Id);
	tree->SetBranchAddress("GenLep4Id", &GenLep4Id);

	TTree* newtree = new TTree("TestTree","");
	newtree->Branch("GenZZMass", &mzz);
	newtree->Branch("GenZ1Mass", &m1);
	newtree->Branch("GenZ2Mass", &m2);
	newtree->Branch("GenhelcosthetaZ1", &h1);
	newtree->Branch("GenhelcosthetaZ2", &h2);
	newtree->Branch("Genhelphi", &phi);
	newtree->Branch("Gencosthetastar", &hs);
	newtree->Branch("GenphistarZ1", &phi1);
	newtree->Branch("GenLep1Id", &GenLep1Id);
	newtree->Branch("GenLep2Id", &GenLep2Id);
	newtree->Branch("GenLep3Id", &GenLep3Id);
	newtree->Branch("GenLep4Id", &GenLep4Id);
	newtree->Branch("bkg_VAMCFM_STU",&bkg_VAMCFM_STU);
	newtree->Branch("bkg_VAMCFM_S",&bkg_VAMCFM_S);
	newtree->Branch("bkg_VAMCFM_TU",&bkg_VAMCFM_TU);
	newtree->Branch("bkg_VAMCFM",&bkg_VAMCFM);

	int nEntries = tree->GetEntries();
	double selfDHvvcoupl[SIZE_HVV][2] = { { 0. } };
	double ggvvcoupl[2]={0,0};
	mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
	for (int ev = 0; ev < nEntries; ev++){
//	for (int ev = 0; ev < 5; ev++){
		tree->GetEntry(ev);

		int lepIdOrdered[4] = { GenLep1Id, GenLep2Id, GenLep3Id, GenLep4Id };
//		int lepIdOrdered[4]={ 11,-11,11,-11 };
		float angularOrdered[8] = { mzz, m1, m2, hs, h1, h2, phi, phi1 };


		mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB);
		bkg_VAMCFM = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |qqZZ|**2

		mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB_STU);
		bkg_VAMCFM_STU = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |qqZZ|**2

		mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB_S);
		bkg_VAMCFM_S = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |qqZZ|**2

		mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB_TU);
		bkg_VAMCFM_TU = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |qqZZ|**2

		newtree->Fill();
	};


	foutput->WriteTObject(newtree);
	foutput->Close();
	finput->Close();
};


void testME_FullMELA_Zgsggs(){
	int erg_tev=8;
	float mPOLE=125.6;
	float wPOLE=4.15e-3;
	char TREE_NAME[] = "GenTree";

//	TVar::VerbosityLevel verbosity = TVar::INFO;

	Mela mela(erg_tev,mPOLE);
	TGraph* vaScale_4e = mela.vaScale_4e;
	TGraph* vaScale_4mu = mela.vaScale_4mu;
	TGraph* vaScale_2e2mu = mela.vaScale_2e2mu;
	TGraph* DggZZ_scalefactor = mela.DggZZ_scalefactor;

	TFile* finput = new TFile("/afs/cern.ch/work/u/usarica/HZZ4l-125p6-FullAnalysis/LHC_8TeV/GenSignal/HZZ4lTree_powheg15jhuGenV3-0PMH125.6_Generated.root","read");
	TFile* foutput = new TFile("HZZ4lTree_powheg15jhuGenV3-0PMH125.6_ZgsggsTest.root","recreate");

	float p0plus_VAJHU;
	float p0hplus_VAJHU;
	float p0minus_VAJHU;
	float p0_g1prime2_VAJHU;
	float pg1g1prime2_VAJHU;
	float pg1g1prime2_pi2_VAJHU;
	float pg1g2_VAJHU;
	float pg1g2_pi2_VAJHU;
	float pg1g4_VAJHU;
	float pg1g4_pi2_VAJHU;
	float p0Zgs_g1prime2_VAJHU;
	float pzzzgs_g1prime2_VAJHU;
	float pzzzgs_g1prime2_pi2_VAJHU;
	float p0_g1prime2_0Zgs_g1prime2_VAJHU;
	float p0_g1prime2_0Zgs_g1prime2_pi2_VAJHU;

	float p0plus_VAJHU_NEW;
	float p0hplus_VAJHU_NEW;
	float p0minus_VAJHU_NEW;
	float p0_g1prime2_VAJHU_NEW;
	float pg1g1prime2_VAJHU_NEW;
	float pg1g1prime2_pi2_VAJHU_NEW;
	float pg1g2_VAJHU_NEW;
	float pg1g2_pi2_VAJHU_NEW;
	float pg1g4_VAJHU_NEW;
	float pg1g4_pi2_VAJHU_NEW;

	float ggHZZ_prob_pure,ggHZZ_prob_int,ggZZ_prob_Total;
	float bkg_VAMCFM,ggzz_VAMCFM;
	float ggHZZ_prob_pure_NEW,ggHZZ_prob_int_NEW,ggZZ_prob_Total_NEW;
	float bkg_VAMCFM_NEW,ggzz_VAMCFM_NEW;
	float bkg_VAMCFM_STU,bkg_VAMCFM_TU,bkg_VAMCFM_S;

	float p0plus_m4l, p0plus_m4l_ScaleUp, p0plus_m4l_ScaleDown, p0plus_m4l_ResUp, p0plus_m4l_ResDown;
	float bkg_m4l, bkg_m4l_ScaleUp, bkg_m4l_ScaleDown, bkg_m4l_ResUp, bkg_m4l_ResDown;

	float mzz = 126.; 
	float m1 = 91.471450;
	float m2 = 12.139782;
	float h1 = 0.2682896;
	float h2 = 0.1679779;
	float phi = 1.5969792;
	float hs = -0.727181;
	float phi1 = 1.8828257;
	float GenLep1Id,GenLep2Id,GenLep3Id,GenLep4Id;

	TTree* tree = (TTree*) finput->Get(TREE_NAME);
	tree->SetBranchAddress("GenZZMass", &mzz);
	tree->SetBranchAddress("GenZ1Mass", &m1);
	tree->SetBranchAddress("GenZ2Mass", &m2);
	tree->SetBranchAddress("GenhelcosthetaZ1", &h1);
	tree->SetBranchAddress("GenhelcosthetaZ2", &h2);
	tree->SetBranchAddress("Genhelphi", &phi);
	tree->SetBranchAddress("Gencosthetastar", &hs);
	tree->SetBranchAddress("GenphistarZ1", &phi1);
	tree->SetBranchAddress("GenLep1Id", &GenLep1Id);
	tree->SetBranchAddress("GenLep2Id", &GenLep2Id);
	tree->SetBranchAddress("GenLep3Id", &GenLep3Id);
	tree->SetBranchAddress("GenLep4Id", &GenLep4Id);

	TTree* newtree = new TTree("TestTree","");
	newtree->Branch("GenZZMass", &mzz);
	newtree->Branch("GenZ1Mass", &m1);
	newtree->Branch("GenZ2Mass", &m2);
	newtree->Branch("GenhelcosthetaZ1", &h1);
	newtree->Branch("GenhelcosthetaZ2", &h2);
	newtree->Branch("Genhelphi", &phi);
	newtree->Branch("Gencosthetastar", &hs);
	newtree->Branch("GenphistarZ1", &phi1);
	newtree->Branch("GenLep1Id", &GenLep1Id);
	newtree->Branch("GenLep2Id", &GenLep2Id);
	newtree->Branch("GenLep3Id", &GenLep3Id);
	newtree->Branch("GenLep4Id", &GenLep4Id);

	newtree->Branch("p0plus_VAJHU",&p0plus_VAJHU);
	newtree->Branch("p0Zgs_g1prime2_VAJHU",&p0Zgs_g1prime2_VAJHU);
	newtree->Branch("p0_g1prime2_VAJHU",&p0_g1prime2_VAJHU);
	newtree->Branch("pzzzgs_g1prime2_VAJHU",&pzzzgs_g1prime2_VAJHU);
	newtree->Branch("pzzzgs_g1prime2_pi2_VAJHU",&pzzzgs_g1prime2_pi2_VAJHU);
	newtree->Branch("pg1g1prime2_VAJHU",&pg1g1prime2_VAJHU);
	newtree->Branch("pg1g1prime2_pi2_VAJHU",&pg1g1prime2_pi2_VAJHU);
	newtree->Branch("p0_g1prime2_0Zgs_g1prime2_VAJHU",&p0_g1prime2_0Zgs_g1prime2_VAJHU);
	newtree->Branch("p0_g1prime2_0Zgs_g1prime2_pi2_VAJHU",&p0_g1prime2_0Zgs_g1prime2_pi2_VAJHU);

	int nEntries = tree->GetEntries();
	double selfDHvvcoupl[SIZE_HVV][2] = { { 0 } };
	double ggvvcoupl[2]={0,0};
	mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
	for (int ev = 0; ev < nEntries; ev++){
//	for (int ev = 0; ev < 100; ev++){
//	for (int ev = 0; ev < 10000; ev++){
		tree->GetEntry(ev);

		int lepIdOrdered[4] = { GenLep1Id, GenLep2Id, GenLep3Id, GenLep4Id };
		float angularOrdered[8] = { mzz, m1, m2, hs, h1, h2, phi, phi1 };

		mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[0][0]=1;
		p0plus_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[30][0]=1;
		p0Zgs_g1prime2_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[30][0]=-7591.914;
		pzzzgs_g1prime2_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[30][1]=-7591.914;
		pzzzgs_g1prime2_pi2_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[11][0]=-12046.01;
		p0_g1prime2_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[11][0]=-12046.01;
		pg1g1prime2_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[11][1]=-12046.01;
		pg1g1prime2_pi2_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[11][0]=-12046.01;
		selfDHvvcoupl[30][0]=-7591.914;
		p0_g1prime2_0Zgs_g1prime2_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[11][0]=-12046.01;
		selfDHvvcoupl[30][1]=-7591.914;
		p0_g1prime2_0Zgs_g1prime2_pi2_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		p0Zgs_g1prime2_VAJHU /= p0plus_VAJHU;
		pzzzgs_g1prime2_VAJHU /= p0plus_VAJHU;
		pzzzgs_g1prime2_pi2_VAJHU /= p0plus_VAJHU;
		p0_g1prime2_VAJHU /= p0plus_VAJHU;
		pg1g1prime2_VAJHU /= p0plus_VAJHU;
		pg1g1prime2_pi2_VAJHU /= p0plus_VAJHU;
		p0_g1prime2_0Zgs_g1prime2_VAJHU /= p0plus_VAJHU;
		p0_g1prime2_0Zgs_g1prime2_pi2_VAJHU /= p0plus_VAJHU;

		newtree->Fill();
	};


	foutput->WriteTObject(newtree);
	foutput->Close();
	finput->Close();
};


void testME_FullMELA_FullSimMCValidation(int flavor=1){
	int erg_tev=8;
	float mPOLE=125.6;
	float wPOLE=4.15e-3;
	char TREE_NAME[] = "SelectedTree";

//	TVar::VerbosityLevel verbosity = TVar::INFO;

	Mela mela(erg_tev,mPOLE);
	TGraph* vaScale_4e = mela.vaScale_4e;
	TGraph* vaScale_4mu = mela.vaScale_4mu;
	TGraph* vaScale_2e2mu = mela.vaScale_2e2mu;
	TGraph* DggZZ_scalefactor = mela.DggZZ_scalefactor;

  TString cinput_main = "/scratch0/hep/ianderso/CJLST/140519/PRODFSR_8TeV";
  TFile* finput = new TFile(Form("%s/%s/HZZ4lTree_powheg15jhuGenV3-0PMH125.6.root", cinput_main.Data(), (flavor==0 ? "2mu2e" : "4e")), "read");
	//TFile* finput = new TFile(Form("/afs/cern.ch/work/u/usarica/HZZ4l-125p6-FullAnalysis/LHC_8TeV/%s/HZZ4lTree_powheg15jhuGenV3-0PMH125.6.root",(flavor==0 ? "2mu2e" : "4e")),"read");
	TFile* foutput = new TFile(Form("HZZ4lTree_powheg15jhuGenV3-0PMH125.6_%s_OriginalMEValidationTestOnly.root",(flavor==0 ? "2mu2e" : "4e")),"recreate");

	float p0plus_VAJHU;
	float p0hplus_VAJHU;
	float p0minus_VAJHU;
	float p0_g1prime2_VAJHU;
	float pg1g1prime2_VAJHU;
	float pg1g1prime2_pi2_VAJHU;
	float pg1g2_VAJHU;
	float pg1g2_pi2_VAJHU;
	float pg1g4_VAJHU;
	float pg1g4_pi2_VAJHU;

	float p0plus_VAJHU_NEW;
	float p0hplus_VAJHU_NEW;
	float p0minus_VAJHU_NEW;
	float p0_g1prime2_VAJHU_NEW;
	float pg1g1prime2_VAJHU_NEW;
	float pg1g1prime2_pi2_VAJHU_NEW;
	float pg1g2_VAJHU_NEW;
	float pg1g2_pi2_VAJHU_NEW;
	float pg1g4_VAJHU_NEW;
	float pg1g4_pi2_VAJHU_NEW;


	float ggHZZ_prob_pure,ggHZZ_prob_int,ggZZ_prob_Total;
	float bkg_VAMCFM,ggzz_VAMCFM;
	float ggHZZ_prob_pure_NEW,ggHZZ_prob_int_NEW,ggZZ_prob_Total_NEW,p0plus_VAMCFM_NEW_BSMOn;
	float bkg_VAMCFM_NEW,ggzz_VAMCFM_NEW;
	float bkg_VAMCFM_STU,bkg_VAMCFM_TU,bkg_VAMCFM_S;

	float p0plus_m4l, p0plus_m4l_ScaleUp, p0plus_m4l_ScaleDown, p0plus_m4l_ResUp, p0plus_m4l_ResDown;
	float bkg_m4l, bkg_m4l_ScaleUp, bkg_m4l_ScaleDown, bkg_m4l_ResUp, bkg_m4l_ResDown;

	float mzz = 126.; 
	float m1 = 91.471450;
	float m2 = 12.139782;
	float h1 = 0.2682896;
	float h2 = 0.1679779;
	float phi = 1.5969792;
	float hs = -0.727181;
	float phi1 = 1.8828257;

	TTree* tree = (TTree*) finput->Get(TREE_NAME);
	tree->SetBranchAddress("ZZMass", &mzz);
	tree->SetBranchAddress("Z1Mass", &m1);
	tree->SetBranchAddress("Z2Mass", &m2);
	tree->SetBranchAddress("helcosthetaZ1", &h1);
	tree->SetBranchAddress("helcosthetaZ2", &h2);
	tree->SetBranchAddress("helphi", &phi);
	tree->SetBranchAddress("costhetastar", &hs);
	tree->SetBranchAddress("phistarZ1", &phi1);
	tree->SetBranchAddress("p0plus_VAJHU",&p0plus_VAJHU);
	tree->SetBranchAddress("p0hplus_VAJHU",&p0hplus_VAJHU);
	tree->SetBranchAddress("p0minus_VAJHU",&p0minus_VAJHU);
	tree->SetBranchAddress("p0_g1prime2_VAJHU",&p0_g1prime2_VAJHU);
	tree->SetBranchAddress("pg1g2_VAJHU",&pg1g2_VAJHU);
	tree->SetBranchAddress("pg1g2_pi2_VAJHU",&pg1g2_pi2_VAJHU);
	tree->SetBranchAddress("pg1g4_VAJHU",&pg1g4_VAJHU);
	tree->SetBranchAddress("pg1g4_pi2_VAJHU",&pg1g4_pi2_VAJHU);
	tree->SetBranchAddress("pg1g1prime2_VAJHU",&pg1g1prime2_VAJHU);
	tree->SetBranchAddress("bkg_VAMCFM",&bkg_VAMCFM);
	tree->SetBranchAddress("ggzz_VAMCFM",&ggzz_VAMCFM);
	tree->SetBranchAddress("p0plus_VAMCFM",&ggHZZ_prob_pure);
	tree->SetBranchAddress("ggzz_p0plus_VAMCFM",&ggHZZ_prob_int);

	TTree* newtree = new TTree("TestTree","");
	newtree->Branch("ZZMass", &mzz);
	newtree->Branch("Z1Mass", &m1);
	newtree->Branch("Z2Mass", &m2);
	newtree->Branch("helcosthetaZ1", &h1);
	newtree->Branch("helcosthetaZ2", &h2);
	newtree->Branch("helphi", &phi);
	newtree->Branch("costhetastar", &hs);
	newtree->Branch("phistarZ1", &phi1);

	newtree->Branch("p0plus_VAJHU",&p0plus_VAJHU);
	newtree->Branch("p0hplus_VAJHU",&p0hplus_VAJHU);
	newtree->Branch("p0minus_VAJHU",&p0minus_VAJHU);
	newtree->Branch("p0_g1prime2_VAJHU",&p0_g1prime2_VAJHU);
	newtree->Branch("pg1g2_VAJHU",&pg1g2_VAJHU);
	newtree->Branch("pg1g2_pi2_VAJHU",&pg1g2_pi2_VAJHU);
	newtree->Branch("pg1g4_VAJHU",&pg1g4_VAJHU);
	newtree->Branch("pg1g4_pi2_VAJHU",&pg1g4_pi2_VAJHU);
	newtree->Branch("pg1g1prime2_VAJHU",&pg1g1prime2_VAJHU);

	newtree->Branch("p0plus_VAJHU_NEW",&p0plus_VAJHU_NEW);
	newtree->Branch("p0hplus_VAJHU_NEW",&p0hplus_VAJHU_NEW);
	newtree->Branch("p0minus_VAJHU_NEW",&p0minus_VAJHU_NEW);
	newtree->Branch("p0_g1prime2_VAJHU_NEW",&p0_g1prime2_VAJHU_NEW);
	newtree->Branch("pg1g2_VAJHU_NEW",&pg1g2_VAJHU_NEW);
	newtree->Branch("pg1g2_pi2_VAJHU_NEW",&pg1g2_pi2_VAJHU_NEW);
	newtree->Branch("pg1g4_VAJHU_NEW",&pg1g4_VAJHU_NEW);
	newtree->Branch("pg1g4_pi2_VAJHU_NEW",&pg1g4_pi2_VAJHU_NEW);
	newtree->Branch("pg1g1prime2_VAJHU_NEW",&pg1g1prime2_VAJHU_NEW);
	newtree->Branch("pg1g1prime2_pi2_VAJHU_NEW",&pg1g1prime2_pi2_VAJHU_NEW);


	newtree->Branch("bkg_VAMCFM",&bkg_VAMCFM);
	newtree->Branch("ggzz_VAMCFM",&ggzz_VAMCFM);
	newtree->Branch("p0plus_VAMCFM",&ggHZZ_prob_pure);
	newtree->Branch("ggzz_p0plus_VAMCFM",&ggHZZ_prob_int);
	newtree->Branch("bkg_VAMCFM_NEW",&bkg_VAMCFM_NEW);
	newtree->Branch("ggzz_VAMCFM_NEW",&ggzz_VAMCFM_NEW);
	newtree->Branch("p0plus_VAMCFM_NEW",&ggHZZ_prob_pure_NEW);
	newtree->Branch("ggzz_p0plus_VAMCFM_NEW",&ggHZZ_prob_int_NEW);
	newtree->Branch("p0plus_VAMCFM_NEW_BSMOn",&p0plus_VAMCFM_NEW_BSMOn);

//	newtree->Branch("bkg_VAMCFM_STU",&bkg_VAMCFM_STU);
//	newtree->Branch("bkg_VAMCFM_S",&bkg_VAMCFM_S);
//	newtree->Branch("bkg_VAMCFM_TU",&bkg_VAMCFM_TU);

	int nEntries = tree->GetEntries();
	double selfDHvvcoupl[SIZE_HVV][2] = { { 0 } };
	double ggvvcoupl[2]={0,0};
	mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
//	for (int ev = 0; ev < nEntries; ev++){
//	for (int ev = 0; ev < 100; ev++){
	for (int ev = 0; ev < 1000; ev++){
		tree->GetEntry(ev);

		int lepIdOrdered[4]={ 11,-11,11,-11 };
		if (flavor == 0){
			lepIdOrdered[0]=13;lepIdOrdered[1]=-13;
		};
		float angularOrdered[8] = { mzz, m1, m2, hs, h1, h2, phi, phi1 };

//		mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
//		selfDHvvcoupl[0][0]=1;
//		mela.setMelaLeptonInterference(TVar::InterfOn);
//		p0plus_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB);
		mela.setMelaLeptonInterference(TVar::InterfOn);
		bkg_VAMCFM_NEW = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl,1); // |qqZZ|**2

		mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZGG);
//		mela.setMelaLeptonInterference(TVar::InterfOn);
		ggzz_VAMCFM_NEW = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |ggZZ|**2

		mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
//		mela.setMelaLeptonInterference(TVar::InterfOn);
		ggHZZ_prob_pure_NEW = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |ggHZZ|**2

		mela.setProcess(TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
//		mela.setMelaLeptonInterference(TVar::InterfOn);
		ggZZ_prob_Total_NEW = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |ggZZ + ggHZZ|**2
		ggHZZ_prob_int_NEW = ggZZ_prob_Total_NEW;

		mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[0][0]=1;
		p0plus_VAMCFM_NEW_BSMOn = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);

		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[1][0]=0;
		selfDHvvcoupl[1][1]=0;
		selfDHvvcoupl[3][0]=0;
		selfDHvvcoupl[3][1]=0;
		selfDHvvcoupl[11][0]=0;
		selfDHvvcoupl[11][1]=0;
		p0plus_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		selfDHvvcoupl[0][0]=0;
		selfDHvvcoupl[1][0]=1;
		selfDHvvcoupl[1][1]=0;
		selfDHvvcoupl[3][0]=0;
		selfDHvvcoupl[3][1]=0;
		selfDHvvcoupl[11][0]=0;
		p0hplus_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		selfDHvvcoupl[0][0]=0;
		selfDHvvcoupl[1][0]=0;
		selfDHvvcoupl[1][1]=0;
		selfDHvvcoupl[3][0]=1;
		selfDHvvcoupl[3][1]=0;
		selfDHvvcoupl[11][0]=0;
		p0minus_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		selfDHvvcoupl[0][0]=0;
		selfDHvvcoupl[1][0]=0;
		selfDHvvcoupl[1][1]=0;
		selfDHvvcoupl[3][0]=0;
		selfDHvvcoupl[3][1]=0;
		selfDHvvcoupl[11][0]=-12046.01;
		p0_g1prime2_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[1][0]=1.638;
		selfDHvvcoupl[1][1]=0;
		selfDHvvcoupl[3][0]=0;
		selfDHvvcoupl[3][1]=0;
		selfDHvvcoupl[11][0]=0;
		pg1g2_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		pg1g2_VAJHU_NEW -= (p0plus_VAJHU_NEW + pow(1.638,2)*p0hplus_VAJHU_NEW);

		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[1][0]=0;
		selfDHvvcoupl[1][1]=1.638;
		selfDHvvcoupl[3][0]=0;
		selfDHvvcoupl[3][1]=0;
		selfDHvvcoupl[11][0]=0;
		pg1g2_pi2_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		pg1g2_pi2_VAJHU_NEW -= (p0plus_VAJHU_NEW + pow(1.638,2)*p0hplus_VAJHU_NEW);

		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[1][0]=0;
		selfDHvvcoupl[1][1]=0;
		selfDHvvcoupl[3][0]=2.521;
		selfDHvvcoupl[3][1]=0;
		selfDHvvcoupl[11][0]=0;
		pg1g4_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		pg1g4_VAJHU_NEW -= (p0plus_VAJHU_NEW + pow(2.521,2)*p0minus_VAJHU_NEW);

		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[1][0]=0;
		selfDHvvcoupl[1][1]=0;
		selfDHvvcoupl[3][0]=0;
		selfDHvvcoupl[3][1]=2.521;
		selfDHvvcoupl[11][0]=0;
		pg1g4_pi2_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		pg1g4_pi2_VAJHU_NEW -= (p0plus_VAJHU_NEW + pow(2.521,2)*p0minus_VAJHU_NEW);

		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[1][0]=0;
		selfDHvvcoupl[1][1]=0;
		selfDHvvcoupl[3][0]=0;
		selfDHvvcoupl[3][1]=0;
		selfDHvvcoupl[11][0]=12046.01;
		pg1g1prime2_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		pg1g1prime2_VAJHU_NEW -= (p0plus_VAJHU_NEW + p0_g1prime2_VAJHU_NEW);

		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[1][0]=0;
		selfDHvvcoupl[1][1]=0;
		selfDHvvcoupl[3][0]=0;
		selfDHvvcoupl[3][1]=0;
		selfDHvvcoupl[11][0]=0;
		selfDHvvcoupl[11][1]=12046.01;
		pg1g1prime2_pi2_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		pg1g1prime2_pi2_VAJHU_NEW -= (p0plus_VAJHU_NEW + p0_g1prime2_VAJHU_NEW);


//		mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB_STU);
//		bkg_VAMCFM_STU = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl,1); // |qqZZ|**2

//		mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB_S);
//		bkg_VAMCFM_S = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl,1); // |qqZZ|**2

//		mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB_TU);
//		bkg_VAMCFM_TU = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl,1); // |qqZZ|**2

		newtree->Fill();
	};


	foutput->WriteTObject(newtree);
	foutput->Close();
	finput->Close();
};

void testME_FullMELA_MCFMBSMHiggs(int flavor=1){
	int erg_tev=8;
	float mPOLE=125.6;
	float wPOLE=4.15e-3;
	char TREE_NAME[] = "GenTree";

//	TVar::VerbosityLevel verbosity = TVar::INFO;

	Mela mela(erg_tev,mPOLE);
	TGraph* vaScale_4e = mela.vaScale_4e;
	TGraph* vaScale_4mu = mela.vaScale_4mu;
	TGraph* vaScale_2e2mu = mela.vaScale_2e2mu;
	TGraph* DggZZ_scalefactor = mela.DggZZ_scalefactor;
//  mela.resetMCFM_EWKParameters(1.16639E-05, 7.81751E-03, 79.9549392, 91.1876);
//  mela.resetMCFM_EWKParameters(1.16639E-05, 7.8125E-03, 79.9549391901877726, 91.1876);
//  mela.resetMCFM_EWKParameters(1.16639E-05, 0.007846559134530596, 80.399, 91.1876);
  mela.resetMCFM_EWKParameters(1.16639E-05, 1./128., 80.399, 91.1876, 0.23119);

	TFile* finput = new TFile("/afs/cern.ch/work/u/usarica/HZZ4l-125p6-FullAnalysis/LHC_8TeV/GenSignal/HZZ4lTree_powheg15jhuGenV3-0PMH125.6_Generated.root","read");
	TFile* foutput = new TFile(Form("HZZ4l-125_6-8TeV-Sig_%s_JHUMCFMTest.root",(flavor==0 ? "2mu2e" : "4e")),"recreate");

	float p0plus_VAJHU;
	float p0plus_VAMCFM;
	float p0hplus_VAJHU;
	float p0minus_VAJHU;
	float p0_g1prime2_VAJHU;
	float pg1g1prime2_VAJHU;
	float pg1g1prime2_pi2_VAJHU;
	float pg1g2_VAJHU;
	float pg1g2_pi2_VAJHU;
	float pg1g4_VAJHU;
	float pg1g4_pi2_VAJHU;

	float p0plus_VAJHU_NEW;
	float p0hplus_VAJHU_NEW;
	float p0minus_VAJHU_NEW;
	float p0_g1prime2_VAJHU_NEW;
	float pg1g1prime2_VAJHU_NEW;
	float pg1g1prime2_pi2_VAJHU_NEW;
	float pg1g2_VAJHU_NEW;
	float pg1g2_pi2_VAJHU_NEW;
	float pg1g4_VAJHU_NEW;
	float pg1g4_pi2_VAJHU_NEW;

	float ggHZZ_prob_pure,ggHZZ_prob_int,ggZZ_prob_Total;
	float bkg_VAMCFM,ggzz_VAMCFM;
	float ggHZZ_prob_pure_NEW,ggHZZ_prob_int_NEW,ggZZ_prob_Total_NEW;
	float bkg_VAMCFM_NEW,ggzz_VAMCFM_NEW;
	float bkg_VAMCFM_STU,bkg_VAMCFM_TU,bkg_VAMCFM_S;

	float p0plus_m4l, p0plus_m4l_ScaleUp, p0plus_m4l_ScaleDown, p0plus_m4l_ResUp, p0plus_m4l_ResDown;
	float bkg_m4l, bkg_m4l_ScaleUp, bkg_m4l_ScaleDown, bkg_m4l_ResUp, bkg_m4l_ResDown;

	float mzz = 126.; 
	float m1 = 91.471450;
	float m2 = 12.139782;
	float h1 = 0.2682896;
	float h2 = 0.1679779;
	float phi = 1.5969792;
	float hs = -0.727181;
	float phi1 = 1.8828257;
	float GenLep1Id,GenLep2Id,GenLep3Id,GenLep4Id;

	TTree* tree = (TTree*) finput->Get(TREE_NAME);
	tree->SetBranchAddress("GenZZMass", &mzz);
	tree->SetBranchAddress("GenZ1Mass", &m1);
	tree->SetBranchAddress("GenZ2Mass", &m2);
	tree->SetBranchAddress("GenhelcosthetaZ1", &h1);
	tree->SetBranchAddress("GenhelcosthetaZ2", &h2);
	tree->SetBranchAddress("Genhelphi", &phi);
	tree->SetBranchAddress("Gencosthetastar", &hs);
	tree->SetBranchAddress("GenphistarZ1", &phi1);
	tree->SetBranchAddress("GenLep1Id", &GenLep1Id);
	tree->SetBranchAddress("GenLep2Id", &GenLep2Id);
	tree->SetBranchAddress("GenLep3Id", &GenLep3Id);
	tree->SetBranchAddress("GenLep4Id", &GenLep4Id);

	TTree* newtree = new TTree("TestTree","");
	newtree->Branch("GenZZMass", &mzz);
	newtree->Branch("GenZ1Mass", &m1);
	newtree->Branch("GenZ2Mass", &m2);
	newtree->Branch("GenhelcosthetaZ1", &h1);
	newtree->Branch("GenhelcosthetaZ2", &h2);
	newtree->Branch("Genhelphi", &phi);
	newtree->Branch("Gencosthetastar", &hs);
	newtree->Branch("GenphistarZ1", &phi1);
	newtree->Branch("GenLep1Id", &GenLep1Id);
	newtree->Branch("GenLep2Id", &GenLep2Id);
	newtree->Branch("GenLep3Id", &GenLep3Id);
	newtree->Branch("GenLep4Id", &GenLep4Id);

	newtree->Branch("p0plus_VAJHU",&p0plus_VAJHU);
	newtree->Branch("p0plus_VAMCFM",&p0plus_VAMCFM);
	newtree->Branch("p0hplus_VAJHU",&p0hplus_VAJHU);
	newtree->Branch("p0minus_VAJHU",&p0minus_VAJHU);
	newtree->Branch("p0_g1prime2_VAJHU",&p0_g1prime2_VAJHU);
	newtree->Branch("pg1g2_VAJHU",&pg1g2_VAJHU);
	newtree->Branch("pg1g2_pi2_VAJHU",&pg1g2_pi2_VAJHU);
	newtree->Branch("pg1g4_VAJHU",&pg1g4_VAJHU);
	newtree->Branch("pg1g4_pi2_VAJHU",&pg1g4_pi2_VAJHU);
	newtree->Branch("pg1g1prime2_VAJHU",&pg1g1prime2_VAJHU);
	newtree->Branch("pg1g1prime2_pi2_VAJHU",&pg1g1prime2_pi2_VAJHU);

	newtree->Branch("p0plus_VAJHU_NEW",&p0plus_VAJHU_NEW);
	newtree->Branch("p0hplus_VAJHU_NEW",&p0hplus_VAJHU_NEW);
	newtree->Branch("p0minus_VAJHU_NEW",&p0minus_VAJHU_NEW);
	newtree->Branch("p0_g1prime2_VAJHU_NEW",&p0_g1prime2_VAJHU_NEW);
	newtree->Branch("pg1g2_VAJHU_NEW",&pg1g2_VAJHU_NEW);
	newtree->Branch("pg1g2_pi2_VAJHU_NEW",&pg1g2_pi2_VAJHU_NEW);
	newtree->Branch("pg1g4_VAJHU_NEW",&pg1g4_VAJHU_NEW);
	newtree->Branch("pg1g4_pi2_VAJHU_NEW",&pg1g4_pi2_VAJHU_NEW);
	newtree->Branch("pg1g1prime2_VAJHU_NEW",&pg1g1prime2_VAJHU_NEW);
	newtree->Branch("pg1g1prime2_pi2_VAJHU_NEW",&pg1g1prime2_pi2_VAJHU_NEW);

	int nEntries = tree->GetEntries();
	double selfDHvvcoupl[SIZE_HVV][2] = { { 0 } };
	double ggvvcoupl[2]={0,0};
	mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
	int ev=0,ctr=0;
//	for (int ev = 0; ev < nEntries; ev++){
	while(ctr<3000){
//	for (int ev = 0; ev < 10000; ev++){
		tree->GetEntry(ev);
		ev++;

		int lepIdOrdered[4] = { GenLep1Id, GenLep2Id, GenLep3Id, GenLep4Id };
		float angularOrdered[8] = { mzz, m1, m2, hs, h1, h2, phi, phi1 };

		if(abs(GenLep1Id)==abs(GenLep3Id) && flavor==0) continue;
		if(abs(GenLep1Id)!=abs(GenLep3Id) && flavor==1) continue;
		ctr++;

		mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[0][0]=1;
		mela.setMelaHiggsWidth(wPOLE);
		p0plus_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[1][0]=1.638;
		mela.setMelaHiggsWidth(wPOLE);
		p0hplus_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[3][0]=2.521;
		mela.setMelaHiggsWidth(wPOLE);
		p0minus_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[11][0]=-12046.01;
		mela.setMelaHiggsWidth(wPOLE);
		p0_g1prime2_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[1][0]=1.638;
		mela.setMelaHiggsWidth(wPOLE);
		pg1g2_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		pg1g2_VAJHU -= (p0plus_VAJHU + p0hplus_VAJHU);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[1][1]=1.638;
		mela.setMelaHiggsWidth(wPOLE);
		pg1g2_pi2_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		pg1g2_pi2_VAJHU -= (p0plus_VAJHU + p0hplus_VAJHU);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[3][0]=2.521;
		mela.setMelaHiggsWidth(wPOLE);
		pg1g4_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		pg1g4_VAJHU -= (p0plus_VAJHU + p0minus_VAJHU);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[3][1]=2.521;
		mela.setMelaHiggsWidth(wPOLE);
		pg1g4_pi2_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		pg1g4_pi2_VAJHU -= (p0plus_VAJHU + p0minus_VAJHU);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[11][0]=12046.01;
		mela.setMelaHiggsWidth(wPOLE);
		pg1g1prime2_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		pg1g1prime2_VAJHU -= (p0plus_VAJHU + p0_g1prime2_VAJHU);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[11][1]=12046.01;
		mela.setMelaHiggsWidth(wPOLE);
		pg1g1prime2_pi2_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		pg1g1prime2_pi2_VAJHU -= (p0plus_VAJHU + p0_g1prime2_VAJHU);


		mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		mela.setMelaLeptonInterference(TVar::InterfOn);
		mela.setMelaHiggsWidth(wPOLE);
		p0plus_VAMCFM = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl);
		selfDHvvcoupl[0][0]=1;
		mela.setMelaHiggsWidth(wPOLE);
		p0plus_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[1][0]=1.638;
		mela.setMelaLeptonInterference(TVar::InterfOn);
		mela.setMelaHiggsWidth(wPOLE);
		p0hplus_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[3][0]=2.521;
		mela.setMelaLeptonInterference(TVar::InterfOn);
		mela.setMelaHiggsWidth(wPOLE);
		p0minus_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[11][0]=-12046.01;
		mela.setMelaLeptonInterference(TVar::InterfOn);
		mela.setMelaHiggsWidth(wPOLE);
		p0_g1prime2_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[1][0]=1.638;
		mela.setMelaLeptonInterference(TVar::InterfOn);
		mela.setMelaHiggsWidth(wPOLE);
		pg1g2_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		pg1g2_VAJHU_NEW -= (p0plus_VAJHU_NEW + p0hplus_VAJHU_NEW);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[1][1]=1.638;
		mela.setMelaLeptonInterference(TVar::InterfOn);
		mela.setMelaHiggsWidth(wPOLE);
		pg1g2_pi2_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		pg1g2_pi2_VAJHU_NEW -= (p0plus_VAJHU_NEW + p0hplus_VAJHU_NEW);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[3][0]=2.521;
		mela.setMelaLeptonInterference(TVar::InterfOn);
		mela.setMelaHiggsWidth(wPOLE);
		pg1g4_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		pg1g4_VAJHU_NEW -= (p0plus_VAJHU_NEW + p0minus_VAJHU_NEW);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[3][1]=2.521;
		mela.setMelaLeptonInterference(TVar::InterfOn);
		mela.setMelaHiggsWidth(wPOLE);
		pg1g4_pi2_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		pg1g4_pi2_VAJHU_NEW -= (p0plus_VAJHU_NEW + p0minus_VAJHU_NEW);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[11][0]=12046.01;
		mela.setMelaLeptonInterference(TVar::InterfOn);
		mela.setMelaHiggsWidth(wPOLE);
		pg1g1prime2_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		pg1g1prime2_VAJHU_NEW -= (p0plus_VAJHU_NEW + p0_g1prime2_VAJHU_NEW);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[11][1]=12046.01;
		mela.setMelaLeptonInterference(TVar::InterfOn);
		mela.setMelaHiggsWidth(wPOLE);
		pg1g1prime2_pi2_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		pg1g1prime2_pi2_VAJHU_NEW -= (p0plus_VAJHU_NEW + p0_g1prime2_VAJHU_NEW);

		p0hplus_VAJHU/=p0plus_VAJHU;
		p0minus_VAJHU/=p0plus_VAJHU;
		p0_g1prime2_VAJHU/=p0plus_VAJHU;
		pg1g1prime2_VAJHU/=p0plus_VAJHU;
		pg1g1prime2_pi2_VAJHU/=p0plus_VAJHU;
		pg1g2_VAJHU/=p0plus_VAJHU;
		pg1g2_pi2_VAJHU/=p0plus_VAJHU;
		pg1g4_VAJHU/=p0plus_VAJHU;
		pg1g4_pi2_VAJHU/=p0plus_VAJHU;

		p0hplus_VAJHU_NEW/=p0plus_VAJHU_NEW;
		p0minus_VAJHU_NEW/=p0plus_VAJHU_NEW;
		p0_g1prime2_VAJHU_NEW/=p0plus_VAJHU_NEW;
		pg1g1prime2_VAJHU_NEW/=p0plus_VAJHU_NEW;
		pg1g1prime2_pi2_VAJHU_NEW/=p0plus_VAJHU_NEW;
		pg1g2_VAJHU_NEW/=p0plus_VAJHU_NEW;
		pg1g2_pi2_VAJHU_NEW/=p0plus_VAJHU_NEW;
		pg1g4_VAJHU_NEW/=p0plus_VAJHU_NEW;
		pg1g4_pi2_VAJHU_NEW/=p0plus_VAJHU_NEW;

		p0plus_VAJHU_NEW/=p0plus_VAJHU;
		p0plus_VAMCFM/=p0plus_VAJHU;

		newtree->Fill();

    if (ctr==1){
      cout << ewscheme_.ewscheme << endl;
      cout
        << "vsq: " << ewcouple_.vevsq << '\t'
        << "GF: " << ewcouple_.Gf << '\t'
        << "xw: " << ewcouple_.xw << '\t'
        << "gW: " << ewcouple_.gw << '\t'
        << "mW: " << masses_mcfm_.wmass << '\t'
        << "mZ: " << masses_mcfm_.zmass << '\t'
        << "mT: " << masses_mcfm_.mt
        << endl;
    }
	}


	foutput->WriteTObject(newtree);
	foutput->Close();
	finput->Close();
}


void testDecayModes(){
  int first=1, second=3;
  __modmisc_MOD_getdecaymodes(&first, &second);
  cout << first << '\t' << second << endl;
  first=1; second=3;
  __modmisc_MOD_setdecaymodes(&first, &second);
  __modmisc_MOD_getdecaymodes(&first, &second);
  cout << first << '\t' << second << endl;
}

//
//void testME_FullMELA_BSMMCFMVBF(float multiplierWW=1, int flavor=0){
//  int erg_tev=8;
//  float mPOLE=125.6;
//  float wPOLE=4.15e-3;
//  char TREE_NAME[] = "SelectedTree";
//
//  Mela mela(erg_tev, mPOLE);
//
//  //TString cinput_main = "/scratch0/hep/ianderso/CJLST/140519/PRODFSR_8TeV";
//  //TFile* finput = new TFile(Form("%s/%s/HZZ4lTree_ZZTo%s.root", cinput_main.Data(), (flavor==0 ? "2mu2e" : "4e"), (flavor==0 ? "2e2mu" : "4e")), "read");
//  TString cinput_main = "/scratch0/hep/usarical/VBF-JHUGen-13TeV/LHC_13TeV_VBFHZZ4l/0-";
//  TFile* finput = new TFile(Form("%s/VBFH125toZZ4l_13TeV.root", cinput_main.Data()), "read");
//  TFile* foutput = new TFile(Form("HZZ4lTree_ZZTo%s_vbfMCFMJHUGenTest.root", (flavor==0 ? "2e2mu" : "4e")), "recreate");
//
//  TLorentzVector nullFourVector(0, 0, 0, 0);
//
//  float pvbf0plus_VAMCFM;
//  float pvbf0hplus_VAMCFM;
//  float pvbf0minus_VAMCFM;
//  float pvbf0_g1prime2_VAMCFM;
//  float pvbfg1g1prime2_VAMCFM;
//  float pvbfg1g1prime2_pi2_VAMCFM;
//  float pvbfg1g2_VAMCFM;
//  float pvbfg1g2_pi2_VAMCFM;
//  float pvbfg1g4_VAMCFM;
//  float pvbfg1g4_pi2_VAMCFM;
//
//  float pFullvbf0plus_VAMCFM[11][11]={ { 0 } };
//  float pFullvbf0hplus_VAMCFM[11][11]={ { 0 } };
//  float pFullvbf0minus_VAMCFM[11][11]={ { 0 } };
//  float pFullvbf0_g1prime2_VAMCFM[11][11]={ { 0 } };
//  float pFullvbfg1g1prime2_VAMCFM[11][11]={ { 0 } };
//  float pFullvbfg1g1prime2_pi2_VAMCFM[11][11]={ { 0 } };
//  float pFullvbfg1g2_VAMCFM[11][11]={ { 0 } };
//  float pFullvbfg1g2_pi2_VAMCFM[11][11]={ { 0 } };
//  float pFullvbfg1g4_VAMCFM[11][11]={ { 0 } };
//  float pFullvbfg1g4_pi2_VAMCFM[11][11]={ { 0 } };
//
//  float p0plus_VAMCFM;
//  float p0hplus_VAMCFM;
//  float p0minus_VAMCFM;
//  float p0_g1prime2_VAMCFM;
//  float pg1g1prime2_VAMCFM;
//  float pg1g1prime2_pi2_VAMCFM;
//  float pg1g2_VAMCFM;
//  float pg1g2_pi2_VAMCFM;
//  float pg1g4_VAMCFM;
//  float pg1g4_pi2_VAMCFM;
//
//  float pvbf0plus_VAJHU;
//  float pvbf0hplus_VAJHU;
//  float pvbf0minus_VAJHU;
//  float pvbf0_g1prime2_VAJHU;
//  float pvbfg1g1prime2_VAJHU;
//  float pvbfg1g1prime2_pi2_VAJHU;
//  float pvbfg1g2_VAJHU;
//  float pvbfg1g2_pi2_VAJHU;
//  float pvbfg1g4_VAJHU;
//  float pvbfg1g4_pi2_VAJHU;
//
//  float pFullvbf0plus_VAJHU[11][11]={ { 0 } };
//  float pFullvbf0hplus_VAJHU[11][11]={ { 0 } };
//  float pFullvbf0minus_VAJHU[11][11]={ { 0 } };
//  float pFullvbf0_g1prime2_VAJHU[11][11]={ { 0 } };
//  float pFullvbfg1g1prime2_VAJHU[11][11]={ { 0 } };
//  float pFullvbfg1g1prime2_pi2_VAJHU[11][11]={ { 0 } };
//  float pFullvbfg1g2_VAJHU[11][11]={ { 0 } };
//  float pFullvbfg1g2_pi2_VAJHU[11][11]={ { 0 } };
//  float pFullvbfg1g4_VAJHU[11][11]={ { 0 } };
//  float pFullvbfg1g4_pi2_VAJHU[11][11]={ { 0 } };
//
//  float p0plus_VAJHU;
//  float p0hplus_VAJHU;
//  float p0minus_VAJHU;
//  float p0_g1prime2_VAJHU;
//  float pg1g1prime2_VAJHU;
//  float pg1g1prime2_pi2_VAJHU;
//  float pg1g2_VAJHU;
//  float pg1g2_pi2_VAJHU;
//  float pg1g4_VAJHU;
//  float pg1g4_pi2_VAJHU;
//
//
//  float jet1Pt, jet2Pt;
//  float jet1px, jet1py, jet1pz, jet1E;
//  float jet2px, jet2py, jet2pz, jet2E;
//  float ZZPx, ZZPy, ZZPz, ZZE, dR;
//  short NJets30;
//  std::vector<double> * JetPt=0;
//  std::vector<double> * JetEta=0;
//  std::vector<double> * JetPhi=0;
//  std::vector<double> * JetMass=0;
//  std::vector<int> * JetId=0;
//  TBranch* bJetPt=0;
//  TBranch* bJetEta=0;
//  TBranch* bJetPhi=0;
//  TBranch* bJetMass=0;
//  TBranch* bJetId=0;
//
//  float mzz = 126.;
//  float m1 = 91.471450;
//  float m2 = 12.139782;
//  float h1 = 0.2682896;
//  float h2 = 0.1679779;
//  float phi = 1.5969792;
//  float hs = -0.727181;
//  float phi1 = 1.8828257;
//  int GenLep1Id, GenLep2Id, GenLep3Id, GenLep4Id;
//  float GenLep1Pt, GenLep2Pt, GenLep3Pt, GenLep4Pt;
//  float GenLep1Eta, GenLep2Eta, GenLep3Eta, GenLep4Eta;
//  float GenLep1Phi, GenLep2Phi, GenLep3Phi, GenLep4Phi;
//  float GenLep1Mass, GenLep2Mass, GenLep3Mass, GenLep4Mass;
//  float ZZPt, ZZPhi, ZZEta, DijetMass;
//
//  TTree* tree = (TTree*)finput->Get(TREE_NAME);
//  tree->SetBranchAddress("GenDijetMass", &DijetMass);
//  tree->SetBranchAddress("GenAssociatedParticleId", &JetId, &bJetId);
//  tree->SetBranchAddress("GenAssociatedParticlePt", &JetPt, &bJetPt);
//  tree->SetBranchAddress("GenAssociatedParticleEta", &JetEta, &bJetEta);
//  tree->SetBranchAddress("GenAssociatedParticlePhi", &JetPhi, &bJetPhi);
//  tree->SetBranchAddress("GenAssociatedParticleMass", &JetMass, &bJetMass);
//  tree->SetBranchAddress("GenLep1Id", &GenLep1Id);
//  tree->SetBranchAddress("GenLep1Mass", &GenLep1Mass);
//  tree->SetBranchAddress("GenLep1Pt", &GenLep1Pt);
//  tree->SetBranchAddress("GenLep1Eta", &GenLep1Eta);
//  tree->SetBranchAddress("GenLep1Phi", &GenLep1Phi);
//  tree->SetBranchAddress("GenLep2Id", &GenLep2Id);
//  tree->SetBranchAddress("GenLep2Mass", &GenLep2Mass);
//  tree->SetBranchAddress("GenLep2Pt", &GenLep2Pt);
//  tree->SetBranchAddress("GenLep2Eta", &GenLep2Eta);
//  tree->SetBranchAddress("GenLep2Phi", &GenLep2Phi);
//  tree->SetBranchAddress("GenLep3Id", &GenLep3Id);
//  tree->SetBranchAddress("GenLep3Mass", &GenLep3Mass);
//  tree->SetBranchAddress("GenLep3Pt", &GenLep3Pt);
//  tree->SetBranchAddress("GenLep3Eta", &GenLep3Eta);
//  tree->SetBranchAddress("GenLep3Phi", &GenLep3Phi);
//  tree->SetBranchAddress("GenLep4Id", &GenLep4Id);
//  tree->SetBranchAddress("GenLep4Mass", &GenLep4Mass);
//  tree->SetBranchAddress("GenLep4Pt", &GenLep4Pt);
//  tree->SetBranchAddress("GenLep4Eta", &GenLep4Eta);
//  tree->SetBranchAddress("GenLep4Phi", &GenLep4Phi);
//  tree->SetBranchAddress("GenHMass", &mzz);
//  tree->SetBranchAddress("GenZ1Mass", &m1);
//  tree->SetBranchAddress("GenZ2Mass", &m2);
//  tree->SetBranchAddress("GenhelcosthetaZ1", &h1);
//  tree->SetBranchAddress("GenhelcosthetaZ2", &h2);
//  tree->SetBranchAddress("Genhelphi", &phi);
//  tree->SetBranchAddress("Gencosthetastar", &hs);
//  tree->SetBranchAddress("GenphistarZ1", &phi1);
//
//  TTree* newtree = new TTree("TestTree", "");
//  newtree->Branch("p0plus_VAJHU", &p0plus_VAJHU);
//  newtree->Branch("p0hplus_VAJHU", &p0hplus_VAJHU);
//  newtree->Branch("p0minus_VAJHU", &p0minus_VAJHU);
//  newtree->Branch("p0_g1prime2_VAJHU", &p0_g1prime2_VAJHU);
//  newtree->Branch("pg1g1prime2_VAJHU", &pg1g1prime2_VAJHU);
//  newtree->Branch("pg1g1prime2_pi2_VAJHU", &pg1g1prime2_pi2_VAJHU);
//  newtree->Branch("pg1g2_VAJHU", &pg1g2_VAJHU);
//  newtree->Branch("pg1g2_pi2_VAJHU", &pg1g2_pi2_VAJHU);
//  newtree->Branch("pg1g4_VAJHU", &pg1g4_VAJHU);
//  newtree->Branch("pg1g4_pi2_VAJHU", &pg1g4_pi2_VAJHU);
//  newtree->Branch("p0plus_VAMCFM", &p0plus_VAMCFM);
//  newtree->Branch("p0hplus_VAMCFM", &p0hplus_VAMCFM);
//  newtree->Branch("p0minus_VAMCFM", &p0minus_VAMCFM);
//  newtree->Branch("p0_g1prime2_VAMCFM", &p0_g1prime2_VAMCFM);
//  newtree->Branch("pg1g1prime2_VAMCFM", &pg1g1prime2_VAMCFM);
//  newtree->Branch("pg1g1prime2_pi2_VAMCFM", &pg1g1prime2_pi2_VAMCFM);
//  newtree->Branch("pg1g2_VAMCFM", &pg1g2_VAMCFM);
//  newtree->Branch("pg1g2_pi2_VAMCFM", &pg1g2_pi2_VAMCFM);
//  newtree->Branch("pg1g4_VAMCFM", &pg1g4_VAMCFM);
//  newtree->Branch("pg1g4_pi2_VAMCFM", &pg1g4_pi2_VAMCFM);
//  newtree->Branch("pvbf0plus_VAJHU", &pvbf0plus_VAJHU);
//  newtree->Branch("pvbf0hplus_VAJHU", &pvbf0hplus_VAJHU);
//  newtree->Branch("pvbf0minus_VAJHU", &pvbf0minus_VAJHU);
//  newtree->Branch("pvbf0_g1prime2_VAJHU", &pvbf0_g1prime2_VAJHU);
//  newtree->Branch("pvbfg1g1prime2_VAJHU", &pvbfg1g1prime2_VAJHU);
//  newtree->Branch("pvbfg1g1prime2_pi2_VAJHU", &pvbfg1g1prime2_pi2_VAJHU);
//  newtree->Branch("pvbfg1g2_VAJHU", &pvbfg1g2_VAJHU);
//  newtree->Branch("pvbfg1g2_pi2_VAJHU", &pvbfg1g2_pi2_VAJHU);
//  newtree->Branch("pvbfg1g4_VAJHU", &pvbfg1g4_VAJHU);
//  newtree->Branch("pvbfg1g4_pi2_VAJHU", &pvbfg1g4_pi2_VAJHU);
//  newtree->Branch("pvbf0plus_VAMCFM", &pvbf0plus_VAMCFM);
//  newtree->Branch("pvbf0hplus_VAMCFM", &pvbf0hplus_VAMCFM);
//  newtree->Branch("pvbf0minus_VAMCFM", &pvbf0minus_VAMCFM);
//  newtree->Branch("pvbf0_g1prime2_VAMCFM", &pvbf0_g1prime2_VAMCFM);
//  newtree->Branch("pvbfg1g1prime2_VAMCFM", &pvbfg1g1prime2_VAMCFM);
//  newtree->Branch("pvbfg1g1prime2_pi2_VAMCFM", &pvbfg1g1prime2_pi2_VAMCFM);
//  newtree->Branch("pvbfg1g2_VAMCFM", &pvbfg1g2_VAMCFM);
//  newtree->Branch("pvbfg1g2_pi2_VAMCFM", &pvbfg1g2_pi2_VAMCFM);
//  newtree->Branch("pvbfg1g4_VAMCFM", &pvbfg1g4_VAMCFM);
//  newtree->Branch("pvbfg1g4_pi2_VAMCFM", &pvbfg1g4_pi2_VAMCFM);
//
//  newtree->Branch("pFullvbf0plus_VAJHU", pFullvbf0plus_VAJHU,"pFullvbf0plus_VAJHU[11][11]/F");
//  newtree->Branch("pFullvbf0hplus_VAJHU", pFullvbf0hplus_VAJHU, "pFullvbf0hplus_VAJHU[11][11]/F");
//  newtree->Branch("pFullvbf0minus_VAJHU", pFullvbf0minus_VAJHU, "pFullvbf0minus_VAJHU[11][11]/F");
//  newtree->Branch("pFullvbf0_g1prime2_VAJHU", pFullvbf0_g1prime2_VAJHU, "pFullvbf0_g1prime2_VAJHU[11][11]/F");
//  newtree->Branch("pFullvbfg1g1prime2_VAJHU", pFullvbfg1g1prime2_VAJHU, "pFullvbfg1g1prime2_VAJHU[11][11]/F");
//  newtree->Branch("pFullvbfg1g1prime2_pi2_VAJHU", pFullvbfg1g1prime2_pi2_VAJHU, "pFullvbfg1g1prime2_pi2_VAJHU[11][11]/F");
//  newtree->Branch("pFullvbfg1g2_VAJHU", pFullvbfg1g2_VAJHU, "pFullvbfg1g2_VAJHU[11][11]/F");
//  newtree->Branch("pFullvbfg1g2_pi2_VAJHU", pFullvbfg1g2_pi2_VAJHU, "pFullvbfg1g2_pi2_VAJHU[11][11]/F");
//  newtree->Branch("pFullvbfg1g4_VAJHU", pFullvbfg1g4_VAJHU, "pFullvbfg1g4_VAJHU[11][11]/F");
//  newtree->Branch("pFullvbfg1g4_pi2_VAJHU", pFullvbfg1g4_pi2_VAJHU, "pFullvbfg1g4_pi2_VAJHU[11][11]/F");
//  newtree->Branch("pFullvbf0plus_VAMCFM", pFullvbf0plus_VAMCFM, "pFullvbf0plus_VAMCFM[11][11]/F");
//  newtree->Branch("pFullvbf0hplus_VAMCFM", pFullvbf0hplus_VAMCFM, "pFullvbf0hplus_VAMCFM[11][11]/F");
//  newtree->Branch("pFullvbf0minus_VAMCFM", pFullvbf0minus_VAMCFM, "pFullvbf0minus_VAMCFM[11][11]/F");
//  newtree->Branch("pFullvbf0_g1prime2_VAMCFM", pFullvbf0_g1prime2_VAMCFM, "pFullvbf0_g1prime2_VAMCFM[11][11]/F");
//  newtree->Branch("pFullvbfg1g1prime2_VAMCFM", pFullvbfg1g1prime2_VAMCFM, "pFullvbfg1g1prime2_VAMCFM[11][11]/F");
//  newtree->Branch("pFullvbfg1g1prime2_pi2_VAMCFM", pFullvbfg1g1prime2_pi2_VAMCFM, "pFullvbfg1g1prime2_pi2_VAMCFM[11][11]/F");
//  newtree->Branch("pFullvbfg1g2_VAMCFM", pFullvbfg1g2_VAMCFM, "pFullvbfg1g2_VAMCFM[11][11]/F");
//  newtree->Branch("pFullvbfg1g2_pi2_VAMCFM", pFullvbfg1g2_pi2_VAMCFM, "pFullvbfg1g2_pi2_VAMCFM[11][11]/F");
//  newtree->Branch("pFullvbfg1g4_VAMCFM", pFullvbfg1g4_VAMCFM, "pFullvbfg1g4_VAMCFM[11][11]/F");
//  newtree->Branch("pFullvbfg1g4_pi2_VAMCFM", pFullvbfg1g4_pi2_VAMCFM, "pFullvbfg1g4_pi2_VAMCFM[11][11]/F");
//
//  newtree->Branch("ZZMass", &mzz);
//  newtree->Branch("DijetMass", &DijetMass);
//
//  int nEntries = tree->GetEntries();
//  double selfDHggcoupl[SIZE_HGG][2] ={ { 0 } };
//  double selfDHvvcoupl[SIZE_HVV_VBF][2] ={ { 0 } };
//  double selfDHwwcoupl[SIZE_HWW_VBF][2] ={ { 0 } };
//  double selfDTotalHvvcoupl[SIZE_HVV][2] ={ { 0 } };
//  double selfDTotalHwwcoupl[SIZE_HVV][2] ={ { 0 } };
//  double ggvvcoupl[2]={ 0, 0 };
//  mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
//  int recorded=0;
//  for (int ev = 0; ev < nEntries; ev++){
//    if (recorded>=1000) break;
//    tree->GetEntry(ev);
//    
//    /*
//    if (flavor == 1){
//      GenLep1Id=11;
//      GenLep2Id=-11;
//      GenLep3Id=11;
//      GenLep4Id=-11;
//    }
//    else{
//      GenLep1Id=13;
//      GenLep2Id=-13;
//      GenLep3Id=11;
//      GenLep4Id=-11;
//    }
//    TLorentzVector jet1(0, 0, 1e-3, 1e-3), jet2(0, 0, 1e-3, 1e-3), higgs(0, 0, 0, 0);
//    higgs.SetPtEtaPhiM(ZZPt, ZZEta, ZZPhi, mzz);
//
//    int lepIdOrdered[4]={ GenLep1Id, GenLep2Id, GenLep3Id, GenLep4Id };
//    float angularOrdered[8] ={ mzz, m1, m2, hs, h1, h2, phi, phi1 };
//    std::vector<TLorentzVector> pZZ4L = mela.calculate4Momentum((double)mzz, (double)m1, (double)m2, (double)acos(hs), (double)acos(h1), (double)acos(h2), (double)phi1, (double)phi);
//    TLorentzVector higgsCM = pZZ4L[0] + pZZ4L[1] + pZZ4L[2] + pZZ4L[3]; // ...because...
//    TLorentzVector higgsCM_boosted = higgsCM;
//    higgsCM_boosted.Boost(-higgsCM.BoostVector());
//    higgsCM_boosted.Boost(higgs.BoostVector());
//    for (int k=0; k<4; k++){
//      pZZ4L[k].Boost(-higgsCM.BoostVector());
//      pZZ4L[k].Boost(higgsCM_boosted.BoostVector());
//    }
//    TLorentzVector pZZ4LOrdered[4] ={ pZZ4L[0], pZZ4L[1], pZZ4L[2], pZZ4L[3] };
//    higgs.SetXYZT(higgsCM_boosted.X(), higgsCM_boosted.Y(), higgsCM_boosted.Z(), higgsCM_boosted.T());
//    */
//
//    if (flavor==1 && GenLep1Id!=GenLep3Id) continue;
//    else if (GenLep1Id==GenLep3Id) continue;
//
//    int lepIdOrdered[4]={ GenLep1Id, GenLep2Id, GenLep3Id, GenLep4Id };
//    float angularOrdered[8] ={ mzz, m1, m2, hs, h1, h2, phi, phi1 };
//
//    TLorentzVector pZZ4LOrdered[4];
//    TLorentzVector pZZ4LOrdered_CM[4];
//    pZZ4LOrdered[0].SetPtEtaPhiM(GenLep1Pt, GenLep1Eta, GenLep1Phi, GenLep1Mass);
//    pZZ4LOrdered[1].SetPtEtaPhiM(GenLep2Pt, GenLep2Eta, GenLep2Phi, GenLep2Mass);
//    pZZ4LOrdered[2].SetPtEtaPhiM(GenLep3Pt, GenLep3Eta, GenLep3Phi, GenLep3Mass);
//    pZZ4LOrdered[3].SetPtEtaPhiM(GenLep4Pt, GenLep4Eta, GenLep4Phi, GenLep4Mass);
//
//    TLorentzVector jet1(0, 0, 0, 0), jet2(0, 0, 0, 0), higgs(0, 0, 0, 0);
//    for (int i=0; i<4; i++) higgs = higgs + pZZ4LOrdered[i];
//
//    pZZ4LOrdered_CM[0].SetPtEtaPhiM(GenLep1Pt, GenLep1Eta, GenLep1Phi, GenLep1Mass);
//    pZZ4LOrdered_CM[1].SetPtEtaPhiM(GenLep2Pt, GenLep2Eta, GenLep2Phi, GenLep2Mass);
//    pZZ4LOrdered_CM[2].SetPtEtaPhiM(GenLep3Pt, GenLep3Eta, GenLep3Phi, GenLep3Mass);
//    pZZ4LOrdered_CM[3].SetPtEtaPhiM(GenLep4Pt, GenLep4Eta, GenLep4Phi, GenLep4Mass);
//    for (int i=0; i<4; i++) pZZ4LOrdered_CM[i].Boost(-higgs.BoostVector());
//    for (int i=0; i<4; i++) {
//      //cout << "pZZ4LOrdered[" << i << "]:\t" << pZZ4LOrdered[i].X() << '\t' << pZZ4LOrdered[i].Y() << '\t' << pZZ4LOrdered[i].Z() << '\t' << pZZ4LOrdered[i].T() << '\t' << pZZ4LOrdered[i].M() << '\n';
//    }
//    for (int i=0; i<4; i++) {
//      //cout << "pZZ4LOrdered_CM[" << i << "]:\t" << pZZ4LOrdered_CM[i].X() << '\t' << pZZ4LOrdered_CM[i].Y() << '\t' << pZZ4LOrdered_CM[i].Z() << '\t' << pZZ4LOrdered_CM[i].T() << '\t' << pZZ4LOrdered_CM[i].M() << '\n';
//    }
//
//    //if (JetPt->size()>=2 && NJets30>=2 && DijetMass>150){
//    if (DijetMass>150){
//      //cout << higgs.M() << '\t' << (pZZ4LOrdered[0]+pZZ4LOrdered[1]+pZZ4LOrdered[2]+pZZ4LOrdered[3]).M() << endl;
//
//      jet1.SetPtEtaPhiM(JetPt->at(0), JetEta->at(0), JetPhi->at(0), JetMass->at(0));
//      jet2.SetPtEtaPhiM(JetPt->at(1), JetEta->at(1), JetPhi->at(1), JetMass->at(1));
//      TLorentzVector myjets[2]={ jet1, jet2 };
//      //int myjetids[2]={ 0 };
//      int myjetids[2]={ JetId->at(0), JetId->at(1) };
//
//      for (int i=0; i<2; i++) {
//        //cout << "myjets[" << i << "]:\t" << myjets[i].X() << '\t' << myjets[i].Y() << '\t' << myjets[i].Z() << '\t' << myjets[i].T() << '\t' << myjets[i].M() << '\n';
//      }
//
//      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
//
//      for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDTotalHvvcoupl[ii][jj] = 0; }
//      selfDTotalHvvcoupl[0][0]=1;
//      mela.setMelaHiggsWidth(wPOLE);
//      p0plus_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDTotalHvvcoupl);
///*
//      if (recorded==0){
//        mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
//        mela.setMelaHiggsWidth(wPOLE);
//        float tmpProb=0;
//        mela.computeP(
//          pZZ4LOrdered_CM[0], lepIdOrdered[0],
//          pZZ4LOrdered_CM[1], lepIdOrdered[1],
//          pZZ4LOrdered_CM[2], lepIdOrdered[2],
//          pZZ4LOrdered_CM[3], lepIdOrdered[3],
//          tmpProb
//          );
//        cout << higgs.Z() << '\t' << p0plus_VAJHU << '\t' << tmpProb << '\t' << p0plus_VAJHU/tmpProb << endl;
//        mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
//      }
//*/
///*
//      for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDTotalHvvcoupl[ii][jj] = 0; }
//      selfDTotalHvvcoupl[1][0]=1.638;
//      mela.setMelaHiggsWidth(wPOLE);
//      p0hplus_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDTotalHvvcoupl);
//
//      for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDTotalHvvcoupl[ii][jj] = 0; }
//      selfDTotalHvvcoupl[3][0]=2.521;
//      mela.setMelaHiggsWidth(wPOLE);
//      p0minus_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDTotalHvvcoupl);
//*/
//      for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDTotalHvvcoupl[ii][jj] = 0; }
//      selfDTotalHvvcoupl[11][0]=-12046.01;
//      mela.setMelaHiggsWidth(wPOLE);
//      p0_g1prime2_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDTotalHvvcoupl);
//
///*
//      for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDTotalHvvcoupl[ii][jj] = 0; }
//      selfDTotalHvvcoupl[0][0]=1;
//      selfDTotalHvvcoupl[1][0]=1.638;
//      mela.setMelaHiggsWidth(wPOLE);
//      pg1g2_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDTotalHvvcoupl);
//
//      for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDTotalHvvcoupl[ii][jj] = 0; }
//      selfDTotalHvvcoupl[0][0]=1;
//      selfDTotalHvvcoupl[1][1]=1.638;
//      mela.setMelaHiggsWidth(wPOLE);
//      pg1g2_pi2_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDTotalHvvcoupl);
//
//      for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDTotalHvvcoupl[ii][jj] = 0; }
//      selfDTotalHvvcoupl[0][0]=1;
//      selfDTotalHvvcoupl[3][0]=2.521;
//      mela.setMelaHiggsWidth(wPOLE);
//      pg1g4_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDTotalHvvcoupl);
//
//      for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDTotalHvvcoupl[ii][jj] = 0; }
//      selfDTotalHvvcoupl[0][0]=1;
//      selfDTotalHvvcoupl[3][1]=2.521;
//      mela.setMelaHiggsWidth(wPOLE);
//      pg1g4_pi2_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDTotalHvvcoupl);
//*/
//
//      for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDTotalHvvcoupl[ii][jj] = 0; }
//      selfDTotalHvvcoupl[0][0]=1;
//      selfDTotalHvvcoupl[11][0]=-12046.01;
//      mela.setMelaHiggsWidth(wPOLE);
//      pg1g1prime2_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDTotalHvvcoupl);
//
//      for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDTotalHvvcoupl[ii][jj] = 0; }
//      selfDTotalHvvcoupl[0][0]=1;
//      selfDTotalHvvcoupl[11][1]=-12046.01;
//      mela.setMelaHiggsWidth(wPOLE);
//      pg1g1prime2_pi2_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDTotalHvvcoupl);
//
//
//// VBF JHUGen
//      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::JJVBF);
//
//      for (int xx = 0; xx < SIZE_HVV_VBF; xx++){ for (int yy = 0; yy < 2; yy++) selfDHvvcoupl[xx][yy] = 0; }
//      for (int xx = 0; xx < SIZE_HWW_VBF; xx++){ for (int yy = 0; yy < 2; yy++) selfDHwwcoupl[xx][yy] = 0; }
//      selfDHvvcoupl[0][0]=1;
//      selfDHwwcoupl[0][0]=1.*multiplierWW;
//      mela.computeProdP(jet1, 0, jet2, 0, higgs, 25, nullFourVector, 0, selfDHggcoupl, selfDHvvcoupl, selfDHwwcoupl, pvbf0plus_VAJHU);
//      pvbf0plus_VAJHU *= p0plus_VAJHU;
//      me_record rcd_pvbf0plus_VAJHU = mela.getFullMERecord();
//      for (int xx=0; xx<11; xx++){
//        for (int yy=0; yy<11; yy++){
//          (rcd_pvbf0plus_VAJHU.MEsq)[xx][yy] *= p0plus_VAJHU;
//          (rcd_pvbf0plus_VAJHU.weightedMEsq)[xx][yy] *= p0plus_VAJHU;
//        }
//      }
//
///*
//      for (int xx = 0; xx < SIZE_HVV_VBF; xx++){ for (int yy = 0; yy < 2; yy++) selfDHvvcoupl[xx][yy] = 0; }
//      for (int xx = 0; xx < SIZE_HWW_VBF; xx++){ for (int yy = 0; yy < 2; yy++) selfDHwwcoupl[xx][yy] = 0; }
//      selfDHvvcoupl[1][0]=1.638;
//      selfDHwwcoupl[1][0]=1.638*multiplierWW;
//      mela.computeProdP(jet1, 0, jet2, 0, higgs, 25, nullFourVector, 0, selfDHggcoupl, selfDHvvcoupl, selfDHwwcoupl, pvbf0hplus_VAJHU);
//      pvbf0hplus_VAJHU *= p0hplus_VAJHU;
//
//      for (int xx = 0; xx < SIZE_HVV_VBF; xx++){ for (int yy = 0; yy < 2; yy++) selfDHvvcoupl[xx][yy] = 0; }
//      for (int xx = 0; xx < SIZE_HWW_VBF; xx++){ for (int yy = 0; yy < 2; yy++) selfDHwwcoupl[xx][yy] = 0; }
//      selfDHvvcoupl[3][0]=2.521;
//      selfDHwwcoupl[3][0]=2.521*multiplierWW;
//      mela.computeProdP(jet1, 0, jet2, 0, higgs, 25, nullFourVector, 0, selfDHggcoupl, selfDHvvcoupl, selfDHwwcoupl, pvbf0minus_VAJHU);
//      pvbf0minus_VAJHU *= p0minus_VAJHU;
//
//*/
//      for (int xx = 0; xx < SIZE_HVV_VBF; xx++){ for (int yy = 0; yy < 2; yy++) selfDHvvcoupl[xx][yy] = 0; }
//      for (int xx = 0; xx < SIZE_HWW_VBF; xx++){ for (int yy = 0; yy < 2; yy++) selfDHwwcoupl[xx][yy] = 0; }
//      selfDHvvcoupl[11][0]=-12046.01;
//      selfDHwwcoupl[11][0]=-12046.01*multiplierWW;
//      mela.computeProdP(jet1, 0, jet2, 0, higgs, 25, nullFourVector, 0, selfDHggcoupl, selfDHvvcoupl, selfDHwwcoupl, pvbf0_g1prime2_VAJHU);
//      pvbf0_g1prime2_VAJHU *= p0_g1prime2_VAJHU;
//      me_record rcd_pvbf0_g1prime2_VAJHU = mela.getFullMERecord();
//      for (int xx=0; xx<11; xx++){
//        for (int yy=0; yy<11; yy++){
//          (rcd_pvbf0_g1prime2_VAJHU.MEsq)[xx][yy] *= p0_g1prime2_VAJHU;
//          (rcd_pvbf0_g1prime2_VAJHU.weightedMEsq)[xx][yy] *= p0_g1prime2_VAJHU;
//        }
//      }
//
//
///*
//      for (int xx = 0; xx < SIZE_HVV_VBF; xx++){ for (int yy = 0; yy < 2; yy++) selfDHvvcoupl[xx][yy] = 0; }
//      for (int xx = 0; xx < SIZE_HWW_VBF; xx++){ for (int yy = 0; yy < 2; yy++) selfDHwwcoupl[xx][yy] = 0; }
//      selfDHvvcoupl[0][0]=1;
//      selfDHwwcoupl[0][0]=1.*multiplierWW;
//      selfDHvvcoupl[1][0]=1.638;
//      selfDHwwcoupl[1][0]=1.638*multiplierWW;
//      mela.computeProdP(jet1, 0, jet2, 0, higgs, 25, nullFourVector, 0, selfDHggcoupl, selfDHvvcoupl, selfDHwwcoupl, pvbfg1g2_VAJHU);
//      pvbfg1g2_VAJHU *= pg1g2_VAJHU;
//
//      for (int xx = 0; xx < SIZE_HVV_VBF; xx++){ for (int yy = 0; yy < 2; yy++) selfDHvvcoupl[xx][yy] = 0; }
//      for (int xx = 0; xx < SIZE_HWW_VBF; xx++){ for (int yy = 0; yy < 2; yy++) selfDHwwcoupl[xx][yy] = 0; }
//      selfDHvvcoupl[0][0]=1;
//      selfDHwwcoupl[0][0]=1.*multiplierWW;
//      selfDHvvcoupl[1][1]=1.638;
//      selfDHwwcoupl[1][1]=1.638*multiplierWW;
//      mela.computeProdP(jet1, 0, jet2, 0, higgs, 25, nullFourVector, 0, selfDHggcoupl, selfDHvvcoupl, selfDHwwcoupl, pvbfg1g2_pi2_VAJHU);
//      pvbfg1g2_pi2_VAJHU *= pg1g2_pi2_VAJHU;
//
//      for (int xx = 0; xx < SIZE_HVV_VBF; xx++){ for (int yy = 0; yy < 2; yy++) selfDHvvcoupl[xx][yy] = 0; }
//      for (int xx = 0; xx < SIZE_HWW_VBF; xx++){ for (int yy = 0; yy < 2; yy++) selfDHwwcoupl[xx][yy] = 0; }
//      selfDHvvcoupl[0][0]=1;
//      selfDHwwcoupl[0][0]=1.*multiplierWW;
//      selfDHvvcoupl[3][0]=2.521;
//      selfDHwwcoupl[3][0]=2.521*multiplierWW;
//      mela.computeProdP(jet1, 0, jet2, 0, higgs, 25, nullFourVector, 0, selfDHggcoupl, selfDHvvcoupl, selfDHwwcoupl, pvbfg1g4_VAJHU);
//      pvbfg1g4_VAJHU *= pg1g4_VAJHU;
//
//      for (int xx = 0; xx < SIZE_HVV_VBF; xx++){ for (int yy = 0; yy < 2; yy++) selfDHvvcoupl[xx][yy] = 0; }
//      for (int xx = 0; xx < SIZE_HWW_VBF; xx++){ for (int yy = 0; yy < 2; yy++) selfDHwwcoupl[xx][yy] = 0; }
//      selfDHvvcoupl[0][0]=1;
//      selfDHwwcoupl[0][0]=1.*multiplierWW;
//      selfDHvvcoupl[3][1]=2.521;
//      selfDHwwcoupl[3][1]=2.521*multiplierWW;
//      mela.computeProdP(jet1, 0, jet2, 0, higgs, 25, nullFourVector, 0, selfDHggcoupl, selfDHvvcoupl, selfDHwwcoupl, pvbfg1g4_pi2_VAJHU);
//      pvbfg1g4_pi2_VAJHU *= pg1g4_pi2_VAJHU;
//
//*/
//      for (int xx = 0; xx < SIZE_HVV_VBF; xx++){ for (int yy = 0; yy < 2; yy++) selfDHvvcoupl[xx][yy] = 0; }
//      for (int xx = 0; xx < SIZE_HWW_VBF; xx++){ for (int yy = 0; yy < 2; yy++) selfDHwwcoupl[xx][yy] = 0; }
//      selfDHvvcoupl[0][0]=1;
//      selfDHwwcoupl[0][0]=1.*multiplierWW;
//      selfDHvvcoupl[11][0]=-12046.01;
//      selfDHwwcoupl[11][0]=-12046.01*multiplierWW;
//      mela.computeProdP(jet1, 0, jet2, 0, higgs, 25, nullFourVector, 0, selfDHggcoupl, selfDHvvcoupl, selfDHwwcoupl, pvbfg1g1prime2_VAJHU);
//      pvbfg1g1prime2_VAJHU *= pg1g1prime2_VAJHU;
//      me_record rcd_pvbfg1g1prime2_VAJHU = mela.getFullMERecord();
//      for (int xx=0; xx<11; xx++){
//        for (int yy=0; yy<11; yy++){
//          (rcd_pvbfg1g1prime2_VAJHU.MEsq)[xx][yy] *= pg1g1prime2_VAJHU;
//          (rcd_pvbfg1g1prime2_VAJHU.weightedMEsq)[xx][yy] *= pg1g1prime2_VAJHU;
//        }
//      }
//
//      for (int xx = 0; xx < SIZE_HVV_VBF; xx++){ for (int yy = 0; yy < 2; yy++) selfDHvvcoupl[xx][yy] = 0; }
//      for (int xx = 0; xx < SIZE_HWW_VBF; xx++){ for (int yy = 0; yy < 2; yy++) selfDHwwcoupl[xx][yy] = 0; }
//      selfDHvvcoupl[0][0]=1;
//      selfDHwwcoupl[0][0]=1.*multiplierWW;
//      selfDHvvcoupl[11][1]=-12046.01;
//      selfDHwwcoupl[11][1]=-12046.01*multiplierWW;
//      mela.computeProdP(jet1, 0, jet2, 0, higgs, 25, nullFourVector, 0, selfDHggcoupl, selfDHvvcoupl, selfDHwwcoupl, pvbfg1g1prime2_pi2_VAJHU);
//      pvbfg1g1prime2_pi2_VAJHU *= pg1g1prime2_pi2_VAJHU;
//      me_record rcd_pvbfg1g1prime2_pi2_VAJHU = mela.getFullMERecord();
//      for (int xx=0; xx<11; xx++){
//        for (int yy=0; yy<11; yy++){
//          (rcd_pvbfg1g1prime2_pi2_VAJHU.MEsq)[xx][yy] *= pg1g1prime2_pi2_VAJHU;
//          (rcd_pvbfg1g1prime2_pi2_VAJHU.weightedMEsq)[xx][yy] *= pg1g1prime2_pi2_VAJHU;
//        }
//      }
//
//
//// VBF MCFM
//      mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::JJVBF);
//
//      if (recorded==0) cout << "SM MCFM" << endl;
//      for (int xx = 0; xx < SIZE_HVV; xx++){ for (int yy = 0; yy < 2; yy++) selfDTotalHvvcoupl[xx][yy] = 0; }
//      for (int xx = 0; xx < SIZE_HVV; xx++){ for (int yy = 0; yy < 2; yy++) selfDTotalHwwcoupl[xx][yy] = 0; }
//      selfDTotalHvvcoupl[0][0]=1;
//      selfDTotalHwwcoupl[0][0]=1.*multiplierWW;
//      mela.computeProdDecP(myjets, pZZ4LOrdered, myjetids, lepIdOrdered, selfDTotalHvvcoupl, selfDTotalHwwcoupl, pvbf0plus_VAMCFM);
//      me_record rcd_pvbf0plus_VAMCFM = mela.getFullMERecord();
//
///*
//      if (recorded==0) cout << "0h+ MCFM" << endl;
//      for (int xx = 0; xx < SIZE_HVV; xx++){ for (int yy = 0; yy < 2; yy++) selfDTotalHvvcoupl[xx][yy] = 0; }
//      for (int xx = 0; xx < SIZE_HVV; xx++){ for (int yy = 0; yy < 2; yy++) selfDTotalHwwcoupl[xx][yy] = 0; }
//      selfDTotalHvvcoupl[1][0]=1.638;
//      selfDTotalHwwcoupl[1][0]=1.638*multiplierWW;
//      mela.computeProdDecP(myjets, pZZ4LOrdered, myjetids, lepIdOrdered, selfDTotalHvvcoupl, selfDTotalHwwcoupl, pvbf0hplus_VAMCFM);
//
//      for (int xx = 0; xx < SIZE_HVV; xx++){ for (int yy = 0; yy < 2; yy++) selfDTotalHvvcoupl[xx][yy] = 0; }
//      for (int xx = 0; xx < SIZE_HVV; xx++){ for (int yy = 0; yy < 2; yy++) selfDTotalHwwcoupl[xx][yy] = 0; }
//      selfDTotalHvvcoupl[3][0]=2.521;
//      selfDTotalHwwcoupl[3][0]=2.521*multiplierWW;
//      mela.computeProdDecP(myjets, pZZ4LOrdered, myjetids, lepIdOrdered, selfDTotalHvvcoupl, selfDTotalHwwcoupl, pvbf0minus_VAMCFM);
//*/
//
//      for (int xx = 0; xx < SIZE_HVV; xx++){ for (int yy = 0; yy < 2; yy++) selfDTotalHvvcoupl[xx][yy] = 0; }
//      for (int xx = 0; xx < SIZE_HVV; xx++){ for (int yy = 0; yy < 2; yy++) selfDTotalHwwcoupl[xx][yy] = 0; }
//      selfDTotalHvvcoupl[11][0]=-12046.01;
//      selfDTotalHwwcoupl[11][0]=-12046.01*multiplierWW;
//      mela.computeProdDecP(myjets, pZZ4LOrdered, myjetids, lepIdOrdered, selfDTotalHvvcoupl, selfDTotalHwwcoupl, pvbf0_g1prime2_VAMCFM);
//      me_record rcd_pvbf0_g1prime2_VAMCFM = mela.getFullMERecord();
//
///*
//      for (int xx = 0; xx < SIZE_HVV; xx++){ for (int yy = 0; yy < 2; yy++) selfDTotalHvvcoupl[xx][yy] = 0; }
//      for (int xx = 0; xx < SIZE_HVV; xx++){ for (int yy = 0; yy < 2; yy++) selfDTotalHwwcoupl[xx][yy] = 0; }
//      selfDTotalHvvcoupl[0][0]=1;
//      selfDTotalHwwcoupl[0][0]=1.*multiplierWW;
//      selfDTotalHvvcoupl[1][0]=1.638;
//      selfDTotalHwwcoupl[1][0]=1.638*multiplierWW;
//      mela.computeProdDecP(myjets, pZZ4LOrdered, myjetids, lepIdOrdered, selfDTotalHvvcoupl, selfDTotalHwwcoupl, pvbfg1g2_VAMCFM);
//
//      for (int xx = 0; xx < SIZE_HVV; xx++){ for (int yy = 0; yy < 2; yy++) selfDTotalHvvcoupl[xx][yy] = 0; }
//      for (int xx = 0; xx < SIZE_HVV; xx++){ for (int yy = 0; yy < 2; yy++) selfDTotalHwwcoupl[xx][yy] = 0; }
//      selfDTotalHvvcoupl[0][0]=1;
//      selfDTotalHwwcoupl[0][0]=1.*multiplierWW;
//      selfDTotalHvvcoupl[1][1]=1.638;
//      selfDTotalHwwcoupl[1][1]=1.638*multiplierWW;
//      mela.computeProdDecP(myjets, pZZ4LOrdered, myjetids, lepIdOrdered, selfDTotalHvvcoupl, selfDTotalHwwcoupl, pvbfg1g2_pi2_VAMCFM);
//
//      for (int xx = 0; xx < SIZE_HVV; xx++){ for (int yy = 0; yy < 2; yy++) selfDTotalHvvcoupl[xx][yy] = 0; }
//      for (int xx = 0; xx < SIZE_HVV; xx++){ for (int yy = 0; yy < 2; yy++) selfDTotalHwwcoupl[xx][yy] = 0; }
//      selfDTotalHvvcoupl[0][0]=1;
//      selfDTotalHwwcoupl[0][0]=1.*multiplierWW;
//      selfDTotalHvvcoupl[3][0]=2.521;
//      selfDTotalHwwcoupl[3][0]=2.521*multiplierWW;
//      mela.computeProdDecP(myjets, pZZ4LOrdered, myjetids, lepIdOrdered, selfDTotalHvvcoupl, selfDTotalHwwcoupl, pvbfg1g4_VAMCFM);
//
//      for (int xx = 0; xx < SIZE_HVV; xx++){ for (int yy = 0; yy < 2; yy++) selfDTotalHvvcoupl[xx][yy] = 0; }
//      for (int xx = 0; xx < SIZE_HVV; xx++){ for (int yy = 0; yy < 2; yy++) selfDTotalHwwcoupl[xx][yy] = 0; }
//      selfDTotalHvvcoupl[0][0]=1;
//      selfDTotalHwwcoupl[0][0]=1.*multiplierWW;
//      selfDTotalHvvcoupl[3][1]=2.521;
//      selfDTotalHwwcoupl[3][1]=2.521*multiplierWW;
//      mela.computeProdDecP(myjets, pZZ4LOrdered, myjetids, lepIdOrdered, selfDTotalHvvcoupl, selfDTotalHwwcoupl, pvbfg1g4_pi2_VAMCFM);
//
//*/
//      for (int xx = 0; xx < SIZE_HVV; xx++){ for (int yy = 0; yy < 2; yy++) selfDTotalHvvcoupl[xx][yy] = 0; }
//      for (int xx = 0; xx < SIZE_HVV; xx++){ for (int yy = 0; yy < 2; yy++) selfDTotalHwwcoupl[xx][yy] = 0; }
//      selfDTotalHvvcoupl[0][0]=1;
//      selfDTotalHwwcoupl[0][0]=1.*multiplierWW;
//      selfDTotalHvvcoupl[11][0]=-12046.01;
//      selfDTotalHwwcoupl[11][0]=-12046.01*multiplierWW;
//      mela.computeProdDecP(myjets, pZZ4LOrdered, myjetids, lepIdOrdered, selfDTotalHvvcoupl, selfDTotalHwwcoupl, pvbfg1g1prime2_VAMCFM);
//      me_record rcd_pvbfg1g1prime2_VAMCFM = mela.getFullMERecord();
//
//      for (int xx = 0; xx < SIZE_HVV; xx++){ for (int yy = 0; yy < 2; yy++) selfDTotalHvvcoupl[xx][yy] = 0; }
//      for (int xx = 0; xx < SIZE_HVV; xx++){ for (int yy = 0; yy < 2; yy++) selfDTotalHwwcoupl[xx][yy] = 0; }
//      selfDTotalHvvcoupl[0][0]=1;
//      selfDTotalHwwcoupl[0][0]=1.*multiplierWW;
//      selfDTotalHvvcoupl[11][1]=-12046.01;
//      selfDTotalHwwcoupl[11][1]=-12046.01*multiplierWW;
//      mela.computeProdDecP(myjets, pZZ4LOrdered, myjetids, lepIdOrdered, selfDTotalHvvcoupl, selfDTotalHwwcoupl, pvbfg1g1prime2_pi2_VAMCFM);
//      me_record rcd_pvbfg1g1prime2_pi2_VAMCFM = mela.getFullMERecord();
//
//
//      mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
//
//      for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDTotalHvvcoupl[ii][jj] = 0; }
//      selfDTotalHvvcoupl[0][0]=1;
//      mela.setMelaHiggsWidth(wPOLE);
//      p0plus_VAMCFM = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDTotalHvvcoupl);
//
///*
//      for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDTotalHvvcoupl[ii][jj] = 0; }
//      selfDTotalHvvcoupl[1][0]=1.638;
//      mela.setMelaHiggsWidth(wPOLE);
//      p0hplus_VAMCFM = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDTotalHvvcoupl);
//
//
//      for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDTotalHvvcoupl[ii][jj] = 0; }
//      selfDTotalHvvcoupl[3][0]=2.521;
//      mela.setMelaHiggsWidth(wPOLE);
//      p0minus_VAMCFM = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDTotalHvvcoupl);
//*/
//
//      for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDTotalHvvcoupl[ii][jj] = 0; }
//      selfDTotalHvvcoupl[11][0]=-12046.01;
//      mela.setMelaHiggsWidth(wPOLE);
//      p0_g1prime2_VAMCFM = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDTotalHvvcoupl);
//
//
///*
//      for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDTotalHvvcoupl[ii][jj] = 0; }
//      selfDTotalHvvcoupl[0][0]=1;
//      selfDTotalHvvcoupl[1][0]=1.638;
//      mela.setMelaHiggsWidth(wPOLE);
//      pg1g2_VAMCFM = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDTotalHvvcoupl);
//
//      for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDTotalHvvcoupl[ii][jj] = 0; }
//      selfDTotalHvvcoupl[0][0]=1;
//      selfDTotalHvvcoupl[1][1]=1.638;
//      mela.setMelaHiggsWidth(wPOLE);
//      pg1g2_pi2_VAMCFM = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDTotalHvvcoupl);
//
//      for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDTotalHvvcoupl[ii][jj] = 0; }
//      selfDTotalHvvcoupl[0][0]=1;
//      selfDTotalHvvcoupl[3][0]=2.521;
//      mela.setMelaHiggsWidth(wPOLE);
//      pg1g4_VAMCFM = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDTotalHvvcoupl);
//
//      for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDTotalHvvcoupl[ii][jj] = 0; }
//      selfDTotalHvvcoupl[0][0]=1;
//      selfDTotalHvvcoupl[3][1]=2.521;
//      mela.setMelaHiggsWidth(wPOLE);
//      pg1g4_pi2_VAMCFM = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDTotalHvvcoupl);
//*/
//
//      for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDTotalHvvcoupl[ii][jj] = 0; }
//      selfDTotalHvvcoupl[0][0]=1;
//      selfDTotalHvvcoupl[11][0]=-12046.01;
//      mela.setMelaHiggsWidth(wPOLE);
//      pg1g1prime2_VAMCFM = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDTotalHvvcoupl);
//
//
//      for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDTotalHvvcoupl[ii][jj] = 0; }
//      selfDTotalHvvcoupl[0][0]=1;
//      selfDTotalHvvcoupl[11][1]=-12046.01;
//      mela.setMelaHiggsWidth(wPOLE);
//      pg1g1prime2_pi2_VAMCFM = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDTotalHvvcoupl);
//
//
//      if (recorded==0){
//
//      }
//
//
//      pg1g2_VAJHU -= (p0plus_VAJHU + p0hplus_VAJHU);
//      pg1g2_pi2_VAJHU -= (p0plus_VAJHU + p0hplus_VAJHU);
//      pg1g4_VAJHU -= (p0plus_VAJHU + p0minus_VAJHU);
//      pg1g4_pi2_VAJHU -= (p0plus_VAJHU + p0minus_VAJHU);
//      pg1g1prime2_VAJHU -= (p0plus_VAJHU + p0_g1prime2_VAJHU);
//      pg1g1prime2_pi2_VAJHU -= (p0plus_VAJHU + p0_g1prime2_VAJHU);
//      pg1g2_VAMCFM -= (p0plus_VAMCFM + p0hplus_VAMCFM);
//      pg1g2_pi2_VAMCFM -= (p0plus_VAMCFM + p0hplus_VAMCFM);
//      pg1g4_VAMCFM -= (p0plus_VAMCFM + p0minus_VAMCFM);
//      pg1g4_pi2_VAMCFM -= (p0plus_VAMCFM + p0minus_VAMCFM);
//      pg1g1prime2_VAMCFM -= (p0plus_VAMCFM + p0_g1prime2_VAMCFM);
//      pg1g1prime2_pi2_VAMCFM -= (p0plus_VAMCFM + p0_g1prime2_VAMCFM);
//      pvbfg1g2_VAJHU -= (pvbf0plus_VAJHU + pvbf0hplus_VAJHU);
//      pvbfg1g2_pi2_VAJHU -= (pvbf0plus_VAJHU + pvbf0hplus_VAJHU);
//      pvbfg1g4_VAJHU -= (pvbf0plus_VAJHU + pvbf0minus_VAJHU);
//      pvbfg1g4_pi2_VAJHU -= (pvbf0plus_VAJHU + pvbf0minus_VAJHU);
//      pvbfg1g1prime2_VAJHU -= (pvbf0plus_VAJHU + pvbf0_g1prime2_VAJHU);
//      pvbfg1g1prime2_pi2_VAJHU -= (pvbf0plus_VAJHU + pvbf0_g1prime2_VAJHU);
//      pvbfg1g2_VAMCFM -= (pvbf0plus_VAMCFM + pvbf0hplus_VAMCFM);
//      pvbfg1g2_pi2_VAMCFM -= (pvbf0plus_VAMCFM + pvbf0hplus_VAMCFM);
//      pvbfg1g4_VAMCFM -= (pvbf0plus_VAMCFM + pvbf0minus_VAMCFM);
//      pvbfg1g4_pi2_VAMCFM -= (pvbf0plus_VAMCFM + pvbf0minus_VAMCFM);
//      pvbfg1g1prime2_VAMCFM -= (pvbf0plus_VAMCFM + pvbf0_g1prime2_VAMCFM);
//      pvbfg1g1prime2_pi2_VAMCFM -= (pvbf0plus_VAMCFM + pvbf0_g1prime2_VAMCFM);
//
//      p0hplus_VAJHU/=p0plus_VAJHU;
//      p0minus_VAJHU/=p0plus_VAJHU;
//      p0_g1prime2_VAJHU/=p0plus_VAJHU;
//      pg1g1prime2_VAJHU/=p0plus_VAJHU;
//      pg1g1prime2_pi2_VAJHU/=p0plus_VAJHU;
//      pg1g2_VAJHU/=p0plus_VAJHU;
//      pg1g2_pi2_VAJHU/=p0plus_VAJHU;
//      pg1g4_VAJHU/=p0plus_VAJHU;
//      pg1g4_pi2_VAJHU/=p0plus_VAJHU;
//      p0hplus_VAMCFM/=p0plus_VAMCFM;
//      p0minus_VAMCFM/=p0plus_VAMCFM;
//      p0_g1prime2_VAMCFM/=p0plus_VAMCFM;
//      pg1g1prime2_VAMCFM/=p0plus_VAMCFM;
//      pg1g1prime2_pi2_VAMCFM/=p0plus_VAMCFM;
//      pg1g2_VAMCFM/=p0plus_VAMCFM;
//      pg1g2_pi2_VAMCFM/=p0plus_VAMCFM;
//      pg1g4_VAMCFM/=p0plus_VAMCFM;
//      pg1g4_pi2_VAMCFM/=p0plus_VAMCFM;
//      pvbf0hplus_VAJHU/=pvbf0plus_VAJHU;
//      pvbf0minus_VAJHU/=pvbf0plus_VAJHU;
//      pvbf0_g1prime2_VAJHU/=pvbf0plus_VAJHU;
//      pvbfg1g1prime2_VAJHU/=pvbf0plus_VAJHU;
//      pvbfg1g1prime2_pi2_VAJHU/=pvbf0plus_VAJHU;
//      pvbfg1g2_VAJHU/=pvbf0plus_VAJHU;
//      pvbfg1g2_pi2_VAJHU/=pvbf0plus_VAJHU;
//      pvbfg1g4_VAJHU/=pvbf0plus_VAJHU;
//      pvbfg1g4_pi2_VAJHU/=pvbf0plus_VAJHU;
//      pvbf0hplus_VAMCFM/=pvbf0plus_VAMCFM;
//      pvbf0minus_VAMCFM/=pvbf0plus_VAMCFM;
//      pvbf0_g1prime2_VAMCFM/=pvbf0plus_VAMCFM;
//      pvbfg1g1prime2_VAMCFM/=pvbf0plus_VAMCFM;
//      pvbfg1g1prime2_pi2_VAMCFM/=pvbf0plus_VAMCFM;
//      pvbfg1g2_VAMCFM/=pvbf0plus_VAMCFM;
//      pvbfg1g2_pi2_VAMCFM/=pvbf0plus_VAMCFM;
//      pvbfg1g4_VAMCFM/=pvbf0plus_VAMCFM;
//      pvbfg1g4_pi2_VAMCFM/=pvbf0plus_VAMCFM;
//
//      for (int xx=0; xx<11; xx++){
//        for (int yy=0; yy<11; yy++){
//          pFullvbf0plus_VAJHU[xx][yy] = (rcd_pvbf0plus_VAJHU.MEsq)[xx][yy];
//          //pFullvbf0hplus_VAJHU[xx][yy] = (rcd_pvbf0hplus_VAJHU.MEsq)[xx][yy];
//          //pFullvbf0minus_VAJHU[xx][yy] = (rcd_pvbf0minus_VAJHU.MEsq)[xx][yy];
//          pFullvbf0_g1prime2_VAJHU[xx][yy] = (rcd_pvbf0_g1prime2_VAJHU.MEsq)[xx][yy];
//          /*pFullvbfg1g2_VAJHU[xx][yy] = (rcd_pvbfg1g2_VAJHU.MEsq)[xx][yy];
//          pFullvbfg1g2_pi2_VAJHU[xx][yy] = (rcd_pvbfg1g2_pi2_VAJHU.MEsq)[xx][yy];
//          pFullvbfg1g4_VAJHU[xx][yy] = (rcd_pvbfg1g4_VAJHU.MEsq)[xx][yy];
//          pFullvbfg1g4_pi2_VAJHU[xx][yy] = (rcd_pvbfg1g4_pi2_VAJHU.MEsq)[xx][yy];
//          */pFullvbfg1g1prime2_VAJHU[xx][yy] = (rcd_pvbfg1g1prime2_VAJHU.MEsq)[xx][yy];
//          pFullvbfg1g1prime2_pi2_VAJHU[xx][yy] = (rcd_pvbfg1g1prime2_pi2_VAJHU.MEsq)[xx][yy];
//          //pFullvbf0plus_VAMCFM[xx][yy] = (rcd_pvbf0plus_VAMCFM.MEsq)[xx][yy];
//          //pFullvbf0hplus_VAMCFM[xx][yy] = (rcd_pvbf0hplus_VAMCFM.MEsq)[xx][yy];
//          //pFullvbf0minus_VAMCFM[xx][yy] = (rcd_pvbf0minus_VAMCFM.MEsq)[xx][yy];
//          pFullvbf0_g1prime2_VAMCFM[xx][yy] = (rcd_pvbf0_g1prime2_VAMCFM.MEsq)[xx][yy];
//          /*pFullvbfg1g2_VAMCFM[xx][yy] = (rcd_pvbfg1g2_VAMCFM.MEsq)[xx][yy];
//          pFullvbfg1g2_pi2_VAMCFM[xx][yy] = (rcd_pvbfg1g2_pi2_VAMCFM.MEsq)[xx][yy];
//          pFullvbfg1g4_VAMCFM[xx][yy] = (rcd_pvbfg1g4_VAMCFM.MEsq)[xx][yy];
//          pFullvbfg1g4_pi2_VAMCFM[xx][yy] = (rcd_pvbfg1g4_pi2_VAMCFM.MEsq)[xx][yy];
//          */pFullvbfg1g1prime2_VAMCFM[xx][yy] = (rcd_pvbfg1g1prime2_VAMCFM.MEsq)[xx][yy];
//          pFullvbfg1g1prime2_pi2_VAMCFM[xx][yy] = (rcd_pvbfg1g1prime2_pi2_VAMCFM.MEsq)[xx][yy];
//
//          /*pFullvbfg1g2_VAJHU[xx][yy] -= (pFullvbf0plus_VAJHU[xx][yy] + pFullvbf0hplus_VAJHU[xx][yy]);
//          pFullvbfg1g2_pi2_VAJHU[xx][yy] -= (pFullvbf0plus_VAJHU[xx][yy] + pFullvbf0hplus_VAJHU[xx][yy]);
//          pFullvbfg1g4_VAJHU[xx][yy] -= (pFullvbf0plus_VAJHU[xx][yy] + pFullvbf0minus_VAJHU[xx][yy]);
//          pFullvbfg1g4_pi2_VAJHU[xx][yy] -= (pFullvbf0plus_VAJHU[xx][yy] + pFullvbf0minus_VAJHU[xx][yy]);
//          */pFullvbfg1g1prime2_VAJHU[xx][yy] -= (pFullvbf0plus_VAJHU[xx][yy] + pFullvbf0_g1prime2_VAJHU[xx][yy]);
//          pFullvbfg1g1prime2_pi2_VAJHU[xx][yy] -= (pFullvbf0plus_VAJHU[xx][yy] + pFullvbf0_g1prime2_VAJHU[xx][yy]);
//          /*pFullvbfg1g2_VAMCFM[xx][yy] -= (pFullvbf0plus_VAMCFM[xx][yy] + pFullvbf0hplus_VAMCFM[xx][yy]);
//          pFullvbfg1g2_pi2_VAMCFM[xx][yy] -= (pFullvbf0plus_VAMCFM[xx][yy] + pFullvbf0hplus_VAMCFM[xx][yy]);
//          pFullvbfg1g4_VAMCFM[xx][yy] -= (pFullvbf0plus_VAMCFM[xx][yy] + pFullvbf0minus_VAMCFM[xx][yy]);
//          pFullvbfg1g4_pi2_VAMCFM[xx][yy] -= (pFullvbf0plus_VAMCFM[xx][yy] + pFullvbf0minus_VAMCFM[xx][yy]);
//          */pFullvbfg1g1prime2_VAMCFM[xx][yy] -= (pFullvbf0plus_VAMCFM[xx][yy] + pFullvbf0_g1prime2_VAMCFM[xx][yy]);
//          pFullvbfg1g1prime2_pi2_VAMCFM[xx][yy] -= (pFullvbf0plus_VAMCFM[xx][yy] + pFullvbf0_g1prime2_VAMCFM[xx][yy]);
//
//          //pFullvbf0hplus_VAJHU[xx][yy]/=pFullvbf0plus_VAJHU[xx][yy];
//          //pFullvbf0minus_VAJHU[xx][yy]/=pFullvbf0plus_VAJHU[xx][yy];
//          pFullvbf0_g1prime2_VAJHU[xx][yy]/=pFullvbf0plus_VAJHU[xx][yy];
//          pFullvbfg1g1prime2_VAJHU[xx][yy]/=pFullvbf0plus_VAJHU[xx][yy];
//          pFullvbfg1g1prime2_pi2_VAJHU[xx][yy]/=pFullvbf0plus_VAJHU[xx][yy];
//          /*pFullvbfg1g2_VAJHU[xx][yy]/=pFullvbf0plus_VAJHU[xx][yy];
//          pFullvbfg1g2_pi2_VAJHU[xx][yy]/=pFullvbf0plus_VAJHU[xx][yy];
//          pFullvbfg1g4_VAJHU[xx][yy]/=pFullvbf0plus_VAJHU[xx][yy];
//          pFullvbfg1g4_pi2_VAJHU[xx][yy]/=pFullvbf0plus_VAJHU[xx][yy];
//          pFullvbf0hplus_VAMCFM[xx][yy]/=pFullvbf0plus_VAMCFM[xx][yy];
//          pFullvbf0minus_VAMCFM[xx][yy]/=pFullvbf0plus_VAMCFM[xx][yy];
//          */pFullvbf0_g1prime2_VAMCFM[xx][yy]/=pFullvbf0plus_VAMCFM[xx][yy];
//          pFullvbfg1g1prime2_VAMCFM[xx][yy]/=pFullvbf0plus_VAMCFM[xx][yy];
//          pFullvbfg1g1prime2_pi2_VAMCFM[xx][yy]/=pFullvbf0plus_VAMCFM[xx][yy];
//          /*pFullvbfg1g2_VAMCFM[xx][yy]/=pFullvbf0plus_VAMCFM[xx][yy];
//          pFullvbfg1g2_pi2_VAMCFM[xx][yy]/=pFullvbf0plus_VAMCFM[xx][yy];
//          pFullvbfg1g4_VAMCFM[xx][yy]/=pFullvbf0plus_VAMCFM[xx][yy];
//          pFullvbfg1g4_pi2_VAMCFM[xx][yy]/=pFullvbf0plus_VAMCFM[xx][yy];
//        */}
//      }
//
//      newtree->Fill();
//      recorded++;
//      //if (recorded==1) break;
//    }
//  }
//
//
//  foutput->WriteTObject(newtree);
//  delete newtree;
//  foutput->Close();
//  finput->Close();
//}
//


void testME_FullMELA_MCFMBSMHiggs_Ping(int flavor=1){
	int erg_tev=8;
	float mPOLE=125.6;
	float wPOLE=4.15e-3;
	char TREE_NAME[] = "GenTree";

//	TVar::VerbosityLevel verbosity = TVar::INFO;

	Mela mela(erg_tev,mPOLE);
	TGraph* vaScale_4e = mela.vaScale_4e;
	TGraph* vaScale_4mu = mela.vaScale_4mu;
	TGraph* vaScale_2e2mu = mela.vaScale_2e2mu;
	TGraph* DggZZ_scalefactor = mela.DggZZ_scalefactor;
	mela.resetMCFM_EWKParameters(1.16639E-05,1./128.,79.9549392,91.1876,0.23119);

	TFile* finput = new TFile("/afs/cern.ch/work/u/usarica/WidthAnalysis/LHC_8TeV/2mu2e/Sig/HZZ4l-125_6-8TeV-Sig.root","read");

	float p0plus_VAJHU;
	float p0hplus_VAJHU;
	float p0minus_VAJHU;
	float p0_g1prime2_VAJHU;
	float pg1g1prime2_VAJHU;
	float pg1g1prime2_pi2_VAJHU;
	float pg1g2_VAJHU;
	float pg1g2_pi2_VAJHU;
	float pg1g4_VAJHU;
	float pg1g4_pi2_VAJHU;

	float p0plus_VAJHU_NEW;
	float p0hplus_VAJHU_NEW;
	float p0minus_VAJHU_NEW;
	float p0_g1prime2_VAJHU_NEW;
	float pg1g1prime2_VAJHU_NEW;
	float pg1g1prime2_pi2_VAJHU_NEW;
	float pg1g2_VAJHU_NEW;
	float pg1g2_pi2_VAJHU_NEW;
	float pg1g4_VAJHU_NEW;
	float pg1g4_pi2_VAJHU_NEW;

	float ggHZZ_prob_pure,ggHZZ_prob_int,ggZZ_prob_Total;
	float bkg_VAMCFM,ggzz_VAMCFM;
	float ggHZZ_prob_pure_NEW,ggHZZ_prob_int_NEW,ggZZ_prob_Total_NEW;
	float bkg_VAMCFM_NEW,ggzz_VAMCFM_NEW;
	float bkg_VAMCFM_STU,bkg_VAMCFM_TU,bkg_VAMCFM_S;

	float p0plus_m4l, p0plus_m4l_ScaleUp, p0plus_m4l_ScaleDown, p0plus_m4l_ResUp, p0plus_m4l_ResDown;
	float bkg_m4l, bkg_m4l_ScaleUp, bkg_m4l_ScaleDown, bkg_m4l_ResUp, bkg_m4l_ResDown;

	float mzz = 126.; 
	float m1 = 91.471450;
	float m2 = 12.139782;
	float h1 = 0.2682896;
	float h2 = 0.1679779;
	float phi = 1.5969792;
	float hs = -0.727181;
	float phi1 = 1.8828257;
	float GenLep1Id,GenLep2Id,GenLep3Id,GenLep4Id;

	TTree* tree = (TTree*) finput->Get(TREE_NAME);
	tree->SetBranchAddress("GenZZMass", &mzz);
	tree->SetBranchAddress("GenZ1Mass", &m1);
	tree->SetBranchAddress("GenZ2Mass", &m2);
	tree->SetBranchAddress("GenhelcosthetaZ1", &h1);
	tree->SetBranchAddress("GenhelcosthetaZ2", &h2);
	tree->SetBranchAddress("Genhelphi", &phi);
	tree->SetBranchAddress("Gencosthetastar", &hs);
	tree->SetBranchAddress("GenphistarZ1", &phi1);
	tree->SetBranchAddress("GenLep1Id", &GenLep1Id);
	tree->SetBranchAddress("GenLep2Id", &GenLep2Id);
	tree->SetBranchAddress("GenLep3Id", &GenLep3Id);
	tree->SetBranchAddress("GenLep4Id", &GenLep4Id);

	int nEntries = tree->GetEntries();
	double selfDHvvcoupl[SIZE_HVV][2] = { { 0 } };
	double ggvvcoupl[2]={0,0};
	mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
	for (int ev = 0; ev < 1; ev++){
		tree->GetEntry(ev);

		if (flavor == 1){
			GenLep1Id=11;
			GenLep2Id=-11;
			GenLep3Id=11;
			GenLep4Id=-11;
		}
		int lepIdOrdered[4] = { GenLep1Id, GenLep2Id, GenLep3Id, GenLep4Id };
		float angularOrdered[8] = { mzz, m1, m2, hs, h1, h2, phi, phi1 };

		mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[0][0]=1;
		cout << "JHUGen g1=1" << endl;
		p0plus_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		cout << endl;
/*
		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[1][0]=1.638;
		p0hplus_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[3][0]=2.521;
		p0minus_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
*/
		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[11][0]=12046.01;
		cout << "JHUGen g1prime2=12046.01i" << endl;
		p0_g1prime2_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[1][0]=1.638;
		cout << "JHUGen g1=1; g2=1.638" << endl;
		pg1g2_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		cout << endl;

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[1][1]=1.638;
		cout << "JHUGen g1=1; g2=1.638i" << endl;
		pg1g2_pi2_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		cout << endl;

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[3][0]=2.521;
		cout << "JHUGen g1=1; g4=2.521" << endl;
		pg1g4_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		cout << endl;

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[3][1]=2.521;
		cout << "JHUGen g1=1; g4=2.521i" << endl;
		pg1g4_pi2_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		cout << endl;

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[11][0]=12046.01;
		cout << "JHUGen g1=1,g1prime2=12046.01" << endl;
		pg1g1prime2_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		cout << endl;

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[11][1]=12046.01;
		cout << "JHUGen g1=1,g1prime2=12046.01i" << endl;
		pg1g1prime2_pi2_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		cout << endl;


		mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[0][0]=1;
		mela.setMelaHiggsWidth(wPOLE);
		cout << "MCFM g1=1" << endl;
		mela.setMelaLeptonInterference(TVar::InterfOn);
		p0plus_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		cout << endl;
/*
		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[1][0]=1.638;
		mela.setMelaHiggsWidth(wPOLE);
		p0hplus_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[3][0]=2.521;
		mela.setMelaHiggsWidth(wPOLE);
		p0minus_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
*/
		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[11][0]=12046.01;
		mela.setMelaHiggsWidth(wPOLE);
		mela.setMelaLeptonInterference(TVar::InterfOn);
		cout << "MCFM g1prime2=12046.01" << endl;
		p0_g1prime2_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		cout << endl;

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[1][0]=1.638;
		mela.setMelaHiggsWidth(wPOLE);
		mela.setMelaLeptonInterference(TVar::InterfOn);
		cout << "MCFM g1=1,g2=1.638" << endl;
		pg1g2_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		cout << endl;

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[1][1]=1.638;
		mela.setMelaHiggsWidth(wPOLE);
		mela.setMelaLeptonInterference(TVar::InterfOn);
		cout << "MCFM g1=1,g2=1.638i" << endl;
		pg1g2_pi2_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		cout << endl;

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[3][0]=2.521;
		mela.setMelaHiggsWidth(wPOLE);
		mela.setMelaLeptonInterference(TVar::InterfOn);
		cout << "MCFM g1=1,g4=2.521" << endl;
		pg1g4_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		cout << endl;

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[3][1]=2.521;
		mela.setMelaHiggsWidth(wPOLE);
		mela.setMelaLeptonInterference(TVar::InterfOn);
		cout << "MCFM g1=1,g4=2.521i" << endl;
		pg1g4_pi2_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		cout << endl;

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[11][0]=12046.01;
		mela.setMelaHiggsWidth(wPOLE);
		mela.setMelaLeptonInterference(TVar::InterfOn);
		cout << "MCFM g1=1,g1prime2=12046.01" << endl;
		pg1g1prime2_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		cout << endl;

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[11][1]=12046.01;
		mela.setMelaHiggsWidth(wPOLE);
		mela.setMelaLeptonInterference(TVar::InterfOn);
		cout << "MCFM g1=1,g1prime2=12046.01i" << endl;
		pg1g1prime2_pi2_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		cout << endl;
	}

	finput->Close();
}


void testME_FullMELA_MCFMBSMHiggs_SignalReweight(int flavor=0){
	int erg_tev=8;
	float mPOLE=125.6;
	float wPOLE=4.15e-3;
	char TREE_NAME[] = "GenTree";

//	TVar::VerbosityLevel verbosity = TVar::INFO;

	Mela mela(erg_tev,mPOLE);
//	mela.resetMCFM_EWKParameters(1.16639E-05,1./128.,79.9549392,91.1876,0.23119);

	TFile* finput = new TFile(Form("/afs/cern.ch/work/u/usarica/WidthAnalysis/LHC_8TeV/%s/Sig/HZZ4l-125_6-8TeV-Sig.root",(flavor==0 ? "2mu2e" : "4e")),"read");
	TFile* foutput = new TFile(Form("HZZ4l-125_6-8TeV-Sig_%s_MCFMSignalReweighting.root",(flavor==0 ? "2mu2e" : "4e")),"recreate");

	float p0plus_VAJHU_NEW;
	float p0plus_VAJHU;
	float p0hplus_VAJHU_NEW;
	float p0minus_VAJHU_NEW;
	float p0_g1prime2_VAJHU_NEW;
	float p0_g1prime4_VAJHU_NEW;
	float p0_g1prime4_VAJHU_NEW2;
	float pg1g1prime2_VAJHU_NEW;
	float pg1g1prime2_pi2_VAJHU_NEW;
	float pg1g2_VAJHU_NEW;
	float pg1g2_pi2_VAJHU_NEW;
	float pg1g4_VAJHU_NEW;
	float pg1g4_pi2_VAJHU_NEW;

	float ggHZZ_prob_pure,ggHZZ_prob_int,ggZZ_prob_Total;
	float bkg_VAMCFM,ggzz_VAMCFM;
	float ggHZZ_prob_pure_NEW,ggHZZ_prob_int_NEW,ggZZ_prob_Total_NEW;
	float bkg_VAMCFM_NEW,ggzz_VAMCFM_NEW;
	float bkg_VAMCFM_STU,bkg_VAMCFM_TU,bkg_VAMCFM_S;

	float p0plus_m4l, p0plus_m4l_ScaleUp, p0plus_m4l_ScaleDown, p0plus_m4l_ResUp, p0plus_m4l_ResDown;
	float bkg_m4l, bkg_m4l_ScaleUp, bkg_m4l_ScaleDown, bkg_m4l_ResUp, bkg_m4l_ResDown;

	float mzz = 126.; 
	float m1 = 91.471450;
	float m2 = 12.139782;
	float h1 = 0.2682896;
	float h2 = 0.1679779;
	float phi = 1.5969792;
	float hs = -0.727181;
	float phi1 = 1.8828257;
	float GenLep1Id,GenLep2Id,GenLep3Id,GenLep4Id;

	TTree* tree = (TTree*) finput->Get(TREE_NAME);
	tree->SetBranchAddress("GenZZMass", &mzz);
	tree->SetBranchAddress("GenZ1Mass", &m1);
	tree->SetBranchAddress("GenZ2Mass", &m2);
	tree->SetBranchAddress("GenhelcosthetaZ1", &h1);
	tree->SetBranchAddress("GenhelcosthetaZ2", &h2);
	tree->SetBranchAddress("Genhelphi", &phi);
	tree->SetBranchAddress("Gencosthetastar", &hs);
	tree->SetBranchAddress("GenphistarZ1", &phi1);
	tree->SetBranchAddress("GenLep1Id", &GenLep1Id);
	tree->SetBranchAddress("GenLep2Id", &GenLep2Id);
	tree->SetBranchAddress("GenLep3Id", &GenLep3Id);
	tree->SetBranchAddress("GenLep4Id", &GenLep4Id);

	TTree* newtree = new TTree("TestTree","");
	newtree->Branch("GenZZMass", &mzz);
	newtree->Branch("GenZ1Mass", &m1);
	newtree->Branch("GenZ2Mass", &m2);
	newtree->Branch("GenhelcosthetaZ1", &h1);
	newtree->Branch("GenhelcosthetaZ2", &h2);
	newtree->Branch("Genhelphi", &phi);
	newtree->Branch("Gencosthetastar", &hs);
	newtree->Branch("GenphistarZ1", &phi1);
	newtree->Branch("GenLep1Id", &GenLep1Id);
	newtree->Branch("GenLep2Id", &GenLep2Id);
	newtree->Branch("GenLep3Id", &GenLep3Id);
	newtree->Branch("GenLep4Id", &GenLep4Id);

	newtree->Branch("p0plus_VAJHU_NEW",&p0plus_VAJHU_NEW);
	newtree->Branch("p0plus_VAJHU",&p0plus_VAJHU);
	newtree->Branch("p0hplus_VAJHU_NEW",&p0hplus_VAJHU_NEW);
	newtree->Branch("p0minus_VAJHU_NEW",&p0minus_VAJHU_NEW);
	newtree->Branch("p0_g1prime2_VAJHU_NEW",&p0_g1prime2_VAJHU_NEW);
	newtree->Branch("p0_g1prime4_VAJHU_NEW",&p0_g1prime4_VAJHU_NEW);
	newtree->Branch("p0_g1prime4_VAJHU_NEW2",&p0_g1prime4_VAJHU_NEW2);

	int nEntries = tree->GetEntries();
	double selfDHvvcoupl[SIZE_HVV][2] = { { 0 } };
	double ggvvcoupl[2]={0,0};
//	for (int ev = 0; ev < nEntries; ev++){
	for (int ev = 0; ev < 1000; ev++){
		tree->GetEntry(ev);

		int lepIdOrdered[4] = { GenLep1Id, GenLep2Id, GenLep3Id, GenLep4Id };
		float angularOrdered[8] = { mzz, m1, m2, hs, h1, h2, phi, phi1 };

		mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[0][0]=1;
		mela.setMelaLeptonInterference(TVar::InterfOn);
		mela.setMelaHiggsWidth(wPOLE);
		p0plus_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[1][0]=1.638;
		mela.setMelaLeptonInterference(TVar::InterfOn);
		mela.setMelaHiggsWidth(wPOLE);
		p0hplus_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[3][0]=2.521;
		mela.setMelaLeptonInterference(TVar::InterfOn);
		mela.setMelaHiggsWidth(wPOLE);
		p0minus_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[11][0]=-12046.01;
		mela.setMelaLeptonInterference(TVar::InterfOn);
		mela.setMelaHiggsWidth(wPOLE);
		p0_g1prime2_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[13][0]=-pow(10000.0/mPOLE,2);
		mela.setMelaLeptonInterference(TVar::InterfOn);
		mela.setMelaHiggsWidth(wPOLE);
		p0_g1prime4_VAJHU_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);


		mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[13][0]=-pow(10000.0/mPOLE,2);
		mela.setMelaLeptonInterference(TVar::InterfOn);
		mela.setMelaHiggsWidth(wPOLE);
		p0_g1prime4_VAJHU_NEW2 = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[0][0]=1;
		mela.setMelaLeptonInterference(TVar::InterfOn);
		mela.setMelaHiggsWidth(wPOLE);
		p0plus_VAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		p0hplus_VAJHU_NEW/=p0plus_VAJHU_NEW;
		p0minus_VAJHU_NEW/=p0plus_VAJHU_NEW;
		p0_g1prime2_VAJHU_NEW/=p0plus_VAJHU_NEW;
		p0_g1prime4_VAJHU_NEW/=p0plus_VAJHU_NEW;
		p0_g1prime4_VAJHU_NEW2/=p0plus_VAJHU;


		newtree->Fill();
	};


	foutput->WriteTObject(newtree);
	foutput->Close();
	finput->Close();
};


void testME_FullMELA_FullSimMCValidation_NonSpin0(int flavor=1){
	int erg_tev=8;
	float mPOLE=125.6;
	float wPOLE=4.15e-3;
	char TREE_NAME[] = "SelectedTree";

//	TVar::VerbosityLevel verbosity = TVar::INFO;

	Mela mela(erg_tev,mPOLE);
	TGraph* vaScale_4e = mela.vaScale_4e;
	TGraph* vaScale_4mu = mela.vaScale_4mu;
	TGraph* vaScale_2e2mu = mela.vaScale_2e2mu;
	TGraph* DggZZ_scalefactor = mela.DggZZ_scalefactor;

  TString cinput_main = "/scratch0/hep/ianderso/CJLST/140519/PRODFSR_8TeV";
  TFile* finput = new TFile(Form("%s/%s/HZZ4lTree_powheg15jhuGenV3-0PMH125.6.root", cinput_main.Data(), (flavor==0 ? "2mu2e" : "4e")), "read");
	//TFile* finput = new TFile(Form("/afs/cern.ch/work/u/usarica/HZZ4l-125p6-FullAnalysis/LHC_8TeV/%s/HZZ4lTree_powheg15jhuGenV3-0PMH125.6.root",(flavor==0 ? "2mu2e" : "4e")),"read");
	TFile* foutput = new TFile(Form("HZZ4lTree_powheg15jhuGenV3-0PMH125.6_%s_OriginalMEValidationTestOnly_NonSpin0.root",(flavor==0 ? "2mu2e" : "4e")),"recreate");

	float p0plus_VAJHU;
	float p0hplus_VAJHU;
	float p0minus_VAJHU;
	float p0_g1prime2_VAJHU;
	float pg1g1prime2_VAJHU;
	float pg1g1prime2_pi2_VAJHU;
	float pg1g2_VAJHU;
	float pg1g2_pi2_VAJHU;
	float pg1g4_VAJHU;
	float pg1g4_pi2_VAJHU;

	float p0plus_VAJHU_NEW;
	float p0hplus_VAJHU_NEW;
	float p0minus_VAJHU_NEW;
	float p0_g1prime2_VAJHU_NEW;
	float pg1g1prime2_VAJHU_NEW;
	float pg1g1prime2_pi2_VAJHU_NEW;
	float pg1g2_VAJHU_NEW;
	float pg1g2_pi2_VAJHU_NEW;
	float pg1g4_VAJHU_NEW;
	float pg1g4_pi2_VAJHU_NEW;


	float ggHZZ_prob_pure,ggHZZ_prob_int,ggZZ_prob_Total;
	float bkg_VAMCFM,ggzz_VAMCFM;
	float ggHZZ_prob_pure_NEW,ggHZZ_prob_int_NEW,ggZZ_prob_Total_NEW,p0plus_VAMCFM_NEW_BSMOn;
	float bkg_VAMCFM_NEW,ggzz_VAMCFM_NEW;
	float bkg_VAMCFM_STU,bkg_VAMCFM_TU,bkg_VAMCFM_S;

	float p0plus_m4l, p0plus_m4l_ScaleUp, p0plus_m4l_ScaleDown, p0plus_m4l_ResUp, p0plus_m4l_ResDown;
	float bkg_m4l, bkg_m4l_ScaleUp, bkg_m4l_ScaleDown, bkg_m4l_ResUp, bkg_m4l_ResDown;

	// probabilities for production independent calculations
	float p1plusProdIndepVA, p1minusProdIndepVA;
	float p1plus_VAJHU,p1_VAJHU;
	float p2minimalVA, p2minimalVA_qqb, p2h2plusVA, p2h2plusVA_qqb;
	float p2h3plusVA, p2h3plusVA_qqb, p2hplusVA, p2hplusVA_qqb;
	float p2bplusVA, p2bplusVA_qqb, p2h6plusVA, p2h6plusVA_qqb, p2h7plusVA, p2h7plusVA_qqb, p2hminusVA, p2hminusVA_qqb, p2h9minusVA, p2h9minusVA_qqb, p2h10minusVA, p2h10minusVA_qqb;
	float p2mProdIndepVA, p2h2plusProdIndepVA, p2h3plusProdIndepVA, p2hplusProdIndepVA, p2bplusProdIndepVA, p2h6plusProdIndepVA, p2h7plusProdIndepVA, p2hminusProdIndepVA, p2h9minusProdIndepVA, p2h10minusProdIndepVA;

	float p1plusProdIndepVA_NEW, p1minusProdIndepVA_NEW;
	float p1plus_VAJHU_NEW,p1_VAJHU_NEW;
	float p2minimalVA_NEW, p2minimalVA_qqb_NEW, p2h2plusVA_NEW, p2h2plusVA_qqb_NEW;
	float p2h3plusVA_NEW, p2h3plusVA_qqb_NEW, p2hplusVA_NEW, p2hplusVA_qqb_NEW;
	float p2bplusVA_NEW, p2bplusVA_qqb_NEW, p2h6plusVA_NEW, p2h6plusVA_qqb_NEW, p2h7plusVA_NEW, p2h7plusVA_qqb_NEW, p2hminusVA_NEW, p2hminusVA_qqb_NEW, p2h9minusVA_NEW, p2h9minusVA_qqb_NEW, p2h10minusVA_NEW, p2h10minusVA_qqb_NEW;
	float p2mProdIndepVA_NEW, p2h2plusProdIndepVA_NEW, p2h3plusProdIndepVA_NEW, p2hplusProdIndepVA_NEW, p2bplusProdIndepVA_NEW, p2h6plusProdIndepVA_NEW, p2h7plusProdIndepVA_NEW, p2hminusProdIndepVA_NEW, p2h9minusProdIndepVA_NEW, p2h10minusProdIndepVA_NEW;

	float mzz = 126.; 
	float m1 = 91.471450;
	float m2 = 12.139782;
	float h1 = 0.2682896;
	float h2 = 0.1679779;
	float phi = 1.5969792;
	float hs = -0.727181;
	float phi1 = 1.8828257;

	TTree* tree = (TTree*) finput->Get(TREE_NAME);
	tree->SetBranchAddress("ZZMass", &mzz);
	tree->SetBranchAddress("Z1Mass", &m1);
	tree->SetBranchAddress("Z2Mass", &m2);
	tree->SetBranchAddress("helcosthetaZ1", &h1);
	tree->SetBranchAddress("helcosthetaZ2", &h2);
	tree->SetBranchAddress("helphi", &phi);
	tree->SetBranchAddress("costhetastar", &hs);
	tree->SetBranchAddress("phistarZ1", &phi1);

	tree->SetBranchAddress("p1plus_prodIndep_VAJHU", &p1plusProdIndepVA);
	tree->SetBranchAddress("p1_prodIndep_VAJHU", &p1minusProdIndepVA);
	tree->SetBranchAddress("p1plus_VAJHU", &p1plus_VAJHU);
	tree->SetBranchAddress("p1_VAJHU", &p1_VAJHU);
	tree->SetBranchAddress("p2_VAJHU", &p2minimalVA);
	tree->SetBranchAddress("p2qqb_VAJHU", &p2minimalVA_qqb);
	tree->SetBranchAddress("p2h2plus_gg_VAJHU", &p2h2plusVA);
	tree->SetBranchAddress("p2h2plus_qqbar_VAJHU", &p2h2plusVA_qqb);
	tree->SetBranchAddress("p2h3plus_gg_VAJHU", &p2h3plusVA);
	tree->SetBranchAddress("p2h3plus_qqbar_VAJHU", &p2h3plusVA_qqb);
	tree->SetBranchAddress("p2hplus_VAJHU", &p2hplusVA);
	tree->SetBranchAddress("p2hplus_qqb_VAJHU", &p2hplusVA_qqb);
	tree->SetBranchAddress("p2bplus_VAJHU", &p2bplusVA);
	tree->SetBranchAddress("p2bplus_qqb_VAJHU", &p2bplusVA_qqb);
	tree->SetBranchAddress("p2h6plus_gg_VAJHU", &p2h6plusVA);
	tree->SetBranchAddress("p2h6plus_qqbar_VAJHU", &p2h6plusVA_qqb);
	tree->SetBranchAddress("p2h7plus_gg_VAJHU", &p2h7plusVA);
	tree->SetBranchAddress("p2h7plus_qqbar_VAJHU", &p2h7plusVA_qqb);
	tree->SetBranchAddress("p2hminus_VAJHU", &p2hminusVA);
	tree->SetBranchAddress("p2hminus_qqb_VAJHU", &p2hminusVA_qqb);
	tree->SetBranchAddress("p2h9minus_gg_VAJHU", &p2h9minusVA);
	tree->SetBranchAddress("p2h9minus_qqbar_VAJHU", &p2h9minusVA_qqb);
	tree->SetBranchAddress("p2h10minus_gg_VAJHU", &p2h10minusVA);
	tree->SetBranchAddress("p2h10minus_qqbar_VAJHU", &p2h10minusVA_qqb);
	tree->SetBranchAddress("p2_prodIndep_VAJHU", &p2mProdIndepVA);
	tree->SetBranchAddress("p2h2plus_prodIndep_VAJHU", &p2h2plusProdIndepVA);
	tree->SetBranchAddress("p2h3plus_prodIndep_VAJHU", &p2h3plusProdIndepVA);
	tree->SetBranchAddress("p2hplus_prodIndep_VAJHU", &p2hplusProdIndepVA);
	tree->SetBranchAddress("p2bplus_prodIndep_VAJHU", &p2bplusProdIndepVA);
	tree->SetBranchAddress("p2h6plus_prodIndep_VAJHU", &p2h6plusProdIndepVA);
	tree->SetBranchAddress("p2h7plus_prodIndep_VAJHU", &p2h7plusProdIndepVA);
	tree->SetBranchAddress("p2hminus_prodIndep_VAJHU", &p2hminusProdIndepVA);
	tree->SetBranchAddress("p2h9minus_prodIndep_VAJHU", &p2h9minusProdIndepVA);
	tree->SetBranchAddress("p2h10minus_prodIndep_VAJHU", &p2h10minusProdIndepVA);


	TTree* newtree = new TTree("TestTree","");
	newtree->Branch("ZZMass", &mzz);
	newtree->Branch("Z1Mass", &m1);
	newtree->Branch("Z2Mass", &m2);
	newtree->Branch("helcosthetaZ1", &h1);
	newtree->Branch("helcosthetaZ2", &h2);
	newtree->Branch("helphi", &phi);
	newtree->Branch("costhetastar", &hs);
	newtree->Branch("phistarZ1", &phi1);

	newtree->Branch("p1plus_prodIndep_VAJHU", &p1plusProdIndepVA);
	newtree->Branch("p1_prodIndep_VAJHU", &p1minusProdIndepVA);
	newtree->Branch("p1plus_VAJHU", &p1plus_VAJHU);
	newtree->Branch("p1_VAJHU", &p1_VAJHU);
	newtree->Branch("p2_VAJHU", &p2minimalVA);
	newtree->Branch("p2qqb_VAJHU", &p2minimalVA_qqb);
	newtree->Branch("p2h2plus_gg_VAJHU", &p2h2plusVA);
	newtree->Branch("p2h2plus_qqbar_VAJHU", &p2h2plusVA_qqb);
	newtree->Branch("p2h3plus_gg_VAJHU", &p2h3plusVA);
	newtree->Branch("p2h3plus_qqbar_VAJHU", &p2h3plusVA_qqb);
	newtree->Branch("p2hplus_VAJHU", &p2hplusVA);
	newtree->Branch("p2hplus_qqb_VAJHU", &p2hplusVA_qqb);
	newtree->Branch("p2bplus_VAJHU", &p2bplusVA);
	newtree->Branch("p2bplus_qqb_VAJHU", &p2bplusVA_qqb);
	newtree->Branch("p2h6plus_gg_VAJHU", &p2h6plusVA);
	newtree->Branch("p2h6plus_qqbar_VAJHU", &p2h6plusVA_qqb);
	newtree->Branch("p2h7plus_gg_VAJHU", &p2h7plusVA);
	newtree->Branch("p2h7plus_qqbar_VAJHU", &p2h7plusVA_qqb);
	newtree->Branch("p2hminus_VAJHU", &p2hminusVA);
	newtree->Branch("p2hminus_qqb_VAJHU", &p2hminusVA_qqb);
	newtree->Branch("p2h9minus_gg_VAJHU", &p2h9minusVA);
	newtree->Branch("p2h9minus_qqbar_VAJHU", &p2h9minusVA_qqb);
	newtree->Branch("p2h10minus_gg_VAJHU", &p2h10minusVA);
	newtree->Branch("p2h10minus_qqbar_VAJHU", &p2h10minusVA_qqb);
	newtree->Branch("p2_prodIndep_VAJHU", &p2mProdIndepVA);
	newtree->Branch("p2h2plus_prodIndep_VAJHU", &p2h2plusProdIndepVA);
	newtree->Branch("p2h3plus_prodIndep_VAJHU", &p2h3plusProdIndepVA);
	newtree->Branch("p2hplus_prodIndep_VAJHU", &p2hplusProdIndepVA);
	newtree->Branch("p2bplus_prodIndep_VAJHU", &p2bplusProdIndepVA);
	newtree->Branch("p2h6plus_prodIndep_VAJHU", &p2h6plusProdIndepVA);
	newtree->Branch("p2h7plus_prodIndep_VAJHU", &p2h7plusProdIndepVA);
	newtree->Branch("p2hminus_prodIndep_VAJHU", &p2hminusProdIndepVA);
	newtree->Branch("p2h9minus_prodIndep_VAJHU", &p2h9minusProdIndepVA);
	newtree->Branch("p2h10minus_prodIndep_VAJHU", &p2h10minusProdIndepVA);

	newtree->Branch("p1plus_prodIndep_VAJHU_NEW", &p1plusProdIndepVA_NEW);
	newtree->Branch("p1_prodIndep_VAJHU_NEW", &p1minusProdIndepVA_NEW);
	newtree->Branch("p1plus_VAJHU_NEW", &p1plus_VAJHU_NEW);
	newtree->Branch("p1_VAJHU_NEW", &p1_VAJHU_NEW);

	newtree->Branch("p2_VAJHU_NEW", &p2minimalVA_NEW);
	newtree->Branch("p2qqb_VAJHU_NEW", &p2minimalVA_qqb_NEW);
	newtree->Branch("p2h2plus_gg_VAJHU_NEW", &p2h2plusVA_NEW);
	newtree->Branch("p2h2plus_qqbar_VAJHU_NEW", &p2h2plusVA_qqb_NEW);
	newtree->Branch("p2h3plus_gg_VAJHU_NEW", &p2h3plusVA_NEW);
	newtree->Branch("p2h3plus_qqbar_VAJHU_NEW", &p2h3plusVA_qqb_NEW);
	newtree->Branch("p2hplus_VAJHU_NEW", &p2hplusVA_NEW);
	newtree->Branch("p2hplus_qqb_VAJHU_NEW", &p2hplusVA_qqb_NEW);
	newtree->Branch("p2bplus_VAJHU_NEW", &p2bplusVA_NEW);
	newtree->Branch("p2bplus_qqb_VAJHU_NEW", &p2bplusVA_qqb_NEW);
	newtree->Branch("p2h6plus_gg_VAJHU_NEW", &p2h6plusVA_NEW);
	newtree->Branch("p2h6plus_qqbar_VAJHU_NEW", &p2h6plusVA_qqb_NEW);
	newtree->Branch("p2h7plus_gg_VAJHU_NEW", &p2h7plusVA_NEW);
	newtree->Branch("p2h7plus_qqbar_VAJHU_NEW", &p2h7plusVA_qqb_NEW);
	newtree->Branch("p2hminus_VAJHU_NEW", &p2hminusVA_NEW);
	newtree->Branch("p2hminus_qqb_VAJHU_NEW", &p2hminusVA_qqb_NEW);
	newtree->Branch("p2h9minus_gg_VAJHU_NEW", &p2h9minusVA_NEW);
	newtree->Branch("p2h9minus_qqbar_VAJHU_NEW", &p2h9minusVA_qqb_NEW);
	newtree->Branch("p2h10minus_gg_VAJHU_NEW", &p2h10minusVA_NEW);
	newtree->Branch("p2h10minus_qqbar_VAJHU_NEW", &p2h10minusVA_qqb_NEW);
	newtree->Branch("p2_prodIndep_VAJHU_NEW", &p2mProdIndepVA_NEW);
	newtree->Branch("p2h2plus_prodIndep_VAJHU_NEW", &p2h2plusProdIndepVA_NEW);
	newtree->Branch("p2h3plus_prodIndep_VAJHU_NEW", &p2h3plusProdIndepVA_NEW);
	newtree->Branch("p2hplus_prodIndep_VAJHU_NEW", &p2hplusProdIndepVA_NEW);
	newtree->Branch("p2bplus_prodIndep_VAJHU_NEW", &p2bplusProdIndepVA_NEW);
	newtree->Branch("p2h6plus_prodIndep_VAJHU_NEW", &p2h6plusProdIndepVA_NEW);
	newtree->Branch("p2h7plus_prodIndep_VAJHU_NEW", &p2h7plusProdIndepVA_NEW);
	newtree->Branch("p2hminus_prodIndep_VAJHU_NEW", &p2hminusProdIndepVA_NEW);
	newtree->Branch("p2h9minus_prodIndep_VAJHU_NEW", &p2h9minusProdIndepVA_NEW);
	newtree->Branch("p2h10minus_prodIndep_VAJHU_NEW", &p2h10minusProdIndepVA_NEW);

	int nEntries = tree->GetEntries();
	double selfDHvvcoupl[SIZE_HVV][2] = { { 0 } };
	double ggvvcoupl[2]={0,0};
	mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
//	for (int ev = 0; ev < nEntries; ev++){
//	for (int ev = 0; ev < 100; ev++){
	for (int ev = 0; ev < 100; ev++){
		tree->GetEntry(ev);

		int lepIdOrdered[4]={ 11,-11,11,-11 };
		if (flavor == 0){
			lepIdOrdered[0]=13;lepIdOrdered[1]=-13;
		};
		float angularOrdered[8] = { mzz, m1, m2, hs, h1, h2, phi, phi1 };

		int myflavor=3;if(flavor==1) myflavor=1;

    // 1-
    mela.setProcess(TVar::H1minus, TVar::JHUGen, TVar::ZZQQB);
    mela.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,myflavor, p1_VAJHU_NEW);
    
    // production independent
    mela.setProcess(TVar::H1minus, TVar::JHUGen, TVar::ZZINDEPENDENT);
    mela.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,myflavor, p1minusProdIndepVA_NEW);
    
    // 1+
    mela.setProcess(TVar::H1plus, TVar::JHUGen, TVar::ZZQQB);
    mela.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,myflavor, p1plus_VAJHU_NEW);
    
    // production independent
    mela.setProcess(TVar::H1plus, TVar::JHUGen, TVar::ZZINDEPENDENT);
    mela.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,myflavor, p1plusProdIndepVA_NEW);

	// gg->2m+
    mela.setProcess(TVar::H2_g1g5, TVar::JHUGen, TVar::ZZGG);
    mela.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,myflavor, p2minimalVA_NEW);
    
    // qqb->2m+
    mela.setProcess(TVar::H2_g1g5, TVar::JHUGen, TVar::ZZQQB);
    mela.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,myflavor, p2minimalVA_qqb_NEW);
    
    // qqb->2m+
    mela.setProcess(TVar::H2_g1g5, TVar::JHUGen, TVar::ZZINDEPENDENT);
    mela.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,myflavor, p2mProdIndepVA_NEW);
    
    // 2h-
    mela.setProcess(TVar::H2_g8, TVar::JHUGen, TVar::ZZGG);
    mela.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,myflavor, p2hminusVA_NEW);
    
    // 2h+
    mela.setProcess(TVar::H2_g4, TVar::JHUGen, TVar::ZZGG);
    mela.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,myflavor, p2hplusVA_NEW);
    
    // 2b+
    mela.setProcess(TVar::H2_g5, TVar::JHUGen, TVar::ZZGG);
    mela.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,myflavor, p2bplusVA_NEW);

    mela.setProcess(TVar::H2_g2, TVar::JHUGen, TVar::ZZGG);
    mela.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,myflavor, p2h2plusVA_NEW);

    mela.setProcess(TVar::H2_g3, TVar::JHUGen, TVar::ZZGG);
    mela.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,myflavor, p2h3plusVA_NEW);

    mela.setProcess(TVar::H2_g6, TVar::JHUGen, TVar::ZZGG);
    mela.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,myflavor, p2h6plusVA_NEW);

    mela.setProcess(TVar::H2_g7, TVar::JHUGen, TVar::ZZGG);
    mela.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,myflavor, p2h7plusVA_NEW);

    mela.setProcess(TVar::H2_g9, TVar::JHUGen, TVar::ZZGG);
    mela.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,myflavor, p2h9minusVA_NEW);

    mela.setProcess(TVar::H2_g10, TVar::JHUGen, TVar::ZZGG);
    mela.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,myflavor, p2h10minusVA_NEW);

    mela.setProcess(TVar::H2_g8, TVar::JHUGen, TVar::ZZQQB);
    mela.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,myflavor, p2hminusVA_qqb_NEW);
    
    // 2h+
    mela.setProcess(TVar::H2_g4, TVar::JHUGen, TVar::ZZQQB);
    mela.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,myflavor, p2hplusVA_qqb_NEW);
    
    // 2b+
    mela.setProcess(TVar::H2_g5, TVar::JHUGen, TVar::ZZQQB);
    mela.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,myflavor, p2bplusVA_qqb_NEW);

    mela.setProcess(TVar::H2_g2, TVar::JHUGen, TVar::ZZQQB);
    mela.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,myflavor, p2h2plusVA_qqb_NEW);

    mela.setProcess(TVar::H2_g3, TVar::JHUGen, TVar::ZZQQB);
    mela.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,myflavor, p2h3plusVA_qqb_NEW);

    mela.setProcess(TVar::H2_g6, TVar::JHUGen, TVar::ZZQQB);
    mela.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,myflavor, p2h6plusVA_qqb_NEW);

    mela.setProcess(TVar::H2_g7, TVar::JHUGen, TVar::ZZQQB);
    mela.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,myflavor, p2h7plusVA_qqb_NEW);

    mela.setProcess(TVar::H2_g9, TVar::JHUGen, TVar::ZZQQB);
    mela.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,myflavor, p2h9minusVA_qqb_NEW);

    mela.setProcess(TVar::H2_g10, TVar::JHUGen, TVar::ZZQQB);
    mela.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,myflavor, p2h10minusVA_qqb_NEW);


    mela.setProcess(TVar::H2_g8, TVar::JHUGen, TVar::ZZINDEPENDENT);
    mela.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,myflavor, p2hminusProdIndepVA_NEW);
    
    // 2h+
    mela.setProcess(TVar::H2_g4, TVar::JHUGen, TVar::ZZINDEPENDENT);
    mela.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,myflavor, p2hplusProdIndepVA_NEW);
    
    // 2b+
    mela.setProcess(TVar::H2_g5, TVar::JHUGen, TVar::ZZINDEPENDENT);
    mela.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,myflavor, p2bplusProdIndepVA_NEW);

    mela.setProcess(TVar::H2_g2, TVar::JHUGen, TVar::ZZINDEPENDENT);
    mela.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,myflavor, p2h2plusProdIndepVA_NEW);

    mela.setProcess(TVar::H2_g3, TVar::JHUGen, TVar::ZZINDEPENDENT);
    mela.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,myflavor, p2h3plusProdIndepVA_NEW);

    mela.setProcess(TVar::H2_g6, TVar::JHUGen, TVar::ZZINDEPENDENT);
    mela.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,myflavor, p2h6plusProdIndepVA_NEW);

    mela.setProcess(TVar::H2_g7, TVar::JHUGen, TVar::ZZINDEPENDENT);
    mela.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,myflavor, p2h7plusProdIndepVA_NEW);

    mela.setProcess(TVar::H2_g9, TVar::JHUGen, TVar::ZZINDEPENDENT);
    mela.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,myflavor, p2h9minusProdIndepVA_NEW);

    mela.setProcess(TVar::H2_g10, TVar::JHUGen, TVar::ZZINDEPENDENT);
    mela.computeP(mzz, m1, m2, 
		    hs,h1,h2,phi,phi1,myflavor, p2h10minusProdIndepVA_NEW);


		newtree->Fill();
	}


	foutput->WriteTObject(newtree);
	foutput->Close();
	finput->Close();
}


void testME_FullMELA_PingWithFourMomenta(int flavor=0){
	int erg_tev=8;
	float mPOLE=125.0;
	float wPOLE=4.165e-3;

//	TVar::VerbosityLevel verbosity = TVar::INFO;

	Mela mela(erg_tev,mPOLE);
	mela.resetMCFM_EWKParameters(1.16639E-05,1./128.,79.9549392,91.1876,0.23119);

	float pVAJHU;
	float pVAMCFM;
	float pVAJHU_SM;
	float pVAMCFM_SM;

	float mzz = 125.; 
	float m1 = 91.471450;
	float m2 = 12.139782;
	float h1 = 0.2682896;
	float h2 = 0.1679779;
	float phi = 1.5969792;
	float hs = -0.727181;
	float phi1 = 1.8828257;
	float GenLep1Id,GenLep2Id,GenLep3Id,GenLep4Id;
/*
	const int nEntries = 3;
	double l1_array[nEntries][4] = {
		{1365.4973807340846,        10.289826593755228,        25.205694382277809,       -1365.2259480507332},
		{238.65751023078761,        9.2808858562825005,        15.827726043466324,       -237.95116187061188},
		{101.52463181523598,        27.359569630718468,      -0.90299073100241323,       -97.764458892691749}
	};
	double l2_array[nEntries][4] = {
		{22.786181013986834,      -0.15136300982222117,      -0.90077551414353962,       -22.767866345236371},
		{101.67043553688544,        2.1169375132239789,       0.77953005873937187,       -101.64540506443268},
		{24.717634703436786,       -1.1722249478288802,       -5.9599387484197646,       -23.959684558009428}
	};
	double l3_array[nEntries][4] = {
		{1895.7562628816693,        25.837804322120054,       -28.821887970086259,       -1895.3610513294620},
		{317.81904277258536,        2.5882005498984775,        21.352807448987718,       -317.09037005377883},
		{180.10885677707822,       -6.7240759244122792,        35.742176497019194,       -176.39865053838915}
	};
	double l4_array[nEntries][4] = {
		{471.71918486784784,       -35.976267906053060,        4.5169691019519895,       -470.32360615864354},
		{95.655512770627581,       -13.986023919404957,       -37.960063551193414,       -86.679881365440792},
		{49.137252081251319,       -19.463268758477309,       -28.879247017597017,       -34.664676589120688}
	};
*/
	const int nEntries = 6;
	double l1_array[nEntries][4] = {
		{51.374202,	25.924766,	12.290178,	42.616376},
		{51.374202,	25.924766,	12.290178,	42.616376},
		{1365.4973807340846,        10.289826593755228,        25.205694382277809,       -1365.2259480507332}, // Massless
		{1365.4973807340846,        10.289826593755228,        25.205694382277809,       -1365.2259480507332}, // Massless
		{1365.4973848483,        10.289826593755228,        25.205694382277809,       -1365.2259480507332}, // Muon via E
		{1365.4973848483,        10.289826593755228,        25.205694382277809,       -1365.2259480507332} // Muon via E
	};
	double l2_array[nEntries][4] = {
		{271.875752,	70.427173,	-11.138146,	261.769598},
		{21.481452,	9.489680,	-9.336587,	16.858699},
		{22.786181013986834,      -0.15136300982222117,      -0.90077551414353962,       -22.767866345236371},
		{471.71918486784784,       -35.976267906053060,        4.5169691019519895,       -470.32360615864354},
		{22.7864275656,      -0.15136300982222117,      -0.90077551414353962,       -22.767866345236371},
		{471.7191967775,       -35.976267906053060,        4.5169691019519895,       -470.32360615864354}
	};
	double l3_array[nEntries][4] = {
		{75.823478,	-16.640412,	23.246999,	70.227220},
		{75.823478,	-16.640412,	23.246999,	70.227220},
		{1895.7562628816693,        25.837804322120054,       -28.821887970086259,       -1895.3610513294620},
		{1895.7562628816693,        25.837804322120054,       -28.821887970086259,       -1895.3610513294620},
		{1895.7562658451,        25.837804322120054,       -28.821887970086259,       -1895.3610513294620},
		{1895.7562658451,        25.837804322120054,       -28.821887970086259,       -1895.3610513294620}
	};
	double l4_array[nEntries][4] = {
		{21.481452,	9.489680,	-9.336587,	16.858699},
		{271.875752,	70.427173,	-11.138146,	261.769598},
		{471.71918486784784,       -35.976267906053060,        4.5169691019519895,       -470.32360615864354},
		{22.786181013986834,      -0.15136300982222117,      -0.90077551414353962,       -22.767866345236371},
		{471.7191967775,       -35.976267906053060,        4.5169691019519895,       -470.32360615864354},
		{22.7864275656,      -0.15136300982222117,      -0.90077551414353962,       -22.767866345236371}
	};

	double selfDHvvcoupl[SIZE_HVV][2] = { { 0 } };
	double ggvvcoupl[2]={0,0};
	mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
	for (int ev = 0; ev < nEntries; ev++){
		if (flavor == 2){
			GenLep1Id=13;
			GenLep2Id=-13;
			GenLep3Id=11;
			GenLep4Id=-11;
		}
		else if (flavor == 1){
			GenLep1Id=11;
			GenLep2Id=-11;
			GenLep3Id=11;
			GenLep4Id=-11;
		}
		else{
			GenLep1Id=13;
			GenLep2Id=-13;
			GenLep3Id=13;
			GenLep4Id=-13;
		}
		int lepIdOrdered[4] = { GenLep1Id, GenLep2Id, GenLep3Id, GenLep4Id };
		TLorentzVector p1(l1_array[ev][1],l1_array[ev][2],l1_array[ev][3],l1_array[ev][0]);
		TLorentzVector p2(l2_array[ev][1],l2_array[ev][2],l2_array[ev][3],l2_array[ev][0]);
		TLorentzVector p3(l3_array[ev][1],l3_array[ev][2],l3_array[ev][3],l3_array[ev][0]);
		TLorentzVector p4(l4_array[ev][1],l4_array[ev][2],l4_array[ev][3],l4_array[ev][0]);
		TLorentzVector pZ1 = p1+p2;
		TLorentzVector pZ2 = p3+p4;
		TLorentzVector pZZ = pZ1+pZ2;
		mzz=pZZ.M();
		m1=pZ1.M();
		m2=pZ2.M();
		mela::computeAngles(p1, GenLep1Id,
			 p2, GenLep2Id,
			 p3, GenLep3Id,
			 p4, GenLep4Id,
			 hs, 
			 h1, 
			 h2, 
			 phi, 
			 phi1);
		float angularOrdered[8] = { mzz, m1, m2, hs, h1, h2, phi, phi1 };

		mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[0][0]=1.0;
//		selfDHvvcoupl[0][1]=8.6;
//		selfDHvvcoupl[1][0]=-3.5;
		selfDHvvcoupl[1][1]=-1.88;
//		selfDHvvcoupl[2][0]=1.0;
//		selfDHvvcoupl[2][1]=-0.88;
//		selfDHvvcoupl[3][0]=1.0;
		selfDHvvcoupl[3][1]=2.6;
		mela.setMelaHiggsWidth(wPOLE);
		mela.setMelaLeptonInterference(TVar::InterfOn);
		pVAMCFM = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[0][0]=1.0;
//		selfDHvvcoupl[0][1]=8.6;
//		selfDHvvcoupl[1][0]=-3.5;
		selfDHvvcoupl[1][1]=-1.88;
//		selfDHvvcoupl[2][0]=1.0;
//		selfDHvvcoupl[2][1]=-0.88;
//		selfDHvvcoupl[3][0]=1.0;
		selfDHvvcoupl[3][1]=2.6;
		mela.setMelaHiggsWidth(wPOLE);
		pVAJHU = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
/*
		cout << "BSM:\n";
		cout << "MCFM:\t";
		cout << pVAMCFM/1.45e-8 << endl;
		cout << "JHUGen:\t";
		cout << pVAJHU/1.45e-8 << endl;
		cout << "JHUGen/MCFM:\t";
		cout << pVAJHU/pVAMCFM << endl;
*/
		mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[0][0]=1.0;
		mela.setMelaHiggsWidth(wPOLE);
		mela.setMelaLeptonInterference(TVar::InterfOn);
		pVAMCFM_SM = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);

		for (int ii = 0; ii < SIZE_HVV; ii++){ for (int jj = 0; jj < 2; jj++) selfDHvvcoupl[ii][jj] = 0; }
		selfDHvvcoupl[0][0]=1.0;
		mela.setMelaHiggsWidth(wPOLE);
		pVAJHU_SM = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		cout << "SM:\n";
		cout << "MCFM:\t";
//		cout << pVAMCFM_SM/1.45e-8 << endl;
		cout << pVAMCFM_SM << endl;
		cout << "JHUGen:\t";
//		cout << pVAJHU_SM/1.45e-8 << endl;
		cout << pVAJHU_SM << endl;
		cout << "JHUGen/MCFM:\t";
		cout << pVAJHU_SM/pVAMCFM_SM << endl;

		cout << "Overall JHUGen/MCFM:\t";
		cout << (pVAJHU/pVAMCFM)/(pVAJHU_SM/pVAMCFM_SM) << endl;
	}
}


void testME_FullMELA_VBF_PingWithFourMomenta(){
  int erg_tev=13;
  float mPOLE=125.0;
  float wPOLE=4.07e-3;

  Mela mela(erg_tev, mPOLE);

  float p0plus_VAJHU;
  float phjj_VAJHU;
  float pvbf_VAJHU;

  float mzz = 125.;
  float m1 = 91.471450;
  float m2 = 12.139782;
  float h1 = 0.2682896;
  float h2 = 0.1679779;
  float phi = 1.5969792;
  float hs = -0.727181;
  float phi1 = 1.8828257;
  float GenLep1Id, GenLep2Id, GenLep3Id, GenLep4Id;

  const int nEntries = 3;
  double event_List[nEntries][6][4]={
    {
      { 44.7736, -0.0781905, 2.39686, -0.00514586 },
      { 22.3954, -0.112428, 1.60894, 0.00596709 },
      { 12.0699, 0.0934713, 0.892571, -0.00316125 },
      { 7.35361, 0.093612, 0.893655, -0.00214606 },
      { 149.721, 0.655096, -0.920623, 19.6351 },
      { 38.854, 0.0861686, 0.913102, 5.54952 }
    },
    {
      { 46.1946, -0.851673, 0.959117, 0.0179481 },
      { 34.8619, 0.280647, -2.61613, -0.00673518 },
      { 13.6921, 0.880684, 2.92821, -0.00127498 },
      { 8.86439, -0.59297, 0.987145, -0.00344459 },
      { 72.5579, -2.15782, -0.701562, 12.7015 },
      { 35.8925, 2.39948, -2.40373, 6.69218 }
    },
    {
      { 31.4753, 2.16809, 1.78843, 0.0276003 },
      { 25.778, 0.267627, 3.08483, -0.00408116 },
      { 33.7643, 0.813051, -1.13545, 0.13957 },
      { 29.3433, 1.20273, 1.57892, 0.13957 },
      { 110.499, 2.08684, 1.80324, 17.2814 },
      { 86.923, -2.20464, -0.956803, 10.2842 }
    }
  };
  int lepIDs[nEntries][4]={
    { 11, -11, 11, -11 },
    { 11, -11, -11, 11 },
    { -11, 11, 13, -13 }
  };

  double selfDHvvcoupl[SIZE_HVV][2] ={ { 0 } };
  double ggvvcoupl[2]={ 0, 0 };
  for (int ev = 0; ev < nEntries; ev++){
    TLorentzVector nullFourVector(0, 0, 0, 0);


    GenLep1Id=lepIDs[ev][0];
    GenLep2Id=lepIDs[ev][1];
    GenLep3Id=lepIDs[ev][2];
    GenLep4Id=lepIDs[ev][3];

    int lepIdOrdered[6] ={ GenLep1Id, GenLep2Id, GenLep3Id, GenLep4Id, 0, 0 };
    TLorentzVector p[6];
    for (int c=0; c<6; c++) p[c].SetPtEtaPhiM(event_List[ev][c][0], event_List[ev][c][1], event_List[ev][c][2], event_List[ev][c][3]);

    TLorentzVector pZ1 = p[0]+p[1];
    TLorentzVector pZ2 = p[2]+p[3];
    TLorentzVector pZZ = pZ1+pZ2;
    mzz=pZZ.M();
    m1=pZ1.M();
    m2=pZ2.M();

    cout << "Passing H mass: " << mzz << endl;

    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    mela.computeP(p[0], lepIdOrdered[0], p[1], lepIdOrdered[1], p[2], lepIdOrdered[2], p[3], lepIdOrdered[3], p0plus_VAJHU);
    mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJGG);
    mela.computeProdP(p[4], lepIdOrdered[4], p[5], lepIdOrdered[5], pZ1, 23, pZ2, 23, phjj_VAJHU);
    mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJVBF);
    mela.computeProdP(p[4], lepIdOrdered[4], p[5], lepIdOrdered[5], pZ1, 23, pZ2, 23, pvbf_VAJHU);

    cout << p0plus_VAJHU << '\t' << pvbf_VAJHU << '\t' << phjj_VAJHU << endl;
  }
}



void testME_FullMELA_MassiveLeptonFourMomenta(int flavor=0){
	int erg_tev=8;
	float mPOLE=125.0;
	float wPOLE=4.165e-3;

//	TVar::VerbosityLevel verbosity = TVar::INFO;

	Mela mela(erg_tev,mPOLE);

	float mzz = 125.; 
	float m1 = 91.471450;
	float m2 = 12.139782;
	float h1 = 0.2682896;
	float h2 = 0.1679779;
	float phi = 1.5969792;
	float hs = -0.727181;
	float phi1 = 1.8828257;
  float GenLep1Id,GenLep2Id,GenLep3Id,GenLep4Id;

  const int nEntries = 2;
	double l1_array[nEntries][4] = {
		{51.374202,	25.924766,	12.290178,	42.616376}, // Massive
		{1365.4973807340846,        10.289826593755228,        25.205694382277809,       -1365.2259480507332} // Massless
	};
	double l2_array[nEntries][4] = {
		{271.875752,	70.427173,	-11.138146,	261.769598},
		{22.786181013986834,      -0.15136300982222117,      -0.90077551414353962,       -22.767866345236371}
	};
	double l3_array[nEntries][4] = {
		{75.823478,	-16.640412,	23.246999,	70.227220},
		{1895.7562628816693,        25.837804322120054,       -28.821887970086259,       -1895.3610513294620}
	};
	double l4_array[nEntries][4] = {
		{21.481452,	9.489680,	-9.336587,	16.858699},
		{471.71918486784784,       -35.976267906053060,        4.5169691019519895,       -470.32360615864354}
	};

	for (int ev = 0; ev < nEntries; ev++){
		if (flavor == 2){
			GenLep1Id=13;
			GenLep2Id=-13;
			GenLep3Id=11;
			GenLep4Id=-11;
		}
		else if (flavor == 1){
			GenLep1Id=11;
			GenLep2Id=-11;
			GenLep3Id=11;
			GenLep4Id=-11;
		}
		else{
			GenLep1Id=13;
			GenLep2Id=-13;
			GenLep3Id=13;
			GenLep4Id=-13;
		}
		int lepIdOrdered[4] = { GenLep1Id, GenLep2Id, GenLep3Id, GenLep4Id };
		TLorentzVector p1(l1_array[ev][1],l1_array[ev][2],l1_array[ev][3],l1_array[ev][0]);
		TLorentzVector p2(l2_array[ev][1],l2_array[ev][2],l2_array[ev][3],l2_array[ev][0]);
		TLorentzVector p3(l3_array[ev][1],l3_array[ev][2],l3_array[ev][3],l3_array[ev][0]);
		TLorentzVector p4(l4_array[ev][1],l4_array[ev][2],l4_array[ev][3],l4_array[ev][0]);
    TLorentzVector pZ1 = p1+p2;
    TLorentzVector pZ2 = p3+p4;

    std::vector<TLorentzVector> pF;
    pF.push_back(p1);
    pF.push_back(p2);
    pF.push_back(p3);
    pF.push_back(p4);

    TLorentzVector p1_new = p1;
    TLorentzVector p2_new = p2;
    TLorentzVector p3_new = p3;
    TLorentzVector p4_new = p4;

    std::vector<TLorentzVector> pF_new;
    TLorentzVector pZ1_new;
    TLorentzVector pZ2_new;
    pF_new.push_back(p1_new);
    pF_new.push_back(p2_new);
    pF_new.push_back(p3_new);
    pF_new.push_back(p4_new);
/*
    for (int v = 0; v < 4; v++){
      double energy = pF_new[v].T();
      double mom = pF_new[v].P();
      double energy_new = (energy + mom) / 2.;
      if (mom != 0){
        double x_new = pF_new[v].X()*energy_new / mom;
        double y_new = pF_new[v].Y()*energy_new / mom;
        double z_new = pF_new[v].Z()*energy_new / mom;

        pF_new[v].SetXYZT(x_new, y_new, z_new, energy_new);
      }
    }
    pZ1_new = pF_new[0]+pF_new[1];
    pZ2_new = pF_new[2]+pF_new[3];
    double delta_m1=pZ1.M()-pZ1_new.M();
    double delta_m2=pZ2.M()-pZ2_new.M();
    TVector3 boostForwardZ1 = -(pZ1_new.BoostVector());
    TVector3 boostForwardZ2 = -(pZ2_new.BoostVector());
    TVector3 boostOldZ1 = -(pZ1.BoostVector());
    TVector3 boostOldZ2 = -(pZ2.BoostVector());
    for (int v = 0; v < 2; v++){
      pF_new[v].Boost(boostForwardZ1);
      double mag = pF_new[v].T();
      double mag_new = mag + delta_m1 / 2.;
      if (mag != 0) pF_new[v] *= mag_new / mag;
      pF_new[v].Boost(-boostOldZ1);
    }
    for (int v = 2; v < 4; v++){
      pF_new[v].Boost(boostForwardZ2);
      double mag = pF_new[v].T();
      double mag_new = mag + delta_m2 / 2.;
      if (mag != 0) pF_new[v] *= mag_new / mag;
      pF_new[v].Boost(-boostOldZ2);
    }
*/

    mela::constrainedRemoveLeptonMass(pF_new[0],pF_new[1]);
    mela::constrainedRemoveLeptonMass(pF_new[2],pF_new[3]);
    if(mela::forbidMassiveLeptons) cout << "Lepton masses already forbidden!" << endl;
    mela::computeAngles(pF[0], GenLep1Id,
			 pF[1], GenLep2Id,
			 pF[2], GenLep3Id,
			 pF[3], GenLep4Id,
			 hs, 
			 h1, 
			 h2, 
			 phi, 
			 phi1);
    mzz = (pF[0]+pF[1]+pF[2]+pF[3]).M();
    m1 = pZ1.M();
    m2 = pZ2.M();

    for (int v = 0; v < 4; v++){
      cout << "Vector " << v << ":\n"
        << "Initial E: " << pF[v].T() << " ... Final E: " << pF_new[v].T() << '\n'
        << "Initial X: " << pF[v].X() << " ... Final X: " << pF_new[v].X() << '\n'
        << "Initial Y: " << pF[v].Y() << " ... Final Y: " << pF_new[v].Y() << '\n'
        << "Initial Z: " << pF[v].Z() << " ... Final Z: " << pF_new[v].Z() << '\n'
        << "Initial M: " << pF[v].M() << " ... Final M: " << pF_new[v].M() << endl;
    }
    cout << "Old mZZ: " << mzz << endl;
    cout << "Old m1: " << m1 << endl;
    cout << "Old m2: " << m2 << endl;
    cout << "Old hs: " << hs << endl;
    cout << "Old h1: " << h1 << endl;
    cout << "Old h2: " << h2 << endl;
    cout << "Old phi1: " << phi1 << endl;
    cout << "Old phi: " << phi << endl;
    mela.setRemoveLeptonMasses(true);
    mela::computeAngles(pF_new[0], GenLep1Id,
			 pF_new[1], GenLep2Id,
			 pF_new[2], GenLep3Id,
			 pF_new[3], GenLep4Id,
			 hs, 
			 h1, 
			 h2, 
			 phi, 
			 phi1);
    mzz = (pF_new[0]+pF_new[1]+pF_new[2]+pF_new[3]).M();
    pZ1_new = pF_new[0]+pF_new[1];
    pZ2_new = pF_new[2]+pF_new[3];
    m1 = pZ1_new.M();
    m2 = pZ2_new.M();
    cout << "New mZZ: " << mzz << endl;
    cout << "New m1: " << m1 << endl;
    cout << "New m2: " << m2 << endl;
    cout << "New hs: " << hs << endl;
    cout << "New h1: " << h1 << endl;
    cout << "New h2: " << h2 << endl;
    cout << "New phi1: " << phi1 << endl;
    cout << "New phi: " << phi << endl;
    mela.setRemoveLeptonMasses(false);
  }

  if(mela::forbidMassiveLeptons) cout << "Lepton masses already forbidden!" << endl;
  else{
    cout << "Lepton masses not forbidden yet!" << endl;
    mela.setRemoveLeptonMasses(true);
    if (mela::forbidMassiveLeptons) cout << "Lepton masses now forbidden!" << endl;
    else cout << "Switch didn't work" << endl;
  }
  mela.setRemoveLeptonMasses(false);
  if (!mela::forbidMassiveLeptons) cout << "Lepton masses not forbidden again." << endl;
  mela::applyLeptonMassCorrection(true);
  if (mela::forbidMassiveLeptons) cout << "Lepton masses forbidden again." << endl;
  mela::applyLeptonMassCorrection(false);
}


void testME_FullMELA_FullSimMC_VBFValidation (int flavor=1){
	int erg_tev=8;
	float mPOLE=125.6;
	float wPOLE=4.15e-3;
	char TREE_NAME[] = "SelectedTree";

//	TVar::VerbosityLevel verbosity = TVar::INFO;

	Mela mela(erg_tev,mPOLE);
	TGraph* vaScale_4e = mela.vaScale_4e;
	TGraph* vaScale_4mu = mela.vaScale_4mu;
	TGraph* vaScale_2e2mu = mela.vaScale_2e2mu;
	TGraph* DggZZ_scalefactor = mela.DggZZ_scalefactor;

  TString cinput_main = "/scratch0/hep/ianderso/CJLST/140519/PRODFSR_8TeV";
  TFile* finput = new TFile(Form("%s/%s/HZZ4lTree_ZZTo%s.root", cinput_main.Data(), (flavor==0 ? "2mu2e" : "4e"), (flavor==0 ? "2e2mu" : "4e")), "read");
	TFile* foutput = new TFile(Form("HZZ4lTree_ZZTo%s_vbfMELATest.root",(flavor==0 ? "2e2mu" : "4e")),"recreate");

	TLorentzVector nullFourVector(0,0,0,0);

	float phjj_VAJHU_old;
	float pvbf_VAJHU_old;
	float phjj_VAJHU_old_NEW;
	float pvbf_VAJHU_old_NEW;
	float phjj_VAJHU_old_NEW_selfD;
	float pvbf_VAJHU_old_NEW_selfD;
	float phjj0minus_VAJHU_old_NEW;
	float pvbf0minus_VAJHU_old_NEW;
	float phjj0minus_VAJHU_old_NEW_selfD;
	float pvbf0minus_VAJHU_old_NEW_selfD;


	float jet1Pt,jet2Pt;
	float jet1px,jet1py,jet1pz,jet1E;
	float jet2px,jet2py,jet2pz,jet2E;
	float ZZPx,ZZPy,ZZPz,ZZE,dR;
	short NJets30;
	std::vector<double> * JetPt=0;
	std::vector<double> * JetEta=0;
	std::vector<double> * JetPhi=0;
	std::vector<double> * JetMass=0;
	std::vector<double> myJetPt;
	std::vector<double> myJetEta;
	std::vector<double> myJetPhi;
	std::vector<double> myJetMass;
	TBranch* bJetPt=0;
	TBranch* bJetEta=0;
	TBranch* bJetPhi=0;
	TBranch* bJetMass=0;

	float mzz = 126.; 
	float m1 = 91.471450;
	float m2 = 12.139782;
	float h1 = 0.2682896;
	float h2 = 0.1679779;
	float phi = 1.5969792;
	float hs = -0.727181;
	float phi1 = 1.8828257;
	int GenLep1Id,GenLep2Id,GenLep3Id,GenLep4Id;
  float ZZPt, ZZPhi, ZZEta;

	TTree* tree = (TTree*) finput->Get(TREE_NAME);
	tree->SetBranchAddress("NJets30",&NJets30);
	tree->SetBranchAddress("JetPt",&JetPt,&bJetPt);
	tree->SetBranchAddress("JetEta", &JetEta,&bJetEta);
	tree->SetBranchAddress("JetPhi", &JetPhi,&bJetPhi);
	tree->SetBranchAddress("JetMass", &JetMass,&bJetMass);
	tree->SetBranchAddress("phjj_VAJHU_old",&phjj_VAJHU_old);
	tree->SetBranchAddress("pvbf_VAJHU_old",&pvbf_VAJHU_old);
	tree->SetBranchAddress("ZZMass", &mzz);
  tree->SetBranchAddress("ZZPt", &ZZPt);
  tree->SetBranchAddress("ZZEta", &ZZEta);
  tree->SetBranchAddress("ZZPhi", &ZZPhi);
  tree->SetBranchAddress("Z1Mass", &m1);
	tree->SetBranchAddress("Z2Mass", &m2);
	tree->SetBranchAddress("helcosthetaZ1", &h1);
	tree->SetBranchAddress("helcosthetaZ2", &h2);
	tree->SetBranchAddress("helphi", &phi);
	tree->SetBranchAddress("costhetastar", &hs);
	tree->SetBranchAddress("phistarZ1", &phi1);

	TTree* newtree = new TTree("TestTree","");
	newtree->Branch("phjj_VAJHU_old",&phjj_VAJHU_old);
	newtree->Branch("pvbf_VAJHU_old",&pvbf_VAJHU_old);
	newtree->Branch("phjj_VAJHU_old_NEW",&phjj_VAJHU_old_NEW);
	newtree->Branch("pvbf_VAJHU_old_NEW",&pvbf_VAJHU_old_NEW);
	newtree->Branch("phjj0minus_VAJHU_old_NEW",&phjj0minus_VAJHU_old_NEW);
	newtree->Branch("pvbf0minus_VAJHU_old_NEW",&pvbf0minus_VAJHU_old_NEW);
	newtree->Branch("phjj_VAJHU_old_NEW_selfD",&phjj_VAJHU_old_NEW_selfD);
	newtree->Branch("pvbf_VAJHU_old_NEW_selfD",&pvbf_VAJHU_old_NEW_selfD);
	newtree->Branch("phjj0minus_VAJHU_old_NEW_selfD",&phjj0minus_VAJHU_old_NEW_selfD);
	newtree->Branch("pvbf0minus_VAJHU_old_NEW_selfD",&pvbf0minus_VAJHU_old_NEW_selfD);
  newtree->Branch("ZZMass", &mzz);

	int nEntries = tree->GetEntries();
	double selfDHggcoupl[SIZE_HGG][2] = { { 0 } };
	double selfDHvvcoupl[SIZE_HVV_VBF][2] = { { 0 } };
	double selfDHwwcoupl[SIZE_HWW_VBF][2] = { { 0 } };
	double ggvvcoupl[2]={0,0};
	mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
	int recorded=0;
  for (int ev = 0; ev < nEntries; ev++){
    if (recorded>=1000) break;
    tree->GetEntry(ev);
    if (flavor == 1){
      GenLep1Id=11;
      GenLep2Id=-11;
      GenLep3Id=11;
      GenLep4Id=-11;
    }
    else{
      GenLep1Id=13;
      GenLep2Id=-13;
      GenLep3Id=11;
      GenLep4Id=-11;
    }
    TLorentzVector jet1(0, 0, 1e-3, 1e-3), jet2(0, 0, 1e-3, 1e-3), higgs(0, 0, 0, 0);
    higgs.SetPtEtaPhiM(ZZPt, ZZEta, ZZPhi, mzz);

    if (JetPt->size()>=2 && NJets30>=2){
      jet1.SetPtEtaPhiM(JetPt->at(0), JetEta->at(0), JetPhi->at(0), JetMass->at(0));
      jet2.SetPtEtaPhiM(JetPt->at(1), JetEta->at(1), JetPhi->at(1), JetMass->at(1));
      mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJVBF);
      mela.computeProdP(jet1, 2, jet2, 2, higgs, 25, nullFourVector, 0, pvbf_VAJHU_old_NEW);
      mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJGG);
      mela.computeProdP(jet1, 2, jet2, 2, higgs, 25, nullFourVector, 0, phjj_VAJHU_old_NEW);
      mela.setProcess(TVar::H0minus, TVar::JHUGen, TVar::JJVBF);
      mela.computeProdP(jet1, 2, jet2, 2, higgs, 25, nullFourVector, 0, pvbf0minus_VAJHU_old_NEW);
      mela.setProcess(TVar::H0minus, TVar::JHUGen, TVar::JJGG);
      mela.computeProdP(jet1, 2, jet2, 2, higgs, 25, nullFourVector, 0, phjj0minus_VAJHU_old_NEW);

      for (int xx = 0; xx < SIZE_HVV_VBF; xx++){ for (int yy = 0; yy < 2; yy++) selfDHvvcoupl[xx][yy] = 0; }
      for (int xx = 0; xx < SIZE_HWW_VBF; xx++){ for (int yy = 0; yy < 2; yy++) selfDHwwcoupl[xx][yy] = 0; }
      for (int xx = 0; xx < SIZE_HGG; xx++){ for (int yy = 0; yy < 2; yy++) selfDHggcoupl[xx][yy] = 0; }
      selfDHvvcoupl[0][0]=1;
      selfDHwwcoupl[0][0]=1;
      selfDHggcoupl[0][0]=1;
      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::JJVBF);
      mela.computeProdP(jet1, 2, jet2, 2, higgs, 25, nullFourVector, 0, selfDHggcoupl, selfDHvvcoupl, selfDHwwcoupl, pvbf_VAJHU_old_NEW_selfD);
      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::JJGG);
      mela.computeProdP(jet1, 2, jet2, 2, higgs, 25, nullFourVector, 0, selfDHggcoupl, selfDHvvcoupl, selfDHwwcoupl, phjj_VAJHU_old_NEW_selfD);

      for (int xx = 0; xx < SIZE_HVV_VBF; xx++){ for (int yy = 0; yy < 2; yy++) selfDHvvcoupl[xx][yy] = 0; }
      for (int xx = 0; xx < SIZE_HWW_VBF; xx++){ for (int yy = 0; yy < 2; yy++) selfDHwwcoupl[xx][yy] = 0; }
      for (int xx = 0; xx < SIZE_HGG; xx++){ for (int yy = 0; yy < 2; yy++) selfDHggcoupl[xx][yy] = 0; }
      selfDHvvcoupl[3][0]=1;
      selfDHwwcoupl[3][0]=1;
      selfDHggcoupl[2][0]=1;
      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::JJVBF);
      mela.computeProdP(jet1, 2, jet2, 2, higgs, 25, nullFourVector, 0, selfDHggcoupl, selfDHvvcoupl, selfDHwwcoupl, pvbf0minus_VAJHU_old_NEW_selfD);
      mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::JJGG);
      mela.computeProdP(jet1, 2, jet2, 2, higgs, 25, nullFourVector, 0, selfDHggcoupl, selfDHvvcoupl, selfDHwwcoupl, phjj0minus_VAJHU_old_NEW_selfD);

      newtree->Fill();
      recorded++;

    }
  }


  foutput->WriteTObject(newtree);
  delete newtree;
	foutput->Close();
  finput->Close();
}

void testME_FullMELA_FullSimMC_HJValidation (int flavor=1){
	int erg_tev=8;
	float mPOLE=125.6;
	float wPOLE=4.15e-3;
	char TREE_NAME[] = "SelectedTree";

//	TVar::VerbosityLevel verbosity = TVar::INFO;

	Mela mela(erg_tev,mPOLE);
	TGraph* vaScale_4e = mela.vaScale_4e;
	TGraph* vaScale_4mu = mela.vaScale_4mu;
	TGraph* vaScale_2e2mu = mela.vaScale_2e2mu;
	TGraph* DggZZ_scalefactor = mela.DggZZ_scalefactor;

  TString cinput_main = "/scratch0/hep/ianderso/CJLST/140519/PRODFSR_8TeV";
  TFile* finput = new TFile(Form("%s/%s/HZZ4lTree_ZZTo%s.root", cinput_main.Data(), (flavor==0 ? "2mu2e" : "4e"), (flavor==0 ? "2e2mu" : "4e")), "read");
//	TFile* finput = new TFile(Form("/afs/cern.ch/work/u/usarica/HZZ4l-125p6-FullAnalysis/LHC_8TeV/%s/HZZ4lTree_ZZTo%s.root",(flavor==0 ? "2mu2e" : "4e"),(flavor==0 ? "2e2mu" : "4e")),"read");
	TFile* foutput = new TFile(Form("HZZ4lTree_ZZTo%s_hjMELATest.root",(flavor==0 ? "2e2mu" : "4e")),"recreate");

	TLorentzVector nullFourVector(0,0,0,0);

	float phj_VAJHU_old_NEW;
	float phj0minus_VAJHU_old_NEW;


	float jet1Pt,jet2Pt;
	float jet1px,jet1py,jet1pz,jet1E;
	float jet2px,jet2py,jet2pz,jet2E;
	float ZZPx,ZZPy,ZZPz,ZZE,dR;
	short NJets30;
	std::vector<double> * JetPt=0;
	std::vector<double> * JetEta=0;
	std::vector<double> * JetPhi=0;
	std::vector<double> * JetMass=0;
	std::vector<double> myJetPt;
	std::vector<double> myJetEta;
	std::vector<double> myJetPhi;
	std::vector<double> myJetMass;
	TBranch* bJetPt=0;
	TBranch* bJetEta=0;
	TBranch* bJetPhi=0;
	TBranch* bJetMass=0;

	float mzz = 126.; 
	float m1 = 91.471450;
	float m2 = 12.139782;
	float h1 = 0.2682896;
	float h2 = 0.1679779;
	float phi = 1.5969792;
	float hs = -0.727181;
	float phi1 = 1.8828257;
	int GenLep1Id,GenLep2Id,GenLep3Id,GenLep4Id;
  float ZZPt, ZZPhi, ZZEta;

	TTree* tree = (TTree*) finput->Get(TREE_NAME);
	tree->SetBranchAddress("NJets30",&NJets30);
	tree->SetBranchAddress("JetPt",&JetPt,&bJetPt);
	tree->SetBranchAddress("JetEta", &JetEta,&bJetEta);
	tree->SetBranchAddress("JetPhi", &JetPhi,&bJetPhi);
	tree->SetBranchAddress("JetMass", &JetMass,&bJetMass);
	tree->SetBranchAddress("ZZMass", &mzz);
  tree->SetBranchAddress("ZZPt", &ZZPt);
  tree->SetBranchAddress("ZZEta", &ZZEta);
  tree->SetBranchAddress("ZZPhi", &ZZPhi);
  tree->SetBranchAddress("Z1Mass", &m1);
	tree->SetBranchAddress("Z2Mass", &m2);
	tree->SetBranchAddress("helcosthetaZ1", &h1);
	tree->SetBranchAddress("helcosthetaZ2", &h2);
	tree->SetBranchAddress("helphi", &phi);
	tree->SetBranchAddress("costhetastar", &hs);
	tree->SetBranchAddress("phistarZ1", &phi1);

	TTree* newtree = new TTree("TestTree","");
	newtree->Branch("phj_VAJHU_old_NEW",&phj_VAJHU_old_NEW);
	newtree->Branch("phj0minus_VAJHU_old_NEW",&phj0minus_VAJHU_old_NEW);
  newtree->Branch("ZZMass", &mzz);


	int nEntries = tree->GetEntries();
	double selfDHggcoupl[SIZE_HGG][2] = { { 0 } };
	double selfDHvvcoupl[SIZE_HVV_VBF][2] = { { 0 } };
	double selfDHwwcoupl[SIZE_HWW_VBF][2] = { { 0 } };
	double ggvvcoupl[2]={0,0};
	mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
	int recorded=0;
  for (int ev = 0; ev < nEntries; ev++){
    if (recorded>=1000) break;
    tree->GetEntry(ev);
    if (flavor == 1){
      GenLep1Id=11;
      GenLep2Id=-11;
      GenLep3Id=11;
      GenLep4Id=-11;
    }
    else{
      GenLep1Id=13;
      GenLep2Id=-13;
      GenLep3Id=11;
      GenLep4Id=-11;
    }
    TLorentzVector jet1(0, 0, 1e-3, 1e-3), jet2(0, 0, 0, 0), higgs(0, 0, 0, 0);
    higgs.SetPtEtaPhiM(ZZPt, ZZEta, ZZPhi, mzz);

    if (JetPt->size()>=1 && NJets30>=1){
      jet1.SetPtEtaPhiM(JetPt->at(0), JetEta->at(0), JetPhi->at(0), JetMass->at(0));
      mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JH);
      mela.computeProdP(jet1, 2, jet2, 0, higgs, 25, nullFourVector, 0, phj_VAJHU_old_NEW);
      mela.setProcess(TVar::H0minus, TVar::JHUGen, TVar::JH);
      mela.computeProdP(jet1, 2, jet2, 0, higgs, 25, nullFourVector, 0, phj0minus_VAJHU_old_NEW);

      newtree->Fill();
      recorded++;
    }
  }


	foutput->WriteTObject(newtree);
	foutput->Close();
	finput->Close();
}

void testME_FullMELA_FullSimMC_VHValidation (int flavor=1){
	int erg_tev=8;
	float mPOLE=125.6;
	float wPOLE=4.15e-3;
	char TREE_NAME[] = "SelectedTree";

//	TVar::VerbosityLevel verbosity = TVar::INFO;

	Mela mela(erg_tev,mPOLE);
	TGraph* vaScale_4e = mela.vaScale_4e;
	TGraph* vaScale_4mu = mela.vaScale_4mu;
	TGraph* vaScale_2e2mu = mela.vaScale_2e2mu;
	TGraph* DggZZ_scalefactor = mela.DggZZ_scalefactor;

  TString cinput_main = "/scratch0/hep/ianderso/CJLST/140519/PRODFSR_8TeV";
  TFile* finput = new TFile(Form("%s/%s/HZZ4lTree_ZZTo%s.root", cinput_main.Data(), (flavor==0 ? "2mu2e" : "4e"), (flavor==0 ? "2e2mu" : "4e")), "read");
//	TFile* finput = new TFile(Form("/afs/cern.ch/work/u/usarica/HZZ4l-125p6-FullAnalysis/LHC_8TeV/%s/HZZ4lTree_ZZTo%s.root",(flavor==0 ? "2mu2e" : "4e"),(flavor==0 ? "2e2mu" : "4e")),"read");
	TFile* foutput = new TFile(Form("HZZ4lTree_ZZTo%s_vhMELATest.root",(flavor==0 ? "2e2mu" : "4e")),"recreate");

	TLorentzVector nullFourVector(0,0,0,0);

	float pzh_VAJHU_old_NEW;
	float pzh0minus_VAJHU_old_NEW;
	float pwh_VAJHU_old_NEW;
	float pwh0minus_VAJHU_old_NEW;


	float jet1Pt,jet2Pt;
	float jet1px,jet1py,jet1pz,jet1E;
	float jet2px,jet2py,jet2pz,jet2E;
	float ZZPx,ZZPy,ZZPz,ZZE,dR;
	short NJets30;
	std::vector<double> * JetPt=0;
	std::vector<double> * JetEta=0;
	std::vector<double> * JetPhi=0;
	std::vector<double> * JetMass=0;
	std::vector<double> myJetPt;
	std::vector<double> myJetEta;
	std::vector<double> myJetPhi;
	std::vector<double> myJetMass;
	TBranch* bJetPt=0;
	TBranch* bJetEta=0;
	TBranch* bJetPhi=0;
	TBranch* bJetMass=0;

	float mzz = 126.; 
	float m1 = 91.471450;
	float m2 = 12.139782;
	float h1 = 0.2682896;
	float h2 = 0.1679779;
	float phi = 1.5969792;
	float hs = -0.727181;
	float phi1 = 1.8828257;
	int GenLep1Id,GenLep2Id,GenLep3Id,GenLep4Id;
  float ZZPt, ZZPhi, ZZEta;

	TTree* tree = (TTree*) finput->Get(TREE_NAME);
	tree->SetBranchAddress("NJets30",&NJets30);
	tree->SetBranchAddress("JetPt",&JetPt,&bJetPt);
	tree->SetBranchAddress("JetEta", &JetEta,&bJetEta);
	tree->SetBranchAddress("JetPhi", &JetPhi,&bJetPhi);
	tree->SetBranchAddress("JetMass", &JetMass,&bJetMass);
	tree->SetBranchAddress("ZZMass", &mzz);
  tree->SetBranchAddress("ZZPt", &ZZPt);
  tree->SetBranchAddress("ZZEta", &ZZEta);
  tree->SetBranchAddress("ZZPhi", &ZZPhi);
  tree->SetBranchAddress("Z1Mass", &m1);
	tree->SetBranchAddress("Z2Mass", &m2);
	tree->SetBranchAddress("helcosthetaZ1", &h1);
	tree->SetBranchAddress("helcosthetaZ2", &h2);
	tree->SetBranchAddress("helphi", &phi);
	tree->SetBranchAddress("costhetastar", &hs);
	tree->SetBranchAddress("phistarZ1", &phi1);

	TTree* newtree = new TTree("TestTree","");
	newtree->Branch("pzh_VAJHU_old_NEW",&pzh_VAJHU_old_NEW);
	newtree->Branch("pzh0minus_VAJHU_old_NEW",&pzh0minus_VAJHU_old_NEW);
	newtree->Branch("pwh_VAJHU_old_NEW",&pwh_VAJHU_old_NEW);
	newtree->Branch("pwh0minus_VAJHU_old_NEW",&pwh0minus_VAJHU_old_NEW);
  newtree->Branch("ZZMass", &mzz);

	int nEntries = tree->GetEntries();
	double selfDHggcoupl[SIZE_HGG][2] = { { 0 } };
	double selfDHvvcoupl[SIZE_HVV_VBF][2] = { { 0 } };
	double selfDHwwcoupl[SIZE_HWW_VBF][2] = { { 0 } };
	double ggvvcoupl[2]={0,0};
	mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
	int recorded=0;
  for (int ev = 0; ev < nEntries; ev++){
    if (recorded>=1000) break;
    tree->GetEntry(ev);
    if (flavor == 1){
      GenLep1Id=11;
      GenLep2Id=-11;
      GenLep3Id=11;
      GenLep4Id=-11;
    }
    else{
      GenLep1Id=13;
      GenLep2Id=-13;
      GenLep3Id=11;
      GenLep4Id=-11;
    }
    int LepIds[4] ={
      GenLep1Id,
      GenLep2Id,
      GenLep3Id,
      GenLep4Id
    };
    TLorentzVector jet1(0, 0, 1e-3, 1e-3), jet2(0, 0, 1e-3, 1e-3), higgs(0, 0, 0, 0);
    higgs.SetPtEtaPhiM(ZZPt, ZZEta, ZZPhi, mzz);

    if (JetPt->size()>=2 && NJets30>=2){
      jet1.SetPtEtaPhiM(JetPt->at(0), JetEta->at(0), JetPhi->at(0), JetMass->at(0));
      jet2.SetPtEtaPhiM(JetPt->at(1), JetEta->at(1), JetPhi->at(1), JetMass->at(1));
      TLorentzVector myjets[2] ={ jet1, jet2 };
      TLorentzVector myHleptons[4] ={ nullFourVector, nullFourVector, nullFourVector, higgs };
      int vdaughters[2] ={ 0, 0 };
      mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZH);
      mela.computeProdP(myjets, myHleptons, vdaughters, LepIds, 0, selfDHvvcoupl, pzh_VAJHU_old_NEW);
      mela.setProcess(TVar::H0minus, TVar::JHUGen, TVar::ZH);
      mela.computeProdP(myjets, myHleptons, vdaughters, LepIds, 0, selfDHvvcoupl, pzh0minus_VAJHU_old_NEW);
      mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::WH);
      mela.computeProdP(myjets, myHleptons, vdaughters, LepIds, 0, selfDHvvcoupl, pwh_VAJHU_old_NEW);
      mela.setProcess(TVar::H0minus, TVar::JHUGen, TVar::WH);
      mela.computeProdP(myjets, myHleptons, vdaughters, LepIds, 0, selfDHvvcoupl, pwh0minus_VAJHU_old_NEW);

      newtree->Fill();
      recorded++;
    }
  }


	foutput->WriteTObject(newtree);
	foutput->Close();
	finput->Close();
}

void testME_FullMELA_FullSimMC_analyticalMELAValidation(int flavor=1){
	int erg_tev=8;
	float mPOLE=125.6;
	float wPOLE=4.15e-3;
	char TREE_NAME[] = "SelectedTree";

//	TVar::VerbosityLevel verbosity = TVar::INFO;

	Mela mela(erg_tev,mPOLE);
	TGraph* vaScale_4e = mela.vaScale_4e;
	TGraph* vaScale_4mu = mela.vaScale_4mu;
	TGraph* vaScale_2e2mu = mela.vaScale_2e2mu;
	TGraph* DggZZ_scalefactor = mela.DggZZ_scalefactor;

	TFile* finput = new TFile(Form("/afs/cern.ch/work/u/usarica/HZZ4l-125p6-FullAnalysis/LHC_8TeV/%s/HZZ4lTree_powheg15jhuGenV3-0PMH125.6.root",(flavor==0 ? "2mu2e" : "4e")),"read");
	TFile* foutput = new TFile(Form("HZZ4lTree_powheg15jhuGenV3-0PMH125.6_%s_Original_analyticalMELAValidationTestOnly.root",(flavor==0 ? "2mu2e" : "4e")),"recreate");

	float p0plus_mela;
	float p0hplus_mela;
	float p0minus_mela;
	float p0_g1prime2_mela;
	float pg1g1prime2_mela;
	float pg1g1prime2_pi2_mela;
	float pg1g2_mela;
	float pg1g2_pi2_mela;
	float pg1g4_mela;
	float pg1g4_pi2_mela;

	float p0plus_mela_NEW;
	float p0hplus_mela_NEW;
	float p0minus_mela_NEW;
	float p0_g1prime2_mela_NEW;
	float pg1g1prime2_mela_NEW;
	float pg1g1prime2_pi2_mela_NEW;
	float pg1g2_mela_NEW;
	float pg1g2_pi2_mela_NEW;
	float pg1g4_mela_NEW;
	float pg1g4_pi2_mela_NEW;


	float ggHZZ_prob_pure,ggHZZ_prob_int,ggZZ_prob_Total;
	float bkg_VAMCFM,ggzz_VAMCFM;
	float ggHZZ_prob_pure_NEW,ggHZZ_prob_int_NEW,ggZZ_prob_Total_NEW,p0plus_VAMCFM_NEW_BSMOn;
	float bkg_VAMCFM_NEW,ggzz_VAMCFM_NEW;
	float bkg_VAMCFM_STU,bkg_VAMCFM_TU,bkg_VAMCFM_S;

	float p0plus_m4l, p0plus_m4l_ScaleUp, p0plus_m4l_ScaleDown, p0plus_m4l_ResUp, p0plus_m4l_ResDown;
	float bkg_m4l, bkg_m4l_ScaleUp, bkg_m4l_ScaleDown, bkg_m4l_ResUp, bkg_m4l_ResDown;

	float mzz = 126.; 
	float m1 = 91.471450;
	float m2 = 12.139782;
	float h1 = 0.2682896;
	float h2 = 0.1679779;
	float phi = 1.5969792;
	float hs = -0.727181;
	float phi1 = 1.8828257;
	int GenLep1Id,GenLep2Id,GenLep3Id,GenLep4Id;

	TTree* tree = (TTree*) finput->Get(TREE_NAME);
	tree->SetBranchAddress("ZZMass", &mzz);
	tree->SetBranchAddress("Z1Mass", &m1);
	tree->SetBranchAddress("Z2Mass", &m2);
	tree->SetBranchAddress("helcosthetaZ1", &h1);
	tree->SetBranchAddress("helcosthetaZ2", &h2);
	tree->SetBranchAddress("helphi", &phi);
	tree->SetBranchAddress("costhetastar", &hs);
	tree->SetBranchAddress("phistarZ1", &phi1);
	tree->SetBranchAddress("Lep1ID", &GenLep1Id);
	tree->SetBranchAddress("Lep2ID", &GenLep2Id);
	tree->SetBranchAddress("Lep3ID", &GenLep3Id);
	tree->SetBranchAddress("Lep4ID", &GenLep4Id);
	tree->SetBranchAddress("pg1g2_mela",&pg1g2_mela);
	tree->SetBranchAddress("pg1g4_mela",&pg1g4_mela);
	tree->SetBranchAddress("pg1g1prime2_mela",&pg1g1prime2_mela);

	TTree* newtree = new TTree("TestTree","");
	newtree->Branch("ZZMass", &mzz);
	newtree->Branch("Z1Mass", &m1);
	newtree->Branch("Z2Mass", &m2);
	newtree->Branch("helcosthetaZ1", &h1);
	newtree->Branch("helcosthetaZ2", &h2);
	newtree->Branch("helphi", &phi);
	newtree->Branch("costhetastar", &hs);
	newtree->Branch("phistarZ1", &phi1);
	newtree->Branch("Lep1Id", &GenLep1Id);
	newtree->Branch("Lep2Id", &GenLep2Id);
	newtree->Branch("Lep3Id", &GenLep3Id);
	newtree->Branch("Lep4Id", &GenLep4Id);

	newtree->Branch("pg1g2_mela",&pg1g2_mela);
	newtree->Branch("pg1g4_mela",&pg1g4_mela);
	newtree->Branch("pg1g1prime2_mela",&pg1g1prime2_mela);

	newtree->Branch("p0plus_mela_NEW",&p0plus_mela_NEW);
	newtree->Branch("p0hplus_mela_NEW",&p0hplus_mela_NEW);
	newtree->Branch("p0minus_mela_NEW",&p0minus_mela_NEW);
	newtree->Branch("p0_g1prime2_mela_NEW",&p0_g1prime2_mela_NEW);
	newtree->Branch("pg1g2_mela_NEW",&pg1g2_mela_NEW);
	newtree->Branch("pg1g4_mela_NEW",&pg1g4_mela_NEW);
	newtree->Branch("pg1g1prime2_mela_NEW",&pg1g1prime2_mela_NEW);

	int nEntries = tree->GetEntries();
	double selfDHvvcoupl[SIZE_HVV][2] = { { 0 } };
	double ggvvcoupl[2]={0,0};
	mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
//	for (int ev = 0; ev < nEntries; ev++){
//	for (int ev = 0; ev < 100; ev++){
	for (int ev = 0; ev < 1000; ev++){
		tree->GetEntry(ev);

		int lepIdOrdered[4]={ 11,-11,11,-11 };
		if (flavor == 0){
			lepIdOrdered[0]=13;lepIdOrdered[1]=-13;
		};
		float angularOrdered[8] = { mzz, m1, m2, hs, h1, h2, phi, phi1 };
//		for(int k=0;k<8;k++) cout << angularOrdered[k] << '\t';
//		cout << endl;

		mela.setProcess(TVar::SelfDefine_spin0, TVar::ANALYTICAL, TVar::ZZGG);
		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[1][0]=0;
		selfDHvvcoupl[1][1]=0;
		selfDHvvcoupl[3][0]=0;
		selfDHvvcoupl[3][1]=0;
		selfDHvvcoupl[11][0]=0;
		selfDHvvcoupl[11][1]=0;
//		cout << "p0plus" << endl;
		p0plus_mela_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		selfDHvvcoupl[0][0]=0;
		selfDHvvcoupl[1][0]=1;
		selfDHvvcoupl[1][1]=0;
		selfDHvvcoupl[3][0]=0;
		selfDHvvcoupl[3][1]=0;
		selfDHvvcoupl[11][0]=0;
//		cout << "p0hplus" << endl;
		p0hplus_mela_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		selfDHvvcoupl[0][0]=0;
		selfDHvvcoupl[1][0]=0;
		selfDHvvcoupl[1][1]=0;
		selfDHvvcoupl[3][0]=1;
		selfDHvvcoupl[3][1]=0;
		selfDHvvcoupl[11][0]=0;
//		cout << "p0minus" << endl;
		p0minus_mela_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		selfDHvvcoupl[0][0]=0;
		selfDHvvcoupl[1][0]=0;
		selfDHvvcoupl[1][1]=0;
		selfDHvvcoupl[3][0]=0;
		selfDHvvcoupl[3][1]=0;
		selfDHvvcoupl[11][0]=-12046.01;
//		cout << "p0g1prime2" << endl;
		p0_g1prime2_mela_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);

		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[1][0]=1.638;
		selfDHvvcoupl[1][1]=0;
		selfDHvvcoupl[3][0]=0;
		selfDHvvcoupl[3][1]=0;
		selfDHvvcoupl[11][0]=0;
//		cout << "pg1g2" << endl;
		pg1g2_mela_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		pg1g2_mela_NEW -= (p0plus_mela_NEW + pow(1.638,2)*p0hplus_mela_NEW);

		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[1][0]=0;
		selfDHvvcoupl[1][1]=0;
		selfDHvvcoupl[3][0]=2.521;
		selfDHvvcoupl[3][1]=0;
		selfDHvvcoupl[11][0]=0;
//		cout << "pg1g4" << endl;
		pg1g4_mela_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		pg1g4_mela_NEW -= (p0plus_mela_NEW + pow(2.521,2)*p0minus_mela_NEW);

		selfDHvvcoupl[0][0]=1;
		selfDHvvcoupl[1][0]=0;
		selfDHvvcoupl[1][1]=0;
		selfDHvvcoupl[3][0]=0;
		selfDHvvcoupl[3][1]=0;
		selfDHvvcoupl[11][0]=12046.01;
//		cout << "pg1g1prime2" << endl;
		pg1g1prime2_mela_NEW = getJHUGenMELAWeight(mela, lepIdOrdered, angularOrdered, selfDHvvcoupl);
		pg1g1prime2_mela_NEW -= (p0plus_mela_NEW + p0_g1prime2_mela_NEW);

		newtree->Fill();
	};


	foutput->WriteTObject(newtree);
	foutput->Close();
	finput->Close();
};

void testME_FullMELA_FullSimMC_ttHValidation(int flavor=1){
  int erg_tev=8;
  float mPOLE=125.6;
  float wPOLE=4.15e-3;
  char TREE_NAME[] = "SelectedTree";

  //	TVar::VerbosityLevel verbosity = TVar::INFO;

  Mela mela(erg_tev, mPOLE);
  TString cinput_main = "/scratch0/hep/ianderso/CJLST/140519/PRODFSR_8TeV";
  TFile* finput = new TFile(Form("%s/%s/HZZ4lTree_ZZTo%s.root", cinput_main.Data(), (flavor==0 ? "2mu2e" : "4e"), (flavor==0 ? "2e2mu" : "4e")), "read");
  TFile* foutput = new TFile(Form("HZZ4lTree_ZZTo%s_tthMELATest.root", (flavor==0 ? "2e2mu" : "4e")), "recreate");

  TLorentzVector nullFourVector(0, 0, 0, 0);

  float ptth_VAJHU;


  float jet1Pt, jet2Pt;
  float jet1px, jet1py, jet1pz, jet1E;
  float jet2px, jet2py, jet2pz, jet2E;
  float ZZPx, ZZPy, ZZPz, ZZE, dR;
  short NJets30;
  std::vector<double> * JetPt=0;
  std::vector<double> * JetEta=0;
  std::vector<double> * JetPhi=0;
  std::vector<double> * JetMass=0;
  std::vector<double> myJetPt;
  std::vector<double> myJetEta;
  std::vector<double> myJetPhi;
  std::vector<double> myJetMass;
  TBranch* bJetPt=0;
  TBranch* bJetEta=0;
  TBranch* bJetPhi=0;
  TBranch* bJetMass=0;

  float mzz = 126.;
  float m1 = 91.471450;
  float m2 = 12.139782;
  float h1 = 0.2682896;
  float h2 = 0.1679779;
  float phi = 1.5969792;
  float hs = -0.727181;
  float phi1 = 1.8828257;
  int GenLep1Id, GenLep2Id, GenLep3Id, GenLep4Id;
  float ZZPt, ZZPhi, ZZEta;

  TTree* tree = (TTree*)finput->Get(TREE_NAME);
  tree->SetBranchAddress("NJets30", &NJets30);
  tree->SetBranchAddress("JetPt", &JetPt, &bJetPt);
  tree->SetBranchAddress("JetEta", &JetEta, &bJetEta);
  tree->SetBranchAddress("JetPhi", &JetPhi, &bJetPhi);
  tree->SetBranchAddress("JetMass", &JetMass, &bJetMass);
  tree->SetBranchAddress("ZZMass", &mzz);
  tree->SetBranchAddress("ZZPt", &ZZPt);
  tree->SetBranchAddress("ZZEta", &ZZEta);
  tree->SetBranchAddress("ZZPhi", &ZZPhi);
  tree->SetBranchAddress("Z1Mass", &m1);
  tree->SetBranchAddress("Z2Mass", &m2);
  tree->SetBranchAddress("helcosthetaZ1", &h1);
  tree->SetBranchAddress("helcosthetaZ2", &h2);
  tree->SetBranchAddress("helphi", &phi);
  tree->SetBranchAddress("costhetastar", &hs);
  tree->SetBranchAddress("phistarZ1", &phi1);

  TTree* newtree = new TTree("TestTree", "");
  newtree->Branch("ptth_VAJHU_new", &ptth_VAJHU);
  newtree->Branch("ZZMass", &mzz);


  int nEntries = tree->GetEntries();
  double selfDHttcoupl[SIZE_TTH][2] ={ { 0 } };
  mela.setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ttH);
  int recorded=0;
  for (int ev = 0; ev < nEntries; ev++){
    if (recorded>=1000) break;
    tree->GetEntry(ev);
    if (flavor == 1){
      GenLep1Id=11;
      GenLep2Id=-11;
      GenLep3Id=11;
      GenLep4Id=-11;
    }
    else{
      GenLep1Id=13;
      GenLep2Id=-13;
      GenLep3Id=11;
      GenLep4Id=-11;
    }
    TLorentzVector jet1(0, 0, 1e-3, 1e-3), jet2(0, 0, 1e-3, 1e-3), higgs(0, 0, 0, 0);
    higgs.SetPtEtaPhiM(ZZPt, ZZEta, ZZPhi, mzz);

    if (JetPt->size()>=2 && NJets30>=2){
      jet1.SetPtEtaPhiM(JetPt->at(0), JetEta->at(0), JetPhi->at(0), JetMass->at(0));
      jet2.SetPtEtaPhiM(JetPt->at(1), JetEta->at(1), JetPhi->at(1), JetMass->at(1));
      mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ttH);
      mela.computeProdP(jet1, 0, jet2, 0, higgs, ptth_VAJHU);

      newtree->Fill();
      recorded++;

    }
  }


  foutput->WriteTObject(newtree);
  delete newtree;
  foutput->Close();
  finput->Close();
}


void testME_FullMELA_FullMC_Parameters(){
  int erg_tev=8;
  float mPOLE=125.6;
  float wPOLE=4.15e-3;
  char TREE_NAME[] = "SelectedTree";

  Mela mela(erg_tev, mPOLE);

  //	char cinput[]="/afs/cern.ch/work/u/usarica/WidthAnalysis/LHC_8TeV/2mu2e/QQZZ/HZZ4l-125_6-8TeV-QQZZ.root";
  //	char cinput[]="/afs/cern.ch/work/u/usarica/WidthAnalysis/LHC_8TeV/2mu2e/Sig/HZZ4l-125_6-8TeV-Sig.root";
  //  char cinput[]="/afs/cern.ch/work/u/usarica/WidthAnalysis/LHC_8TeV/4e/Sig/HZZ4l-125_6-8TeV-Sig.root";
  TString cinput_main = "/scratch0/hep/ianderso/CJLST/140519/PRODFSR";
  if (erg_tev==8) cinput_main.Append("_8TeV");
//  TString cinput = Form("%s/4e/HZZ4lTree_ggTo4e_Contin-MCFM67.root", cinput_main.Data());
  TString cinput = Form("%s/4e/HZZ4lTree_ggTo4e_SMHContinInterf-MCFM67_H125.6.root", cinput_main.Data());
  char coutput[]="HZZ4l-125_6-8TeV-Sig_ParametersTest.root";
  TFile* finput = new TFile(cinput, "read");
  TFile* foutput = new TFile(coutput, "recreate");

  double p0plus_VAJHU, ggHZZ_prob_pure, ggHZZ_prob_int, ggZZ_prob_Total;
  double bkg_VAMCFM, ggzz_VAMCFM;
  double p0plus_VAJHU_NEW, ggHZZ_prob_pure_NEW, ggHZZ_prob_int_NEW, ggZZ_prob_Total_NEW;
  double bkg_VAMCFM_NEW, ggzz_VAMCFM_NEW;
  double bkg_VAMCFM_STU, bkg_VAMCFM_TU, bkg_VAMCFM_S;
  float Dgg10, Dgg10_NEW;

  double p0plus_m4l, p0plus_m4l_ScaleUp, p0plus_m4l_ScaleDown, p0plus_m4l_ResUp, p0plus_m4l_ResDown;
  double bkg_m4l, bkg_m4l_ScaleUp, bkg_m4l_ScaleDown, bkg_m4l_ResUp, bkg_m4l_ResDown;

  float mzz = 126.;
  float m1 = 91.471450;
  float m2 = 12.139782;
  float h1 = 0.2682896;
  float h2 = 0.1679779;
  float phi = 1.5969792;
  float hs = -0.727181;
  float phi1 = 1.8828257;
  float GenLep1Id, GenLep2Id, GenLep3Id, GenLep4Id;

  TTree* tree = (TTree*)finput->Get(TREE_NAME);
  tree->SetBranchAddress("ZZMass", &mzz);
  tree->SetBranchAddress("Z1Mass", &m1);
  tree->SetBranchAddress("Z2Mass", &m2);
  tree->SetBranchAddress("helcosthetaZ1", &h1);
  tree->SetBranchAddress("helcosthetaZ2", &h2);
  tree->SetBranchAddress("helphi", &phi);
  tree->SetBranchAddress("costhetastar", &hs);
  tree->SetBranchAddress("phistarZ1", &phi1);


  TTree* newtree = new TTree("TestTree", "");
  newtree->Branch("ZZMass", &mzz);
  newtree->Branch("Z1Mass", &m1);
  newtree->Branch("Z2Mass", &m2);
  newtree->Branch("helcosthetaZ1", &h1);
  newtree->Branch("helcosthetaZ2", &h2);
  newtree->Branch("helphi", &phi);
  newtree->Branch("costhetastar", &hs);
  newtree->Branch("phistarZ1", &phi1);

  newtree->Branch("bkg_VAMCFM", &bkg_VAMCFM);
  newtree->Branch("ggzz_VAMCFM", &ggzz_VAMCFM);
  newtree->Branch("ggHZZ_prob_pure", &ggHZZ_prob_pure);
  newtree->Branch("ggHZZ_prob_int", &ggHZZ_prob_int);
  newtree->Branch("Dgg10", &Dgg10);

  newtree->Branch("bkg_VAMCFM_NEW", &bkg_VAMCFM_NEW);
  newtree->Branch("ggzz_VAMCFM_NEW", &ggzz_VAMCFM_NEW);
  newtree->Branch("ggHZZ_prob_pure_NEW", &ggHZZ_prob_pure_NEW);
  newtree->Branch("ggHZZ_prob_int_NEW", &ggHZZ_prob_int_NEW);
  newtree->Branch("Dgg10_NEW", &Dgg10_NEW);

  int nEntries = tree->GetEntries();
  double selfDHvvcoupl[SIZE_HVV][2] ={ { 0. } };
  double ggvvcoupl[2]={ 0, 0 };
  //	for (int ev = 0; ev < nEntries; ev++){
  int ctr=0;
  int ev = 0;
  while (ctr<3000){
    if (ev==tree->GetEntries()) break;
    tree->GetEntry(ev);
    ev++;
//    if (mzz>200) ctr++;
    ctr++;
    int lepIdOrdered[4]={ 11, -11, 11, -11 };
    float angularOrdered[8] ={ mzz, m1, m2, hs, h1, h2, phi, phi1 };

    mela.resetMCFM_EWKParameters(1.16639E-05, 7.562468901984759e-3, 80.385, 91.1876, 0.2228972225239183);


    mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB);
    bkg_VAMCFM = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |qqZZ|**2

    mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZGG);
    ggzz_VAMCFM = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |ggZZ|**2

    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    ggHZZ_prob_pure = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |ggHZZ|**2

    mela.setProcess(TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
    ggZZ_prob_Total = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |ggZZ + ggHZZ|**2
    ggHZZ_prob_int = ggZZ_prob_Total - ggHZZ_prob_pure - ggzz_VAMCFM;

    mela.setProcess(TVar::D_gg10, TVar::MCFM, TVar::ZZGG);
    mela.computeD_gg(angularOrdered[0], angularOrdered[1], angularOrdered[2], angularOrdered[3], angularOrdered[4], angularOrdered[5], angularOrdered[6], angularOrdered[7], 1, TVar::MCFM, TVar::D_gg10, Dgg10);

    mela.resetMCFM_EWKParameters(1.16639E-05, 1./128., 80.399, 91.1876, 0.23119);


    mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB);
    bkg_VAMCFM_NEW = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |qqZZ|**2

    mela.setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZGG);
    ggzz_VAMCFM_NEW = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |ggZZ|**2

    mela.setProcess(TVar::HSMHiggs, TVar::MCFM, TVar::ZZGG);
    ggHZZ_prob_pure_NEW = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |ggHZZ|**2

    mela.setProcess(TVar::bkgZZ_SMHiggs, TVar::MCFM, TVar::ZZGG);
    ggZZ_prob_Total_NEW = getMCFMMELAWeight(mela, lepIdOrdered, angularOrdered, ggvvcoupl); // |ggZZ + ggHZZ|**2
    ggHZZ_prob_int_NEW = ggZZ_prob_Total_NEW - ggHZZ_prob_pure_NEW - ggzz_VAMCFM_NEW;

    mela.setProcess(TVar::D_gg10, TVar::MCFM, TVar::ZZGG);
    mela.computeD_gg(angularOrdered[0], angularOrdered[1], angularOrdered[2], angularOrdered[3], angularOrdered[4], angularOrdered[5], angularOrdered[6], angularOrdered[7], 1, TVar::MCFM, TVar::D_gg10, Dgg10_NEW);

    newtree->Fill();
  };


  foutput->WriteTObject(newtree);
  foutput->Close();
  finput->Close();
};


