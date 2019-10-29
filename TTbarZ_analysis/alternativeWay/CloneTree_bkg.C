//**********************************************************************************************************************************
// Remove some branches + selects the events + add variables -- for muonic channel
//***************************************** To Compile******************************************************************************
// g++ -g -std=c++11 -Wl,--no-as-needed `root-config --cflags` `root-config --libs` -lMinuit CloneTree.C -o CloneTree.exe
//**********************************************************************************************************************************

#ifndef __CINT__
#include "RooGlobalFunc.h"
//------------------------------------------------
     
#endif
#include "RooMCStudy.h"
#include "RooFitResult.h"
#include "RooStats/SPlot.h"
#include <vector>
#include <string>
#include <iostream>
#include "RooRandom.h"
#include "RooMinuit.h"
#include "TRandom3.h"
#include <time.h>
#include <TROOT.h>
#include <TH2.h>
#include <TF1.h>
#include <TTree.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TString.h>
#include <TTimeStamp.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <iostream>
#include <TMath.h>
#include "TH1D.h"
#include "TH2.h"
#include "RooFormulaVar.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooArgusBG.h"
#include "TString.h"
#include "RooFormulaVar.h"
#include "RooAddPdf.h"
#include "RooLandau.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooMappedCategory.h"
#include "RooCmdArg.h"
#include "RooChebychev.h"
#include "RooUnblindUniform.h"
#include "RooUnblindPrecision.h"
#include "RooExponential.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooSimWSTool.h"
#include "RooWorkspace.h"
#include <TLatex.h>
#include "RooFit.h"
#include "RooConstVar.h"
#include "RooSimPdfBuilder.h"
#include "RooStringVar.h"
#include "TText.h"
#include "TPaveText.h"
#include "RooHist.h"
#include "TLorentzVector.h"


using namespace std;
using namespace RooFit;
using namespace RooStats;


int main(){

  //Get old file, old tree and set top branch address
  //TFile *oldfile = new TFile("/afs/cern.ch/work/o/okara/TTbarZ/TTbarZ_analysis/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8_Summer16.root");
  //TFile *oldfile = new TFile("/user/moanwar/Run2016/CMSSW_8_0_29/src/HNL/HeavyNeutralLeptonAnalysis/test/signal/HNL_M5_mu_2.95.root");
  

  TChain *fChain = new TChain("blackJackAndHookers/blackJackAndHookersTree","");

  string inputPath = "/eos/user/o/okara/ntuples_ttbar/ntuples_ttV_2016/";
  string tree = "blackJackAndHookers/blackJackAndHookersTree";



 //string fileName = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Summer16.root" ;
  //string fileName = "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Summer16.root" ;
  //string fileName = "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Summer16.root" ;
  //string fileName = "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Summer16.root" ;
  //string fileName = "GluGluHToWWTo2L2Nu_M125_13TeV_amcatnloFXFX_pythia8_Summer16.root" ;
  //string fileName = "GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_Summer16.root" ;
  //string fileName = "GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8_Summer16.root" ;
  //string fileName = "GluGluToContinToZZTo2e2nu_13TeV_MCFM701_pythia8_Summer16.root" ;
  //string fileName = "GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8_Summer16.root" ;
  //string fileName = "GluGluToContinToZZTo2mu2nu_13TeV_MCFM701_pythia8_Summer16.root" ;
  //string fileName = "GluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8_Summer16.root" ;
  //string fileName = "GluGluToContinToZZTo4e_13TeV_MCFM701_pythia8_Summer16.root" ;
  //string fileName = "GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8_Summer16.root" ;
  //string fileName = "GluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8_Summer16.root" ;
 //string fileName = "GluGluWWTo2L2Nu_MCFM_13TeV_Summer16.root" ;
 //string fileName = "GluGluZH_HToWW_M125_13TeV_powheg_pythia8_Summer16.root" ;
  //string fileName = "ST_s-channel_4f_InclusiveDecays_13TeV-amcatnlo-pythia8_Summer16.root" ;
 //string fileName = "ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1_Summer16.root" ;
  //string fileName = "ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1_Summer16.root" ;
  //string fileName = "ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1_Summer16.root" ;
  //string fileName = "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_Summer16.root" ;
  //string fileName = "ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_Summer16.root" ;
  //string fileName = "ST_tWll_5f_LO_13TeV-MadGraph-pythia8_Summer16.root" ;
  //string fileName = "ST_tWll_5f_LO_13TeV_MadGraph_pythia8_Summer16.root" ;
  //string fileName = "ST_tWnunu_5f_LO_13TeV-MadGraph-pythia8_Summer16.root" ;
  //string fileName = "TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8_Summer16.root" ;
  //string fileName = "THQ_Hincl_13TeV-madgraph-pythia8_TuneCUETP8M1_Summer16.root" ;
  //string fileName = "THW_Hincl_13TeV-madgraph-pythia8_TuneCUETP8M1_Summer16.root" ;
  //string fileName = "TTGG_0Jets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8_Summer16.root" ;
  //string fileName = "TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8_Summer16.root" ;
  //string fileName = "TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Summer16.root" ;
  //string fileName = "TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Summer16.root" ;
  //string fileName = "TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Summer16.root" ;
  //string fileName = "TTTT_TuneCUETP8M1_13TeV-amcatnlo-pythia8_Summer16.root" ;
  //string fileName = "TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8_Summer16.root" ;
  //string fileName = "TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8_Summer16.root" ; 
  //string fileName = "TTWW_TuneCUETP8M2T4_13TeV-madgraph-pythia8_Summer16.root" ;
  //string fileName = "TTWZ_TuneCUETP8M2T4_13TeV-madgraph-pythia8_Summer16.root" ;
  // string fileName = "TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8_Summer16.root" ;
  //string fileName = "TTZToLL_M-1to10_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Summer16.root" ;
  //string fileName = "TTZZ_TuneCUETP8M2T4_13TeV-madgraph-pythia8_Summer16.root" ;
  


  
//////////////////////////////////////////////////////////////////////////////////////////////////////
////string fileName = "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_Summer16.root" ; that one gave error/// 
////string fileName = "TT_TuneCUETP8M2T4_GluonMoveCRTune_13TeV-powheg-pythia8_Summer16.root" ;////// 
////string fileName = "TT_TuneCUETP8M2T4_QCDbasedCRTune_erdON_13TeV-powheg-pythia8_Summer16.root" ;///////
///////////////////////////////////////////////////////////////////////////////////////////////////



  
  
  //string fileName = "VBFHToWWTo2L2Nu_M125_13TeV_amcatnlo_pythia8_Summer16.root" ;
  //string fileName = "VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_Summer16.root" ;
  //string fileName = "VHToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8_Summer16.root" ;
  //string fileName = "WGGJets_TuneCUETP8M1_13TeV_madgraphMLM_pythia8_Summer16.root" ;
  //string fileName = "WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Summer16.root" ;
  //string fileName = "WHiggs0PHToWW_2LFilter_M-125_13TeV-JHUGenV6_pythia8_Summer16.root" ;
  //string fileName = "WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Summer16.root" ;
  //string fileName = "WWG_TuneCUETP8M1_13TeV-amcatnlo-pythia8_Summer16.root" ;
  //string fileName = "WWTo2L2Nu_13TeV-powheg_Summer16.root" ;
  //string fileName = "WWTo2L2Nu_DoubleScattering_13TeV-pythia8_Summer16.root" ;
  //string fileName = "WWToLNuQQ_13TeV-powheg_Summer16.root" ;
  //string fileName = "WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8_Summer16.root" ;
  //string fileName = "WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8_Summer16.root" ;
  //string fileName = "WZG_TuneCUETP8M1_13TeV-amcatnlo-pythia8_Summer16.root" ;
  //string fileName = "WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8_Summer16.root" ;
  //string fileName = "WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8_Summer16.root" ;
  //string fileName = "WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_Summer16.root" ;
  //string fileName = "WZTo3LNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Summer16.root" ;
  //string fileName = "WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8_Summer16.root" ;
  //string fileName = "WZTo3LNu_mllmin01_13TeV-powheg-pythia8_ext1_Summer16.root" ;
  //string fileName = "WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8_Summer16.root" ;
  //string fileName = "WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8_Summer16.root" ;
  //string fileName = "WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8_Summer16.root" ;
  //string fileName = "ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Summer16.root" ;
  //string fileName = "ZGToLLG_01J_5f_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Summer16.root" ;
  //string fileName = "ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUgenV6_pythia8_Summer16.root" ;
  //string fileName = "ZHiggs0PHToWW_2LFilter_M-125_13TeV-JHUGenV6_pythia8_Summer16.root" ;
  //string fileName = "ZZJJTo4L_EWK_13TeV-madgraph-pythia8_Summer16.root" ;
  //string fileName = "ZZTo2L2Nu_13TeV_powheg_pythia8_Summer16.root" ;
  //string fileName = "ZZTo2L2Nu_13TeV_powheg_pythia8_ext1_Summer16.root" ;
  //string fileName = "ZZTo2L2Q_13TeV_powheg_pythia8_Summer16.root" ;
  string fileName = "ZZTo4L_13TeV_powheg_pythia8_Summer16.root" ;
  //string fileName = "ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8_Summer16.root" ;
  //string fileName = "tZq_ll_4f_13TeV-amcatnlo-herwigpp_Summer16.root" ;
 //string fileName = "tZq_ll_4f_13TeV-amcatnlo-pythia8_Summer16.root" ;
 // string fileName = "tZq_ll_4f_ckm_NLO_13TeV-amcatnlo-herwigpp_Summer16.root" ;
  //string fileName = "ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8_Summer16.root" ;
  


  fChain->Add(Form("%s/%s/%s",inputPath.c_str(),fileName.c_str(),tree.c_str()));

  string outputPath ="output";
  
  TFile *newfile = new TFile(Form("%s/%s",outputPath.c_str(),fileName.c_str()),"recreate");

  //////selection cuts

  Float_t isoCut = 0.15;
  bool is2017 = false;
  bool isMC = true;
  float Zmass = 90;

  //if I want to use a TChain.....
  //cout<< "starting..."<<endl;
  //TChain * fChain = new TChain("blackJackAndHookers/blackJackAndHookersTree",""); //okara

  //oldtree->Add("/eos/cms/store/group/phys_exotica/HNL/Data/SingleMuon/crab_Run_2016B-v3_SingleMuon/Data_Analysis10.root/HeavyNeutralLepton/tree_");

//=============================================================================================//


  //TFile *newfile = new TFile("skimmedSignale.root","recreate");
  //Create a new file + a clone of old tree in new file 
  //TTree *newtree = oldtree->CloneTree(0); 
  TTree *newtree1  = new TTree("tree_4lep","Analysis Tree");

  //cout<<"cloning done"<<endl;

  // Long64_t nentries = oldtree->GetEntries();
  Long64_t nentries = fChain->GetEntries();
  cout  <<nentries<<endl;

//======================= Old Tree Variables ==========================================// 
  // These are the variables I cut on 
   // OZGUN ADD it 


   ULong64_t       _runNb;
   ULong64_t       _lumiBlock;
   ULong64_t       _eventNb;
   UChar_t         _nVertex;
   Double_t        _met;
   Double_t        _metJECDown;
   Double_t        _metJECUp;
   Double_t        _metJetResDown;
   Double_t        _metJetResUp;
   Double_t        _metUnclDown;
   Double_t        _metUnclUp;
   Double_t        _metPhi;
   Double_t        _metPhiJECDown;
   Double_t        _metPhiJECUp;
   Double_t        _metPhiJetResDown;
   Double_t        _metPhiJetResUp;
   Double_t        _metPhiUnclDown;
   Double_t        _metPhiUnclUp;
   Bool_t          _2016_FR;
   Bool_t          _2017_FR;
   Bool_t          _HLT_Mu3_PFJet40;
   Int_t           _HLT_Mu3_PFJet40_prescale;
   Bool_t          _HLT_Mu8;
   Int_t           _HLT_Mu8_prescale;
   Bool_t          _HLT_Mu17;
   Int_t           _HLT_Mu17_prescale;
   Bool_t          _HLT_Mu27;
   Int_t           _HLT_Mu27_prescale;
   Bool_t          _HLT_Ele8_CaloIdM_TrackIdM_PFJet30;
   Int_t           _HLT_Ele8_CaloIdM_TrackIdM_PFJet30_prescale;
   Bool_t          _HLT_Ele12_CaloIdM_TrackIdM_PFJet30;
   Int_t           _HLT_Ele12_CaloIdM_TrackIdM_PFJet30_prescale;
   Bool_t          _HLT_Ele17_CaloIdM_TrackIdM_PFJet30;
   Int_t           _HLT_Ele17_CaloIdM_TrackIdM_PFJet30_prescale;
   Bool_t          _HLT_Ele23_CaloIdM_TrackIdM_PFJet30;
   Int_t           _HLT_Ele23_CaloIdM_TrackIdM_PFJet30_prescale;
   Bool_t          _passTrigger_e;
   Bool_t          _passTrigger_m;
   Bool_t          _passTrigger_ee;
   Bool_t          _passTrigger_em;
   Bool_t          _passTrigger_mm;
   Bool_t          _passTrigger_eee;
   Bool_t          _passTrigger_eem;
   Bool_t          _passTrigger_emm;
   Bool_t          _passTrigger_mmm;
   Bool_t          _passMETFilters;
   UChar_t         _nL;
   UChar_t         _nMu;
   UChar_t         _nEle;
   UChar_t         _nLight;
   Double_t        _lPt[11];   //[_nL]
   Double_t        _lEta[11];   //[_nL]
   Double_t        _lEtaSC[11];   //[_nLight]
   Double_t        _lPhi[11];   //[_nL]
   Double_t        _lE[11];   //[_nL]
   UInt_t          _lFlavor[11];   //[_nL]
   Int_t           _lCharge[11];   //[_nL]
   Double_t        _dxy[11];   //[_nL]
   Double_t        _dz[11];   //[_nL]
   Double_t        _3dIP[11];   //[_nL]
   Double_t        _3dIPSig[11];   //[_nL]
   Float_t         _lElectronMva[11];   //[_nLight]
   Float_t         _lElectronMvaHZZ[11];   //[_nLight]
   Float_t         _lElectronMvaFall17Iso[11];   //[_nLight]
   Float_t         _lElectronMvaFall17NoIso[11];   //[_nLight]
   Bool_t          _lElectronPassEmu[11];   //[_nLight]
   Bool_t          _lElectronPassConvVeto[11];   //[_nLight]
   Bool_t          _lElectronChargeConst[11];   //[_nLight]
   UInt_t          _lElectronMissingHits[11];   //[_nLight]
   Double_t        _leptonMvaSUSY[11];   //[_nLight]
   Double_t        _leptonMvaTTH[11];   //[_nLight]
   Double_t        _leptonMvatZqTTV[11];   //[_nLight]
   Bool_t          _lPOGTight[11];   //[_nL]
   Bool_t          _lPOGLoose[11];   //[_nL]
   Bool_t          _lPOGMedium[11];   //[_nL]
   Bool_t          _lPOGLooseWOIso[11];   //[_nL]
   Bool_t          _lPOGMediumWOIso[11];   //[_nL]
   Bool_t          _lPOGTightWOIso[11];   //[_nL]
   Double_t        _relIso[11];   //[_nLight]
   Double_t        _relIso0p4Mu[8];   //[_nMu]
   Double_t        _relIso0p4[11];   //[_nLight]
   Double_t        _relIso0p6[11];   //[_nLight]
   Double_t        _relIso0p8[11];   //[_nLight]
   Double_t        _relIso1p0[11];   //[_nLight]
   Double_t        _miniIso[11];   //[_nLight]
   Double_t        _miniIsoCharged[11];   //[_nLight]
   Double_t        _ptRel[11];   //[_nLight]
   Double_t        _ptRatio[11];   //[_nLight]
   Double_t        _closestJetCsvV2[11];   //[_nLight]
   Double_t        _closestJetDeepCsv_b[11];   //[_nLight]
   Double_t        _closestJetDeepCsv_bb[11];   //[_nLight]
   UInt_t          _selectedTrackMult[11];   //[_nLight]
   Double_t        _lMuonSegComp[8];   //[_nMu]
   Double_t        _lMuonTrackPt[8];   //[_nMu]
   Double_t        _lMuonTrackPtErr[8];   //[_nMu]
   UChar_t         _nJets;
   Double_t        _jetPt[20];   //[_nJets]
   Double_t        _jetSmearedPt[20];   //[_nJets]
   Double_t        _jetSmearedPt_JECDown[20];   //[_nJets]
   Double_t        _jetSmearedPt_JECUp[20];   //[_nJets]
   Double_t        _jetSmearedPt_JERDown[20];   //[_nJets]
   Double_t        _jetSmearedPt_JERUp[20];   //[_nJets]
   Double_t        _jetPt_JECUp[20];   //[_nJets]
   Double_t        _jetPt_JECDown[20];   //[_nJets]
   Double_t        _jetPt_JERUp[20];   //[_nJets]
   Double_t        _jetPt_JERDown[20];   //[_nJets]
   Double_t        _jetEta[20];   //[_nJets]
   Double_t        _jetPhi[20];   //[_nJets]
   Double_t        _jetE[20];   //[_nJets]
   Double_t        _jetCsvV2[20];   //[_nJets]
   Double_t        _jetDeepCsv_udsg[20];   //[_nJets]
   Double_t        _jetDeepCsv_b[20];   //[_nJets]
   Double_t        _jetDeepCsv_c[20];   //[_nJets]
   Double_t        _jetDeepCsv_bb[20];   //[_nJets]
   UInt_t          _jetHadronFlavor[20];   //[_nJets]
   Bool_t          _jetIsTight[20];   //[_nJets]
   Bool_t          _jetIsTightLepVeto[20];   //[_nJets]
   UChar_t         _nLheWeights;
   Double_t        _lheWeight[110];   //[_nLheWeights]
   UChar_t         _nPsWeights;
   Double_t        _psWeight[1];   //[_nPsWeights]
   Bool_t          _lIsPrompt[11];   //[_nL]
   Int_t           _lMatchPdgId[11];   //[_nL]
   Int_t           _lMomPdgId[11];   //[_nL]
   Double_t        _weight;
   Float_t         _nTrueInt;
   Double_t        _gen_met;
   Double_t        _gen_metPhi;
   UChar_t         _gen_nL;
   Double_t        _gen_lPt[20];   //[_gen_nL]
   Double_t        _gen_lEta[20];   //[_gen_nL]
   Double_t        _gen_lPhi[20];   //[_gen_nL]
   Double_t        _gen_lE[20];   //[_gen_nL]
   UInt_t          _gen_lFlavor[20];   //[_gen_nL]
   Int_t           _gen_lCharge[20];   //[_gen_nL]
   Int_t           _gen_lMomPdg[20];   //[_gen_nL]
   Bool_t          _gen_lIsPrompt[20];   //[_gen_nL]
   Double_t        _gen_partonPt[20];   //[_gen_nL]
   UInt_t          _lProvenance[11];   //[_nL]
   UInt_t          _lProvenanceCompressed[11];   //[_nL]



   fChain->SetBranchAddress("_runNb", &_runNb);
   fChain->SetBranchAddress("_lumiBlock", &_lumiBlock);
   fChain->SetBranchAddress("_eventNb", &_eventNb);
   fChain->SetBranchAddress("_nVertex", &_nVertex);
   fChain->SetBranchAddress("_met", &_met);
   fChain->SetBranchAddress("_metJECDown", &_metJECDown);
   fChain->SetBranchAddress("_metJECUp", &_metJECUp);
   fChain->SetBranchAddress("_metJetResDown", &_metJetResDown);
   fChain->SetBranchAddress("_metJetResUp", &_metJetResUp);
   fChain->SetBranchAddress("_metUnclDown", &_metUnclDown);
   fChain->SetBranchAddress("_metUnclUp", &_metUnclUp);
   fChain->SetBranchAddress("_metPhi", &_metPhi);
   fChain->SetBranchAddress("_metPhiJECDown", &_metPhiJECDown);
   fChain->SetBranchAddress("_metPhiJECUp", &_metPhiJECUp);
   fChain->SetBranchAddress("_metPhiJetResDown", &_metPhiJetResDown);
   fChain->SetBranchAddress("_metPhiJetResUp", &_metPhiJetResUp);
   fChain->SetBranchAddress("_metPhiUnclDown", &_metPhiUnclDown);
   fChain->SetBranchAddress("_metPhiUnclUp", &_metPhiUnclUp);
   fChain->SetBranchAddress("_2016_FR", &_2016_FR);
   fChain->SetBranchAddress("_2017_FR", &_2017_FR);
   fChain->SetBranchAddress("_HLT_Mu3_PFJet40", &_HLT_Mu3_PFJet40);
   fChain->SetBranchAddress("_HLT_Mu3_PFJet40_prescale", &_HLT_Mu3_PFJet40_prescale);
   fChain->SetBranchAddress("_HLT_Mu8", &_HLT_Mu8);
   fChain->SetBranchAddress("_HLT_Mu8_prescale", &_HLT_Mu8_prescale);
   fChain->SetBranchAddress("_HLT_Mu17", &_HLT_Mu17);
   fChain->SetBranchAddress("_HLT_Mu17_prescale", &_HLT_Mu17_prescale);
   fChain->SetBranchAddress("_HLT_Mu27", &_HLT_Mu27);
   fChain->SetBranchAddress("_HLT_Mu27_prescale", &_HLT_Mu27_prescale);
   fChain->SetBranchAddress("_HLT_Ele8_CaloIdM_TrackIdM_PFJet30", &_HLT_Ele8_CaloIdM_TrackIdM_PFJet30);
   fChain->SetBranchAddress("_HLT_Ele8_CaloIdM_TrackIdM_PFJet30_prescale", &_HLT_Ele8_CaloIdM_TrackIdM_PFJet30_prescale);
   fChain->SetBranchAddress("_HLT_Ele12_CaloIdM_TrackIdM_PFJet30", &_HLT_Ele12_CaloIdM_TrackIdM_PFJet30);
   fChain->SetBranchAddress("_HLT_Ele12_CaloIdM_TrackIdM_PFJet30_prescale", &_HLT_Ele12_CaloIdM_TrackIdM_PFJet30_prescale);
   fChain->SetBranchAddress("_HLT_Ele17_CaloIdM_TrackIdM_PFJet30", &_HLT_Ele17_CaloIdM_TrackIdM_PFJet30);
   fChain->SetBranchAddress("_HLT_Ele17_CaloIdM_TrackIdM_PFJet30_prescale", &_HLT_Ele17_CaloIdM_TrackIdM_PFJet30_prescale);
   fChain->SetBranchAddress("_HLT_Ele23_CaloIdM_TrackIdM_PFJet30", &_HLT_Ele23_CaloIdM_TrackIdM_PFJet30);
   fChain->SetBranchAddress("_HLT_Ele23_CaloIdM_TrackIdM_PFJet30_prescale", &_HLT_Ele23_CaloIdM_TrackIdM_PFJet30_prescale);
   fChain->SetBranchAddress("_passTrigger_e", &_passTrigger_e);
   fChain->SetBranchAddress("_passTrigger_m", &_passTrigger_m);
   fChain->SetBranchAddress("_passTrigger_ee", &_passTrigger_ee);
   fChain->SetBranchAddress("_passTrigger_em", &_passTrigger_em);
   fChain->SetBranchAddress("_passTrigger_mm", &_passTrigger_mm);
   fChain->SetBranchAddress("_passTrigger_eee", &_passTrigger_eee);
   fChain->SetBranchAddress("_passTrigger_eem", &_passTrigger_eem);
   fChain->SetBranchAddress("_passTrigger_emm", &_passTrigger_emm);
   fChain->SetBranchAddress("_passTrigger_mmm", &_passTrigger_mmm);
   fChain->SetBranchAddress("_passMETFilters", &_passMETFilters);
   fChain->SetBranchAddress("_nL", &_nL);
   fChain->SetBranchAddress("_nMu", &_nMu);
   fChain->SetBranchAddress("_nEle", &_nEle);
   fChain->SetBranchAddress("_nLight", &_nLight);
   fChain->SetBranchAddress("_lPt", &_lPt);
   fChain->SetBranchAddress("_lEta", &_lEta);
   fChain->SetBranchAddress("_lEtaSC", &_lEtaSC);
   fChain->SetBranchAddress("_lPhi", &_lPhi);
   fChain->SetBranchAddress("_lE", &_lE);
   fChain->SetBranchAddress("_lFlavor", &_lFlavor);
   fChain->SetBranchAddress("_lCharge", &_lCharge);
   fChain->SetBranchAddress("_dxy", &_dxy);
   fChain->SetBranchAddress("_dz", &_dz);
   fChain->SetBranchAddress("_3dIP", &_3dIP);
   fChain->SetBranchAddress("_3dIPSig", &_3dIPSig);
   fChain->SetBranchAddress("_lElectronMva", &_lElectronMva);
   fChain->SetBranchAddress("_lElectronMvaHZZ", &_lElectronMvaHZZ);
   fChain->SetBranchAddress("_lElectronMvaFall17Iso", &_lElectronMvaFall17Iso);
   fChain->SetBranchAddress("_lElectronMvaFall17NoIso", &_lElectronMvaFall17NoIso);
   fChain->SetBranchAddress("_lElectronPassEmu", &_lElectronPassEmu);
   fChain->SetBranchAddress("_lElectronPassConvVeto", &_lElectronPassConvVeto);
   fChain->SetBranchAddress("_lElectronChargeConst", &_lElectronChargeConst);
   fChain->SetBranchAddress("_lElectronMissingHits", &_lElectronMissingHits);
   fChain->SetBranchAddress("_leptonMvaSUSY", &_leptonMvaSUSY);
   fChain->SetBranchAddress("_leptonMvaTTH", &_leptonMvaTTH);
   fChain->SetBranchAddress("_leptonMvatZqTTV", &_leptonMvatZqTTV);
   fChain->SetBranchAddress("_lPOGTight", &_lPOGTight);
   fChain->SetBranchAddress("_lPOGLoose", &_lPOGLoose);
   fChain->SetBranchAddress("_lPOGMedium", &_lPOGMedium);
   fChain->SetBranchAddress("_lPOGLooseWOIso", &_lPOGLooseWOIso);
   fChain->SetBranchAddress("_lPOGMediumWOIso", &_lPOGMediumWOIso);
   fChain->SetBranchAddress("_lPOGTightWOIso", &_lPOGTightWOIso);
   fChain->SetBranchAddress("_relIso", &_relIso);
   fChain->SetBranchAddress("_relIso0p4Mu", &_relIso0p4Mu);
   fChain->SetBranchAddress("_relIso0p4", &_relIso0p4);
   fChain->SetBranchAddress("_relIso0p6", &_relIso0p6);
   fChain->SetBranchAddress("_relIso0p8", &_relIso0p8);
   fChain->SetBranchAddress("_relIso1p0", &_relIso1p0);
   fChain->SetBranchAddress("_miniIso", &_miniIso);
   fChain->SetBranchAddress("_miniIsoCharged", &_miniIsoCharged);
   fChain->SetBranchAddress("_ptRel", &_ptRel);
   fChain->SetBranchAddress("_ptRatio", &_ptRatio);
   fChain->SetBranchAddress("_closestJetCsvV2", &_closestJetCsvV2);
   fChain->SetBranchAddress("_closestJetDeepCsv_b", &_closestJetDeepCsv_b);
   fChain->SetBranchAddress("_closestJetDeepCsv_bb", &_closestJetDeepCsv_bb);
   fChain->SetBranchAddress("_selectedTrackMult", &_selectedTrackMult);
   fChain->SetBranchAddress("_lMuonSegComp", &_lMuonSegComp);
   fChain->SetBranchAddress("_lMuonTrackPt", &_lMuonTrackPt);
   fChain->SetBranchAddress("_lMuonTrackPtErr", &_lMuonTrackPtErr);
   fChain->SetBranchAddress("_nJets", &_nJets);
   fChain->SetBranchAddress("_jetPt", &_jetPt);
   fChain->SetBranchAddress("_jetSmearedPt", &_jetSmearedPt);
   fChain->SetBranchAddress("_jetSmearedPt_JECDown", &_jetSmearedPt_JECDown);
   fChain->SetBranchAddress("_jetSmearedPt_JECUp", &_jetSmearedPt_JECUp);
   fChain->SetBranchAddress("_jetSmearedPt_JERDown", &_jetSmearedPt_JERDown);
   fChain->SetBranchAddress("_jetSmearedPt_JERUp", &_jetSmearedPt_JERUp);
   fChain->SetBranchAddress("_jetPt_JECUp", &_jetPt_JECUp);
   fChain->SetBranchAddress("_jetPt_JECDown", &_jetPt_JECDown);
   fChain->SetBranchAddress("_jetPt_JERUp", &_jetPt_JERUp);
   fChain->SetBranchAddress("_jetPt_JERDown", &_jetPt_JERDown);
   fChain->SetBranchAddress("_jetEta", &_jetEta);
   fChain->SetBranchAddress("_jetPhi", &_jetPhi);
   fChain->SetBranchAddress("_jetE", &_jetE);
   fChain->SetBranchAddress("_jetCsvV2", &_jetCsvV2);
   fChain->SetBranchAddress("_jetDeepCsv_udsg", &_jetDeepCsv_udsg);
   fChain->SetBranchAddress("_jetDeepCsv_b", &_jetDeepCsv_b);
   fChain->SetBranchAddress("_jetDeepCsv_c", &_jetDeepCsv_c);
   fChain->SetBranchAddress("_jetDeepCsv_bb", &_jetDeepCsv_bb);
   fChain->SetBranchAddress("_jetHadronFlavor", &_jetHadronFlavor);
   fChain->SetBranchAddress("_jetIsTight", &_jetIsTight);
   fChain->SetBranchAddress("_jetIsTightLepVeto", &_jetIsTightLepVeto);
   fChain->SetBranchAddress("_nLheWeights", &_nLheWeights);
   fChain->SetBranchAddress("_lheWeight", &_lheWeight);
   fChain->SetBranchAddress("_nPsWeights", &_nPsWeights);
   fChain->SetBranchAddress("_psWeight", &_psWeight);
   fChain->SetBranchAddress("_lIsPrompt", &_lIsPrompt);
   fChain->SetBranchAddress("_lMatchPdgId", &_lMatchPdgId);
   fChain->SetBranchAddress("_lMomPdgId", &_lMomPdgId);
  /*
   fChain->SetBranchAddress("_weight", &_weight, &b__weight);
   fChain->SetBranchAddress("_nTrueInt", &b__nTrueInt);
   fChain->SetBranchAddress("_gen_met", &_gen_met, &b__gen_met);
   fChain->SetBranchAddress("_gen_metPhi", &_gen_metPhi, &b__gen_metPhi);
   fChain->SetBranchAddress("_gen_nL", &_gen_nL, &b__gen_nL);
   fChain->SetBranchAddress("_gen_lPt", _gen_lPt, &b__gen_lPt);
   fChain->SetBranchAddress("_gen_lEta", _gen_lEta, &b__gen_lEta);
   fChain->SetBranchAddress("_gen_lPhi", _gen_lPhi, &b__gen_lPhi);
   fChain->SetBranchAddress("_gen_lE", _gen_lE, &b__gen_lE);
   fChain->SetBranchAddress("_gen_lFlavor", _gen_lFlavor, &b__gen_lFlavor);
   fChain->SetBranchAddress("_gen_lCharge", _gen_lCharge, &b__gen_lCharge);
   fChain->SetBranchAddress("_gen_lMomPdg", _gen_lMomPdg, &b__gen_lMomPdg);
   fChain->SetBranchAddress("_gen_lIsPrompt", &b__gen_lIsPrompt);
   fChain->SetBranchAddress("_gen_partonPt", &b__gen_partonPt);
   fChain->SetBranchAddress"_lProvenanceCompressed", _lProvenanceCompressed, &b__lProvenanceCompressed);

 */
   //===============================Identify the new vairable okara =========================//

 // =================================  mu branches ===========================================//
   Float_t  lep_1stPt,  lep_1stEta, lep_1stPhi , lep_1stE , lep_1stCharge, lep_1stFlavor;

   TBranch* branch_lep_1stPt_tree1      = newtree1->Branch("lep_1stPt",&lep_1stPt, "lep_1stPt/F");
   TBranch* branch_lep_1stEta_tree1     = newtree1->Branch("lep_1stEta" ,&lep_1stEta ,"lep_1stEta/F");
   TBranch* branch_lep_1stPhi_tree1     = newtree1->Branch("lep_1stPhi" ,&lep_1stPhi ,"lep_1stPhi/F");
   TBranch* branch_lep_1stE_tree1       = newtree1->Branch("lep_1stE" ,&lep_1stE ,"lep_1stE/F");
   TBranch* branch_lep_1stCharge_tree1  = newtree1->Branch("lep_1stCharge" ,&lep_1stCharge ,"lep_1stCharge/F");
   TBranch* branch_lep_1stFlavor_tree1  = newtree1->Branch("lep_1stFlavor" ,&lep_1stFlavor ,"lep_1stFlavor/F");

   Float_t  lep_2ndPt,  lep_2ndEta, lep_2ndPhi , lep_2ndE , lep_2ndCharge, lep_2ndFlavor;

   TBranch* branch_lep_2ndPt_tree1      = newtree1->Branch("lep_2ndPt",&lep_2ndPt, "lep_2ndPt/F");
   TBranch* branch_lep_2ndEta_tree1     = newtree1->Branch("lep_2ndEta" ,&lep_2ndEta ,"lep_2ndEta/F");
   TBranch* branch_lep_2ndPhi_tree1     = newtree1->Branch("lep_2ndPhi" ,&lep_2ndPhi ,"lep_2ndPhi/F");
   TBranch* branch_lep_2ndE_tree1       = newtree1->Branch("lep_2ndE" ,&lep_2ndE ,"lep_2ndE/F");
   TBranch* branch_lep_2ndCharge_tree1  = newtree1->Branch("lep_2ndCharge" ,&lep_2ndCharge ,"lep_2ndCharge/F");
   TBranch* branch_lep_2ndFlavor_tree1  = newtree1->Branch("lep_2ndFlavor" ,&lep_2ndFlavor ,"lep_2ndFlavor/F");

   Float_t  lep_3rdPt,  lep_3rdEta, lep_3rdPhi , lep_3rdE , lep_3rdCharge, lep_3rdFlavor;

   TBranch* branch_lep_3rdPt_tree1      = newtree1->Branch("lep_3rdPt",&lep_3rdPt, "lep_3rdPt/F");
   TBranch* branch_lep_3rdEta_tree1     = newtree1->Branch("lep_3rdEta" ,&lep_3rdEta ,"lep_3rdEta/F");
   TBranch* branch_lep_3rdPhi_tree1     = newtree1->Branch("lep_3rdPhi" ,&lep_3rdPhi ,"lep_3rdPhi/F");
   TBranch* branch_lep_3rdE_tree1       = newtree1->Branch("lep_3rdE" ,&lep_3rdE ,"lep_3rdE/F");
   TBranch* branch_lep_3rdCharge_tree1  = newtree1->Branch("lep_3rdCharge" ,&lep_3rdCharge ,"lep_3rdCharge/F");
   TBranch* branch_lep_3rdFlavor_tree1  = newtree1->Branch("lep_3rdFlavor" ,&lep_3rdFlavor ,"lep_3rdFlavor/F");

   Float_t  lep_4thPt,  lep_4thEta, lep_4thPhi , lep_4thE , lep_4thCharge, lep_4thFlavor;

   TBranch* branch_lep_4thPt_tree1      = newtree1->Branch("lep_4thPt",&lep_4thPt, "lep_rthPt/F");
   TBranch* branch_lep_4thEta_tree1     = newtree1->Branch("lep_4thEta" ,&lep_4thEta ,"lep_4thEta/F");
   TBranch* branch_lep_4thPhi_tree1     = newtree1->Branch("lep_4thPhi" ,&lep_4thPhi ,"lep_4thPhi/F");
   TBranch* branch_lep_4thE_tree1       = newtree1->Branch("lep_4thE" ,&lep_4thE ,"lep_3rdE/F");
   TBranch* branch_lep_4thCharge_tree1  = newtree1->Branch("lep_4thCharge" ,&lep_4thCharge ,"lep_4thCharge/F");
   TBranch* branch_lep_4thFlavor_tree1  = newtree1->Branch("lep_4thFlavor" ,&lep_4thFlavor ,"lep_4thFalvor/F");

 // =================================  DiLeptons branches ===========================================//   

   Float_t DilepMass1, DilepMass2, DilepMass3, DilepMass4, DilepMass5, DilepMass6;
   Int_t Jet_Count,bJet_Count;

 TBranch* branch_DilepMass1_tree1 = newtree1->Branch("DilepMass1",&DilepMass1 ,"DilepMass1/F");
 TBranch* branch_DilepMass2_tree1 = newtree1->Branch("DilepMass2",&DilepMass2 ,"DilepMass2/F");
 TBranch* branch_DilepMass3_tree1 = newtree1->Branch("DilepMass3",&DilepMass3 ,"DilepMass3/F");
 TBranch* branch_DilepMass4_tree1 = newtree1->Branch("DilepMass4",&DilepMass4 ,"DilepMass4/F");
 TBranch* branch_DilepMass5_tree1 = newtree1->Branch("DilepMass5",&DilepMass5 ,"DilepMass5/F");
 TBranch* branch_DielpMass6_tree1 = newtree1->Branch("DilepMass6",&DilepMass6 ,"DilepMass6/F");

 TBranch* branch_Jet_Count = newtree1->Branch("Jet_Count",&Jet_Count ,"Jet_Count/I"); 
 TBranch* branch_bJet_Count = newtree1->Branch("bJet_Count",&bJet_Count ,"bJet_Count/I");

//======================= Start the running over input branches ==========================================//
 for (int i=0;i<fChain->GetEntries(); i++) {
 //for (int i=0;i<100000; i++) {

   if (i%10000==0)       cout<<i<<endl;
    fChain->GetEntry(i);
    //    if (passIsoMu24All==0 && passIsoMu27All == 0) continue;  // cut on the trigger!


    Float_t   minPt1 = -1000;
    unsigned  firstlep = -1;
    for(unsigned i=0; i < _nL ; i++){
      if(_leptonMvatZqTTV[i] < -0.4) continue;
      if (!(_lPt[i] > 40) ) continue;  
      if (_lPt[i] > minPt1){
	minPt1= _lPt[i];
	firstlep = i;	  
      }
    }

    //if( f1stlep != -1)   cout <<"1stPt =  "<<_lPt[f1stlep]<<endl;
	 
    Float_t   minPt2 = -1000;
    unsigned  secondlep = -1;
    for(unsigned i=0; i < _nL ; i++){
      if(_leptonMvatZqTTV[i] < -0.4) continue;
      if( firstlep == -1)  continue;
      if( _lPt[firstlep] > _lPt[i] and  _lPt[i] > 10) { 
	if (_lPt[i] > minPt2){
	  minPt2= _lPt[i];
	  secondlep = i;
	}
      }
    }

    //if( seclep != -1)   cout <<"2ndtPt =  "<<_lPt[seclep]<<endl;

    Float_t   minPt3 = -1000;
    unsigned  thirdlep = -1;
    for(unsigned i=0; i < _nL ; i++){
      if(_leptonMvatZqTTV[i] < -0.4) continue;
      if( secondlep == -1) continue;  
      if( _lPt[secondlep] > _lPt[i] and _lPt[i] > 10) {
	if (_lPt[i] > minPt3){
	  minPt3= _lPt[i];
	  thirdlep = i;
	}
      }
    }

    Float_t   minPt4 = -1000;
    unsigned  fourthlep = -1;
    for(unsigned i=0; i < _nL ; i++){
      if(_leptonMvatZqTTV[i] < -0.4) continue;
      if( thirdlep == -1) continue;
      if( _lPt[thirdlep] > _lPt[i] and _lPt[i] > 10) {
        if (_lPt[i] > minPt4){
          minPt4= _lPt[i];
          fourthlep = i;
        }
      }
    }

///////////////////////jet counting////////////////////////////////////

        Int_t  jetcount = 0;
	Int_t bjet = 0;
	Float_t  wp = (is2017) ? 0.4941 : 0.6324;
        for(unsigned i=0; i < _nJets ; i++){
          if(_jetSmearedPt[i] < 30 && abs(_jetEta[i]) > 2.4 && !_jetIsTight[i]) continue;
	  jetcount++;
	 if(_jetDeepCsv_b[i] > wp ) bjet++;
     }



    //if( thirdlep != -1)   cout <<"3rdPt =  "<<_lPt[thirdlep]<<endl;
    //cout<<"==================== end of the event ================"<<endl;


      TLorentzVector Lep1;
      TLorentzVector Lep2;
      TLorentzVector Lep3;
      TLorentzVector Lep4;

      float Lep_1stPt       = (firstlep != -1)  ?  _lPt[firstlep]         : -999 ;
      float Lep_1stEta      = (firstlep != -1)  ?  _lEta[firstlep]        : -999 ;
      float Lep_1stPhi      = (firstlep != -1)  ?  _lPhi[firstlep]        : -999 ;
      float Lep_1stE        = (firstlep != -1)  ?  _lE[firstlep]          : -999 ;
      float Lep_1stCharge   = (firstlep != -1)  ?  _lCharge[firstlep]     : -999 ;
      float Lep_1stFlavor   = (firstlep != -1)  ?  _lFlavor[firstlep]     : -999 ;

      float Lep_2ndPt       = (secondlep != -1) ?  _lPt[secondlep]        : -999 ;
      float Lep_2ndEta      = (secondlep != -1) ?  _lEta[secondlep]       : -999 ;
      float Lep_2ndPhi      = (secondlep != -1) ?  _lPhi[secondlep]       : -999 ;
      float Lep_2ndE        = (secondlep != -1) ?  _lE[secondlep]         : -999 ;
      float Lep_2ndCharge   = (secondlep != -1) ?  _lCharge[secondlep]    : -999 ;
      float Lep_2ndFlavor   = (secondlep != -1) ?  _lFlavor[secondlep]    : -999 ;

      float Lep_3rdPt       = (thirdlep != -1)  ? _lPt[thirdlep]          : -999 ;
      float Lep_3rdEta      = (thirdlep != -1)  ? _lEta[thirdlep]         : -999 ;
      float Lep_3rdPhi      = (thirdlep != -1)  ? _lPhi[thirdlep]         : -999 ;
      float Lep_3rdE        = (thirdlep != -1)  ? _lE[thirdlep]           : -999 ;
      float Lep_3rdCharge   = (thirdlep != -1)  ? _lCharge[thirdlep]      : -999 ;
      float Lep_3rdFlavor   = (thirdlep != -1)  ? _lFlavor[thirdlep]      : -999 ;
 
      float Lep_4thPt       = (fourthlep != -1) ? _lPt[fourthlep]         : -999 ;
      float Lep_4thEta      = (fourthlep != -1) ? _lEta[fourthlep]        : -999 ;
      float Lep_4thPhi      = (fourthlep != -1) ? _lPhi[fourthlep]        : -999 ;
      float Lep_4thE        = (fourthlep != -1) ? _lE[fourthlep]          : -999 ;
      float Lep_4thCharge   = (fourthlep != -1) ? _lCharge[fourthlep]     : -999 ;
      float Lep_4thFlavor   = (fourthlep != -1) ? _lFlavor[fourthlep]     : -999 ;

      Lep1.SetPtEtaPhiE(Lep_1stPt,Lep_1stEta,Lep_1stPhi,Lep_1stE);
      Lep2.SetPtEtaPhiE(Lep_2ndPt,Lep_2ndEta,Lep_2ndPhi,Lep_2ndE);
      Lep3.SetPtEtaPhiE(Lep_3rdPt,Lep_3rdEta,Lep_3rdPhi,Lep_3rdE);
      Lep4.SetPtEtaPhiE(Lep_4thPt,Lep_4thEta,Lep_4thPhi,Lep_4thE);

      float DiLepMass1  = (Lep1 + Lep2).M();
      float DiLepMass2  = (Lep1 + Lep3).M();
      float DiLepMass3  = (Lep1 + Lep4).M();
      float DiLepMass4  = (Lep2 + Lep3).M();
      float DiLepMass5  = (Lep2 + Lep4).M();
      float DiLepMass6  = (Lep3 + Lep4).M();


      if(firstlep != -1 && secondlep != -1 && thirdlep != -1 && fourthlep != -1){

      if(   ((Lep_1stCharge + Lep_2ndCharge) == 0  && abs(DiLepMass1 - Zmass) < 20 && Lep_1stFlavor == Lep_2ndFlavor  and 
	     (Lep_2ndCharge + Lep_3rdCharge) == 0  && abs(DiLepMass4 - Zmass) < 20 && Lep_2ndFlavor == Lep_3rdFlavor) 


	 or ((Lep_1stCharge + Lep_3rdCharge) == 0  && abs(DiLepMass2 - Zmass) < 20 && Lep_1stFlavor == Lep_3rdFlavor and 
	     (Lep_2ndCharge + Lep_4thCharge) == 0  && abs(DiLepMass5 - Zmass) < 20 && Lep_2ndFlavor == Lep_4thFlavor )	 


	 or ((Lep_1stCharge + Lep_4thCharge) == 0  && abs(DiLepMass3 - Zmass) < 20   && Lep_1stFlavor == Lep_4thFlavor and
	     (Lep_3rdCharge + Lep_4thCharge) == 0  && abs(DiLepMass6 - Zmass) < 20 ) && Lep_4thFlavor == Lep_3rdFlavor)   


      {

	lep_1stPt     = _lPt[firstlep];
        lep_1stEta    = _lEta[firstlep] ;
        lep_1stPhi    = _lPhi[firstlep];
        lep_1stE      = _lE[firstlep];
	lep_1stCharge = _lCharge[firstlep] ;
        lep_1stFlavor = _lFlavor[firstlep] ;


	lep_2ndPt     = _lPt[secondlep];
        lep_2ndEta    = _lEta[secondlep] ;
        lep_2ndPhi    = _lPhi[secondlep];
        lep_2ndE      = _lE[secondlep] ;
        lep_2ndCharge = _lCharge[secondlep] ;
        lep_2ndFlavor = _lFlavor[secondlep] ;

	lep_3rdPt     = _lPt[thirdlep];
	lep_3rdEta    = _lEta[thirdlep] ;
	lep_3rdPhi    = _lPhi[thirdlep];
	lep_3rdE      = _lE[thirdlep] ;
        lep_3rdCharge = _lCharge[thirdlep] ;
        lep_3rdFlavor = _lFlavor[thirdlep] ;

        lep_4thPt     = _lPt[fourthlep] ;
        lep_4thEta    = _lEta[fourthlep] ;
        lep_4thPhi    = _lPhi[fourthlep];
        lep_4thE      = _lE[fourthlep] ;
        lep_4thCharge = _lCharge[fourthlep] ;
        lep_4thFlavor = _lFlavor[fourthlep] ;


	DilepMass1 = DiLepMass1;
	DilepMass2 = DiLepMass2 ;
	DilepMass3 = DiLepMass3 ;
	DilepMass4 = DiLepMass4;
        DilepMass5 = DiLepMass5 ;
	DilepMass6 = DiLepMass6 ;

	Jet_Count = jetcount;
	bJet_Count = bjet;
	newtree1->Fill();
      }
      }
 }

 newfile->Write();  
 newtree1->Print();
 newtree1->AutoSave();

  delete newfile;

  return 0;
}
