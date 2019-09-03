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
  TFile *oldfile = new TFile("/eos/user/m/moanwar/Okara/CMSSW_9_4_13/src/ana/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8_Summer16.root");
  //TFile *oldfile = new TFile("/user/moanwar/Run2016/CMSSW_8_0_29/src/HNL/HeavyNeutralLeptonAnalysis/test/signal/HNL_M5_mu_2.95.root");
  TTree *fChain = (TTree*)oldfile->Get("blackJackAndHookers/blackJackAndHookersTree");


  //////selection cuts

  Float_t isoCut = 0.15;
  bool isMC = true;
  float Zmass = 90;

  //if I want to use a TChain.....
  //cout<< "starting..."<<endl;
  //TChain * fChain = new TChain("blackJackAndHookers/blackJackAndHookersTree",""); //okara

  //oldtree->Add("/eos/cms/store/group/phys_exotica/HNL/Data/SingleMuon/crab_Run_2016B-v3_SingleMuon/Data_Analysis10.root/HeavyNeutralLepton/tree_");

//=============================================================================================//


  char fileName[256];
  cout<<"Please enter the name of the output root file you want to create (yyy.root) : "<<endl;
  cin.getline(fileName,256);
  TFile *newfile = new TFile(fileName,"recreate");

  //TFile *newfile = new TFile("skimmedSignale.root","recreate");
  //Create a new file + a clone of old tree in new file 
  //TTree *newtree = oldtree->CloneTree(0); 
  TTree *newtree1  = new TTree("tree_4mu","Analysis Tree");
  TTree *newtree2  = new TTree("tree_3mu1ele","Analysis Tree");
  TTree *newtree3  = new TTree("tree_2mu2ele","Analysis Tree");
  TTree *newtree4  = new TTree("tree_1mu3ele","Analysis Tree");
  TTree *newtree5  = new TTree("tree_4ele","Analysis Tree");

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
 Float_t  mu_1stPt,  mu_1stEta, mu_1stPhi , mu_1stE , mu_1stCharge;

 TBranch* branch_mu_1stPt_tree1      = newtree1->Branch("mu_1stPt",&mu_1stPt, "mu_1stPt/F");
 TBranch* branch_mu_1stEta_tree1     = newtree1->Branch("mu_1stEta" ,&mu_1stEta ,"mu_1stEta/F");
 TBranch* branch_mu_1stPhi_tree1     = newtree1->Branch("mu_1stPhi" ,&mu_1stPhi ,"mu_1stPhi/F");
 TBranch* branch_mu_1stE_tree1       = newtree1->Branch("mu_1stE" ,&mu_1stE ,"mu_1stE/F");
 TBranch* branch_mu_1stCharge_tree1  = newtree1->Branch("mu_1stCharge" ,&mu_1stCharge ,"mu_1stCharge/F");

 TBranch* branch_mu_1stPt_tree2      = newtree2->Branch("mu_1stPt",&mu_1stPt, "mu_1stPt/F");
 TBranch* branch_mu_1stEta_tree2     = newtree2->Branch("mu_1stEta" ,&mu_1stEta ,"mu_1stEta/F");
 TBranch* branch_mu_1stPhi_tree2     = newtree2->Branch("mu_1stPhi" ,&mu_1stPhi ,"mu_1stPhi/F");
 TBranch* branch_mu_1stE_tree2       = newtree2->Branch("mu_1stE" ,&mu_1stE ,"mu_1stE/F");
 TBranch* branch_mu_1stCharge_tree2  = newtree2->Branch("mu_1stCharge" ,&mu_1stCharge ,"mu_1stCharge/F");

 TBranch* branch_mu_1stPt_tree3      = newtree3->Branch("mu_1stPt",&mu_1stPt, "mu_1stPt/F");
 TBranch* branch_mu_1stEta_tree3     = newtree3->Branch("mu_1stEta" ,&mu_1stEta ,"mu_1stEta/F");
 TBranch* branch_mu_1stPhi_tree3     = newtree3->Branch("mu_1stPhi" ,&mu_1stPhi ,"mu_1stPhi/F");
 TBranch* branch_mu_1stE_tree3       = newtree3->Branch("mu_1stE" ,&mu_1stE ,"mu_1stE/F");
 TBranch* branch_mu_1stCharge_tree3  = newtree3->Branch("mu_1stCharge" ,&mu_1stCharge ,"mu_1stCharge/F");

 TBranch* branch_mu_1stPt_tree4      = newtree4->Branch("mu_1stPt",&mu_1stPt, "mu_1stPt/F");
 TBranch* branch_mu_1stEta_tree4     = newtree4->Branch("mu_1stEta" ,&mu_1stEta ,"mu_1stEta/F");
 TBranch* branch_mu_1stPhi_tree4     = newtree4->Branch("mu_1stPhi" ,&mu_1stPhi ,"mu_1stPhi/F");
 TBranch* branch_mu_1stE_tree4       = newtree4->Branch("mu_1stE" ,&mu_1stE ,"mu_1stE/F");
 TBranch* branch_mu_1stCharge_tree4  = newtree4->Branch("mu_1stCharge" ,&mu_1stCharge ,"mu_1stCharge/F");


 Float_t  mu_2ndPt,  mu_2ndEta, mu_2ndPhi , mu_2ndE , mu_2ndCharge;

 TBranch* branch_mu_2ndPt_tree1      = newtree1->Branch("mu_2ndPt",&mu_2ndPt, "mu_2ndPt/F");
 TBranch* branch_mu_2ndEta_tree1     = newtree1->Branch("mu_2ndEta" ,&mu_2ndEta ,"mu_2ndEta/F");
 TBranch* branch_mu_2ndPhi_tree1     = newtree1->Branch("mu_2ndPhi" ,&mu_2ndPhi ,"mu_2ndPhi/F");
 TBranch* branch_mu_2ndE_tree1       = newtree1->Branch("mu_2ndE" ,&mu_2ndE ,"mu_2ndE/F");
 TBranch* branch_mu_2ndCharge_tree1  = newtree1->Branch("mu_2ndCharge" ,&mu_2ndCharge ,"mu_2ndCharge/F");

 TBranch* branch_mu_2ndPt_tree2       = newtree2->Branch("mu_2ndPt",&mu_2ndPt, "mu_2ndPt/F");
 TBranch* branch_mu_2ndEta_tree2      = newtree2->Branch("mu_2ndEta" ,&mu_2ndEta ,"mu_2ndEta/F");
 TBranch* branch_mu_2ndPhi_tree2      = newtree2->Branch("mu_2ndPhi" ,&mu_2ndPhi ,"mu_2ndPhi/F");
 TBranch* branch_mu_2ndE_tree2        = newtree2->Branch("mu_2ndE" ,&mu_2ndE ,"mu_2ndE/F");
 TBranch* branch_mu_2ndCharge_tree2   = newtree2->Branch("mu_2ndCharge" ,&mu_2ndCharge ,"mu_2ndCharge/F");

 TBranch* branch_mu_2ndPt_tree3      = newtree3->Branch("mu_2ndPt",&mu_2ndPt, "mu_2ndPt/F");
 TBranch* branch_mu_2ndEta_tree3     = newtree3->Branch("mu_2ndEta" ,&mu_2ndEta ,"mu_2ndEta/F");
 TBranch* branch_mu_2ndPhi_tree3     = newtree3->Branch("mu_2ndPhi" ,&mu_2ndPhi ,"mu_2ndPhi/F");
 TBranch* branch_mu_2ndE_tree3       = newtree3->Branch("mu_2ndE" ,&mu_2ndE ,"mu_2ndE/F");
 TBranch* branch_mu_2ndCharge_tree3  = newtree3->Branch("mu_2ndCharge" ,&mu_2ndCharge ,"mu_2ndCharge/F");


 Float_t  mu_3rdPt,  mu_3rdEta, mu_3rdPhi , mu_3rdE , mu_3rdCharge;
 
 TBranch* branch_mu_3rdPt_tree1      = newtree1->Branch("mu_3rdPt",&mu_3rdPt, "mu_3rdPt/F");
 TBranch* branch_mu_3rdEta_tree1     = newtree1->Branch("mu_3rdEta" ,&mu_3rdEta ,"mu_3rdEta/F");
 TBranch* branch_mu_3rdPhi_tree1     = newtree1->Branch("mu_3rdPhi" ,&mu_3rdPhi ,"mu_3rdPhi/F");
 TBranch* branch_mu_3rdE_tree1       = newtree1->Branch("mu_3rdE" ,&mu_3rdE ,"mu_3rdE/F");
 TBranch* branch_mu_3rdCharge_tree1  = newtree1->Branch("mu_3rdCharge" ,&mu_3rdCharge ,"mu_3rdCharge/F");

 TBranch* branch_mu_3rdPt_tree2       = newtree2->Branch("mu_3rdPt",&mu_3rdPt, "mu_3rdPt/F");
 TBranch* branch_mu_3rdEta_tree2      = newtree2->Branch("mu_3rdEta" ,&mu_3rdEta ,"mu_3rdEta/F");
 TBranch* branch_mu_3rdPhi_tree2      = newtree2->Branch("mu_3rdPhi" ,&mu_3rdPhi ,"mu_3rdPhi/F");
 TBranch* branch_mu_3rdE_tree2        = newtree2->Branch("mu_3rdE" ,&mu_3rdE ,"mu_3rdE/F");
 TBranch* branch_mu_3rdCharge_tree2   = newtree2->Branch("mu_3rdCharge" ,&mu_3rdCharge ,"mu_3rdCharge/F");

 Float_t  mu_4thPt,  mu_4thEta, mu_4thPhi , mu_4thE , mu_4thCharge;

 TBranch* branch_mu_4thPt_tree1      = newtree1->Branch("mu_4thPt",&mu_4thPt, "mu_rthPt/F");
 TBranch* branch_mu_4thEta_tree1     = newtree1->Branch("mu_4thEta" ,&mu_4thEta ,"mu_4thEta/F");
 TBranch* branch_mu_4thPhi_tree1     = newtree1->Branch("mu_4thPhi" ,&mu_4thPhi ,"mu_4thPhi/F");
 TBranch* branch_mu_4thE_tree1       = newtree1->Branch("mu_4thE" ,&mu_4thE ,"mu_3rdE/F");
 TBranch* branch_mu_4thCharge_tree1  = newtree1->Branch("mu_4thCharge" ,&mu_4thCharge ,"mu_4thCharge/F");

 // =================================  ele branches ===========================================//   

 Float_t  ele_1stPt,  ele_1stEta, ele_1stPhi , ele_1stE , ele_1stCharge;

 TBranch* branch_ele_1stPt_tree2      = newtree2->Branch("ele_1stPt",&ele_1stPt, "ele_1stPt/F");
 TBranch* branch_ele_1stEta_tree2     = newtree2->Branch("ele_1stEta" ,&ele_1stEta ,"ele_1stEta/F");
 TBranch* branch_ele_1stPhi_tree2     = newtree2->Branch("ele_1stPhi" ,&ele_1stPhi ,"ele_1stPhi/F");
 TBranch* branch_ele_1stE_tree2       = newtree2->Branch("ele_1stE" ,&ele_1stE ,"ele_1stE/F");
 TBranch* branch_ele_1stCharge_tree2  = newtree2->Branch("ele_1stCharge" ,&ele_1stCharge ,"ele_1stCharge/F");

 TBranch* branch_ele_1stPt_tree3      = newtree3->Branch("ele_1stPt",&ele_1stPt, "ele_1stPt/F");
 TBranch* branch_ele_1stEta_tree3     = newtree3->Branch("ele_1stEta" ,&ele_1stEta ,"ele_1stEta/F");
 TBranch* branch_ele_1stPhi_tree3     = newtree3->Branch("ele_1stPhi" ,&ele_1stPhi ,"ele_1stPhi/F");
 TBranch* branch_ele_1stE_tree3       = newtree3->Branch("ele_1stE" ,&ele_1stE ,"ele_1stE/F");
 TBranch* branch_ele_1stCharge_tree3  = newtree3->Branch("ele_1stCharge" ,&ele_1stCharge ,"ele_1stCharge/F");

 TBranch* branch_ele_1stPt_tree4      = newtree4->Branch("ele_1stPt",&ele_1stPt, "ele_1stPt/F");
 TBranch* branch_ele_1stEta_tree4     = newtree4->Branch("ele_1stEta" ,&ele_1stEta ,"ele_1stEta/F");
 TBranch* branch_ele_1stPhi_tree4     = newtree4->Branch("ele_1stPhi" ,&ele_1stPhi ,"ele_1stPhi/F");
 TBranch* branch_ele_1stE_tree4       = newtree4->Branch("ele_1stE" ,&ele_1stE ,"ele_1stE/F");
 TBranch* branch_ele_1stCharge_tree4  = newtree4->Branch("ele_1stCharge" ,&ele_1stCharge ,"ele_1stCharge/F");

 TBranch* branch_ele_1stPt_tree5      = newtree5->Branch("ele_1stPt",&ele_1stPt, "ele_1stPt/F");
 TBranch* branch_ele_1stEta_tree5     = newtree5->Branch("ele_1stEta" ,&ele_1stEta ,"ele_1stEta/F");
 TBranch* branch_ele_1stPhi_tree5     = newtree5->Branch("ele_1stPhi" ,&ele_1stPhi ,"ele_1stPhi/F");
 TBranch* branch_ele_1stE_tree5       = newtree5->Branch("ele_1stE" ,&ele_1stE ,"ele_1stE/F");
 TBranch* branch_ele_1stCharge_tree5  = newtree5->Branch("ele_1stCharge" ,&ele_1stCharge ,"ele_1stCharge/F");


 Float_t  ele_2ndPt,  ele_2ndEta, ele_2ndPhi , ele_2ndE , ele_2ndCharge;

 TBranch* branch_ele_2ndPt_tree3      = newtree3->Branch("ele_2ndPt",&ele_2ndPt, "ele_2ndPt/F");
 TBranch* branch_ele_2ndEta_tree3     = newtree3->Branch("ele_2ndEta" ,&ele_2ndEta ,"ele_2ndEta/F");
 TBranch* branch_ele_2ndPhi_tree3     = newtree3->Branch("ele_2ndPhi" ,&ele_2ndPhi ,"ele_2ndPhi/F");
 TBranch* branch_ele_2ndE_tree3       = newtree3->Branch("ele_2ndE" ,&ele_2ndE ,"ele_2ndE/F");
 TBranch* branch_ele_2ndCharge_tree3  = newtree3->Branch("ele_2ndCharge" ,&ele_2ndCharge ,"ele_2ndCharge/F");

 TBranch* branch_ele_2ndPt_tree4       = newtree4->Branch("ele_2ndPt",&ele_2ndPt, "ele_2ndPt/F");
 TBranch* branch_ele_2ndEta_tree4      = newtree4->Branch("ele_2ndEta" ,&ele_2ndEta ,"ele_2ndEta/F");
 TBranch* branch_ele_2ndPhi_tree4      = newtree4->Branch("ele_2ndPhi" ,&ele_2ndPhi ,"ele_2ndPhi/F");
 TBranch* branch_ele_2ndE_tree4        = newtree4->Branch("ele_2ndE" ,&ele_2ndE ,"ele_2ndE/F");
 TBranch* branch_ele_2ndCharge_tree4   = newtree4->Branch("ele_2ndCharge" ,&ele_2ndCharge ,"ele_2ndCharge/F");

 TBranch* branch_ele_2ndPt_tree5      = newtree5->Branch("ele_2ndPt",&ele_2ndPt, "ele_2ndPt/F");
 TBranch* branch_ele_2ndEta_tree5     = newtree5->Branch("ele_2ndEta" ,&ele_2ndEta ,"ele_2ndEta/F");
 TBranch* branch_ele_2ndPhi_tree5     = newtree5->Branch("ele_2ndPhi" ,&ele_2ndPhi ,"ele_2ndPhi/F");
 TBranch* branch_ele_2ndE_tree5       = newtree5->Branch("ele_2ndE" ,&ele_2ndE ,"ele_2ndE/F");
 TBranch* branch_ele_2ndCharge_tree5  = newtree5->Branch("ele_2ndCharge" ,&ele_2ndCharge ,"ele_2ndCharge/F");


 Float_t  ele_3rdPt,  ele_3rdEta, ele_3rdPhi , ele_3rdE , ele_3rdCharge;
 
 TBranch* branch_ele_3rdPt_tree4      = newtree4->Branch("ele_3rdPt",&ele_3rdPt, "ele_3rdPt/F");
 TBranch* branch_ele_3rdEta_tree4     = newtree4->Branch("ele_3rdEta" ,&ele_3rdEta ,"ele_3rdEta/F");
 TBranch* branch_ele_3rdPhi_tree4     = newtree4->Branch("ele_3rdPhi" ,&ele_3rdPhi ,"ele_3rdPhi/F");
 TBranch* branch_ele_3rdE_tree4       = newtree4->Branch("ele_3rdE" ,&ele_3rdE ,"ele_3rdE/F");
 TBranch* branch_ele_3rdCharge_tree4  = newtree4->Branch("ele_3rdCharge" ,&ele_3rdCharge ,"ele_3rdCharge/F");

 TBranch* branch_ele_3rdPt_tree5       = newtree5->Branch("ele_3rdPt",&ele_3rdPt, "ele_3rdPt/F");
 TBranch* branch_ele_3rdEta_tree5      = newtree5->Branch("ele_3rdEta" ,&ele_3rdEta ,"ele_3rdEta/F");
 TBranch* branch_ele_3rdPhi_tree5      = newtree5->Branch("ele_3rdPhi" ,&ele_3rdPhi ,"ele_3rdPhi/F");
 TBranch* branch_ele_3rdE_tree5        = newtree5->Branch("ele_3rdE" ,&ele_3rdE ,"ele_3rdE/F");
 TBranch* branch_ele_3rdCharge_tree5   = newtree5->Branch("ele_3rdCharge" ,&ele_3rdCharge ,"ele_3rdCharge/F");

 Float_t  ele_4thPt,  ele_4thEta, ele_4thPhi , ele_4thE , ele_4thCharge;

 TBranch* branch_ele_4thPt_tree5      = newtree5->Branch("ele_4thPt",&ele_4thPt, "ele_rthPt/F");
 TBranch* branch_ele_4thEta_tree5     = newtree5->Branch("ele_4thEta" ,&ele_4thEta ,"ele_4thEta/F");
 TBranch* branch_ele_4thPhi_tree5     = newtree5->Branch("ele_4thPhi" ,&ele_4thPhi ,"ele_4thPhi/F");
 TBranch* branch_ele_4thE_tree5       = newtree5->Branch("ele_4thE" ,&ele_4thE ,"ele_3rdE/F");
 TBranch* branch_ele_4thCharge_tree5  = newtree5->Branch("ele_4thCharge" ,&ele_4thCharge ,"ele_4thCharge/F");

 // =================================  DiLeptons branches ===========================================//   

 Float_t DimuMass1, DimuMass2, DimuMass3, DimuMass4, DimuMass5, DimuMass6;

 Float_t DieleMass1, DieleMass2, DieleMass3, DieleMass4, DieleMass5, DieleMass6;

 TBranch* branch_DimuMass1_tree1 = newtree1->Branch("DimuMass1",&DimuMass1 ,"DimuMass1/F");
 TBranch* branch_DimuMass2_tree1 = newtree1->Branch("DimuMass2",&DimuMass2 ,"DimuMass2/F");
 TBranch* branch_DimuMass3_tree1 = newtree1->Branch("DimuMass3",&DimuMass3 ,"DimuMass3/F");
 TBranch* branch_DimuMass4_tree1 = newtree1->Branch("DimuMass4",&DimuMass4 ,"DimuMass4/F");
 TBranch* branch_DimuMass5_tree1 = newtree1->Branch("DimuMass5",&DimuMass5 ,"DimuMass5/F");
 TBranch* branch_DimuMass6_tree1 = newtree1->Branch("DimuMass6",&DimuMass6 ,"DimuMass6/F");

 TBranch* branch_DimuMass1_tree2 = newtree2->Branch("DimuMass1",&DimuMass1 ,"DimuMass1/F");
 TBranch* branch_DimuMass2_tree2 = newtree2->Branch("DimuMass2",&DimuMass2 ,"DimuMass2/F");
 TBranch* branch_DimuMass3_tree2 = newtree2->Branch("DimuMass3",&DimuMass3 ,"DimuMass3/F");

 TBranch* branch_DimuMass1_tree3  = newtree3->Branch("DimuMass1",&DimuMass1 ,"DimuMass1/F");
 TBranch* branch_DieleMass1_tree3 = newtree3->Branch("DieleMass1",&DieleMass1 ,"DieleMass1/F");

 TBranch* branch_DieleMass1_tree4 = newtree4->Branch("DieleMass1",&DieleMass1 ,"DieleMass1/F");
 TBranch* branch_DieleMass2_tree4 = newtree4->Branch("DieleMass2",&DieleMass2 ,"DieleMass2/F");
 TBranch* branch_DieleMass3_tree4 = newtree4->Branch("DieleMass3",&DieleMass3 ,"DieleMass3/F");

 TBranch* branch_DieleMass1_tree5 = newtree5->Branch("DieleMass1",&DieleMass1 ,"DieleMass1/F");
 TBranch* branch_DieleMass2_tree5 = newtree5->Branch("DieleMass2",&DieleMass2 ,"DieleMass2/F");
 TBranch* branch_DieleMass3_tree5 = newtree5->Branch("DieleMass3",&DieleMass3 ,"DieleMass3/F");
 TBranch* branch_DieleMass4_tree5 = newtree5->Branch("DieleMass4",&DieleMass4 ,"DieleMass4/F");
 TBranch* branch_DieleMass5_tree5 = newtree5->Branch("DieleMass5",&DieleMass5 ,"DieleMass5/F");
 TBranch* branch_DieleMass6_tree5 = newtree5->Branch("DieleMass6",&DieleMass6 ,"DieleMass6/F");

//======================= Start the running over input branches ==========================================//
 //for (int i=0;i<fChain->GetEntries(); i++) {
 for (int i=0;i<100000; i++) {

   if (i%10000==0)       cout<<i<<endl;
    fChain->GetEntry(i);
    //    if (passIsoMu24All==0 && passIsoMu27All == 0) continue;  // cut on the trigger!


    Float_t   minPt1 = -1000;
    unsigned  firstMu = -1;
    for(unsigned i=0; i < _nL ; i++){
      if(_lFlavor[i] != 1 or _leptonMvatZqTTV[i] < -0.4) continue;
      if (!(_lPt[i] > 40) ) continue;  
      if (_lPt[i] > minPt1){
	minPt1= _lPt[i];
	firstMu = i;	  
      }
    }

    //if( f1stlep != -1)   cout <<"1stPt =  "<<_lPt[f1stlep]<<endl;
	 
    Float_t   minPt2 = -1000;
    unsigned  secondMu = -1;
    for(unsigned i=0; i < _nL ; i++){
      if(_lFlavor[i] != 1 or _leptonMvatZqTTV[i] < -0.4) continue;
      if( firstMu == -1)  continue;
      if( _lPt[firstMu] > _lPt[i] and  _lPt[i] > 10) { 
	if (_lPt[i] > minPt2){
	  minPt2= _lPt[i];
	  secondMu = i;
	}
      }
    }

    //if( seclep != -1)   cout <<"2ndtPt =  "<<_lPt[seclep]<<endl;

    Float_t   minPt3 = -1000;
    unsigned  thirdMu = -1;
    for(unsigned i=0; i < _nL ; i++){
      if(_lFlavor[i] != 1 or _leptonMvatZqTTV[i] < -0.4) continue;
      if( secondMu == -1) continue;  
      if( _lPt[secondMu] > _lPt[i] and _lPt[i] > 10) {
	if (_lPt[i] > minPt3){
	  minPt3= _lPt[i];
	  thirdMu = i;
	}
      }
    }

    Float_t   minPt4 = -1000;
    unsigned  fourthMu = -1;
    for(unsigned i=0; i < _nL ; i++){
      if(_lFlavor[i] != 1 or _leptonMvatZqTTV[i] < -0.4) continue;
      if( thirdMu == -1) continue;
      if( _lPt[thirdMu] > _lPt[i] and _lPt[i] > 10) {
        if (_lPt[i] > minPt4){
          minPt4= _lPt[i];
          fourthMu = i;
        }
      }
    }

    Float_t   minPt5 = -1000;
    unsigned  firstEle = -1;
    for(unsigned i=0; i < _nL ; i++){
      if(_lFlavor[i] != 0 or _leptonMvatZqTTV[i] < -0.4) continue;
      if (!(_lPt[i] > 40) ) continue;
      if (_lPt[i] > minPt5){
	minPt5= _lPt[i];
        firstEle = i;
      }
    }

    Float_t   minPt6 = -1000;
    unsigned  secondEle = -1;
    for(unsigned i=0; i < _nL ; i++){
      if(_lFlavor[i] != 0 or _leptonMvatZqTTV[i] < -0.4) continue;
      if( firstEle == -1)  continue;
      if( _lPt[firstEle] > _lPt[i] and  _lPt[i] > 10) {
        if (_lPt[i] > minPt6){
          minPt6= _lPt[i];
          secondEle = i;
        }
      }
    }

    Float_t   minPt7 = -1000;
    unsigned  thirdEle = -1;
    for(unsigned i=0; i < _nL ; i++){
      if(_lFlavor[i] != 0 or _leptonMvatZqTTV[i] < -0.4) continue;
      if( secondEle == -1)  continue;
      if( _lPt[secondEle] > _lPt[i] and  _lPt[i] > 10) {
        if (_lPt[i] > minPt7){
          minPt7= _lPt[i];
          thirdEle = i;
        }
      }
    }


    Float_t   minPt8 = -1000;
    unsigned  fourthEle = -1;
    for(unsigned i=0; i < _nL ; i++){
      if(_lFlavor[i] != 0 or _leptonMvatZqTTV[i] < -0.4) continue;
      if( thirdEle == -1)  continue;
      if( _lPt[thirdEle] > _lPt[i] and  _lPt[i] > 10) {
        if (_lPt[i] > minPt8){
          minPt8= _lPt[i];
          fourthEle = i;
        }
      }
    }



    //if( thirdlep != -1)   cout <<"3rdPt =  "<<_lPt[thirdlep]<<endl;
    //cout<<"==================== end of the event ================"<<endl;


      TLorentzVector Mu1;
      TLorentzVector Mu2;
      TLorentzVector Mu3;
      TLorentzVector Mu4;

      float Mu_1stPt       = (firstMu != -1)  ?  _lPt[firstMu]         : -999 ;
      float Mu_1stEta      = (firstMu != -1)  ?  _lEta[firstMu]        : -999 ;
      float Mu_1stPhi      = (firstMu != -1)  ?  _lPhi[firstMu]        : -999 ;
      float Mu_1stE        = (firstMu != -1)  ?  _lE[firstMu]          : -999 ;
      float Mu_1stCharge   = (firstMu != -1)  ?  _lCharge[firstMu]     : -999 ;

      float Mu_2ndPt       = (secondMu != -1) ?  _lPt[secondMu]        : -999 ;
      float Mu_2ndEta      = (secondMu != -1) ?  _lEta[secondMu]       : -999 ;
      float Mu_2ndPhi      = (secondMu != -1) ?  _lPhi[secondMu]       : -999 ;
      float Mu_2ndE        = (secondMu != -1) ?  _lE[secondMu]         : -999 ;
      float Mu_2ndCharge   = (secondMu != -1) ?  _lCharge[secondMu]    : -999 ;

      float Mu_3rdPt       = (thirdMu != -1)  ? _lPt[thirdMu]          : -999 ;
      float Mu_3rdEta      = (thirdMu != -1)  ? _lEta[thirdMu]         : -999 ;
      float Mu_3rdPhi      = (thirdMu != -1)  ? _lPhi[thirdMu]         : -999 ;
      float Mu_3rdE        = (thirdMu != -1)  ? _lE[thirdMu]           : -999 ;
      float Mu_3rdCharge   = (thirdMu != -1)  ? _lCharge[thirdMu]      : -999 ;
 
      float Mu_4thPt       = (fourthMu != -1) ? _lPt[fourthMu]         : -999 ;
      float Mu_4thEta      = (fourthMu != -1) ? _lEta[fourthMu]        : -999 ;
      float Mu_4thPhi      = (fourthMu != -1) ? _lPhi[fourthMu]        : -999 ;
      float Mu_4thE        = (fourthMu != -1) ? _lE[fourthMu]          : -999 ;
      float Mu_4thCharge   = (fourthMu != -1) ? _lCharge[fourthMu]     : -999 ;

      Mu1.SetPtEtaPhiE(Mu_1stPt,Mu_1stEta,Mu_1stPhi,Mu_1stE);
      Mu2.SetPtEtaPhiE(Mu_2ndPt,Mu_2ndEta,Mu_2ndPhi,Mu_2ndE);
      Mu3.SetPtEtaPhiE(Mu_3rdPt,Mu_3rdEta,Mu_3rdPhi,Mu_3rdE);
      Mu4.SetPtEtaPhiE(Mu_4thPt,Mu_4thEta,Mu_4thPhi,Mu_4thE);

      float DiMuMass1  = (Mu1 + Mu2).M();
      float DiMuMass2  = (Mu1 + Mu3).M();
      float DiMuMass3  = (Mu1 + Mu4).M();
      float DiMuMass4  = (Mu2 + Mu3).M();
      float DiMuMass5  = (Mu2 + Mu4).M();
      float DiMuMass6  = (Mu3 + Mu4).M();

      TLorentzVector Ele1;
      TLorentzVector Ele2;
      TLorentzVector Ele3;
      TLorentzVector Ele4;

      float Ele_1stPt       = (firstEle != -1)  ?  _lPt[firstEle]         : -999 ;
      float Ele_1stEta      = (firstEle != -1)  ?  _lEta[firstEle]        : -999 ;
      float Ele_1stPhi      = (firstEle != -1)  ?  _lPhi[firstEle]        : -999 ;
      float Ele_1stE        = (firstEle != -1)  ?  _lE[firstEle]          : -999 ;
      float Ele_1stCharge   = (firstEle != -1)  ?  _lCharge[firstEle]     : -999 ;

      float Ele_2ndPt       = (secondEle != -1) ?  _lPt[secondEle]        : -999 ;
      float Ele_2ndEta      = (secondEle != -1) ?  _lEta[secondEle]       : -999 ;
      float Ele_2ndPhi      = (secondEle != -1) ?  _lPhi[secondEle]       : -999 ;
      float Ele_2ndE        = (secondEle != -1) ?  _lE[secondEle]         : -999 ;
      float Ele_2ndCharge   = (secondEle != -1) ?  _lCharge[secondEle]    : -999 ;

      float Ele_3rdPt       = (thirdEle != -1)  ?  _lPt[thirdEle]         : -999 ;
      float Ele_3rdEta      = (thirdEle != -1)  ?  _lEta[thirdEle]        : -999 ;
      float Ele_3rdPhi      = (thirdEle != -1)  ?  _lPhi[thirdEle]        : -999 ;
      float Ele_3rdE        = (thirdEle != -1)  ?  _lE[thirdEle]          : -999 ;
      float Ele_3rdCharge   = (thirdEle != -1)  ?  _lCharge[thirdEle]     : -999 ;

      float Ele_4thPt       = (fourthEle != -1) ?  _lPt[fourthEle]        : -999 ;
      float Ele_4thEta      = (fourthEle != -1) ?  _lEta[fourthEle]       : -999 ;
      float Ele_4thPhi      = (fourthEle != -1) ?  _lPhi[fourthEle]       : -999 ;
      float Ele_4thE        = (fourthEle != -1) ?  _lE[fourthEle]         : -999 ;
      float Ele_4thCharge   = (fourthEle != -1) ?  _lCharge[fourthEle]    : -999 ;

      Ele1.SetPtEtaPhiE(Ele_1stPt,Ele_1stEta,Ele_1stPhi,Ele_1stE);
      Ele2.SetPtEtaPhiE(Ele_2ndPt,Ele_2ndEta,Ele_2ndPhi,Ele_2ndE);
      Ele3.SetPtEtaPhiE(Ele_3rdPt,Ele_3rdEta,Ele_3rdPhi,Ele_3rdE);
      Ele4.SetPtEtaPhiE(Ele_4thPt,Ele_4thEta,Ele_4thPhi,Ele_4thE);

      float DiEleMass1  = (Ele1 + Ele2).M();
      float DiEleMass2  = (Ele1 + Ele3).M();
      float DiEleMass3  = (Ele1 + Ele4).M();
      float DiEleMass4  = (Ele2 + Ele3).M();
      float DiEleMass5  = (Ele2 + Ele4).M();
      float DiEleMass6  = (Ele3 + Ele4).M();



      if(firstMu != -1 && secondMu != -1 && thirdMu != -1 && fourthMu != -1){

      if(   ((Mu_1stCharge + Mu_2ndCharge) == 0  && abs(DiMuMass1 - Zmass) < 20 ) 
	 or ((Mu_1stCharge + Mu_3rdCharge) == 0  && abs(DiMuMass2 - Zmass) < 20 )	 
	 or ((Mu_1stCharge + Mu_4thCharge) == 0  && abs(DiMuMass3 - Zmass) < 20 )   
	 or ((Mu_2ndCharge + Mu_3rdCharge) == 0  && abs(DiMuMass4 - Zmass) < 20 )
	 or ((Mu_2ndCharge + Mu_4thCharge) == 0  && abs(DiMuMass5 - Zmass) < 20 )
	    or ((Mu_3rdCharge + Mu_4thCharge) == 0  && abs(DiMuMass6 - Zmass) < 20 ) )


      {

	mu_1stPt     = _lPt[firstMu];
        mu_1stEta    = _lEta[firstMu] ;
        mu_1stPhi    = _lPhi[firstMu];
        mu_1stE      = _lE[firstMu];
	mu_1stCharge = _lCharge[firstMu] ;


	mu_2ndPt     = _lPt[secondMu];
        mu_2ndEta    = _lEta[secondMu] ;
        mu_2ndPhi    = _lPhi[secondMu];
        mu_2ndE      = _lE[secondMu] ;
        mu_2ndCharge = _lCharge[secondMu] ;

	mu_3rdPt     = _lPt[thirdMu];
	mu_3rdEta    = _lEta[thirdMu] ;
	mu_3rdPhi    = _lPhi[thirdMu];
	mu_3rdE      = _lE[thirdMu] ;
        mu_3rdCharge = _lCharge[thirdMu] ;

        mu_4thPt     = _lPt[fourthMu] ;
        mu_4thEta    = _lEta[fourthMu] ;
        mu_4thPhi    = _lPhi[fourthMu];
        mu_4thE      = _lE[fourthMu] ;
        mu_4thCharge = _lCharge[fourthMu] ;


	DimuMass1 = DiMuMass1;
	DimuMass2 = DiMuMass2 ;
	DimuMass3 = DiMuMass3 ;
	DimuMass4 = DiMuMass4;
        DimuMass5 = DiMuMass5 ;
	DimuMass6 = DiMuMass6 ;

	newtree1->Fill();
      }
    }


      if(firstMu != -1 && secondMu != -1 && thirdMu != -1 && firstEle != -1){

	if(   ((Mu_1stCharge + Mu_2ndCharge) == 0  && abs(DiMuMass1 - Zmass) < 20 )
	      or ((Mu_1stCharge + Mu_3rdCharge) == 0  && abs(DiMuMass2 - Zmass) < 20 )
	      or ((Mu_2ndCharge + Mu_3rdCharge) == 0  && abs(DiMuMass4 - Zmass) < 20 ))
	  {
	    ////////////////////////


	    mu_1stPt     = _lPt[firstMu];
	    mu_1stEta    = _lEta[firstMu] ;
	    mu_1stPhi    = _lPhi[firstMu];
	    mu_1stE      = _lE[firstMu];
	    mu_1stCharge = _lCharge[firstMu] ;


	    mu_2ndPt     = _lPt[secondMu];
	    mu_2ndEta    = _lEta[secondMu] ;
	    mu_2ndPhi    = _lPhi[secondMu];
	    mu_2ndE      = _lE[secondMu] ;
	    mu_2ndCharge = _lCharge[secondMu] ;

	    mu_3rdPt     = _lPt[thirdMu];
	    mu_3rdEta    = _lEta[thirdMu] ;
	    mu_3rdPhi    = _lPhi[thirdMu];
	    mu_3rdE      = _lE[thirdMu] ;
	    mu_3rdCharge = _lCharge[thirdMu] ;

            ele_1stPt     = _lPt[firstEle];
            ele_1stEta    = _lEta[firstEle] ;
            ele_1stPhi    = _lPhi[firstEle];
            ele_1stE      = _lE[firstEle];
            ele_1stCharge = _lCharge[firstEle] ;


	    DimuMass1 = DiMuMass1;
	    DimuMass2 = DiMuMass2 ;
	    DimuMass4 = DiMuMass4;
  
	    newtree2->Fill();

	  }
      }



      if(firstMu != -1 && secondMu != -1 && firstEle != -1 && secondEle != -1){

        if(   ((Mu_1stCharge + Mu_2ndCharge) == 0  && abs(DiMuMass1 - Zmass) < 20 )
              or ((Ele_1stCharge + Ele_2ndCharge) == 0  && abs(DiEleMass2 - Zmass) < 20 ))
          {
	    //////////////////////////////

	    mu_1stPt     = _lPt[firstMu];
	    mu_1stEta    = _lEta[firstMu] ;
	    mu_1stPhi    = _lPhi[firstMu];
	    mu_1stE      = _lE[firstMu];
	    mu_1stCharge = _lCharge[firstMu] ;


	    mu_2ndPt     = _lPt[secondMu];
	    mu_2ndEta    = _lEta[secondMu] ;
	    mu_2ndPhi    = _lPhi[secondMu];
	    mu_2ndE      = _lE[secondMu] ;
	    mu_2ndCharge = _lCharge[secondMu] ;

            ele_1stPt     = _lPt[firstEle];
            ele_1stEta    = _lEta[firstEle] ;
            ele_1stPhi    = _lPhi[firstEle];
            ele_1stE      = _lE[firstEle];
            ele_1stCharge = _lCharge[firstEle] ;


            ele_2ndPt     = _lPt[secondEle];
            ele_2ndEta    = _lEta[secondEle] ;
            ele_2ndPhi    = _lPhi[secondEle];
            ele_2ndE      = _lE[secondEle] ;
            ele_2ndCharge = _lCharge[secondEle] ;

	    DimuMass1 = DiMuMass1;
            DieleMass1 = DiEleMass1;

	    newtree3->Fill();

          }
      }


      if(firstMu != -1 && firstEle != -1 && secondEle != -1 && thirdEle != -1){

	if(   ((Ele_1stCharge + Ele_2ndCharge) == 0  && abs(DiEleMass1 - Zmass) < 20 )
              or ((Ele_1stCharge + Ele_3rdCharge) == 0  && abs(DiEleMass2 - Zmass) < 20 )
              or ((Ele_2ndCharge + Ele_3rdCharge) == 0  && abs(DiEleMass4 - Zmass) < 20 ))
          {


	    mu_1stPt     = _lPt[firstMu];
	    mu_1stEta    = _lEta[firstMu] ;
	    mu_1stPhi    = _lPhi[firstMu];
	    mu_1stE      = _lE[firstMu];
	    mu_1stCharge = _lCharge[firstMu] ;


            ele_1stPt     = _lPt[firstEle];
            ele_1stEta    = _lEta[firstEle] ;
            ele_1stPhi    = _lPhi[firstEle];
            ele_1stE      = _lE[firstEle];
            ele_1stCharge = _lCharge[firstEle] ;


            ele_2ndPt     = _lPt[secondEle];
            ele_2ndEta    = _lEta[secondEle] ;
            ele_2ndPhi    = _lPhi[secondEle];
            ele_2ndE      = _lE[secondEle] ;
            ele_2ndCharge = _lCharge[secondEle] ;

            ele_3rdPt     = _lPt[thirdEle];
            ele_3rdEta    = _lEta[thirdEle] ;
            ele_3rdPhi    = _lPhi[thirdEle];
            ele_3rdE      = _lE[thirdEle] ;
            ele_3rdCharge = _lCharge[thirdEle] ;


            DieleMass1 = DiEleMass1;
            DieleMass2 = DiEleMass2 ;
            DieleMass4 = DiEleMass4;

	    newtree4->Fill();

	    /////////////////////////
          }
      }

      if(firstEle != -1 && secondEle != -1 && thirdEle != -1 && fourthEle != -1){

	if(   ((Ele_1stCharge + Ele_2ndCharge) == 0  && abs(DiEleMass1 - Zmass) < 20 )
	      or ((Ele_1stCharge + Ele_3rdCharge) == 0  && abs(DiEleMass2 - Zmass) < 20 )
	      or ((Ele_1stCharge + Ele_4thCharge) == 0  && abs(DiEleMass3 - Zmass) < 20 )
	      or ((Ele_2ndCharge + Ele_3rdCharge) == 0  && abs(DiEleMass4 - Zmass) < 20 )
	      or ((Ele_2ndCharge + Ele_4thCharge) == 0  && abs(DiEleMass5 - Zmass) < 20 )
	      or ((Ele_3rdCharge + Ele_4thCharge) == 0  && abs(DiEleMass6 - Zmass) < 20 ) )

	  {

	    ele_1stPt     = _lPt[firstEle];
	    ele_1stEta    = _lEta[firstEle] ;
	    ele_1stPhi    = _lPhi[firstEle];
	    ele_1stE      = _lE[firstEle];
	    ele_1stCharge = _lCharge[firstEle] ;


	    ele_2ndPt     = _lPt[secondEle];
	    ele_2ndEta    = _lEta[secondEle] ;
	    ele_2ndPhi    = _lPhi[secondEle];
	    ele_2ndE      = _lE[secondEle] ;
	    ele_2ndCharge = _lCharge[secondEle] ;

	    ele_3rdPt     = _lPt[thirdEle];
	    ele_3rdEta    = _lEta[thirdEle] ;
	    ele_3rdPhi    = _lPhi[thirdEle];
	    ele_3rdE      = _lE[thirdEle] ;
	    ele_3rdCharge = _lCharge[thirdEle] ;

	    ele_4thPt     = _lPt[fourthEle] ;
	    ele_4thEta    = _lEta[fourthEle] ;
	    ele_4thPhi    = _lPhi[fourthEle];
	    ele_4thE      = _lE[fourthEle] ;
	    ele_4thCharge = _lCharge[fourthEle] ;


	    DieleMass1 = DiEleMass1;
	    DieleMass2 = DiEleMass2 ;
	    DieleMass3 = DiEleMass3 ;
	    DieleMass4 = DiEleMass1;
	    DieleMass5 = DiEleMass2 ;
	    DieleMass6 = DiEleMass3 ;

	    newtree5->Fill();

	    //////////////////////////
	  }
      }


  }


 newfile->Write();  
 newtree1->Print();
 newtree2->Print();
 newtree3->Print();
 newtree4->Print();
 newtree5->Print();

 newtree1->AutoSave();
 newtree2->AutoSave();
 newtree3->AutoSave();
 newtree4->AutoSave();
 newtree5->AutoSave();

  delete newfile;

  return 0;
}
