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
#include "../interface/treeReader.h"

using namespace std;
using namespace RooFit;
using namespace RooStats;

//=============================================================================================//

int main(){
  
  //Get old file, old tree and set top branch address
  TFile *oldfile = new TFile("/afs/cern.ch/work/o/okara/TTbarZ/TTbarZ_analysis/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8_Summer16.root");
  TTree *fChain = (TTree*)oldfile->Get("blackJackAndHookers/blackJackAndHookersTree");

//=============================================================================================//

  char fileName[256];
  cout<<"Please enter the name of the output root file you want to create (yyy.root) : "<<endl;
  cin.getline(fileName,256);
  TFile *newfile = new TFile(fileName,"recreate");

  Long64_t nentries = fChain->GetEntries();
  cout  <<nentries<<endl;

 //Histograms

 TH1D *lepton_1stPt = new TH1D("lepton_1stPt", "", 30, 0, 300);
 TH1D *lepton_2ndPt = new TH1D("lepton_2ndPt", "", 30, 0, 300);
 TH1D *lepton_3rdPt = new TH1D("lepton_3rdPt", "", 30, 0, 300);
 TH1D *lepton_4thPt = new TH1D("lepton_4thPt", "", 30, 0, 300);

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

  //fChain start here 

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

//======================= Start the running over input branches ==========================================//

for (int i=0;i<100000; i++) {

	if (i%10000==0)       cout<<i<<endl;
    	fChain->GetEntry(i);
    	//    if (passIsoMu24All==0 && passIsoMu27All == 0) continue;  // cut on the trigger!


  	std::vector<unsigned> leptonIndex;
	unsigned lCount = 0;
	
	for(unsigned l = 0; l < _nLight; l++){
		if(_lPt[l] > 10){
			lCount++;
			leptonIndex.push_back(l);		
		}
	}
  
	if (leptonIndex.size() == 4)
	{
		for(std::size_t l=0; l<leptonIndex.size(); ++l){
			
			lepton_1stPt->Fill(_lPt[leptonIndex[0]]);
                        lepton_2ndPt->Fill(_lPt[leptonIndex[1]]);
                        lepton_3rdPt->Fill(_lPt[leptonIndex[2]]);
                        lepton_4thPt->Fill(_lPt[leptonIndex[3]]);


		}
	}
}


//histogram styles //

//gStyle->SetOptStat(0);

lepton_1stPt->SetXTitle("p_{T} [GeV]");
lepton_1stPt->SetYTitle("Number Of Event");
lepton_2ndPt->SetXTitle("p_{T} [GeV]");
lepton_2ndPt->SetYTitle("Number Of Event");
lepton_3rdPt->SetXTitle("p_{T} [GeV]");
lepton_3rdPt->SetYTitle("Number Of Event");
lepton_4thPt->SetXTitle("p_{T} [GeV]");
lepton_4thPt->SetYTitle("Number Of Event");

newfile->Write();

return 0;
}
