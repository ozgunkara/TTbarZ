//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jul  8 14:37:31 2019 by ROOT version 6.16/00
// from TTree blackJackAndHookersTree/blackJackAndHookersTree
// found on file: /afs/cern.ch/user/o/okara/public/forsoner/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8_Summer16.root
//////////////////////////////////////////////////////////

#ifndef okara_h
#define okara_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class okara {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
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

   // List of branches
   TBranch        *b__runNb;   //!
   TBranch        *b__lumiBlock;   //!
   TBranch        *b__eventNb;   //!
   TBranch        *b__nVertex;   //!
   TBranch        *b__met;   //!
   TBranch        *b__metJECDown;   //!
   TBranch        *b__metJECUp;   //!
   TBranch        *b__metJetResDown;   //!
   TBranch        *b__metJetResUp;   //!
   TBranch        *b__metUnclDown;   //!
   TBranch        *b__metUnclUp;   //!
   TBranch        *b__metPhi;   //!
   TBranch        *b__metPhiJECDown;   //!
   TBranch        *b__metPhiJECUp;   //!
   TBranch        *b__metPhiJetResDown;   //!
   TBranch        *b__metPhiJetResUp;   //!
   TBranch        *b__metPhiUnclDown;   //!
   TBranch        *b__metPhiUnclUp;   //!
   TBranch        *b__2016_FR;   //!
   TBranch        *b__2017_FR;   //!
   TBranch        *b__HLT_Mu3_PFJet40;   //!
   TBranch        *b__HLT_Mu3_PFJet40_prescale;   //!
   TBranch        *b__HLT_Mu8;   //!
   TBranch        *b__HLT_Mu8_prescale;   //!
   TBranch        *b__HLT_Mu17;   //!
   TBranch        *b__HLT_Mu17_prescale;   //!
   TBranch        *b__HLT_Mu27;   //!
   TBranch        *b__HLT_Mu27_prescale;   //!
   TBranch        *b__HLT_Ele8_CaloIdM_TrackIdM_PFJet30;   //!
   TBranch        *b__HLT_Ele8_CaloIdM_TrackIdM_PFJet30_prescale;   //!
   TBranch        *b__HLT_Ele12_CaloIdM_TrackIdM_PFJet30;   //!
   TBranch        *b__HLT_Ele12_CaloIdM_TrackIdM_PFJet30_prescale;   //!
   TBranch        *b__HLT_Ele17_CaloIdM_TrackIdM_PFJet30;   //!
   TBranch        *b__HLT_Ele17_CaloIdM_TrackIdM_PFJet30_prescale;   //!
   TBranch        *b__HLT_Ele23_CaloIdM_TrackIdM_PFJet30;   //!
   TBranch        *b__HLT_Ele23_CaloIdM_TrackIdM_PFJet30_prescale;   //!
   TBranch        *b__passTrigger_e;   //!
   TBranch        *b__passTrigger_m;   //!
   TBranch        *b__passTrigger_ee;   //!
   TBranch        *b__passTrigger_em;   //!
   TBranch        *b__passTrigger_mm;   //!
   TBranch        *b__passTrigger_eee;   //!
   TBranch        *b__passTrigger_eem;   //!
   TBranch        *b__passTrigger_emm;   //!
   TBranch        *b__passTrigger_mmm;   //!
   TBranch        *b__passMETFilters;   //!
   TBranch        *b__nL;   //!
   TBranch        *b__nMu;   //!
   TBranch        *b__nEle;   //!
   TBranch        *b__nLight;   //!
   TBranch        *b__lPt;   //!
   TBranch        *b__lEta;   //!
   TBranch        *b__lEtaSC;   //!
   TBranch        *b__lPhi;   //!
   TBranch        *b__lE;   //!
   TBranch        *b__lFlavor;   //!
   TBranch        *b__lCharge;   //!
   TBranch        *b__dxy;   //!
   TBranch        *b__dz;   //!
   TBranch        *b__3dIP;   //!
   TBranch        *b__3dIPSig;   //!
   TBranch        *b__lElectronMva;   //!
   TBranch        *b__lElectronMvaHZZ;   //!
   TBranch        *b__lElectronMvaFall17Iso;   //!
   TBranch        *b__lElectronMvaFall17NoIso;   //!
   TBranch        *b__lElectronPassEmu;   //!
   TBranch        *b__lElectronPassConvVeto;   //!
   TBranch        *b__lElectronChargeConst;   //!
   TBranch        *b__lElectronMissingHits;   //!
   TBranch        *b__leptonMvaSUSY;   //!
   TBranch        *b__leptonMvaTTH;   //!
   TBranch        *b__leptonMvatZqTTV;   //!
   TBranch        *b__lPOGTight;   //!
   TBranch        *b__lPOGLoose;   //!
   TBranch        *b__lPOGMedium;   //!
   TBranch        *b__lPOGLooseWOIso;   //!
   TBranch        *b__lPOGMediumWOIso;   //!
   TBranch        *b__lPOGTightWOIso;   //!
   TBranch        *b__relIso;   //!
   TBranch        *b__relIso0p4Mu;   //!
   TBranch        *b__relIso0p4;   //!
   TBranch        *b__relIso0p6;   //!
   TBranch        *b__relIso0p8;   //!
   TBranch        *b__relIso1p0;   //!
   TBranch        *b__miniIso;   //!
   TBranch        *b__miniIsoCharged;   //!
   TBranch        *b__ptRel;   //!
   TBranch        *b__ptRatio;   //!
   TBranch        *b__closestJetCsvV2;   //!
   TBranch        *b__closestJetDeepCsv_b;   //!
   TBranch        *b__closestJetDeepCsv_bb;   //!
   TBranch        *b__selectedTrackMult;   //!
   TBranch        *b__lMuonSegComp;   //!
   TBranch        *b__lMuonTrackPt;   //!
   TBranch        *b__lMuonTrackPtErr;   //!
   TBranch        *b__nJets;   //!
   TBranch        *b__jetPt;   //!
   TBranch        *b__jetSmearedPt;   //!
   TBranch        *b__jetSmearedPt_JECDown;   //!
   TBranch        *b__jetSmearedPt_JECUp;   //!
   TBranch        *b__jetSmearedPt_JERDown;   //!
   TBranch        *b__jetSmearedPt_JERUp;   //!
   TBranch        *b__jetPt_JECUp;   //!
   TBranch        *b__jetPt_JECDown;   //!
   TBranch        *b__jetPt_JERUp;   //!
   TBranch        *b__jetPt_JERDown;   //!
   TBranch        *b__jetEta;   //!
   TBranch        *b__jetPhi;   //!
   TBranch        *b__jetE;   //!
   TBranch        *b__jetCsvV2;   //!
   TBranch        *b__jetDeepCsv_udsg;   //!
   TBranch        *b__jetDeepCsv_b;   //!
   TBranch        *b__jetDeepCsv_c;   //!
   TBranch        *b__jetDeepCsv_bb;   //!
   TBranch        *b__jetHadronFlavor;   //!
   TBranch        *b__jetIsTight;   //!
   TBranch        *b__jetIsTightLepVeto;   //!
   TBranch        *b__nLheWeights;   //!
   TBranch        *b__lheWeight;   //!
   TBranch        *b__nPsWeights;   //!
   TBranch        *b__psWeight;   //!
   TBranch        *b__lIsPrompt;   //!
   TBranch        *b__lMatchPdgId;   //!
   TBranch        *b__lMomPdgId;   //!
   TBranch        *b__weight;   //!
   TBranch        *b__nTrueInt;   //!
   TBranch        *b__gen_met;   //!
   TBranch        *b__gen_metPhi;   //!
   TBranch        *b__gen_nL;   //!
   TBranch        *b__gen_lPt;   //!
   TBranch        *b__gen_lEta;   //!
   TBranch        *b__gen_lPhi;   //!
   TBranch        *b__gen_lE;   //!
   TBranch        *b__gen_lFlavor;   //!
   TBranch        *b__gen_lCharge;   //!
   TBranch        *b__gen_lMomPdg;   //!
   TBranch        *b__gen_lIsPrompt;   //!
   TBranch        *b__gen_partonPt;   //!
   TBranch        *b__lProvenance;   //!
   TBranch        *b__lProvenanceCompressed;   //!

   okara(TTree *tree=0);
   virtual ~okara();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef okara_cxx
okara::okara(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/afs/cern.ch/user/o/okara/public/forsoner/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8_Summer16.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/afs/cern.ch/user/o/okara/public/forsoner/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8_Summer16.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/afs/cern.ch/user/o/okara/public/forsoner/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8_Summer16.root:/blackJackAndHookers");
      dir->GetObject("blackJackAndHookersTree",tree);

   }
   Init(tree);
}

okara::~okara()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t okara::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t okara::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void okara::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("_runNb", &_runNb, &b__runNb);
   fChain->SetBranchAddress("_lumiBlock", &_lumiBlock, &b__lumiBlock);
   fChain->SetBranchAddress("_eventNb", &_eventNb, &b__eventNb);
   fChain->SetBranchAddress("_nVertex", &_nVertex, &b__nVertex);
   fChain->SetBranchAddress("_met", &_met, &b__met);
   fChain->SetBranchAddress("_metJECDown", &_metJECDown, &b__metJECDown);
   fChain->SetBranchAddress("_metJECUp", &_metJECUp, &b__metJECUp);
   fChain->SetBranchAddress("_metJetResDown", &_metJetResDown, &b__metJetResDown);
   fChain->SetBranchAddress("_metJetResUp", &_metJetResUp, &b__metJetResUp);
   fChain->SetBranchAddress("_metUnclDown", &_metUnclDown, &b__metUnclDown);
   fChain->SetBranchAddress("_metUnclUp", &_metUnclUp, &b__metUnclUp);
   fChain->SetBranchAddress("_metPhi", &_metPhi, &b__metPhi);
   fChain->SetBranchAddress("_metPhiJECDown", &_metPhiJECDown, &b__metPhiJECDown);
   fChain->SetBranchAddress("_metPhiJECUp", &_metPhiJECUp, &b__metPhiJECUp);
   fChain->SetBranchAddress("_metPhiJetResDown", &_metPhiJetResDown, &b__metPhiJetResDown);
   fChain->SetBranchAddress("_metPhiJetResUp", &_metPhiJetResUp, &b__metPhiJetResUp);
   fChain->SetBranchAddress("_metPhiUnclDown", &_metPhiUnclDown, &b__metPhiUnclDown);
   fChain->SetBranchAddress("_metPhiUnclUp", &_metPhiUnclUp, &b__metPhiUnclUp);
   fChain->SetBranchAddress("_2016_FR", &_2016_FR, &b__2016_FR);
   fChain->SetBranchAddress("_2017_FR", &_2017_FR, &b__2017_FR);
   fChain->SetBranchAddress("_HLT_Mu3_PFJet40", &_HLT_Mu3_PFJet40, &b__HLT_Mu3_PFJet40);
   fChain->SetBranchAddress("_HLT_Mu3_PFJet40_prescale", &_HLT_Mu3_PFJet40_prescale, &b__HLT_Mu3_PFJet40_prescale);
   fChain->SetBranchAddress("_HLT_Mu8", &_HLT_Mu8, &b__HLT_Mu8);
   fChain->SetBranchAddress("_HLT_Mu8_prescale", &_HLT_Mu8_prescale, &b__HLT_Mu8_prescale);
   fChain->SetBranchAddress("_HLT_Mu17", &_HLT_Mu17, &b__HLT_Mu17);
   fChain->SetBranchAddress("_HLT_Mu17_prescale", &_HLT_Mu17_prescale, &b__HLT_Mu17_prescale);
   fChain->SetBranchAddress("_HLT_Mu27", &_HLT_Mu27, &b__HLT_Mu27);
   fChain->SetBranchAddress("_HLT_Mu27_prescale", &_HLT_Mu27_prescale, &b__HLT_Mu27_prescale);
   fChain->SetBranchAddress("_HLT_Ele8_CaloIdM_TrackIdM_PFJet30", &_HLT_Ele8_CaloIdM_TrackIdM_PFJet30, &b__HLT_Ele8_CaloIdM_TrackIdM_PFJet30);
   fChain->SetBranchAddress("_HLT_Ele8_CaloIdM_TrackIdM_PFJet30_prescale", &_HLT_Ele8_CaloIdM_TrackIdM_PFJet30_prescale, &b__HLT_Ele8_CaloIdM_TrackIdM_PFJet30_prescale);
   fChain->SetBranchAddress("_HLT_Ele12_CaloIdM_TrackIdM_PFJet30", &_HLT_Ele12_CaloIdM_TrackIdM_PFJet30, &b__HLT_Ele12_CaloIdM_TrackIdM_PFJet30);
   fChain->SetBranchAddress("_HLT_Ele12_CaloIdM_TrackIdM_PFJet30_prescale", &_HLT_Ele12_CaloIdM_TrackIdM_PFJet30_prescale, &b__HLT_Ele12_CaloIdM_TrackIdM_PFJet30_prescale);
   fChain->SetBranchAddress("_HLT_Ele17_CaloIdM_TrackIdM_PFJet30", &_HLT_Ele17_CaloIdM_TrackIdM_PFJet30, &b__HLT_Ele17_CaloIdM_TrackIdM_PFJet30);
   fChain->SetBranchAddress("_HLT_Ele17_CaloIdM_TrackIdM_PFJet30_prescale", &_HLT_Ele17_CaloIdM_TrackIdM_PFJet30_prescale, &b__HLT_Ele17_CaloIdM_TrackIdM_PFJet30_prescale);
   fChain->SetBranchAddress("_HLT_Ele23_CaloIdM_TrackIdM_PFJet30", &_HLT_Ele23_CaloIdM_TrackIdM_PFJet30, &b__HLT_Ele23_CaloIdM_TrackIdM_PFJet30);
   fChain->SetBranchAddress("_HLT_Ele23_CaloIdM_TrackIdM_PFJet30_prescale", &_HLT_Ele23_CaloIdM_TrackIdM_PFJet30_prescale, &b__HLT_Ele23_CaloIdM_TrackIdM_PFJet30_prescale);
   fChain->SetBranchAddress("_passTrigger_e", &_passTrigger_e, &b__passTrigger_e);
   fChain->SetBranchAddress("_passTrigger_m", &_passTrigger_m, &b__passTrigger_m);
   fChain->SetBranchAddress("_passTrigger_ee", &_passTrigger_ee, &b__passTrigger_ee);
   fChain->SetBranchAddress("_passTrigger_em", &_passTrigger_em, &b__passTrigger_em);
   fChain->SetBranchAddress("_passTrigger_mm", &_passTrigger_mm, &b__passTrigger_mm);
   fChain->SetBranchAddress("_passTrigger_eee", &_passTrigger_eee, &b__passTrigger_eee);
   fChain->SetBranchAddress("_passTrigger_eem", &_passTrigger_eem, &b__passTrigger_eem);
   fChain->SetBranchAddress("_passTrigger_emm", &_passTrigger_emm, &b__passTrigger_emm);
   fChain->SetBranchAddress("_passTrigger_mmm", &_passTrigger_mmm, &b__passTrigger_mmm);
   fChain->SetBranchAddress("_passMETFilters", &_passMETFilters, &b__passMETFilters);
   fChain->SetBranchAddress("_nL", &_nL, &b__nL);
   fChain->SetBranchAddress("_nMu", &_nMu, &b__nMu);
   fChain->SetBranchAddress("_nEle", &_nEle, &b__nEle);
   fChain->SetBranchAddress("_nLight", &_nLight, &b__nLight);
   fChain->SetBranchAddress("_lPt", _lPt, &b__lPt);
   fChain->SetBranchAddress("_lEta", _lEta, &b__lEta);
   fChain->SetBranchAddress("_lEtaSC", _lEtaSC, &b__lEtaSC);
   fChain->SetBranchAddress("_lPhi", _lPhi, &b__lPhi);
   fChain->SetBranchAddress("_lE", _lE, &b__lE);
   fChain->SetBranchAddress("_lFlavor", _lFlavor, &b__lFlavor);
   fChain->SetBranchAddress("_lCharge", _lCharge, &b__lCharge);
   fChain->SetBranchAddress("_dxy", _dxy, &b__dxy);
   fChain->SetBranchAddress("_dz", _dz, &b__dz);
   fChain->SetBranchAddress("_3dIP", _3dIP, &b__3dIP);
   fChain->SetBranchAddress("_3dIPSig", _3dIPSig, &b__3dIPSig);
   fChain->SetBranchAddress("_lElectronMva", _lElectronMva, &b__lElectronMva);
   fChain->SetBranchAddress("_lElectronMvaHZZ", _lElectronMvaHZZ, &b__lElectronMvaHZZ);
   fChain->SetBranchAddress("_lElectronMvaFall17Iso", _lElectronMvaFall17Iso, &b__lElectronMvaFall17Iso);
   fChain->SetBranchAddress("_lElectronMvaFall17NoIso", _lElectronMvaFall17NoIso, &b__lElectronMvaFall17NoIso);
   fChain->SetBranchAddress("_lElectronPassEmu", _lElectronPassEmu, &b__lElectronPassEmu);
   fChain->SetBranchAddress("_lElectronPassConvVeto", _lElectronPassConvVeto, &b__lElectronPassConvVeto);
   fChain->SetBranchAddress("_lElectronChargeConst", _lElectronChargeConst, &b__lElectronChargeConst);
   fChain->SetBranchAddress("_lElectronMissingHits", _lElectronMissingHits, &b__lElectronMissingHits);
   fChain->SetBranchAddress("_leptonMvaSUSY", _leptonMvaSUSY, &b__leptonMvaSUSY);
   fChain->SetBranchAddress("_leptonMvaTTH", _leptonMvaTTH, &b__leptonMvaTTH);
   fChain->SetBranchAddress("_leptonMvatZqTTV", _leptonMvatZqTTV, &b__leptonMvatZqTTV);
   fChain->SetBranchAddress("_lPOGTight", _lPOGTight, &b__lPOGTight);
   fChain->SetBranchAddress("_lPOGLoose", _lPOGLoose, &b__lPOGLoose);
   fChain->SetBranchAddress("_lPOGMedium", _lPOGMedium, &b__lPOGMedium);
   fChain->SetBranchAddress("_lPOGLooseWOIso", _lPOGLooseWOIso, &b__lPOGLooseWOIso);
   fChain->SetBranchAddress("_lPOGMediumWOIso", _lPOGMediumWOIso, &b__lPOGMediumWOIso);
   fChain->SetBranchAddress("_lPOGTightWOIso", _lPOGTightWOIso, &b__lPOGTightWOIso);
   fChain->SetBranchAddress("_relIso", _relIso, &b__relIso);
   fChain->SetBranchAddress("_relIso0p4Mu", _relIso0p4Mu, &b__relIso0p4Mu);
   fChain->SetBranchAddress("_relIso0p4", _relIso0p4, &b__relIso0p4);
   fChain->SetBranchAddress("_relIso0p6", _relIso0p6, &b__relIso0p6);
   fChain->SetBranchAddress("_relIso0p8", _relIso0p8, &b__relIso0p8);
   fChain->SetBranchAddress("_relIso1p0", _relIso1p0, &b__relIso1p0);
   fChain->SetBranchAddress("_miniIso", _miniIso, &b__miniIso);
   fChain->SetBranchAddress("_miniIsoCharged", _miniIsoCharged, &b__miniIsoCharged);
   fChain->SetBranchAddress("_ptRel", _ptRel, &b__ptRel);
   fChain->SetBranchAddress("_ptRatio", _ptRatio, &b__ptRatio);
   fChain->SetBranchAddress("_closestJetCsvV2", _closestJetCsvV2, &b__closestJetCsvV2);
   fChain->SetBranchAddress("_closestJetDeepCsv_b", _closestJetDeepCsv_b, &b__closestJetDeepCsv_b);
   fChain->SetBranchAddress("_closestJetDeepCsv_bb", _closestJetDeepCsv_bb, &b__closestJetDeepCsv_bb);
   fChain->SetBranchAddress("_selectedTrackMult", _selectedTrackMult, &b__selectedTrackMult);
   fChain->SetBranchAddress("_lMuonSegComp", _lMuonSegComp, &b__lMuonSegComp);
   fChain->SetBranchAddress("_lMuonTrackPt", _lMuonTrackPt, &b__lMuonTrackPt);
   fChain->SetBranchAddress("_lMuonTrackPtErr", _lMuonTrackPtErr, &b__lMuonTrackPtErr);
   fChain->SetBranchAddress("_nJets", &_nJets, &b__nJets);
   fChain->SetBranchAddress("_jetPt", _jetPt, &b__jetPt);
   fChain->SetBranchAddress("_jetSmearedPt", _jetSmearedPt, &b__jetSmearedPt);
   fChain->SetBranchAddress("_jetSmearedPt_JECDown", _jetSmearedPt_JECDown, &b__jetSmearedPt_JECDown);
   fChain->SetBranchAddress("_jetSmearedPt_JECUp", _jetSmearedPt_JECUp, &b__jetSmearedPt_JECUp);
   fChain->SetBranchAddress("_jetSmearedPt_JERDown", _jetSmearedPt_JERDown, &b__jetSmearedPt_JERDown);
   fChain->SetBranchAddress("_jetSmearedPt_JERUp", _jetSmearedPt_JERUp, &b__jetSmearedPt_JERUp);
   fChain->SetBranchAddress("_jetPt_JECUp", _jetPt_JECUp, &b__jetPt_JECUp);
   fChain->SetBranchAddress("_jetPt_JECDown", _jetPt_JECDown, &b__jetPt_JECDown);
   fChain->SetBranchAddress("_jetPt_JERUp", _jetPt_JERUp, &b__jetPt_JERUp);
   fChain->SetBranchAddress("_jetPt_JERDown", _jetPt_JERDown, &b__jetPt_JERDown);
   fChain->SetBranchAddress("_jetEta", _jetEta, &b__jetEta);
   fChain->SetBranchAddress("_jetPhi", _jetPhi, &b__jetPhi);
   fChain->SetBranchAddress("_jetE", _jetE, &b__jetE);
   fChain->SetBranchAddress("_jetCsvV2", _jetCsvV2, &b__jetCsvV2);
   fChain->SetBranchAddress("_jetDeepCsv_udsg", _jetDeepCsv_udsg, &b__jetDeepCsv_udsg);
   fChain->SetBranchAddress("_jetDeepCsv_b", _jetDeepCsv_b, &b__jetDeepCsv_b);
   fChain->SetBranchAddress("_jetDeepCsv_c", _jetDeepCsv_c, &b__jetDeepCsv_c);
   fChain->SetBranchAddress("_jetDeepCsv_bb", _jetDeepCsv_bb, &b__jetDeepCsv_bb);
   fChain->SetBranchAddress("_jetHadronFlavor", _jetHadronFlavor, &b__jetHadronFlavor);
   fChain->SetBranchAddress("_jetIsTight", _jetIsTight, &b__jetIsTight);
   fChain->SetBranchAddress("_jetIsTightLepVeto", _jetIsTightLepVeto, &b__jetIsTightLepVeto);
   fChain->SetBranchAddress("_nLheWeights", &_nLheWeights, &b__nLheWeights);
   fChain->SetBranchAddress("_lheWeight", _lheWeight, &b__lheWeight);
   fChain->SetBranchAddress("_nPsWeights", &_nPsWeights, &b__nPsWeights);
   fChain->SetBranchAddress("_psWeight", _psWeight, &b__psWeight);
   fChain->SetBranchAddress("_lIsPrompt", _lIsPrompt, &b__lIsPrompt);
   fChain->SetBranchAddress("_lMatchPdgId", _lMatchPdgId, &b__lMatchPdgId);
   fChain->SetBranchAddress("_lMomPdgId", _lMomPdgId, &b__lMomPdgId);
   fChain->SetBranchAddress("_weight", &_weight, &b__weight);
   fChain->SetBranchAddress("_nTrueInt", &_nTrueInt, &b__nTrueInt);
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
   fChain->SetBranchAddress("_gen_lIsPrompt", _gen_lIsPrompt, &b__gen_lIsPrompt);
   fChain->SetBranchAddress("_gen_partonPt", _gen_partonPt, &b__gen_partonPt);
   fChain->SetBranchAddress("_lProvenance", _lProvenance, &b__lProvenance);
   fChain->SetBranchAddress("_lProvenanceCompressed", _lProvenanceCompressed, &b__lProvenanceCompressed);
   Notify();
}

Bool_t okara::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void okara::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t okara::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef okara_cxx
