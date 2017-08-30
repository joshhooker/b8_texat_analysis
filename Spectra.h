//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Aug 10 09:32:04 2017 by ROOT version 6.08/06
// from TTree mfmData/Experimental Event Data
// found on file: run1_0.root
//////////////////////////////////////////////////////////

#ifndef Spectra_h
#define Spectra_h

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TROOT.h>
#include <TStyle.h>

#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

// Header file for the classes stored in the TTree if any.

class Spectra {
public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

  // Declaration of leaf types
  Int_t           mmMul;
  Int_t           mmCobo[2000];   //[mmMul]
  Int_t           mmAsad[2000];   //[mmMul]
  Int_t           mmAget[2000];   //[mmMul]
  Int_t           mmChan[2000];   //[mmMul]
  Float_t         mmTime[2000];   //[mmMul]
  Float_t         mmEnergy[2000];   //[mmMul]

  // List of branches
  TBranch        *b_mmMul;   //!
  TBranch        *b_mmCobo;   //!
  TBranch        *b_mmAsad;   //!
  TBranch        *b_mmAget;   //!
  TBranch        *b_mmChan;   //!
  TBranch        *b_mmTime;   //!
  TBranch        *b_mmEnergy;   //!

  Spectra(TTree *tree=0);
  virtual ~Spectra();
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);

// Channel Map
private:
  void InitChannelMap();

  map<int,pair<int,int> > siForwardMap;
  int siForwardChannel[10][4];

  map<int,pair<int,int> > siLeftMap;
  int siLeftChannel[6][4];

  map<int,int> CsIForwardMap;
  int CsIForwardChannel[10];

  map<int,int> CsILeftMap;
  int CsILeftChannel[6];

  map<int, int> Aget_Map;

  map<int,pair<int,int> > MM_Map_Asad0_Aget0;
  map<int,pair<int,int> > MM_Map_Asad0_Aget1;
  map<int,pair<int,int> > MM_Map_Asad0_Aget2;
  map<int,pair<int,int> > MM_Map_Asad0_Aget3;

  map<int,pair<int,int> > MM_Map_Asad1_Aget0;
  map<int,pair<int,int> > MM_Map_Asad1_Aget1;
  map<int,pair<int,int> > MM_Map_Asad1_Aget2;
  map<int,pair<int,int> > MM_Map_Asad1_Aget3;

  map<int,int> MM_Map_Asad2_Aget0;
  map<int,int> MM_Map_Asad2_Aget1;
  map<int,pair<int,int> > MM_Map_Asad2_Aget2;
  map<int,pair<int,int> > MM_Map_Asad2_Aget3;

  map<int,int> MM_Map_Asad3_Aget0;
  map<int,int> MM_Map_Asad3_Aget1;
  map<int,pair<int,int> > MM_Map_Asad3_Aget2;
  map<int,pair<int,int> > MM_Map_Asad3_Aget3;

// Histograms
private:
  TH1F* hSiEForward[10][4];
  TH1F* hSiTForward[10][4];

  TH1F* hSiELeft[8][4];
  TH1F* hSiTLeft[8][4];

  TH1F* hCsIEForward[10];
  TH1F* hCsITForward[10];

  TH1F* hCsIELeft[8];
  TH1F* hCsITLeft[8];
};

#endif

#ifdef Spectra_cxx

inline Spectra::Spectra(TTree *tree) : fChain(0) {
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("run1_0.root");
    if (!f || !f->IsOpen()) {
      f = new TFile("run1_0.root");
    }
    f->GetObject("mfmData",tree);
  }
  Init(tree);
}

inline Spectra::~Spectra() {
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

inline Int_t Spectra::GetEntry(Long64_t entry) {
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

inline Long64_t Spectra::LoadTree(Long64_t entry) {
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

inline void Spectra::Init(TTree *tree) {
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

  fChain->SetBranchAddress("mmMul", &mmMul, &b_mmMul);
  fChain->SetBranchAddress("mmCobo", mmCobo, &b_mmCobo);
  fChain->SetBranchAddress("mmAsad", mmAsad, &b_mmAsad);
  fChain->SetBranchAddress("mmAget", mmAget, &b_mmAget);
  fChain->SetBranchAddress("mmChan", mmChan, &b_mmChan);
  fChain->SetBranchAddress("mmTime", mmTime, &b_mmTime);
  fChain->SetBranchAddress("mmEnergy", mmEnergy, &b_mmEnergy);
  Notify();
}

Bool_t Spectra::Notify() {
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

inline void Spectra::Show(Long64_t entry) {
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}

inline void Spectra::InitChannelMap() {
  //////////////////////////
  // Forward Si Detectors //
  //////////////////////////

  // Detector 1
  siForwardMap[5] = make_pair(1, 1);
  siForwardMap[12] = make_pair(1, 2);
  siForwardMap[18] = make_pair(1, 3);
  siForwardMap[25] = make_pair(1, 4);

  // Detector 2
  siForwardMap[31] = make_pair(2, 1);
  siForwardMap[37] = make_pair(2, 2);
  siForwardMap[43] = make_pair(2, 3);
  siForwardMap[50] = make_pair(2, 4);

  // Detector 3
  siForwardMap[54] = make_pair(3, 1);
  siForwardMap[59] = make_pair(3, 2);
  siForwardMap[63] = make_pair(3, 3);
  siForwardMap[67] = make_pair(3, 4);

  // Detector 4
  siForwardMap[23] = make_pair(4, 1);
  siForwardMap[3] = make_pair(4, 2);
  siForwardMap[9] = make_pair(4, 3);
  siForwardMap[16] = make_pair(4, 4);

  // Detector 5
  siForwardMap[29] = make_pair(5, 1);
  siForwardMap[35] = make_pair(5, 2);
  siForwardMap[41] = make_pair(5, 3);
  siForwardMap[48] = make_pair(5, 4);

  // Detector 6
  siForwardMap[52] = make_pair(6, 1);
  siForwardMap[57] = make_pair(6, 2);
  siForwardMap[61] = make_pair(6, 3);
  siForwardMap[65] = make_pair(6, 4);

  // Detector 7
  siForwardMap[1] = make_pair(7, 1);
  siForwardMap[7] = make_pair(7, 2);
  siForwardMap[14] = make_pair(7, 3);
  siForwardMap[20] = make_pair(7, 4);

  // Detector 8
  siForwardMap[46] = make_pair(8, 1);
  siForwardMap[27] = make_pair(8, 2);
  siForwardMap[33] = make_pair(8, 3);
  siForwardMap[39] = make_pair(8, 4);

  // Detector 9
  siForwardMap[4] = make_pair(9, 1);
  siForwardMap[10] = make_pair(9, 2);
  siForwardMap[17] = make_pair(9, 3);
  siForwardMap[24] = make_pair(9, 4);

  // Detector 10
  siForwardMap[30] = make_pair(10, 1);
  siForwardMap[36] = make_pair(10, 2);
  siForwardMap[42] = make_pair(10, 3);
  siForwardMap[49] = make_pair(10, 4);

  siForwardChannel[0][0] = 5;
  siForwardChannel[0][1] = 12;
  siForwardChannel[0][2] = 18;
  siForwardChannel[0][3] = 25;
  siForwardChannel[1][0] = 31;
  siForwardChannel[1][1] = 37;
  siForwardChannel[1][2] = 43;
  siForwardChannel[1][3] = 50;
  siForwardChannel[2][0] = 54;
  siForwardChannel[2][1] = 59;
  siForwardChannel[2][2] = 63;
  siForwardChannel[2][3] = 67;
  siForwardChannel[3][0] = 23;
  siForwardChannel[3][1] = 3;
  siForwardChannel[3][2] = 9;
  siForwardChannel[3][3] = 16;
  siForwardChannel[4][0] = 29;
  siForwardChannel[4][1] = 35;
  siForwardChannel[4][2] = 41;
  siForwardChannel[4][3] = 48;
  siForwardChannel[5][0] = 52;
  siForwardChannel[5][1] = 57;
  siForwardChannel[5][2] = 61;
  siForwardChannel[5][3] = 65;
  siForwardChannel[6][0] = 1;
  siForwardChannel[6][1] = 7;
  siForwardChannel[6][2] = 14;
  siForwardChannel[6][3] = 20;
  siForwardChannel[7][0] = 46;
  siForwardChannel[7][1] = 27;
  siForwardChannel[7][2] = 33;
  siForwardChannel[7][3] = 39;
  siForwardChannel[8][0] = 4;
  siForwardChannel[8][1] = 10;
  siForwardChannel[8][2] = 17;
  siForwardChannel[8][3] = 24;
  siForwardChannel[9][0] = 30;
  siForwardChannel[9][1] = 36;
  siForwardChannel[9][2] = 42;
  siForwardChannel[9][3] = 49;

  ////////////////////////////
  // Beam Left Si Detectors //
  ////////////////////////////

  // Detector 1
  siLeftMap[5] = make_pair(1, 1);
  siLeftMap[12] = make_pair(1, 2);
  siLeftMap[18] = make_pair(1, 3);
  siLeftMap[25] = make_pair(1, 4);

  // Detector 2
  siLeftMap[50] = make_pair(2, 1);
  siLeftMap[31] = make_pair(2, 2);
  siLeftMap[37] = make_pair(2, 3);
  siLeftMap[43] = make_pair(2, 4);

  // Detector 3
  siLeftMap[54] = make_pair(3, 1);
  siLeftMap[59] = make_pair(3, 2);
  siLeftMap[63] = make_pair(3, 3);
  siLeftMap[67] = make_pair(3, 4);

  // Detector 4
  siLeftMap[23] = make_pair(4, 1);
  siLeftMap[3] = make_pair(4, 2);
  siLeftMap[9] = make_pair(4, 3);
  siLeftMap[16] = make_pair(4, 4);

  // Detector 5
  siLeftMap[29] = make_pair(5, 1);
  siLeftMap[35] = make_pair(5, 2);
  siLeftMap[41] = make_pair(5, 3);
  siLeftMap[48] = make_pair(5, 4);

  // Detector 6
  siLeftMap[65] = make_pair(6, 1);
  siLeftMap[52] = make_pair(6, 2);
  siLeftMap[57] = make_pair(6, 3);
  siLeftMap[61] = make_pair(6, 4);

  siLeftChannel[0][0] = 5;
  siLeftChannel[0][1] = 12;
  siLeftChannel[0][2] = 18;
  siLeftChannel[0][3] = 25;
  siLeftChannel[1][0] = 50;
  siLeftChannel[1][1] = 31;
  siLeftChannel[1][2] = 37;
  siLeftChannel[1][3] = 43;
  siLeftChannel[2][0] = 54;
  siLeftChannel[2][1] = 59;
  siLeftChannel[2][2] = 63;
  siLeftChannel[2][3] = 67;
  siLeftChannel[3][0] = 23;
  siLeftChannel[3][1] = 3;
  siLeftChannel[3][2] = 9;
  siLeftChannel[3][3] = 16;
  siLeftChannel[4][0] = 29;
  siLeftChannel[4][1] = 35;
  siLeftChannel[4][2] = 41;
  siLeftChannel[4][3] = 48;
  siLeftChannel[5][0] = 65;
  siLeftChannel[5][1] = 52;
  siLeftChannel[5][2] = 57;
  siLeftChannel[5][3] = 61;

  ///////////////////////////
  // Forward CsI Detectors //
  ///////////////////////////

  CsIForwardMap[2] = 1;
  CsIForwardMap[7] = 2;
  CsIForwardMap[10] = 3;
  CsIForwardMap[16] = 4;
  CsIForwardMap[19] = 5;
  CsIForwardMap[25] = 6;
  CsIForwardMap[28] = 7;
  CsIForwardMap[33] = 8;
  CsIForwardMap[36] = 9;
  CsIForwardMap[41] = 10;

  CsIForwardChannel[0] = 2;
  CsIForwardChannel[1] = 7;
  CsIForwardChannel[2] = 10;
  CsIForwardChannel[3] = 16;
  CsIForwardChannel[4] = 19;
  CsIForwardChannel[5] = 25;
  CsIForwardChannel[6] = 28;
  CsIForwardChannel[7] = 33;
  CsIForwardChannel[8] = 36;
  CsIForwardChannel[9] = 41;

  /////////////////////////////
  // Beam Left CsI Detectors //
  /////////////////////////////

  CsILeftMap[44] = 1;
  CsILeftMap[50] = 2;
  CsILeftMap[53] = 3;
  CsILeftMap[59] = 4;
  CsILeftMap[62] = 5;
  CsILeftMap[67] = 6;

  CsILeftChannel[0] = 44;
  CsILeftChannel[1] = 50;
  CsILeftChannel[2] = 53;
  CsILeftChannel[3] = 59;
  CsILeftChannel[4] = 62;
  CsILeftChannel[5] = 67;

  ////////////////////
  // Micromegas Map //
  ////////////////////

  // Asad 0 and Asad 1 are all central pads
  // Asad 2 Aget 3+4 and Asad 3 Aget 3+4 are also central pads (outside)
  // Asad 2 Aget 1+2 are beam left
  // Asad 3 Aget 1+2 are beam right

  // Aget_Map
  int j=0;
  for(int i=0; i<68; i++) {
    if(i==11 || i==22 || i==45 || i==56) continue; // FPN Channels
    Aget_Map[j] = i;
    j++;
  }

  // Asad0_Aget0
  j=31;
  for(int i=0; i<64; i++) {
    if(i%4==0 && i!=0) j--;
    int remain = i%4;
    MM_Map_Asad0_Aget0[Aget_Map[i]] = make_pair(4-remain, j);
  }

  // Asad0_Aget1
  j=63;
  for(int i=0; i<64; i++) {
    if(i%4==0 && i!=0) j--;
    int remain = i%4;
    MM_Map_Asad0_Aget1[Aget_Map[i]] = make_pair(4-remain, j);
  }

  // Asad0_Aget2
  j=95;
  for(int i=0; i<64; i++) {
    if(i%4==0 && i!=0) j--;
    int remain = i%4;
    MM_Map_Asad0_Aget2[Aget_Map[i]] = make_pair(1+remain, j);
  }

  // Asad0_Aget3
  j=127;
  for(int i=0; i<64; i++) {
    if(i%4==0 && i!=0) j--;
    int remain = i%4;
    MM_Map_Asad0_Aget3[Aget_Map[i]] = make_pair(1+remain, j);
  }

  // Asad1_Aget0
  j=15;
  for(int i=0; i<64; i++) {
    if(i%4==0 && i!=0) j--;
    int remain = i%4;
    MM_Map_Asad1_Aget0[Aget_Map[i]] = make_pair(1+remain, j);
  }

  // Asad1_Aget1
  j=47;
  for(int i=0; i<64; i++) {
    if(i%4==0 && i!=0) j--;
    int remain = i%4;
    MM_Map_Asad1_Aget1[Aget_Map[i]] = make_pair(1+remain, j);
  }

  // Asad1_Aget2
  j=79;
  for(int i=0; i<64; i++) {
    if(i%4==0 && i!=0) j--;
    int remain = i%4;
    MM_Map_Asad1_Aget2[Aget_Map[i]] = make_pair(4-remain, j);
  }

  // Asad1_Aget3
  j=111;
  for(int i=0; i<64; i++) {
    if(i%4==0 && i!=0) j--;
    int remain = i%4;
    MM_Map_Asad1_Aget3[Aget_Map[i]] = make_pair(4-remain, j);
  }

  // Asad2_Aget0
  // Strips part (Beam left)
  j=31;
  for(int i=0; i<32; i++) {
    MM_Map_Asad2_Aget0[Aget_Map[i]] = j;
    j--;
  }
  // Chains part (Beam left)
  j=32;
  for(int i=32; i<64; i++) {
    MM_Map_Asad2_Aget0[Aget_Map[i]] = j;
    j++;
  }

  // Asad2_Aget1
  // Strips part (Beam left)
  j=63;
  for(int i=0; i<32; i++) {
    MM_Map_Asad2_Aget1[Aget_Map[i]] = j;
    j--;
  }
  // Chains part (Beam left)
  j=0;
  for(int i=32; i<64; i++) {
    MM_Map_Asad2_Aget1[Aget_Map[i]] = j;
    j++;
  }

  // Asad2_Aget2
  j=63;
  for(int i=0; i<64; i++) {
    MM_Map_Asad2_Aget2[Aget_Map[i]] = make_pair(0,j);
    j--;
  }

  // Asad2_Aget3
  map<int, pair<int, int> > MM_Map_Asad2_Aget3;
  j=127;
  for(int i=0; i<64; i++) {
    MM_Map_Asad2_Aget3[Aget_Map[i]] = make_pair(0,j);
    j--;
  }

  // Asad3_Aget0
  j=31;
  for(int i=0; i<32; i++) {
    MM_Map_Asad3_Aget0[Aget_Map[i]] = j;
    j--;
  }
  // Chains part (Beam right)
  j=32;
  for(int i=32; i<64; i++) {
    MM_Map_Asad3_Aget0[Aget_Map[i]] = j;
    j++;
  }

  // Asad3_Aget1
  j=63;
  for(int i=0; i<32; i++) {
    MM_Map_Asad3_Aget1[Aget_Map[i]] = j;
    j--;
  }
  // Chains part (Beam right)
  j=0;
  for(int i=32; i<64; i++) {
    MM_Map_Asad3_Aget1[Aget_Map[i]] = j;
    j++;
  }

  // Asad3_Aget2
  j=63;
  for(int i=0; i<64; i++) {
    MM_Map_Asad3_Aget2[Aget_Map[i]] = make_pair(5,j);
    j--;
  }

  // Asad3_Aget3
  map<int, pair<int, int> > MM_Map_Asad3_Aget3;
  j=127;
  for(int i=0; i<64; i++) {
    MM_Map_Asad3_Aget3[Aget_Map[i]] = make_pair(5,j);
    j--;
  }

  // Asad3_Aget3

}

#endif // #ifdef Spectra_cxx
