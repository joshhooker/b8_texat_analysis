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
#include <TCutG.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>

#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

#include "CubicSpline.h"
#include "EnergyLoss.h"

typedef struct siDetect {
  int detect;
  int quad;
  double energy;
  double time;
} siHit;

typedef struct csiDetect {
  int detect;
  double energy;
  double time;
} csiHit;

typedef struct mmCenter {
  int column;
  int row;
  double energy;
  double time;
} mmCenter;

typedef struct mmStrip {
  int row;
  double energy;
  double time;
} mmStrip;

typedef struct mmChain {
  int column;
  double energy;
  double time;
} mmChain;

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
  ~Spectra();
  Int_t GetEntry(Long64_t entry);
  Long64_t LoadTree(Long64_t entry);
  void Init(TTree *tree);
  void Loop();
  Bool_t Notify();
  void Show(Long64_t entry = -1);

// Channel Map
private:
  void InitChannelMap();

  std::map<int, std::pair<int, int> > siForwardMap;
  int siForwardChannel[10][4];

  std::map<int, std::pair<int, int> > siLeftMap;
  int siLeftChannel[6][4];

  std::map<int, int> csiForwardMap;
  int csiForwardChannel[10];

  std::map<int, int> csiLeftMap;
  int csiLeftChannel[6];

  std::map<int, int> Aget_Map;

  std::map<int, std::pair<int, int> > MM_Map_Asad0_Aget0;
  std::map<int, std::pair<int, int> > MM_Map_Asad0_Aget1;
  std::map<int, std::pair<int, int> > MM_Map_Asad0_Aget2;
  std::map<int, std::pair<int, int> > MM_Map_Asad0_Aget3;

  std::map<int, std::pair<int, int> > MM_Map_Asad1_Aget0;
  std::map<int, std::pair<int, int> > MM_Map_Asad1_Aget1;
  std::map<int, std::pair<int, int> > MM_Map_Asad1_Aget2;
  std::map<int, std::pair<int, int> > MM_Map_Asad1_Aget3;

  std::map<int, int> MM_Map_Asad2_Aget0;
  std::map<int, int> MM_Map_Asad2_Aget1;
  std::map<int, std::pair<int, int> > MM_Map_Asad2_Aget2;
  std::map<int, std::pair<int, int> > MM_Map_Asad2_Aget3;

  std::map<int, int> MM_Map_Asad3_Aget0;
  std::map<int, int> MM_Map_Asad3_Aget1;
  std::map<int, std::pair<int, int> > MM_Map_Asad3_Aget2;
  std::map<int, std::pair<int, int> > MM_Map_Asad3_Aget3;

// Histograms
private:
  TH1F* hSiEForwardTotal[10];
  TH1F* hSiTForwardTotal[10];
  TH1F* hSiEForward[10][4];
  TH1F* hSiEForwardCal[10][4];
  TH1F* hSiTForward[10][4];

  TH1F* hSiELeftTotal[8];
  TH1F* hSiELeft[8][4];
  TH1F* hSiELeftCal[8][4];
  TH1F* hSiTLeft[8][4];

  TH1F* hCsIEForward[10];
  TH1F* hCsITForward[10];

  TH1F* hCsIELeft[8];
  TH1F* hCsITLeft[8];

  TH1F* hCSDet5;
  TH1F* hCSDet5Counts;

// General Variables
private:
  void Initialize();
  Double_t m1;
  Double_t m2;
  Double_t beamEnergy;
  Double_t density;
  Double_t distanceHavarToSilicon;
  Double_t numberB8;
  EnergyLoss *boronMethane;
  EnergyLoss *protonMethane;

// Cross Section
private:
  void SimpleCrossSection(TH1F*);
  void SimpleSolidAngleDet5(TH1F*);
  Double_t CalcSimpleSolidAngleDet5(Double_t);
  void DivideTargetThickness(TH1F*);

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
  siForwardMap[5] = std::make_pair(1, 1);
  siForwardMap[12] = std::make_pair(1, 2);
  siForwardMap[18] = std::make_pair(1, 3);
  siForwardMap[25] = std::make_pair(1, 4);

  // Detector 2
  siForwardMap[31] = std::make_pair(2, 1);
  siForwardMap[37] = std::make_pair(2, 2);
  siForwardMap[43] = std::make_pair(2, 3);
  siForwardMap[50] = std::make_pair(2, 4);

  // Detector 3
  siForwardMap[54] = std::make_pair(3, 1);
  siForwardMap[59] = std::make_pair(3, 2);
  siForwardMap[63] = std::make_pair(3, 3);
  siForwardMap[67] = std::make_pair(3, 4);

  // Detector 4
  siForwardMap[23] = std::make_pair(4, 1);
  siForwardMap[3] = std::make_pair(4, 2);
  siForwardMap[9] = std::make_pair(4, 3);
  siForwardMap[16] = std::make_pair(4, 4);

  // Detector 5
  siForwardMap[29] = std::make_pair(5, 1);
  siForwardMap[35] = std::make_pair(5, 2);
  siForwardMap[41] = std::make_pair(5, 3);
  siForwardMap[48] = std::make_pair(5, 4);

  // Detector 6
  siForwardMap[52] = std::make_pair(6, 1);
  siForwardMap[57] = std::make_pair(6, 2);
  siForwardMap[61] = std::make_pair(6, 3);
  siForwardMap[65] = std::make_pair(6, 4);

  // Detector 7
  siForwardMap[1] = std::make_pair(7, 1);
  siForwardMap[7] = std::make_pair(7, 2);
  siForwardMap[14] = std::make_pair(7, 3);
  siForwardMap[20] = std::make_pair(7, 4);

  // Detector 8
  siForwardMap[46] = std::make_pair(8, 1);
  siForwardMap[27] = std::make_pair(8, 2);
  siForwardMap[33] = std::make_pair(8, 3);
  siForwardMap[39] = std::make_pair(8, 4);

  // Detector 9
  siForwardMap[4] = std::make_pair(9, 1);
  siForwardMap[10] = std::make_pair(9, 2);
  siForwardMap[17] = std::make_pair(9, 3);
  siForwardMap[24] = std::make_pair(9, 4);

  // Detector 10
  siForwardMap[30] = std::make_pair(10, 1);
  siForwardMap[36] = std::make_pair(10, 2);
  siForwardMap[42] = std::make_pair(10, 3);
  siForwardMap[49] = std::make_pair(10, 4);

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
  siLeftMap[5] = std::make_pair(1, 1);
  siLeftMap[12] = std::make_pair(1, 2);
  siLeftMap[18] = std::make_pair(1, 3);
  siLeftMap[25] = std::make_pair(1, 4);

  // Detector 2
  siLeftMap[50] = std::make_pair(2, 1);
  siLeftMap[31] = std::make_pair(2, 2);
  siLeftMap[37] = std::make_pair(2, 3);
  siLeftMap[43] = std::make_pair(2, 4);

  // Detector 3
  siLeftMap[54] = std::make_pair(3, 1);
  siLeftMap[59] = std::make_pair(3, 2);
  siLeftMap[63] = std::make_pair(3, 3);
  siLeftMap[67] = std::make_pair(3, 4);

  // Detector 4
  siLeftMap[23] = std::make_pair(4, 1);
  siLeftMap[3] = std::make_pair(4, 2);
  siLeftMap[9] = std::make_pair(4, 3);
  siLeftMap[16] = std::make_pair(4, 4);

  // Detector 5
  siLeftMap[29] = std::make_pair(5, 1);
  siLeftMap[35] = std::make_pair(5, 2);
  siLeftMap[41] = std::make_pair(5, 3);
  siLeftMap[48] = std::make_pair(5, 4);

  // Detector 6
  siLeftMap[65] = std::make_pair(6, 1);
  siLeftMap[52] = std::make_pair(6, 2);
  siLeftMap[57] = std::make_pair(6, 3);
  siLeftMap[61] = std::make_pair(6, 4);

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

  csiForwardMap[2] = 1;
  csiForwardMap[7] = 2;
  csiForwardMap[10] = 3;
  csiForwardMap[16] = 4;
  csiForwardMap[19] = 5;
  csiForwardMap[25] = 6;
  csiForwardMap[28] = 7;
  csiForwardMap[33] = 8;
  csiForwardMap[36] = 9;
  csiForwardMap[41] = 10;

  csiForwardChannel[0] = 2;
  csiForwardChannel[1] = 7;
  csiForwardChannel[2] = 10;
  csiForwardChannel[3] = 16;
  csiForwardChannel[4] = 19;
  csiForwardChannel[5] = 25;
  csiForwardChannel[6] = 28;
  csiForwardChannel[7] = 33;
  csiForwardChannel[8] = 36;
  csiForwardChannel[9] = 41;

  /////////////////////////////
  // Beam Left CsI Detectors //
  /////////////////////////////

  csiLeftMap[44] = 1;
  csiLeftMap[50] = 2;
  csiLeftMap[53] = 3;
  csiLeftMap[59] = 4;
  csiLeftMap[62] = 5;
  csiLeftMap[67] = 6;

  csiLeftChannel[0] = 44;
  csiLeftChannel[1] = 50;
  csiLeftChannel[2] = 53;
  csiLeftChannel[3] = 59;
  csiLeftChannel[4] = 62;
  csiLeftChannel[5] = 67;

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
    MM_Map_Asad0_Aget0[Aget_Map[i]] = std::make_pair(4-remain, j);
  }

  // Asad0_Aget1
  j=63;
  for(int i=0; i<64; i++) {
    if(i%4==0 && i!=0) j--;
    int remain = i%4;
    MM_Map_Asad0_Aget1[Aget_Map[i]] = std::make_pair(4-remain, j);
  }

  // Asad0_Aget2
  j=95;
  for(int i=0; i<64; i++) {
    if(i%4==0 && i!=0) j--;
    int remain = i%4;
    MM_Map_Asad0_Aget2[Aget_Map[i]] = std::make_pair(1+remain, j);
  }

  // Asad0_Aget3
  j=127;
  for(int i=0; i<64; i++) {
    if(i%4==0 && i!=0) j--;
    int remain = i%4;
    MM_Map_Asad0_Aget3[Aget_Map[i]] = std::make_pair(1+remain, j);
  }

  // Asad1_Aget0
  j=15;
  for(int i=0; i<64; i++) {
    if(i%4==0 && i!=0) j--;
    int remain = i%4;
    MM_Map_Asad1_Aget0[Aget_Map[i]] = std::make_pair(1+remain, j);
  }

  // Asad1_Aget1
  j=47;
  for(int i=0; i<64; i++) {
    if(i%4==0 && i!=0) j--;
    int remain = i%4;
    MM_Map_Asad1_Aget1[Aget_Map[i]] = std::make_pair(1+remain, j);
  }

  // Asad1_Aget2
  j=79;
  for(int i=0; i<64; i++) {
    if(i%4==0 && i!=0) j--;
    int remain = i%4;
    MM_Map_Asad1_Aget2[Aget_Map[i]] = std::make_pair(4-remain, j);
  }

  // Asad1_Aget3
  j=111;
  for(int i=0; i<64; i++) {
    if(i%4==0 && i!=0) j--;
    int remain = i%4;
    MM_Map_Asad1_Aget3[Aget_Map[i]] = std::make_pair(4-remain, j);
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
    MM_Map_Asad2_Aget2[Aget_Map[i]] = std::make_pair(0,j);
    j--;
  }

  // Asad2_Aget3
  j=127;
  for(int i=0; i<64; i++) {
    MM_Map_Asad2_Aget3[Aget_Map[i]] = std::make_pair(0,j);
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
    MM_Map_Asad3_Aget2[Aget_Map[i]] = std::make_pair(5,j);
    j--;
  }

  // Asad3_Aget3
  j=127;
  for(int i=0; i<64; i++) {
    MM_Map_Asad3_Aget3[Aget_Map[i]] = std::make_pair(5,j);
    j--;
  }

  // Asad3_Aget3

}

#endif // #ifdef Spectra_cxx
