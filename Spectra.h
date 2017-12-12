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
#include <TF2.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <TPolyLine3D.h>
#include <TROOT.h>
#include <TStyle.h>

#include <Fit/Fitter.h>

#include <Math/Functor.h>
#include <Math/Vector3D.h>

#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

#include "CubicSpline.h"
#include "EnergyLoss.h"

// Define the parameteric line equation
void line(double t, const double *p, double &x, double &y, double &z) {
  // a parameteric line is define from 6 parameters but 4 are independent
  // x0,y0,z0,z1,y1,z1 which are the coordinates of two points on the line
  // can choose z0 = 0 if line not parallel to x-y plane and z1 = 1
  x = p[0] + p[1]*t;
  y = p[2] + p[3]*t;
  z = t;
}

// Calculate distance line-point
double distance2(double x, double y, double z, const double *p) {
  // distance line point is D = | (xp-x0) cross ux |
  // where ux is direction of line and x0 is a point in the line (like t = 0)
  ROOT::Math::XYZVector xp(x, y, z);
  ROOT::Math::XYZVector x0(p[0], p[2], 0);
  ROOT::Math::XYZVector x1(p[0] + p[1], p[2] + p[3], 1.);
  ROOT::Math::XYZVector u = (x1-x0).Unit();
  double d2 = ((xp-x0).Cross(u)).Mag2();
  return d2;
}

// function Object to be minimized
struct SumDistance2 {
  // the TGraph is a data memeber of the object
  TGraph2D *fGraph;
  SumDistance2(TGraph2D *g) : fGraph(g) {}

  // implementation of the function to be minimized
  double operator() (const double *par) {
    assert(fGraph != 0);
    double *x = fGraph->GetX();
    double *y = fGraph->GetY();
    double *z = fGraph->GetZ();
    int npoints = fGraph->GetN();
    double sum = 0;
    for(int i = 0; i < npoints ; i++) {
      double d = distance2(x[i], y[i], z[i], par);
      sum += d;
    }
    return sum;
  }
};


typedef struct siDetect {
  int detect;
  int quad;
  int channel;
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

typedef struct mmCenterTrack {
  double position;
  int row;
  double time;
  double energy;
  double height;
  int total;
} mmCenterTrack;


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
  void InitHistograms();

  TH1F* hSiEForwardTotal[10];
  TH1F* hSiEForwardTotalCal[10];
  TH1F* hSiTForwardTotal[10];

  TH1F* hSiEForward[10][4];
  TH1F* hSiEForwardCal[10][4];
  TH1F* hSiTForward[10][4];

  TH1F* hSiELeftTotal[6];
  TH1F* hSiELeftTotalCal[6];
  TH1F* hSiTLeftTotal[6];

  TH1F* hSiELeft[6][4];
  TH1F* hSiELeftCal[6][4];
  TH1F* hSiTLeft[6][4];

  TH1F* hCsIEForward[10];
  TH1F* hCsIEForwardCal[10];
  TH1F* hCsITForward[10];

  TH1F* hCsIELeft[6];
  TH1F* hCsIELeftCal[6];
  TH1F* hCsITLeft[6];

  TH1F* hCSDet5;
  TH1F* hCSDet5Counts;

// Silicon Energy Calibration
private:
  void InitSiEForwardCalibration();
  std::pair<double, double> siEForwardCalibration[10][4] = {std::make_pair(0., 0.)};

// Center Pad Gain Match
private:
  void InitCentralPadGainMatch();
  double scale[6][128];

// Beam Average Central Pads
  void InitAverageBeamEnergy();
  double averageBeamEnergy[128];

// Cross Section
private:
  void SimpleCrossSection(TH1F*);
  void SimpleSolidAngleDet5(TH1F*);
  Double_t CalcSimpleSolidAngleDet5(Double_t);
  void DivideTargetThickness(TH1F*);

// General Variables
private:
  void InitVariables();
  Double_t m1;
  Double_t m2;
  Double_t zeroTime;
  Double_t driftVelocity;
  Double_t beamEnergy;
  Double_t density;
  Double_t distanceHavarToSilicon;
  Double_t numberB8;
  EnergyLoss *boronMethane;
  EnergyLoss *protonMethane;

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
  printf("Initializing Channel Map\n");

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
    MM_Map_Asad2_Aget2[Aget_Map[i]] = std::make_pair(0, j);
    j--;
  }

  // Asad2_Aget3
  j=127;
  for(int i=0; i<64; i++) {
    MM_Map_Asad2_Aget3[Aget_Map[i]] = std::make_pair(0, j);
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
    MM_Map_Asad3_Aget2[Aget_Map[i]] = std::make_pair(5, j);
    j--;
  }

  // Asad3_Aget3
  j=127;
  for(int i=0; i<64; i++) {
    MM_Map_Asad3_Aget3[Aget_Map[i]] = std::make_pair(5, j);
    j--;
  }

}

inline void Spectra::InitHistograms() {
  printf("Initializing Histograms\n");

  // Histograms for the Forward Si Detectors
  for(unsigned int i=0; i<10; i++) {
    TString name = Form("siEForward_d%d", i + 1);
    hSiEForwardTotal[i] = new TH1F(name, name, 500, 0, 4000);
    hSiEForwardTotal[i]->GetXaxis()->SetTitle("Channels"); hSiEForwardTotal[i]->GetXaxis()->CenterTitle();
    hSiEForwardTotal[i]->GetYaxis()->SetTitle("Counts"); hSiEForwardTotal[i]->GetYaxis()->CenterTitle();
    hSiEForwardTotal[i]->GetYaxis()->SetTitleOffset(1.4);

    name = Form("siEForwardCal_d%d", i + 1);
    hSiEForwardTotalCal[i] = new TH1F(name, name, 2000, 0, 16000);
    hSiEForwardTotalCal[i]->GetXaxis()->SetTitle("Energy [keV]"); hSiEForwardTotalCal[i]->GetXaxis()->CenterTitle();
    hSiEForwardTotalCal[i]->GetYaxis()->SetTitle("Counts"); hSiEForwardTotalCal[i]->GetYaxis()->CenterTitle();
    hSiEForwardTotalCal[i]->GetYaxis()->SetTitleOffset(1.4);

    name = Form("siTForward_d%d", i + 1);
    hSiTForwardTotal[i] = new TH1F(name, name, 2500, 0, 20000);
    hSiTForwardTotal[i]->GetXaxis()->SetTitle("Channels"); hSiTForwardTotal[i]->GetXaxis()->CenterTitle();
    hSiTForwardTotal[i]->GetYaxis()->SetTitle("Counts"); hSiTForwardTotal[i]->GetYaxis()->CenterTitle();
    hSiTForwardTotal[i]->GetYaxis()->SetTitleOffset(1.4);

    for(unsigned int j=0; j<4; j++) {
      TString name = Form("siEForward_d%d_q%d_ch_%d", i + 1, j + 1, siForwardChannel[i][j]);
      hSiEForward[i][j] = new TH1F(name, name, 500, 0, 4000);
      hSiEForward[i][j]->GetXaxis()->SetTitle("Channels"); hSiEForward[i][j]->GetXaxis()->CenterTitle();
      hSiEForward[i][j]->GetYaxis()->SetTitle("Counts"); hSiEForward[i][j]->GetYaxis()->CenterTitle();
      hSiEForward[i][j]->GetYaxis()->SetTitleOffset(1.4);

      name = Form("siEForwardCal_d%d_q%d_ch_%d", i + 1, j + 1, siForwardChannel[i][j]);
      hSiEForwardCal[i][j] = new TH1F(name, name, 500, 0, 4000);
      hSiEForwardCal[i][j]->GetXaxis()->SetTitle("Channels"); hSiEForwardCal[i][j]->GetXaxis()->CenterTitle();
      hSiEForwardCal[i][j]->GetYaxis()->SetTitle("Counts"); hSiEForwardCal[i][j]->GetYaxis()->CenterTitle();
      hSiEForwardCal[i][j]->GetYaxis()->SetTitleOffset(1.4);

      name = Form("siTForward_d%d_q%d_ch_%d", i + 1, j + 1, siForwardChannel[i][j]);
      hSiTForward[i][j] = new TH1F(name, name, 500, 0, 4000);
      hSiTForward[i][j]->GetXaxis()->SetTitle("Channels"); hSiTForward[i][j]->GetXaxis()->CenterTitle();
      hSiTForward[i][j]->GetYaxis()->SetTitle("Counts"); hSiTForward[i][j]->GetYaxis()->CenterTitle();
      hSiTForward[i][j]->GetYaxis()->SetTitleOffset(1.4);
    }
  }

  // Histograms for the Beam Left Si Detectors
  for(unsigned int i=0; i<6; i++) {
    if(i > 7) break;
    TString name = Form("siELeft_d%d", i + 1);
    hSiELeftTotal[i] = new TH1F(name, name, 500, 0, 4000);
    hSiELeftTotal[i]->GetXaxis()->SetTitle("Channels"); hSiELeftTotal[i]->GetXaxis()->CenterTitle();
    hSiELeftTotal[i]->GetYaxis()->SetTitle("Counts"); hSiELeftTotal[i]->GetYaxis()->CenterTitle();
    hSiELeftTotal[i]->GetYaxis()->SetTitleOffset(1.4);

    name = Form("siELeftCal_d%d", i + 1);
    hSiELeftTotalCal[i] = new TH1F(name, name, 2000, 0, 16000);
    hSiELeftTotalCal[i]->GetXaxis()->SetTitle("Energy [keV]"); hSiELeftTotalCal[i]->GetXaxis()->CenterTitle();
    hSiELeftTotalCal[i]->GetYaxis()->SetTitle("Counts"); hSiELeftTotalCal[i]->GetYaxis()->CenterTitle();
    hSiELeftTotalCal[i]->GetYaxis()->SetTitleOffset(1.4);

    name = Form("siTLeft_d%d", i + 1);
    hSiTLeftTotal[i] = new TH1F(name, name, 2500, 0, 20000);
    hSiTLeftTotal[i]->GetXaxis()->SetTitle("Channels"); hSiTLeftTotal[i]->GetXaxis()->CenterTitle();
    hSiTLeftTotal[i]->GetYaxis()->SetTitle("Counts"); hSiTLeftTotal[i]->GetYaxis()->CenterTitle();
    hSiTLeftTotal[i]->GetYaxis()->SetTitleOffset(1.4);

    for(unsigned int j=0; j<4; j++) {
      TString name = Form("siELeft_d%d_q%d_ch_%d", i + 1, j + 1, siLeftChannel[i][j]);
      hSiELeft[i][j] = new TH1F(name, name, 500, 0, 4000);
      hSiELeft[i][j]->GetXaxis()->SetTitle("Channels"); hSiELeft[i][j]->GetXaxis()->CenterTitle();
      hSiELeft[i][j]->GetYaxis()->SetTitle("Counts"); hSiELeft[i][j]->GetYaxis()->CenterTitle();
      hSiELeft[i][j]->GetYaxis()->SetTitleOffset(1.4);

      name = Form("siELeftCal_d%d_q%d_ch_%d", i + 1, j + 1, siLeftChannel[i][j]);
      hSiELeftCal[i][j] = new TH1F(name, name, 500, 0, 4000);
      hSiELeftCal[i][j]->GetXaxis()->SetTitle("Channels"); hSiELeftCal[i][j]->GetXaxis()->CenterTitle();
      hSiELeftCal[i][j]->GetYaxis()->SetTitle("Counts"); hSiELeftCal[i][j]->GetYaxis()->CenterTitle();
      hSiELeftCal[i][j]->GetYaxis()->SetTitleOffset(1.4);

      name = Form("siTLeft_d%d_q%d_ch_%d", i + 1, j + 1, siLeftChannel[i][j]);
      hSiTLeft[i][j] = new TH1F(name, name, 500, 0, 4000);
      hSiTLeft[i][j]->GetXaxis()->SetTitle("Channels"); hSiTLeft[i][j]->GetXaxis()->CenterTitle();
      hSiTLeft[i][j]->GetYaxis()->SetTitle("Counts"); hSiTLeft[i][j]->GetYaxis()->CenterTitle();
      hSiTLeft[i][j]->GetYaxis()->SetTitleOffset(1.4);
    }
  }

  // Histograms for the Forward CsI Detectors
  for(unsigned int i=0; i<10; i++) {
    TString name = Form("CsIEForward_d%d_ch%d", i + 1, csiForwardChannel[i]);
    hCsIEForward[i] = new TH1F(name, name, 1000, 0, 4000);
    hCsIEForward[i]->GetXaxis()->SetTitle("Channels"); hCsIEForward[i]->GetXaxis()->CenterTitle();
    hCsIEForward[i]->GetYaxis()->SetTitle("Counts"); hCsIEForward[i]->GetYaxis()->CenterTitle();
    hCsIEForward[i]->GetYaxis()->SetTitleOffset(1.4);

    name = Form("CsIEForwardCal_d%d_ch%d", i + 1, csiForwardChannel[i]);
    hCsIEForwardCal[i] = new TH1F(name, name, 1000, 0, 4000);
    hCsIEForwardCal[i]->GetXaxis()->SetTitle("Channels"); hCsIEForwardCal[i]->GetXaxis()->CenterTitle();
    hCsIEForwardCal[i]->GetYaxis()->SetTitle("Counts"); hCsIEForwardCal[i]->GetYaxis()->CenterTitle();
    hCsIEForwardCal[i]->GetYaxis()->SetTitleOffset(1.4);

    name = Form("CsITForward_d%d_ch%d", i + 1, csiForwardChannel[i]);
    hCsITForward[i] = new TH1F(name, name, 1000, 0, 20000);
    hCsITForward[i]->GetXaxis()->SetTitle("Channels"); hCsITForward[i]->GetXaxis()->CenterTitle();
    hCsITForward[i]->GetYaxis()->SetTitle("Counts"); hCsITForward[i]->GetYaxis()->CenterTitle();
    hCsITForward[i]->GetYaxis()->SetTitleOffset(1.4);
  }

  // Histograms for the Beam Left CsI Detectors
  for(unsigned int i=0; i<6; i++) {
    TString name = Form("CsIELeft_d%d_ch%d", i + 1, csiLeftChannel[i]);
    hCsIELeft[i] = new TH1F(name,name,1000,0,4000);
    hCsIELeft[i]->GetXaxis()->SetTitle("Channels"); hCsIELeft[i]->GetXaxis()->CenterTitle();
    hCsIELeft[i]->GetYaxis()->SetTitle("Counts"); hCsIELeft[i]->GetYaxis()->CenterTitle();
    hCsIELeft[i]->GetYaxis()->SetTitleOffset(1.4);

    name = Form("CsIELeftCal_d%d_ch%d", i + 1, csiLeftChannel[i]);
    hCsIELeftCal[i] = new TH1F(name,name,1000,0,4000);
    hCsIELeftCal[i]->GetXaxis()->SetTitle("Channels"); hCsIELeftCal[i]->GetXaxis()->CenterTitle();
    hCsIELeftCal[i]->GetYaxis()->SetTitle("Counts"); hCsIELeftCal[i]->GetYaxis()->CenterTitle();
    hCsIELeftCal[i]->GetYaxis()->SetTitleOffset(1.4);

    name = Form("CsITLeft_d%d_ch%d", i + 1, csiLeftChannel[i]);
    hCsITLeft[i] = new TH1F(name,name,1000,0,20000);
    hCsITLeft[i]->GetXaxis()->SetTitle("Channels"); hCsITLeft[i]->GetXaxis()->CenterTitle();
    hCsITLeft[i]->GetYaxis()->SetTitle("Counts"); hCsITLeft[i]->GetYaxis()->CenterTitle();
    hCsITLeft[i]->GetYaxis()->SetTitleOffset(1.4);
  }
}

inline void Spectra::InitSiEForwardCalibration() {
  printf("Reading Si Energy Forward Calibrations File\n");
  std::ifstream inSiCalFile("siCalibration.dat");
  assert(inSiCalFile.is_open());
  int var1, var2, var3;
  double slope, intercept;
  while(inSiCalFile >> var1 >> var2 >> var3 >> slope >> intercept) {
    siEForwardCalibration[var1-1][var2-1] = std::make_pair(slope, intercept);
  }
  inSiCalFile.close();
}

inline void Spectra::InitCentralPadGainMatch() {
  printf("Reading Central Pad Gain Matching File\n");
  std::ifstream inGainFile("gainFile.dat");
  assert(inGainFile.is_open());
  int varI, varJ;
  double varScale;
  while(inGainFile >> varI >> varJ >> varScale) {
    scale[varI][varJ] = varScale;
  }
  inGainFile.close();
}

inline void Spectra::InitAverageBeamEnergy() {
  printf("Reading Average Beam Energy File\n");
  std::ifstream inBeamFile("averageBeamEnergy.out");
  assert(inBeamFile.is_open());
  int varJ;
  double varE;
  while(inBeamFile >> varJ >> varE) {
    averageBeamEnergy[varJ] = varE;
  }
}

inline void Spectra::InitVariables() {
  printf("Initializing Variables\n");

  m1 = 8.; // AMU of projectile
  m2 = 1.; // AMU of target

  // zeroTime = 6081.81; // Time of 0 height in chamber (in ns)
  zeroTime = 6081.81 - 352.457; // Time of 0 height in chamber (in ns)
  driftVelocity = 0.05674449; // in mm/ns

  beamEnergy = 56.; // In MeV, after havar window
  density = 0.00038175; // in g/cm3, from LISE++ (Methane at 435 torr)
  numberB8 = 174809089.;

  distanceHavarToSilicon = 544.07; // Distance from Havar to Forward Silicon in mm

  // Initialize EnergyLoss
  boronMethane = new EnergyLoss("b8_methane.dat");
  protonMethane = new EnergyLoss("proton_methane.dat");
}

#endif // #ifdef Spectra_cxx
