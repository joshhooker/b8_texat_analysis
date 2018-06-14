//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Aug 10 09:32:04 2017 by ROOT version 6.08/06
// from TTree mfmData/Experimental Event Data
// found on file: run1_0.root
//////////////////////////////////////////////////////////

#ifndef SPECTRA_H
#define SPECTRA_H

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
#include <TVector3.h>

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
#include "FitTrack.h"
#include "Hough2D.h"
#include "TypeDef.h"

// Define the parameteric line equation
void line(Double_t t, std::vector<Double_t> p, Double_t &x, Double_t &y, Double_t &z) {
  x = p[0] + p[1]*t;
  y = t;
  z = p[2] + p[3]*t;
}

struct sortByRowMMTrack {
  inline Bool_t operator() (const mmTrack& struct1, const mmTrack& struct2) {
    return (struct1.row < struct2.row);
  }
};

struct sortByRowMMStripChain {
  inline Bool_t operator() (const mmStripChain& struct1, const mmStripChain& struct2) {
    return (struct1.row < struct2.row);
  }
};

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

  std::map<Int_t, std::pair<Int_t, Int_t> > siForwardMap;
  Int_t siForwardChannel[10][4];

  std::map<Int_t, std::pair<Int_t, Int_t> > siLeftMap;
  Int_t siLeftChannel[6][4];

  std::map<Int_t, Int_t> csiForwardMap;
  Int_t csiForwardChannel[10];

  std::map<Int_t, Int_t> csiLeftMap;
  Int_t csiLeftChannel[6];

  std::map<Int_t, Int_t> Aget_Map;

  std::map<Int_t, std::pair<Int_t, Int_t> > MM_Map_Asad0_Aget0;
  std::map<Int_t, std::pair<Int_t, Int_t> > MM_Map_Asad0_Aget1;
  std::map<Int_t, std::pair<Int_t, Int_t> > MM_Map_Asad0_Aget2;
  std::map<Int_t, std::pair<Int_t, Int_t> > MM_Map_Asad0_Aget3;

  std::map<Int_t, std::pair<Int_t, Int_t> > MM_Map_Asad1_Aget0;
  std::map<Int_t, std::pair<Int_t, Int_t> > MM_Map_Asad1_Aget1;
  std::map<Int_t, std::pair<Int_t, Int_t> > MM_Map_Asad1_Aget2;
  std::map<Int_t, std::pair<Int_t, Int_t> > MM_Map_Asad1_Aget3;

  std::map<Int_t, Int_t> MM_Map_Asad2_Aget0;
  std::map<Int_t, Int_t> MM_Map_Asad2_Aget1;
  std::map<Int_t, std::pair<Int_t, Int_t> > MM_Map_Asad2_Aget2;
  std::map<Int_t, std::pair<Int_t, Int_t> > MM_Map_Asad2_Aget3;

  std::map<Int_t, Int_t> MM_Map_Asad3_Aget0;
  std::map<Int_t, Int_t> MM_Map_Asad3_Aget1;
  std::map<Int_t, std::pair<Int_t, Int_t> > MM_Map_Asad3_Aget2;
  std::map<Int_t, std::pair<Int_t, Int_t> > MM_Map_Asad3_Aget3;

// Histograms
private:
  void InitHistograms();

  // Ionization Chamber
  TH1F* hIonizationChamberE;
  TH1F* hIonizationChamberT;

  // Micromegas - Central
  TH2F* hMicroMegasCenterCumulative;
  TH2F* hMicroMegasCenterCumulativePosition;
  TH2F* hMicroMegasCenterEnergyCumulative;
  TH2F* hMicroMegasCenterEnergyAverage;
  TH2F* hMicroMegasCenterEnergyAverageScaled;
  TH2F* hMicroMegasCenterTimeAverage;

  // Micromegas - Left
  TH1F* hMicroMegasStripLeftCumulative;
  TH1F* hMicroMegasChainLeftCumulative;

  // Micromegas - Right
  TH1F* hMicroMegasStripRightCumulative;
  TH1F* hMicroMegasChainRightCumulative;

  // Forward Si Detectors
  TH1F* hSiEForwardDet[10];
  TH1F* hSiEForwardDetCal[10];
  TH1F* hSiTForwardDet[10];
  TH1F* hSiEForward[10][4];
  TH1F* hSiEForwardCal[10][4];
  TH1F* hSiTForward[10][4];

  // Left Si Detectors
  TH1F* hSiELeftDet[6];
  TH1F* hSiELeftDetCal[6];
  TH1F* hSiTLeftDet[6];
  TH1F* hSiELeft[6][4];
  TH1F* hSiELeftCal[6][4];
  TH1F* hSiTLeft[6][4];

  // Forward CsI Detectors
  TH1F* hCsIEForward[10];
  TH1F* hCsIEForwardCal[10];
  TH1F* hCsITForward[10];
  TH2F* hCsIETForward[10];

  // Left CsI Detectors
  TH1F* hCsIELeft[6];
  TH1F* hCsIELeftCal[6];
  TH1F* hCsITLeft[6];
  TH2F* hCsIETLeft[6];

  // Si vs CsI
  TH2F* hSiCsIEForwardDet[10];
  TH2F* hSiCsIEForwardDetCal[10];
  TH2F* hSiCsIEForward[10][4];
  TH2F* hSiCsIEForwardCal[10][4];
  TH2F* hSiCsIELeft[6];
  TH2F* hSiCsIELeftCal[6];

  // dE vs E Forward Detectors
  TH2F* hdEEForward[10];
  TH2F* hdEEForwardCal[10];
  TH2F* hdEEForwardCalTotal[10];

  // dE vs E Left Detectors
  TH2F* hdEELeft[6];
  TH2F* hdEELeftCal[6];
  TH2F* hdEELeftCalTotal[6];

  // Angle vs E Forward Detectors
  TH2F* hAngleEForward[10];
  TH2F* hAngleEForwardCal[10];
  TH2F* hAngleEForwardCalTotal[10];

  // Angle vs E Left Detectors
  TH2F* hAngleELeft[10];
  TH2F* hAngleELeftCal[10];
  TH2F* hAngleELeftCalTotal[10];

  // Vertex vs E Forward Detectors
  TH2F* hVertexSiEForward[10];
  TH2F* hVertexSiEForwardCal[10];
  TH2F* hVertexSiEForwardCalTotal[10];

  // Vertex vs E Left Detectors
  TH2F* hVertexSiELeft[6];
  TH2F* hVertexSiELeftCal[6];
  TH2F* hVertexSiELeftCalTotal[6];

  // Vertex vs Angle Forward Detectors
  TH2F* hVertexAngleForward[10];

  // Cross Section
  TH1F* s1;

// Silicon Energy Calibration
private:
  void InitSiEForwardCalibration();
  std::pair<Float_t, Float_t> siEForwardCalibration[10][4] = {std::make_pair(0., 0.)};
  std::pair<Float_t, Float_t> siELeftCalibration[6][4] = {std::make_pair(0., 0.)};

  void InitCsIECalibration();
  std::pair<Float_t, Float_t> csiEForwardCalibration[10] = {std::make_pair(0., 0.)};

// Center Pad Gain Match
private:
  void InitCentralPadGainMatch();
  Double_t scale[6][128];

// Beam Average Central Pads
  void InitAverageBeamEnergy();
  Double_t averageBeamEnergy[128];

// Strip and Chain Matching
  void StripChainMatch(std::vector<mmTrack> &stripChainMatched, std::vector<mmTrack> &stripChainRaw, std::vector<mmStripChain> chain_,
                       std::vector<mmStripChain> strip_, Bool_t leftSide, Double_t siTime);
  size_t StripChainTime0TimeBuckets(std::vector<mmTrack> matched);
  size_t StripChainNumberTimeBuckets(std::vector<mmStripChain> chain, std::vector<mmStripChain> strip);
  void StripChainMatchingOutward(std::vector<mmTrack> &stripChainMatched, std::vector<mmStripChain> chain,
                                 std::vector<mmStripChain> strip, Bool_t leftSide, Double_t siTime);
  void StripChainMatchingBox(std::vector<mmTrack> &stripChainMatched, std::vector<mmStripChain> chain,
                             std::vector<mmStripChain> strip, Bool_t leftSide, Double_t siTime);
  void StripChainMatchingBoxTime0(std::vector<mmTrack> &stripChainMatched, std::vector<mmTrack> time0);
  void StripChainMatchingTime(std::vector<mmTrack> &stripChainMatched, std::vector<mmStripChain> chain,
                              std::vector<mmStripChain> strip, Bool_t leftSide, Double_t siTime, Int_t timeWindow);
  void StripChainMatchingTimeSlopeFit(std::vector<mmTrack> &stripChainMatched, std::vector<mmStripChain> chain,
                                   std::vector<mmStripChain> strip, Bool_t leftSide, Double_t siTime,
                                   Double_t timeWindow);
  void StripChainMatchingTimeSlopeHough(std::vector<mmTrack> &stripChainMatched, std::vector<mmStripChain> chain,
                                   std::vector<mmStripChain> strip, Bool_t leftSide, Double_t siTime,
                                   Double_t timeWindow);

// Cross Section
private:
  void DivideTargetThickness(TH1F *f);
  void ReadSolidAngle();
  void SolidAngle(TH1F *f);
  void WriteSpectrumToFile(TH1F *f, Int_t region);
  CubicSpline reg3SA;
  CubicSpline reg3CMAngle;

// General Variables
private:
  void InitVariables();
  TFile *file;
  Float_t m1;
  Float_t m2;
  Float_t rowConversion;
  Float_t rowConversionOffset;
  Float_t heightOffset;
  Float_t driftVelocity;
  Float_t timeResolution;
  Long64_t entry;
  Float_t beamEnergy;
  Float_t density;
  Float_t distanceHavarToSilicon;
  Float_t numberB8;
  Float_t siXPosForward[10][4];

  EnergyLoss *boronMethane;
  EnergyLoss *protonMethane;

// Tree Variables
  void InitTree();
  void FillTree();
  void WriteTree();
  TTree *outTree;
  Int_t siDet;
  Int_t siQuad;
  Int_t siChannel;
  Float_t siEnergy;
  Float_t siEnergyCal;
  Float_t siTime;
  Float_t csiEnergy;
  Float_t csiEnergyCal;
  Float_t csiTime;
  Float_t totalEnergy;
  Bool_t punchthrough;
  Float_t dE;
  Float_t vertexPositionX;
  Float_t vertexPositionY;
  Float_t vertexPositionZ;
  Float_t angle;
  Float_t cmEnergy;

};
#endif

#ifdef Spectra_cxx

inline Spectra::Spectra(TTree *tree) : fChain(0) {
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if(tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("run1_0.root");
    if(!f || !f->IsOpen()) {
      f = new TFile("run1_0.root");
    }
    f->GetObject("mfmData",tree);
  }
  Init(tree);
}

inline Spectra::~Spectra() {
  if(!fChain) return;
  delete fChain->GetCurrentFile();
}

inline Int_t Spectra::GetEntry(Long64_t entry) {
  // Read contents of entry.
  if(!fChain) return 0;
  return fChain->GetEntry(entry);
}

inline Long64_t Spectra::LoadTree(Long64_t entry) {
// Set the environment to read one entry
  if(!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if(centry < 0) return centry;
  if(fChain->GetTreeNumber() != fCurrent) {
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
  if(!tree) return;
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
  if(!fChain) return;
  fChain->Show(entry);
}

inline void Spectra::InitChannelMap() {
  printf("Initializing Channel Map\n");

  //////////////////////////
  // Forward Si Detectors //
  //////////////////////////

  // Detector 1
  siForwardMap[5]  = std::make_pair(0, 0);
  siForwardMap[12] = std::make_pair(0, 1);
  siForwardMap[18] = std::make_pair(0, 2);
  siForwardMap[25] = std::make_pair(0, 3);

  // Detector 2
  siForwardMap[31] = std::make_pair(1, 0);
  siForwardMap[37] = std::make_pair(1, 1);
  siForwardMap[43] = std::make_pair(1, 2);
  siForwardMap[50] = std::make_pair(1, 3);

  // Detector 3
  siForwardMap[54] = std::make_pair(2, 0);
  siForwardMap[59] = std::make_pair(2, 1);
  siForwardMap[63] = std::make_pair(2, 2);
  siForwardMap[67] = std::make_pair(2, 3);

  // Detector 4
  siForwardMap[23] = std::make_pair(3, 0);
  siForwardMap[3]  = std::make_pair(3, 1);
  siForwardMap[9]  = std::make_pair(3, 2);
  siForwardMap[16] = std::make_pair(3, 3);

  // Detector 5
  siForwardMap[29] = std::make_pair(4, 0);
  siForwardMap[35] = std::make_pair(4, 1);
  siForwardMap[41] = std::make_pair(4, 2);
  siForwardMap[48] = std::make_pair(4, 3);

  // Detector 6
  siForwardMap[52] = std::make_pair(5, 0);
  siForwardMap[57] = std::make_pair(5, 1);
  siForwardMap[61] = std::make_pair(5, 2);
  siForwardMap[65] = std::make_pair(5, 3);

  // Detector 7
  siForwardMap[1]  = std::make_pair(6, 0);
  siForwardMap[7]  = std::make_pair(6, 1);
  siForwardMap[14] = std::make_pair(6, 2);
  siForwardMap[20] = std::make_pair(6, 3);

  // Detector 8
  siForwardMap[46] = std::make_pair(7, 0);
  siForwardMap[27] = std::make_pair(7, 1);
  siForwardMap[33] = std::make_pair(7, 2);
  siForwardMap[39] = std::make_pair(7, 3);

  // Detector 9
  siForwardMap[4]  = std::make_pair(8, 0);
  siForwardMap[10] = std::make_pair(8, 1);
  siForwardMap[17] = std::make_pair(8, 2);
  siForwardMap[24] = std::make_pair(8, 3);

  // Detector 10
  siForwardMap[30] = std::make_pair(9, 0);
  siForwardMap[36] = std::make_pair(9, 1);
  siForwardMap[42] = std::make_pair(9, 2);
  siForwardMap[49] = std::make_pair(9, 3);

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
  siLeftMap[5]  = std::make_pair(0, 0);
  siLeftMap[12] = std::make_pair(0, 1);
  siLeftMap[18] = std::make_pair(0, 2);
  siLeftMap[25] = std::make_pair(0, 3);

  // Detector 2
  siLeftMap[50] = std::make_pair(1, 0);
  siLeftMap[31] = std::make_pair(1, 1);
  siLeftMap[37] = std::make_pair(1, 2);
  siLeftMap[43] = std::make_pair(1, 3);

  // Detector 3
  siLeftMap[54] = std::make_pair(2, 0);
  siLeftMap[59] = std::make_pair(2, 1);
  siLeftMap[63] = std::make_pair(2, 2);
  siLeftMap[67] = std::make_pair(2, 3);

  // Detector 4
  siLeftMap[23] = std::make_pair(3, 0);
  siLeftMap[3]  = std::make_pair(3, 1);
  siLeftMap[9]  = std::make_pair(3, 2);
  siLeftMap[16] = std::make_pair(3, 3);

  // Detector 5
  siLeftMap[29] = std::make_pair(4, 0);
  siLeftMap[35] = std::make_pair(4, 1);
  siLeftMap[41] = std::make_pair(4, 2);
  siLeftMap[48] = std::make_pair(4, 3);

  // Detector 6
  siLeftMap[65] = std::make_pair(5, 0);
  siLeftMap[52] = std::make_pair(5, 1);
  siLeftMap[57] = std::make_pair(5, 2);
  siLeftMap[61] = std::make_pair(5, 3);

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

  csiForwardMap[2]  = 0;
  csiForwardMap[7]  = 1;
  csiForwardMap[10] = 2;
  csiForwardMap[16] = 3;
  csiForwardMap[19] = 4;
  csiForwardMap[25] = 5;
  csiForwardMap[28] = 6;
  csiForwardMap[33] = 7;
  csiForwardMap[36] = 8;
  csiForwardMap[41] = 9;

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

  csiLeftMap[44] = 0;
  csiLeftMap[50] = 1;
  csiLeftMap[53] = 2;
  csiLeftMap[59] = 3;
  csiLeftMap[62] = 4;
  csiLeftMap[67] = 5;

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
  Int_t j = 0;
  for(Int_t i = 0; i < 68; i++) {
    if(i == 11 || i == 22 || i == 45 || i == 56) continue; // FPN Channels
    Aget_Map[j] = i;
    j++;
  }

  // Asad0_Aget0
  j = 31;
  for(Int_t i = 0; i < 64; i++) {
    if(i % 4 == 0 && i != 0) j--;
    Int_t remain = i % 4;
    MM_Map_Asad0_Aget0[Aget_Map[i]] = std::make_pair(4 - remain, j);
  }

  // Asad0_Aget1
  j = 63;
  for(Int_t i = 0; i < 64; i++) {
    if(i % 4 == 0 && i != 0) j--;
    Int_t remain = i % 4;
    MM_Map_Asad0_Aget1[Aget_Map[i]] = std::make_pair(4 - remain, j);
  }

  // Asad0_Aget2
  j = 95;
  for(Int_t i = 0; i < 64; i++) {
    if(i % 4 == 0 && i != 0) j--;
    Int_t remain = i % 4;
    MM_Map_Asad0_Aget2[Aget_Map[i]] = std::make_pair(1 + remain, j);
  }

  // Asad0_Aget3
  j = 127;
  for(Int_t i = 0; i < 64; i++) {
    if(i % 4 == 0 && i != 0) j--;
    Int_t remain = i % 4;
    MM_Map_Asad0_Aget3[Aget_Map[i]] = std::make_pair(1 + remain, j);
  }

  // Asad1_Aget0
  j = 15;
  for(Int_t i = 0; i < 64; i++) {
    if(i % 4 == 0 && i != 0) j--;
    Int_t remain = i % 4;
    MM_Map_Asad1_Aget0[Aget_Map[i]] = std::make_pair(1 + remain, j);
  }

  // Asad1_Aget1
  j = 47;
  for(Int_t i = 0; i < 64; i++) {
    if(i % 4 == 0 && i != 0) j--;
    Int_t remain = i % 4;
    MM_Map_Asad1_Aget1[Aget_Map[i]] = std::make_pair(1 + remain, j);
  }

  // Asad1_Aget2
  j = 79;
  for(Int_t i = 0; i < 64; i++) {
    if(i % 4 == 0 && i != 0) j--;
    Int_t remain = i % 4;
    MM_Map_Asad1_Aget2[Aget_Map[i]] = std::make_pair(4 - remain, j);
  }

  // Asad1_Aget3
  j = 111;
  for(Int_t i = 0; i < 64; i++) {
    if(i % 4 == 0 && i != 0) j--;
    Int_t remain = i % 4;
    MM_Map_Asad1_Aget3[Aget_Map[i]] = std::make_pair(4 - remain, j);
  }

  // Asad2_Aget0 (Beam left)
  // Strips
  MM_Map_Asad2_Aget0[Aget_Map[28]] = 0;
  MM_Map_Asad2_Aget0[Aget_Map[29]] = 1;
  MM_Map_Asad2_Aget0[Aget_Map[30]] = 2;
  MM_Map_Asad2_Aget0[Aget_Map[31]] = 3;
  MM_Map_Asad2_Aget0[Aget_Map[24]] = 4;
  MM_Map_Asad2_Aget0[Aget_Map[25]] = 5;
  MM_Map_Asad2_Aget0[Aget_Map[26]] = 6;
  MM_Map_Asad2_Aget0[Aget_Map[27]] = 7;
  MM_Map_Asad2_Aget0[Aget_Map[20]] = 8;
  MM_Map_Asad2_Aget0[Aget_Map[21]] = 9;
  MM_Map_Asad2_Aget0[Aget_Map[22]] = 10;
  MM_Map_Asad2_Aget0[Aget_Map[23]] = 11;
  MM_Map_Asad2_Aget0[Aget_Map[16]] = 12;
  MM_Map_Asad2_Aget0[Aget_Map[17]] = 13;
  MM_Map_Asad2_Aget0[Aget_Map[18]] = 14;
  MM_Map_Asad2_Aget0[Aget_Map[19]] = 15;
  MM_Map_Asad2_Aget0[Aget_Map[12]] = 16;
  MM_Map_Asad2_Aget0[Aget_Map[13]] = 17;
  MM_Map_Asad2_Aget0[Aget_Map[14]] = 18;
  MM_Map_Asad2_Aget0[Aget_Map[15]] = 19;
  MM_Map_Asad2_Aget0[Aget_Map[8]] = 20;
  MM_Map_Asad2_Aget0[Aget_Map[9]] = 21;
  MM_Map_Asad2_Aget0[Aget_Map[10]] = 22;
  MM_Map_Asad2_Aget0[Aget_Map[11]] = 23;
  MM_Map_Asad2_Aget0[Aget_Map[4]] = 24;
  MM_Map_Asad2_Aget0[Aget_Map[5]] = 25;
  MM_Map_Asad2_Aget0[Aget_Map[6]] = 26;
  MM_Map_Asad2_Aget0[Aget_Map[7]] = 27;
  MM_Map_Asad2_Aget0[Aget_Map[0]] = 28;
  MM_Map_Asad2_Aget0[Aget_Map[1]] = 29;
  MM_Map_Asad2_Aget0[Aget_Map[2]] = 30;
  MM_Map_Asad2_Aget0[Aget_Map[3]] = 31;
  // j = 31;
  // for(Int_t i = 0; i < 32; i++) {
  //   MM_Map_Asad2_Aget0[Aget_Map[i]] = j;
  //   j--;
  // }
  // Chains
  MM_Map_Asad2_Aget0[Aget_Map[35]] = 32;
  MM_Map_Asad2_Aget0[Aget_Map[34]] = 33;
  MM_Map_Asad2_Aget0[Aget_Map[33]] = 34;
  MM_Map_Asad2_Aget0[Aget_Map[32]] = 35;
  MM_Map_Asad2_Aget0[Aget_Map[39]] = 36;
  MM_Map_Asad2_Aget0[Aget_Map[38]] = 37;
  MM_Map_Asad2_Aget0[Aget_Map[37]] = 38;
  MM_Map_Asad2_Aget0[Aget_Map[36]] = 39;
  MM_Map_Asad2_Aget0[Aget_Map[43]] = 40;
  MM_Map_Asad2_Aget0[Aget_Map[42]] = 41;
  MM_Map_Asad2_Aget0[Aget_Map[41]] = 42;
  MM_Map_Asad2_Aget0[Aget_Map[40]] = 43;
  MM_Map_Asad2_Aget0[Aget_Map[47]] = 44;
  MM_Map_Asad2_Aget0[Aget_Map[46]] = 45;
  MM_Map_Asad2_Aget0[Aget_Map[45]] = 46;
  MM_Map_Asad2_Aget0[Aget_Map[44]] = 47;
  MM_Map_Asad2_Aget0[Aget_Map[51]] = 48;
  MM_Map_Asad2_Aget0[Aget_Map[50]] = 49;
  MM_Map_Asad2_Aget0[Aget_Map[49]] = 50;
  MM_Map_Asad2_Aget0[Aget_Map[48]] = 51;
  MM_Map_Asad2_Aget0[Aget_Map[55]] = 52;
  MM_Map_Asad2_Aget0[Aget_Map[54]] = 53;
  MM_Map_Asad2_Aget0[Aget_Map[53]] = 54;
  MM_Map_Asad2_Aget0[Aget_Map[52]] = 55;
  MM_Map_Asad2_Aget0[Aget_Map[59]] = 56;
  MM_Map_Asad2_Aget0[Aget_Map[58]] = 57;
  MM_Map_Asad2_Aget0[Aget_Map[57]] = 58;
  MM_Map_Asad2_Aget0[Aget_Map[56]] = 59;
  MM_Map_Asad2_Aget0[Aget_Map[63]] = 60;
  MM_Map_Asad2_Aget0[Aget_Map[62]] = 61;
  MM_Map_Asad2_Aget0[Aget_Map[61]] = 62;
  MM_Map_Asad2_Aget0[Aget_Map[60]] = 63;
  // j = 32;
  // for(Int_t i = 32; i < 64; i++) {
  //   MM_Map_Asad2_Aget0[Aget_Map[i]] = j;
  //   j++;
  // }

  // Asad2_Aget1 (Beam left)
  // Strips
  MM_Map_Asad2_Aget1[Aget_Map[28]] = 32;
  MM_Map_Asad2_Aget1[Aget_Map[29]] = 33;
  MM_Map_Asad2_Aget1[Aget_Map[30]] = 34;
  MM_Map_Asad2_Aget1[Aget_Map[31]] = 35;
  MM_Map_Asad2_Aget1[Aget_Map[24]] = 36;
  MM_Map_Asad2_Aget1[Aget_Map[25]] = 37;
  MM_Map_Asad2_Aget1[Aget_Map[26]] = 38;
  MM_Map_Asad2_Aget1[Aget_Map[27]] = 39;
  MM_Map_Asad2_Aget1[Aget_Map[20]] = 40;
  MM_Map_Asad2_Aget1[Aget_Map[21]] = 41;
  MM_Map_Asad2_Aget1[Aget_Map[22]] = 42;
  MM_Map_Asad2_Aget1[Aget_Map[23]] = 43;
  MM_Map_Asad2_Aget1[Aget_Map[16]] = 44;
  MM_Map_Asad2_Aget1[Aget_Map[17]] = 45;
  MM_Map_Asad2_Aget1[Aget_Map[18]] = 46;
  MM_Map_Asad2_Aget1[Aget_Map[19]] = 47;
  MM_Map_Asad2_Aget1[Aget_Map[12]] = 48;
  MM_Map_Asad2_Aget1[Aget_Map[13]] = 49;
  MM_Map_Asad2_Aget1[Aget_Map[14]] = 50;
  MM_Map_Asad2_Aget1[Aget_Map[15]] = 51;
  MM_Map_Asad2_Aget1[Aget_Map[8]] = 52;
  MM_Map_Asad2_Aget1[Aget_Map[9]] = 53;
  MM_Map_Asad2_Aget1[Aget_Map[10]] = 54;
  MM_Map_Asad2_Aget1[Aget_Map[11]] = 55;
  MM_Map_Asad2_Aget1[Aget_Map[4]] = 56;
  MM_Map_Asad2_Aget1[Aget_Map[5]] = 57;
  MM_Map_Asad2_Aget1[Aget_Map[6]] = 58;
  MM_Map_Asad2_Aget1[Aget_Map[7]] = 59;
  MM_Map_Asad2_Aget1[Aget_Map[0]] = 60;
  MM_Map_Asad2_Aget1[Aget_Map[1]] = 61;
  MM_Map_Asad2_Aget1[Aget_Map[2]] = 62;
  MM_Map_Asad2_Aget1[Aget_Map[3]] = 63;
  // j = 63;
  // for(Int_t i = 0; i < 32; i++) {
  //   MM_Map_Asad2_Aget1[Aget_Map[i]] = j;
  //   j--;
  // }
  // Chains
  MM_Map_Asad2_Aget1[Aget_Map[35]] = 0;
  MM_Map_Asad2_Aget1[Aget_Map[34]] = 1;
  MM_Map_Asad2_Aget1[Aget_Map[33]] = 2;
  MM_Map_Asad2_Aget1[Aget_Map[32]] = 3;
  MM_Map_Asad2_Aget1[Aget_Map[39]] = 4;
  MM_Map_Asad2_Aget1[Aget_Map[38]] = 5;
  MM_Map_Asad2_Aget1[Aget_Map[37]] = 6;
  MM_Map_Asad2_Aget1[Aget_Map[36]] = 7;
  MM_Map_Asad2_Aget1[Aget_Map[43]] = 8;
  MM_Map_Asad2_Aget1[Aget_Map[42]] = 9;
  MM_Map_Asad2_Aget1[Aget_Map[41]] = 10;
  MM_Map_Asad2_Aget1[Aget_Map[40]] = 11;
  MM_Map_Asad2_Aget1[Aget_Map[47]] = 12;
  MM_Map_Asad2_Aget1[Aget_Map[46]] = 13;
  MM_Map_Asad2_Aget1[Aget_Map[45]] = 14;
  MM_Map_Asad2_Aget1[Aget_Map[44]] = 15;
  MM_Map_Asad2_Aget1[Aget_Map[51]] = 16;
  MM_Map_Asad2_Aget1[Aget_Map[50]] = 17;
  MM_Map_Asad2_Aget1[Aget_Map[49]] = 18;
  MM_Map_Asad2_Aget1[Aget_Map[48]] = 19;
  MM_Map_Asad2_Aget1[Aget_Map[55]] = 20;
  MM_Map_Asad2_Aget1[Aget_Map[54]] = 21;
  MM_Map_Asad2_Aget1[Aget_Map[53]] = 22;
  MM_Map_Asad2_Aget1[Aget_Map[52]] = 23;
  MM_Map_Asad2_Aget1[Aget_Map[59]] = 24;
  MM_Map_Asad2_Aget1[Aget_Map[58]] = 25;
  MM_Map_Asad2_Aget1[Aget_Map[57]] = 26;
  MM_Map_Asad2_Aget1[Aget_Map[56]] = 27;
  MM_Map_Asad2_Aget1[Aget_Map[63]] = 28;
  MM_Map_Asad2_Aget1[Aget_Map[62]] = 29;
  MM_Map_Asad2_Aget1[Aget_Map[61]] = 30;
  MM_Map_Asad2_Aget1[Aget_Map[60]] = 31;
  // j = 0;
  // for(Int_t i = 32; i < 64; i++) {
  //   MM_Map_Asad2_Aget1[Aget_Map[i]] = j;
  //   j++;
  // }

  // Asad2_Aget2
  j = 63;
  for(Int_t i = 0; i < 64; i++) {
    MM_Map_Asad2_Aget2[Aget_Map[i]] = std::make_pair(0, j);
    j--;
  }

  // Asad2_Aget3
  j = 127;
  for(Int_t i = 0; i < 64; i++) {
    MM_Map_Asad2_Aget3[Aget_Map[i]] = std::make_pair(0, j);
    j--;
  }

  // Asad3_Aget0 (Beam right)
  // Strips
  j = 31;
  for(Int_t i = 0; i < 32; i++) {
    MM_Map_Asad3_Aget0[Aget_Map[i]] = j;
    j--;
  }
  // Chains
  j = 32;
  for(Int_t i = 32; i < 64; i++) {
    MM_Map_Asad3_Aget0[Aget_Map[i]] = j;
    j++;
  }

  // Asad3_Aget1 (Beam right)
  // Strips
  j = 63;
  for(Int_t i = 0; i < 32; i++) {
    MM_Map_Asad3_Aget1[Aget_Map[i]] = j;
    j--;
  }
  // Chains
  j = 0;
  for(Int_t i = 32; i < 64; i++) {
    MM_Map_Asad3_Aget1[Aget_Map[i]] = j;
    j++;
  }

  // Asad3_Aget2
  j = 63;
  for(Int_t i = 0; i < 64; i++) {
    MM_Map_Asad3_Aget2[Aget_Map[i]] = std::make_pair(5, j);
    j--;
  }

  // Asad3_Aget3
  j = 127;
  for(Int_t i = 0; i < 64; i++) {
    MM_Map_Asad3_Aget3[Aget_Map[i]] = std::make_pair(5, j);
    j--;
  }

}

inline void Spectra::InitHistograms() {
  printf("Initializing Histograms\n");

  // Ionization Chamber
  hIonizationChamberE = new TH1F("icE", "icE", 500, 0, 4000);
  hIonizationChamberE->GetXaxis()->SetTitle("Energy [channels]"); hIonizationChamberE->GetXaxis()->CenterTitle();
  hIonizationChamberE->GetYaxis()->SetTitle("Counts"); hIonizationChamberE->GetYaxis()->CenterTitle();
  hIonizationChamberE->GetYaxis()->SetTitleOffset(1.4);

  hIonizationChamberT = new TH1F("icT", "icT", 500, 0, 20000);
  hIonizationChamberT->GetXaxis()->SetTitle("Time [ns]"); hIonizationChamberT->GetXaxis()->CenterTitle();
  hIonizationChamberT->GetYaxis()->SetTitle("Counts"); hIonizationChamberT->GetYaxis()->CenterTitle();
  hIonizationChamberT->GetYaxis()->SetTitleOffset(1.4);

  // Micromegas - Central
  hMicroMegasCenterCumulative = new TH2F("MM_Center_Cumulative", "MM_Center_Cumulative", 82, -41, 41, 128, 0, 128);
  hMicroMegasCenterCumulativePosition = new TH2F("MM_Center_Cumulative_Position", "MM_Center_Cumulative_Position", 82, -41, 41, 128, 0, 250);
  hMicroMegasCenterEnergyCumulative = new TH2F("MM_Center_Energy_Cumulative", "MM_Center_Energy_Cumulative", 20, -10, 10, 128, 0, 128);
  hMicroMegasCenterEnergyAverage = new TH2F("MM_Center_Energy_Average", "MM_Center_Energy_Average", 20, -10, 10, 128, 0, 128);
  hMicroMegasCenterEnergyAverageScaled = new TH2F("MM_Center_Energy_Average_Scaled", "MM_Center_Energy_Average_Scaled", 20, -10, 10, 128, 0, 128);
  hMicroMegasCenterTimeAverage = new TH2F("MM_Center_Time_Average", "MM_Center_Time_Average", 20, -10, 10, 128, 0, 128);

  // Micromegas - Left
  hMicroMegasStripLeftCumulative = new TH1F("MM_Strip_Left_Cumulative", "MM_Strip_Left_Cumulative", 64, 0, 64);
  hMicroMegasChainLeftCumulative = new TH1F("MM_Chain_Left_Cumulative", "MM_Chain_Left_Cumulative", 64, 0, 64);

  // Micromegas - Right
  hMicroMegasStripRightCumulative = new TH1F("MM_Strip_Right_Cumulative", "MM_Strip_Right_Cumulative", 64, 0, 64);
  hMicroMegasChainRightCumulative = new TH1F("MM_Chain_Right_Cumulative", "MM_Chain_Right_Cumulative", 64, 0, 64);

  // Histograms for the Forward Si Detectors
  for(UInt_t i = 0; i < 10; i++) {
    TString name = Form("siEForward_d%d", i);
    hSiEForwardDet[i] = new TH1F(name, name, 500, 0, 4000);
    hSiEForwardDet[i]->GetXaxis()->SetTitle("Energy [channels]"); hSiEForwardDet[i]->GetXaxis()->CenterTitle();
    hSiEForwardDet[i]->GetYaxis()->SetTitle("Counts"); hSiEForwardDet[i]->GetYaxis()->CenterTitle();
    hSiEForwardDet[i]->GetYaxis()->SetTitleOffset(1.4);

    name = Form("siEForwardCal_d%d", i);
    hSiEForwardDetCal[i] = new TH1F(name, name, 2000, 0, 16000);
    hSiEForwardDetCal[i]->GetXaxis()->SetTitle("Energy [keV]"); hSiEForwardDetCal[i]->GetXaxis()->CenterTitle();
    hSiEForwardDetCal[i]->GetYaxis()->SetTitle("Counts"); hSiEForwardDetCal[i]->GetYaxis()->CenterTitle();
    hSiEForwardDetCal[i]->GetYaxis()->SetTitleOffset(1.4);

    name = Form("siTForward_d%d", i);
    hSiTForwardDet[i] = new TH1F(name, name, 2500, 0, 20000);
    hSiTForwardDet[i]->GetXaxis()->SetTitle("Time [ns]"); hSiTForwardDet[i]->GetXaxis()->CenterTitle();
    hSiTForwardDet[i]->GetYaxis()->SetTitle("Counts"); hSiTForwardDet[i]->GetYaxis()->CenterTitle();
    hSiTForwardDet[i]->GetYaxis()->SetTitleOffset(1.4);

    for(UInt_t j = 0; j < 4; j++) {
      TString name = Form("siEForward_d%d_q%d", i, j);
      hSiEForward[i][j] = new TH1F(name, name, 1000, 0, 1000);
      hSiEForward[i][j]->GetXaxis()->SetTitle("Energy [channels]"); hSiEForward[i][j]->GetXaxis()->CenterTitle();
      hSiEForward[i][j]->GetYaxis()->SetTitle("Counts"); hSiEForward[i][j]->GetYaxis()->CenterTitle();
      hSiEForward[i][j]->GetYaxis()->SetTitleOffset(1.4);

      name = Form("siEForwardCal_d%d_q%d", i, j);
      hSiEForwardCal[i][j] = new TH1F(name, name, 500, 0, 4000);
      hSiEForwardCal[i][j]->GetXaxis()->SetTitle("Energy [keV]"); hSiEForwardCal[i][j]->GetXaxis()->CenterTitle();
      hSiEForwardCal[i][j]->GetYaxis()->SetTitle("Counts"); hSiEForwardCal[i][j]->GetYaxis()->CenterTitle();
      hSiEForwardCal[i][j]->GetYaxis()->SetTitleOffset(1.4);

      name = Form("siTForward_d%d_q%d", i, j);
      hSiTForward[i][j] = new TH1F(name, name, 500, 0, 4000);
      hSiTForward[i][j]->GetXaxis()->SetTitle("Time [ns]"); hSiTForward[i][j]->GetXaxis()->CenterTitle();
      hSiTForward[i][j]->GetYaxis()->SetTitle("Counts"); hSiTForward[i][j]->GetYaxis()->CenterTitle();
      hSiTForward[i][j]->GetYaxis()->SetTitleOffset(1.4);
    }
  }

  // Histograms for the Beam Left Si Detectors
  for(UInt_t i = 0; i < 6; i++) {
    if(i > 7) break;
    TString name = Form("siELeft_d%d", i);
    hSiELeftDet[i] = new TH1F(name, name, 500, 0, 4000);
    hSiELeftDet[i]->GetXaxis()->SetTitle("Energy [channels]"); hSiELeftDet[i]->GetXaxis()->CenterTitle();
    hSiELeftDet[i]->GetYaxis()->SetTitle("Counts"); hSiELeftDet[i]->GetYaxis()->CenterTitle();
    hSiELeftDet[i]->GetYaxis()->SetTitleOffset(1.4);

    name = Form("siELeftCal_d%d", i);
    hSiELeftDetCal[i] = new TH1F(name, name, 2000, 0, 16000);
    hSiELeftDetCal[i]->GetXaxis()->SetTitle("Energy [keV]"); hSiELeftDetCal[i]->GetXaxis()->CenterTitle();
    hSiELeftDetCal[i]->GetYaxis()->SetTitle("Counts"); hSiELeftDetCal[i]->GetYaxis()->CenterTitle();
    hSiELeftDetCal[i]->GetYaxis()->SetTitleOffset(1.4);

    name = Form("siTLeft_d%d", i);
    hSiTLeftDet[i] = new TH1F(name, name, 2500, 0, 20000);
    hSiTLeftDet[i]->GetXaxis()->SetTitle("Time [ns]"); hSiTLeftDet[i]->GetXaxis()->CenterTitle();
    hSiTLeftDet[i]->GetYaxis()->SetTitle("Counts"); hSiTLeftDet[i]->GetYaxis()->CenterTitle();
    hSiTLeftDet[i]->GetYaxis()->SetTitleOffset(1.4);

    for(UInt_t j = 0; j < 4; j++) {
      TString name = Form("siELeft_d%d_q%d", i, j);
      hSiELeft[i][j] = new TH1F(name, name, 1000, 0, 1000);
      hSiELeft[i][j]->GetXaxis()->SetTitle("Energy [channels]"); hSiELeft[i][j]->GetXaxis()->CenterTitle();
      hSiELeft[i][j]->GetYaxis()->SetTitle("Counts"); hSiELeft[i][j]->GetYaxis()->CenterTitle();
      hSiELeft[i][j]->GetYaxis()->SetTitleOffset(1.4);

      name = Form("siELeftCal_d%d_q%d", i, j);
      hSiELeftCal[i][j] = new TH1F(name, name, 500, 0, 4000);
      hSiELeftCal[i][j]->GetXaxis()->SetTitle("Energy [keV]"); hSiELeftCal[i][j]->GetXaxis()->CenterTitle();
      hSiELeftCal[i][j]->GetYaxis()->SetTitle("Counts"); hSiELeftCal[i][j]->GetYaxis()->CenterTitle();
      hSiELeftCal[i][j]->GetYaxis()->SetTitleOffset(1.4);

      name = Form("siTLeft_d%d_q%d", i, j);
      hSiTLeft[i][j] = new TH1F(name, name, 500, 0, 4000);
      hSiTLeft[i][j]->GetXaxis()->SetTitle("Time [ns]"); hSiTLeft[i][j]->GetXaxis()->CenterTitle();
      hSiTLeft[i][j]->GetYaxis()->SetTitle("Counts"); hSiTLeft[i][j]->GetYaxis()->CenterTitle();
      hSiTLeft[i][j]->GetYaxis()->SetTitleOffset(1.4);
    }
  }

  // Histograms for the Forward CsI Detectors
  for(UInt_t i = 0; i < 10; i++) {
    TString name = Form("CsIEForward_d%d", i);
    hCsIEForward[i] = new TH1F(name, name, 500, 0, 4000);
    hCsIEForward[i]->GetXaxis()->SetTitle("Energy [channels]"); hCsIEForward[i]->GetXaxis()->CenterTitle();
    hCsIEForward[i]->GetYaxis()->SetTitle("Counts"); hCsIEForward[i]->GetYaxis()->CenterTitle();
    hCsIEForward[i]->GetYaxis()->SetTitleOffset(1.4);

    name = Form("CsIEForwardCal_d%d", i);
    hCsIEForwardCal[i] = new TH1F(name, name, 500, 0, 4000);
    hCsIEForwardCal[i]->GetXaxis()->SetTitle("Energy [keV]"); hCsIEForwardCal[i]->GetXaxis()->CenterTitle();
    hCsIEForwardCal[i]->GetYaxis()->SetTitle("Counts"); hCsIEForwardCal[i]->GetYaxis()->CenterTitle();
    hCsIEForwardCal[i]->GetYaxis()->SetTitleOffset(1.4);

    name = Form("CsITForward_d%d", i);
    hCsITForward[i] = new TH1F(name, name, 500, 0, 20000);
    hCsITForward[i]->GetXaxis()->SetTitle("Time [ns]"); hCsITForward[i]->GetXaxis()->CenterTitle();
    hCsITForward[i]->GetYaxis()->SetTitle("Counts"); hCsITForward[i]->GetYaxis()->CenterTitle();
    hCsITForward[i]->GetYaxis()->SetTitleOffset(1.4);

    name = Form("CsIETForward_d%d", i);
    hCsIETForward[i] = new TH2F(name, name, 500, 0, 10000, 500, 0, 20000);
    hCsIETForward[i]->GetXaxis()->SetTitle("Energy [Channels]"); hCsIETForward[i]->GetXaxis()->CenterTitle();
    hCsIETForward[i]->GetYaxis()->SetTitle("Time [ns]"); hCsIETForward[i]->GetYaxis()->CenterTitle();
    hCsIETForward[i]->GetYaxis()->SetTitleOffset(1.4);
  }

  // Histograms for the Beam Left CsI Detectors
  for(UInt_t i = 0; i < 6; i++) {
    TString name = Form("CsIELeft_d%d", i);
    hCsIELeft[i] = new TH1F(name,name,1000,0,4000);
    hCsIELeft[i]->GetXaxis()->SetTitle("Energy [channels]"); hCsIELeft[i]->GetXaxis()->CenterTitle();
    hCsIELeft[i]->GetYaxis()->SetTitle("Counts"); hCsIELeft[i]->GetYaxis()->CenterTitle();
    hCsIELeft[i]->GetYaxis()->SetTitleOffset(1.4);

    name = Form("CsIELeftCal_d%d", i);
    hCsIELeftCal[i] = new TH1F(name,name,1000,0,4000);
    hCsIELeftCal[i]->GetXaxis()->SetTitle("Energy [keV]"); hCsIELeftCal[i]->GetXaxis()->CenterTitle();
    hCsIELeftCal[i]->GetYaxis()->SetTitle("Counts"); hCsIELeftCal[i]->GetYaxis()->CenterTitle();
    hCsIELeftCal[i]->GetYaxis()->SetTitleOffset(1.4);

    name = Form("CsITLeft_d%d", i);
    hCsITLeft[i] = new TH1F(name,name,1000,0,20000);
    hCsITLeft[i]->GetXaxis()->SetTitle("Time [ns]"); hCsITLeft[i]->GetXaxis()->CenterTitle();
    hCsITLeft[i]->GetYaxis()->SetTitle("Counts"); hCsITLeft[i]->GetYaxis()->CenterTitle();
    hCsITLeft[i]->GetYaxis()->SetTitleOffset(1.4);

    name = Form("CsIETLeft_d%d", i);
    hCsIETLeft[i] = new TH2F(name, name, 500, 0, 10000, 500, 0, 20000);
    hCsIETLeft[i]->GetXaxis()->SetTitle("Energy [Channels]"); hCsIETLeft[i]->GetXaxis()->CenterTitle();
    hCsIETLeft[i]->GetYaxis()->SetTitle("Time [ns]"); hCsIETLeft[i]->GetYaxis()->CenterTitle();
    hCsIETLeft[i]->GetYaxis()->SetTitleOffset(1.4);
  }

  // Si vs CsI Forward Detectors
  for(UInt_t i = 0; i < 10; i++) {
    TString name = Form("siCsIEForward_d%d", i);
    hSiCsIEForwardDet[i] = new TH2F(name, name, 500, 0, 4000, 500, 0, 10000);
    hSiCsIEForwardDet[i]->GetXaxis()->SetTitle("Si Energy [channels]"); hSiCsIEForwardDet[i]->GetXaxis()->CenterTitle();
    hSiCsIEForwardDet[i]->GetYaxis()->SetTitle("CsI Energy [channels]"); hSiCsIEForwardDet[i]->GetYaxis()->CenterTitle();
    hSiCsIEForwardDet[i]->GetYaxis()->SetTitleOffset(1.4);

    name = Form("siCsIEForwardCal_d%d", i);
    hSiCsIEForwardDetCal[i] = new TH2F(name, name, 500, 0, 14000, 500, 0, 10000);
    hSiCsIEForwardDetCal[i]->GetXaxis()->SetTitle("Si Energy [keV]"); hSiCsIEForwardDetCal[i]->GetXaxis()->CenterTitle();
    hSiCsIEForwardDetCal[i]->GetYaxis()->SetTitle("CsI Energy [channels]"); hSiCsIEForwardDetCal[i]->GetYaxis()->CenterTitle();
    hSiCsIEForwardDetCal[i]->GetYaxis()->SetTitleOffset(1.4);

    for(UInt_t j = 0; j < 4; j++) {
      name = Form("siCsIEForward_d%d_q%d", i, j);
      hSiCsIEForward[i][j] = new TH2F(name, name, 500, 0, 4000, 500, 0, 10000);
      hSiCsIEForward[i][j]->GetXaxis()->SetTitle("Si Energy [channels]"); hSiCsIEForward[i][j]->GetXaxis()->CenterTitle();
      hSiCsIEForward[i][j]->GetYaxis()->SetTitle("CsI Energy [channels]"); hSiCsIEForward[i][j]->GetYaxis()->CenterTitle();
      hSiCsIEForward[i][j]->GetYaxis()->SetTitleOffset(1.4);

      name = Form("siCsIEForwardCal_d%d_q%d", i, j);
      hSiCsIEForwardCal[i][j] = new TH2F(name, name, 500, 0, 14000, 500, 0, 10000);
      hSiCsIEForwardCal[i][j]->GetXaxis()->SetTitle("Si Energy [keV]"); hSiCsIEForwardCal[i][j]->GetXaxis()->CenterTitle();
      hSiCsIEForwardCal[i][j]->GetYaxis()->SetTitle("CsI Energy [channels]"); hSiCsIEForwardCal[i][j]->GetYaxis()->CenterTitle();
      hSiCsIEForwardCal[i][j]->GetYaxis()->SetTitleOffset(1.4);
    }
  }

  // dE vs E Histograms
  for(UInt_t i = 0; i < 10; i++) {
    TString name = Form("dEEForward_d%d", i);
    hdEEForward[i] = new TH2F(name, name, 500, 0, 4000, 500, 0, 4000);
    hdEEForward[i]->GetXaxis()->SetTitle("Si Energy [channels]"); hdEEForward[i]->GetXaxis()->CenterTitle();
    hdEEForward[i]->GetYaxis()->SetTitle("dE [channels]"); hdEEForward[i]->GetYaxis()->CenterTitle();
    hdEEForward[i]->GetYaxis()->SetTitleOffset(1.4);

    name = Form("dEEForwardCal_d%d", i);
    hdEEForwardCal[i] = new TH2F(name, name, 500, 0, 14000, 500, 0, 4000);
    hdEEForwardCal[i]->GetXaxis()->SetTitle("Si Energy [keV]"); hdEEForwardCal[i]->GetXaxis()->CenterTitle();
    hdEEForwardCal[i]->GetYaxis()->SetTitle("dE [channels]"); hdEEForwardCal[i]->GetYaxis()->CenterTitle();
    hdEEForwardCal[i]->GetYaxis()->SetTitleOffset(1.4);

    name = Form("dEEForwardCalTotal_d%d", i);
    hdEEForwardCalTotal[i] = new TH2F(name, name, 500, 0, 25000, 500, 0, 4000);
    hdEEForwardCalTotal[i]->GetXaxis()->SetTitle("Si Energy [keV]"); hdEEForwardCalTotal[i]->GetXaxis()->CenterTitle();
    hdEEForwardCalTotal[i]->GetYaxis()->SetTitle("dE [channels]"); hdEEForwardCalTotal[i]->GetYaxis()->CenterTitle();
    hdEEForwardCalTotal[i]->GetYaxis()->SetTitleOffset(1.4);
  }

  for(UInt_t i = 0; i < 6; i++) {
    TString name = Form("dEELeft_d%d", i);
    hdEELeft[i] = new TH2F(name, name, 500, 0, 4000, 500, 0, 4000);
    hdEELeft[i]->GetXaxis()->SetTitle("Si Energy [channels]"); hdEELeft[i]->GetXaxis()->CenterTitle();
    hdEELeft[i]->GetYaxis()->SetTitle("dE [channels]"); hdEELeft[i]->GetYaxis()->CenterTitle();
    hdEELeft[i]->GetYaxis()->SetTitleOffset(1.4);
  }

  // Vertex vs E Histograms
  for(UInt_t i = 0; i < 10; i++) {
    TString name = Form("vertexSiEForward_d%d", i);
    hVertexSiEForward[i] = new TH2F(name, name, 500, 0, 4000, 500, -250, 300);
    hVertexSiEForward[i]->GetXaxis()->SetTitle("Si Energy [channels]"); hVertexSiEForward[i]->GetXaxis()->CenterTitle();
    hVertexSiEForward[i]->GetYaxis()->SetTitle("Vertex [mm]"); hVertexSiEForward[i]->GetYaxis()->CenterTitle();
    hVertexSiEForward[i]->GetYaxis()->SetTitleOffset(1.4);

    name = Form("vertexSiEForwardCal_d%d", i);
    hVertexSiEForwardCal[i] = new TH2F(name, name, 500, 0, 14000, 500, -250, 300);
    hVertexSiEForwardCal[i]->GetXaxis()->SetTitle("Si Energy [keV]"); hVertexSiEForwardCal[i]->GetXaxis()->CenterTitle();
    hVertexSiEForwardCal[i]->GetYaxis()->SetTitle("Vertex [mm]"); hVertexSiEForwardCal[i]->GetYaxis()->CenterTitle();
    hVertexSiEForwardCal[i]->GetYaxis()->SetTitleOffset(1.4);
  }

  for(UInt_t i = 0; i < 6; i++) {
    TString name = Form("vertexSiELeft_d%d", i);
    hVertexSiELeft[i] = new TH2F(name, name, 500, 0, 4000, 500, -250, 300);
    hVertexSiELeft[i]->GetXaxis()->SetTitle("Si Energy [channels]"); hVertexSiELeft[i]->GetXaxis()->CenterTitle();
    hVertexSiELeft[i]->GetYaxis()->SetTitle("Vertex [mm]"); hVertexSiELeft[i]->GetYaxis()->CenterTitle();
    hVertexSiELeft[i]->GetYaxis()->SetTitleOffset(1.4);
  }

  // Angle vs E Histograms hAngleEForward[10]
  for(UInt_t i = 0; i < 10; i++) {
    TString name = Form("angleSiEForward_d%d", i);
    hAngleEForward[i] = new TH2F(name, name, 200, 0, 4000, 200, 0, 1.5);
    hAngleEForward[i]->GetXaxis()->SetTitle("Si Energy [channels]"); hAngleEForward[i]->GetXaxis()->CenterTitle();
    hAngleEForward[i]->GetYaxis()->SetTitle("Vertex [mm]"); hAngleEForward[i]->GetYaxis()->CenterTitle();
    hAngleEForward[i]->GetYaxis()->SetTitleOffset(1.4);

    name = Form("angleSiEForwardCal_d%d", i);
    hAngleEForwardCal[i] = new TH2F(name, name, 200, 0, 20000, 200, 0, 1.5);
    hAngleEForwardCal[i]->GetXaxis()->SetTitle("Si Energy [keV]"); hAngleEForwardCal[i]->GetXaxis()->CenterTitle();
    hAngleEForwardCal[i]->GetYaxis()->SetTitle("Vertex [mm]"); hAngleEForwardCal[i]->GetYaxis()->CenterTitle();
    hAngleEForwardCal[i]->GetYaxis()->SetTitleOffset(1.4);

    name = Form("angleSiEForwardCalTotal_d%d", i);
    hAngleEForwardCalTotal[i] = new TH2F(name, name, 200, 0, 20000, 200, 0, 1.5);
    hAngleEForwardCalTotal[i]->GetXaxis()->SetTitle("Si Energy [keV]"); hAngleEForwardCalTotal[i]->GetXaxis()->CenterTitle();
    hAngleEForwardCalTotal[i]->GetYaxis()->SetTitle("Vertex [mm]"); hAngleEForwardCalTotal[i]->GetYaxis()->CenterTitle();
    hAngleEForwardCalTotal[i]->GetYaxis()->SetTitleOffset(1.4);
  }

  // Vertex vs Angle Histograms Forward Si
  for(UInt_t i = 0; i < 10; i++) {
    TString name = Form("vertexAngleForward_d%d", i);
    hVertexAngleForward[i] = new TH2F(name, name, 500, -300, 300, 500, 0, 3);
  }

  // Cross Section Histograms
  s1 = new TH1F("s1", "Outside Forward", 70, 0, 6);
  s1->Sumw2();
  s1->GetXaxis()->SetTitle("Center of Mass Energy [MeV]"); s1->GetXaxis()->CenterTitle();
  s1->GetYaxis()->SetTitle("Cross Section [b/sr]"); s1->GetYaxis()->CenterTitle();
  s1->GetYaxis()->SetTitleOffset(1.2);
}

inline void Spectra::InitSiEForwardCalibration() {
  printf("Reading Si Energy Forward Calibrations File\n");
  std::ifstream inSiCalFile("siCalibration.dat");
  assert(inSiCalFile.is_open());
  Int_t var1, var2, var3;
  Double_t slope, intercept;
  while(inSiCalFile >> var1 >> var2 >> var3 >> slope >> intercept) {
    if(var1 < 10) siEForwardCalibration[var1][var2] = std::make_pair(slope, intercept);
    else siELeftCalibration[var1 - 10][var2] = std::make_pair(slope, intercept);
  }
  inSiCalFile.close();
}

inline void Spectra::InitCsIECalibration() {
  printf("Reading CsI Energy Calibrations File\n");
  std::ifstream inCsICalFile("csiCalibration.dat");
  assert(inCsICalFile.is_open());
  Int_t var1, var2;
  Double_t slope, intercept;
  while(inCsICalFile >> var1 >> var2 >> slope >> intercept) {
    if(var1 < 10) csiEForwardCalibration[var1] = std::make_pair(slope, intercept);
    else csiEForwardCalibration[var1 - 10] = std::make_pair(slope, intercept);
  }
  inCsICalFile.close();
}

inline void Spectra::InitCentralPadGainMatch() {
  printf("Reading Central Pad Gain Matching File\n");
  std::ifstream inGainFile("gainFile.dat");
  assert(inGainFile.is_open());
  Int_t varI, varJ;
  Double_t varScale;
  while(inGainFile >> varI >> varJ >> varScale) {
    scale[varI][varJ] = varScale;
  }
  inGainFile.close();
}

inline void Spectra::InitAverageBeamEnergy() {
  printf("Reading Average Beam Energy File\n");
  std::ifstream inBeamFile("averageBeamEnergy.out");
  assert(inBeamFile.is_open());
  Int_t varJ;
  Double_t varE;
  while(inBeamFile >> varJ >> varE) {
    averageBeamEnergy[varJ] = varE;
  }
}

inline void Spectra::InitVariables() {
  printf("Initializing Variables\n");

  // Make output file
  file = new TFile("spectra.root", "recreate");

  m1 = 8.; // AMU of projectile
  m2 = 1.; // AMU of target

  heightOffset = 6081.81 - 352.457; // Time of 0 height in chamber (in ns)
  driftVelocity = 0.05674449; // in mm/ns

  rowConversion = 1.75; // Side of row in mm
  rowConversionOffset = 1.75/2.; // Put the position in middle of row

  timeResolution = 40; // Buckets to time (40 ns time buckets)

  beamEnergy = 56.; // In MeV, after havar window
  density = 0.00038175; // in g/cm3, from LISE++ (Methane at 435 torr)
  numberB8 = 174809089.;

  distanceHavarToSilicon = 544.07; // Distance from Havar to Forward Silicon in mm

  // Initialize EnergyLoss
  boronMethane = new EnergyLoss("b8_methane.dat");
  protonMethane = new EnergyLoss("proton_methane.dat");

  // X positions of Forward Si Detector Quadrants
  siXPosForward[0][0] = 124. + 12.5;
  siXPosForward[0][1] = 124. - 12.5;
  siXPosForward[0][2] = 124. - 12.5;
  siXPosForward[0][3] = 124. + 12.5;
  siXPosForward[1][0] = 124. + 12.5;
  siXPosForward[1][1] = 124. - 12.5;
  siXPosForward[1][2] = 124. - 12.5;
  siXPosForward[1][3] = 124. + 12.5;
  siXPosForward[2][0] = 62.2 + 12.5;
  siXPosForward[2][1] = 62.2 + 12.5;
  siXPosForward[2][2] = 62.2 - 12.5;
  siXPosForward[2][3] = 62.2 - 12.5;
  siXPosForward[3][0] = 62.2 + 12.5;
  siXPosForward[3][1] = 62.2 - 12.5;
  siXPosForward[3][2] = 62.2 - 12.5;
  siXPosForward[3][3] = 62.2 + 12.5;
  siXPosForward[4][0] = 0. + 12.5;
  siXPosForward[4][1] = 0. + 12.5;
  siXPosForward[4][2] = 0. + 12.5;
  siXPosForward[4][3] = 0. + 12.5;
  siXPosForward[5][0] = 0. + 12.5;
  siXPosForward[5][1] = 0. + 12.5;
  siXPosForward[5][2] = 0. + 12.5;
  siXPosForward[5][3] = 0. + 12.5;
  siXPosForward[6][0] = 62.2 - 12.5;
  siXPosForward[6][1] = 62.2 - 12.5;
  siXPosForward[6][2] = 62.2 + 12.5;
  siXPosForward[6][3] = 62.2 + 12.5;
  siXPosForward[7][0] = 62.2 - 12.5;
  siXPosForward[7][1] = 62.2 + 12.5;
  siXPosForward[7][2] = 62.2 + 12.5;
  siXPosForward[7][3] = 62.2 - 12.5;
  siXPosForward[8][0] = 124. - 12.5;
  siXPosForward[8][1] = 124. + 12.5;
  siXPosForward[8][2] = 124. + 12.5;
  siXPosForward[8][3] = 124. - 12.5;
  siXPosForward[9][0] = 124. - 12.5;
  siXPosForward[9][1] = 124. + 12.5;
  siXPosForward[9][2] = 124. + 12.5;
  siXPosForward[9][3] = 124. - 12.5;
}

inline void Spectra::InitTree() {
  outTree = new TTree("outTree","Events from Digitizer");

  outTree->Branch("siDet", &siDet);
  outTree->Branch("siQuad", &siQuad);
  outTree->Branch("siChannel", &siChannel);
  outTree->Branch("siEnergy", &siEnergy);
  outTree->Branch("siEnergyCal", &siEnergyCal);
  outTree->Branch("siTime", &siTime);
  outTree->Branch("csiEnergy", &csiEnergy);
  outTree->Branch("csiEnergyCal", &csiEnergyCal);
  outTree->Branch("csiTime", &csiTime);
  outTree->Branch("totalEnergy", &totalEnergy);
  outTree->Branch("punchthrough", &punchthrough);
  outTree->Branch("dE", &dE);
  outTree->Branch("vertexPositionX", &vertexPositionX);
  outTree->Branch("vertexPositionY", &vertexPositionY);
  outTree->Branch("vertexPositionZ", &vertexPositionZ);
  outTree->Branch("angle", &angle);
  outTree->Branch("cmEnergy", &cmEnergy);
  return;
}

inline void Spectra::FillTree() {
  outTree->Fill();
}

inline void Spectra::WriteTree() {
  outTree->Write();
}

#endif // #ifdef Spectra_cxx
