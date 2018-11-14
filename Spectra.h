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
#include <math.h>
#include <sstream>

#include "CubicSpline.h"
#include "EnergyLoss.h"
#include "Hough2D.h"
#include "HoughTrack.h"
#include "TypeDef.h"

// Define the parameteric line equation
void line(Double_t t, std::vector<Double_t> p, Double_t &x, Double_t &y, Double_t &z) {
  x = p[0] + p[1]*t;
  y = t;
  z = p[2] + p[3]*t;
}

struct sortByRowMMChainStrip {
  inline Bool_t operator() (const mmChainStrip& struct1, const mmChainStrip& struct2) {
    return (struct1.row < struct2.row);
  }
};

struct sortByRowMMCenter {
  inline Bool_t operator() (const mmCenter& struct1, const mmCenter& struct2) {
    return (struct1.row < struct2.row);
  }
};

struct sortByRowMMTrack {
  inline Bool_t operator() (const mmTrack& struct1, const mmTrack& struct2) {
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
  Float_t         mmEnergy[2000]; //[mmMul]
  Float_t         mmPa[2000][5];  //[mmMul][5]

  // List of branches
  TBranch        *b_mmMul;    //!
  TBranch        *b_mmCobo;   //!
  TBranch        *b_mmAsad;   //!
  TBranch        *b_mmAget;   //!
  TBranch        *b_mmChan;   //!
  TBranch        *b_mmTime;   //!
  TBranch        *b_mmEnergy; //!
  TBranch        *b_mmPa;     //!

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

// Cuts
  TCutG *dEEForwardCut[10];

  TCutG *angleTotEnergyCut[10];

  TCutG *timeChainForwardCut[10][4];
  TCutG *timeStripForwardCut[10][4];

  TCutG *cwtE_CentralCut;
  TCutG *cwtE_CentralProtonCut;

// Histograms
private:
  void InitHistograms();

  // Ionization Chamber
  TH1F* hIonizationChamberE;
  TH1F* hIonizationChamberT;

  // Micromegas - Central
  TH2F* hMicroMegasCenterCumulative;
  TH2F* hMicroMegasCenterCumulativePosition;
  TH2F* hMicroMegasCenterCumulativePositionRaw;
  TH2F* hMicroMegasCenterEnergyCumulative;
  TH2F* hMicroMegasCenterEnergyAverage;
  TH2F* hMicroMegasCenterEnergyAverageScaled;
  TH2F* hMicroMegasCenterTime;
  TH2F* hMicroMegasCenterHeight;

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

  // Forward CsI Detectors
  TH1F* hCsIEForward[10];
  TH1F* hCsIEForwardCal[10];
  TH1F* hCsITForward[10];
  TH2F* hCsIETForward[10];

  // Si vs CsI
  TH2F* hSiCsIEForwardDet[10];
  TH2F* hSiCsIEForwardDetCal[10];
  TH2F* hSiCsIEForward[10][4];
  TH2F* hSiCsIEForwardCal[10][4];

  // Si + CsI vs Si (raw)
  TH2F* hSumSiEForwardDet[10];
  TH2F* hSumSiEForward[10][4];

  // Si + CsI vs CsI (raw)
  TH2F* hSumCsIEForwardDet[10];
  TH2F* hSumCsIEForward[10][4];

  // dE vs E Forward Detectors
  TH2F* hdEEForward[10];
  TH2F* hdEEForwardCal[10];
  TH2F* hdEEForwardCalTotal[10];

  // Hough Angle Forward Detectors
  TH1F* hHoughAngle[10];

  // Angle vs E Forward Detectors
  TH2F* hAngleEForward[10];
  TH2F* hAngleEForwardCal[10];
  TH2F* hAngleEForwardCalTotal[10];
  TH2F* hAngleEForwardProtonEnergy[10];
  TH2F* hAngleEForwardCMEnergy[10];

  // Vertex vs E Forward Detectors
  TH2F* hVertexSiEForward[10];
  TH2F* hVertexSiEForwardCal[10];
  TH2F* hVertexSiEForwardCalTotal[10];
  TH2F* hVertexSiETotalRegion3;
  TH2F* hVertexCMERegion3;

  // Vertex vs Angle Forward Detectors
  TH2F* hVertexAngleForward[10];

  // Time vs Column Number Forward Detectors
  TH2F* hTimeChainForward[10][4];

  // Time vs Strip Number Forward Detectors
  TH2F* hTimeStripForward[10][4];

  // Time vs Central Region Forward Detectors
  TH2F* hTimeCentralForward[10];

  // Forward Wall XZ Hit Positions
  TH2F* hHitPositionsXZForward;
  TH2F* hHitPositionsXZForwardInd[10];

  // Central Micromegas Energy vs Pa
  TH2F* hCWTECentral;

  // Silicon Energy vs Pa
  TH2F* hCWTSiE;

  // Max Peak Location in Central Pad vs Si E Forward Detectors
  TH2F* hMaxPeakSiE[10];

  // Max Peak Location in Central Pad vs Average Peak Energy Forward Detectors
  TH2F* hMaxPeakAvgPeakE[10];

  // Max Peak Energy in Central Pad vs Average Peak Energy Forward Detectors
  TH2F* hMaxPeakEAvgPeakE[10];

  // Average Peak Energy vs Si E Forward Detectors
  TH2F* hAvgPeakESiE[10];

  // Max Peak Energy in Central Pad vs Derivative around Peak Forward Detectors
  TH2F* hMaxPeakEDerivPeak[10];

  // Max Peak Location in Central Pad vs Difference between Peak Location and Derivative Max
  TH2F* hMaxPeakDerivDiff[10];

  // Derivative Location in Central Pad vs Difference between Peak Location and Derivative Max
  TH2F* hDerivPeakDerivDiff[10];

  // Cross Section
  TH1F* s1;

  void WriteHistograms();

// TCanvas

  void InitCanvas();

  // For central pads
  void DrawCenterEnergyCanvas(Int_t count, std::vector<mmCenter> centerMatched_, std::vector<mmTrack> centerBeamTotal_);
  Int_t totalCenterEnergyCanvas;
  Int_t centerEnergyCanvasNum, centerEnergyCanvasXYNum, centerEnergyCanvasXNum, centerEnergyCanvasYNum;
  TCanvas* centerEnergyCanvas[5];

  // Beam track
  void DrawCenterBeamCanvas(Int_t count, std::vector<mmTrack> centerBeamTrack_, std::vector<Double_t> pars,
      Int_t lastRow);
  Int_t totalCenterBeamCanvas;
  Int_t centerBeamCanvasNum, centerBeamCanvasXYNum, centerBeamCanvasXNum, centerBeamCanvasYNum;
  TCanvas* centerBeamCanvas[5];

  void DrawCenterNoiseCanvas(Int_t count, std::vector<mmCenter> centerReduced_,
                             std::vector<mmCenter> centerReducedNoise_);
  Int_t totalCenterNoiseCanvas;
  Int_t centerNoiseCanvasNum, centerNoiseCanvasXYNum, centerNoiseCanvasXNum, centerNoiseCanvasYNum;
  TCanvas* centerNoiseCanvas[10];

  // Draw central pad raw energy vs running average energy
  void DrawCenterEnergyRunningAverageCanvas(Int_t count, std::vector<mmTrack> rawTrack_,
                                             std::vector<mmTrack> averageTrack_);
  Int_t totalCenterEnergyRunningCanvas;
  Int_t centerEnergyRunningCanvasNum, centerEnergyRunningCanvasXYNum, centerEnergyRunningCanvasXNum, centerEnergyRunningCanvasYNum;
  TCanvas* centerEnergyRunningCanvas[5];

  // Draw three and five point derivatives with central region with respect to energy
  void DrawCenterEnergyDeriv(Int_t count, std::vector<centerDeriv> threePoint_, std::vector<centerDeriv> fivePoint_);
  Int_t totalCenterEnergyDerivCanvas;
  Int_t centerEnergyDerivCanvasNum, centerEnergyDerivCanvasXYNum, centerEnergyDerivCanvasXNum, centerEnergyDerivCanvasYNum;
  TCanvas* centerEnergyDerivCanvas[5];

  void WriteCanvas();

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
  Double_t averageBeamEnergy[6][128];

// General Methods
  Bool_t AnalysisForwardCentral(std::vector<mmCenter> centerMatched_, std::vector<mmTrack> centerBeamTotal_,
                              std::vector<mmTrack> centerProton_, std::map<Int_t, Double_t> centralPadTotalEnergy);
  Bool_t AnalysisForwardSide(std::vector<mmCenter> centerMatched_, std::vector<mmTrack> centerBeamTotal_,
                           std::vector<mmChainStrip> leftChain_, std::vector<mmChainStrip> leftStrip_,
                           std::vector<mmChainStrip> rightChain_,
                           std::vector<mmChainStrip> rightStrip_);
  Bool_t AnalysisLeftSide();

// MicroMegas Functions

  // Center Functions
  std::vector<mmCenter> CenterReduceNoise(std::vector<mmCenter> center);
  void CorrectCenterEnergy(std::vector<mmTrack> &centerMatched_, std::vector<Double_t> parsBeam, Int_t lastRow);
  Double_t GaussianCDF(Double_t x, Double_t mean, Double_t sigma);
  void FindMaxCentralEnergy(std::vector<mmTrack> centerMatched_, Int_t &maxEnergyRow, Double_t &maxEnergy,
                            Double_t &averageMaxEnergy, Double_t &maxDeriv);
  std::vector<mmTrack> GetRunningEnergyAverageThree(std::vector<mmTrack> centerMatched_);
  std::vector<mmTrack> GetRunningEnergyAverageFive(std::vector<mmTrack> centerMatched_);
  Bool_t CenterOnlyOneColumn(std::vector<mmTrack> centerMatched_);
  std::vector<centerDeriv> CenterEnergyThreePointDeriv(std::vector<mmTrack> centerMatched_);
  std::vector<centerDeriv> CenterEnergyFivePointDeriv(std::vector<mmTrack> centerMatched_);
  std::pair<Int_t, Int_t> CenterGetDerivMax(std::vector<centerDeriv> threePoint_, std::vector<centerDeriv> fivePoint_);

  // Side Functions
  void ChainStripMatch(std::vector<mmTrack> &chainStripMatched, std::vector<mmTrack> &chainStripRaw,
                       std::vector<mmChainStrip> chain_,
                       std::vector<mmChainStrip> strip_, Bool_t leftSide, Double_t siTime);
  size_t ChainStripTime0TimeBuckets(std::vector<mmTrack> matched);
  size_t ChainStripNumberTimeBuckets(std::vector<mmChainStrip> chain, std::vector<mmChainStrip> strip);
  void ChainStripMatchingOutward(std::vector<mmTrack> &chainStripMatched, std::vector<mmChainStrip> chain,
                                 std::vector<mmChainStrip> strip, Bool_t leftSide, Double_t siTime);
  void ChainStripMatchingBox(std::vector<mmTrack> &chainStripMatched, std::vector<mmChainStrip> chain,
                             std::vector<mmChainStrip> strip, Bool_t leftSide, Double_t siTime);
  void ChainStripMatchingBoxTime0(std::vector<mmTrack> &chainStripMatched, std::vector<mmTrack> time0);
  void ChainStripMatchingTime(std::vector<mmTrack> &chainStripMatched, std::vector<mmChainStrip> chain,
                              std::vector<mmChainStrip> strip, Bool_t leftSide, Double_t siTime, Int_t timeWindow);
  void ChainStripMatchingTimeSlopeFit(std::vector<mmTrack> &chainStripMatched, std::vector<mmChainStrip> chain,
                                      std::vector<mmChainStrip> strip, Bool_t leftSide, Double_t siTime,
                                      Double_t timeWindow);
  void ChainStripMatchingTimeSlopeHough(std::vector<mmTrack> &chainStripMatched, std::vector<mmChainStrip> chain,
                                        std::vector<mmChainStrip> strip, Bool_t leftSide, Double_t siTime,
                                        Double_t timeWindow);

// Visualize Hough Transform
  void GetMinMaxD(std::vector<mmTrack> initPoints, Int_t &minXY, Int_t &maxXY, Int_t &minYZ, Int_t &maxYZ);
  void VisualizeHough(std::vector<mmTrack> initPoints, TH2I* fXY, TH2I* fYZ);
  void GetMinMaxDRestricted(std::vector<mmTrack> initPoints, Int_t &minXY, Int_t &maxXY, Int_t &minYZ, Int_t &maxYZ, Int_t siDet);
  void VisualizeHoughRestricted(std::vector<mmTrack> initPoints, TH2I* fXY, TH2I* fYZ, Int_t siDet);
  void GetHoughStdDevXYRestricted(std::vector<mmTrack> initPoints, std::vector<Double_t> &angle_, std::vector<Double_t> &stdDev_, Int_t siDet);

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
  Double_t m1;
  Double_t m2;
  Double_t rowConversion;
  Double_t rowConversionOffset;
  Double_t heightOffset;
  Double_t driftVelocity;
  Double_t timeResolution;
  Long64_t entry;
  Double_t beamEnergy;
  Double_t density;
  Double_t distanceHavarToSilicon;
  Double_t distanceHavarToMM;
  Double_t numberB8;
  Double_t siXPosForward[10][4];
  Double_t siYPosForward;
  Double_t gasPositionResolution;
  std::pair<Double_t, Double_t> mmColumnSize[6];

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
  Double_t siEnergy;
  Double_t siEnergyCal;
  Double_t siTime;
  Double_t csiEnergy;
  Double_t csiEnergyCal;
  Double_t csiTime;
  Double_t totalEnergy;
  Bool_t punchthrough;
  Double_t dE;
  Double_t vertexPositionX;
  Double_t vertexPositionY;
  Double_t vertexPositionZ;
  Double_t angle;
  Double_t cmEnergy;
  Double_t siPosX;
  Double_t siPosY;
  Double_t siPosZ;

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
  fChain->SetBranchAddress("mmPa", mmPa, &b_mmPa);
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
  hMicroMegasCenterCumulative = new TH2F("MM_Center_Cumulative", "MM_Center_Cumulative", 10, -5, 5, 128, 0, 128);
  hMicroMegasCenterCumulativePosition = new TH2F("MM_Center_Cumulative_Position", "MM_Center_Cumulative_Position", 24, -21, 21, 143, 0, 250.25);
  hMicroMegasCenterCumulativePositionRaw = new TH2F("MM_Center_Cumulative_Position_Raw", "MM_Center_Cumulative_Position_Raw", 24, -21, 21, 143, 0, 250.25);
  hMicroMegasCenterEnergyCumulative = new TH2F("MM_Center_Energy_Cumulative", "MM_Center_Energy_Cumulative", 20, -10, 10, 128, 0, 128);
  hMicroMegasCenterEnergyAverage = new TH2F("MM_Center_Energy_Average", "MM_Center_Energy_Average", 20, -10, 10, 128, 0, 128);
  hMicroMegasCenterEnergyAverageScaled = new TH2F("MM_Center_Energy_Average_Scaled", "MM_Center_Energy_Average_Scaled", 20, -10, 10, 128, 0, 128);
  hMicroMegasCenterTime = new TH2F("MM_Center_Time", "MM_Center_Time", 130, -1, 129, 250, 0, 10000);
  hMicroMegasCenterHeight = new TH2F("MM_Center_Height", "MM_Center_Height", 130, -1, 129, 160, -200, 200);

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
    hSiEForwardDetCal[i] = new TH1F(name, name, 500, 0, 16000);
    hSiEForwardDetCal[i]->GetXaxis()->SetTitle("Energy [keV]"); hSiEForwardDetCal[i]->GetXaxis()->CenterTitle();
    hSiEForwardDetCal[i]->GetYaxis()->SetTitle("Counts"); hSiEForwardDetCal[i]->GetYaxis()->CenterTitle();
    hSiEForwardDetCal[i]->GetYaxis()->SetTitleOffset(1.4);

    name = Form("siTForward_d%d", i);
    hSiTForwardDet[i] = new TH1F(name, name, 500, 0, 20000);
    hSiTForwardDet[i]->GetXaxis()->SetTitle("Time [ns]"); hSiTForwardDet[i]->GetXaxis()->CenterTitle();
    hSiTForwardDet[i]->GetYaxis()->SetTitle("Counts"); hSiTForwardDet[i]->GetYaxis()->CenterTitle();
    hSiTForwardDet[i]->GetYaxis()->SetTitleOffset(1.4);

    for(UInt_t j = 0; j < 4; j++) {
      TString name = Form("siEForward_d%d_q%d", i, j);
      hSiEForward[i][j] = new TH1F(name, name, 500, 0, 4000);
      hSiEForward[i][j]->GetXaxis()->SetTitle("Energy [channels]"); hSiEForward[i][j]->GetXaxis()->CenterTitle();
      hSiEForward[i][j]->GetYaxis()->SetTitle("Counts"); hSiEForward[i][j]->GetYaxis()->CenterTitle();
      hSiEForward[i][j]->GetYaxis()->SetTitleOffset(1.4);

      name = Form("siEForwardCal_d%d_q%d", i, j);
      hSiEForwardCal[i][j] = new TH1F(name, name, 500, 0, 16000);
      hSiEForwardCal[i][j]->GetXaxis()->SetTitle("Energy [keV]"); hSiEForwardCal[i][j]->GetXaxis()->CenterTitle();
      hSiEForwardCal[i][j]->GetYaxis()->SetTitle("Counts"); hSiEForwardCal[i][j]->GetYaxis()->CenterTitle();
      hSiEForwardCal[i][j]->GetYaxis()->SetTitleOffset(1.4);

      name = Form("siTForward_d%d_q%d", i, j);
      hSiTForward[i][j] = new TH1F(name, name, 500, 0, 20000);
      hSiTForward[i][j]->GetXaxis()->SetTitle("Time [ns]"); hSiTForward[i][j]->GetXaxis()->CenterTitle();
      hSiTForward[i][j]->GetYaxis()->SetTitle("Counts"); hSiTForward[i][j]->GetYaxis()->CenterTitle();
      hSiTForward[i][j]->GetYaxis()->SetTitleOffset(1.4);
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
    hCsIEForwardCal[i] = new TH1F(name, name, 500, 0, 16000);
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

  // Si vs CsI Forward Detectors
  for(UInt_t i = 0; i < 10; i++) {
    TString name = Form("siCsIEForward_d%d", i);
    hSiCsIEForwardDet[i] = new TH2F(name, name, 500, 0, 4000, 500, 0, 10000);
    hSiCsIEForwardDet[i]->GetXaxis()->SetTitle("Si Energy [channels]"); hSiCsIEForwardDet[i]->GetXaxis()->CenterTitle();
    hSiCsIEForwardDet[i]->GetYaxis()->SetTitle("CsI Energy [channels]"); hSiCsIEForwardDet[i]->GetYaxis()->CenterTitle();
    hSiCsIEForwardDet[i]->GetYaxis()->SetTitleOffset(1.4);

    name = Form("siCsIEForwardCal_d%d", i);
    hSiCsIEForwardDetCal[i] = new TH2F(name, name, 500, 0, 16000, 500, 0, 10000);
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
      hSiCsIEForwardCal[i][j] = new TH2F(name, name, 500, 0, 16000, 500, 0, 10000);
      hSiCsIEForwardCal[i][j]->GetXaxis()->SetTitle("Si Energy [keV]"); hSiCsIEForwardCal[i][j]->GetXaxis()->CenterTitle();
      hSiCsIEForwardCal[i][j]->GetYaxis()->SetTitle("CsI Energy [channels]"); hSiCsIEForwardCal[i][j]->GetYaxis()->CenterTitle();
      hSiCsIEForwardCal[i][j]->GetYaxis()->SetTitleOffset(1.4);
    }
  }

  // Si + CsI (raw) vs Si (raw) Forward Detectors
  for(UInt_t i = 0; i < 10; i++) {
    TString name = Form("hSumSiEForward_d%d", i);
    hSumSiEForwardDet[i] = new TH2F(name, name, 500, 0, 8000, 500, 0, 4000);
    hSumSiEForwardDet[i]->GetXaxis()->SetTitle("Si + CsI [channels]"); hSumSiEForwardDet[i]->GetXaxis()->CenterTitle();
    hSumSiEForwardDet[i]->GetYaxis()->SetTitle("Si [channels]"); hSumSiEForwardDet[i]->GetYaxis()->CenterTitle();
    hSumSiEForwardDet[i]->GetYaxis()->SetTitleOffset(1.4);
    for(UInt_t j = 0; j < 4; j++) {
      name = Form("hSumSiEForward_d%d_q%d", i, j);
      hSumSiEForward[i][j] = new TH2F(name, name, 500, 0, 8000, 500, 0, 4000);
      hSumSiEForward[i][j]->GetXaxis()->SetTitle("Si + CsI [channels]"); hSumSiEForward[i][j]->GetXaxis()->CenterTitle();
      hSumSiEForward[i][j]->GetYaxis()->SetTitle("Si [channels]"); hSumSiEForward[i][j]->GetYaxis()->CenterTitle();
      hSumSiEForward[i][j]->GetYaxis()->SetTitleOffset(1.4);
    }
  }

  // Si + CsI (raw) vs CsI (raw) Forward Detectors
  for(UInt_t i = 0; i < 10; i++) {
    TString name = Form("hSumCsIEForward_d%d", i);
    hSumCsIEForwardDet[i] = new TH2F(name, name, 500, 0, 8000, 500, 0, 4000);
    hSumCsIEForwardDet[i]->GetXaxis()->SetTitle("Si + CsI [channels]"); hSumCsIEForwardDet[i]->GetXaxis()->CenterTitle();
    hSumCsIEForwardDet[i]->GetYaxis()->SetTitle("CsI [channels]"); hSumCsIEForwardDet[i]->GetYaxis()->CenterTitle();
    hSumCsIEForwardDet[i]->GetYaxis()->SetTitleOffset(1.4);
    for(UInt_t j = 0; j < 4; j++) {
      name = Form("hSumCsIEForward_d%d_q%d", i, j);
      hSumCsIEForward[i][j] = new TH2F(name, name, 500, 0, 8000, 500, 0, 4000);
      hSumCsIEForward[i][j]->GetXaxis()->SetTitle("Si + CsI [channels]"); hSumCsIEForward[i][j]->GetXaxis()->CenterTitle();
      hSumCsIEForward[i][j]->GetYaxis()->SetTitle("CsI [channels]"); hSumCsIEForward[i][j]->GetYaxis()->CenterTitle();
      hSumCsIEForward[i][j]->GetYaxis()->SetTitleOffset(1.4);
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
    hdEEForwardCalTotal[i]->GetXaxis()->SetTitle("Total Energy [keV]"); hdEEForwardCalTotal[i]->GetXaxis()->CenterTitle();
    hdEEForwardCalTotal[i]->GetYaxis()->SetTitle("dE [channels]"); hdEEForwardCalTotal[i]->GetYaxis()->CenterTitle();
    hdEEForwardCalTotal[i]->GetYaxis()->SetTitleOffset(1.4);
  }

  // Hough Angle - Forward Detectors
  for(UInt_t i = 0; i < 10; i++) {
    TString name = Form("houghAngleSiForward_d%d", i);
    hHoughAngle[i] = new TH1F(name, name, 180, 0, 180);
    hHoughAngle[i]->GetXaxis()->SetTitle("Hough Angle"); hHoughAngle[i]->GetXaxis()->CenterTitle();
  }

  // Angle vs E Histograms hAngleEForward[10]
  for(UInt_t i = 0; i < 10; i++) {
    TString name = Form("angleSiEForward_d%d", i);
    hAngleEForward[i] = new TH2F(name, name, 200, 0, 4000, 75, 0, 1.7);
    hAngleEForward[i]->GetXaxis()->SetTitle("Si Energy [channels]"); hAngleEForward[i]->GetXaxis()->CenterTitle();
    hAngleEForward[i]->GetYaxis()->SetTitle("Angle [rad]"); hAngleEForward[i]->GetYaxis()->CenterTitle();
    hAngleEForward[i]->GetYaxis()->SetTitleOffset(1.4);

    name = Form("angleSiEForwardCal_d%d", i);
    hAngleEForwardCal[i] = new TH2F(name, name, 200, 0, 20000, 75, 0, 1.7);
    hAngleEForwardCal[i]->GetXaxis()->SetTitle("Si Energy [keV]"); hAngleEForwardCal[i]->GetXaxis()->CenterTitle();
    hAngleEForwardCal[i]->GetYaxis()->SetTitle("Angle [rad]"); hAngleEForwardCal[i]->GetYaxis()->CenterTitle();
    hAngleEForwardCal[i]->GetYaxis()->SetTitleOffset(1.4);

    name = Form("angleSiEForwardCalTotal_d%d", i);
    hAngleEForwardCalTotal[i] = new TH2F(name, name, 100, 0, 20000, 75, 0, 1.7);
    hAngleEForwardCalTotal[i]->GetXaxis()->SetTitle("Total Energy [keV]"); hAngleEForwardCalTotal[i]->GetXaxis()->CenterTitle();
    hAngleEForwardCalTotal[i]->GetYaxis()->SetTitle("Angle [rad]"); hAngleEForwardCalTotal[i]->GetYaxis()->CenterTitle();
    hAngleEForwardCalTotal[i]->GetYaxis()->SetTitleOffset(1.4);

    name = Form("angleProtonEForward_d%d", i);
    hAngleEForwardProtonEnergy[i] = new TH2F(name, name, 75, 0, 20, 50, 0, 1.7);
    hAngleEForwardProtonEnergy[i]->GetXaxis()->SetTitle("Proton Energy at Vertex [MeV]"); hAngleEForwardProtonEnergy[i]->GetXaxis()->CenterTitle();
    hAngleEForwardProtonEnergy[i]->GetYaxis()->SetTitle("Angle [rad]"); hAngleEForwardProtonEnergy[i]->GetYaxis()->CenterTitle();
    hAngleEForwardProtonEnergy[i]->GetYaxis()->SetTitleOffset(1.4);

    name = Form("angleCMEForward_d%d", i);
    hAngleEForwardCMEnergy[i] = new TH2F(name, name, 100, 0, 8, 75, 0, 1.7);
    hAngleEForwardCMEnergy[i]->GetXaxis()->SetTitle("CM Energy [MeV]"); hAngleEForwardCMEnergy[i]->GetXaxis()->CenterTitle();
    hAngleEForwardCMEnergy[i]->GetYaxis()->SetTitle("Angle [rad]"); hAngleEForwardCMEnergy[i]->GetYaxis()->CenterTitle();
    hAngleEForwardCMEnergy[i]->GetYaxis()->SetTitleOffset(1.4);
  }

  // Vertex vs E Histograms
  for(UInt_t i = 0; i < 10; i++) {
    TString name = Form("vertexSiEForward_d%d", i);
    hVertexSiEForward[i] = new TH2F(name, name, 500, 0, 4000, 500, -450, 300);
    hVertexSiEForward[i]->GetXaxis()->SetTitle("Si Energy [channels]"); hVertexSiEForward[i]->GetXaxis()->CenterTitle();
    hVertexSiEForward[i]->GetYaxis()->SetTitle("Vertex [mm]"); hVertexSiEForward[i]->GetYaxis()->CenterTitle();
    hVertexSiEForward[i]->GetYaxis()->SetTitleOffset(1.4); hVertexSiEForward[i]->SetStats(false);

    name = Form("vertexSiEForwardCal_d%d", i);
    hVertexSiEForwardCal[i] = new TH2F(name, name, 500, 0, 16000, 500, -450, 300);
    hVertexSiEForwardCal[i]->GetXaxis()->SetTitle("Si Energy [keV]"); hVertexSiEForwardCal[i]->GetXaxis()->CenterTitle();
    hVertexSiEForwardCal[i]->GetYaxis()->SetTitle("Vertex [mm]"); hVertexSiEForwardCal[i]->GetYaxis()->CenterTitle();
    hVertexSiEForwardCal[i]->GetYaxis()->SetTitleOffset(1.4); hVertexSiEForwardCal[i]->SetStats(false);

    name = Form("vertexSiEForwardCalTotal_d%d", i);
    hVertexSiEForwardCalTotal[i] = new TH2F(name, name, 500, 0, 20000, 500, -450, 300);
    hVertexSiEForwardCalTotal[i]->GetXaxis()->SetTitle("Total Energy [keV]"); hVertexSiEForwardCalTotal[i]->GetXaxis()->CenterTitle();
    hVertexSiEForwardCalTotal[i]->GetYaxis()->SetTitle("Vertex [mm]"); hVertexSiEForwardCalTotal[i]->GetYaxis()->CenterTitle();
    hVertexSiEForwardCalTotal[i]->GetYaxis()->SetTitleOffset(1.4); hVertexSiEForwardCalTotal[i]->SetStats(false);
  }

  hVertexSiETotalRegion3 = new TH2F("vertexSiERegion3", "vertexSiERegion3", 500, 0, 20000, 500, -450, 300);
  hVertexSiETotalRegion3->GetXaxis()->SetTitle("Total Energy [keV]"); hVertexSiETotalRegion3->GetXaxis()->CenterTitle();
  hVertexSiETotalRegion3->GetYaxis()->SetTitle("Vertex [mm]"); hVertexSiETotalRegion3->GetYaxis()->CenterTitle();
  hVertexSiETotalRegion3->GetYaxis()->SetTitleOffset(1.4); hVertexSiETotalRegion3->SetStats(false);

  hVertexCMERegion3 = new TH2F("vertexCMERegion3", "vertexCMERegion3", 500, 0, 6, 500, -450, 300);
  hVertexCMERegion3->GetXaxis()->SetTitle("CM Energy [MeV]"); hVertexCMERegion3->GetXaxis()->CenterTitle();
  hVertexCMERegion3->GetYaxis()->SetTitle("Vertex [mm]"); hVertexCMERegion3->GetYaxis()->CenterTitle();
  hVertexCMERegion3->GetYaxis()->SetTitleOffset(1.4); hVertexCMERegion3->SetStats(false);

  // Vertex vs Angle Histograms Forward Si
  for(UInt_t i = 0; i < 10; i++) {
    TString name = Form("vertexAngleForward_d%d", i);
    hVertexAngleForward[i] = new TH2F(name, name, 500, -300, 300, 500, 0, 3);
  }

  // Time vs Column Number Forward Detectors
  for(UInt_t i = 0; i < 10; i++) {
    for(UInt_t j = 0; j < 4; j++) {
      TString name = Form("timeColumnForward_d%d_q%d", i, j);
      hTimeChainForward[i][j] = new TH2F(name, name, 68, -2, 66, 75, 0, 3000);
      hTimeChainForward[i][j]->GetXaxis()->SetTitle("Chain #"); hTimeChainForward[i][j]->GetXaxis()->CenterTitle();
      hTimeChainForward[i][j]->GetYaxis()->SetTitle("Time [ns]"); hTimeChainForward[i][j]->GetYaxis()->CenterTitle();
      hTimeChainForward[i][j]->GetYaxis()->SetTitleOffset(1.4);
      hTimeChainForward[i][j]->SetStats(false);
    }
  }

  // Time vs Strip Number Forward Detectors
  for(UInt_t i = 0; i < 10; i++) {
    for(UInt_t j = 0; j < 4; j++) {
      TString name = Form("timeStripForward_d%d_q%d", i, j);
      hTimeStripForward[i][j] = new TH2F(name, name, 68, -2, 66, 75, 0, 3000);
      hTimeStripForward[i][j]->GetXaxis()->SetTitle("Strip #"); hTimeStripForward[i][j]->GetXaxis()->CenterTitle();
      hTimeStripForward[i][j]->GetYaxis()->SetTitle("Time [ns]"); hTimeStripForward[i][j]->GetYaxis()->CenterTitle();
      hTimeStripForward[i][j]->GetYaxis()->SetTitleOffset(1.4);
      hTimeStripForward[i][j]->SetStats(false);
    }
  }

  // Time vs Central Region Forward Detectors
  for(UInt_t i = 0; i < 10; i++) {
    TString name = Form("timeCentralForward_d%d", i);
    hTimeCentralForward[i] = new TH2F(name, name, 130, -2, 128, 75, 0, 3000);
    hTimeCentralForward[i]->GetXaxis()->SetTitle("Row #"); hTimeCentralForward[i]->GetXaxis()->CenterTitle();
    hTimeCentralForward[i]->GetYaxis()->SetTitle("Time [ns]"); hTimeCentralForward[i]->GetYaxis()->CenterTitle();
    hTimeCentralForward[i]->GetYaxis()->SetTitleOffset(1.4);
    hTimeCentralForward[i]->SetStats(false);
  }

  // Forward Wall XZ Hit Positions
  hHitPositionsXZForward = new TH2F("hitPositionXZForward", "hitPositionXZForward", 200, -200, 200, 100, -100, 100);
  hHitPositionsXZForward->GetXaxis()->SetTitle("X [mm]"); hHitPositionsXZForward->GetXaxis()->CenterTitle();
  hHitPositionsXZForward->GetYaxis()->SetTitle("Z [mm]"); hHitPositionsXZForward->GetYaxis()->CenterTitle();
  hHitPositionsXZForward->GetYaxis()->SetTitleOffset(1.4);
  hHitPositionsXZForward->SetStats(false);

  for(UInt_t i = 0; i < 10; i++) {
    TString name = Form("hitPositionXZForward_d%d", i);
    hHitPositionsXZForwardInd[i] = new TH2F(name, name, 200, -200, 200, 150, -150, 150);
    hHitPositionsXZForwardInd[i]->GetXaxis()->SetTitle("X [mm]"); hHitPositionsXZForwardInd[i]->GetXaxis()->CenterTitle();
    hHitPositionsXZForwardInd[i]->GetYaxis()->SetTitle("Z [mm]"); hHitPositionsXZForwardInd[i]->GetYaxis()->CenterTitle();
    hHitPositionsXZForwardInd[i]->GetYaxis()->SetTitleOffset(1.4);
    hHitPositionsXZForwardInd[i]->SetStats(false);
  }

  hCWTECentral = new TH2F("cwtECentral", "cwtECentral", 500, 0, 4000, 1000, 0, 0.1);
  hCWTECentral->GetXaxis()->SetTitle("Energy [channels]"); hCWTECentral->GetXaxis()->CenterTitle();
  hCWTECentral->GetYaxis()->SetTitle("CWT"); hCWTECentral->GetYaxis()->CenterTitle();
  hCWTECentral->GetYaxis()->SetTitleOffset(1.4);
  hCWTECentral->SetStats(false);

  hCWTSiE = new TH2F("cwtSiE", "cwtSiE", 500, 0, 4000, 1000, 0, 0.1);
  hCWTSiE->GetXaxis()->SetTitle("Energy [channels]"); hCWTSiE->GetXaxis()->CenterTitle();
  hCWTSiE->GetYaxis()->SetTitle("CWT"); hCWTSiE->GetYaxis()->CenterTitle();
  hCWTSiE->GetYaxis()->SetTitleOffset(1.4);
  hCWTSiE->SetStats(false);

  // Max Peak Location Central Pad vs Si E Forward Wall
  for(UInt_t i = 0; i < 10; i++) {
    TString name = Form("maxPeakSiE_d%d", i);
    hMaxPeakSiE[i] = new TH2F(name, name, 128, 0, 128, 500, 0, 15000);
    hMaxPeakSiE[i]->GetXaxis()->SetTitle("Max Peak [row #]"); hMaxPeakSiE[i]->GetXaxis()->CenterTitle();
    hMaxPeakSiE[i]->GetYaxis()->SetTitle("Si Energy [keV]"); hMaxPeakSiE[i]->GetYaxis()->CenterTitle();
    hMaxPeakSiE[i]->GetYaxis()->SetTitleOffset(1.4);
    hMaxPeakSiE[i]->SetStats(false);
  }

  // Max Peak Location Central Pad vs Average E Forward Wall
  for(UInt_t i = 0; i < 10; i++) {
    TString name = Form("maxPeakAvgE_d%d", i);
    hMaxPeakAvgPeakE[i] = new TH2F(name, name, 128, 0, 128, 128, 0, 4096);
    hMaxPeakAvgPeakE[i]->GetXaxis()->SetTitle("Max Peak [row #]"); hMaxPeakAvgPeakE[i]->GetXaxis()->CenterTitle();
    hMaxPeakAvgPeakE[i]->GetYaxis()->SetTitle("Average Peak Energy [channels]"); hMaxPeakAvgPeakE[i]->GetYaxis()->CenterTitle();
    hMaxPeakAvgPeakE[i]->GetYaxis()->SetTitleOffset(1.4);
    hMaxPeakAvgPeakE[i]->SetStats(false);
  }

  // Max Peak Energy in Central Pad vs Average Peak Energy Forward Detectors
  for(UInt_t i = 0; i < 10; i++) {
    TString name = Form("maxPeakEAvgE_d%d", i);
    hMaxPeakEAvgPeakE[i] = new TH2F(name, name, 256, 0, 4096, 128, 0, 4096);
    hMaxPeakEAvgPeakE[i]->GetXaxis()->SetTitle("Max Peak Energy [channels]"); hMaxPeakEAvgPeakE[i]->GetXaxis()->CenterTitle();
    hMaxPeakEAvgPeakE[i]->GetYaxis()->SetTitle("Average Peak Energy [channels]"); hMaxPeakEAvgPeakE[i]->GetYaxis()->CenterTitle();
    hMaxPeakEAvgPeakE[i]->GetYaxis()->SetTitleOffset(1.4);
    hMaxPeakEAvgPeakE[i]->SetStats(false);
  }

  // Average Peak Energy vs Si E Forward Detectors
  for(UInt_t i = 0; i < 10; i++) {
    TString name = Form("avgPeakESiE_d%d", i);
    hAvgPeakESiE[i] = new TH2F(name, name, 128, 0, 4096, 500, 0, 15000);
    hAvgPeakESiE[i]->GetXaxis()->SetTitle("Average Peak Energy [channels]"); hAvgPeakESiE[i]->GetXaxis()->CenterTitle();
    hAvgPeakESiE[i]->GetYaxis()->SetTitle("Si Energy [keV]"); hAvgPeakESiE[i]->GetYaxis()->CenterTitle();
    hAvgPeakESiE[i]->GetYaxis()->SetTitleOffset(1.4);
    hAvgPeakESiE[i]->SetStats(false);
  }

  // Max Peak Energy in Central Pad vs Derivative around Peak Forward Detectors
  for(UInt_t i = 0; i < 10; i++) {
    TString name = Form("maxPeakEDeriv_d%d", i);
    hMaxPeakEDerivPeak[i] = new TH2F(name, name, 256, 0, 4096, 128, 0, 2048);
    hMaxPeakEDerivPeak[i]->GetXaxis()->SetTitle("Max Peak Energy [channels]"); hMaxPeakEDerivPeak[i]->GetXaxis()->CenterTitle();
    hMaxPeakEDerivPeak[i]->GetYaxis()->SetTitle("Derivative at Peak"); hMaxPeakEDerivPeak[i]->GetYaxis()->CenterTitle();
    hMaxPeakEDerivPeak[i]->GetYaxis()->SetTitleOffset(1.4);
    hMaxPeakEDerivPeak[i]->SetStats(false);
  }

  // Max Peak Location in Central Pad vs Difference between Peak Location and Derivative Max
  for(UInt_t i = 0; i < 10; i++) {
    TString name = Form("maxPeakDerivDiff_d%d", i);
    hMaxPeakDerivDiff[i] = new TH2F(name, name, 128, 0, 128, 128, 0, 128);
    hMaxPeakDerivDiff[i]->GetXaxis()->SetTitle("Max Peak Location [row]"); hMaxPeakDerivDiff[i]->GetXaxis()->CenterTitle();
    hMaxPeakDerivDiff[i]->GetYaxis()->SetTitle("Diff in Derivative and Max Peak [row]"); hMaxPeakDerivDiff[i]->GetYaxis()->CenterTitle();
    hMaxPeakDerivDiff[i]->GetYaxis()->SetTitleOffset(1.4);
    hMaxPeakDerivDiff[i]->SetStats(false);
  }

  // Derivative Location in Central Pad vs Difference between Peak Location and Derivative Max
  for(UInt_t i = 0; i < 10; i++) {
    TString name = Form("derivPeakDerivDiff_d%d", i);
    hDerivPeakDerivDiff[i] = new TH2F(name, name, 128, 0, 128, 128, 0, 128);
    hDerivPeakDerivDiff[i]->GetXaxis()->SetTitle("Derivative Peak Location [row]"); hDerivPeakDerivDiff[i]->GetXaxis()->CenterTitle();
    hDerivPeakDerivDiff[i]->GetYaxis()->SetTitle("Diff in Derivative and Max Peak [row]"); hDerivPeakDerivDiff[i]->GetYaxis()->CenterTitle();
    hDerivPeakDerivDiff[i]->GetYaxis()->SetTitleOffset(1.4);
    hDerivPeakDerivDiff[i]->SetStats(false);
  }

  // Cross Section Histograms
  s1 = new TH1F("s1", "Outside Forward", 70, 0, 6);
  s1->Sumw2();
  s1->GetXaxis()->SetTitle("Center of Mass Energy [MeV]"); s1->GetXaxis()->CenterTitle();
  s1->GetYaxis()->SetTitle("Cross Section [b/sr]"); s1->GetYaxis()->CenterTitle();
  s1->GetYaxis()->SetTitleOffset(1.2);
}

inline void Spectra::InitCanvas() {
  totalCenterEnergyCanvas = 0;
  centerEnergyCanvasNum = static_cast<Int_t>(sizeof(centerEnergyCanvas)/sizeof(centerEnergyCanvas[0]));
  centerEnergyCanvasXNum = 4;
  centerEnergyCanvasYNum = 4;
  centerEnergyCanvasXYNum = centerEnergyCanvasXNum*centerEnergyCanvasYNum;
  for(Int_t i = 0; i < centerEnergyCanvasNum; i++) {
    TString name = Form("centerEnergy%d", i + 1);
    centerEnergyCanvas[i] = new TCanvas(name, name, 1600, 1200);
    centerEnergyCanvas[i]->Divide(centerEnergyCanvasXNum, centerEnergyCanvasYNum);
    centerEnergyCanvas[i]->Update();
  }

  totalCenterBeamCanvas = 0;
  centerBeamCanvasNum = static_cast<Int_t>(sizeof(centerBeamCanvas)/sizeof(centerBeamCanvas[0]));
  centerBeamCanvasXNum = 4;
  centerBeamCanvasYNum = 4;
  centerBeamCanvasXYNum = centerBeamCanvasXNum*centerBeamCanvasYNum;
  for(Int_t i = 0; i < centerBeamCanvasNum; i++) {
    TString name = Form("centerBeam%d", i + 1);
    centerBeamCanvas[i] = new TCanvas(name, name, 1600, 1200);
    centerBeamCanvas[i]->Divide(centerBeamCanvasXNum, centerBeamCanvasYNum);
    centerBeamCanvas[i]->Update();
  }

  totalCenterNoiseCanvas = 0;
  centerNoiseCanvasNum = static_cast<Int_t>(sizeof(centerNoiseCanvas)/sizeof(centerNoiseCanvas[0]));
  centerNoiseCanvasXNum = 4;
  centerNoiseCanvasYNum = 4;
  centerNoiseCanvasXYNum = centerNoiseCanvasXNum*centerNoiseCanvasYNum;
  for(Int_t i = 0; i < centerNoiseCanvasNum; i++) {
    TString name = Form("centerNoise%d", i + 1);
    centerNoiseCanvas[i] = new TCanvas(name, name, 1600, 1200);
    centerNoiseCanvas[i]->Divide(centerNoiseCanvasXNum, centerNoiseCanvasYNum);
    centerNoiseCanvas[i]->Update();
  }

  totalCenterEnergyRunningCanvas = 0;
  centerEnergyRunningCanvasNum = static_cast<Int_t>(sizeof(centerEnergyRunningCanvas)/sizeof(centerEnergyRunningCanvas[0]));
  centerEnergyRunningCanvasXNum = 4;
  centerEnergyRunningCanvasYNum = 4;
  centerEnergyRunningCanvasXYNum = centerEnergyRunningCanvasXNum*centerEnergyRunningCanvasYNum;
  for(Int_t i = 0; i < centerEnergyRunningCanvasNum; i++) {
    TString name = Form("centerEnergyRunning%d", i + 1);
    centerEnergyRunningCanvas[i] = new TCanvas(name, name, 1600, 1200);
    centerEnergyRunningCanvas[i]->Divide(centerEnergyRunningCanvasXNum, centerEnergyRunningCanvasYNum);
    centerEnergyRunningCanvas[i]->Update();
  }

  totalCenterEnergyDerivCanvas = 0;
  centerEnergyDerivCanvasNum = static_cast<Int_t>(sizeof(centerEnergyDerivCanvas)/sizeof(centerEnergyDerivCanvas[0]));
  centerEnergyDerivCanvasXNum = 4;
  centerEnergyDerivCanvasYNum = 4;
  centerEnergyDerivCanvasXYNum = centerEnergyDerivCanvasXNum*centerEnergyDerivCanvasYNum;
  for(Int_t i = 0; i < centerEnergyDerivCanvasNum; i++) {
    TString name = Form("centerEnergyDeriv%d", i + 1);
    centerEnergyDerivCanvas[i] = new TCanvas(name, name, 1600, 1200);
    centerEnergyDerivCanvas[i]->Divide(centerEnergyDerivCanvasXNum, centerEnergyDerivCanvasYNum);
    centerEnergyDerivCanvas[i]->Update();
  }
}

inline void Spectra::DrawCenterEnergyCanvas(Int_t count, std::vector<mmCenter> centerMatched_,
                                            std::vector<mmTrack> centerBeamTotal_) {
  // Make maps of different central columns
  std::map<Int_t, std::map<Int_t, Double_t> > centralEnergyMap;
  for(auto mm : centerMatched_) {
    if(mm.row > 111) continue;
    centralEnergyMap[mm.column][mm.row] = mm.energy;
  }

  TMultiGraph* mgColumn = new TMultiGraph();
  TGraph* graphColumn0 = new TGraph();
  TGraph* graphColumn1 = new TGraph();
  TGraph* graphColumn2 = new TGraph();
  TGraph* graphColumn3 = new TGraph();
  TGraph* graphColumn4 = new TGraph();
  TGraph* graphColumn5 = new TGraph();
  Int_t i = 0;
  for(auto map : centralEnergyMap[0]) {
    graphColumn0->SetPoint(i, map.first, map.second);
    i++;
  }
  i = 0;
  for(auto map : centralEnergyMap[1]) {
    graphColumn1->SetPoint(i, map.first, map.second);
    i++;
  }
  i = 0;
  for(auto map : centralEnergyMap[2]) {
    graphColumn2->SetPoint(i, map.first, map.second);
    i++;
  }
  i = 0;
  for(auto map : centralEnergyMap[3]) {
    graphColumn3->SetPoint(i, map.first, map.second);
    i++;
  }
  i = 0;
  for(auto map : centralEnergyMap[4]) {
    graphColumn4->SetPoint(i, map.first, map.second);
    i++;
  }
  i = 0;
  for(auto map : centralEnergyMap[5]) {
    graphColumn5->SetPoint(i, map.first, map.second);
    i++;
  }

  TGraph* graphColumnTot = new TGraph();
  i = 0;
  for(auto mm : centerBeamTotal_) {
    graphColumnTot->SetPoint(i, mm.row, mm.energy);
    i++;
  }

  graphColumn0->SetLineColor(28);
  graphColumn1->SetLineColor(3);
  graphColumn2->SetLineColor(4);
  graphColumn3->SetLineColor(6);
  graphColumn4->SetLineColor(7);
  graphColumn5->SetLineColor(9);
  graphColumnTot->SetLineColor(1);

  if(count < centerEnergyCanvasNum*centerEnergyCanvasXYNum) {
    TString name = Form("Event_%lld", entry);
    mgColumn->SetTitle(name);
    Int_t histNum = count/centerEnergyCanvasXYNum;
    centerEnergyCanvas[histNum]->cd(count + 1 - histNum*centerEnergyCanvasXYNum);
    if(graphColumn0->GetN() > 0) mgColumn->Add(graphColumn0);
    if(graphColumn1->GetN() > 0) mgColumn->Add(graphColumn1);
    if(graphColumn2->GetN() > 0) mgColumn->Add(graphColumn2);
    if(graphColumn3->GetN() > 0) mgColumn->Add(graphColumn3);
    if(graphColumn4->GetN() > 0) mgColumn->Add(graphColumn4);
    if(graphColumn5->GetN() > 0) mgColumn->Add(graphColumn5);
    if(graphColumnTot->GetN() > 0) mgColumn->Add(graphColumnTot);
    mgColumn->GetXaxis()->SetLimits(0, 128);
    mgColumn->SetMinimum(0);
    mgColumn->Draw("a");
    centerEnergyCanvas[histNum]->Update();
  }
}

inline void Spectra::DrawCenterBeamCanvas(Int_t count, std::vector<mmTrack> centerBeamTrack_,
    std::vector<Double_t> pars, Int_t lastRow) {
  auto* beamMG = new TMultiGraph();
  auto* beamGraph = new TGraph();
  auto* beamGraphFit = new TGraph();

  // Draw beam
  beamGraph->SetMarkerStyle(8);
  Int_t i = 0;
  for(auto mm : centerBeamTrack_) {
    beamGraph->SetPoint(i, mm.xPosition, mm.yPosition/1.75);
    i++;
  }

  // Draw fit
  beamGraphFit->SetMarkerStyle(7);
  beamGraphFit->SetMarkerColor(2);
  Int_t numPoints = 1000;
  Double_t x, y, z;
  Double_t y_begin = centerBeamTrack_[0].yPosition;
  Double_t y_end = centerBeamTrack_[centerBeamTrack_.size() - 1].yPosition;
  Double_t y_step = (y_end - y_begin)/static_cast<Double_t>(numPoints);
  for(Int_t i = 0; i < numPoints; i++) {
    line(i*y_step + y_begin, pars, x, y, z);
    if(y/1.75 > lastRow) continue;
    beamGraphFit->SetPoint(i, x, y/1.75);
  }

  // Draw boundary
  auto *boundGraph = new TGraph();
  boundGraph->SetMarkerStyle(8);
  boundGraph->SetMarkerColor(0);
  boundGraph->SetPoint(0, -20, 0);
  boundGraph->SetPoint(1, 20, 130);

  if(count < centerBeamCanvasNum*centerBeamCanvasXYNum) {
    TString name = Form("Event_%lld", entry);
    beamMG->SetTitle(name);
    Int_t histNum = count/centerBeamCanvasXYNum;
    centerBeamCanvas[histNum]->cd(count + 1 - histNum*centerBeamCanvasXYNum);
    beamMG->Add(beamGraph);
    beamMG->Add(beamGraphFit);
    beamMG->Add(boundGraph);
    beamMG->Draw("ap");
    centerBeamCanvas[histNum]->Update();
  }
}

inline void Spectra::DrawCenterNoiseCanvas(Int_t count, std::vector<mmCenter> centerReduced_,
                                           std::vector<mmCenter> centerReducedNoise_) {
  TMultiGraph* mgCenter = new TMultiGraph();
  TGraph* graphReduced = new TGraph();
  TGraph* graphReducedNoise = new TGraph();

  Int_t i = 0;
  for(auto mm : centerReduced_) {
    graphReduced->SetPoint(i, mm.column - 3, mm.row);
    i++;
  }
  i = 0;
  for(auto mm : centerReducedNoise_) {
    graphReducedNoise->SetPoint(i, mm.column - 3, mm.row);
    i++;
  }

  graphReduced->SetMarkerColor(1);
  graphReduced->SetMarkerStyle(8);
  graphReducedNoise->SetMarkerColor(2);
  graphReducedNoise->SetMarkerStyle(7);

  if(count < centerNoiseCanvasNum*centerNoiseCanvasXYNum) {
    Int_t histNum = count/centerNoiseCanvasXYNum;
    TString name = Form("Event_%lld", entry);
    mgCenter->SetTitle(name);
    centerNoiseCanvas[histNum]->cd(count + 1 - histNum*centerNoiseCanvasXYNum);
    if(graphReduced->GetN() > 0) mgCenter->Add(graphReduced);
    if(graphReducedNoise->GetN() > 0) mgCenter->Add(graphReducedNoise);
    mgCenter->GetXaxis()->SetLimits(-5, 5);
    mgCenter->SetMinimum(-1);
    mgCenter->SetMaximum(130);
    mgCenter->Draw("ap");
    centerNoiseCanvas[histNum]->Update();
  }
}

inline void Spectra::DrawCenterEnergyRunningAverageCanvas(Int_t count, std::vector<mmTrack> rawTrack_,
                                                          std::vector<mmTrack> averageTrack_) {
  TMultiGraph* mgColumn = new TMultiGraph();
  TGraph* graphRaw = new TGraph();
  TGraph* graphAverage = new TGraph();

  Int_t i = 0;
  for(auto mm : rawTrack_) {
    graphRaw->SetPoint(i, mm.row, mm.energy);
    i++;
  }

  i = 0;
  for(auto mm : averageTrack_) {
    graphAverage->SetPoint(i, mm.row, mm.energy);
    i++;
  }

  graphRaw->SetLineColor(1);
  graphAverage->SetLineColor(2);

  if(count < centerEnergyRunningCanvasNum*centerEnergyRunningCanvasXYNum) {
    Int_t histNum = count/centerEnergyRunningCanvasXYNum;
    TString name = Form("Event_%lld", entry);
    mgColumn->SetTitle(name);
    centerEnergyRunningCanvas[histNum]->cd(count + 1 - histNum*centerEnergyRunningCanvasXYNum);
    if(graphRaw->GetN() > 0) mgColumn->Add(graphRaw);
    if(graphAverage->GetN() > 0) mgColumn->Add(graphAverage);
    mgColumn->GetXaxis()->SetLimits(0, 128);
    mgColumn->SetMinimum(0);
    mgColumn->Draw("a");
    centerEnergyRunningCanvas[histNum]->Update();
  }
}

inline void Spectra::DrawCenterEnergyDeriv(Int_t count, std::vector<centerDeriv> threePoint_,
                                           std::vector<centerDeriv> fivePoint_) {
  TMultiGraph* mgColumn = new TMultiGraph();
  TGraph* graphThree = new TGraph();
  TGraph* graphFive = new TGraph();

  Int_t i = 0;
  for(auto mm : threePoint_) {
    graphThree->SetPoint(i, mm.row, mm.deriv);
    i++;
  }

  i = 0;
  for(auto mm : fivePoint_) {
    graphFive->SetPoint(i, mm.row, mm.deriv);
    i++;
  }

  graphThree->SetLineColor(2);
  graphFive->SetLineColor(4);

  if(count < centerEnergyDerivCanvasNum*centerEnergyDerivCanvasXYNum) {
    Int_t histNum = count/centerEnergyRunningCanvasXYNum;
    TString name = Form("Event_%lld", entry);
    mgColumn->SetTitle(name);
    centerEnergyDerivCanvas[histNum]->cd(count + 1 - histNum*centerEnergyDerivCanvasXYNum);
    if(graphThree->GetN() > 0) mgColumn->Add(graphThree);
    if(graphFive->GetN() > 0) mgColumn->Add(graphFive);
    mgColumn->GetXaxis()->SetLimits(0, 128);
    mgColumn->Draw("a");
    centerEnergyDerivCanvas[histNum]->Update();
  }
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
  std::ifstream inGainFile("centerGain.dat");
  assert(inGainFile.is_open());
  Int_t varI, varJ;
  Double_t varScale;
  while(inGainFile >> varI >> varJ >> varScale) {
    scale[varJ][varI] = varScale;
  }
  inGainFile.close();
}

inline void Spectra::InitAverageBeamEnergy() {
  printf("Reading Average Beam Energy File\n");
  std::ifstream inBeamFile("averageBeamEnergy.dat");
  assert(inBeamFile.is_open());
  Int_t varI, varJ;
  Double_t varE;
  while(inBeamFile >> varI >> varJ >> varE) {
    averageBeamEnergy[varJ][varI] = varE;
  }
}

inline void Spectra::InitVariables() {
  printf("Initializing Variables\n");

  // Make output file
  file = new TFile("spectra.root", "recreate");

  m1 = 8.; // AMU of projectile
  m2 = 1.; // AMU of target

  heightOffset = 70; // Time of 0 height in chamber (in ns)
  driftVelocity = 0.05674449; // in mm/ns

  rowConversion = 1.75; // Side of row in mm
  rowConversionOffset = 1.75/2.; // Put the position in middle of row

  timeResolution = 40; // Buckets to time (40 ns time buckets)

  beamEnergy = 56.; // In MeV, after havar window
  density = 0.00038175; // in g/cm3, from LISE++ (Methane at 435 torr)
  numberB8 = 174809089.;

  distanceHavarToSilicon = 544.07; // Distance from Havar to Forward Silicon in mm
  distanceHavarToMM = 268.73; // Distance from Havar to beginning of MM in mm

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

  // Y position of Forward Si Detectors
  siYPosForward = 275.34;

  // Position resolution of gas (-2000 V of 435 torr Methane)
  gasPositionResolution = 0.029364;

  // MM Column Size
  mmColumnSize[0] = std::make_pair(-10.5, -7.0);
  mmColumnSize[1] = std::make_pair(-7.0, -3.5);
  mmColumnSize[2] = std::make_pair(-3.5, 0);
  mmColumnSize[3] = std::make_pair(0, 3.5);
  mmColumnSize[4] = std::make_pair(3.5, 7.0);
  mmColumnSize[5] = std::make_pair(7.0, 10.5);
}

inline void Spectra::WriteHistograms() {
  // hIonizationChamberE->Write();
  // hIonizationChamberT->Write();

  hMicroMegasCenterCumulative->Write();
  hMicroMegasCenterCumulativePosition->Write();
  hMicroMegasCenterCumulativePositionRaw->Write();
  hMicroMegasCenterTime->Write();
  hMicroMegasCenterHeight->Write();

  // Forward Si Detectors
  // for(UInt_t i = 0; i < 10; i++) {
    // hSiEForwardDet[i]->Write();
    // hSiTForwardDet[i]->Write();
    // hSiEForwardDetCal[i]->Write();
    // for (int j = 0; j < 4; j++) {
      // hSiEForward[i][j]->Write();
      // hSiEForwardCal[i][j]->Write();
      // hSiTForward[i][j]->Write();
    // }
    // hCsIEForward[i]->Write();
    // hCsITForward[i]->Write();
  // }

  // Forward CsI Energy vs Time
  // for(UInt_t i = 0; i < 10; i++) {
    // hCsIETForward[i]->Write();
  // }

  // Forward Si Energy vs CsI Energy
  // for(UInt_t i = 0; i < 10; i++) {
    // hSiCsIEForwardDet[i]->Write();
    // hSiCsIEForwardDetCal[i]->Write();
    //  for(UInt_t j = 0; j < 4; j++) {
      //  hSiCsIEForward[i][j]->Write();
      //  hSiCsIEForwardCal[i][j]->Write();
    //  }
  // }

  // Forward dE vs Si Energy
  for(UInt_t i = 0; i < 10; i++) {
  //  hdEEForward[i]->Write();
  //  hdEEForwardCal[i]->Write();
    hdEEForwardCalTotal[i]->Write();
  }

  // Forward Hough Angle
  // for(UInt_t i = 0; i < 10; i++) {
    // hHoughAngle[i]->Write();
  // }

  // Forward Vertex vs Si Energy
  // for(UInt_t i = 0; i < 10; i++) {
    // hVertexSiEForward[i]->Write();
    // hVertexSiEForwardCal[i]->Write();
    // hVertexSiEForwardCalTotal[i]->Write();
  // }
  // hVertexSiETotalRegion3->Write();
  // hVertexCMERegion3->Write();

  // Forward Angle vs Si Energy
  // for(UInt_t i = 0; i < 10; i++) {
    // hAngleEForward[i]->Write();
    // hAngleEForwardCal[i]->Write();
    // hAngleEForwardCalTotal[i]->Write();
    // hAngleEForwardProtonEnergy[i]->Write();
    // hAngleEForwardCMEnergy[i]->Write();
  // }

  // Forward Vertex vs Angle
  // for(UInt_t i = 0; i < 10; i++) {
    // hVertexAngleForward[i]->Write();
  // }

  // Time vs Column/Strip Number Forward Detectors
  // for(UInt_t i = 0; i < 10; i++) {
    // for(UInt_t j = 0; j < 4; j++) {
      // hTimeChainForward[i][j]->Write();
      // hTimeStripForward[i][j]->Write();
    // }
  // }

  // Time vs Central Row Forward Detectors
  // for(UInt_t i = 0; i < 10; i++) {
    // hTimeCentralForward[i]->Write();
  // }

  // Forward Wall XZ Hit Positions
  // hHitPositionsXZForward->Write();
  // for(UInt_t i = 0; i < 10; i++) {
    // hHitPositionsXZForwardInd[i]->Write();
  // }

  // Max Peak Location Central Pad vs Si E Forward Wall
  // for(UInt_t i = 0; i < 10; i++) {
  for(UInt_t i = 4; i < 6; i++) {
    hMaxPeakSiE[i]->Write();
  }

  // Max Peak Location Central Pad vs Average E Forward Wall
  // for(UInt_t i = 0; i < 10; i++) {
  for(UInt_t i = 4; i < 6; i++) {
    hMaxPeakAvgPeakE[i]->Write();
  }

  // Max Peak Energy in Central Pad vs Average Peak Energy Forward Detectors
  // for(UInt_t i = 0; i < 10; i++) {
  for(UInt_t i = 4; i < 6; i++) {
    hMaxPeakEAvgPeakE[i]->Write();
  }

  // Average Peak Energy vs Si E Forward Detectors
  // for(UInt_t i = 0; i < 10; i++) {
  for(UInt_t i = 4; i < 6; i++) {
    hAvgPeakESiE[i]->Write();
  }

  // Max Peak Energy in Central Pad vs Derivative around Peak Forward Detectors
  // for(UInt_t i = 0; i < 10; i++) {
  for(UInt_t i = 4; i < 6; i++) {
    hMaxPeakEDerivPeak[i]->Write();
  }

  // Max Peak Location in Central Pad vs Difference between Peak Location and Derivative Max
  // for(UInt_t i = 0; i < 10; i++) {
  for(UInt_t i = 4; i < 6; i++) {
    hMaxPeakDerivDiff[i]->Write();
  }

  // Derivative Location in Central Pad vs Difference between Peak Location and Derivative Max
  // for(UInt_t i = 0; i < 10; i++) {
  for(UInt_t i = 4; i < 6; i++) {
    hDerivPeakDerivDiff[i]->Write();
  }

  hCWTECentral->Write();
  hCWTSiE->Write();
}

inline void Spectra::WriteCanvas() {
  for(Int_t i = 0; i < centerEnergyCanvasNum; i++) {
    centerEnergyCanvas[i]->Write();
  }

  for(Int_t i = 0; i < centerBeamCanvasNum; i++) {
    centerBeamCanvas[i]->Write();
  }

  for(Int_t i = 0; i < centerEnergyRunningCanvasNum; i++) {
    centerEnergyRunningCanvas[i]->Write();
  }

  for(Int_t i = 0; i < centerEnergyDerivCanvasNum; i++) {
    centerEnergyDerivCanvas[i]->Write();
  }
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
  outTree->Branch("siPosX", &siPosX);
  outTree->Branch("siPosZ", &siPosZ);
  return;
}

inline void Spectra::FillTree() {
  outTree->Fill();
}

inline void Spectra::WriteTree() {
  outTree->Write();
}

#endif // #ifdef Spectra_cxx
