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
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <TPolyLine3D.h>
#include <TRandom3.h>
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
void line(double t, std::vector<double> p, double &x, double &y, double &z) {
  x = p[0] + p[1]*t;
  y = t;
  z = p[2] + p[3]*t;
}

struct sortByRowMMChainStrip {
  inline bool operator() (const mmChainStrip& struct1, const mmChainStrip& struct2) {
    return (struct1.row < struct2.row);
  }
};

struct sortByRowMMCenter {
  inline bool operator() (const mmCenter& struct1, const mmCenter& struct2) {
    return (struct1.row < struct2.row);
  }
};

struct sortByRowMMTrack {
  inline bool operator() (const mmTrack& struct1, const mmTrack& struct2) {
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

  std::map<Int_t, std::pair<int, int> > siForwardMap;
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

// Cuts
  TCutG *dEEForwardCut[10];

  TCutG *angleTotEnergyCut[10];

  TCutG *timeChainForwardCut[10][4];
  TCutG *timeStripForwardCut[10][4];

  TCutG *cwtE_CentralCut;
  TCutG *cwtE_CentralProtonCut;

  TCutG *siCsiEForwardCut[10];

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
  TH1F* hMicroMegasChainLeftCumulative;
  TH1F* hMicroMegasStripLeftCumulative;

  // Micromegas - Right
  TH1F* hMicroMegasChainRightCumulative;
  TH1F* hMicroMegasStripRightCumulative;

  // Number of Si detectors fired
  TH1I* hSiFired;

  // Forward Si Detectors
  TH1F* hSiEForwardDet[10];
  TH1F* hSiEForwardDetCal[10];
  TH1F* hSiTForwardDet[10];
  TH1F* hSiEForward[10][4];
  TH1F* hSiEForwardCal[10][4];
  TH1F* hSiTForward[10][4];

  // Forward CsI Detectors
  TH1F* hCsiEForward[10];
  TH1F* hCsiEForwardCal[10];
  TH1F* hCsiTForward[10];
  TH2F* hCsiETForward[10];

  // Forward Total Energy
  TH1F* hTotalEForward[10];

  // Si vs CsI
  TH2F* hSiCsiEForwardDet[10];
  TH2F* hSiCsiEForwardDetCal[10];
  TH2F* hSiCsiEForward[10][4];
  TH2F* hSiCsiEForwardCal[10][4];

  // Si + CsI vs Si (raw)
  TH2F* hSumSiEForwardDet[10];
  TH2F* hSumSiEForward[10][4];

  // Si + CsI vs CsI (raw)
  TH2F* hSumCsiEForwardDet[10];
  TH2F* hSumCsiEForward[10][4];

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
  TH2F* hVertexSiETotalRegion1;
  TH2F* hVertexSiETotalRegion3;
  TH2F* hVertexCMERegion3;

  // Vertex vs Angle Forward Detectors
  TH2F* hVertexAngleForward[10];

  // Time vs Column Number Forward Detectors
  TH2F* hTimeChainForward[10][4];
  TH2F* hTimeChainForwardCumulative;

  // Time vs Strip Number Forward Detectors
  TH2F* hTimeStripForward[10][4];
  TH2F* hTimeStripForwardCumulative;

  // Time vs Central Region Forward Detectors
  TH2F* hTimeCentralForward[10];

  // Time vs dE Row Number Forward Detectors
  TH2F* hTimeCentraldEForward[10];

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
  TH1F* s2;
  TH1F* s3;

  void WriteHistograms();

// TCanvas

  void InitCanvas();

  // Draw Central Pads Energy
  void DrawCenterEnergyCanvas(int count, std::vector<mmCenter> centerMatched_, std::vector<mmTrack> centerBeamTotal_);
  int totalCenterEnergyCanvas;
  int centerEnergyCanvasNum, centerEnergyCanvasXYNum, centerEnergyCanvasXNum, centerEnergyCanvasYNum;
  TCanvas* centerEnergyCanvas[5];

  // Draw Beam Track
  void DrawCenterBeamCanvas(int count, std::vector<mmTrack> centerBeamTrack_, std::vector<double> pars,
      int lastRow);
  int totalCenterBeamCanvas;
  int centerBeamCanvasNum, centerBeamCanvasXYNum, centerBeamCanvasXNum, centerBeamCanvasYNum;
  TCanvas* centerBeamCanvas[5];

  void DrawCenterNoiseCanvas(int count, std::vector<mmCenter> centerReduced_,
                             std::vector<mmCenter> centerReducedNoise_);
  int totalCenterNoiseCanvas;
  int centerNoiseCanvasNum, centerNoiseCanvasXYNum, centerNoiseCanvasXNum, centerNoiseCanvasYNum;
  TCanvas* centerNoiseCanvas[10];

  // Draw Central Pad Raw Energy vs Running Average Energy
  void DrawCenterEnergyRunningAverageCanvas(int count, std::vector<mmTrack> rawTrack_,
                                             std::vector<mmTrack> averageTrack_);
  int totalCenterEnergyRunningCanvas;
  int centerEnergyRunningCanvasNum, centerEnergyRunningCanvasXYNum, centerEnergyRunningCanvasXNum, centerEnergyRunningCanvasYNum;
  TCanvas* centerEnergyRunningCanvas[5];

  // Draw Three and Five Point Derivatives with Central Region with Respect to Energy
  void DrawCenterEnergyDerivCanvas(int count, std::vector<centerDeriv> threePoint_, std::vector<centerDeriv> fivePoint_);
  int totalCenterEnergyDerivCanvas;
  int centerEnergyDerivCanvasNum, centerEnergyDerivCanvasXYNum, centerEnergyDerivCanvasXNum, centerEnergyDerivCanvasYNum;
  TCanvas* centerEnergyDerivCanvas[5];

  // Draw Linear Fit of Central Regions
  void DrawCenterBeamLinearCanvas(int count, std::vector<mmTrack> centerMatched_, std::vector<double> parsBeam, std::vector<double> parsRecoil,
      int vertexRow);
  int totalCenterBeamLinearCanvas;
  int centerBeamLinearCanvasNum, centerBeamLinearCanvasXYNum, centerBeamLinearCanvasXNum, centerBeamLinearCanvasYNum;
  TCanvas* centerBeamLinearCanvas[5];

  // Draw Event Track
  void DrawEventTrackSideCanvas(int count, std::vector<mmTrack> center_, std::vector<mmTrack> centerProton_, std::vector<mmTrack> left_, std::vector<mmTrack> leftRaw_,
      std::vector<mmTrack> right_, std::vector<mmTrack> rightRaw_, std::vector<double> parProton);
  int totalEventTrackSideCanvas;
  int eventTrackSideCanvasNum, eventTrackSideCanvasXYNum, eventTrackSideCanvasXNum, eventTrackSideCanvasYNum;
  TCanvas* eventTrackSideCanvas[5];

  // Draw dE Event Track for Side Region
  void DrawEventTrackSidedECanvas(int count, std::vector<mmTrack> leftRaw_, std::vector<mmTrack> rightRaw_, double leftAngle, double rightAngle);
  int totalEventTrackSidedECanvas;
  int eventTrackSidedECanvasNum, eventTrackSidedECanvasXYNum, eventTrackSidedECanvasXNum, eventTrackSidedECanvasYNum;
  TCanvas* eventTrackSidedECanvas[5];

  void WriteCanvas();

// Silicon Energy Calibration
private:
  void InitSiEForwardCalibration();
  std::pair<double, double> siEForwardCalibration[10][4] = {std::make_pair(0., 0.)};
  std::pair<double, double> siELeftCalibration[6][4] = {std::make_pair(0., 0.)};

  void InitCsIECalibration();
  std::pair<double, double> csiEForwardCalibration[10] = {std::make_pair(0., 0.)};

// Center Pad Gain Match
private:
  void InitCentralPadGainMatch();
  double scale[6][128];

// Beam Average Central Pads
  void InitAverageBeamEnergy();
  double averageBeamEnergy[6][128];

// Central Max Peak
  void InitMaxPeak();
  std::map<int, double> maxPeakMap;

// General Methods
  bool AnalysisForwardCentral(std::vector<mmCenter> centerMatched_, std::vector<mmTrack> centerBeamTotal_,
                              std::vector<mmTrack> centerProton_, std::map<int, double> centralPadTotalEnergy,
                              std::vector<mmChainStrip> leftChain_, std::vector<mmChainStrip> leftStrip_,
                              std::vector<mmChainStrip> rightChain_, std::vector<mmChainStrip> rightStrip_);
  bool AnalysisForwardSide(std::vector<mmCenter> centerMatched_, std::vector<mmTrack> centerBeamTotal_,
                           std::vector<mmTrack> centerProton_,
                           std::vector<mmChainStrip> leftChain_, std::vector<mmChainStrip> leftStrip_,
                           std::vector<mmChainStrip> rightChain_, std::vector<mmChainStrip> rightStrip_);
  bool AnalysisLeftSide();

// MicroMegas Functions

  // Center Functions
  std::vector<mmCenter> CenterReduceNoise(std::vector<mmCenter> center);
  void CorrectCenterEnergy(std::vector<mmTrack> &centerMatched_, std::vector<double> parsBeam, int lastRow);
  double GaussianCDF(double x, double mean, double sigma);
  void FindMaxCentralEnergy(std::vector<mmTrack> centerMatched_, int &maxEnergyRow, double &maxEnergy,
                            double &averageMaxEnergy, double &maxDeriv);
  std::vector<mmTrack> GetRunningEnergyAverageThree(std::vector<mmTrack> centerMatched_);
  std::vector<mmTrack> GetRunningEnergyAverageFive(std::vector<mmTrack> centerMatched_);
  bool CenterOnlyOneColumn(std::vector<mmTrack> centerMatched_);
  std::vector<centerDeriv> CenterEnergyThreePointDeriv(std::vector<mmTrack> centerMatched_);
  std::vector<centerDeriv> CenterEnergyFivePointDeriv(std::vector<mmTrack> centerMatched_);
  std::pair<int, int> CenterGetDerivMax(std::vector<centerDeriv> threePoint_, std::vector<centerDeriv> fivePoint_);
  std::vector<double> FitCenterPadsLinear(std::vector<mmTrack> centerMatched_, int vertexRow);

  // Side Functions
  void CleanLeftStrip(std::vector<mmChainStrip> &cleaned, std::vector<mmChainStrip> left);
  void CleanRightStrip(std::vector<mmChainStrip> &cleaned, std::vector<mmChainStrip> rightStrip);
  void ChainStripMatch(std::vector<mmTrack> &chainStripMatched, std::vector<mmTrack> &chainStripRaw,
                       std::vector<mmChainStrip> chain_,
                       std::vector<mmChainStrip> strip_, bool leftSide, double siTime);
  void ChainStripMatchingTime(std::vector<mmTrack> &chainStripMatched, std::vector<mmChainStrip> chain,
                              std::vector<mmChainStrip> strip, bool leftSide, double siTime, int timeWindow);
  size_t ChainStripTime0NumTimeBuckets(std::vector<mmTrack> matched);
  void ChainStripMatchingBoxTime0(std::vector<mmTrack> &chainStripMatched, std::vector<mmTrack> time0);
  double ChainStripSize(std::vector<mmTrack> chainStripMatched);
  void CleanSideTracks(std::vector<mmTrack> &cleaned, std::vector<mmTrack> &cleanedRaw, std::vector<mmTrack> track,
                       double angle, int &changes);

  void SideVertexFinder(std::vector<mmTrack> centerBeam_, std::vector<mmTrack> protonTrack_, double &vertexPositionX,
                        double &vertexPositionY, double &vertexPositionZ, std::vector<double> &parsProton);
  void SideVertexFinderSingleHelp(std::vector<mmTrack> centerBeam_, std::vector<mmTrack> protonTrack_, double &vertexPositionX,
                        double &vertexPositionY, double &vertexPositionZ, std::vector<double> &parsProton);
  void SideVertexFinderHelp(std::vector<mmTrack> centerBeam_, std::vector<mmTrack> protonTrack_, double &vertexPositionX,
                        double &vertexPositionY, double &vertexPositionZ, std::vector<double> &parsProton);

  // Visualize Hough Transform
  void GetMinMaxD(std::vector<mmTrack> initPoints, int &minXY, int &maxXY, int &minYZ, int &maxYZ);
  void VisualizeHough(std::vector<mmTrack> initPoints, TH2I* fXY, TH2I* fYZ);
  void GetMinMaxDRestricted(std::vector<mmTrack> initPoints, int &minXY, int &maxXY, int &minYZ, int &maxYZ, int siDet);
  void VisualizeHoughRestricted(std::vector<mmTrack> initPoints, TH2I* fXY, TH2I* fYZ, int siDet);
  void GetHoughStdDevXYRestricted(std::vector<mmTrack> initPoints, std::vector<double> &angle_, std::vector<double> &stdDev_, int siDet);

// Cross Section
private:
  void DivideTargetThickness(TH1F *f);
  void ReadSolidAngle();
  void SolidAngle(TH1F *f, int region);
  void WriteSpectrumToFile(TH1F *f, int region);
  CubicSpline reg1SA;
  CubicSpline reg1CMAngle;
  CubicSpline reg2SA;
  CubicSpline reg2CMAngle;
  CubicSpline reg3SA;
  CubicSpline reg3CMAngle;

// General Variables
private:
  void InitVariables();
  TFile *file;
  double m1;
  double m2;
  double rowConversion;
  double rowConversionOffset;
  double heightOffset;
  double driftVelocity;
  double timeResolution;
  long entry;
  double beamEnergy;
  double density;
  double distanceHavarToSilicon;
  double distanceHavarToMM;
  double numberB8;
  double siXPosForward[10][4];
  double siYPosForward;
  double gasPositionResolution;
  std::pair<double, double> mmColumnSize[6];
  double padError;
  double padErrorTwoColumns;
  std::vector<int> leftChainHits;
  std::vector<int> leftStripHits;
  std::vector<int> rightChainHits;
  std::vector<int> rightStripHits;

  EnergyLoss *boronMethane;
  EnergyLoss *protonMethane;

// Tree Variables
  void InitTree();
  void FillTree();
  void WriteTree();
  TTree *outTree;
  int siDet;
  int siQuad;
  int siChannel;
  double siEnergy;
  double siEnergyCal;
  double siTime;
  double csiEnergy;
  double csiEnergyCal;
  double csiTime;
  double totalEnergy;
  bool punchthrough;
  double dE;
  double vertexPositionX;
  double vertexPositionY;
  double vertexPositionZ;
  double angle;
  double cmEnergy;
  double siPosX;
  double siPosY;
  double siPosZ;
  TRandom3 *rng;

};
#endif

// #ifdef Spectra_cxx

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
  int j = 0;
  for(int i = 0; i < 68; i++) {
    if(i == 11 || i == 22 || i == 45 || i == 56) continue; // FPN Channels
    Aget_Map[j] = i;
    j++;
  }

  // Asad0_Aget0
  j = 31;
  for(int i = 0; i < 64; i++) {
    if(i % 4 == 0 && i != 0) j--;
    int remain = i % 4;
    MM_Map_Asad0_Aget0[Aget_Map[i]] = std::make_pair(4 - remain, j);
  }

  // Asad0_Aget1
  j = 63;
  for(int i = 0; i < 64; i++) {
    if(i % 4 == 0 && i != 0) j--;
    int remain = i % 4;
    MM_Map_Asad0_Aget1[Aget_Map[i]] = std::make_pair(4 - remain, j);
  }

  // Asad0_Aget2
  j = 95;
  for(int i = 0; i < 64; i++) {
    if(i % 4 == 0 && i != 0) j--;
    int remain = i % 4;
    MM_Map_Asad0_Aget2[Aget_Map[i]] = std::make_pair(1 + remain, j);
  }

  // Asad0_Aget3
  j = 127;
  for(int i = 0; i < 64; i++) {
    if(i % 4 == 0 && i != 0) j--;
    int remain = i % 4;
    MM_Map_Asad0_Aget3[Aget_Map[i]] = std::make_pair(1 + remain, j);
  }

  // Asad1_Aget0
  j = 15;
  for(int i = 0; i < 64; i++) {
    if(i % 4 == 0 && i != 0) j--;
    int remain = i % 4;
    MM_Map_Asad1_Aget0[Aget_Map[i]] = std::make_pair(1 + remain, j);
  }

  // Asad1_Aget1
  j = 47;
  for(int i = 0; i < 64; i++) {
    if(i % 4 == 0 && i != 0) j--;
    int remain = i % 4;
    MM_Map_Asad1_Aget1[Aget_Map[i]] = std::make_pair(1 + remain, j);
  }

  // Asad1_Aget2
  j = 79;
  for(int i = 0; i < 64; i++) {
    if(i % 4 == 0 && i != 0) j--;
    int remain = i % 4;
    MM_Map_Asad1_Aget2[Aget_Map[i]] = std::make_pair(4 - remain, j);
  }

  // Asad1_Aget3
  j = 111;
  for(int i = 0; i < 64; i++) {
    if(i % 4 == 0 && i != 0) j--;
    int remain = i % 4;
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

  // Asad2_Aget2
  j = 63;
  for(int i = 0; i < 64; i++) {
    MM_Map_Asad2_Aget2[Aget_Map[i]] = std::make_pair(0, j);
    j--;
  }

  // Asad2_Aget3
  j = 127;
  for(int i = 0; i < 64; i++) {
    MM_Map_Asad2_Aget3[Aget_Map[i]] = std::make_pair(0, j);
    j--;
  }

  // Asad3_Aget0 (Beam right)
  // Strips
  j = 31;
  for(int i = 0; i < 32; i++) {
    MM_Map_Asad3_Aget0[Aget_Map[i]] = j;
    j--;
  }
  // Chains
  j = 32;
  for(int i = 32; i < 64; i++) {
    MM_Map_Asad3_Aget0[Aget_Map[i]] = j;
    j++;
  }

  // Asad3_Aget1 (Beam right)
  // Strips
  j = 63;
  for(int i = 0; i < 32; i++) {
    MM_Map_Asad3_Aget1[Aget_Map[i]] = j;
    j--;
  }
  // Chains
  j = 0;
  for(int i = 32; i < 64; i++) {
    MM_Map_Asad3_Aget1[Aget_Map[i]] = j;
    j++;
  }

  // Asad3_Aget2
  j = 63;
  for(int i = 0; i < 64; i++) {
    MM_Map_Asad3_Aget2[Aget_Map[i]] = std::make_pair(5, j);
    j--;
  }

  // Asad3_Aget3
  j = 127;
  for(int i = 0; i < 64; i++) {
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
  hMicroMegasCenterCumulative = new TH2F("MM_Center_Cumulative", "MM_Center_Cumulative", 20, -10, 10, 128, 0, 128);
  hMicroMegasCenterCumulativePosition = new TH2F("MM_Center_Cumulative_Position", "MM_Center_Cumulative_Position", 24, -21, 21, 143, 0, 250.25);
  hMicroMegasCenterCumulativePositionRaw = new TH2F("MM_Center_Cumulative_Position_Raw", "MM_Center_Cumulative_Position_Raw", 24, -21, 21, 143, 0, 250.25);
  hMicroMegasCenterEnergyCumulative = new TH2F("MM_Center_Energy_Cumulative", "MM_Center_Energy_Cumulative", 20, -10, 10, 128, 0, 128);
  hMicroMegasCenterEnergyAverage = new TH2F("MM_Center_Energy_Average", "MM_Center_Energy_Average", 20, -10, 10, 128, 0, 128);
  hMicroMegasCenterEnergyAverageScaled = new TH2F("MM_Center_Energy_Average_Scaled", "MM_Center_Energy_Average_Scaled", 20, -10, 10, 128, 0, 128);
  hMicroMegasCenterTime = new TH2F("MM_Center_Time", "MM_Center_Time", 130, -1, 129, 250, 0, 10000);
  hMicroMegasCenterHeight = new TH2F("MM_Center_Height", "MM_Center_Height", 130, -1, 129, 160, -200, 200);
}

inline void Spectra::InitCanvas() {
  totalCenterEnergyCanvas = 0;
  centerEnergyCanvasNum = static_cast<int>(sizeof(centerEnergyCanvas)/sizeof(centerEnergyCanvas[0]));
  centerEnergyCanvasXNum = 4;
  centerEnergyCanvasYNum = 4;
  centerEnergyCanvasXYNum = centerEnergyCanvasXNum*centerEnergyCanvasYNum;
  for(int i = 0; i < centerEnergyCanvasNum; i++) {
    TString name = Form("centerEnergy%d", i + 1);
    centerEnergyCanvas[i] = new TCanvas(name, name, 1600, 1200);
    centerEnergyCanvas[i]->Divide(centerEnergyCanvasXNum, centerEnergyCanvasYNum);
    centerEnergyCanvas[i]->Update();
  }

  totalCenterBeamCanvas = 0;
  centerBeamCanvasNum = static_cast<int>(sizeof(centerBeamCanvas)/sizeof(centerBeamCanvas[0]));
  centerBeamCanvasXNum = 4;
  centerBeamCanvasYNum = 4;
  centerBeamCanvasXYNum = centerBeamCanvasXNum*centerBeamCanvasYNum;
  for(int i = 0; i < centerBeamCanvasNum; i++) {
    TString name = Form("centerBeam%d", i + 1);
    centerBeamCanvas[i] = new TCanvas(name, name, 1600, 1200);
    centerBeamCanvas[i]->Divide(centerBeamCanvasXNum, centerBeamCanvasYNum);
    centerBeamCanvas[i]->Update();
  }

  totalCenterNoiseCanvas = 0;
  centerNoiseCanvasNum = static_cast<int>(sizeof(centerNoiseCanvas)/sizeof(centerNoiseCanvas[0]));
  centerNoiseCanvasXNum = 4;
  centerNoiseCanvasYNum = 4;
  centerNoiseCanvasXYNum = centerNoiseCanvasXNum*centerNoiseCanvasYNum;
  for(int i = 0; i < centerNoiseCanvasNum; i++) {
    TString name = Form("centerNoise%d", i + 1);
    centerNoiseCanvas[i] = new TCanvas(name, name, 1600, 1200);
    centerNoiseCanvas[i]->Divide(centerNoiseCanvasXNum, centerNoiseCanvasYNum);
    centerNoiseCanvas[i]->Update();
  }

  totalCenterEnergyRunningCanvas = 0;
  centerEnergyRunningCanvasNum = static_cast<int>(sizeof(centerEnergyRunningCanvas)/sizeof(centerEnergyRunningCanvas[0]));
  centerEnergyRunningCanvasXNum = 4;
  centerEnergyRunningCanvasYNum = 4;
  centerEnergyRunningCanvasXYNum = centerEnergyRunningCanvasXNum*centerEnergyRunningCanvasYNum;
  for(int i = 0; i < centerEnergyRunningCanvasNum; i++) {
    TString name = Form("centerEnergyRunning%d", i + 1);
    centerEnergyRunningCanvas[i] = new TCanvas(name, name, 1600, 1200);
    centerEnergyRunningCanvas[i]->Divide(centerEnergyRunningCanvasXNum, centerEnergyRunningCanvasYNum);
    centerEnergyRunningCanvas[i]->Update();
  }

  totalCenterEnergyDerivCanvas = 0;
  centerEnergyDerivCanvasNum = static_cast<int>(sizeof(centerEnergyDerivCanvas)/sizeof(centerEnergyDerivCanvas[0]));
  centerEnergyDerivCanvasXNum = 4;
  centerEnergyDerivCanvasYNum = 4;
  centerEnergyDerivCanvasXYNum = centerEnergyDerivCanvasXNum*centerEnergyDerivCanvasYNum;
  for(int i = 0; i < centerEnergyDerivCanvasNum; i++) {
    TString name = Form("centerEnergyDeriv%d", i + 1);
    centerEnergyDerivCanvas[i] = new TCanvas(name, name, 1600, 1200);
    centerEnergyDerivCanvas[i]->Divide(centerEnergyDerivCanvasXNum, centerEnergyDerivCanvasYNum);
    centerEnergyDerivCanvas[i]->Update();
  }

  totalCenterBeamLinearCanvas = 0;
  centerBeamLinearCanvasNum = static_cast<int>(sizeof(centerBeamLinearCanvas)/sizeof(centerBeamLinearCanvas[0]));
  centerBeamLinearCanvasXNum = 4;
  centerBeamLinearCanvasYNum = 4;
  centerBeamLinearCanvasXYNum = centerBeamLinearCanvasXNum*centerBeamLinearCanvasYNum;
  for(int i = 0; i < centerBeamLinearCanvasNum; i++) {
    TString name = Form("centerBeamLinear%d", i + 1);
    centerBeamLinearCanvas[i] = new TCanvas(name, name, 1600, 1200);
    centerBeamLinearCanvas[i]->Divide(centerBeamLinearCanvasXNum, centerBeamLinearCanvasYNum);
    centerBeamLinearCanvas[i]->Update();
  }

  totalEventTrackSideCanvas = 0;
  eventTrackSideCanvasNum = static_cast<int>(sizeof(eventTrackSideCanvas)/sizeof(eventTrackSideCanvas[0]));
  eventTrackSideCanvasXNum = 3;
  eventTrackSideCanvasYNum = 3;
  eventTrackSideCanvasXYNum = eventTrackSideCanvasXNum*eventTrackSideCanvasYNum;
  for(int i = 0; i < eventTrackSideCanvasNum; i++) {
    TString name = Form("eventTrackSide%d", i + 1);
    eventTrackSideCanvas[i] = new TCanvas(name, name, 1600, 1200);
    eventTrackSideCanvas[i]->Divide(eventTrackSideCanvasXNum, eventTrackSideCanvasYNum);
    eventTrackSideCanvas[i]->Update();
  }

  totalEventTrackSidedECanvas = 0;
  eventTrackSidedECanvasNum = static_cast<int>(sizeof(eventTrackSidedECanvas)/sizeof(eventTrackSidedECanvas[0]));
  eventTrackSidedECanvasXNum = 3;
  eventTrackSidedECanvasYNum = 3;
  eventTrackSidedECanvasXYNum = eventTrackSidedECanvasXNum*eventTrackSidedECanvasYNum;
  for(int i = 0; i < eventTrackSidedECanvasNum; i++) {
    TString name = Form("eventTrackSidedE%d", i + 1);
    eventTrackSidedECanvas[i] = new TCanvas(name, name, 1600, 1200);
    eventTrackSidedECanvas[i]->Divide(eventTrackSidedECanvasXNum, eventTrackSidedECanvasYNum);
    eventTrackSidedECanvas[i]->Update();
  }
}

inline void Spectra::DrawCenterEnergyCanvas(int count, std::vector<mmCenter> centerMatched_,
                                            std::vector<mmTrack> centerBeamTotal_) {
  // Make maps of different central columns
  std::map<int, std::map<int, double> > centralEnergyMap;
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
  int i = 0;
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
    TString name = Form("Event_%ld", entry);
    mgColumn->SetTitle(name);
    int histNum = count/centerEnergyCanvasXYNum;
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

inline void Spectra::DrawCenterBeamCanvas(int count, std::vector<mmTrack> centerBeamTrack_,
    std::vector<double> pars, int lastRow) {
  auto* beamMG = new TMultiGraph();
  auto* beamGraph = new TGraph();
  auto* beamGraphFit = new TGraph();

  // Draw beam
  beamGraph->SetMarkerStyle(8);
  int i = 0;
  for(auto mm : centerBeamTrack_) {
    beamGraph->SetPoint(i, mm.xPosition, mm.yPosition/1.75);
    i++;
  }

  // Draw fit
  beamGraphFit->SetMarkerStyle(7);
  beamGraphFit->SetMarkerColor(2);
  int numPoints = 1000;
  double x, y, z;
  double y_begin = centerBeamTrack_[0].yPosition;
  double y_end = centerBeamTrack_[centerBeamTrack_.size() - 1].yPosition;
  double y_step = (y_end - y_begin)/static_cast<double>(numPoints);
  for(int i = 0; i < numPoints; i++) {
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
    TString name = Form("Event_%ld", entry);
    beamMG->SetTitle(name);
    int histNum = count/centerBeamCanvasXYNum;
    centerBeamCanvas[histNum]->cd(count + 1 - histNum*centerBeamCanvasXYNum);
    beamMG->Add(beamGraph);
    beamMG->Add(beamGraphFit);
    beamMG->Add(boundGraph);
    beamMG->Draw("ap");
    centerBeamCanvas[histNum]->Update();
  }
}

inline void Spectra::DrawCenterNoiseCanvas(int count, std::vector<mmCenter> centerReduced_,
                                           std::vector<mmCenter> centerReducedNoise_) {
  TMultiGraph* mgCenter = new TMultiGraph();
  TGraph* graphReduced = new TGraph();
  TGraph* graphReducedNoise = new TGraph();

  int i = 0;
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
    int histNum = count/centerNoiseCanvasXYNum;
    TString name = Form("Event_%ld", entry);
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

inline void Spectra::DrawCenterEnergyRunningAverageCanvas(int count, std::vector<mmTrack> rawTrack_,
                                                          std::vector<mmTrack> averageTrack_) {
  TMultiGraph* mgColumn = new TMultiGraph();
  TGraph* graphRaw = new TGraph();
  TGraph* graphAverage = new TGraph();

  int i = 0;
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
    int histNum = count/centerEnergyRunningCanvasXYNum;
    TString name = Form("Event_%ld", entry);
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

inline void Spectra::DrawCenterEnergyDerivCanvas(int count, std::vector<centerDeriv> threePoint_,
                                           std::vector<centerDeriv> fivePoint_) {
  TMultiGraph* mgColumn = new TMultiGraph();
  TGraph* graphThree = new TGraph();
  TGraph* graphFive = new TGraph();

  int i = 0;
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
    int histNum = count/centerEnergyRunningCanvasXYNum;
    TString name = Form("Event_%ld", entry);
    mgColumn->SetTitle(name);
    centerEnergyDerivCanvas[histNum]->cd(count + 1 - histNum*centerEnergyDerivCanvasXYNum);
    if(graphThree->GetN() > 0) mgColumn->Add(graphThree);
    if(graphFive->GetN() > 0) mgColumn->Add(graphFive);
    mgColumn->GetXaxis()->SetLimits(0, 128);
    mgColumn->Draw("a");
    centerEnergyDerivCanvas[histNum]->Update();
  }
}

inline void Spectra::DrawCenterBeamLinearCanvas(int count, std::vector<mmTrack> center_, std::vector<double> parsBeam,
                                                std::vector<double> parsRecoil, int vertexRow) {
  TMultiGraph* mgCenter = new TMultiGraph();
  TGraphErrors* graphCenter = new TGraphErrors();
  TGraph* graphBeam = new TGraph();
  TGraph* graphRecoil = new TGraph();

  int i = 0;
  for(auto mm : center_) {
    if(mm.total == 1) {
      graphCenter->SetPoint(i, mm.row, mm.xPosition);
      graphCenter->SetPointError(i, 0, padError);
    }
    else if(mm.total == 2) {
      graphCenter->SetPoint(i, mm.row, mm.xPosition);
      graphCenter->SetPointError(i, 0, padErrorTwoColumns);
    }
    i++;
  }

  if(count < centerBeamLinearCanvasNum*centerBeamLinearCanvasXYNum) {
    int histNum = count/centerBeamCanvasXYNum;
    TString name = Form("Event_%ld", entry);
    mgCenter->SetTitle(name);
    centerBeamLinearCanvas[histNum]->cd(count + 1 - histNum*centerBeamLinearCanvasXYNum);
    if(graphCenter->GetN() > 0) mgCenter->Add(graphCenter);
    mgCenter->GetXaxis()->SetLimits(0, 128);
    mgCenter->SetMinimum(-15);
    mgCenter->SetMaximum(15);
    mgCenter->Draw("a");
    centerBeamLinearCanvas[histNum]->Update();
  }
}

inline void Spectra::DrawEventTrackSideCanvas(int count, std::vector<mmTrack> center_, std::vector<mmTrack> centerProton_,
    std::vector<mmTrack> proton_, std::vector<mmTrack> protonRaw_, std::vector<mmTrack> other_, std::vector<mmTrack> otherRaw_,
                                              std::vector<double> parProton) {
  TMultiGraph* mgEvent = new TMultiGraph();
  TGraph* graphCenter = new TGraph();
  TGraph* graphCenterProton = new TGraph();
  TGraph* graphProton = new TGraph();
  TGraph* graphProtonRaw = new TGraph();
  TGraph* graphOther = new TGraph();
  TGraph* graphOtherRaw = new TGraph();

  TGraph* graphProtonFit = new TGraph();

  int i = 0;
  for(auto mm : center_) {
    graphCenter->SetPoint(i, mm.xPosition, mm.yPosition);
    i++;
  }

  i = 0;
  for(auto mm : centerProton_) {
    graphCenterProton->SetPoint(i, mm.xPosition, mm.yPosition);
    i++;
  }

  i = 0;
  for(auto mm : proton_) {
    graphProton->SetPoint(i, mm.xPosition, mm.yPosition);
    i++;
  }

  i = 0;
  for(auto mm : protonRaw_) {
    graphProtonRaw->SetPoint(i, mm.xPosition, mm.yPosition);
    i++;
  }

  i = 0;
  for(auto mm : other_) {
    graphOther->SetPoint(i, mm.xPosition, mm.yPosition);
    i++;
  }

  i = 0;
  for(auto mm : otherRaw_) {
    graphOtherRaw->SetPoint(i, mm.xPosition, mm.yPosition);
    i++;
  }

  i = 0;
  double x, y, z;
  for(double point = 0; point < 250; point += 0.1) {
    line(point, parProton, x, y, z);
    graphProtonFit->SetPoint(i, x, y);
    i++;
  }

  graphCenter->SetMarkerColor(1);
  graphCenter->SetMarkerStyle(8);

  graphCenterProton->SetMarkerColor(2);
  graphCenterProton->SetMarkerStyle(8);

  graphProton->SetMarkerColor(2);
  graphProton->SetMarkerStyle(8);
  graphProtonRaw->SetMarkerColor(3);
  graphProtonRaw->SetMarkerStyle(7);

  graphOther->SetMarkerColor(2);
  graphOther->SetMarkerStyle(8);
  graphOtherRaw->SetMarkerColor(4);
  graphOtherRaw->SetMarkerStyle(7);

  graphProtonFit->SetMarkerColor(6);
  graphProtonFit->SetMarkerStyle(7);

  if(count < eventTrackSideCanvasNum*eventTrackSideCanvasXYNum) {
    int histNum = count/eventTrackSideCanvasXYNum;
    TString name = Form("Event_%ld", entry);
    mgEvent->SetTitle(name);
    eventTrackSideCanvas[histNum]->cd(count + 1 - histNum*eventTrackSideCanvasXYNum);
    if(graphCenter->GetN() > 0) mgEvent->Add(graphCenter);
    if(graphCenterProton->GetN() > 0) mgEvent->Add(graphCenterProton);
    if(graphProton->GetN() > 0) mgEvent->Add(graphProton);
    if(graphProtonRaw->GetN() > 0) mgEvent->Add(graphProtonRaw);
    if(graphOther->GetN() > 0) mgEvent->Add(graphOther);
    if(graphOtherRaw->GetN() > 0) mgEvent->Add(graphOtherRaw);
    if(graphProtonFit->GetN() > 0) mgEvent->Add(graphProtonFit);
    mgEvent->GetXaxis()->SetLimits(-200, 200);
    mgEvent->SetMinimum(0);
    mgEvent->SetMaximum(300);
    mgEvent->Draw("ap");
    eventTrackSideCanvas[histNum]->Update();
  }
}

inline void Spectra::DrawEventTrackSidedECanvas(int count, std::vector<mmTrack> protonRaw_, std::vector<mmTrack> otherRaw_,
                                                double protonAngle, double otherAngle) {
  TMultiGraph* mgEvent = new TMultiGraph();
  TGraph* graphProtonRaw = new TGraph();
  TGraph* graphOtherRaw = new TGraph();

  int i = 0;
  for(auto mm : protonRaw_) {
    graphProtonRaw->SetPoint(i, mm.yPosition, mm.energy*cos(protonAngle));
    i++;
  }

  i = 0;
  for(auto mm : otherRaw_) {
    graphOtherRaw->SetPoint(i, mm.yPosition, mm.energy*cos(otherAngle));
    i++;
  }

  graphProtonRaw->SetLineColor(3);
  graphOtherRaw->SetLineColor(4);

  if(count < eventTrackSidedECanvasNum*eventTrackSideCanvasXYNum) {
    int histNum = count/eventTrackSidedECanvasXYNum;
    TString name = Form("Event_%ld", entry);
    mgEvent->SetTitle(name);
    eventTrackSidedECanvas[histNum]->cd(count + 1 - histNum*eventTrackSidedECanvasXYNum);
    if(graphProtonRaw->GetN() > 0) mgEvent->Add(graphProtonRaw);
    if(graphOtherRaw->GetN() > 0) mgEvent->Add(graphOtherRaw);
    mgEvent->GetXaxis()->SetLimits(0, 300);
    mgEvent->SetMinimum(0);
    mgEvent->SetMaximum(4000);
    mgEvent->Draw("a");
    eventTrackSidedECanvas[histNum]->Update();
  }
}

inline void Spectra::InitSiEForwardCalibration() {
  printf("Reading Si Energy Forward Calibrations File\n");
  std::ifstream inSiCalFile("siCalibration.dat");
  assert(inSiCalFile.is_open());
  int var1, var2, var3;
  double slope, intercept;
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
  int var1, var2;
  double slope, intercept;
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
  int varI, varJ;
  double varScale;
  while(inGainFile >> varI >> varJ >> varScale) {
    scale[varJ][varI] = varScale;
  }
  inGainFile.close();
}

inline void Spectra::InitAverageBeamEnergy() {
  printf("Reading Average Beam Energy File\n");
  std::ifstream inBeamFile("averageBeamEnergy.dat");
  assert(inBeamFile.is_open());
  int varI, varJ;
  double varE;
  while(inBeamFile >> varI >> varJ >> varE) {
    averageBeamEnergy[varJ][varI] = varE;
  }
}

inline void Spectra::InitMaxPeak() {
  std::ifstream inMaxPeakFile("centralPeakVertex.out");
  assert(inMaxPeakFile.is_open());
  double var1, var2, var3;
  int var4;
  while(inMaxPeakFile >> var1 >> var2 >> var3 >> var4) {
    maxPeakMap[var4] = var3;

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

  timeResolution = 40.; // Buckets to time (40 ns time buckets)

  beamEnergy = 56.; // In MeV, after havar window
  density = 0.00038175; // in g/cm3, from LISE++ (Methane at 435 torr)
  numberB8 = 174809089.;

  distanceHavarToSilicon = 544.; // Distance from Havar to Forward Silicon in mm
  distanceHavarToMM = 270.; // Distance from Havar to beginning of MM in mm

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
  siYPosForward = 274.;

  // Position resolution of gas (-2000 V of 435 torr Methane)
  gasPositionResolution = 0.029364;

  // MM Column Size
  mmColumnSize[0] = std::make_pair(-10.5, -7.0);
  mmColumnSize[1] = std::make_pair(-7.0, -3.5);
  mmColumnSize[2] = std::make_pair(-3.5, 0);
  mmColumnSize[3] = std::make_pair(0, 3.5);
  mmColumnSize[4] = std::make_pair(3.5, 7.0);
  mmColumnSize[5] = std::make_pair(7.0, 10.5);

  // Pad Error for Fitting of Central Pads
  padError = 1.75; // in mm
  padErrorTwoColumns = 0.5; // in mm

  // Read in hits on leftChain, leftStrip, rightChain, rightStrip
  std::ifstream inLeftChain("leftChainHits.out");
  int leftChain;
  while(inLeftChain >> leftChain) {
    leftChainHits.push_back(leftChain);
  }
  inLeftChain.close();

  std::ifstream inLeftStrip("leftStripHits.out");
  int leftStrip;
  while(inLeftStrip >> leftStrip) {
    leftStripHits.push_back(leftStrip);
  }
  inLeftStrip.close();

  std::ifstream inRightChain("rightChainHits.out");
  int rightChain;
  while(inRightChain >> rightChain) {
    rightChainHits.push_back(rightChain);
  }
  inRightChain.close();

  std::ifstream inRightStrip("rightStripHits.out");
  int rightStrip;
  while(inRightStrip >> rightStrip) {
    rightStripHits.push_back(rightStrip);
  }
  inRightStrip.close();

  rng = new TRandom3();
}

inline void Spectra::WriteHistograms() {
  hIonizationChamberE->Write();
  hIonizationChamberT->Write();

  hMicroMegasCenterCumulative->Write();
  hMicroMegasCenterCumulativePosition->Write();
  hMicroMegasCenterCumulativePositionRaw->Write();
  hMicroMegasCenterTime->Write();
  hMicroMegasCenterEnergyAverageScaled->Write();

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

// #endif // #ifdef Spectra_cxx
