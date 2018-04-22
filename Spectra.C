#define Spectra_cxx
#include "Spectra.h"

TChain* MakeChain();

int main() {

  TChain *chain = MakeChain();

  Spectra t(chain);
  t.Loop();

  return 0;
}

TChain* MakeChain() {
  TChain *chain = new TChain("mfmData");

  // Home
  // TString PathToFiles = "/hd3/research/data/run0817a/rootM2R-MaxModule/run";

  // Laptop
  TString PathToFiles = "/Users/joshhooker/Desktop/data/run0817a/run";

  // chain->Add(PathToFiles+"004.root");
  chain->Add(PathToFiles+"175.root");
  // chain->Add(PathToFiles+"178.root");
  // chain->Add(PathToFiles+"180.root");
  // chain->Add(PathToFiles+"181.root");
  // chain->Add(PathToFiles+"182.root");
  // chain->Add(PathToFiles+"183.root");
  // chain->Add(PathToFiles+"184.root");
  // chain->Add(PathToFiles+"185.root");
  // chain->Add(PathToFiles+"186.root");
  // chain->Add(PathToFiles+"187.root");
  // chain->Add(PathToFiles+"188.root");
  // chain->Add(PathToFiles+"189.root");
  // chain->Add(PathToFiles+"190.root");
  // chain->Add(PathToFiles+"191.root");
  // chain->Add(PathToFiles+"192.root");
  // chain->Add(PathToFiles+"193.root");
  // chain->Add(PathToFiles+"195.root");
  // chain->Add(PathToFiles+"196.root");
  // chain->Add(PathToFiles+"197.root");
  // chain->Add(PathToFiles+"198.root");
  // chain->Add(PathToFiles+"199.root");
  // chain->Add(PathToFiles+"200.root");
  // chain->Add(PathToFiles+"201.root");
  // chain->Add(PathToFiles+"202.root");
  // chain->Add(PathToFiles+"203.root");
  // chain->Add(PathToFiles+"204.root");
  // chain->Add(PathToFiles+"205.root");
  // chain->Add(PathToFiles+"206.root");
  // chain->Add(PathToFiles+"207.root");
  // chain->Add(PathToFiles+"208.root");
  // chain->Add(PathToFiles+"210.root");
  // chain->Add(PathToFiles+"211.root");
  // chain->Add(PathToFiles+"212.root");
  // chain->Add(PathToFiles+"213.root");
  // chain->Add(PathToFiles+"214.root");
  // chain->Add(PathToFiles+"215.root");
  // chain->Add(PathToFiles+"216.root");
  // chain->Add(PathToFiles+"217.root");
  // chain->Add(PathToFiles+"218.root");
  // chain->Add(PathToFiles+"219.root");
  // chain->Add(PathToFiles+"220.root");
  // chain->Add(PathToFiles+"221.root");
  // chain->Add(PathToFiles+"223.root");
  // chain->Add(PathToFiles+"224.root");
  // chain->Add(PathToFiles+"232.root");

  return chain;
}

void Spectra::Loop() {

  if (fChain == 0) return;

  InitChannelMap();
  InitHistograms();
  InitVariables();

  InitSiEForwardCalibration();

  InitCentralPadGainMatch();
  InitAverageBeamEnergy();

  Bool_t individualMMHistograms = false;

  Long64_t nentries = fChain->GetEntriesFast();

  /////////////////
  // Set up cuts //
  /////////////////
  // TFile* cutFile = TFile::Open("cuts.root");
  // TCutG* det5_dEE_noPunchthroughCut = (TCutG*)cutFile->Get("det5_dEECut_noPunchthrough");
  // TCutG* det56_dEECut = (TCutG*)cutFile->Get("det56_dEECut");
  // cutFile->Close();

  ////////////////
  // Histograms //
  ////////////////

  TH2F* hMicroMegasCenterCumulative = new TH2F("MM_Center_Cumulative", "MM_Center_Cumulative", 82, -41, 41, 128, 0, 128);
  TH2F* hMicroMegasCenterEnergyCumulative = new TH2F("MM_Center_Energy_Cumulative", "MM_Center_Energy_Cumulative", 20, -10, 10, 128, 0, 128);
  TH2F* hMicroMegasCenterEnergyAverage = new TH2F("MM_Center_Energy_Average", "MM_Center_Energy_Average", 20, -10, 10, 128, 0, 128);
  TH2F* hMicroMegasCenterEnergyAverageScaled = new TH2F("MM_Center_Energy_Average_Scaled", "MM_Center_Energy_Average_Scaled", 20, -10, 10, 128, 0, 128);
  TH2F* hMicroMegasCenterTimeAverage = new TH2F("MM_Center_Time_Average", "MM_Center_Time_Average", 20, -10, 10, 128, 0, 128);

  TH1F* hMicroMegasStripLeftCumulative = new TH1F("MM_Strip_Left_Cumulative", "MM_Strip_Left_Cumulative", 64, 0, 64);
  TH1F* hMicroMegasStripRightCumulative = new TH1F("MM_Strip_Right_Cumulative", "MM_Strip_Right_Cumulative", 64, 0, 64);

  TH1F* hMicroMegasChainLeftCumulative = new TH1F("MM_Chain_Left_Cumulative", "MM_Chain_Left_Cumulative", 64, 0, 64);
  TH1F* hMicroMegasChainRightCumulative = new TH1F("MM_Chain_Right_Cumulative", "MM_Chain_Right_Cumulative", 64, 0, 64);

  TH1F* hIonizationChamberE = new TH1F("icE", "icE", 500, 0, 4000);
  TH1F* hIonizationChamberT = new TH1F("icT", "icT", 500, 0, 20000);

  TH1F* hRow116TimeDet5 = new TH1F("row116_time_det5", "row116_time_det5", 100, 0, 10000);
  TH1F* hRow116TimeDet6 = new TH1F("row116_time_det6", "row116_time_det6", 100, 0, 10000);
  TH1F* hRow123TimeDet5 = new TH1F("row123_time_det5", "row123_time_det5", 100, 0, 10000);
  TH1F* hRow123TimeDet6 = new TH1F("row123_time_det6", "row123_time_det6", 100, 0, 10000);

  // dE plots
  TH2F* hdEECenterPad = new TH2F("dEECenterPad", "dEECenterPad", 500, 0, 4000, 500, 0, 25000);
  hdEECenterPad->GetXaxis()->SetTitle("Energy [Channels]"); hdEECenterPad->GetXaxis()->CenterTitle();
  hdEECenterPad->GetYaxis()->SetTitle("dE [Channels]"); hdEECenterPad->GetYaxis()->CenterTitle(); hdEECenterPad->GetYaxis()->SetTitleOffset(1.4);
  TH2F* hdEECenterPadPunchthrough = new TH2F("dEECenterPadPunchthrough", "dEECenterPadPunchthrough", 500, 0, 4000, 500, 0, 25000);
  hdEECenterPadPunchthrough->GetXaxis()->SetTitle("Energy [Channels]"); hdEECenterPadPunchthrough->GetXaxis()->SetTitleSize(0.05); hdEECenterPadPunchthrough->GetXaxis()->CenterTitle();
  hdEECenterPadPunchthrough->GetYaxis()->SetTitle("dE [Channels]"); hdEECenterPadPunchthrough->GetYaxis()->SetTitleSize(0.05); hdEECenterPadPunchthrough->GetYaxis()->CenterTitle();

  TH2F* hdEECenterPadCal = new TH2F("dEECenterPadCal", "dEECenterPadCal", 500, 0, 20000, 500, 0, 25000);
  hdEECenterPadCal->GetXaxis()->SetTitle("Energy [keV]"); hdEECenterPadCal->GetXaxis()->SetTitleSize(0.05); hdEECenterPadCal->GetXaxis()->CenterTitle();
  hdEECenterPadCal->GetYaxis()->SetTitle("dE [Channels]"); hdEECenterPadCal->GetYaxis()->SetTitleSize(0.05); hdEECenterPadCal->GetYaxis()->CenterTitle(); hdEECenterPadCal->GetYaxis()->SetTitleOffset(1.);
  TH2F* hdEECenterPadPunchthroughCal = new TH2F("dEECenterPadPunchthroughCal", "dEECenterPadPunchthroughCal", 500, 0, 20000, 500, 0, 25000);
  hdEECenterPadPunchthroughCal->GetXaxis()->SetTitle("Energy [keV]"); hdEECenterPadPunchthroughCal->GetXaxis()->SetTitleSize(0.05); hdEECenterPadPunchthroughCal->GetXaxis()->CenterTitle();
  hdEECenterPadPunchthroughCal->GetYaxis()->SetTitle("dE [Channels]"); hdEECenterPadPunchthroughCal->GetYaxis()->SetTitleSize(0.05); hdEECenterPadPunchthroughCal->GetYaxis()->CenterTitle(); hdEECenterPadPunchthroughCal->GetYaxis()->SetTitleOffset(1.);

  TH1F* hCenterPadTime = new TH1F("centerPadTime", "centerPadTime", 1000, 0, 20480);

  // Vertex vs E plots
  TH2F* hVertexE = new TH2F("vertexE", "vertexE", 128, 0, 128, 500, 0, 4000);
  hVertexE->GetXaxis()->SetTitle("Vertex [pad #]"); hVertexE->GetXaxis()->SetTitleSize(0.05); hVertexE->GetXaxis()->CenterTitle();
  hVertexE->GetYaxis()->SetTitle("Energy [Channels]"); hVertexE->GetYaxis()->SetTitleSize(0.05); hVertexE->GetYaxis()->CenterTitle(); hVertexE->GetYaxis()->SetTitleOffset(1.);
  TH2F* hVertexEPunchthrough = new TH2F("vertexEPunchthrough", "vertexEPunchthrough", 128, 0, 128, 500, 0, 4000);
  hVertexEPunchthrough->GetXaxis()->SetTitle("Vertex [pad #]"); hVertexEPunchthrough->GetXaxis()->SetTitleSize(0.05); hVertexEPunchthrough->GetXaxis()->CenterTitle();
  hVertexEPunchthrough->GetYaxis()->SetTitle("Energy [Channels]"); hVertexEPunchthrough->GetYaxis()->SetTitleSize(0.05); hVertexEPunchthrough->GetYaxis()->CenterTitle(); hVertexEPunchthrough->GetYaxis()->SetTitleOffset(1.);
  TH2F* hVertexECal = new TH2F("vertexECal", "vertexECal", 128, 0, 128, 500, 0, 20000);
  hVertexECal->GetXaxis()->SetTitle("Vertex [pad #]"); hVertexECal->GetXaxis()->SetTitleSize(0.05); hVertexECal->GetXaxis()->CenterTitle();
  hVertexECal->GetYaxis()->SetTitle("Energy [keV]"); hVertexECal->GetYaxis()->SetTitleSize(0.05); hVertexECal->GetYaxis()->CenterTitle(); hVertexECal->GetYaxis()->SetTitleOffset(1.);
  TH2F* hVertexEPunchthroughCal = new TH2F("vertexEPunchthroughCal", "vertexEPunchthroughCal", 128, 0, 128, 500, 0, 20000);
  hVertexEPunchthroughCal->GetXaxis()->SetTitle("Vertex [pad #]"); hVertexEPunchthroughCal->GetXaxis()->SetTitleSize(0.05); hVertexEPunchthroughCal->GetXaxis()->CenterTitle();
  hVertexEPunchthroughCal->GetYaxis()->SetTitle("Energy [keV]"); hVertexEPunchthroughCal->GetYaxis()->SetTitleSize(0.05); hVertexEPunchthroughCal->GetYaxis()->CenterTitle(); hVertexEPunchthroughCal->GetYaxis()->SetTitleOffset(1.);

  hCSDet5 = new TH1F("CSDet5", "CSDDet5", 40, 0, 5); hCSDet5->Sumw2();
  hCSDet5->GetXaxis()->SetTitle("Center of Mass Energy [MeV]"); hCSDet5->GetXaxis()->SetTitleSize(0.05); hCSDet5->GetXaxis()->CenterTitle(); hCSDet5->GetXaxis()->SetTitleOffset(0.95);
  hCSDet5->GetYaxis()->SetTitle("Cross Section [b/sr]"); hCSDet5->GetYaxis()->SetTitleSize(0.05); hCSDet5->GetYaxis()->CenterTitle(); hCSDet5->GetYaxis()->SetTitleOffset(0.95);

  hCSDet5Counts = new TH1F("CSDet5Counts", "CSDDet5Counts", 40, 0, 5); hCSDet5Counts->Sumw2();
  hCSDet5Counts->GetXaxis()->SetTitle("Center of Mass Energy [MeV]"); hCSDet5Counts->GetXaxis()->SetTitleSize(0.05); hCSDet5Counts->GetXaxis()->CenterTitle();
  hCSDet5Counts->GetYaxis()->SetTitle("Counts"); hCSDet5Counts->GetYaxis()->SetTitleSize(0.05); hCSDet5Counts->GetYaxis()->CenterTitle(); hCSDet5Counts->GetYaxis()->SetTitleOffset(1.);

  // For average energy deposited
  std::vector<Double_t> padCumulative[6][128];
  std::vector<Double_t> padCumulativeScaled[6][128];
  std::vector<Double_t> rowCumulativeScaled[128];
  std::vector<Double_t> padCumulativeTime[6][128];

  printf("Starting Main Loop\n");

  Long64_t nbytes = 0, nb = 0;
  // for(Long64_t jentry = 0; jentry < 400; jentry++) {
  for(Long64_t jentry = 0; jentry < nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if(ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;

    if(jentry != 0 && jentry % 200 == 0) printf("Processed %lld events\n",jentry);

    // Find hit in Si detectors
    std::vector<siDetect> siDetect_;

    // Find hit in CsI Detectors
    std::vector<csiDetect> csiDetect_;

    // Central Pads in Micromegas
    std::vector<mmCenter> mmCenter_;
    std::vector<mmCenter> mmCenterMatched_;

    // Beam Left Strips in Micromegas
    std::vector<mmStripChain> mmLeftStrip_;

    // Beam Left Chains in Micromegas
    std::vector<mmStripChain> mmLeftChain_;

    // Beam Right Strips in Micromegas
    std::vector<mmStripChain> mmRightStrip_;

    // Beam Right Chains in Micromegas
    std::vector<mmStripChain> mmRightChain_;

    // IC
    Int_t icE = 0.;
    Int_t icT = 0.;

    for(Int_t i = 0; i < mmMul; i++) {
      if(mmChan[i] == 11 || mmChan[i] == 22 || mmChan[i] == 45 || mmChan[i] == 56) continue; // Skip FPN Channels

      // Micromegas
      if(mmCobo[i] == 0) {
        // Asad0 - All Center Pads
        if(mmAsad[i] == 0) {
          // Aget0
          if(mmAget[i] == 0) {
            std::pair<Int_t, Int_t> pad = MM_Map_Asad0_Aget0[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, mmEnergy[i], mmTime[i]};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, mmEnergy[i]*scale[pad.first][pad.second], mmTime[i]};
            mmCenterMatched_.push_back(mmHitMatched);
          }
          // Aget1
          else if(mmAget[i] == 1) {
            std::pair<Int_t, Int_t> pad = MM_Map_Asad0_Aget1[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, mmEnergy[i], mmTime[i]};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, mmEnergy[i]*scale[pad.first][pad.second], mmTime[i]};
            mmCenterMatched_.push_back(mmHitMatched);
          }
          // Aget2
          else if(mmAget[i] == 2) {
            std::pair<Int_t, Int_t> pad = MM_Map_Asad0_Aget2[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, mmEnergy[i], mmTime[i]};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, mmEnergy[i]*scale[pad.first][pad.second], mmTime[i]};
            mmCenterMatched_.push_back(mmHitMatched);
          }
          // Aget3
          else if(mmAget[i] == 3) {
            std::pair<Int_t, Int_t> pad = MM_Map_Asad0_Aget3[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, mmEnergy[i], mmTime[i]};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, mmEnergy[i]*scale[pad.first][pad.second], mmTime[i]};
            mmCenterMatched_.push_back(mmHitMatched);
          }
        }
        // Asad1
        else if(mmAsad[i] == 1) {
          // Aget0
          if(mmAget[i] == 0) {
            std::pair<Int_t, Int_t> pad = MM_Map_Asad1_Aget0[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, mmEnergy[i], mmTime[i]};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, mmEnergy[i]*scale[pad.first][pad.second], mmTime[i]};
            mmCenterMatched_.push_back(mmHitMatched);
          }
          // Aget1
          else if(mmAget[i] == 1) {
            std::pair<Int_t, Int_t> pad = MM_Map_Asad1_Aget1[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, mmEnergy[i], mmTime[i]};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, mmEnergy[i]*scale[pad.first][pad.second], mmTime[i]};
            mmCenterMatched_.push_back(mmHitMatched);
          }
          // Aget2
          else if(mmAget[i] == 2) {
            std::pair<Int_t, Int_t> pad = MM_Map_Asad1_Aget2[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, mmEnergy[i], mmTime[i]};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, mmEnergy[i]*scale[pad.first][pad.second], mmTime[i]};
            mmCenterMatched_.push_back(mmHitMatched);
          }
          // Aget3
          else if(mmAget[i] == 3) {
            std::pair<Int_t, Int_t> pad = MM_Map_Asad1_Aget3[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, mmEnergy[i], mmTime[i]};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, mmEnergy[i]*scale[pad.first][pad.second], mmTime[i]};
            mmCenterMatched_.push_back(mmHitMatched);
          }
        }
        // Asad2
        else if(mmAsad[i] == 2) {
          // Aget0 - Strips and Chains (Beam left)
          if(mmAget[i] == 0) {
            Int_t bin = MM_Map_Asad2_Aget0[mmChan[i]];
            if(mmChan[i] < 34) { // Strips
              mmStripChain mmHit = {bin, mmEnergy[i], mmTime[i]};
              mmLeftStrip_.push_back(mmHit);
            }
            else if(mmChan[i] > 33 && mmChan[i] < 68) { // Chains
              mmStripChain mmHit = {bin, mmEnergy[i], mmTime[i]};
              mmLeftChain_.push_back(mmHit);
            }
          }
          // Aget1 - Strips and Chains (Beam left)
          else if(mmAget[i] == 1) {
            Int_t bin = MM_Map_Asad2_Aget1[mmChan[i]];
            if(mmChan[i] < 34) { // Strips
              mmStripChain mmHit = {bin, mmEnergy[i], mmTime[i]};
              mmLeftStrip_.push_back(mmHit);
            }
            else if(mmChan[i] > 33 && mmChan[i] < 68) { // Chains
              mmStripChain mmHit = {bin, mmEnergy[i], mmTime[i]};
              mmLeftChain_.push_back(mmHit);
            }
          }
          // Aget2 - Outside Central Pads
          else if(mmAget[i] == 2) {
            std::pair<Int_t, Int_t> pad = MM_Map_Asad2_Aget2[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, mmEnergy[i], mmTime[i]};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, mmEnergy[i]*scale[pad.first][pad.second], mmTime[i]};
            mmCenterMatched_.push_back(mmHitMatched);
          }
          // Aget3 - Outside Central Pads
          else if(mmAget[i] == 3) {
            std::pair<Int_t, Int_t> pad = MM_Map_Asad2_Aget3[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, mmEnergy[i], mmTime[i]};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, mmEnergy[i]*scale[pad.first][pad.second], mmTime[i]};
            mmCenterMatched_.push_back(mmHitMatched);
          }
        }
        // Asad3
        else if(mmAsad[i] == 3) {
          // Aget0 - Strips and Chains (Beam right)
          if(mmAget[i] == 0) {
            Int_t bin = MM_Map_Asad3_Aget0[mmChan[i]];
            if(mmChan[i] < 34) { // Strips
              mmStripChain mmHit = {bin, mmEnergy[i], mmTime[i]};
              mmRightStrip_.push_back(mmHit);
            }
            else if(mmChan[i] > 33 && mmChan[i] < 68) { // Chains
              mmStripChain mmHit = {bin, mmEnergy[i], mmTime[i]};
              mmRightChain_.push_back(mmHit);
            }
          }
          // Aget1 - Strips and Chains (Beam right)
          else if(mmAget[i] == 1) {
            Int_t bin = MM_Map_Asad3_Aget1[mmChan[i]];
            if(mmChan[i] < 34) { // Strips
              mmStripChain mmHit = {bin, mmEnergy[i], mmTime[i]};
              mmRightStrip_.push_back(mmHit);
            }
            else if(mmChan[i] > 33 && mmChan[i] < 68) { // Chains
              mmStripChain mmHit = {bin, mmEnergy[i], mmTime[i]};
              mmRightChain_.push_back(mmHit);
            }
          }
          // Aget2
          else if(mmAget[i] == 2) {
            std::pair<Int_t, Int_t> pad = MM_Map_Asad3_Aget2[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, mmEnergy[i], mmTime[i]};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, mmEnergy[i]*scale[pad.first][pad.second], mmTime[i]};
            mmCenterMatched_.push_back(mmHitMatched);
          }
          // Aget3
          else if(mmAget[i] == 3) {
            std::pair<Int_t, Int_t> pad = MM_Map_Asad3_Aget3[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, mmEnergy[i], mmTime[i]};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, mmEnergy[i]*scale[pad.first][pad.second], mmTime[i]};
            mmCenterMatched_.push_back(mmHitMatched);
          }
        }
      }

      else if(mmCobo[i] == 1) {
        // Si Detectors
        if(mmAsad[i] == 0) {
          // Forward Detectors
          if(mmAget[i] == 0) {
            if(siForwardMap.count(mmChan[i]) == 0) continue;
            Int_t detect = siForwardMap[mmChan[i]].first;
            Int_t quad = siForwardMap[mmChan[i]].second;
            Int_t channel = mmChan[i];
            siDetect siHit = {detect, quad, channel, mmEnergy[i], mmTime[i]};
            siDetect_.push_back(siHit);
          }
          // Beam Left Detectors
          else if(mmAget[i] == 3) {
            if(siLeftMap.count(mmChan[i]) == 0) continue;
            Int_t detect = siLeftMap[mmChan[i]].first;
            Int_t quad = siLeftMap[mmChan[i]].second;
            Int_t channel = mmChan[i];
            siDetect siHit = {detect+10, quad, channel, mmEnergy[i], mmTime[i]};
            siDetect_.push_back(siHit);
          }
        }

        // IC/CsI
        else if(mmAsad[i] == 1) {
          // Ionization Chamber
          if(mmAget[i] == 0 && mmChan[i] == 10) {
            icE = mmEnergy[i];
            icT = mmTime[i];
          }
          // CsI Detectors
          else if(mmAget[i] == 3) {
            if(mmChan[i] < 43) {
              csiDetect csiHit = {csiForwardMap[mmChan[i]], mmEnergy[i], mmTime[i]};
              csiDetect_.push_back(csiHit);
            }
            else {
              csiDetect csiHit = {csiLeftMap[mmChan[i]], mmEnergy[i], mmTime[i]};
              csiDetect_.push_back(csiHit);
            }
          }
        }
      }
    }

    // Only one silicon should fire
    if(siDetect_.size() != 1) continue;

    hIonizationChamberE->Fill(icE);
    hIonizationChamberT->Fill(icT);

    // if(icE<1250) continue;
    if(icE < 1400 || icE > 1900) continue;
    if(icT < 5500 || icT > 6500) continue;

    siDet = siDetect_[0].detect;
    siQuad = siDetect_[0].quad;
    siChannel = siDetect_[0].channel;
    siEnergy = siDetect_[0].energy;
    siEnergyCal = 0.;
    if(siDet < 10) {
      siEnergyCal = siEnergy*siEForwardCalibration[siDet][siQuad].first + siEForwardCalibration[siDet][siQuad].second;
    }
    siTime = siDetect_[0].time;

    if(siDet < 10) {
      hSiEForwardTotal[siDet]->Fill(siEnergy);
      hSiTForwardTotal[siDet]->Fill(siTime);
      hSiEForward[siDet][siQuad]->Fill(siEnergy);
    }
    else {
      hSiELeftTotal[siDet - 10]->Fill(siEnergy);
      hSiTLeftTotal[siDet - 10]->Fill(siTime);
      hSiELeft[siDet - 10][siQuad]->Fill(siEnergy);
    }

    // Find if CsI behind Si fired
    punchthrough = false;
    for(UInt_t i = 0; i < csiDetect_.size(); i++) {
      if(csiDetect_[i].detect == siDet) {
        punchthrough = true;
        hSiCsIForward[siDet]->Fill(siEnergy, csiDetect_[i].energy);
      }
    }

    // Reduce MM Center to one entry per row
    std::map<Int_t, Double_t> centralPadPosition;
    std::map<Int_t, Double_t> centralPadTotalEnergy;
    std::map<Int_t, Double_t> centralPadTime;
    std::map<Int_t, Int_t> centralPadTotal;
    for(UInt_t i = 0; i < 128; i++) {
      centralPadTotal[i] = 0;
    }
    for(UInt_t i = 0; i < mmCenterMatched_.size(); i++) {
      if(centralPadTotal[mmCenterMatched_[i].row] == 0) {
        centralPadPosition[mmCenterMatched_[i].row] = (mmCenterMatched_[i].column*3.5 - 8.75)*mmCenterMatched_[i].energy;
        centralPadTotalEnergy[mmCenterMatched_[i].row] = mmCenterMatched_[i].energy;
        centralPadTime[mmCenterMatched_[i].row] = mmCenterMatched_[i].time;
      }
      else {
        centralPadPosition[mmCenterMatched_[i].row] += (mmCenterMatched_[i].column*3.5 - 8.75)*mmCenterMatched_[i].energy;
        centralPadTotalEnergy[mmCenterMatched_[i].row] += mmCenterMatched_[i].energy;
        centralPadTime[mmCenterMatched_[i].row] += mmCenterMatched_[i].time;
      }
      centralPadTotal[mmCenterMatched_[i].row]++;
    }

    std::vector<mmTrack> mmCenterBeamTotal;
    std::map<Int_t, Double_t>::iterator it;
    for(it = centralPadPosition.begin(); it != centralPadPosition.end(); it++) {
      Int_t row = it->first;
      mmTrack hit = {0, 0., 0., 0., 0., 0., 0};
      hit.row = row;
      hit.xPosition = centralPadPosition[row]/centralPadTotalEnergy[row];
      hit.yPosition = row*rowConversion + rowConversionOffset;
      hit.time = centralPadTime[row]/static_cast<Double_t>(centralPadTotal[row]) - siTime;
      hit.energy = centralPadTotalEnergy[row];
      hit.height = heightOffset - hit.time*driftVelocity;
      hit.total = centralPadTotal[row];
      mmCenterBeamTotal.push_back(hit);
    }

    // Sorting mmCenterBeamTotal
    if(mmCenterBeamTotal.size() > 0) {
      std::sort(mmCenterBeamTotal.begin(), mmCenterBeamTotal.end(), sortByRowMMTrack());
    }

    // Sorting mmLeftChain
    if(mmLeftChain_.size() > 0) {
      std::sort(mmLeftChain_.begin(), mmLeftChain_.end(), sortByRowMMStripChain());
    }

    // Sorting mmLeftStrip
    if(mmLeftStrip_.size() > 0) {
      std::sort(mmLeftStrip_.begin(), mmLeftStrip_.end(), sortByRowMMStripChain());
    }

    // Sorting mmRightChain
    if(mmRightChain_.size() > 0) {
      std::sort(mmRightChain_.begin(), mmRightChain_.end(), sortByRowMMStripChain());
    }

    // Sorting mmRightStrip
    if(mmRightStrip_.size() > 0) {
      std::sort(mmRightStrip_.begin(), mmRightStrip_.end(), sortByRowMMStripChain());
    }

    std::vector<mmTrack> stripChainMatched;
    std::vector<mmTrack> stripChainRaw;
    if(mmLeftChain_.size() > 0 && mmLeftStrip_.size() > 0) {
      StripChainMatch(stripChainMatched, stripChainRaw, mmLeftChain_, mmLeftStrip_, true, siTime);
    }
    if(mmRightChain_.size() > 0 && mmRightStrip_.size() > 0) {
      StripChainMatch(stripChainMatched, stripChainRaw, mmRightChain_, mmRightStrip_, false, siTime);
    }

    //  ** End of event by event analysis ** //
  }

  for(int i=0; i<10 ; i++) {
    hSiEForwardTotal[i]->Write();
    // hSiTForwardTotal[i]->Write();
    // for(int j=0; j<4; j++) {
      // hSiEForward[i][j]->Write();
      // hSiEForwardCal[i][j]->Write();
      // hSiTForward[i][j]->Write();
    // }
    // hCsIEForward[i]->Write();
    // hCsITForward[i]->Write();
    // hSiCsIForward[i]->Write();
  }

  for(int i=0; i<6 ; i++) {
    hSiELeftTotal[i]->Write();
    // for(int j=0; j<4; j++) {
      // hSiELeft[i][j]->Write();
      // hSiTLeft[i][j]->Write();
    // }
    // hCsIELeft[i]->Write();
    // hCsITLeft[i]->Write();
  }

  hIonizationChamberE->Write();
  hIonizationChamberT->Write();

  // hMicroMegasCenterCumulative->Write();
  // hMicroMegasCenterEnergyCumulative->Write();
  // hMicroMegasCenterEnergyAverage->Write();
  // hMicroMegasCenterEnergyAverageScaled->Write();
  // hMicroMegasCenterTimeAverage->Write();

  // hdEECenterPad->Write();
  // hdEECenterPadPunchthrough->Write();
  // hdEECenterPadCal->Write();
  // hdEECenterPadPunchthroughCal->Write();
  // hCenterPadTime->Write();

  // hVertexE->Write();
  // hVertexEPunchthrough->Write();
  // hVertexECal->Write();
  // hVertexEPunchthroughCal->Write();

  // hCSDet5->Write();
  // hCSDet5Counts->Write();

  // hRow116TimeDet5->Write();
  // hRow116TimeDet6->Write();
  // hRow123TimeDet5->Write();
  // hRow123TimeDet6->Write();

  file->Close();
}

void Spectra::StripChainMatch(std::vector<mmTrack> &stripChainMatched, std::vector<mmTrack> &stripChainRaw, std::vector<mmStripChain> chain_,
                              std::vector<mmStripChain> strip_, Bool_t leftSide, Double_t siTime) {
  // stripChainMatched.clear();
  // stripChainRaw.clear();
  std::vector<mmTrack> totalTime0;
  StripChainMatchingTime(totalTime0, chain_, strip_, leftSide, siTime, 0);
  StripChainMatchingTime(stripChainRaw, chain_, strip_, leftSide, siTime, 0);
  size_t numTimeBuckets = StripChainTime0TimeBuckets(totalTime0);
  if(numTimeBuckets < 4) {
    StripChainMatchingBoxTime0(stripChainMatched, totalTime0);
  }
  else if(chain_.size() > 8 && strip_.size() > 8) {
    StripChainMatchingTimeSlopeHough(stripChainMatched, chain_, strip_, leftSide, siTime, 10.);
  }
  else {
    StripChainMatchingTime(stripChainMatched, chain_, strip_, leftSide, siTime, 0);
  }
}

size_t Spectra::StripChainTime0TimeBuckets(std::vector<mmTrack> matched) {
  // Function that finds the number of different time buckets for the time0 algorithm
  // This is used if you want to change the strip/chain matching algorithm (if all on the same plane)

  std::map<Int_t, Int_t> timeMap;

  for(UInt_t i = 0; i < matched.size(); i++) {
    timeMap[matched[i].time] = 1;
  }

  return timeMap.size();
}

size_t Spectra::StripChainNumberTimeBuckets(std::vector<mmStripChain> chain, std::vector<mmStripChain> strip) {
  // Function that finds the number of different time buckets for strips and chains
  // This is used if you want to change the strip/chain matching algorithm (if all on the same plane)

  std::map<Int_t, Int_t> individualTimeStripMap;
  std::map<Int_t, Int_t> individualTimeChainMap;

  for(UInt_t i = 0; i < chain.size(); i++) {
    individualTimeChainMap[chain[i].time] = 1;
  }
  for(UInt_t i = 0; i < strip.size(); i++) {
    individualTimeStripMap[strip[i].time] = 1;
  }

  size_t sizeTimeChainMap = individualTimeChainMap.size();
  size_t sizeTimeStripMap = individualTimeStripMap.size();

  return std::min(sizeTimeChainMap, sizeTimeStripMap);
}

void Spectra::StripChainMatchingOutward(std::vector<mmTrack> &stripChainMatched, std::vector<mmStripChain> chain,
                                        std::vector<mmStripChain> strip, Bool_t leftSide, Double_t siTime) {
  // Function to match strips and chains outward. Meaning that the matched strips and chains must be moving away from the central region
  // Strip and chains need to be sorted for this
  // Position transform = 10.5 + 1.75/2 + 1.75*chain (chain starting at 0)
  // Row transform = (strip - 1)*2
  // Row position transform = Row * rowConversion + rowConversionOffset

  std::vector<mmStripChain> stripCopy = strip;
  for(UInt_t i = 0; i < chain.size(); i++) {
    Double_t timeChain = chain[i].time;
    // Look in the strip vector for this time
    std::vector<mmStripChain>::iterator it = std::find_if(stripCopy.begin(), stripCopy.end(),
      [timeChain] (const mmStripChain& d) { return d.time == timeChain; });
    if(it != stripCopy.end()) {
      Double_t position = 10.5 + 1.75/2 + 1.75*chain[i].row;
      if(leftSide) position = -position;
      Int_t row = (*it).row*2;
      mmTrack hit = {0, 0., 0., 0., 0., 0., 0};
      hit.row = row;
      hit.xPosition = position;
      hit.yPosition = hit.row*rowConversion + rowConversionOffset;
      hit.time = chain[i].time - siTime;
      hit.energy = (*it).energy;
      hit.height = heightOffset - hit.time*driftVelocity;
      hit.total = 1;
      // if(fabs(hit.position) > 150) std::cout << chain[i].row << '\t' << (*it).row << '\t' << hit.time/1000. << std::endl;
      stripChainMatched.push_back(hit);
      stripCopy.erase(stripCopy.begin(), stripCopy.begin() + std::distance(stripCopy.begin(), it) + 1);
    }
    else continue;
  }
}

void Spectra::StripChainMatchingBox(std::vector<mmTrack> &stripChainMatched, std::vector<mmStripChain> chain,
                                    std::vector<mmStripChain> strip, Bool_t leftSide, Double_t siTime) {
  // Function to match strips and chains. This is the simplest algorithm, making two points in the side
  // region which is where the particle entered and exited.
  // Strips and chains need to be sorted for this
  // Position transform = 10.5 + 1.75/2 + 1.75*chain (chain starting at 0)
  // Row transform = (strip - 1)*2
  // Row position transform = Row * rowConversion + rowConversionOffset
  Size_t chainSize = chain.size();
  Size_t stripSize = strip.size();

  Double_t position0 = 10.5 + 1.75/2 + 1.75*chain[0].row;
  if(leftSide) position0 = -position0;
  Int_t row0 = strip[0].row*2;
  mmTrack hit0 = {0, 0., 0., 0., 0., 0., 0};
  hit0.row = row0;
  hit0.xPosition = position0;
  hit0.yPosition = hit0.row*rowConversion + rowConversionOffset;
  hit0.time = (chain[0].time + strip[0].time)/2. - siTime;
  hit0.energy = strip[0].energy;
  hit0.height = heightOffset - hit0.time*driftVelocity;
  hit0.total = 1;

  Double_t position1 = 10.5 + 1.75/2 + 1.75*chain[chainSize - 1].row;
  if(leftSide) position1 = -position1;
  Int_t row1 = strip[stripSize - 1].row*2;
  mmTrack hit1 = {0, 0., 0., 0., 0., 0., 0};
  hit1.row = row1;
  hit1.xPosition = position1;
  hit1.yPosition = hit1.row*rowConversion + rowConversionOffset;
  hit1.time = (chain[chainSize - 1].time + strip[stripSize - 1].time)/2. - siTime;
  hit1.energy = strip[stripSize - 1].energy;
  hit1.height = heightOffset - hit1.time*driftVelocity;
  hit1.total = 1;

  stripChainMatched.push_back(hit0);
  stripChainMatched.push_back(hit1);
}

void Spectra::StripChainMatchingBoxTime0(std::vector<mmTrack> &stripChainMatched, std::vector<mmTrack> time0) {
  UInt_t size = time0.size();
  if(size == 0) return;
  stripChainMatched.push_back(time0[0]);
  stripChainMatched.push_back(time0[size - 1]);
}

void Spectra::StripChainMatchingTime(std::vector<mmTrack> &stripChainMatched, std::vector<mmStripChain> chain,
                                     std::vector<mmStripChain> strip, Bool_t leftSide, Double_t siTime,
                                     Int_t timeWindow) {
  // Function to match strips and chains. Chains and strips are matched if their time is within the timeWindow parameter
  // The timeWindow parameter is the window of time buckets to look for using the timeResolution parameter
  // which is the ns -> bucket
  // Position transform = 10.5 + 1.75/2 + 1.75*chain (chain starting at 0)
  // Row transform = strip*2
  // Row position transform = Row * rowConversion + rowConversionOffset

  // Remove all strips and chains with time outside boundary
  std::vector<mmStripChain> stripReduced;
  std::vector<mmStripChain> chainReduced;
  for(UInt_t i = 0; i < strip.size(); i++) {
    if(strip[i].time < 8000 || strip[i].time > 11000) continue;
    stripReduced.push_back(strip[i]);
  }
  for(UInt_t i = 0; i < chain.size(); i++) {
    if(chain[i].time < 8000 || chain[i].time > 11000) continue;
    chainReduced.push_back(chain[i]);
  }

  for(UInt_t i = 0; i < chainReduced.size(); i++) {
    Double_t timeChain = chainReduced[i].time;
    // Loop over strips and find strips when the time is within the timeWindow
    for(UInt_t j = 0; j < stripReduced.size(); j++) {
      if((stripReduced[j].time > timeChain - timeWindow*timeResolution - timeResolution/2.) &&
          (stripReduced[j].time < timeChain + timeWindow*timeResolution + timeResolution/2.)) {
        mmTrack hit = {0, 0., 0., 0., 0., 0., 0};
        Double_t position = 10.5 + 1.75/2 + 1.75*chainReduced[i].row;
        if(leftSide) position = -position;
        Int_t row = stripReduced[j].row*2;
        hit.row = row;
        hit.xPosition = position;
        hit.yPosition = hit.row*rowConversion + rowConversionOffset;
        hit.time = chainReduced[i].time - siTime;
        hit.energy = stripReduced[j].energy;
        hit.height = heightOffset - hit.time*driftVelocity;
        hit.total = 1;
        stripChainMatched.push_back(hit);
      }
    }
  }
}

void Spectra::StripChainMatchingTimeSlopeFit(std::vector<mmTrack> &stripChainMatched, std::vector<mmStripChain> chain,
                                          std::vector<mmStripChain> strip, Bool_t leftSide, Double_t siTime,
                                          Double_t timeWindow) {
  // Function to match strips and chains by fitting with time individually
  // Then is able to correlate strips and chains by time. Requires multiple time steps in both strips and chains
  // After the fit, it goes through each chain and finds strips based off of the fit when match within the timeWindow input parameter
  // Position transform = 10.5 + 1.75/2 + 1.75*chain (chain starting at 0)
  // Row transform = (strip - 1)*2
  // Row position transform = Row * rowConversion + rowConversionOffset

  TGraph* hStrip = new TGraph();
  TGraph* hChain = new TGraph();

  // Fill hStrip
  for(UInt_t i = 0; i < strip.size(); i++) {
    hStrip->SetPoint(i, strip[i].row, strip[i].time);
  }

  // Fill hChain
  for(UInt_t i = 0; i < chain.size(); i++) {
    hChain->SetPoint(i, chain[i].row, chain[i].time);
  }

  // Fit each with a 1-D polynomial
  hStrip->Fit("pol1", "Q");
  hChain->Fit("pol1", "Q");
  TF1 *fitStrip = hStrip->GetFunction("pol1");
  TF1 *fitChain = hChain->GetFunction("pol1");

  hStrip->SetName(Form("hStrip_%lld", entry));
  hChain->SetName(Form("hChain_%lld", entry));
//  hStrip->Write();
//  hChain->Write();

  Double_t stripP0 = fitStrip->GetParameter(0);
  Double_t stripP1 = fitStrip->GetParameter(1);
  Double_t chainP0 = fitChain->GetParameter(0);
  Double_t chainP1 = fitChain->GetParameter(1);

  // Loop through chains and find strips within timeWindow
  for(UInt_t i = 0; i < chain.size(); i++) {
    Double_t chainTime = chain[i].row*chainP1 + chainP0;
    for(UInt_t j = 0; j < strip.size(); j++) {
      Double_t stripTime = strip[j].row*stripP1 + stripP0;
      if((stripTime > chainTime - timeWindow) && (stripTime < chainTime + timeWindow)) {
        // Check if chains and strips are actually within a timebucket
        if((chain[i].time < strip[j].time - timeResolution) ||
           (chain[i].time > strip[j].time + timeResolution)) continue;
        mmTrack hit = {0, 0., 0., 0., 0., 0., 0};
        Double_t position = 10.5 + 1.75/2 + 1.75*chain[i].row;
        if(leftSide) position = -position;
        Int_t row = strip[j].row*2;
        hit.row = row;
        hit.xPosition = position;
        hit.yPosition = hit.row*rowConversion + rowConversionOffset;
        hit.time = chain[i].time - siTime;
        hit.energy = strip[j].energy;
        hit.height = heightOffset - hit.time*driftVelocity;
        hit.total = 1;
        stripChainMatched.push_back(hit);
      }
    }
  }

  delete fitStrip;
  delete fitChain;

  delete hStrip;
  delete hChain;
}

void Spectra::StripChainMatchingTimeSlopeHough(std::vector<mmTrack> &stripChainMatched, std::vector<mmStripChain> chain,
                                             std::vector<mmStripChain> strip, Bool_t leftSide, Double_t siTime,
                                             Double_t timeWindow) {
  // Function to match strips and chains by fitting with time individually (using Hough 2D)
  // Then is able to correlate strips and chains by time. Requires multiple time steps in both strips and chains
  // After the fit, it goes through each chain and finds strips based off of the fit when match within the timeWindow input parameter
  // Position transform = 10.5 + 1.75/2 + 1.75*chain (chain starting at 0)
  // Row transform = (strip - 1)*2
  // Row position transform = Row * rowConversion + rowConversionOffset

  std::vector<mmTrack> newChain;
  for(UInt_t i = 0; i < chain.size(); i++) {
    mmTrack hit = {0, 0., 0., 0., 0., 0., 0};
    hit.xPosition = chain[i].row;
    hit.yPosition = chain[i].time/timeResolution; // Convert to time buckets
    newChain.push_back(hit);
  }

  std::vector<mmTrack> newStrip;
  for(UInt_t i = 0; i < strip.size(); i++) {
    mmTrack hit = {0, 0., 0., 0., 0., 0., 0};
    hit.xPosition = strip[i].row;
    hit.yPosition = strip[i].time/timeResolution; // Convert to time buckets
    newStrip.push_back(hit);
  }

  Hough2D *houghChain = new Hough2D();
  houghChain->SetPoints(newChain);
  houghChain->CalculateHoughXY();
  Double_t chainTheta = houghChain->GetMaxThetaXY();
  Double_t chainD = houghChain->GetMaxDXY();
  delete houghChain;

  Hough2D *houghStrip = new Hough2D();
  houghStrip->SetPoints(newStrip);
  houghStrip->CalculateHoughXY();
  Double_t stripTheta = houghStrip->GetMaxThetaXY();
  Double_t stripD = houghStrip->GetMaxDXY();
  delete houghStrip;

  // Loop through chains and find strips within timeWindow
  for(UInt_t i = 0; i < chain.size(); i++) {
    Double_t chainTime = -(cos(chainTheta*M_PI/180.)/sin(chainTheta*M_PI/180.))*chain[i].row +
        chainD/sin(chainTheta*M_PI/180.);
    chainTime *= 40.;
    for(UInt_t j = 0; j < strip.size(); j++) {
//      Double_t stripTime = strip[j].row*stripP1 + stripP0;
      Double_t stripTime = -(cos(stripTheta*M_PI/180.)/sin(stripTheta*M_PI/180.))*strip[j].row +
          stripD/sin(stripTheta*M_PI/180.);
      stripTime *= 40.;
      if((stripTime > chainTime - timeWindow) && (stripTime < chainTime + timeWindow)) {
        // Check if chains and strips are actually within a timebucket
        if((chain[i].time < strip[j].time - timeResolution) ||
           (chain[i].time > strip[j].time + timeResolution)) continue;
        mmTrack hit = {0, 0., 0., 0., 0., 0., 0};
        Double_t position = 10.5 + 1.75/2 + 1.75*chain[i].row;
        if(leftSide) position = -position;
        Int_t row = strip[j].row*2;
        hit.row = row;
        hit.xPosition = position;
        hit.yPosition = hit.row*rowConversion + rowConversionOffset;
        hit.time = chain[i].time - siTime;
        hit.energy = strip[j].energy;
        hit.height = heightOffset - hit.time*driftVelocity;
        hit.total = 1;
        stripChainMatched.push_back(hit);
      }
    }
  }
}
