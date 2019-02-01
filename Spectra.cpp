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
  auto *chain = new TChain("mfmData");

  // Home
  // TString PathToFiles = "/hd/research/data/run0817a/rootM2R-WaveformReduced/run";

  // Laptop
  TString PathToFiles = "/Users/joshhooker/Desktop/data/run0817a/run";

  // Alpha source test in gas
  // chain->Add(PathToFiles+"004.root");

  // 12C Runs
  // chain->Add(PathToFiles + "111.root");
  // chain->Add(PathToFiles + "113.root");
  // chain->Add(PathToFiles + "114.root");
  // chain->Add(PathToFiles + "116.root");
  // chain->Add(PathToFiles + "117.root");
  // chain->Add(PathToFiles + "118.root");
  // chain->Add(PathToFiles + "120.root");
  // chain->Add(PathToFiles + "121.root");
  // chain->Add(PathToFiles + "122.root");
  // chain->Add(PathToFiles + "127.root");
  // chain->Add(PathToFiles + "128.root");
  // chain->Add(PathToFiles + "129.root");
  // chain->Add(PathToFiles + "136.root");
  // chain->Add(PathToFiles + "138.root");
  // chain->Add(PathToFiles + "139.root");
  // chain->Add(PathToFiles + "140.root");
  // chain->Add(PathToFiles + "141.root");
  // chain->Add(PathToFiles + "142.root");
  // chain->Add(PathToFiles + "144.root");

  // 8B beam with IC
  chain->Add(PathToFiles + "158.root");

  // 8B Runs
  // chain->Add(PathToFiles + "175.root");
  // chain->Add(PathToFiles + "178.root");
  // chain->Add(PathToFiles + "180.root");
  // chain->Add(PathToFiles + "181.root");
  // chain->Add(PathToFiles + "182.root");
  // chain->Add(PathToFiles + "183.root");
  // chain->Add(PathToFiles + "184.root");
  // chain->Add(PathToFiles + "185.root");
  // chain->Add(PathToFiles + "186.root");
  // chain->Add(PathToFiles + "187.root");
  // chain->Add(PathToFiles + "188.root");
  // chain->Add(PathToFiles + "189.root");
  // chain->Add(PathToFiles + "190.root");
  // chain->Add(PathToFiles + "191.root");
  // chain->Add(PathToFiles + "192.root");
  // chain->Add(PathToFiles + "193.root");
  // chain->Add(PathToFiles + "195.root");
  // chain->Add(PathToFiles + "196.root");
  // chain->Add(PathToFiles + "197.root");
  // chain->Add(PathToFiles + "198.root");
  // chain->Add(PathToFiles + "199.root");
  // chain->Add(PathToFiles + "200.root");
  // chain->Add(PathToFiles + "201.root");
  // chain->Add(PathToFiles + "202.root");
  // chain->Add(PathToFiles + "203.root");
  // chain->Add(PathToFiles + "204.root");
  // chain->Add(PathToFiles + "205.root");
  // chain->Add(PathToFiles + "206.root");
  // chain->Add(PathToFiles + "207.root");
  // chain->Add(PathToFiles + "208.root");
  // chain->Add(PathToFiles + "210.root");
  // chain->Add(PathToFiles + "211.root");
  // chain->Add(PathToFiles + "212.root");
  // chain->Add(PathToFiles + "213.root");
  // chain->Add(PathToFiles + "214.root");
  // chain->Add(PathToFiles + "215.root");
  // chain->Add(PathToFiles + "216.root");
  // chain->Add(PathToFiles + "217.root");
  // chain->Add(PathToFiles + "218.root");
  // chain->Add(PathToFiles + "219.root");
  // chain->Add(PathToFiles + "220.root");
  // chain->Add(PathToFiles + "221.root");
  // chain->Add(PathToFiles + "223.root");
  // chain->Add(PathToFiles + "224.root");

  // Alpha source test in vacuum

  // Alpha source test in gas
  // chain->Add(PathToFiles + "232.root");

  return chain;
}

void Spectra::Loop() {

  if (fChain == nullptr) return;

  /////////////////
  // Set up cuts //
  /////////////////
  TFile *cutFile = TFile::Open("cuts.root");

  // dE cuts (Calibrated Si energy vs dE)
  for(uint i = 0; i < 10; i++) {
    if(i == 8) continue;
    dEEForwardCut[i] = static_cast<TCutG*>(cutFile->Get(Form("dEEForward_d%dCut", i)));
  }

  // Angle vs Energy Cut
  angleTotEnergyCut[0] = static_cast<TCutG*>(cutFile->Get("angleTotEnergy_d0Cut"));
  angleTotEnergyCut[1] = static_cast<TCutG*>(cutFile->Get("angleTotEnergy_d1Cut"));
  angleTotEnergyCut[9] = static_cast<TCutG*>(cutFile->Get("angleTotEnergy_d9Cut"));

  // Chain/Strip vs Time Cuts
  for(uint i = 0; i < 10; i++) {
    if(i == 4 || i == 5) continue;
    for(uint j = 0; j < 4; j++) {
      if(i == 7 && j == 0) continue;
      timeChainForwardCut[i][j] = static_cast<TCutG*>(cutFile->Get(Form("timeChainForward_d%d_q%dCut", i, j)));
      timeStripForwardCut[i][j] = static_cast<TCutG*>(cutFile->Get(Form("timeStripForward_d%d_q%dCut", i, j)));
    }
  }

  cwtE_CentralCut = static_cast<TCutG*>((cutFile->Get("cwtE_CentralCut")));
  cwtE_CentralProtonCut = static_cast<TCutG*>(cutFile->Get("cwtE_CentralProtonCut"));

  // Cut on raw Si E and CsI E (channels)
  siCsiEForwardCut[0] = static_cast<TCutG*>(cutFile->Get("siCsiEForward_d0Cut"));
  siCsiEForwardCut[1] = static_cast<TCutG*>(cutFile->Get("siCsiEForward_d1Cut"));
  siCsiEForwardCut[4] = static_cast<TCutG*>(cutFile->Get("siCsiEForward_d4Cut"));
  siCsiEForwardCut[5] = static_cast<TCutG*>(cutFile->Get("siCsiEForward_d5Cut"));
  siCsiEForwardCut[9] = static_cast<TCutG*>(cutFile->Get("siCsiEForward_d9Cut"));

  cutFile->Close();

  InitChannelMap();
  InitHistograms();
  InitCanvas();
  InitVariables();

  InitSiEForwardCalibration();
  InitCsIECalibration();

  InitCentralPadGainMatch();
  InitAverageBeamEnergy();

  InitMaxPeak();

  InitTree();

  long nentries = fChain->GetEntriesFast();

  double padEnergy[6][128];
  int padTotal[6][128];

  for(Int_t i = 0; i < 6; i++) {
    for(Int_t j = 0; j < 128; j++) {
      padEnergy[i][j] = 0.;
      padTotal[i][j] = 0;
    }
  }

  printf("Starting Main Loop\n");

  long nbytes = 0, nb = 0;
  // for(long jentry = 0; jentry < 5000; jentry++) {
  // for(long jentry = 4747; jentry < 4748; jentry++) {
  for(long jentry = 0; jentry < nentries; jentry++) {
    long ientry = LoadTree(jentry);
    if(ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;

    entry = jentry;

    if(jentry != 0 && jentry % 10000 == 0) printf("Processed %ld events\n", jentry);

    // Find hit in Si detectors
    std::vector<siDetect> siDetect_;

    // Find hit in CsI Detectors
    std::vector<csiDetect> csiDetect_;

    // Central Pads in Micromegas
    std::vector<mmCenter> mmCenter_;
    std::vector<mmCenter> mmCenterMatched_;

    // Beam Left Strips in Micromegas
    std::vector<mmChainStrip> mmLeftStrip_;

    // Beam Left Chains in Micromegas
    std::vector<mmChainStrip> mmLeftChain_;

    // Beam Right Strips in Micromegas
    std::vector<mmChainStrip> mmRightStrip_;

    // Beam Right Chains in Micromegas
    std::vector<mmChainStrip> mmRightChain_;

    // IC
    double icE = 0.;
    double icT = 0.;

    for(int i = 0; i < mmMul; i++) {
      if(mmChan[i] == 11 || mmChan[i] == 22 || mmChan[i] == 45 || mmChan[i] == 56) continue; // Skip FPN Channels

      if(mmEnergy[i] < 0) continue;

      // Micromegas
      if(mmCobo[i] == 0) {
        // Asad0 - All Center Pads
        if(mmAsad[i] == 0) {
          // Aget0
          if(mmAget[i] == 0) {
            std::pair<int, int> pad = MM_Map_Asad0_Aget0[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, static_cast<double>(mmEnergy[i]), static_cast<double>(mmTime[i]), static_cast<double>(mmPa[i][3]/mmPa[i][1]),
                             static_cast<int>(floor((static_cast<double>(mmTime[i]) + 0.01)/timeResolution))};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, static_cast<double>(mmEnergy[i])*scale[pad.first][pad.second], static_cast<double>(mmTime[i]),
                                     static_cast<double>(mmPa[i][3]/mmPa[i][1]), static_cast<int>(floor((static_cast<double>(mmTime[i]) + 0.01)/timeResolution))};
            mmCenterMatched_.push_back(mmHitMatched);
          }
          // Aget1
          else if(mmAget[i] == 1) {
            std::pair<int, int> pad = MM_Map_Asad0_Aget1[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, static_cast<double>(mmEnergy[i]),static_cast<double>(mmTime[i]), static_cast<double>(mmPa[i][3]/mmPa[i][1]),
                             static_cast<int>(floor((static_cast<double>(mmTime[i]) + 0.01)/timeResolution))};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, static_cast<double>(mmEnergy[i])*scale[pad.first][pad.second], static_cast<double>(mmTime[i]),
                                     static_cast<double>(mmPa[i][3]/mmPa[i][1]), static_cast<int>(floor((static_cast<double>(mmTime[i]) + 0.01)/timeResolution))};
            mmCenterMatched_.push_back(mmHitMatched);
          }
          // Aget2
          else if(mmAget[i] == 2) {
            std::pair<int, int> pad = MM_Map_Asad0_Aget2[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, static_cast<double>(mmEnergy[i]), static_cast<double>(mmTime[i]), static_cast<double>(mmPa[i][3]/mmPa[i][1]),
                             static_cast<int>(floor((static_cast<double>(mmTime[i]) + 0.01)/timeResolution))};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, static_cast<double>(mmEnergy[i])*scale[pad.first][pad.second], static_cast<double>(mmTime[i]),
                                     static_cast<double>(mmPa[i][3]/mmPa[i][1]), static_cast<int>(floor((static_cast<double>(mmTime[i]) + 0.01)/timeResolution))};
            mmCenterMatched_.push_back(mmHitMatched);
          }
          // Aget3
          else if(mmAget[i] == 3) {
            std::pair<int, int> pad = MM_Map_Asad0_Aget3[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, static_cast<double>(mmEnergy[i]), static_cast<double>(mmTime[i]), static_cast<double>(mmPa[i][3]/mmPa[i][1]),
                             static_cast<int>(floor((static_cast<double>(mmTime[i]) + 0.01)/timeResolution))};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, static_cast<double>(mmEnergy[i])*scale[pad.first][pad.second], static_cast<double>(mmTime[i]),
                                     static_cast<double>(mmPa[i][3]/mmPa[i][1]), static_cast<int>(floor((static_cast<double>(mmTime[i]) + 0.01)/timeResolution))};
            mmCenterMatched_.push_back(mmHitMatched);
          }
        }
        // Asad1
        else if(mmAsad[i] == 1) {
          // Aget0
          if(mmAget[i] == 0) {
            std::pair<int, int> pad = MM_Map_Asad1_Aget0[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, static_cast<double>(mmEnergy[i]), static_cast<double>(mmTime[i]), static_cast<double>(mmPa[i][3]/mmPa[i][1]),
                             static_cast<int>(floor((static_cast<double>(mmTime[i]) + 0.01)/timeResolution))};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, static_cast<double>(mmEnergy[i])*scale[pad.first][pad.second], static_cast<double>(mmTime[i]),
                                     static_cast<double>(mmPa[i][3]/mmPa[i][1]), static_cast<int>(floor((static_cast<double>(mmTime[i]) + 0.01)/timeResolution))};
            mmCenterMatched_.push_back(mmHitMatched);
          }
          // Aget1
          else if(mmAget[i] == 1) {
            std::pair<int, int> pad = MM_Map_Asad1_Aget1[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, static_cast<double>(mmEnergy[i]), static_cast<double>(mmTime[i]), static_cast<double>(mmPa[i][3]/mmPa[i][1]),
                             static_cast<int>(floor((static_cast<double>(mmTime[i]) + 0.01)/timeResolution))};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, static_cast<double>(mmEnergy[i])*scale[pad.first][pad.second], static_cast<double>(mmTime[i]),
                                     static_cast<double>(mmPa[i][3]/mmPa[i][1]), static_cast<int>(floor((static_cast<double>(mmTime[i]) + 0.01)/timeResolution))};
            mmCenterMatched_.push_back(mmHitMatched);
          }
          // Aget2
          else if(mmAget[i] == 2) {
            std::pair<int, int> pad = MM_Map_Asad1_Aget2[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, static_cast<double>(mmEnergy[i]), static_cast<double>(mmTime[i]), static_cast<double>(mmPa[i][3]/mmPa[i][1]),
                             static_cast<int>(floor((static_cast<double>(mmTime[i]) + 0.01)/timeResolution))};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, static_cast<double>(mmEnergy[i])*scale[pad.first][pad.second], static_cast<double>(mmTime[i]),
                                     static_cast<double>(mmPa[i][3]/mmPa[i][1]), static_cast<int>(floor((static_cast<double>(mmTime[i]) + 0.01)/timeResolution))};
            mmCenterMatched_.push_back(mmHitMatched);
          }
          // Aget3
          else if(mmAget[i] == 3) {
            std::pair<int, int> pad = MM_Map_Asad1_Aget3[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, static_cast<double>(mmEnergy[i]), static_cast<double>(mmTime[i]), static_cast<double>(mmPa[i][3]/mmPa[i][1]),
                             static_cast<int>(floor((static_cast<double>(mmTime[i]) + 0.01)/timeResolution))};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, static_cast<double>(mmEnergy[i])*scale[pad.first][pad.second], static_cast<double>(mmTime[i]),
                                     static_cast<double>(mmPa[i][3]/mmPa[i][1]), static_cast<int>(floor((static_cast<double>(mmTime[i]) + 0.01)/timeResolution))};
            mmCenterMatched_.push_back(mmHitMatched);
          }
        }
        // Asad2
        else if(mmAsad[i] == 2) {
          // Aget0 - Strips and Chains (Beam left)
          if(mmAget[i] == 0) {
            int bin = MM_Map_Asad2_Aget0[mmChan[i]];
            if(mmChan[i] < 34) { // Strips
              mmChainStrip mmHit = {bin, static_cast<double>(mmEnergy[i]), static_cast<double>(mmTime[i]), static_cast<int>(floor((static_cast<double>(mmTime[i]) + 0.01)/timeResolution))};
              mmLeftStrip_.push_back(mmHit);
            }
            else if(mmChan[i] > 33 && mmChan[i] < 68) { // Chains
              mmChainStrip mmHit = {bin, static_cast<double>(mmEnergy[i]), static_cast<double>(mmTime[i]), static_cast<int>(floor((static_cast<double>(mmTime[i]) + 0.01)/timeResolution))};
              mmLeftChain_.push_back(mmHit);
            }
          }
          // Aget1 - Strips and Chains (Beam left)
          else if(mmAget[i] == 1) {
            int bin = MM_Map_Asad2_Aget1[mmChan[i]];
            if(mmChan[i] < 34) { // Strips
              mmChainStrip mmHit = {bin, static_cast<double>(mmEnergy[i]), static_cast<double>(mmTime[i]), static_cast<int>(floor((static_cast<double>(mmTime[i]) + 0.01)/timeResolution))};
              mmLeftStrip_.push_back(mmHit);
            }
            else if(mmChan[i] > 33 && mmChan[i] < 68) { // Chains
              mmChainStrip mmHit = {bin, static_cast<double>(mmEnergy[i]), static_cast<double>(mmTime[i]), static_cast<int>(floor((static_cast<double>(mmTime[i]) + 0.01)/timeResolution))};
              mmLeftChain_.push_back(mmHit);
            }
          }
          // Aget2 - Outside Central Pads
          else if(mmAget[i] == 2) {
            std::pair<int, int> pad = MM_Map_Asad2_Aget2[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, static_cast<double>(mmEnergy[i]), static_cast<double>(mmTime[i]), static_cast<double>(mmPa[i][3]/mmPa[i][1]),
                             static_cast<int>(floor((static_cast<double>(mmTime[i]) + 0.01)/timeResolution))};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, static_cast<double>(mmEnergy[i])*scale[pad.first][pad.second], static_cast<double>(mmTime[i]),
                                     static_cast<double>(mmPa[i][3]/mmPa[i][1]), static_cast<int>(floor((static_cast<double>(mmTime[i]) + 0.01)/timeResolution))};
            mmCenterMatched_.push_back(mmHitMatched);
          }
          // Aget3 - Outside Central Pads
          else if(mmAget[i] == 3) {
            std::pair<int, int> pad = MM_Map_Asad2_Aget3[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, static_cast<double>(mmEnergy[i]), static_cast<double>(mmTime[i]), static_cast<double>(mmPa[i][3]/mmPa[i][1]),
                             static_cast<int>(floor((static_cast<double>(mmTime[i]) + 0.01)/timeResolution))};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, static_cast<double>(mmEnergy[i])*scale[pad.first][pad.second], static_cast<double>(mmTime[i]),
                                     static_cast<double>(mmPa[i][3]/mmPa[i][1]), static_cast<int>(floor((static_cast<double>(mmTime[i]) + 0.01)/timeResolution))};
            mmCenterMatched_.push_back(mmHitMatched);
          }
        }
        // Asad3
        else if(mmAsad[i] == 3) {
          // Aget0 - Strips and Chains (Beam right)
          if(mmAget[i] == 0) {
            int bin = MM_Map_Asad3_Aget0[mmChan[i]];
            if(mmChan[i] < 34) { // Strips
              mmChainStrip mmHit = {bin, static_cast<double>(mmEnergy[i]), static_cast<double>(mmTime[i]), static_cast<int>(floor((static_cast<double>(mmTime[i]) + 0.01)/timeResolution))};
              mmRightStrip_.push_back(mmHit);
            }
            else if(mmChan[i] > 33 && mmChan[i] < 68) { // Chains
              mmChainStrip mmHit = {bin, static_cast<double>(mmEnergy[i]), static_cast<double>(mmTime[i]), static_cast<int>(floor((static_cast<double>(mmTime[i]) + 0.01)/timeResolution))};
              mmRightChain_.push_back(mmHit);
            }
          }
          // Aget1 - Strips and Chains (Beam right)
          else if(mmAget[i] == 1) {
            int bin = MM_Map_Asad3_Aget1[mmChan[i]];
            if(mmChan[i] < 34) { // Strips
              mmChainStrip mmHit = {bin, static_cast<double>(mmEnergy[i]), static_cast<double>(mmTime[i]), static_cast<int>(floor((static_cast<double>(mmTime[i]) + 0.01)/timeResolution))};
              mmRightStrip_.push_back(mmHit);
            }
            else if(mmChan[i] > 33 && mmChan[i] < 68) { // Chains
              mmChainStrip mmHit = {bin, static_cast<double>(mmEnergy[i]), static_cast<double>(mmTime[i]), static_cast<int>(floor((static_cast<double>(mmTime[i]) + 0.01)/timeResolution))};
              mmRightChain_.push_back(mmHit);
            }
          }
          // Aget2
          else if(mmAget[i] == 2) {
            std::pair<int, int> pad = MM_Map_Asad3_Aget2[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, static_cast<double>(mmEnergy[i]), static_cast<double>(mmTime[i]), static_cast<double>(mmPa[i][3]/mmPa[i][1]),
                              static_cast<int>(floor((static_cast<double>(mmTime[i]) + 0.01)/timeResolution))};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, static_cast<double>(mmEnergy[i])*scale[pad.first][pad.second], static_cast<double>(mmTime[i]),
                                     static_cast<double>(mmPa[i][3]/mmPa[i][1]), static_cast<int>(floor((static_cast<double>(mmTime[i]) + 0.01)/timeResolution))};
            mmCenterMatched_.push_back(mmHitMatched);
          }
          // Aget3
          else if(mmAget[i] == 3) {
            std::pair<int, int> pad = MM_Map_Asad3_Aget3[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, static_cast<double>(mmEnergy[i]), static_cast<double>(mmTime[i]), static_cast<double>(mmPa[i][3]/mmPa[i][1]),
                             static_cast<int>(floor((static_cast<double>(mmTime[i]) + 0.01)/timeResolution))};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, static_cast<double>(mmEnergy[i])*scale[pad.first][pad.second], static_cast<double>(mmTime[i]),
                                     static_cast<double>(mmPa[i][3]/mmPa[i][1]), static_cast<int>(floor((static_cast<double>(mmTime[i]) + 0.01)/timeResolution))};
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
            int detect = siForwardMap[mmChan[i]].first;
            int quad = siForwardMap[mmChan[i]].second;
            int channel = mmChan[i];
            siDetect siHit = {detect, quad, channel, static_cast<double>(mmEnergy[i]), static_cast<double>(mmTime[i]), static_cast<double>(mmPa[i][3]/mmPa[i][1])};
            siDetect_.push_back(siHit);
          }
          // Beam Left Detectors
          else if(mmAget[i] == 3) {
            if(siLeftMap.count(mmChan[i]) == 0) continue;
            int detect = siLeftMap[mmChan[i]].first;
            int quad = siLeftMap[mmChan[i]].second;
            int channel = mmChan[i];
            siDetect siHit = {detect+10, quad, channel, static_cast<double>(mmEnergy[i]), static_cast<double>(mmTime[i]), static_cast<double>(mmPa[i][3]/mmPa[i][1])};
            siDetect_.push_back(siHit);
          }
        }

        // IC/CsI
        else if(mmAsad[i] == 1) {
          // Ionization Chamber
          if(mmAget[i] == 0 && mmChan[i] == 10) {
            icE = static_cast<double>(mmEnergy[i]);
            icT = static_cast<double>(mmTime[i]);
          }
          // CsI Detectors
          else if(mmAget[i] == 3) {
            if(mmChan[i] < 43) {
              csiDetect csiHit = {csiForwardMap[mmChan[i]], static_cast<double>(mmEnergy[i]), static_cast<double>(mmTime[i])};
              csiDetect_.push_back(csiHit);
            }
            else {
              csiDetect csiHit = {csiLeftMap[mmChan[i]], static_cast<double>(mmEnergy[i]), static_cast<double>(mmTime[i])};
              csiDetect_.push_back(csiHit);
            }
          }
        }
      }
    }

    hIonizationChamberE->Fill(icE);
    hIonizationChamberT->Fill(icT);

    if(icE < 1400 || icE > 1900) continue;
    if(icT < 5500 || icT > 6500) continue;

    //************//
    // Micromegas //
    //************//

    // Gate on time for central pads
    std::vector<mmCenter> mmCenterMatchedReduced_;
//    for(auto mm : mmCenterMatched_) {
    for(auto mm : mmCenter_) {
      double time = mm.time;
      if(mm.row < 2) continue;
      if(mm.row < 112) {
        if(time < 4700 || time > 6300) continue;
        if(mm.energy < 200) continue;
        hCWTECentral->Fill(mm.energy, mm.cwt);
        if(!cwtE_CentralCut->IsInside(mm.energy, mm.cwt)) continue;
        mmCenterMatchedReduced_.push_back(mm);
      }
      // These should not fire
      else if(mm.row > 111 && (mm.column == 0 || mm.column == 5)) {
        if(time < 4700 || time > 6300) continue;
        if(mm.energy < 200) continue;
        hCWTECentral->Fill(mm.energy, mm.cwt);
        if(!cwtE_CentralCut->IsInside(mm.energy, mm.cwt)) continue;
        mmCenterMatchedReduced_.push_back(mm);
      }
      else {
        if(time < 4700 || time > 6300) continue;
        hCWTECentral->Fill(mm.energy, mm.cwt);
        // if(!cwtE_CentralProtonCut->IsInside(mm.energy, mm.cwt)) continue;
        mmCenterMatchedReduced_.push_back(mm);
      }
    }

    // Sorting mmCenterMatchedReduced_ by row and reduce noise
     std::vector<mmCenter> mmCenterMatchedReducedNoise_;
     if(!mmCenterMatchedReduced_.empty()) {
      std::sort(mmCenterMatchedReduced_.begin(), mmCenterMatchedReduced_.end(), sortByRowMMCenter());
      mmCenterMatchedReducedNoise_ = CenterReduceNoise(mmCenterMatchedReduced_);
    }

    // Reduce mmCenter to one entry per row for beam and protons
    std::map<int, double> centralPadPosition;
    std::map<int, double> centralPadTotalEnergy;
    std::map<int, double> centralPadTime;
    std::map<int, int> centralPadTotal;
    for(int i = 0; i < 128; i++) {
      centralPadTotal[i] = 0;
    }
    for(auto mm : mmCenterMatchedReducedNoise_) {
      double xPosition = mm.column*3.5 - 8.75;
      double yPosition = mm.row*rowConversion + rowConversionOffset;
      hMicroMegasCenterCumulative->Fill(mm.column - 3, mm.row);
      hMicroMegasCenterCumulativePositionRaw->Fill(xPosition, yPosition);
      hMicroMegasCenterTime->Fill(mm.row, mm.time);
      padEnergy[mm.column][mm.row] += mm.energy;
      padTotal[mm.column][mm.row]++;
      if(centralPadTotal[mm.row] == 0) {
        if(mm.energy < 0) continue;
        centralPadPosition[mm.row] = (mm.column*3.5 - 8.75)*mm.energy;
        centralPadTotalEnergy[mm.row] = mm.energy;
        centralPadTime[mm.row] = mm.time;
      }
      else {
        if(mm.energy < 0) continue;
        centralPadPosition[mm.row] += (mm.column*3.5 - 8.75)*mm.energy;
        centralPadTotalEnergy[mm.row] += mm.energy;
        centralPadTime[mm.row] += mm.time;
      }
      centralPadTotal[mm.row]++;
    }

    // Separate beam and proton in central Micromegas
    std::vector<mmTrack> mmCenterBeamTotal_;
    std::vector<mmTrack> mmCenterProton_;
    std::map<int, double>::iterator it;
    for(it = centralPadPosition.begin(); it != centralPadPosition.end(); it++) {
      int row = it->first;
      mmTrack hit = {0, 0., 0., 0., 0., 0., 0, 0};
      if(centralPadTotalEnergy[row] < 1) continue;
      if(centralPadTotal[row] == 0) continue;
      hit.row = row;
      hit.xPosition = centralPadPosition[row]/centralPadTotalEnergy[row];
      hit.yPosition = row*rowConversion + rowConversionOffset;
      hit.time = centralPadTime[row]/static_cast<double>(centralPadTotal[row]) - siTime;
      hit.energy = centralPadTotalEnergy[row];
      hit.height = heightOffset - hit.time*driftVelocity;
      hit.timeBucket = static_cast<int>(floor((centralPadTime[row]/static_cast<double>(centralPadTotal[row]) + 0.01)/timeResolution));
      hit.total = centralPadTotal[row];
      hMicroMegasCenterCumulativePosition->Fill(hit.xPosition, hit.yPosition);
      if(row < 112) {
        mmCenterBeamTotal_.push_back(hit);
      }
      else if(row > 111) {
        mmCenterProton_.push_back(hit);
      }
    }

  }

  double padAverage[6][128];
  for(Int_t i = 0; i < 6; i++) {
    for(Int_t j = 0; j < 128; j++) {
      if(padEnergy[i][j] < 1) continue;
      padAverage[i][j] = padEnergy[i][j]/static_cast<double>(padTotal[i][j]);
      int padAvgNum = static_cast<int>(padAverage[i][j]);
      for(Int_t k = 0; k < padAvgNum; k++) {
        hMicroMegasCenterEnergyAverageScaled->Fill(i - 3, j);
      }
      std::cout << i << '\t' << j << '\t' << padAverage[i][j] << '\t' << padAvgNum << std::endl;
    }
  }

  WriteHistograms();

  WriteTree();

  file->Close();
}

std::vector<mmCenter> Spectra::CenterReduceNoise(std::vector<mmCenter> center) {
  std::vector<mmCenter> mmNew;

  std::vector<int> rows[128];

  for(auto mm : center) {
    if(mm.row > 111 && mm.column != 0 && mm.column != 5) {
      mmNew.push_back(mm);
    }
    else {
      if(mm.row < 4) {
        rows[mm.row - 2].push_back(mm.column);
        mmNew.push_back(mm);
      }
      else {
        // Check previous 2 rows if they have any columns
        if(rows[mm.row - 4].empty() && rows[mm.row - 3].empty()) continue;
        if(!rows[mm.row - 3].empty()) {
          bool found = false;
          for(auto column : rows[mm.row - 3]) {
            if(column == mm.column - 1 || column == mm.column || column == mm.column + 1) found = true;
          }
          if(found) {
            bool inRow = false;
            for(auto column : rows[mm.row - 2]) {
              if(column == mm.column) inRow = true;
            }
            if(!inRow) {
              rows[mm.row - 2].push_back(mm.column);
              mmNew.push_back(mm);
            }
          }
        }
        else if(!rows[mm.row - 4].empty()) {
          bool found = false;
          for(auto column : rows[mm.row - 4]) {
            if(column == mm.column - 1 || column == mm.column || column == mm.column + 1) found = true;
          }
          if(found) {
            bool inRow = false;
            for(auto column : rows[mm.row - 2]) {
              if(column == mm.column) inRow = true;
            }
            if(!inRow) {
              rows[mm.row - 2].push_back(mm.column);
              mmNew.push_back(mm);
            }
          }
        }
      }
    }
  }

  return mmNew;
}
