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
  TString PathToFiles = "/hd/research/data/run0817a/rootM2R-WaveformReduced/run";

  // Laptop
  // TString PathToFiles = "/Users/joshhooker/Desktop/data/run0817a/run";

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
  // chain->Add(PathToFiles + "158.root");

  // 8B Runs
  chain->Add(PathToFiles + "175.root");
  chain->Add(PathToFiles + "178.root");
  chain->Add(PathToFiles + "180.root");
  chain->Add(PathToFiles + "181.root");
  chain->Add(PathToFiles + "182.root");
  chain->Add(PathToFiles + "183.root");
  chain->Add(PathToFiles + "184.root");
  chain->Add(PathToFiles + "185.root");
  chain->Add(PathToFiles + "186.root");
  chain->Add(PathToFiles + "187.root");
  chain->Add(PathToFiles + "188.root");
  chain->Add(PathToFiles + "189.root");
  chain->Add(PathToFiles + "190.root");
  chain->Add(PathToFiles + "191.root");
  chain->Add(PathToFiles + "192.root");
  chain->Add(PathToFiles + "193.root");
  chain->Add(PathToFiles + "195.root");
  chain->Add(PathToFiles + "196.root");
  chain->Add(PathToFiles + "197.root");
  chain->Add(PathToFiles + "198.root");
  chain->Add(PathToFiles + "199.root");
  chain->Add(PathToFiles + "200.root");
  chain->Add(PathToFiles + "201.root");
  chain->Add(PathToFiles + "202.root");
  chain->Add(PathToFiles + "203.root");
  chain->Add(PathToFiles + "204.root");
  chain->Add(PathToFiles + "205.root");
  chain->Add(PathToFiles + "206.root");
  chain->Add(PathToFiles + "207.root");
  chain->Add(PathToFiles + "208.root");
  chain->Add(PathToFiles + "210.root");
  chain->Add(PathToFiles + "211.root");
  chain->Add(PathToFiles + "212.root");
  chain->Add(PathToFiles + "213.root");
  chain->Add(PathToFiles + "214.root");
  chain->Add(PathToFiles + "215.root");
  chain->Add(PathToFiles + "216.root");
  chain->Add(PathToFiles + "217.root");
  chain->Add(PathToFiles + "218.root");
  chain->Add(PathToFiles + "219.root");
  chain->Add(PathToFiles + "220.root");
  chain->Add(PathToFiles + "221.root");
  chain->Add(PathToFiles + "223.root");
  chain->Add(PathToFiles + "224.root");

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

    bool central = false;
    bool left = false;
    bool right = false;

    //*********//
    // Silicon //
    //*********//

    // Look at all Si that fired and gate
    std::vector<siDetect> siDetectReduced_;
    for(auto si : siDetect_) {
      hCWTSiE->Fill(si.energy, si.cwt);
      if(si.time > 3600 && si.time < 5200) siDetectReduced_.push_back(si);
    }

    std::vector<protonDetect> protonDet_;
    for(auto si : siDetectReduced_) {
      // Skip bad Kiev detectors (side)
      if(si.detect > 9) continue;
      // Skip bad Kiev detectors (forward)
      if(si.detect == 2 || si.detect == 3 || si.detect == 6 || si.detect == 7) continue;
      protonDetect hit = {si.detect, si.quad, si.channel, si.energy, 0., si.time, 0., 0., 0., 0., false};
      protonDet_.push_back(hit);
    }

    for(auto &det : protonDet_) {
      double siEnergyCal = det.siEnergy*siEForwardCalibration[det.det][det.quad].first +
          siEForwardCalibration[det.det][det.quad].second;
      det.siEnergyCal = siEnergyCal;
      det.totalEnergy = siEnergyCal;
    }

    // Remove events only in bad detectors
    if(protonDet_.empty()) continue;

    //*****//
    // CsI //
    //*****//

    for(auto csi : csiDetect_) {
      if(csi.time < 11000 || csi.time > 12250) continue;
      for(auto &det : protonDet_) {
        if(det.det != csi.detect) continue;
        if(!siCsiEForwardCut[det.det]->IsInside(det.siEnergy, csi.energy)) continue;
        double csiEnergyCal = csi.energy*csiEForwardCalibration[det.det].first + csiEForwardCalibration[det.det].second;
        det.csiEnergy = csi.energy;
        det.csiEnergyCal = csiEnergyCal;
        det.csiTime = csi.time;
        det.totalEnergy += csiEnergyCal;
        det.punchthrough = true;
      }
    }

    //****************//
    // Plot that shit //
    //****************//

    // hSiFired->Fill(protonDet_.size());

    for(auto det : protonDet_) {
      hSiEForward[det.det][det.quad]->Fill(det.siEnergy);
      hSiEForwardDet[det.det]->Fill(det.siEnergy);

      hSiEForwardCal[det.det][det.quad]->Fill(det.siEnergyCal);
      hSiEForwardDetCal[det.det]->Fill(det.siEnergyCal);

      hSiTForward[det.det][det.quad]->Fill(det.siTime);
      hSiTForwardDet[det.det]->Fill(det.siTime);

      if(det.punchthrough) {
        hCsiEForward[det.det]->Fill(det.csiEnergy);
        hCsiEForwardCal[det.det]->Fill(det.csiEnergyCal);
        hCsiTForward[det.det]->Fill(det.csiTime);

        hCsiETForward[det.det]->Fill(det.csiEnergy, det.csiTime);

        hSiCsiEForward[det.det][det.quad]->Fill(det.siEnergy, det.csiEnergy);
        hSiCsiEForwardDet[det.det]->Fill(det.siEnergy, det.csiEnergy);

        hSiCsiEForwardCal[det.det][det.quad]->Fill(det.siEnergyCal, det.csiEnergy);
        hSiCsiEForwardDetCal[det.det]->Fill(det.siEnergyCal, det.csiEnergy);
      }

      hTotalEForward[det.det]->Fill(det.totalEnergy);

      hSumSiEForward[det.det][det.quad]->Fill(det.totalEnergy, det.siEnergy);
      hSumSiEForwardDet[det.det]->Fill(det.totalEnergy, det.siEnergy);

      hSumCsiEForward[det.det][det.quad]->Fill(det.totalEnergy, det.csiEnergy);
      hSumCsiEForwardDet[det.det]->Fill(det.totalEnergy, det.csiEnergy);
    }

    //************//
    // Micromegas //
    //************//

    // Gate on time for central pads
    std::vector<mmCenter> mmCenterMatchedReduced_;
    for(auto mm : mmCenterMatched_) {
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
        if(!cwtE_CentralProtonCut->IsInside(mm.energy, mm.cwt)) continue;
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

    // Sort Chains/Strips and Clean Strips
    std::vector<mmChainStrip> mmLeftStripCleaned_;
    std::vector<mmChainStrip> mmRightStripCleaned_;

    // Sorting mmLeftChain_
    if(!mmLeftChain_.empty()) {
      std::sort(mmLeftChain_.begin(), mmLeftChain_.end(), sortByRowMMChainStrip());
    }

    // Sorting mmLeftStrip_
    if(!mmLeftStrip_.empty()) {
      std::sort(mmLeftStrip_.begin(), mmLeftStrip_.end(), sortByRowMMChainStrip());
      // CleanLeftStrip(mmLeftStripCleaned_, mmLeftStrip_);
    }

    // Sorting mmRightChain_
    if(!mmRightChain_.empty()) {
      std::sort(mmRightChain_.begin(), mmRightChain_.end(), sortByRowMMChainStrip());
    }

    // Sorting mmRightStrip_
    if(!mmRightStrip_.empty()) {
      std::sort(mmRightStrip_.begin(), mmRightStrip_.end(), sortByRowMMChainStrip());
      // CleanRightStrip(mmRightStripCleaned_, mmRightStrip_);
    }

    //*************************************//
    // Get dE for side and central regions //
    //*************************************//

    // Central region
    double dECentral = 0.;
    int totalRowsCentral = 0;
    for(auto mm : mmCenterProton_) {
      if(mm.row > 115 && mm.row < 124 && mm.row != 117) {
        dECentral += mm.energy;
        totalRowsCentral++;
      }
    }
    if(totalRowsCentral > 0) {
      dECentral /= static_cast<double>(totalRowsCentral);
    }

    // Side Region
    double dELeft = 0.;
    int totalRowsLeft = 0;
    for(auto mm : mmLeftStrip_) {
      if(mm.row < 50 || mm.row > 62) continue;
      dELeft += mm.energy;
      totalRowsLeft++;
    }
    if(totalRowsLeft > 0) {
      dELeft /= static_cast<double>(totalRowsLeft);
    }
    double dERight = 0.;
    int totalRowsRight = 0;
    for(auto mm : mmRightStrip_) {
      if(mm.row < 50 || mm.row > 62) continue;
      dERight += mm.energy;
      totalRowsRight++;
    }
    if(totalRowsRight > 0) {
      dERight /= static_cast<double>(totalRowsRight);
    }

    std::vector<protonDetect> protonDetReduced_;
    for(auto det : protonDet_) {
      double dERegion = 0.;
      if(det.det < 4) dERegion = dERight;
      else if(det.det > 5) dERegion = dELeft;
      else dERegion = dECentral;
      hdEEForward[det.det]->Fill(det.siEnergy, dERegion);
      hdEEForwardCal[det.det]->Fill(det.siEnergyCal, dERegion);
      hdEEForwardCalTotal[det.det]->Fill(det.totalEnergy, dERegion);

      if(!dEEForwardCut[det.det]->IsInside(det.siEnergyCal, dERegion)) continue;
      protonDetReduced_.push_back(det);
    }

    if(protonDetReduced_.empty()) continue;

    hSiFired->Fill(protonDetReduced_.size());

    if(protonDetReduced_.size() != 1) continue;

    siDet = protonDetReduced_[0].det;
    siQuad = protonDetReduced_[0].quad;
    siChannel = protonDetReduced_[0].siChannel;
    siEnergy = protonDetReduced_[0].siEnergy;
    siEnergyCal = protonDetReduced_[0].siEnergyCal;
    siTime = protonDetReduced_[0].siTime;
    csiEnergy = protonDetReduced_[0].csiEnergy;
    csiEnergyCal = protonDetReduced_[0].csiEnergyCal;
    csiTime = protonDetReduced_[0].csiTime;
    totalEnergy = protonDetReduced_[0].totalEnergy;
    punchthrough = protonDetReduced_[0].punchthrough;
    if(siDet < 4) {
      dE = static_cast<double>(dERight);
      right = true;
    }
    else if(siDet > 5) {
      dE = static_cast<double>(dELeft);
      left = true;
    }
    else {
      dE = static_cast<double>(dECentral);
      central = true;
    }


    bool event = false;
    if(central) {
      event = AnalysisForwardCentral(mmCenterMatchedReducedNoise_, mmCenterBeamTotal_, mmCenterProton_, centralPadTotalEnergy,
                                     mmLeftChain_, mmLeftStrip_, mmRightChain_, mmRightStrip_);
    }
//    if(left || right) {
//      event = AnalysisForwardSide(mmCenterMatchedReducedNoise_, mmCenterBeamTotal_, mmCenterProton_,
//                                  mmLeftChain_, mmLeftStrip_, mmRightChain_, mmRightStrip_);
//    }

    if(!event) continue;

    //  ** End of event by event analysis ** //
    FillTree();
  }

  ReadSolidAngle();

  DivideTargetThickness(s1);
  SolidAngle(s1, 1);
  s1->Scale(1./numberB8);
  s1->Write();

  DivideTargetThickness(s3);
  SolidAngle(s3, 3);
  s3->Scale(1./numberB8);
  s3->Write();

  WriteHistograms();

  WriteCanvas();

  WriteTree();

  file->Close();
}

bool Spectra::AnalysisForwardCentral(std::vector<mmCenter> centerMatched_, std::vector<mmTrack> centerBeamTotal_,
                                     std::vector<mmTrack> centerProton_, std::map<int, double> centralPadTotalEnergy,
                                     std::vector<mmChainStrip> leftChain_, std::vector<mmChainStrip> leftStrip_,
                                     std::vector<mmChainStrip> rightChain_, std::vector<mmChainStrip> rightStrip_) {

  // Find dE for center region
  dE = 0.;
  int totalRows = 0;
  for(auto mm : centerProton_) {
    if(mm.row > 115 && mm.row < 124 && mm.row != 117) {
      if(mm.energy < 50) continue;
      dE += mm.energy;
      totalRows++;
    }
  }
  dE /= static_cast<double>(totalRows);

  for(auto mm : centerMatched_) {
    hTimeCentraldEForward[siDet]->Fill(mm.row, mm.time - siTime);
  }

  hdEEForward[siDet]->Fill(siEnergy, dE);
  hdEEForwardCal[siDet]->Fill(siEnergyCal, dE);
  hdEEForwardCalTotal[siDet]->Fill(totalEnergy, dE);

  if(!dEEForwardCut[siDet]->IsInside(siEnergyCal, dE)) return false;

//  // Match strips with chains
//  std::vector<mmTrack> chainStripMatchedLeft_;
//  std::vector<mmTrack> chainStripMatchedRight_;
//  std::vector<mmTrack> chainStripRawLeft_;
//  std::vector<mmTrack> chainStripRawRight_;
//  if(!leftChain_.empty() && !leftStrip_.empty()) {
//    ChainStripMatch(chainStripMatchedLeft_, chainStripRawLeft_, leftChain_, leftStrip_, true, siTime);
//  }
//  if(!rightChain_.empty() && !rightStrip_.empty()) {
//    ChainStripMatch(chainStripMatchedRight_, chainStripRawRight_, rightChain_, rightStrip_, false,
//                    siTime);
//  }

  bool oneProtonEvent = true;

  if(siDet != 4) return false;

//  // Clean left side track
//  double leftTrackAngle = 0;
//  if(chainStripMatchedLeft_.size() > 1) {
//    auto *fitLeft = new HoughTrack();
//    fitLeft->AddTrack(chainStripMatchedLeft_, siDet, siQuad);
//    double minDistOther = fitLeft->Fit();
//    std::vector<double> parsOther = fitLeft->GetPars();
//    delete fitLeft;
//    double m_xComponentOther = fabs(parsOther[1]);
//    leftTrackAngle = atan(m_xComponentOther);

//    int changes;
//    std::vector<mmTrack> leftTrackCleaned_;
//    std::vector<mmTrack> leftTrackCleanedRaw_;
//    bool change = true;
//    while(change) {
//      CleanSideTracks(leftTrackCleaned_, leftTrackCleanedRaw_, chainStripRawLeft_, leftTrackAngle, changes);
//      chainStripMatchedLeft_ = leftTrackCleaned_;
//      chainStripRawLeft_ = leftTrackCleanedRaw_;
//      if(changes == 0) change = false;

//      if(leftTrackCleanedRaw_.size() < 2 || chainStripMatchedLeft_.size() < 2) change = false;

//      auto *fitLeft = new HoughTrack();
//      fitLeft->AddTrack(chainStripMatchedLeft_, siDet, siQuad);
//      double minDistOther = fitLeft->Fit();
//      std::vector<double> parsOther = fitLeft->GetPars();
//      delete fitLeft;
//      double m_xComponentOther = fabs(parsOther[1]);
//      leftTrackAngle = atan(m_xComponentOther);
//    }
//  }

//  // Check dE of left track and see if proton
//  if(chainStripMatchedLeft_.size() > 1) {
//    double dESide = 0.;
//    int totalSide = 0;
//    for(auto mm : chainStripMatchedLeft_) {
//      dESide += mm.energy*cos(leftTrackAngle);
//      totalSide++;
//    }
//    dESide /= static_cast<double>(totalSide);

//    if(dESide > 200 && dESide < 1700) oneProtonEvent = false;
//  }

//  if(!oneProtonEvent) return false;

//  // Clean right side track
//  double rightTrackAngle = 0;
//  if(chainStripMatchedRight_.size() > 1) {
//    auto *fitRight = new HoughTrack();
//    fitRight->AddTrack(chainStripMatchedRight_, siDet, siQuad);
//    double minDistOther = fitRight->Fit();
//    std::vector<double> parsOther = fitRight->GetPars();
//    delete fitRight;
//    double m_xComponentOther = fabs(parsOther[1]);
//    rightTrackAngle = atan(m_xComponentOther);

//    int changes;
//    std::vector<mmTrack> rightTrackCleaned_;
//    std::vector<mmTrack> rightTrackCleanedRaw_;
//    bool change = true;
//    while(change) {
//      CleanSideTracks(rightTrackCleaned_, rightTrackCleanedRaw_, chainStripRawRight_, rightTrackAngle, changes);
//      chainStripMatchedRight_ = rightTrackCleaned_;
//      chainStripRawRight_ = rightTrackCleanedRaw_;
//      if(changes == 0) change = false;

//      if(rightTrackCleanedRaw_.size() < 2 || chainStripMatchedRight_.size() < 2) change = false;

//      auto *fitRight = new HoughTrack();
//      fitRight->AddTrack(chainStripMatchedRight_, siDet, siQuad);
//      double minDistOther = fitRight->Fit();
//      std::vector<double> parsOther = fitRight->GetPars();
//      delete fitRight;
//      double m_xComponentOther = fabs(parsOther[1]);
//      rightTrackAngle = atan(m_xComponentOther);
//    }
//  }

//  // Check dE of right track and see if proton
//  if(chainStripMatchedRight_.size() > 1) {
//    double dESide = 0.;
//    int totalSide = 0;
//    for(auto mm : chainStripMatchedRight_) {
//      dESide += mm.energy*cos(rightTrackAngle);
//      totalSide++;
//    }
//    dESide /= static_cast<double>(totalSide);

//    if(dESide > 200 && dESide < 1700) oneProtonEvent = false;
//  }

//  if(!oneProtonEvent) return false;

  // Remove events with no or very little counts for beam
  if(centerBeamTotal_.size() < 5) return false;

  int maxPeakRow;
  double maxPeakEnergy, avgPeakEnergy, maxPeakDeriv;
  FindMaxCentralEnergy(centerBeamTotal_, maxPeakRow, maxPeakEnergy, avgPeakEnergy, maxPeakDeriv);

  if(maxPeakRow < 5) return false;

  // Find maxPeakRow in maxPeakMap
  if(maxPeakMap.find(maxPeakRow) == maxPeakMap.end()) return false;

  double vertex = maxPeakMap[maxPeakRow];
  int vertexRow = 0;
  if(vertex > 0) {
    vertexRow = static_cast<int>(vertex/1.75);
  }

  // std::cout << maxPeakRow << '\t' << vertex << '\t' << vertexRow << std::endl;

  mmTrack centerVertexPoint;
  if(vertexRow == 0) {
    mmTrack hit = {0, 0., vertex, 0, 0, 0, 0, 0};
    centerVertexPoint = hit;
  }
  else {
    // Get vertex from centralBeam
    double closestVertexValue = 1000;
    mmTrack closestVertex;
    for(auto mm : centerBeamTotal_) {
      double diff = fabs(mm.yPosition - vertex);
      if(diff < closestVertexValue) {
        closestVertexValue = diff;
        closestVertex = mm;
      }
    }
    centerVertexPoint = closestVertex;
  }

  centerProton_.push_back(centerVertexPoint);

  auto *fitProton = new HoughTrack();
  fitProton->AddTrack(centerProton_, siDet, siQuad);
  double minDistProton = fitProton->Fit();
  std::vector<double> parsProton = fitProton->GetPars();
  delete fitProton;

  // Plot XZ hit position of forward non-central detectors
  double x, y, z;
  line(siYPosForward, parsProton, x, y, z);
  siPosX = x;
  siPosY = siYPosForward;
  siPosZ = z;
  hHitPositionsXZForward->Fill(x, z);
  hHitPositionsXZForwardInd[siDet]->Fill(x, z);

  double vertexToSi = siYPosForward - vertex;

  double siX, siY, siZ;
  line(siYPosForward, parsProton, siX, siY, siZ);

  TVector3 v1(siX, vertexToSi, siZ);
  TVector3 v2(0, vertexToSi, 0);
  angle = v1.Angle(v2);
  double cosAngle = cos(angle);

  double pathLength = vertexToSi/cosAngle;

  if(pathLength < 0) pathLength = 200.;
  if(pathLength > 700) pathLength = 700.;

  double protonEnergy = protonMethane->AddBack(totalEnergy/1000., pathLength);
  double beamEnergy = protonEnergy*(m1 + m2)*(m1 + m2)/(4.*m1*m2*cosAngle*cosAngle);

  cmEnergy = beamEnergy*m2/(m1 + m2);

  hVertexSiEForward[siDet]->Fill(siEnergy, vertexPositionY);
  hVertexSiEForwardCal[siDet]->Fill(siEnergyCal, vertexPositionY);
  hVertexSiEForwardCalTotal[siDet]->Fill(totalEnergy, vertexPositionY);

  hAngleEForward[siDet]->Fill(siEnergy, angle);
  hAngleEForwardCal[siDet]->Fill(siEnergyCal, angle);
  hAngleEForwardCalTotal[siDet]->Fill(totalEnergy, angle);
  hAngleEForwardProtonEnergy[siDet]->Fill(protonEnergy, angle);

  hVertexAngleForward[siDet]->Fill(vertexPositionY, angle);

  if(siDet == 4 && cmEnergy > 0.9) {
    s1->Fill(cmEnergy);
  }

  return true;
}

bool Spectra::AnalysisForwardSide(std::vector<mmCenter> centerMatched_, std::vector<mmTrack> centerBeamTotal_,
                                  std::vector<mmTrack> centerProton_,
                                  std::vector<mmChainStrip> leftChain_, std::vector<mmChainStrip> leftStrip_,
                                  std::vector<mmChainStrip> rightChain_, std::vector<mmChainStrip> rightStrip_) {
  // Sorting mmCenterBeamTotal_ by row
  if(!centerBeamTotal_.empty()) {
    std::sort(centerBeamTotal_.begin(), centerBeamTotal_.end(), sortByRowMMTrack());
  }

  bool left = false;
  bool right = false;
  if(siDet < 4) right = true;
  else if(siDet > 5) left = true;

  // Left Chain
  for(auto mm : leftChain_) {
    hTimeChainForward[siDet][siQuad]->Fill(mm.row, mm.time - siTime);
    hTimeChainForwardCumulative->Fill(mm.row, mm.time - siTime);
    hMicroMegasChainLeftCumulative->Fill(mm.row);
  }
  // Left Strip
  for(auto mm : leftStrip_) {
    hTimeStripForward[siDet][siQuad]->Fill(mm.row, mm.time - siTime);
    hTimeStripForwardCumulative->Fill(mm.row, mm.time - siTime);
    hMicroMegasStripLeftCumulative->Fill(mm.row);
  }
  // Right Chain
  for(auto mm : rightChain_) {
    hTimeChainForward[siDet][siQuad]->Fill(mm.row, mm.time - siTime);
    hTimeChainForwardCumulative->Fill(mm.row, mm.time - siTime);
    hMicroMegasChainRightCumulative->Fill(mm.row);
  }
  // Right Strip
  for(auto mm : rightStrip_) {
    hTimeStripForward[siDet][siQuad]->Fill(mm.row, mm.time - siTime);
    hTimeStripForwardCumulative->Fill(mm.row, mm.time - siTime);
    hMicroMegasStripRightCumulative->Fill(mm.row);
  }

  // Gate on time vs chain/strip
  std::vector<mmChainStrip> leftChainReduced_;
  std::vector<mmChainStrip> leftStripReduced_;
  std::vector<mmChainStrip> rightChainReduced_;
  std::vector<mmChainStrip> rightStripReduced_;
  if(left) {
    // Left Chain
    for(auto mm : leftChain_) {
      if(timeChainForwardCut[siDet][siQuad]->IsInside(mm.row, mm.time - siTime)) {
        leftChainReduced_.push_back(mm);
      }
    }
    // Left Strip
    for(auto mm : leftStrip_) {
      if(timeStripForwardCut[siDet][siQuad]->IsInside(mm.row, mm.time - siTime)) {
        leftStripReduced_.push_back(mm);
      }
    }
    for(auto mm : rightChain_) {
      rightChainReduced_.push_back(mm);
    }
    for(auto mm : rightStrip_) {
      rightStripReduced_.push_back(mm);
    }
  }
  else if(right) {
    // Right Chain
    for(auto mm : rightChain_) {
      if(timeChainForwardCut[siDet][siQuad]->IsInside(mm.row, mm.time - siTime)) {
        rightChainReduced_.push_back(mm);
      }
    }
    // Right Strip
    for(auto mm : rightStrip_) {
      if(timeStripForwardCut[siDet][siQuad]->IsInside(mm.row, mm.time - siTime)) {
        rightStripReduced_.push_back(mm);
      }
    }
    for(auto mm : leftChain_) {
      leftChainReduced_.push_back(mm);
    }
    for(auto mm : leftStrip_) {
      leftStripReduced_.push_back(mm);
    }
  }

  for(auto mm : centerMatched_) {
    hTimeCentralForward[siDet]->Fill(mm.row, mm.time - siTime);
  }

  // Match strips with chains
  std::vector<mmTrack> chainStripMatchedLeft_;
  std::vector<mmTrack> chainStripMatchedRight_;
  std::vector<mmTrack> chainStripRawLeft_;
  std::vector<mmTrack> chainStripRawRight_;
  if(!leftChainReduced_.empty() && !leftStripReduced_.empty()) {
    ChainStripMatch(chainStripMatchedLeft_, chainStripRawLeft_, leftChainReduced_, leftStripReduced_, true, siTime);
  }
  if(!rightChainReduced_.empty() && !rightStripReduced_.empty()) {
    ChainStripMatch(chainStripMatchedRight_, chainStripRawRight_, rightChainReduced_, rightStripReduced_, false,
                    siTime);
  }

  // Assign proton track for side regions
  std::vector<mmTrack> protonTrack_;
  std::vector<mmTrack> protonTrackRaw_;
  std::vector<mmTrack> otherTrack_;
  std::vector<mmTrack> otherTrackRaw_;

  double leftTrackLength = 0.;
  double rightTrackLength = 0.;
  if(!chainStripMatchedLeft_.empty()) leftTrackLength = ChainStripSize(chainStripMatchedLeft_);
  if(!chainStripMatchedRight_.empty()) rightTrackLength = ChainStripSize(chainStripMatchedRight_);

  double protonTrackLength = 0.;
  double otherTrackLength = 0.;
  double protonTrackAngle = 0.;
  double otherTrackAngle = 0.;
  if(right) {
    protonTrack_ = chainStripMatchedRight_;
    protonTrackRaw_ = chainStripRawRight_;
    otherTrack_ = chainStripMatchedLeft_;
    otherTrackRaw_ = chainStripRawLeft_;

    protonTrackLength = rightTrackLength;
    otherTrackLength = leftTrackLength;
  }
  else if(left) {
    protonTrack_ = chainStripMatchedLeft_;
    protonTrackRaw_ = chainStripRawLeft_;
    otherTrack_ = chainStripMatchedRight_;
    otherTrackRaw_ = chainStripRawRight_;

    protonTrackLength = leftTrackLength;
    otherTrackLength = rightTrackLength;
  }

  // If proton length is very short or 0, event is bad (no matching?)
  if(protonTrackLength < 0.1 || protonTrack_.empty()) return false;

  bool oneProtonEvent = true;
  bool goodEvent = true;

  // Clean other track
  if(otherTrack_.size() > 1) {
    auto *fitOther = new HoughTrack();
    fitOther->AddTrack(otherTrack_, siDet, siQuad);
    double minDistOther = fitOther->Fit();
    std::vector<double> parsOther = fitOther->GetPars();
    delete fitOther;
    double m_xComponentOther = fabs(parsOther[1]);
    otherTrackAngle = atan(m_xComponentOther);

    int changes;
    std::vector<mmTrack> otherTrackCleaned_;
    std::vector<mmTrack> otherTrackCleanedRaw_;
    bool change = true;
    while(change) {
      CleanSideTracks(otherTrackCleaned_, otherTrackCleanedRaw_, otherTrackRaw_, otherTrackAngle, changes);
      otherTrack_ = otherTrackCleaned_;
      otherTrackRaw_ = otherTrackCleanedRaw_;
      if(changes == 0) change = false;

      if(otherTrackCleanedRaw_.size() < 2 || otherTrack_.size() < 2) change = false;

      auto *fitOther = new HoughTrack();
      fitOther->AddTrack(otherTrack_, siDet, siQuad);
      double minDistOther = fitOther->Fit();
      std::vector<double> parsOther = fitOther->GetPars();
      delete fitOther;
      double m_xComponentOther = fabs(parsOther[1]);
      otherTrackAngle = atan(m_xComponentOther);
    }
  }

  // Check dE of other track and see if proton
  if(otherTrack_.size() > 1) {
    double dESide = 0.;
    int totalSide = 0;
    for(auto mm : otherTrackRaw_) {
      dESide += mm.energy*cos(otherTrackAngle);
      totalSide++;
    }
    dESide /= static_cast<double>(totalSide);

    if(dESide > 200 && dESide < 1700) oneProtonEvent = false;
  }

  if(!oneProtonEvent) return false;

  // Check if dE is proton (dE average is between 200 - 2700)
  if(!centerProton_.empty()) {
    // Find dE for center region
    double dECenter = 0.;
    int totalRows = 0;
    for(auto mm : centerProton_) {
      if(mm.row > 115 && mm.row < 124 && mm.row != 117) {
        if(mm.energy < 50) continue;
        dECenter += mm.energy;
        totalRows++;
      }
    }
    dECenter /= static_cast<double>(totalRows);

    if(dECenter > 200 && dECenter < 2700) oneProtonEvent = false;
  }

  if(!oneProtonEvent) return false;

  // Initial fit of proton track
  auto *fitProton = new HoughTrack();
  fitProton->AddTrack(protonTrack_, siDet, siQuad);
  double minDistProton = fitProton->FitRestricted();
  std::vector<double> parsProton = fitProton->GetPars();
  double houghAngleXY = fitProton->GetHoughAngleXY();
  double houghAngleYZ = fitProton->GetHoughAngleYZ();
  hHoughAngle[siDet]->Fill(houghAngleXY);
  delete fitProton;
  double m_xComponentProton = fabs(parsProton[1]);
  protonTrackAngle = atan(m_xComponentProton);

  // Clean proton track
  int changes;
  std::vector<mmTrack> protonTrackCleaned_;
  std::vector<mmTrack> protonTrackCleanedRaw_;
  bool change = true;
  while(change) {
    CleanSideTracks(protonTrackCleaned_, protonTrackCleanedRaw_, protonTrackRaw_, protonTrackAngle, changes);
    protonTrack_ = protonTrackCleaned_;
    protonTrackRaw_ = protonTrackCleanedRaw_;
    if(changes == 0) change = false;

    if(protonTrackRaw_.size() < 2 || protonTrack_.size() < 2) change = false;

    auto *fitProton = new HoughTrack();
    fitProton->AddTrack(protonTrack_, siDet, siQuad);
    double minDistProton = fitProton->FitRestricted();
    std::vector<double> parsProton = fitProton->GetPars();
    delete fitProton;
    double m_xComponentProton = fabs(parsProton[1]);
    protonTrackAngle = atan(m_xComponentProton);
  }

  if(protonTrack_.empty()) return false;

  // Remove really small tracks when no beam track
  if(centerBeamTotal_.empty()) {
    double protonTrackLength = ChainStripSize(protonTrack_);
    if(protonTrackLength < 10) goodEvent = false;
  }

  if(!oneProtonEvent) return false;
  if(!goodEvent) return false;

  // Find vertex using the proton track in side regions
  vertexPositionX = 0.;
  vertexPositionY = -400.;
  vertexPositionZ = 0.;

  SideVertexFinder(centerBeamTotal_, protonTrack_, vertexPositionX, vertexPositionY, vertexPositionZ, parsProton);
  if(vertexPositionY < -390) {
    SideVertexFinderSingleHelp(centerBeamTotal_, protonTrack_, vertexPositionX, vertexPositionY, vertexPositionZ, parsProton);
  }
  if(vertexPositionY < -390) {
    SideVertexFinderHelp(centerBeamTotal_, protonTrack_, vertexPositionX, vertexPositionY, vertexPositionZ, parsProton);
  }

  if(vertexPositionY < -390) {
    DrawEventTrackSideCanvas(totalEventTrackSideCanvas, centerBeamTotal_, centerProton_, protonTrack_, protonTrackRaw_,
              otherTrack_, otherTrackRaw_, parsProton);
    totalEventTrackSideCanvas++;
  }

  // Plot XZ hit position of forward non-central detectors
  double x, y, z;
  line(siYPosForward, parsProton, x, y, z);
  siPosX = x;
  siPosY = siYPosForward;
  siPosZ = z;
  hHitPositionsXZForward->Fill(x, z);
  hHitPositionsXZForwardInd[siDet]->Fill(x, z);

  double vertexToSi = siYPosForward - vertexPositionY;

  double siX, siY, siZ;
  line(siYPosForward, parsProton, siX, siY, siZ);

  TVector3 v1(siX, vertexToSi, siZ);
  TVector3 v2(0, vertexToSi, 0);
  angle = v1.Angle(v2);
  double cosAngle = cos(angle);

  double pathLength = vertexToSi/cosAngle;

  if(pathLength < 0) pathLength = 200.;
  if(pathLength > 700) pathLength = 700.;

  double protonEnergy = protonMethane->AddBack(totalEnergy/1000., pathLength);
  double beamEnergy = protonEnergy*(m1 + m2)*(m1 + m2)/(4.*m1*m2*cosAngle*cosAngle);

  cmEnergy = beamEnergy*m2/(m1 + m2);

  hVertexSiEForward[siDet]->Fill(siEnergy, vertexPositionY);
  hVertexSiEForwardCal[siDet]->Fill(siEnergyCal, vertexPositionY);
  hVertexSiEForwardCalTotal[siDet]->Fill(totalEnergy, vertexPositionY);

  hAngleEForward[siDet]->Fill(siEnergy, angle);
  hAngleEForwardCal[siDet]->Fill(siEnergyCal, angle);
  hAngleEForwardCalTotal[siDet]->Fill(totalEnergy, angle);
  hAngleEForwardProtonEnergy[siDet]->Fill(protonEnergy, angle);

  hVertexAngleForward[siDet]->Fill(vertexPositionY, angle);

  if(siDet == 0 || siDet == 1 || siDet == 9) {
//    hVertexSiETotalRegion3->Fill(totalEnergy, vertexPositionY);
//    s3->Fill(cmEnergy);
    if(angleTotEnergyCut[siDet]->IsInside(protonEnergy, angle)) {
      hVertexSiETotalRegion3->Fill(totalEnergy, vertexPositionY);
      hVertexCMERegion3->Fill(cmEnergy, vertexPositionY);
      s3->Fill(cmEnergy);
    }
  }

  return true;
}

void Spectra::SideVertexFinder(std::vector<mmTrack> centerBeam_, std::vector<mmTrack> protonTrack_, double &vertexPositionX,
                               double &vertexPositionY, double &vertexPositionZ, std::vector<double> &parsProton) {
  vertexPositionX = 0.;
  vertexPositionY = -400.;
  vertexPositionZ = 0.;

  auto *fitProton = new HoughTrack();
  fitProton->AddTrack(protonTrack_, siDet, siQuad);
  double minDistProton = fitProton->FitRestricted();
  parsProton = fitProton->GetPars();
  delete fitProton;

  if(centerBeam_.empty()) {
    double beamX_old = 0.;
    double beamY_old = 250.;
    double beamZ_old = 0.;
    double beamX = beamX_old;
    double beamY = beamY_old;
    double beamZ = beamZ_old;
    double x, y, z;
    line(beamY, parsProton, x, y, z);
    double protonX_old = x;
    double protonX = protonX_old;
    double xDiff_old = beamX - protonX_old;
    double xDiff = xDiff_old;
    double yPos = beamY_old;
    while(yPos > -350.) {
      beamY = yPos;
      line(beamY, parsProton, x, y, z);
      protonX = x;
      xDiff = beamX - protonX;
      if(xDiff_old*xDiff < 0) {
        double m = (beamY - beamY_old)/(xDiff - xDiff_old);
        double b = beamY - m*xDiff;
        vertexPositionX = beamX;
        vertexPositionY = beamY;
        vertexPositionZ = beamZ;
        break;
      }
      beamY_old = beamY;
      protonX_old = protonX;
      xDiff_old = xDiff;
      yPos -= 0.2;
    }
  }
  else {
    // Loop through beam in central region starting from last
    double beamX_old = centerBeam_[0].xPosition;
    double beamY_old = centerBeam_[0].yPosition;
    double beamZ_old = centerBeam_[0].height;
    double beamX = beamX_old;
    double beamY = beamY_old;
    double beamZ = beamZ_old;
    double x, y, z;
    line(beamY, parsProton, x, y, z);
    double protonX_old = x;
    double protonX = protonX_old;
    double xDiff_old = beamX - protonX_old;
    double xDiff = xDiff_old;
    Bool_t foundVertex = false;
    for(auto mm : centerBeam_) {
      beamX = mm.xPosition;
      beamY = mm.yPosition;
      beamZ = mm.height;
      line(beamY, parsProton, x, y, z);
      protonX = x;
      xDiff = beamX - protonX;
      if(xDiff_old*xDiff < 0) {
        double m = (beamY - beamY_old)/(xDiff - xDiff_old);
        double b = beamY - m*xDiff;
        vertexPositionX = beamX;
        vertexPositionY = b;
        vertexPositionZ = beamZ;
        foundVertex = true;
        break;
      }
      beamX_old = beamX;
      beamY_old = beamY;
      beamZ_old = beamZ;
      protonX_old = protonX;
      xDiff_old = xDiff;
    }
    if(!foundVertex) {
      beamX = centerBeam_[0].xPosition;
      beamY = 250.;
      while(beamY > -300.) {
        line(beamY, parsProton, x, y, z);
        protonX = x;
        xDiff = beamX - x;
        if(xDiff*xDiff_old < 0) {
          double m = (beamY - beamY_old)/(xDiff - xDiff_old);
          double b = beamY - m*xDiff;
          vertexPositionX = beamX;
          vertexPositionY = b;
          vertexPositionZ = 0.;
          break;
        }
        beamY_old = beamY;
        protonX_old = protonX;
        xDiff_old = xDiff;
        beamY -= 1.;
      }
    }
  }
}

void Spectra::SideVertexFinderSingleHelp(std::vector<mmTrack> centerBeam_, std::vector<mmTrack> protonTrack_, double &vertexPositionX,
                               double &vertexPositionY, double &vertexPositionZ, std::vector<double> &parsProton) {
  vertexPositionX = 0.;
  vertexPositionY = -400.;
  vertexPositionZ = 0.;

  auto *fitProton = new HoughTrack();
  fitProton->AddTrack(protonTrack_, siDet, siQuad);
  double minDistProton = fitProton->FitRestrictedSingleHelp();
  parsProton = fitProton->GetPars();
  delete fitProton;

  if(centerBeam_.empty()) {
    double beamX_old = 0.;
    double beamY_old = 250.;
    double beamZ_old = 0.;
    double beamX = beamX_old;
    double beamY = beamY_old;
    double beamZ = beamZ_old;
    double x, y, z;
    line(beamY, parsProton, x, y, z);
    double protonX_old = x;
    double protonX = protonX_old;
    double xDiff_old = beamX - protonX_old;
    double xDiff = xDiff_old;
    double yPos = beamY_old;
    while(yPos > -350.) {
      beamY = yPos;
      line(beamY, parsProton, x, y, z);
      protonX = x;
      xDiff = beamX - protonX;
      if(xDiff_old*xDiff < 0) {
        double m = (beamY - beamY_old)/(xDiff - xDiff_old);
        double b = beamY - m*xDiff;
        vertexPositionX = beamX;
        vertexPositionY = beamY;
        vertexPositionZ = beamZ;
        break;
      }
      beamY_old = beamY;
      protonX_old = protonX;
      xDiff_old = xDiff;
      yPos -= 0.2;
    }
  }
  else {
    // Loop through beam in central region starting from last
    double beamX_old = centerBeam_[0].xPosition;
    double beamY_old = centerBeam_[0].yPosition;
    double beamZ_old = centerBeam_[0].height;
    double beamX = beamX_old;
    double beamY = beamY_old;
    double beamZ = beamZ_old;
    double x, y, z;
    line(beamY, parsProton, x, y, z);
    double protonX_old = x;
    double protonX = protonX_old;
    double xDiff_old = beamX - protonX_old;
    double xDiff = xDiff_old;
    Bool_t foundVertex = false;
    for(auto mm : centerBeam_) {
      beamX = mm.xPosition;
      beamY = mm.yPosition;
      beamZ = mm.height;
      line(beamY, parsProton, x, y, z);
      protonX = x;
      xDiff = beamX - protonX;
      if(xDiff_old*xDiff < 0) {
        double m = (beamY - beamY_old)/(xDiff - xDiff_old);
        double b = beamY - m*xDiff;
        vertexPositionX = beamX;
        vertexPositionY = b;
        vertexPositionZ = beamZ;
        foundVertex = true;
        break;
      }
      beamX_old = beamX;
      beamY_old = beamY;
      beamZ_old = beamZ;
      protonX_old = protonX;
      xDiff_old = xDiff;
    }
    if(!foundVertex) {
      beamX = centerBeam_[0].xPosition;
      beamY = 250.;
      while(beamY > -300.) {
        line(beamY, parsProton, x, y, z);
        protonX = x;
        xDiff = beamX - x;
        if(xDiff*xDiff_old < 0) {
          double m = (beamY - beamY_old)/(xDiff - xDiff_old);
          double b = beamY - m*xDiff;
          vertexPositionX = beamX;
          vertexPositionY = b;
          vertexPositionZ = 0.;
          break;
        }
        beamY_old = beamY;
        protonX_old = protonX;
        xDiff_old = xDiff;
        beamY -= 1.;
      }
    }
  }
  if(vertexPositionY < -390) {
    double m_xcomponent = fabs(parsProton[1]);
    double xAngle = atan(m_xcomponent);
    double yDist = siXPosForward[siDet][siQuad]/fabs(tan(xAngle));
    if(275. - yDist < -400) {
      vertexPositionY = -400.;
    }
    else {
      vertexPositionY = 275. - yDist;
    }
  }
}

void Spectra::SideVertexFinderHelp(std::vector<mmTrack> centerBeam_, std::vector<mmTrack> protonTrack_, double &vertexPositionX,
                               double &vertexPositionY, double &vertexPositionZ, std::vector<double> &parsProton) {
  vertexPositionX = 0.;
  vertexPositionY = -400.;
  vertexPositionZ = 0.;

  auto *fitProton = new HoughTrack();
  fitProton->AddTrack(protonTrack_, siDet, siQuad);
  double minDistProton = fitProton->FitRestrictedHelp();
  parsProton = fitProton->GetPars();
  delete fitProton;

  if(centerBeam_.empty()) {
    double beamX_old = 0.;
    double beamY_old = 250.;
    double beamZ_old = 0.;
    double beamX = beamX_old;
    double beamY = beamY_old;
    double beamZ = beamZ_old;
    double x, y, z;
    line(beamY, parsProton, x, y, z);
    double protonX_old = x;
    double protonX = protonX_old;
    double xDiff_old = beamX - protonX_old;
    double xDiff = xDiff_old;
    double yPos = beamY_old;
    while(yPos > -350.) {
      beamY = yPos;
      line(beamY, parsProton, x, y, z);
      protonX = x;
      xDiff = beamX - protonX;
      if(xDiff_old*xDiff < 0) {
        double m = (beamY - beamY_old)/(xDiff - xDiff_old);
        double b = beamY - m*xDiff;
        vertexPositionX = beamX;
        vertexPositionY = beamY;
        vertexPositionZ = beamZ;
        break;
      }
      beamY_old = beamY;
      protonX_old = protonX;
      xDiff_old = xDiff;
      yPos -= 0.2;
    }
  }
  else {
    // Loop through beam in central region starting from last
    double beamX_old = centerBeam_[0].xPosition;
    double beamY_old = centerBeam_[0].yPosition;
    double beamZ_old = centerBeam_[0].height;
    double beamX = beamX_old;
    double beamY = beamY_old;
    double beamZ = beamZ_old;
    double x, y, z;
    line(beamY, parsProton, x, y, z);
    double protonX_old = x;
    double protonX = protonX_old;
    double xDiff_old = beamX - protonX_old;
    double xDiff = xDiff_old;
    Bool_t foundVertex = false;
    for(auto mm : centerBeam_) {
      beamX = mm.xPosition;
      beamY = mm.yPosition;
      beamZ = mm.height;
      line(beamY, parsProton, x, y, z);
      protonX = x;
      xDiff = beamX - protonX;
      if(xDiff_old*xDiff < 0) {
        double m = (beamY - beamY_old)/(xDiff - xDiff_old);
        double b = beamY - m*xDiff;
        vertexPositionX = beamX;
        vertexPositionY = b;
        vertexPositionZ = beamZ;
        foundVertex = true;
        break;
      }
      beamX_old = beamX;
      beamY_old = beamY;
      beamZ_old = beamZ;
      protonX_old = protonX;
      xDiff_old = xDiff;
    }
    if(!foundVertex) {
      beamX = centerBeam_[0].xPosition;
      beamY = 250.;
      while(beamY > -300.) {
        line(beamY, parsProton, x, y, z);
        protonX = x;
        xDiff = beamX - x;
        if(xDiff*xDiff_old < 0) {
          double m = (beamY - beamY_old)/(xDiff - xDiff_old);
          double b = beamY - m*xDiff;
          vertexPositionX = beamX;
          vertexPositionY = b;
          vertexPositionZ = 0.;
          break;
        }
        beamY_old = beamY;
        protonX_old = protonX;
        xDiff_old = xDiff;
        beamY -= 1.;
      }
    }
  }
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

void Spectra::CorrectCenterEnergy(std::vector<mmTrack> &centerBeam_, std::vector<double> parsBeam, int lastRow) {
  double x, y, z;

  // Correct energy
  // for(auto mm : centerBeam_) {
  for(uint i = 0; i < centerBeam_.size(); i++) {
    if(centerBeam_[i].row > lastRow) continue;

    // if(centerBeam_[i].total != 1) continue;
    line(centerBeam_[i].yPosition, parsBeam, x, y, z);
    double distToMesh = 6.54 + centerBeam_[i].height/10.;
    double gasSigma = sqrt(distToMesh)*gasPositionResolution*10.;
    // double gausSigma = gasSigma;
    double gausSigma = sqrt(gasSigma*gasSigma + 1.5*1.5);
    double gausMean = x;

    // Find pad where the hit fell into
    double padFirst = 0;
    double padSecond = 0;

    if(centerBeam_[i].total == 1) {
      for(int j = 0; j < 6; j++) {
        if((centerBeam_[i].xPosition > mmColumnSize[j].first) &&
            (centerBeam_[i].xPosition < mmColumnSize[j].second)) {
          padFirst = mmColumnSize[j].first;
          padSecond = mmColumnSize[j].second;
          break;
        }
      }
//      std::cout << entry << '\t' << centerBeam_[i].row << '\t' << centerBeam_[i].xPosition << '\t' << padFirst << '\t' << padSecond << std::endl;
    }
    else if(centerBeam_[i].total == 2) {
      for(int j = 1; j < 6; j++) {
        if((centerBeam_[i].xPosition > mmColumnSize[j - 1].first) &&
           (centerBeam_[i].xPosition < mmColumnSize[j].second) &&
           fabs(centerBeam_[i].xPosition - mmColumnSize[j - 1].first) < 5) {
          padFirst = mmColumnSize[j - 1].first;
          padSecond = mmColumnSize[j].second;
          break;
        }
      }
      // std::cout << entry << '\t' << centerBeam_[i].xPosition << '\t' << padFirst << '\t' << padSecond << std::endl;
    }
    else continue;

    double gausCDFFirst = GaussianCDF(padFirst, gausMean, gausSigma);
    double gausCDFSecond = GaussianCDF(padSecond, gausMean, gausSigma);
    double gausPercentage = gausCDFSecond - gausCDFFirst;

    if(gausPercentage > 0.5) {
      centerBeam_[i].energy = centerBeam_[i].energy/gausPercentage;
    }
    else {
//      centerBeam_[i].energy = 0.;
      centerBeam_[i].energy = centerBeam_[i].energy/0.5;
    }
  }
}

double Spectra::GaussianCDF(double x, double mean, double sigma) {

  return 0.5*(1. + erf((x - mean)/(sqrt(2.)*sigma)));
}

void Spectra::FindMaxCentralEnergy(std::vector<mmTrack> centerMatched_, int &maxEnergyRow, double &maxEnergy,
                                   double &averageMaxEnergy, double &maxDeriv) {
  maxEnergy = -1000.;
  int maxRow = -1;
  int maxRowBin = 0;
  double totalEnergy = 0.;

  for(uint i = 0; i < centerMatched_.size(); i++) {
    if(centerMatched_[i].energy > maxEnergy) {
      maxEnergy = centerMatched_[i].energy;
      maxRow = centerMatched_[i].row;
      maxRowBin = static_cast<int>(i);
    }
  }

  int minRowBin;
  if(maxRowBin < 5) {
    minRowBin = 0;
  }
  else {
    minRowBin = maxRowBin - 4;
  }

  double derivative;
  if(maxRowBin > 1) {
    derivative = (centerMatched_[maxRowBin].energy - centerMatched_[maxRowBin - 2].energy)/2.;
  }
  else {
    derivative = 0.;
  }

  int count = 0;
  for(int i = minRowBin; i < maxRowBin + 1; i++) {
    totalEnergy += centerMatched_[i].energy;
    count++;
  }

  maxEnergyRow = maxRow;
  averageMaxEnergy = totalEnergy/static_cast<double>(count);
  maxDeriv = derivative;
}

std::vector<mmTrack> Spectra::GetRunningEnergyAverageThree(std::vector<mmTrack> centerMatched_) {
  std::vector<mmTrack> newTrack_;
  newTrack_.push_back(centerMatched_[0]);
  for(uint i = 1; i < centerMatched_.size() - 2; i++) {
    double avgEnergy = (centerMatched_[i - 1].energy + centerMatched_[i].energy + centerMatched_[i + 1].energy)/3.;
    newTrack_.push_back(centerMatched_[i]);
    newTrack_[i].energy = avgEnergy;
  }
  newTrack_.push_back(centerMatched_[centerMatched_.size() - 1]);

  return newTrack_;
}

std::vector<mmTrack> Spectra::GetRunningEnergyAverageFive(std::vector<mmTrack> centerMatched_) {
  std::vector<mmTrack> newTrack_;
  newTrack_.push_back(centerMatched_[0]);
  newTrack_.push_back(centerMatched_[1]);
  for(uint i = 2; i < centerMatched_.size() - 3; i++) {
    double avgEnergy = (centerMatched_[i - 2].energy + centerMatched_[i - 1].energy +
        centerMatched_[i].energy + centerMatched_[i + 1].energy + centerMatched_[i + 2].energy)/5.;
    newTrack_.push_back(centerMatched_[i]);
    newTrack_[i].energy = avgEnergy;
  }
  newTrack_.push_back(centerMatched_[centerMatched_.size() - 2]);
  newTrack_.push_back(centerMatched_[centerMatched_.size() - 1]);

  return newTrack_;
}

bool Spectra::CenterOnlyOneColumn(std::vector<mmTrack> centerMatched_) {
  bool singleColumn = true;

  for(auto mm : centerMatched_) {
    if(mm.total > 1) singleColumn = false;
  }

  return singleColumn;
}

std::vector<centerDeriv> Spectra::CenterEnergyThreePointDeriv(std::vector<mmTrack> center_) {
  std::vector<centerDeriv> mmDeriv;
  for(uint i = 1; i < center_.size() - 1; i++) {
    double deriv = (center_[i + 1].energy - center_[i - 1].energy)/2.;
    centerDeriv hit = {center_[i].row, deriv};
    mmDeriv.push_back(hit);
  }
  return mmDeriv;
}

std::vector<centerDeriv> Spectra::CenterEnergyFivePointDeriv(std::vector<mmTrack> center_) {
  std::vector<centerDeriv> mmDeriv;
  for(uint i = 2; i < center_.size() - 2; i++) {
    double deriv = (-center_[i + 2].energy + 8.*center_[i + 1].energy -
                      8.*center_[i - 1].energy + center_[i - 2].energy)/12.;
    centerDeriv hit = {center_[i].row, deriv};
    mmDeriv.push_back(hit);
  }
  return mmDeriv;
}

std::pair<int, int> Spectra::CenterGetDerivMax(std::vector<centerDeriv> threePoint_, std::vector<centerDeriv> fivePoint_) {
  int maxThreePointRow, maxFivePointRow;

  std::pair<int, int> pair;

  double maxDeriv = -1000.;
  for(auto mm : threePoint_) {
    if(mm.deriv > maxDeriv) {
      maxDeriv = mm.deriv;
      maxThreePointRow = mm.row;
    }
  }

  maxDeriv = -1000.;
  for(auto mm : fivePoint_) {
    if(mm.deriv > maxDeriv) {
      maxDeriv = mm.deriv;
      maxFivePointRow = mm.row;
    }
  }

  pair = std::make_pair(maxThreePointRow, maxFivePointRow);

  return pair;
}

void Spectra::CleanLeftStrip(std::vector<mmChainStrip> &cleaned, std::vector<mmChainStrip> left) {
  cleaned.clear();

  if(left.size() < 2) return;

  for(uint i = 0; i < left.size(); i++) {
    bool validStrip = false;

    std::vector<int>::iterator it = std::find(leftStripHits.begin(), leftStripHits.end(), left[i].row);
    if(it == leftStripHits.end()) continue;
    int listIndex = std::distance(leftStripHits.begin(), it);

    if(!validStrip) return;

    cleaned.push_back(left[i]);
    // std::cout << mm.row << '\t' << listIndex << std::endl;
  }

  return;
}

void Spectra::CleanRightStrip(std::vector<mmChainStrip> &cleaned, std::vector<mmChainStrip> right) {
  cleaned.clear();

  if(right.size() < 2) return;

  for(auto mm : right) {
    std::vector<int>::iterator it = std::find(rightStripHits.begin(), rightStripHits.end(), mm.row);
    if(it == rightStripHits.end()) continue;
    int listIndex = std::distance(rightStripHits.begin(), it);
    // std::cout << mm.row << '\t' << listIndex << std::endl;
  }

  return;
}

void Spectra::ChainStripMatch(std::vector<mmTrack> &chainStripMatched, std::vector<mmTrack> &chainStripRaw,
                              std::vector<mmChainStrip> chain_,
                              std::vector<mmChainStrip> strip_, bool leftSide, double siTime) {
  // Function that matches chains and strips together
  // If number of different time buckets < 4, uses box method
  // If > 3, uses the same time to match

  chainStripMatched.clear();
  chainStripRaw.clear();
  std::vector<mmTrack> totalTime0;
  std::vector<mmTrack> totalTime1;
  ChainStripMatchingTime(totalTime0, chain_, strip_, leftSide, siTime, 0);
  ChainStripMatchingTime(totalTime1, chain_, strip_, leftSide, siTime, 1);

  unsigned long numTimeBuckets = ChainStripTime0NumTimeBuckets(totalTime0);
  if(numTimeBuckets < 4) {
    ChainStripMatchingBoxTime0(chainStripMatched, totalTime0);
    // ChainStripMatchingBoxTime0(chainStripMatched, totalTime1);
  }
  else {
    ChainStripMatchingTime(chainStripMatched, chain_, strip_, leftSide, siTime, 0);
  }

  chainStripRaw = totalTime0;
  // chainStripRaw = totalTime1;
}

void Spectra::ChainStripMatchingTime(std::vector<mmTrack> &chainStripMatched, std::vector<mmChainStrip> chain,
                                     std::vector<mmChainStrip> strip, bool leftSide, double siTime,
                                     int timeWindow) {
  // Function to match strips and chains. Chains and strips are matched if their time is within the timeWindow parameter
  // The timeWindow parameter is the window of time buckets to look for using the timeResolution parameter
  // which is the ns -> bucket
  // Position transform = 10.5 + 1.75/2 + 1.75*chain (chain starting at 0)
  // Row transform = strip*2
  // Row position transform = Row * rowConversion + rowConversionOffset

  for(auto chainVec : chain) {
    int timeChain = chainVec.timeBucket;
    // Loop over strips and find strips when the time is within the timeWindow
    for(auto stripVec : strip) {
      int timeStrip = stripVec.timeBucket;
      if((timeStrip > timeChain - timeWindow - 1) &&
          (timeStrip < timeChain + timeWindow + 1)) {
        mmTrack hit = {0, 0., 0., 0., 0., 0., 0, 0};
        double position = 10.5 + 1.75/2 + 1.75*chainVec.row;
        if(leftSide) position = -position;
        int row = stripVec.row*2;
        hit.row = row;
        hit.xPosition = position;
        hit.yPosition = hit.row*rowConversion + rowConversionOffset;
        hit.time = stripVec.time - siTime;
        hit.energy = stripVec.energy;
        hit.height = heightOffset - hit.time*driftVelocity;
        hit.timeBucket = timeStrip;
        hit.total = 1;
        chainStripMatched.push_back(hit);
      }
    }
  }
}

size_t Spectra::ChainStripTime0NumTimeBuckets(std::vector<mmTrack> matched) {
  // Function that finds the number of different time buckets for the time0 algorithm
  // This is used if you want to change the strip/chain matching algorithm (if all on the same plane)

  std::map<int, int> timeMap;

  for(auto mm : matched) {
    timeMap[mm.timeBucket] = 1;
  }

  return timeMap.size();
}

void Spectra::ChainStripMatchingBoxTime0(std::vector<mmTrack> &chainStripMatched, std::vector<mmTrack> time0) {
  // Function to quickly get the box points after matching using the same time
  // Finds the points closest and furthest from the origin for box
  chainStripMatched.clear();
  if(time0.empty()) return;

  mmTrack closestHit;
  double closestDistance = 5000.;
  mmTrack furthestHit;
  double furthestDistance = 0.;

  for(auto mm : time0) {
    double distance = sqrt(mm.xPosition*mm.xPosition + mm.yPosition*mm.yPosition);
    if(distance < closestDistance) {
      closestHit = mm;
      closestDistance = distance;
    }
    if(distance > furthestDistance) {
      furthestHit = mm;
      furthestDistance = distance;
    }
  }

  chainStripMatched.push_back(closestHit);
  chainStripMatched.push_back(furthestHit);
}

double Spectra::ChainStripSize(std::vector<mmTrack> chainStripMatched) {
  double initX = chainStripMatched[0].xPosition;
  double initY = chainStripMatched[0].yPosition;
  double finalX = chainStripMatched[chainStripMatched.size() - 1].xPosition;
  double finalY = chainStripMatched[chainStripMatched.size() - 1].yPosition;

  double diffX = fabs(finalX - initX);
  double diffY = fabs(finalY - initY);

  double distance = sqrt(diffX*diffX + diffY*diffY);

  return distance;
}

void Spectra::CleanSideTracks(std::vector<mmTrack> &cleaned_, std::vector<mmTrack> &cleanedRaw_, std::vector<mmTrack> trackRaw_,
                              double angle, int &changes) {
  // Get dE average and compare, remove hits with difference between average > 400
  changes = 0;
  double maxDiff = 400.;

  double averageEnergy = 0.;
  int totalRows = 0;
  for(auto mm : trackRaw_) {
    if(mm.energy > 2700) continue;
    averageEnergy += mm.energy*cos(angle);
    totalRows++;
  }
  averageEnergy /= static_cast<double>(totalRows);

  cleaned_.clear();
  cleanedRaw_.clear();
  for(auto mm : trackRaw_) {
    if(fabs(mm.energy*cos(angle) - averageEnergy) > maxDiff) {
      changes++;
      continue;
    }
    cleanedRaw_.push_back(mm);
  }

  if(cleanedRaw_.size() < 2) return;

  unsigned long numTimeBuckets = ChainStripTime0NumTimeBuckets(cleanedRaw_);
  if(numTimeBuckets < 4) {
    ChainStripMatchingBoxTime0(cleaned_, cleanedRaw_);
  }
  else {
    cleaned_ = cleanedRaw_;
  }
  return;
}

void Spectra::DivideTargetThickness(TH1F *f) {
  int i_size = f->GetSize();

  TAxis *x_axis = f->GetXaxis();
  for(int i = 1; i <= i_size - 1; i++) {
    double binLowEdge = x_axis->GetBinLowEdge(i);
    double binUpEdge = x_axis->GetBinUpEdge(i);
    if(binLowEdge == 0.) binLowEdge += 0.001;
    binLowEdge *= (m1 + m2)/m2;
    binUpEdge *= (m1 + m2)/m2;
    double binContent = f->GetBinContent(i);
    double binError = f->GetBinError(i);
    double delta_x = boronMethane->CalcRange(binUpEdge, binLowEdge);
    delta_x /= 10.;
    double molarMassMethane = 0.01604;
    double factor = 4.e-27*density*delta_x*TMath::Na()/molarMassMethane;
    binContent /= factor;
    binError /= factor;
    f->SetBinContent(i, binContent);
    f->SetBinError(i, binError);
  }
}

void Spectra::ReadSolidAngle() {

  double var1, var2, var3, var4;

  std::ifstream inReg1("csReg1_4.out");
  assert(inReg1.is_open());
  std::vector<double> cmEnergyReg1_, solidAngleReg1_, labAngleReg1_, cmAngleReg1_;
  while(inReg1 >> var1 >> var2 >> var3 >> var4) {
    cmEnergyReg1_.push_back(var1);
    solidAngleReg1_.push_back(var2);
    labAngleReg1_.push_back(var3);
    cmAngleReg1_.push_back(var4);
  }
  reg1SA.SetPoints(cmEnergyReg1_, solidAngleReg1_);
  reg1CMAngle.SetPoints(cmEnergyReg1_, cmAngleReg1_);
  inReg1.close();

  std::ifstream inReg3("csReg3.out");
  assert(inReg3.is_open());
  std::vector<double> cmEnergyReg3_, solidAngleReg3_, labAngleReg3_, cmAngleReg3_;
  while(inReg3 >> var1 >> var2 >> var3 >> var4) {
    cmEnergyReg3_.push_back(var1);
    solidAngleReg3_.push_back(var2);
    labAngleReg3_.push_back(var3);
    cmAngleReg3_.push_back(var4);
  }
  reg3SA.SetPoints(cmEnergyReg3_, solidAngleReg3_);
  reg3CMAngle.SetPoints(cmEnergyReg3_, cmAngleReg3_);
  inReg3.close();
}

void Spectra::SolidAngle(TH1F *f, int region) {
  int i_size = f->GetSize();
  TAxis *x_axis = f->GetXaxis();

  for(int i = 1; i <= i_size - 1; i++) {
    double binCenter = x_axis->GetBinCenter(i);
    double binContent = f->GetBinContent(i);
    double binError = f->GetBinError(i);

    double simCS = 0;
    if(region == 1) {
      simCS = reg1SA(binCenter);
    }
    if(region == 3) {
      simCS = reg3SA(binCenter);
    }

    if(simCS < 0.001) {
      f->SetBinContent(i, 0.);
      f->SetBinError(i, 0.);
    }
    else {
      f->SetBinContent(i, binContent/simCS);
      f->SetBinError(i, binError/simCS);
    }
  }
}

void Spectra::WriteSpectrumToFile(TH1F *f, int region) {
  FILE* spectrumFile = fopen(Form("spectrum_reg%d.out", region), "w");

  int i_size = f->GetSize();
  TAxis *x_axis = f->GetXaxis();

  for(int i = 1; i <= i_size - 1; i++) {
    double binCenter = x_axis->GetBinCenter(i);
    double binContent = f->GetBinContent(i);
    double binError = f->GetBinError(i);

    double cmAngle = reg3CMAngle(binCenter);

    fprintf(spectrumFile, "%f %f %f %f\n", binCenter, binContent, binError, cmAngle);
  }

  fflush(spectrumFile);
  fclose(spectrumFile);
}

void Spectra::GetMinMaxD(std::vector<mmTrack> initPoints, int &minXY, int &maxXY, int &minYZ, int &maxYZ) {
  std::vector<xy> pointsXY;
  std::vector<yz> pointsYZ;
  for(uint i = 0; i < initPoints.size(); i++) {
    xy initPointsXY = {initPoints[i].xPosition, initPoints[i].yPosition};
    yz initPointsYZ = {initPoints[i].yPosition, initPoints[i].height};
    pointsXY.push_back(initPointsXY);
    pointsYZ.push_back(initPointsYZ);
  }

  double minDHistogramXY = 1000;
  double maxDHistogramXY = -1000;
  double minDHistogramYZ = 1000;
  double maxDHistogramYZ = -1000;
  for(int j = 0; j < 180; j++) {
    double cosj = cos(j*M_PI/180.);
    double sinj = sin(j*M_PI/180.);
    for(uint i = 0; i < pointsXY.size(); i++) {
      double d = pointsXY[i].x*cosj + pointsXY[i].y*sinj;
      if(d < minDHistogramXY) minDHistogramXY = d;
      if(d > maxDHistogramXY) maxDHistogramXY = d;
    }
    for(uint i = 0; i < pointsYZ.size(); i++) {
      double d = pointsYZ[i].y*cosj + pointsYZ[i].z*sinj;
      if(d < minDHistogramYZ) minDHistogramYZ = d;
      if(d > maxDHistogramYZ) maxDHistogramYZ = d;
    }
  }

  minXY = static_cast<int>(minDHistogramXY);
  maxXY = static_cast<int>(maxDHistogramXY);
  minYZ = static_cast<int>(minDHistogramYZ);
  maxYZ = static_cast<int>(maxDHistogramYZ);
}

void Spectra::VisualizeHough(std::vector<mmTrack> initPoints, TH2I* fXY, TH2I* fYZ) {
  int nBinsX = 500;
  int nBinsY = 500;
  std::vector<xy> pointsXY;
  std::vector<yz> pointsYZ;
  for(uint i = 0; i < initPoints.size(); i++) {
    xy initPointsXY = {initPoints[i].xPosition, initPoints[i].yPosition};
    yz initPointsYZ = {initPoints[i].yPosition, initPoints[i].height};
    pointsXY.push_back(initPointsXY);
    pointsYZ.push_back(initPointsYZ);
  }

  double thetaStepXY = 180./static_cast<double>(nBinsX);

  // Fill Hough Matrix
  for(double j = 0; j < 180; j += thetaStepXY) {
    if(j > 89 && j < 91) continue;
    double cosj = cos(j*M_PI/180.);
    double sinj = sin(j*M_PI/180.);
    for(uint i = 0; i < pointsXY.size(); i++) {
      double d = pointsXY[i].x*cosj + pointsXY[i].y*sinj;
      fXY->Fill(j, d);
    }
  }
}

void Spectra::GetMinMaxDRestricted(std::vector<mmTrack> initPoints, int &minXY, int &maxXY, int &minYZ, int &maxYZ, int siDet) {
  int minAngle = 0;
  int maxAngle = 0;
  if(siDet < 4) {
    minAngle = 91;
    maxAngle = 179;
  }
  else if(siDet > 5 && siDet < 10) {
    minAngle = 1;
    maxAngle = 89;
  }

  std::vector<xy> pointsXY;
  std::vector<yz> pointsYZ;
  for(uint i = 0; i < initPoints.size(); i++) {
    xy initPointsXY = {initPoints[i].xPosition, initPoints[i].yPosition};
    yz initPointsYZ = {initPoints[i].yPosition, initPoints[i].height};
    pointsXY.push_back(initPointsXY);
    pointsYZ.push_back(initPointsYZ);
  }

  double minDHistogramXY = 1000;
  double maxDHistogramXY = -1000;
  double minDHistogramYZ = 1000;
  double maxDHistogramYZ = -1000;
  for(int j = minAngle; j < maxAngle; j++) {
    double cosj = cos(j*M_PI/180.);
    double sinj = sin(j*M_PI/180.);
    for(uint i = 0; i < pointsXY.size(); i++) {
      double d = pointsXY[i].x*cosj + pointsXY[i].y*sinj;
      if(d < minDHistogramXY) minDHistogramXY = d;
      if(d > maxDHistogramXY) maxDHistogramXY = d;
    }
    for(uint i = 0; i < pointsYZ.size(); i++) {
      double d = pointsYZ[i].y*cosj + pointsYZ[i].z*sinj;
      if(d < minDHistogramYZ) minDHistogramYZ = d;
      if(d > maxDHistogramYZ) maxDHistogramYZ = d;
    }
  }

  minXY = static_cast<int>(minDHistogramXY);
  maxXY = static_cast<int>(maxDHistogramXY);
  minYZ = static_cast<int>(minDHistogramYZ);
  maxYZ = static_cast<int>(maxDHistogramYZ);
}

void Spectra::VisualizeHoughRestricted(std::vector<mmTrack> initPoints, TH2I* fXY, TH2I* fYZ, int siDet) {
  int minAngle = 0;
  int maxAngle = 180;
  if(siDet < 4) {
    minAngle = 91;
    maxAngle = 179;
  }
  else if(siDet > 5 && siDet < 10) {
    minAngle = 1;
    maxAngle = 89;
  }

  int nBinsX = 360;
  int nBinsY = 360;
  std::vector<xy> pointsXY;
  std::vector<yz> pointsYZ;
  for(auto point : initPoints) {
    xy initPointsXY = {point.xPosition, point.yPosition};
    yz initPointsYZ = {point.yPosition, point.height};
    pointsXY.push_back(initPointsXY);
    pointsYZ.push_back(initPointsYZ);
  }

  double thetaStepXY = 180./static_cast<double>(nBinsX);

  // Fill Hough Matrix
  for(double j = minAngle; j < maxAngle; j += thetaStepXY) {
    if(j > 89 && j < 91) continue;
    double cosj = cos(j*M_PI/180.);
    double sinj = sin(j*M_PI/180.);
    for(auto point : pointsXY) {
      double d = point.x*cosj + point.y*sinj;
      fXY->Fill(j, d);
    }
  }
}

void Spectra::GetHoughStdDevXYRestricted(std::vector<mmTrack> initPoints, std::vector<double> &angle_, std::vector<double> &stdDev_, int siDet) {

  angle_.clear();
  stdDev_.clear();

  int minAngle = 0;
  int maxAngle = 180;
  if(siDet < 4) {
    minAngle = 91;
    maxAngle = 179;
  }
  else if(siDet > 5 && siDet < 10) {
    minAngle = 1;
    maxAngle = 89;
  }

  int nBinsX = 360;
  int nBinsY = 360;
  std::vector<xy> pointsXY;
  std::vector<yz> pointsYZ;
  for(uint i = 0; i < initPoints.size(); i++) {
    xy initPointsXY = {initPoints[i].xPosition, initPoints[i].yPosition};
    yz initPointsYZ = {initPoints[i].yPosition, initPoints[i].height};
    pointsXY.push_back(initPointsXY);
    pointsYZ.push_back(initPointsYZ);
  }

  double thetaStepXY = 180./static_cast<double>(nBinsX);

  // Fill Hough Matrix
  for(double j = minAngle; j < maxAngle; j += thetaStepXY) {
    if(j > 89 && j < 91) continue;
    double cosj = cos(j*M_PI/180.);
    double sinj = sin(j*M_PI/180.);
    std::vector<double> dVector;
    double mean = 0.;
    for(uint i = 0; i < pointsXY.size(); i++) {
      double d = pointsXY[i].x*cosj + pointsXY[i].y*sinj;
      mean += d;
      dVector.push_back(d);
    }
    if(dVector.empty()) continue;

    mean /= static_cast<double>(pointsXY.size());
    double stdDev = 0.;
    for(uint i = 0; i < pointsXY.size(); i++) {
      double d = pointsXY[i].x*cosj + pointsXY[i].y*sinj;
      stdDev += (d - mean)*(d - mean);
    }
    stdDev /= static_cast<double>(pointsXY.size() - 1.);

    angle_.push_back(j);
    stdDev_.push_back(stdDev);

  }

}
