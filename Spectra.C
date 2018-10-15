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
  TString PathToFiles = "/hd3/research/data/run0817a/rootM2R-WaveformReduced/run";

  // Mac Laptop
  // TString PathToFiles = "/Users/joshhooker/Desktop/data/run0817a/run";

  // Linux Laptop
  // TString PathToFiles = "/home/joshhooker/Desktop/data/run0817a/run";

  // Alpha source test in gas
  chain->Add(PathToFiles+"004.root");

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
//  chain->Add(PathToFiles + "175.root");
//  chain->Add(PathToFiles + "178.root");
//  chain->Add(PathToFiles + "180.root");
//  chain->Add(PathToFiles + "181.root");
//  chain->Add(PathToFiles + "182.root");
//  chain->Add(PathToFiles + "183.root");
//  chain->Add(PathToFiles + "184.root");
//  chain->Add(PathToFiles + "185.root");
//  chain->Add(PathToFiles + "186.root");
//  chain->Add(PathToFiles + "187.root");
//  chain->Add(PathToFiles + "188.root");
//  chain->Add(PathToFiles + "189.root");
//  chain->Add(PathToFiles + "190.root");
//  chain->Add(PathToFiles + "191.root");
//  chain->Add(PathToFiles + "192.root");
//  chain->Add(PathToFiles + "193.root");
//  chain->Add(PathToFiles + "195.root");
//  chain->Add(PathToFiles + "196.root");
//  chain->Add(PathToFiles + "197.root");
//  chain->Add(PathToFiles + "198.root");
//  chain->Add(PathToFiles + "199.root");
//  chain->Add(PathToFiles + "200.root");
//  chain->Add(PathToFiles + "201.root");
//  chain->Add(PathToFiles + "202.root");
//  chain->Add(PathToFiles + "203.root");
//  chain->Add(PathToFiles + "204.root");
//  chain->Add(PathToFiles + "205.root");
//  chain->Add(PathToFiles + "206.root");
//  chain->Add(PathToFiles + "207.root");
//  chain->Add(PathToFiles + "208.root");
//  chain->Add(PathToFiles + "210.root");
//  chain->Add(PathToFiles + "211.root");
//  chain->Add(PathToFiles + "212.root");
//  chain->Add(PathToFiles + "213.root");
//  chain->Add(PathToFiles + "214.root");
//  chain->Add(PathToFiles + "215.root");
//  chain->Add(PathToFiles + "216.root");
//  chain->Add(PathToFiles + "217.root");
//  chain->Add(PathToFiles + "218.root");
//  chain->Add(PathToFiles + "219.root");
//  chain->Add(PathToFiles + "220.root");
//  chain->Add(PathToFiles + "221.root");
//  chain->Add(PathToFiles + "223.root");
//  chain->Add(PathToFiles + "224.root");

  // Alpha source test in vacuum

  // Alpha source test in gas
  // chain->Add(PathToFiles + "232.root");

  return chain;
}

void Spectra::Loop() {

  if (fChain == nullptr) return;

  InitChannelMap();
  InitHistograms();
  InitVariables();

  InitSiEForwardCalibration();
  InitCsIECalibration();

  InitCentralPadGainMatch();
  InitAverageBeamEnergy();

  InitTree();

  Long64_t nentries = fChain->GetEntriesFast();

  printf("Starting Main Loop\n");

  Double_t totalEnergy[128][6];
  Int_t totalFired[128][6];
  for(UInt_t i = 0; i < 128; i++) {
    for(UInt_t j = 0; j < 6; j++) {
      totalFired[i][j] = 0;
      totalEnergy[i][j] = 0.;
    }
  }

  Long64_t nbytes = 0, nb = 0;
  for(Long64_t jentry = 0; jentry < nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if(ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;

    if(jentry != 0 && jentry % 10000 == 0) printf("Processed %lld events\n", jentry);

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
    Double_t icE = 0.;
    Double_t icT = 0.;

    for(Int_t i = 0; i < mmMul; i++) {
      if(mmChan[i] == 11 || mmChan[i] == 22 || mmChan[i] == 45 || mmChan[i] == 56) continue; // Skip FPN Channels

      if(mmEnergy[i] < 0) continue;

      // Micromegas
      if(mmCobo[i] == 0) {
        // Asad0 - All Center Pads
        if(mmAsad[i] == 0) {
          // Aget0
          if(mmAget[i] == 0) {
            std::pair<Int_t, Int_t> pad = MM_Map_Asad0_Aget0[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, mmEnergy[i], mmTime[i], mmPa[i][3]/mmPa[i][1]};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, mmEnergy[i]*scale[pad.first][pad.second], mmTime[i],
                                     mmPa[i][3]/mmPa[i][1]};
            mmCenterMatched_.push_back(mmHitMatched);
          }
          // Aget1
          else if(mmAget[i] == 1) {
            std::pair<Int_t, Int_t> pad = MM_Map_Asad0_Aget1[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, mmEnergy[i], mmTime[i], mmPa[i][3]/mmPa[i][1]};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, mmEnergy[i]*scale[pad.first][pad.second], mmTime[i],
                                     mmPa[i][3]/mmPa[i][1]};
            mmCenterMatched_.push_back(mmHitMatched);
          }
          // Aget2
          else if(mmAget[i] == 2) {
            std::pair<Int_t, Int_t> pad = MM_Map_Asad0_Aget2[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, mmEnergy[i], mmTime[i], mmPa[i][3]/mmPa[i][1]};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, mmEnergy[i]*scale[pad.first][pad.second], mmTime[i],
                                     mmPa[i][3]/mmPa[i][1]};
            mmCenterMatched_.push_back(mmHitMatched);
          }
          // Aget3
          else if(mmAget[i] == 3) {
            std::pair<Int_t, Int_t> pad = MM_Map_Asad0_Aget3[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, mmEnergy[i], mmTime[i], mmPa[i][3]/mmPa[i][1]};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, mmEnergy[i]*scale[pad.first][pad.second], mmTime[i],
                                     mmPa[i][3]/mmPa[i][1]};
            mmCenterMatched_.push_back(mmHitMatched);
          }
        }
        // Asad1
        else if(mmAsad[i] == 1) {
          // Aget0
          if(mmAget[i] == 0) {
            std::pair<Int_t, Int_t> pad = MM_Map_Asad1_Aget0[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, mmEnergy[i], mmTime[i], mmPa[i][3]/mmPa[i][1]};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, mmEnergy[i]*scale[pad.first][pad.second], mmTime[i],
                                     mmPa[i][3]/mmPa[i][1]};
            mmCenterMatched_.push_back(mmHitMatched);
          }
          // Aget1
          else if(mmAget[i] == 1) {
            std::pair<Int_t, Int_t> pad = MM_Map_Asad1_Aget1[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, mmEnergy[i], mmTime[i], mmPa[i][3]/mmPa[i][1]};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, mmEnergy[i]*scale[pad.first][pad.second], mmTime[i],
                                     mmPa[i][3]/mmPa[i][1]};
            mmCenterMatched_.push_back(mmHitMatched);
          }
          // Aget2
          else if(mmAget[i] == 2) {
            std::pair<Int_t, Int_t> pad = MM_Map_Asad1_Aget2[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, mmEnergy[i], mmTime[i], mmPa[i][3]/mmPa[i][1]};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, mmEnergy[i]*scale[pad.first][pad.second], mmTime[i],
                                     mmPa[i][3]/mmPa[i][1]};
            mmCenterMatched_.push_back(mmHitMatched);
          }
          // Aget3
          else if(mmAget[i] == 3) {
            std::pair<Int_t, Int_t> pad = MM_Map_Asad1_Aget3[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, mmEnergy[i], mmTime[i], mmPa[i][3]/mmPa[i][1]};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, mmEnergy[i]*scale[pad.first][pad.second], mmTime[i],
                                     mmPa[i][3]/mmPa[i][1]};
            mmCenterMatched_.push_back(mmHitMatched);
          }
        }
        // Asad2
        else if(mmAsad[i] == 2) {
          // Aget0 - Strips and Chains (Beam left)
          if(mmAget[i] == 0) {
            Int_t bin = MM_Map_Asad2_Aget0[mmChan[i]];
            if(mmChan[i] < 34) { // Strips
              mmChainStrip mmHit = {bin, mmEnergy[i], mmTime[i]};
              mmLeftStrip_.push_back(mmHit);
            }
            else if(mmChan[i] > 33 && mmChan[i] < 68) { // Chains
              mmChainStrip mmHit = {bin, mmEnergy[i], mmTime[i]};
              mmLeftChain_.push_back(mmHit);
            }
          }
          // Aget1 - Strips and Chains (Beam left)
          else if(mmAget[i] == 1) {
            Int_t bin = MM_Map_Asad2_Aget1[mmChan[i]];
            if(mmChan[i] < 34) { // Strips
              mmChainStrip mmHit = {bin, mmEnergy[i], mmTime[i]};
              mmLeftStrip_.push_back(mmHit);
            }
            else if(mmChan[i] > 33 && mmChan[i] < 68) { // Chains
              mmChainStrip mmHit = {bin, mmEnergy[i], mmTime[i]};
              mmLeftChain_.push_back(mmHit);
            }
          }
          // Aget2 - Outside Central Pads
          else if(mmAget[i] == 2) {
            std::pair<Int_t, Int_t> pad = MM_Map_Asad2_Aget2[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, mmEnergy[i], mmTime[i], mmPa[i][3]/mmPa[i][1]};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, mmEnergy[i]*scale[pad.first][pad.second], mmTime[i],
                                     mmPa[i][3]/mmPa[i][1]};
            mmCenterMatched_.push_back(mmHitMatched);
          }
          // Aget3 - Outside Central Pads
          else if(mmAget[i] == 3) {
            std::pair<Int_t, Int_t> pad = MM_Map_Asad2_Aget3[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, mmEnergy[i], mmTime[i], mmPa[i][3]/mmPa[i][1]};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, mmEnergy[i]*scale[pad.first][pad.second], mmTime[i],
                                     mmPa[i][3]/mmPa[i][1]};
            mmCenterMatched_.push_back(mmHitMatched);
          }
        }
        // Asad3
        else if(mmAsad[i] == 3) {
          // Aget0 - Strips and Chains (Beam right)
          if(mmAget[i] == 0) {
            Int_t bin = MM_Map_Asad3_Aget0[mmChan[i]];
            if(mmChan[i] < 34) { // Strips
              mmChainStrip mmHit = {bin, mmEnergy[i], mmTime[i]};
              mmRightStrip_.push_back(mmHit);
            }
            else if(mmChan[i] > 33 && mmChan[i] < 68) { // Chains
              mmChainStrip mmHit = {bin, mmEnergy[i], mmTime[i]};
              mmRightChain_.push_back(mmHit);
            }
          }
          // Aget1 - Strips and Chains (Beam right)
          else if(mmAget[i] == 1) {
            Int_t bin = MM_Map_Asad3_Aget1[mmChan[i]];
            if(mmChan[i] < 34) { // Strips
              mmChainStrip mmHit = {bin, mmEnergy[i], mmTime[i]};
              mmRightStrip_.push_back(mmHit);
            }
            else if(mmChan[i] > 33 && mmChan[i] < 68) { // Chains
              mmChainStrip mmHit = {bin, mmEnergy[i], mmTime[i]};
              mmRightChain_.push_back(mmHit);
            }
          }
          // Aget2
          else if(mmAget[i] == 2) {
            std::pair<Int_t, Int_t> pad = MM_Map_Asad3_Aget2[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, mmEnergy[i], mmTime[i], mmPa[i][3]/mmPa[i][1]};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, mmEnergy[i]*scale[pad.first][pad.second], mmTime[i],
                                     mmPa[i][3]/mmPa[i][1]};
            mmCenterMatched_.push_back(mmHitMatched);
          }
          // Aget3
          else if(mmAget[i] == 3) {
            std::pair<Int_t, Int_t> pad = MM_Map_Asad3_Aget3[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, mmEnergy[i], mmTime[i], mmPa[i][3]/mmPa[i][1]};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, mmEnergy[i]*scale[pad.first][pad.second], mmTime[i],
                                     mmPa[i][3]/mmPa[i][1]};
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

    //*********//
    // Silicon //
    //*********//

    // Look at all Si that fired and gate
    std::vector<siDetect> siDetectNew_;
    for(auto si : siDetect_) {
      siDetectNew_.push_back(si);
      // if(si.time > 3600 && si.time < 5200) siDetectNew_.push_back(si);
    }

    // If more than 1 Si fired, get highest energy (good enough for now for (p, p))
    siDetect highSi;
    Double_t maxEnergy = 0.;
    for(auto si : siDetectNew_) {
      if(si.energy > maxEnergy) {
        maxEnergy = si.energy;
        highSi = si;
      }
    }

    siDet = highSi.detect;
    siQuad = highSi.quad;
    siChannel = highSi.channel;
    siEnergy = highSi.energy;
    siEnergyCal = 0.;
    if(siDet < 10) {
      siEnergyCal = siEnergy*siEForwardCalibration[siDet][siQuad].first + siEForwardCalibration[siDet][siQuad].second;
    }
    else {
      siEnergyCal = siEnergy*siELeftCalibration[siDet - 10][siQuad].first + siELeftCalibration[siDet - 10][siQuad].second;
    }
    siTime = highSi.time;

    if(siEnergyCal < 1.) continue;

    if(siDet < 10) {
      hSiEForwardDet[siDet]->Fill(siEnergy);
      hSiEForwardDetCal[siDet]->Fill(siEnergyCal);
      hSiTForwardDet[siDet]->Fill(siTime);
    }

    if(siDet > 5) continue;
    if(siDet < 2) continue;

    if(siDet == 5 && (siEnergyCal < 3100 || siEnergyCal > 3600)) continue;
    if(siDet == 4 && (siEnergyCal < 4000 || siEnergyCal > 4400)) continue;
    if(siDet == 3 && (siEnergyCal < 3750 || siEnergyCal > 4150)) continue;
    if(siDet == 2 && (siEnergyCal < 3800 || siEnergyCal > 4100)) continue;

    //************//
    // Micromegas //
    //************//

    // Loop through Micromegas and get total number fired
    std::map<Int_t, Int_t> centralPadTotal;
    for(UInt_t i = 0; i < 128; i++) {
      centralPadTotal[i] = 0;
    }
    for(auto mm : mmCenterMatched_) {
      centralPadTotal[mm.row]++;
      hMicroMegasCenterTime->Fill(mm.row, mm.time - siTime);
    }

    // Make sure only one pad fired
    for(auto mm : mmCenterMatched_) {
      if(centralPadTotal[mm.row] != 1) continue;
      hMicroMegasCenterCumulative->Fill(mm.column - 3, mm.row);
      totalFired[mm.row][mm.column]++;
      totalEnergy[mm.row][mm.column] += mm.energy;
    }

    //  ** End of event by event analysis ** //
    FillTree();

  }

  Double_t averageEnergy[128][6];
  for(Int_t i = 0; i < 128; i++) {
    for(Int_t j = 0; j < 6; j++) {
      averageEnergy[i][j] = 0.;
      if(totalFired[i][j] == 0) continue;
      averageEnergy[i][j] = totalEnergy[i][j]/static_cast<Double_t>(totalFired[i][j]);
      for(UInt_t k = 0; k < static_cast<UInt_t>(averageEnergy[i][j]); k++) {
        hMicroMegasCenterEnergyAverage->Fill(j - 3, i);
      }
    }
  }

//  // Gain match all pads
//  FILE* gainFile = fopen("centerGain.dat", "w");
//  Double_t padBest = averageEnergy[64][2];
//  Double_t rowEnergySlope = 0.01338555;
//  Double_t rowEnergyIntercept = 7.7997787;
//  Double_t padBestEnergy = 8.85;
//  Double_t gainScale[128][6];
//  for(Int_t i = 0; i < 128; i++) {
//    for(Int_t j = 0; j < 6; j++) {
//      Double_t energyRow = i*rowEnergySlope + rowEnergyIntercept;
//      gainScale[i][j] = padBest*energyRow/(averageEnergy[i][j]*padBestEnergy);
//      fprintf(gainFile, "%d %d %f\n", i, j, gainScale[i][j]);
//    }
//  }
//  fflush(gainFile);
//  fclose(gainFile);

  hMicroMegasCenterCumulative->Write();
  hMicroMegasCenterTime->Write();
  hMicroMegasCenterEnergyAverage->Write();

  // Forward Si Detectors
  for(UInt_t i = 0; i < 10; i++) {
//    hSiEForwardDet[i]->Write();
    hSiTForwardDet[i]->Write();
    hSiEForwardDetCal[i]->Write();
  }

  hCWTECentral->Write();

  WriteTree();

  file->Close();
}