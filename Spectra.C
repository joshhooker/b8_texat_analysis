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

  /////////////////
  // Set up cuts //
  /////////////////
  TFile *cutFile = TFile::Open("cuts.root");

  // dE cuts
  for(UInt_t i = 0; i < 10; i++) {
    if(i == 8) continue;
    dEEForwardCut[i] = (TCutG*)cutFile->Get(Form("dEEForward_d%dCut", i));
  }

  // Angle vs Energy Cut
  angleTotEnergyCut[0] = (TCutG*)cutFile->Get("angleTotEnergy_d0Cut");
  angleTotEnergyCut[1] = (TCutG*)cutFile->Get("angleTotEnergy_d1Cut");
  angleTotEnergyCut[9] = (TCutG*)cutFile->Get("angleTotEnergy_d9Cut");

  // Chain/Strip vs Time Cuts
  for(UInt_t i = 0; i < 10; i++) {
    if(i == 4 || i == 5) continue;
    for(UInt_t j = 0; j < 4; j++) {
      if(i == 7 && j == 0) continue;
      timeChainForwardCut[i][j] = (TCutG*)cutFile->Get(Form("timeChainForward_d%d_q%dCut", i, j));
      timeStripForwardCut[i][j] = (TCutG*)cutFile->Get(Form("timeStripForward_d%d_q%dCut", i, j));
    }
  }

  cutFile->Close();

  InitChannelMap();
  InitHistograms();
  InitCanvas();
  InitVariables();

  InitSiEForwardCalibration();
  InitCsIECalibration();

  InitCentralPadGainMatch();
  InitAverageBeamEnergy();

  InitTree();

  Long64_t nentries = fChain->GetEntriesFast();

  printf("Starting Main Loop\n");

  Long64_t nbytes = 0, nb = 0;
//  for(Long64_t jentry = 0; jentry < 50; jentry++) {
//  for(Long64_t jentry = 49983; jentry < 49984; jentry++) {
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

    hIonizationChamberE->Fill(icE);
    hIonizationChamberT->Fill(icT);

    if(icE < 1400 || icE > 1900) continue;
    if(icT < 5500 || icT > 6500) continue;

    //*********//
    // Silicon //
    //*********//

    // Look at all Si that fired and gate
    std::vector<siDetect> siDetectNew_;
    for(auto si : siDetect_) {
      if(si.time > 3600 && si.time < 5200) siDetectNew_.push_back(si);
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
      hSiEForward[siDet][siQuad]->Fill(siEnergy);
      hSiEForwardCal[siDet][siQuad]->Fill(siEnergyCal);
    }

    Bool_t left = false;
    Bool_t central = false;
    Bool_t right = false;
    if(siDet == 4 || siDet == 5) { // Central detectors
      central = true;
    }
    else if(siDet < 4) { // Right detectors
      right = true;
    }
    else if(siDet > 5) { // Left detectors
      left = true;
    }

    Bool_t sideDet = false;
    if(siDet > 9) sideDet = true;

    if(sideDet) continue;

    // Skip bad Kiev detectors
//    if(siDet == 2 || siDet == 3 || siDet == 6 || siDet == 7) continue;

    //*****//
    // CsI //
    //*****//

    // Find if CsI behind Si fired
    punchthrough = false;
    csiEnergy = 0.;
    csiEnergyCal = 0.;
    csiTime = 0.;
    for(auto csi : csiDetect_) {
      if(csi.detect == siDet) {
        csiEnergy = csi.energy;
        csiTime = csi.time;
        if(csiTime < 11000 || csiTime > 12250) continue;
        punchthrough = true;
        if(!sideDet) {
          hCsIETForward[siDet]->Fill(csiEnergy, csiTime);
          csiEnergyCal = csiEnergy*csiEForwardCalibration[siDet].first + csiEForwardCalibration[siDet].second;
          hSiCsIEForwardDet[siDet]->Fill(siEnergy, csiEnergy);
          hSiCsIEForwardDetCal[siDet]->Fill(siEnergyCal, csiEnergy);
          hSiCsIEForward[siDet][siQuad]->Fill(siEnergy, csiEnergy);
          hSiCsIEForwardCal[siDet][siQuad]->Fill(siEnergyCal, csiEnergy);
        }
      }
    }

    if(punchthrough) {
      hSumSiEForwardDet[siDet]->Fill(siEnergy + csiEnergy, siEnergy);
      hSumSiEForward[siDet][siQuad]->Fill(siEnergy + csiEnergy, siEnergy);

      hSumCsIEForwardDet[siDet]->Fill(siEnergy + csiEnergy, csiEnergy);
      hSumCsIEForward[siDet][siQuad]->Fill(siEnergy + csiEnergy, csiEnergy);
    }

    totalEnergy = siEnergyCal + csiEnergyCal;

    //************//
    // Micromegas //
    //************//

    // Gate on time for central pads
    std::vector<mmCenter> mmCenterMatchedReduced_;
    for(auto mm : mmCenterMatched_) {
      Double_t time = mm.time - siTime;
      if(mm.row < 2) continue;
      if(mm.row < 112) {
        if(time < 900 || time > 1600) continue;
        mmCenterMatchedReduced_.push_back(mm);
      }
      else {
        if(time > 2500) continue;
        mmCenterMatchedReduced_.push_back(mm);
      }
    }

    // Sorting mmCenterMatchedReduced_ by row and reduce noise
    std::vector<mmCenter> mmCenterMatchedReducedNoise_;
    if(!mmCenterMatchedReduced_.empty()) {
      std::sort(mmCenterMatchedReduced_.begin(), mmCenterMatchedReduced_.end(), sortByRowMMCenter());
      mmCenterMatchedReducedNoise_ = CenterReduceNoise(mmCenterMatchedReduced_);
    }

    // Reduce mmCenter to one entry per row
    std::map<Int_t, Double_t> centralPadPosition;
    std::map<Int_t, Double_t> centralPadTotalEnergy;
    std::map<Int_t, Double_t> centralPadTime;
    std::map<Int_t, Int_t> centralPadColumn;
    std::map<Int_t, Int_t> centralPadTotal;
    for(UInt_t i = 0; i < 128; i++) {
      centralPadTotal[i] = 0;
    }
    for(auto mm : mmCenterMatchedReducedNoise_) {
      Double_t xPosition = mm.column*3.5 - 8.75;
      Double_t yPosition = mm.row*rowConversion + rowConversionOffset;
      hMicroMegasCenterCumulative->Fill(mm.column - 3, mm.row);
      hMicroMegasCenterCumulativePositionRaw->Fill(xPosition, yPosition);
      hMicroMegasCenterTime->Fill(mm.row, mm.time - siTime);
      hMicroMegasCenterHeight->Fill(mm.row, heightOffset - (mm.time - siTime)*driftVelocity);
      if(centralPadTotal[mm.row] == 0) {
        if(mm.energy < 0) continue;
        centralPadPosition[mm.row] = (mm.column*3.5 - 8.75)*mm.energy;
        centralPadTotalEnergy[mm.row] = mm.energy;
        centralPadTime[mm.row] = mm.time;
        centralPadColumn[mm.row] = mm.column;
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
    std::map<Int_t, Double_t>::iterator it;
    for(it = centralPadPosition.begin(); it != centralPadPosition.end(); it++) {
      Int_t row = it->first;
      mmTrack hit = {0, 0., 0., 0., 0., 0., 0};
      if(centralPadTotalEnergy[row] == 0) continue;
      if(centralPadTotal[row] == 0) continue;
      hit.row = row;
      hit.xPosition = centralPadPosition[row]/centralPadTotalEnergy[row];
      hit.yPosition = row*rowConversion + rowConversionOffset;
      hit.time = centralPadTime[row]/static_cast<Double_t>(centralPadTotal[row]) - siTime;
      hit.energy = centralPadTotalEnergy[row];
      hit.height = heightOffset - hit.time*driftVelocity;
      hit.total = centralPadTotal[row];
//      hMicroMegasCenterCumulative->Fill(centralPadColumn[row] - 3, row);
      hMicroMegasCenterCumulativePosition->Fill(hit.xPosition, hit.yPosition);
      if(row < 112) {
        mmCenterBeamTotal_.push_back(hit);
      }
      else if(row > 111) {
        mmCenterProton_.push_back(hit);
      }
    }

    Bool_t event = false;
    if(central) {
      event = AnalysisForwardCentral(mmCenterMatched_, mmCenterBeamTotal_, mmCenterProton_, centralPadTotalEnergy);
    }
    if((left || right) && !sideDet) {
      event = AnalysisForwardSide(mmCenterMatched_, mmCenterBeamTotal_, mmLeftChain_, mmLeftStrip_, mmRightChain_,
          mmRightStrip_);
    }

    if(!event) continue;

    //  ** End of event by event analysis ** //
    FillTree();

  }

  // TGraph boundaries
  auto *h_track_bound = new TGraph();
  h_track_bound->SetMarkerStyle(8);
  h_track_bound->SetMarkerColor(5);
  h_track_bound->SetPoint(0, -150, -300);
  h_track_bound->SetPoint(1, 150, 300);
  h_track_bound->SetName("Bound");
  h_track_bound->Write();
  delete h_track_bound;

//  hIonizationChamberE->Write();
//  hIonizationChamberT->Write();

  hMicroMegasCenterCumulative->Write();
  hMicroMegasCenterCumulativePosition->Write();
  hMicroMegasCenterCumulativePositionRaw->Write();
  hMicroMegasCenterTime->Write();
  hMicroMegasCenterHeight->Write();

  // Forward Si Detectors
//  for(UInt_t i = 0; i < 10; i++) {
//    hSiEForwardDet[i]->Write();
//    hSiTForwardDet[i]->Write();
//    hSiEForwardDetCal[i]->Write();
//    for (int j = 0; j < 4; j++) {
//      hSiEForward[i][j]->Write();
//      hSiEForwardCal[i][j]->Write();
//      hSiTForward[i][j]->Write();
//    }
//    hCsIEForward[i]->Write();
//    hCsITForward[i]->Write();
//  }

  // Left Si Detectors
//  for(UInt_t i = 0; i < 6; i++) {
//    hSiELeftDet[i]->Write();
//    hSiTLeftDet[i]->Write();
//    hSiELeftDetCal[i]->Write();
//    for(UInt_t j = 0; j < 4; j++) {
//      hSiELeft[i][j]->Write();
//      hSiTLeft[i][j]->Write();
//    }
//    hCsIELeft[i]->Write();
//    hCsITLeft[i]->Write();
//    hCsIETLeft[i]->Write();
//  }

//   hMicroMegasCenterCumulative->Write();
//   hMicroMegasCenterEnergyCumulative->Write();
//   hMicroMegasCenterEnergyAverage->Write();
//   hMicroMegasCenterEnergyAverageScaled->Write();
//   hMicroMegasCenterTimeAverage->Write();

  // Forward CsI Energy vs Time
//  for(UInt_t i = 0; i < 10; i++) {
//    hCsIETForward[i]->Write();
//  }

//   Forward Si Energy vs CsI Energy
//   for(UInt_t i = 0; i < 10; i++) {
//     hSiCsIEForwardDet[i]->Write();
//     hSiCsIEForwardDetCal[i]->Write();
//      for(UInt_t j = 0; j < 4; j++) {
//        hSiCsIEForward[i][j]->Write();
//        hSiCsIEForwardCal[i][j]->Write();
//      }
//   }

  // Forward Si + CsI vs Si (raw)
  // for(UInt_t i = 0; i < 10; i++) {
  //   hSumSiEForwardDet[i]->Write();
  //   for(UInt_t j = 0; j < 4; j++) {
  //     hSumSiEForward[i][j]->Write();
  //   }
  // }

  // Forward Si + CsI vs CsI (raw)
  // for(UInt_t i = 0; i < 10; i++) {
  //   hSumCsIEForwardDet[i]->Write();
  //   for(UInt_t j = 0; j < 4; j++) {
  //     hSumCsIEForward[i][j]->Write();
  //   }
  // }

  // Forward dE vs Si Energy
//  for(UInt_t i = 0; i < 10; i++) {
//    hdEEForward[i]->Write();
//    hdEEForwardCal[i]->Write();
//    hdEEForwardCalTotal[i]->Write();
//  }

  // Forward Hough Angle
//  for(UInt_t i = 0; i < 10; i++) {
//    hHoughAngle[i]->Write();
//  }

  // Forward Vertex vs Si Energy
//  for(UInt_t i = 0; i < 10; i++) {
//    hVertexSiEForward[i]->Write();
//    hVertexSiEForwardCal[i]->Write();
//    hVertexSiEForwardCalTotal[i]->Write();
//  }
//  hVertexSiETotalRegion3->Write();
//  hVertexCMERegion3->Write();

  // Forward Angle vs Si Energy
//  for(UInt_t i = 0; i < 10; i++) {
//    hAngleEForward[i]->Write();
//    hAngleEForwardCal[i]->Write();
//    hAngleEForwardCalTotal[i]->Write();
//    hAngleEForwardProtonEnergy[i]->Write();
//    hAngleEForwardCMEnergy[i]->Write();
//  }

  // Forward Vertex vs Angle
//  for(UInt_t i = 0; i < 10; i++) {
//    hVertexAngleForward[i]->Write();
//  }

  // Time vs Column/Strip Number Forward Detectors
//  for(UInt_t i = 0; i < 10; i++) {
//    for(UInt_t j = 0; j < 4; j++) {
//      hTimeChainForward[i][j]->Write();
//      hTimeStripForward[i][j]->Write();
//    }
//  }

  // Time vs Central Row Forward Detectors
//  for(UInt_t i = 0; i < 10; i++) {
//    hTimeCentralForward[i]->Write();
//  }

  // Forward Wall XZ Hit Positions
//  hHitPositionsXZForward->Write();
//  for(UInt_t i = 0; i < 10; i++) {
//    hHitPositionsXZForwardInd[i]->Write();
//  }

  hCWTECentral->Write();

  DivideTargetThickness(s1);
  ReadSolidAngle();
  SolidAngle(s1);
  s1->Scale(1./numberB8);
  s1->Write();
  WriteSpectrumToFile(s1, 3);

  for(Int_t i = 0; i < centerEnergyCanvasNum; i++) {
    centerEnergyCanvas[i]->Write();
  }

  WriteTree();

  file->Close();
}

Bool_t Spectra::AnalysisForwardCentral(std::vector<mmCenter> centerMatched_, std::vector<mmTrack> centerBeamTotal_,
                                     std::vector<mmTrack> centerProton_, std::map<Int_t, Double_t> centralPadTotalEnergy) {

  // Find dE for center region
  dE = 0.;
  Int_t totalRows = 0;
  for(Int_t i = 116; i < 124; i++) {
    if(i == 117) continue;
    auto iterator = centralPadTotalEnergy.find(i);
    if(iterator != centralPadTotalEnergy.end()) {
      dE += centralPadTotalEnergy[i];
      totalRows++;
    }
    dE /= static_cast<Double_t>(totalRows);
  }

  hdEEForward[siDet]->Fill(siEnergy, dE);
  hdEEForwardCal[siDet]->Fill(siEnergyCal, dE);
  hdEEForwardCalTotal[siDet]->Fill(totalEnergy, dE);

  if(!dEEForwardCut[siDet]->IsInside(siEnergy, dE)) return false;

  for(auto mm : centerMatched_) {
    hTimeCentralForward[siDet]->Fill(mm.row, mm.time - siTime);
  }

  if(!centerMatched_.empty()) {
    std::sort(centerMatched_.begin(), centerMatched_.end(), sortByRowMMCenter());
  }

  return true;
}

Bool_t Spectra::AnalysisForwardSide(std::vector<mmCenter> centerMatched_, std::vector<mmTrack> centerBeamTotal_,
                                  std::vector<mmChainStrip> leftChain_, std::vector<mmChainStrip> leftStrip_,
                                  std::vector<mmChainStrip> rightChain_, std::vector<mmChainStrip> rightStrip_) {
  // Sorting mmCenterBeamTotal_ by row
  if(!centerBeamTotal_.empty()) {
    std::sort(centerBeamTotal_.begin(), centerBeamTotal_.end(), sortByRowMMTrack());
  }

  // Sorting mmLeftChain_
  if(!leftChain_.empty()) {
    std::sort(leftChain_.begin(), leftChain_.end(), sortByRowMMChainStrip());
  }

  // Sorting mmLeftStrip_
  if(!leftStrip_.empty()) {
    std::sort(leftStrip_.begin(), leftStrip_.end(), sortByRowMMChainStrip());
  }

  // Sorting mmRightChain_
  if(!rightChain_.empty()) {
    std::sort(rightChain_.begin(), rightChain_.end(), sortByRowMMChainStrip());
  }

  // Sorting mmRightStrip_
  if(!rightStrip_.empty()) {
    std::sort(rightStrip_.begin(), rightStrip_.end(), sortByRowMMChainStrip());
  }

  Bool_t left = false;
  Bool_t right = false;
  if(siDet < 4) right = true;
  else if(siDet > 5) left = true;

  // Look at time vs chain/strip
  if(left) {
    // Chain
    for(auto mm : leftChain_) {
      hTimeChainForward[siDet][siQuad]->Fill(mm.row, mm.time - siTime);
    }
    // Strip
    for(auto mm : leftStrip_) {
      hTimeStripForward[siDet][siQuad]->Fill(mm.row, mm.time - siTime);
    }
  }
  else if(right) {
    // Chain
    for(auto mm : rightChain_) {
      hTimeChainForward[siDet][siQuad]->Fill(mm.row, mm.time - siTime);
    }
    for(auto mm : rightStrip_) {
      hTimeStripForward[siDet][siQuad]->Fill(mm.row, mm.time - siTime);
    }
  }

  // Gate on time vs chain/strip
  std::vector<mmChainStrip> leftChainReduced_;
  std::vector<mmChainStrip> leftStripReduced_;
  std::vector<mmChainStrip> rightChainReduced_;
  std::vector<mmChainStrip> rightStripReduced_;
  if(left) {
    // Chain
    for(auto mm : leftChain_) {
      if(timeChainForwardCut[siDet][siQuad]->IsInside(mm.row, mm.time - siTime)) {
        leftChainReduced_.push_back(mm);
      }
    }
    // Strip
    for(auto mm : leftStrip_) {
      if(timeStripForwardCut[siDet][siQuad]->IsInside(mm.row, mm.time - siTime)) {
        leftStripReduced_.push_back(mm);
      }
    }
  }
  else if(right) {
    // Chain
    for(auto mm : rightChain_) {
      if(timeChainForwardCut[siDet][siQuad]->IsInside(mm.row, mm.time - siTime)) {
        rightChainReduced_.push_back(mm);
      }
    }
    // Strip
    for(auto mm : rightStrip_) {
      if(timeStripForwardCut[siDet][siQuad]->IsInside(mm.row, mm.time - siTime)) {
        rightStripReduced_.push_back(mm);
      }
    }
  }

  // Find dE for left and right regions
  dE = 0.;
  if(left) {
    for(auto mm: leftStripReduced_) {
      dE += mm.energy;
    }
    dE /= static_cast<Double_t>(leftStripReduced_.size());
  }
  else if(right) {
    for(auto mm : rightStripReduced_) {
      dE += mm.energy;
    }
    dE /= static_cast<Double_t>(rightStripReduced_.size());
  }

  hdEEForward[siDet]->Fill(siEnergy, dE);
  hdEEForwardCal[siDet]->Fill(siEnergyCal, dE);
  hdEEForwardCalTotal[siDet]->Fill(totalEnergy, dE);

  if(!dEEForwardCut[siDet]->IsInside(totalEnergy, dE)) return false;

  for(auto mm : centerMatched_) {
    hTimeCentralForward[siDet]->Fill(mm.row, mm.time - siTime);
  }

  // Match strips with chains
  std::vector<mmTrack> chainStripMatchedLeft;
  std::vector<mmTrack> chainStripMatchedRight;
  std::vector<mmTrack> chainStripRawLeft;
  std::vector<mmTrack> chainStripRawRight;
  if(!leftChainReduced_.empty() && !leftStripReduced_.empty()) {
    ChainStripMatch(chainStripMatchedLeft, chainStripRawLeft, leftChainReduced_, leftStripReduced_, true, siTime);
  }
  if(!rightChainReduced_.empty() && !rightStripReduced_.empty()) {
    ChainStripMatch(chainStripMatchedRight, chainStripRawRight, rightChainReduced_, rightStripReduced_, false,
                    siTime);
  }

  // Assign proton track for side regions
  std::vector<mmTrack> protonTrack;
  if(left) {
    protonTrack = chainStripMatchedLeft;
  }
  else if(right) {
    protonTrack = chainStripMatchedRight;
  }

  // Find vertex using the proton track in side regions
  vertexPositionX = 0.;
  vertexPositionY = -400.;
  vertexPositionZ = 0.;

  Double_t protonEnergy = 0.;

  if(protonTrack.empty()) return false;

  auto *fitProton = new HoughTrack();
  fitProton->AddTrack(protonTrack, siDet, siQuad);
  Double_t minDist = fitProton->FitRestricted();
  std::vector<Double_t> parsProton = fitProton->GetPars();
  Double_t houghAngleXY = fitProton->GetHoughAngleXY();
  Double_t houghAngleYZ = fitProton->GetHoughAngleYZ();
  hHoughAngle[siDet]->Fill(houghAngleXY);

  delete fitProton;

  ULong_t mmCenterSize = centerBeamTotal_.size();

  if(mmCenterSize == 0) {
    Double_t beamX_old = 0.;
    Double_t beamY_old = 250.;
    Double_t beamZ_old = 0.;
    Double_t beamX = beamX_old;
    Double_t beamY = beamY_old;
    Double_t beamZ = beamZ_old;
    Double_t x, y, z;
    line(beamY, parsProton, x, y, z);
    Double_t protonX_old = x;
    Double_t protonX = protonX_old;
    Double_t xDiff_old = beamX - protonX_old;
    Double_t xDiff = xDiff_old;
    Double_t yPos = beamY_old;
    while(yPos > -300.) {
      beamY = yPos;
      line(beamY, parsProton, x, y, z);
      protonX = x;
      xDiff = beamX - protonX;
      if(xDiff_old*xDiff < 0) {
        Double_t m = (beamY - beamY_old)/(xDiff - xDiff_old);
        Double_t b = beamY - m*xDiff;
        vertexPositionX = beamX;
        vertexPositionY = beamY;
        vertexPositionZ = beamZ;
        break;
      }
      beamY_old = beamY;
      protonX_old = protonX;
      xDiff_old = xDiff;
      yPos -= 1.;
    }
    // Did not fit well
    if(vertexPositionY == -400) {
      Double_t m_xcomponent = fabs(parsProton[1]);
      Double_t xAngle = atan(m_xcomponent);

      Double_t yDist = siXPosForward[siDet][siQuad]/fabs(tan(xAngle));

      if(250. - yDist < -400) {
        vertexPositionY = -400.;
      }
      else {
        vertexPositionY = 250. - yDist;
      }
    }
  }
  else {
    // Loop through beam in central region starting from last
    Double_t beamX_old = centerBeamTotal_[mmCenterSize - 1].xPosition;
    Double_t beamY_old = centerBeamTotal_[mmCenterSize - 1].yPosition;
    Double_t beamZ_old = centerBeamTotal_[mmCenterSize - 1].height;
    Double_t beamX = beamX_old;
    Double_t beamY = beamY_old;
    Double_t beamZ = beamZ_old;
    Double_t x, y, z;
    line(beamY, parsProton, x, y, z);
    Double_t protonX_old = x;
    Double_t protonX = protonX_old;
    Double_t xDiff_old = beamX - protonX_old;
    Double_t xDiff = xDiff_old;
    Bool_t foundVertex = false;
    for(ULong_t i = mmCenterSize - 1; i > -1; i--) {
      beamX = centerBeamTotal_[i].xPosition;
      beamY = centerBeamTotal_[i].yPosition;
      beamZ = centerBeamTotal_[i].height;
      line(beamY, parsProton, x, y, z);
      protonX = x;
      xDiff = beamX - protonX;
      if(xDiff_old*xDiff < 0) {
        Double_t m = (beamY - beamY_old)/(xDiff - xDiff_old);
        Double_t b = beamY - m*xDiff;
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
      Double_t beamY = -1.;
      while(beamY > -300.) {
        line(beamY, parsProton, x, y, z);
        protonX = x;
        xDiff = beamX - x;
        if(xDiff*xDiff_old < 0) {
          Double_t m = (beamY - beamY_old)/(xDiff - xDiff_old);
          Double_t b = beamY - m*xDiff;
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
    // Did not fit well
    if(vertexPositionY == -400) {
      Double_t m_xcomponent = fabs(parsProton[1]);
      Double_t xAngle = atan(m_xcomponent);

      Double_t yDist = siXPosForward[siDet][siQuad]/fabs(tan(xAngle));

      if(250. - yDist < -400) {
        vertexPositionY = -400.;
      }
      else {
        vertexPositionY = 250. - yDist;
      }
    }
  }

  // Plot XZ hit position of forward non-central detectors
  Double_t x, y, z;
  line(siYPosForward, parsProton, x, y, z);
  siPosX = x;
  siPosY = siYPosForward;
  siPosZ = z;
  hHitPositionsXZForward->Fill(x, z);
  hHitPositionsXZForwardInd[siDet]->Fill(x, z);

  Double_t vertexToSi = siYPosForward - vertexPositionY;

  Double_t siX, siY, siZ;
  line(siYPosForward, parsProton, siX, siY, siZ);

  // Double_t angleX = atan(fabs(siX)/(vertexToSi));
  // Double_t angleZ = atan(fabs(siZ)/(vertexToSi));
  // Double_t cosAngle = cos(angleX)*cos(angleZ);
  // angle = acos(cosAngle);

  TVector3 v1(siX, vertexToSi, siZ);
  TVector3 v2(0, vertexToSi, 0);
  angle = v1.Angle(v2);
  Double_t cosAngle = cos(angle);

  Double_t pathLength = vertexToSi/cosAngle;

  if(pathLength < 0) pathLength = 200.;
  if(pathLength > 700) pathLength = 700.;

  protonEnergy = protonMethane->AddBack(totalEnergy/1000., pathLength);
  Double_t beamEnergy = protonEnergy*(m1 + m2)*(m1 + m2)/(4.*m1*m2*cosAngle*cosAngle);

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
    if(angleTotEnergyCut[siDet]->IsInside(protonEnergy, angle)) {
      hVertexSiETotalRegion3->Fill(totalEnergy, vertexPositionY);
      hVertexCMERegion3->Fill(cmEnergy, vertexPositionY);
      s1->Fill(cmEnergy);
    }
  }

  return true;
}

std::vector<mmCenter> Spectra::CenterReduceNoise(std::vector<mmCenter> center) {
  std::vector<mmCenter> mmNew;

  std::vector<Int_t> rows[128];

  for(auto mm : center) {
    if(mm.row > 111) {
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
          Bool_t found = false;
          for(auto column : rows[mm.row - 3]) {
            if(column == mm.column - 1 || column == mm.column || column == mm.column + 1) found = true;
          }
          if(found) {
            Bool_t inRow = false;
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
          Bool_t found = false;
          for(auto column : rows[mm.row - 4]) {
            if(column == mm.column - 1 || column == mm.column || column == mm.column + 1) found = true;
          }
          if(found) {
            Bool_t inRow = false;
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

void Spectra::ChainStripMatch(std::vector<mmTrack> &chainStripMatched, std::vector<mmTrack> &chainStripRaw,
                              std::vector<mmChainStrip> chain_,
                              std::vector<mmChainStrip> strip_, Bool_t leftSide, Double_t siTime) {
  chainStripMatched.clear();
  chainStripRaw.clear();
  std::vector<mmTrack> totalTime0;
  ChainStripMatchingTime(totalTime0, chain_, strip_, leftSide, siTime, 0);
  ChainStripMatchingTime(chainStripRaw, chain_, strip_, leftSide, siTime, 0);

  size_t numTimeBuckets = ChainStripTime0TimeBuckets(totalTime0);
  // std::cout << "Time buckets: " << numTimeBuckets << std::endl;
  if(numTimeBuckets < 4) {
    ChainStripMatchingBoxTime0(chainStripMatched, totalTime0);
    // std::cout << "BoxTime0: " << chainStripMatched.size() << std::endl;
  }
  // else if(chain_.size() > 12 && strip_.size() > 12) {
  //   ChainStripMatchingTimeSlopeHough(chainStripMatched, chain_, strip_, leftSide, siTime, 10.);
  //   if(chainStripMatched.size() < 2) {
  //     chainStripMatched.clear();
  //     ChainStripMatchingBoxTime0(chainStripMatched, totalTime0);
  //   }
  //   // std::cout << "SlopeHough0: " << chainStripMatched.size() << std::endl;
  // }
  else {
    ChainStripMatchingTime(chainStripMatched, chain_, strip_, leftSide, siTime, 0);
    // std::cout << "MatchingTime: " << chainStripMatched.size() << std::endl;
  }
  // std::cout << chain_.size() << '\t' << strip_.size() << '\t' << chainStripMatched.size() << std::endl;
}

size_t Spectra::ChainStripTime0TimeBuckets(std::vector<mmTrack> matched) {
  // Function that finds the number of different time buckets for the time0 algorithm
  // This is used if you want to change the strip/chain matching algorithm (if all on the same plane)

  std::map<Int_t, Int_t> timeMap;

  for(auto mm : matched) {
    timeMap[mm.time] = 1;
  }

  return timeMap.size();
}

size_t Spectra::ChainStripNumberTimeBuckets(std::vector<mmChainStrip> chain, std::vector<mmChainStrip> strip) {
  // Function that finds the number of different time buckets for strips and chains
  // This is used if you want to change the strip/chain matching algorithm (if all on the same plane)

  std::map<Int_t, Int_t> individualTimeStripMap;
  std::map<Int_t, Int_t> individualTimeChainMap;

  for(auto mm : chain) {
    individualTimeChainMap[mm.time] = 1;
  }
  for(auto mm : strip) {
    individualTimeStripMap[mm.time] = 1;
  }

  size_t sizeTimeChainMap = individualTimeChainMap.size();
  size_t sizeTimeStripMap = individualTimeStripMap.size();

  return std::min(sizeTimeChainMap, sizeTimeStripMap);
}

void Spectra::ChainStripMatchingOutward(std::vector<mmTrack> &chainStripMatched, std::vector<mmChainStrip> chain,
                                        std::vector<mmChainStrip> strip, Bool_t leftSide, Double_t siTime) {
  // Function to match strips and chains outward. Meaning that the matched strips and chains must be moving away from the central region
  // Strip and chains need to be sorted for this
  // Position transform = 10.5 + 1.75/2 + 1.75*chain (chain starting at 0)
  // Row transform = (strip - 1)*2
  // Row position transform = Row * rowConversion + rowConversionOffset

  std::vector<mmChainStrip> stripCopy = strip;
  for(UInt_t i = 0; i < chain.size(); i++) {
    Double_t timeChain = chain[i].time;
    // Look in the strip vector for this time
    auto it = std::find_if(stripCopy.begin(), stripCopy.end(),
      [timeChain] (const mmChainStrip& d) { return d.time == timeChain; });
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
      chainStripMatched.push_back(hit);
      stripCopy.erase(stripCopy.begin(), stripCopy.begin() + std::distance(stripCopy.begin(), it) + 1);
    }
    else continue;
  }
}

void Spectra::ChainStripMatchingBox(std::vector<mmTrack> &chainStripMatched, std::vector<mmChainStrip> chain,
                                    std::vector<mmChainStrip> strip, Bool_t leftSide, Double_t siTime) {
  // Function to match strips and chains. This is the simplest algorithm, making two points in the side
  // region which is where the particle entered and exited.
  // Strips and chains need to be sorted for this
  // Position transform = 10.5 + 1.75/2 + 1.75*chain (chain starting at 0)
  // Row transform = (strip - 1)*2
  // Row position transform = Row * rowConversion + rowConversionOffset
  Size_t chainSize = chain.size();
  Size_t stripSize = strip.size();

  Double_t position0 = 10.5 + 1.75/2 + 1.75*chain[0].row;
  if(leftSide) position0 = -1.*position0;
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
  if(leftSide) position1 = -1.*position1;
  Int_t row1 = strip[stripSize - 1].row*2;
  mmTrack hit1 = {0, 0., 0., 0., 0., 0., 0};
  hit1.row = row1;
  hit1.xPosition = position1;
  hit1.yPosition = hit1.row*rowConversion + rowConversionOffset;
  hit1.time = (chain[chainSize - 1].time + strip[stripSize - 1].time)/2. - siTime;
  hit1.energy = strip[stripSize - 1].energy;
  hit1.height = heightOffset - hit1.time*driftVelocity;
  hit1.total = 1;

  chainStripMatched.push_back(hit0);
  chainStripMatched.push_back(hit1);
}

void Spectra::ChainStripMatchingBoxTime0(std::vector<mmTrack> &chainStripMatched, std::vector<mmTrack> time0) {
  UInt_t size = time0.size();
  if(size == 0) return;
  chainStripMatched.push_back(time0[0]);
  chainStripMatched.push_back(time0[size - 1]);
}

void Spectra::ChainStripMatchingTime(std::vector<mmTrack> &chainStripMatched, std::vector<mmChainStrip> chain,
                                     std::vector<mmChainStrip> strip, Bool_t leftSide, Double_t siTime,
                                     Int_t timeWindow) {
  // Function to match strips and chains. Chains and strips are matched if their time is within the timeWindow parameter
  // The timeWindow parameter is the window of time buckets to look for using the timeResolution parameter
  // which is the ns -> bucket
  // Position transform = 10.5 + 1.75/2 + 1.75*chain (chain starting at 0)
  // Row transform = strip*2
  // Row position transform = Row * rowConversion + rowConversionOffset

  for(auto chainVec : chain) {
    Double_t timeChain = chainVec.time;
    // Loop over strips and find strips when the time is within the timeWindow
    for(auto stripVec : strip) {
      if((stripVec.time > timeChain - timeWindow*timeResolution - timeResolution/2.) &&
          (stripVec.time < timeChain + timeWindow*timeResolution + timeResolution/2.)) {
        mmTrack hit = {0, 0., 0., 0., 0., 0., 0};
        Double_t position = 10.5 + 1.75/2 + 1.75*chainVec.row;
        if(leftSide) position = -position;
        Int_t row = stripVec.row*2;
        hit.row = row;
        hit.xPosition = position;
        hit.yPosition = hit.row*rowConversion + rowConversionOffset;
        hit.time = chainVec.time - siTime;
        hit.energy = stripVec.energy;
        hit.height = heightOffset - hit.time*driftVelocity;
        hit.total = 1;
        chainStripMatched.push_back(hit);
      }
    }
  }
}

void Spectra::ChainStripMatchingTimeSlopeFit(std::vector<mmTrack> &chainStripMatched, std::vector<mmChainStrip> chain,
                                             std::vector<mmChainStrip> strip, Bool_t leftSide, Double_t siTime,
                                             Double_t timeWindow) {
  // Function to match strips and chains by fitting with time individually
  // Then is able to correlate strips and chains by time. Requires multiple time steps in both strips and chains
  // After the fit, it goes through each chain and finds strips based off of the fit when match within the timeWindow input parameter
  // Position transform = 10.5 + 1.75/2 + 1.75*chain (chain starting at 0)
  // Row transform = (strip - 1)*2
  // Row position transform = Row * rowConversion + rowConversionOffset

  auto *hChain = new TGraph();
  auto *hStrip = new TGraph();

  // Fill hChain
  Int_t i = 0;
  for(auto mm : chain) {
    hChain->SetPoint(i, mm.row, mm.time);
    i++;
  }

  // Fill hStrip
  i = 0;
  for(auto mm : strip) {
    hStrip->SetPoint(i, mm.row, mm.time);
    i++;
  }

  // Fit each with a 1-D polynomial
  hChain->Fit("pol1", "Q");
  hStrip->Fit("pol1", "Q");
  TF1 *fitChain = hChain->GetFunction("pol1");
  TF1 *fitStrip = hStrip->GetFunction("pol1");

  hChain->SetName(Form("hChain_%lld", entry));
  hStrip->SetName(Form("hStrip_%lld", entry));
//  hChain->Write();
//  hStrip->Write();

  Double_t chainP0 = fitChain->GetParameter(0);
  Double_t chainP1 = fitChain->GetParameter(1);
  Double_t stripP0 = fitStrip->GetParameter(0);
  Double_t stripP1 = fitStrip->GetParameter(1);

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
        chainStripMatched.push_back(hit);
      }
    }
  }

  delete fitStrip;
  delete fitChain;

  delete hStrip;
  delete hChain;
}

void Spectra::ChainStripMatchingTimeSlopeHough(std::vector<mmTrack> &chainStripMatched, std::vector<mmChainStrip> chain,
                                               std::vector<mmChainStrip> strip, Bool_t leftSide, Double_t siTime,
                                               Double_t timeWindow) {
  // Function to match strips and chains by fitting with time individually (using Hough 2D)
  // Then is able to correlate strips and chains by time. Requires multiple time steps in both strips and chains
  // After the fit, it goes through each chain and finds strips based off of the fit when match within the timeWindow input parameter
  // Position transform = 10.5 + 1.75/2 + 1.75*chain (chain starting at 0)
  // Row transform = (strip - 1)*2
  // Row position transform = Row * rowConversion + rowConversionOffset

  std::vector<mmTrack> newChain;
  for(auto mm : chain) {
    mmTrack hit = {0, 0., 0., 0., 0., 0., 0};
    hit.xPosition = mm.row;
    hit.yPosition = mm.time/timeResolution; // Convert to time buckets
    newChain.push_back(hit);
  }

  std::vector<mmTrack> newStrip;
  for(auto mm : strip) {
    mmTrack hit = {0, 0., 0., 0., 0., 0., 0};
    hit.xPosition = mm.row;
    hit.yPosition = mm.time/timeResolution; // Convert to time buckets
    newStrip.push_back(hit);
  }

  auto *houghChain = new Hough2D();
  houghChain->SetPoints(newChain);
  houghChain->CalculateHoughXY();
  Double_t chainTheta = houghChain->GetMaxThetaXY();
  Double_t chainD = houghChain->GetMaxDXY();
  delete houghChain;

  auto *houghStrip = new Hough2D();
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
        chainStripMatched.push_back(hit);
      }
    }
  }
}

void Spectra::DivideTargetThickness(TH1F *f) {
  Int_t i_size = f->GetSize();

  TAxis *x_axis = f->GetXaxis();
  for(Int_t i = 1; i <= i_size - 1; i++) {
    Double_t binLowEdge = x_axis->GetBinLowEdge(i);
    Double_t binUpEdge = x_axis->GetBinUpEdge(i);
    if(binLowEdge == 0.) binLowEdge += 0.001;
    binLowEdge *= (m1 + m2)/m2;
    binUpEdge *= (m1 + m2)/m2;
    Double_t binContent = f->GetBinContent(i);
    Double_t binError = f->GetBinError(i);
    Double_t delta_x = boronMethane->CalcRange(binUpEdge, binLowEdge);
    delta_x /= 10.;
    Double_t molarMassMethane = 0.01604;
    Double_t factor = 4.e-27*density*delta_x*TMath::Na()/molarMassMethane;
    binContent /= factor;
    binError /= factor;
    f->SetBinContent(i, binContent);
    f->SetBinError(i, binError);
  }
}

void Spectra::ReadSolidAngle() {
  std::ifstream in("csReg3.out");
  assert(in.is_open());

  Double_t var1, var2, var3, var4;
  std::vector<Double_t> cmEnergy_, solidAngle_, labAngle_, cmAngle_;
  while(in >> var1 >> var2 >> var3 >> var3) {
    cmEnergy_.push_back(var1);
    solidAngle_.push_back(var2);
    labAngle_.push_back(var3);
    cmAngle_.push_back(var4);
  }

  reg3SA.SetPoints(cmEnergy_, solidAngle_);
  reg3CMAngle.SetPoints(cmEnergy_, cmAngle_);

  in.close();
}

void Spectra::SolidAngle(TH1F *f) {
  Int_t i_size = f->GetSize();
  TAxis *x_axis = f->GetXaxis();

  for(Int_t i = 1; i <= i_size - 1; i++) {
    Double_t binCenter = x_axis->GetBinCenter(i);
    Double_t binContent = f->GetBinContent(i);
    Double_t binError = f->GetBinError(i);

    Double_t simCS = reg3SA(binCenter);
    // std::cout << binCenter << '\t' << simCS << std::endl;

    if(simCS == 0) {
      f->SetBinContent(i, 0.);
      f->SetBinError(i, 0.);
    }
    else {
      f->SetBinContent(i, binContent/simCS);
      f->SetBinError(i, binError/simCS);
    }
  }
}

void Spectra::WriteSpectrumToFile(TH1F *f, Int_t region) {
  FILE* spectrumFile = fopen(Form("spectrum_reg%d.out", region), "w");

  Int_t i_size = f->GetSize();
  TAxis *x_axis = f->GetXaxis();

  for(Int_t i = 1; i <= i_size - 1; i++) {
    Double_t binCenter = x_axis->GetBinCenter(i);
    Double_t binContent = f->GetBinContent(i);
    Double_t binError = f->GetBinError(i);

    Double_t cmAngle = reg3CMAngle(binCenter);

    fprintf(spectrumFile, "%f %f %f %f\n", binCenter, binContent, binError, cmAngle);
  }

  fflush(spectrumFile);
  fclose(spectrumFile);
}

void Spectra::GetMinMaxD(std::vector<mmTrack> initPoints, Int_t &minXY, Int_t &maxXY, Int_t &minYZ, Int_t &maxYZ) {
  std::vector<xy> pointsXY;
  std::vector<yz> pointsYZ;
  for(UInt_t i = 0; i < initPoints.size(); i++) {
    xy initPointsXY = {initPoints[i].xPosition, initPoints[i].yPosition};
    yz initPointsYZ = {initPoints[i].yPosition, initPoints[i].height};
    pointsXY.push_back(initPointsXY);
    pointsYZ.push_back(initPointsYZ);
  }

  Double_t minDHistogramXY = 1000;
  Double_t maxDHistogramXY = -1000;
  Double_t minDHistogramYZ = 1000;
  Double_t maxDHistogramYZ = -1000;
  for(Int_t j = 0; j < 180; j++) {
    Double_t cosj = cos(j*M_PI/180.);
    Double_t sinj = sin(j*M_PI/180.);
    for(UInt_t i = 0; i < pointsXY.size(); i++) {
      Double_t d = pointsXY[i].x*cosj + pointsXY[i].y*sinj;
      if(d < minDHistogramXY) minDHistogramXY = d;
      if(d > maxDHistogramXY) maxDHistogramXY = d;
    }
    for(UInt_t i = 0; i < pointsYZ.size(); i++) {
      Double_t d = pointsYZ[i].y*cosj + pointsYZ[i].z*sinj;
      if(d < minDHistogramYZ) minDHistogramYZ = d;
      if(d > maxDHistogramYZ) maxDHistogramYZ = d;
    }
  }

  minXY = static_cast<Int_t>(minDHistogramXY);
  maxXY = static_cast<Int_t>(maxDHistogramXY);
  minYZ = static_cast<Int_t>(minDHistogramYZ);
  maxYZ = static_cast<Int_t>(maxDHistogramYZ);
}

void Spectra::VisualizeHough(std::vector<mmTrack> initPoints, TH2I* fXY, TH2I* fYZ) {
  Int_t nBinsX = 500;
  Int_t nBinsY = 500;
  std::vector<xy> pointsXY;
  std::vector<yz> pointsYZ;
  for(UInt_t i = 0; i < initPoints.size(); i++) {
    xy initPointsXY = {initPoints[i].xPosition, initPoints[i].yPosition};
    yz initPointsYZ = {initPoints[i].yPosition, initPoints[i].height};
    pointsXY.push_back(initPointsXY);
    pointsYZ.push_back(initPointsYZ);
  }

  Double_t thetaStepXY = 180./static_cast<Double_t>(nBinsX);

  // Fill Hough Matrix
  for(Double_t j = 0; j < 180; j += thetaStepXY) {
    if(j > 89 && j < 91) continue;
    Double_t cosj = cos(j*M_PI/180.);
    Double_t sinj = sin(j*M_PI/180.);
    for(UInt_t i = 0; i < pointsXY.size(); i++) {
      Double_t d = pointsXY[i].x*cosj + pointsXY[i].y*sinj;
      fXY->Fill(j, d);
    }
  }
}

void Spectra::GetMinMaxDRestricted(std::vector<mmTrack> initPoints, Int_t &minXY, Int_t &maxXY, Int_t &minYZ, Int_t &maxYZ, Int_t siDet) {
  Int_t minAngle = 0;
  Int_t maxAngle = 0;
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
  for(UInt_t i = 0; i < initPoints.size(); i++) {
    xy initPointsXY = {initPoints[i].xPosition, initPoints[i].yPosition};
    yz initPointsYZ = {initPoints[i].yPosition, initPoints[i].height};
    pointsXY.push_back(initPointsXY);
    pointsYZ.push_back(initPointsYZ);
  }

  Double_t minDHistogramXY = 1000;
  Double_t maxDHistogramXY = -1000;
  Double_t minDHistogramYZ = 1000;
  Double_t maxDHistogramYZ = -1000;
  for(Int_t j = minAngle; j < maxAngle; j++) {
    Double_t cosj = cos(j*M_PI/180.);
    Double_t sinj = sin(j*M_PI/180.);
    for(UInt_t i = 0; i < pointsXY.size(); i++) {
      Double_t d = pointsXY[i].x*cosj + pointsXY[i].y*sinj;
      if(d < minDHistogramXY) minDHistogramXY = d;
      if(d > maxDHistogramXY) maxDHistogramXY = d;
    }
    for(UInt_t i = 0; i < pointsYZ.size(); i++) {
      Double_t d = pointsYZ[i].y*cosj + pointsYZ[i].z*sinj;
      if(d < minDHistogramYZ) minDHistogramYZ = d;
      if(d > maxDHistogramYZ) maxDHistogramYZ = d;
    }
  }

  minXY = static_cast<Int_t>(minDHistogramXY);
  maxXY = static_cast<Int_t>(maxDHistogramXY);
  minYZ = static_cast<Int_t>(minDHistogramYZ);
  maxYZ = static_cast<Int_t>(maxDHistogramYZ);
}

void Spectra::VisualizeHoughRestricted(std::vector<mmTrack> initPoints, TH2I* fXY, TH2I* fYZ, Int_t siDet) {
  Int_t minAngle = 0;
  Int_t maxAngle = 180;
  if(siDet < 4) {
    minAngle = 91;
    maxAngle = 179;
  }
  else if(siDet > 5 && siDet < 10) {
    minAngle = 1;
    maxAngle = 89;
  }

  Int_t nBinsX = 360;
  Int_t nBinsY = 360;
  std::vector<xy> pointsXY;
  std::vector<yz> pointsYZ;
  for(auto point : initPoints) {
    xy initPointsXY = {point.xPosition, point.yPosition};
    yz initPointsYZ = {point.yPosition, point.height};
    pointsXY.push_back(initPointsXY);
    pointsYZ.push_back(initPointsYZ);
  }

  Double_t thetaStepXY = 180./static_cast<Double_t>(nBinsX);

  // Fill Hough Matrix
  for(Double_t j = minAngle; j < maxAngle; j += thetaStepXY) {
    if(j > 89 && j < 91) continue;
    Double_t cosj = cos(j*M_PI/180.);
    Double_t sinj = sin(j*M_PI/180.);
    for(auto point : pointsXY) {
      Double_t d = point.x*cosj + point.y*sinj;
      fXY->Fill(j, d);
    }
  }
}

void Spectra::GetHoughStdDevXYRestricted(std::vector<mmTrack> initPoints, std::vector<Double_t> &angle_, std::vector<Double_t> &stdDev_, Int_t siDet) {

  angle_.clear();
  stdDev_.clear();

  Int_t minAngle = 0;
  Int_t maxAngle = 180;
  if(siDet < 4) {
    minAngle = 91;
    maxAngle = 179;
  }
  else if(siDet > 5 && siDet < 10) {
    minAngle = 1;
    maxAngle = 89;
  }

  Int_t nBinsX = 360;
  Int_t nBinsY = 360;
  std::vector<xy> pointsXY;
  std::vector<yz> pointsYZ;
  for(UInt_t i = 0; i < initPoints.size(); i++) {
    xy initPointsXY = {initPoints[i].xPosition, initPoints[i].yPosition};
    yz initPointsYZ = {initPoints[i].yPosition, initPoints[i].height};
    pointsXY.push_back(initPointsXY);
    pointsYZ.push_back(initPointsYZ);
  }

  Double_t thetaStepXY = 180./static_cast<Double_t>(nBinsX);

  // Fill Hough Matrix
  for(Double_t j = minAngle; j < maxAngle; j += thetaStepXY) {
    if(j > 89 && j < 91) continue;
    Double_t cosj = cos(j*M_PI/180.);
    Double_t sinj = sin(j*M_PI/180.);
    std::vector<Double_t> dVector;
    Double_t mean = 0.;
    for(UInt_t i = 0; i < pointsXY.size(); i++) {
      Double_t d = pointsXY[i].x*cosj + pointsXY[i].y*sinj;
      mean += d;
      dVector.push_back(d);
    }
    if(dVector.empty()) continue;

    mean /= static_cast<Double_t>(pointsXY.size());
    Double_t stdDev = 0.;
    for(UInt_t i = 0; i < pointsXY.size(); i++) {
      Double_t d = pointsXY[i].x*cosj + pointsXY[i].y*sinj;
      stdDev += (d - mean)*(d - mean);
    }
    stdDev /= static_cast<Double_t>(pointsXY.size() - 1.);

    angle_.push_back(j);
    stdDev_.push_back(stdDev);

  }

}