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
  TString PathToFiles = "/hd3/research/data/run0817a/rootM2R-WaveformReduced/run";

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

  if (fChain == 0) return;

  /////////////////
  // Set up cuts //
  /////////////////
  TFile *cutFile = TFile::Open("cuts.root");
  TCutG *dEEForwardCut[10];
  for(UInt_t i = 0; i < 10; i++) {
    if(i == 8) continue;
    dEEForwardCut[i] = (TCutG*)cutFile->Get(Form("dEEForward_d%dCut", i));
  }
  TCutG *angleTotEnergyCut_d0 = (TCutG*)cutFile->Get("angleTotEnergy_d0Cut");
  TCutG *angleTotEnergyCut_d1 = (TCutG*)cutFile->Get("angleTotEnergy_d1Cut");
  TCutG *angleTotEnergyCut_d9 = (TCutG*)cutFile->Get("angleTotEnergy_d9Cut");
  TCutG *angleTotEnergyCut_d019 = (TCutG*)cutFile->Get("angleTotEnergy_d019Cut");
  cutFile->Close();

  InitChannelMap();
  InitHistograms();
  InitVariables();

  InitSiEForwardCalibration();
  InitCsIECalibration();

  InitCentralPadGainMatch();
  InitAverageBeamEnergy();

  Bool_t individualMMHistograms = false;

  InitTree();

  Long64_t nentries = fChain->GetEntriesFast();

  printf("Starting Main Loop\n");

  Long64_t nbytes = 0, nb = 0;
  // for(Long64_t jentry = 0; jentry < 4000; jentry++) {
  // for(Long64_t jentry = 3431; jentry < 3432; jentry++) {
  for(Long64_t jentry = 0; jentry < nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if(ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;

    if(jentry != 0 && jentry % 2000 == 0) printf("Processed %lld events\n",jentry);

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

    hIonizationChamberE->Fill(icE);
    hIonizationChamberT->Fill(icT);

    // if(icE < 1400 || icE > 1900) continue;
    if(icT < 5500 || icT > 6500) continue;

    if(siDetect_.size() != 1) continue;

    siDet = siDetect_[0].detect;
    siQuad = siDetect_[0].quad;
    siChannel = siDetect_[0].channel;
    siEnergy = siDetect_[0].energy;
    siEnergyCal = 0.;
    if(siDet < 10) {
      siEnergyCal = siEnergy*siEForwardCalibration[siDet][siQuad].first + siEForwardCalibration[siDet][siQuad].second;
    }
    else {
      siEnergyCal = siEnergy*siELeftCalibration[siDet - 10][siQuad].first + siELeftCalibration[siDet - 10][siQuad].second;
    }
    siTime = siDetect_[0].time;

    if(siEnergyCal < 1) continue;

    if(siDet < 10) {
      hSiEForwardDet[siDet]->Fill(siEnergy);
      hSiEForwardDetCal[siDet]->Fill(siEnergyCal);
      hSiTForwardDet[siDet]->Fill(siTime);
      hSiEForward[siDet][siQuad]->Fill(siEnergy);
      hSiEForwardCal[siDet][siQuad]->Fill(siEnergyCal);
    }
    else {
      hSiELeftDet[siDet - 10]->Fill(siEnergy);
      hSiELeftDetCal[siDet - 10]->Fill(siEnergyCal);
      hSiTLeftDet[siDet - 10]->Fill(siTime);
      hSiELeft[siDet - 10][siQuad]->Fill(siEnergy);
      hSiELeftCal[siDet - 10][siQuad]->Fill(siEnergyCal);
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

    Bool_t forwardDet = false;
    Bool_t sideDet = false;
    if(siDet < 10) forwardDet = true;
    else sideDet = true;

    if(sideDet) continue;

    // Find if CsI behind Si fired
    punchthrough = false;
    csiEnergy = 0.;
    csiEnergyCal = 0.;
    csiTime = 0.;
    for(UInt_t i = 0; i < csiDetect_.size(); i++) {
      if(csiDetect_[i].detect == siDet) {
        csiEnergy = csiDetect_[i].energy;
        csiTime = csiDetect_[i].time;
        if(siDet < 10) hCsIETForward[siDet]->Fill(csiEnergy, csiTime);
        else hCsIETLeft[siDet - 10]->Fill(csiEnergy, csiTime);
        if(csiTime < 11000 || csiTime > 12250) continue;
        punchthrough = true;
        if(siDet < 10) csiEnergyCal = csiEnergy*csiEForwardCalibration[siDet].first + csiEForwardCalibration[siDet].second;
        hSiCsIEForwardDet[siDet]->Fill(siEnergy, csiEnergy);
        hSiCsIEForwardDetCal[siDet]->Fill(siEnergyCal, csiEnergy);

        hSiCsIEForward[siDet][siQuad]->Fill(siEnergy, csiEnergy);
        hSiCsIEForwardCal[siDet][siQuad]->Fill(siEnergyCal, csiEnergy);
      }
    }

    if(punchthrough) {
      hSumSiEForwardDet[siDet]->Fill(siEnergy + csiEnergy, siEnergy);
      hSumSiEForward[siDet][siQuad]->Fill(siEnergy + csiEnergy, siEnergy);

      hSumCsIEForwardDet[siDet]->Fill(siEnergy + csiEnergy, csiEnergy);
      hSumCsIEForward[siDet][siQuad]->Fill(siEnergy + csiEnergy, csiEnergy);
    }

    totalEnergy = siEnergyCal + csiEnergyCal;

    // Reduce MM Center to one entry per row
    std::map<Int_t, Double_t> centralPadPosition;
    std::map<Int_t, Double_t> centralPadTotalEnergy;
    std::map<Int_t, Double_t> centralPadTime;
    std::map<Int_t, Int_t> centralPadColumn;
    std::map<Int_t, Int_t> centralPadTotal;
    for(UInt_t i = 0; i < 128; i++) {
      centralPadTotal[i] = 0;
    }
    for(UInt_t i = 0; i < mmCenterMatched_.size(); i++) {
      if(centralPadTotal[mmCenterMatched_[i].row] == 0) {
        centralPadPosition[mmCenterMatched_[i].row] = (mmCenterMatched_[i].column*3.5 - 8.75)*mmCenterMatched_[i].energy;
        centralPadTotalEnergy[mmCenterMatched_[i].row] = mmCenterMatched_[i].energy;
        centralPadTime[mmCenterMatched_[i].row] = mmCenterMatched_[i].time;
        centralPadColumn[mmCenterMatched_[i].row] = mmCenterMatched_[i].column;
      }
      else {
        centralPadPosition[mmCenterMatched_[i].row] += (mmCenterMatched_[i].column*3.5 - 8.75)*mmCenterMatched_[i].energy;
        centralPadTotalEnergy[mmCenterMatched_[i].row] += mmCenterMatched_[i].energy;
        centralPadTime[mmCenterMatched_[i].row] += mmCenterMatched_[i].time;
      }
      centralPadTotal[mmCenterMatched_[i].row]++;
    }

    std::vector<mmTrack> mmCenterBeamTotal;
    std::vector<mmTrack> mmCenterProton;
    std::map<Int_t, Double_t>::iterator it;
    for(it = centralPadPosition.begin(); it != centralPadPosition.end(); it++) {
      Int_t row = it->first;
      mmTrack hit = {0, 0., 0., 0., 0., 0., 0};
      if(centralPadTotalEnergy[row] == 0) continue;
      hit.row = row;
      hit.xPosition = centralPadPosition[row]/centralPadTotalEnergy[row];
      hit.yPosition = row*rowConversion + rowConversionOffset;
      hit.time = centralPadTime[row]/static_cast<Double_t>(centralPadTotal[row]) - siTime;
      hit.energy = centralPadTotalEnergy[row];
      hit.height = heightOffset - hit.time*driftVelocity;
      hit.total = centralPadTotal[row];
      hMicroMegasCenterCumulative->Fill(centralPadColumn[row], row);
      hMicroMegasCenterCumulativePosition->Fill(centralPadColumn[row], hit.yPosition);
      if(row < 112) {
        mmCenterBeamTotal.push_back(hit);
      }
      else if(row > 111) {
        mmCenterProton.push_back(hit);
      }
    }

    // Sorting mmCenterBeamTotal by row
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

    std::vector<mmTrack> stripChainMatchedLeft;
    std::vector<mmTrack> stripChainMatchedRight;
    std::vector<mmTrack> stripChainRawLeft;
    std::vector<mmTrack> stripChainRawRight;
    if(mmLeftChain_.size() > 0 && mmLeftStrip_.size() > 0) {
      StripChainMatch(stripChainMatchedLeft, stripChainRawLeft, mmLeftChain_, mmLeftStrip_, true, siTime);
    }
    if(mmRightChain_.size() > 0 && mmRightStrip_.size() > 0) {
      StripChainMatch(stripChainMatchedRight, stripChainRawRight, mmRightChain_, mmRightStrip_, false, siTime);
    }

    // Find dE
    dE = 0.;
    if(central) { // Forward central detectors
      Int_t totalRows = 0;
      for(Int_t i = 116; i < 124; i++) {
        if(i == 117) continue;
        auto it = centralPadTotalEnergy.find(i);
        if(it != centralPadTotalEnergy.end()) {
          dE += centralPadTotalEnergy[i];
          totalRows++;
        }
      }
      dE /= totalRows;
    }
    else if(right) { // Forward beam right detectors
      for(UInt_t i = 0; i < mmRightStrip_.size(); i++) {
        dE += mmRightStrip_[i].energy;
      }
      dE /= mmRightStrip_.size();
    }
    else if(left && siDet < 10) { // Forward beam left detectors
      for(UInt_t i = 0; i < mmLeftStrip_.size(); i++) {
        dE += mmLeftStrip_[i].energy;
      }
      dE /= mmLeftStrip_.size();
    }
    else if(left && siDet > 9) { // Left wall detectors
      for(UInt_t i = 0; i < mmLeftStrip_.size(); i++) {
        dE += mmLeftStrip_[i].energy;
      }
      dE /= mmLeftStrip_.size();
    }

    hdEEForward[siDet]->Fill(siEnergy, dE);
    hdEEForwardCal[siDet]->Fill(siEnergyCal, dE);
    hdEEForwardCalTotal[siDet]->Fill(siEnergyCal + csiEnergyCal, dE);

    if(!dEEForwardCut[siDet]->IsInside(siEnergy, dE)) continue;

    // Assign proton track for side regions
    std::vector<mmTrack> protonTrack;
    if(left) {
      protonTrack = stripChainMatchedLeft;
    }
    else if(right) {
      protonTrack = stripChainMatchedRight;
    }

    // Find vertex using the proton track in side regions
    vertexPositionX = 0.;
    vertexPositionY = -400.;
    vertexPositionZ = 0.;

    Double_t houghAngleXY = 89.;

    if(!central) {

      if(protonTrack.size() == 0) {
        std::cout << siEnergy << '\t' << dE << std::endl;
        continue;
      }
      // for(UInt_t i = 0; i < protonTrack.size(); i++) {
        // printf("%f %f %f\n", protonTrack[i].xPosition, protonTrack[i].yPosition, protonTrack[i].height);
      // }

      // std::cout << jentry << '\t' << siDet << '\t' << mmRightChain_.size() << '\t' << mmRightStrip_.size() << '\t' << protonTrack.size() << std::endl;


      mmTrack detPoint1, detPoint2, detPoint3;
      if(siDet < 4) {
        detPoint1 = {0, siXPosForward[siDet][siQuad] + 12.5, 275.34, 0, 0, 0, 0};
        detPoint2 = {0, siXPosForward[siDet][siQuad], 275.34, 0, 0, 0, 0};
        detPoint3 = {0, siXPosForward[siDet][siQuad] - 12.5, 275.34, 0, 0, 0, 0};
        for(UInt_t j = 0; j < 1; j++) {
          protonTrack.push_back(detPoint1);
          protonTrack.push_back(detPoint2);
          protonTrack.push_back(detPoint3);
        }
      }
      else if(siDet > 5 && siDet < 10) {
        detPoint1 = {0, -siXPosForward[siDet][siQuad] + 12.5, 275.34, 0, 0, 0, 0};
        detPoint2 = {0, -siXPosForward[siDet][siQuad], 275.34, 0, 0, 0, 0};
        detPoint3 = {0, -siXPosForward[siDet][siQuad] - 12.5, 275.34, 0, 0, 0, 0};
        for(UInt_t j = 0; j < 1; j++) {
          protonTrack.push_back(detPoint1);
          protonTrack.push_back(detPoint2);
          protonTrack.push_back(detPoint3);
        }
      }

      HoughTrack* fitProton = new HoughTrack();
      fitProton->AddTrack(protonTrack, siDet, siQuad);
      Double_t minDist = fitProton->FitRestricted();
      if(minDist > 3.e+07) {
        for(UInt_t j = 0; j < 34; j++) {
          protonTrack.push_back(detPoint1);
          protonTrack.push_back(detPoint2);
          protonTrack.push_back(detPoint3);
        }
        fitProton->AddTrack(protonTrack, siDet, siQuad);
        minDist = fitProton->FitRestricted();
      }
      std::vector<Double_t> parsProton = fitProton->GetPars();
      houghAngleXY = fitProton->GetHoughAngleXY();
      hHoughAngle[siDet]->Fill(houghAngleXY);

      delete fitProton;

      // Check there were beam ions
      UInt_t mmCenterSize = mmCenterBeamTotal.size();

      // std::cout << parsProton[0] << '\t' << parsProton[1] << std::endl;
      // std::cout << parsProton[2] << '\t' << parsProton[3] << std::endl;
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
        Double_t xDiff = xDiff;
        for(Double_t yPos = beamY_old; yPos > -300; yPos -= 1.) {
          beamY = yPos;
          line(beamY, parsProton, x, y, z);
          protonX = x;
          xDiff = beamX - protonX;
          if(xDiff_old * xDiff_old < 0) {
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
      else{
        // Loop through beam in central region starting from last
        Double_t beamX_old = mmCenterBeamTotal[mmCenterSize - 1].xPosition;
        Double_t beamY_old = mmCenterBeamTotal[mmCenterSize - 1].yPosition;
        Double_t beamZ_old = mmCenterBeamTotal[mmCenterSize - 1].height;
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
        // std::cout << jentry << '\t' << mmCenterBeamTotal.size() << '\t' << mmCenterSize - 1 << std::endl;
        for(Int_t i = mmCenterSize - 1; i > -1; i--) {
          beamX = mmCenterBeamTotal[i].xPosition;
          beamY = mmCenterBeamTotal[i].yPosition;
          beamZ = mmCenterBeamTotal[i].height;
          line(beamY, parsProton, x, y, z);
          protonX = x;
          xDiff = beamX - protonX;
          // printf("%lld %f %f %f %f %f %f\n", jentry, beamX_old, beamY_old, beamX, beamY, xDiff, xDiff_old);
          if(xDiff_old * xDiff < 0) {
            Double_t m = (beamY - beamY_old)/(xDiff - xDiff_old);
            Double_t b = beamY - m*xDiff;
            // printf("%lld %f %f %f %f %f %f %f\n", jentry, beamX_old, beamY_old, beamX, beamY, protonX, protonX_old, b);
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
          for(Double_t beamY = -1; beamY > -300; beamY -= 1) {
            line(beamY, parsProton, x, y, z);
            protonX = x;
            xDiff = beamX - x;
            // printf("%lld %f %f %f %f %f %f\n", jentry, beamX_old, beamY_old, beamX, beamY, xDiff, xDiff_old);
            if(xDiff * xDiff_old < 0) {
              Double_t m = (beamY - beamY_old)/(xDiff - xDiff_old);
              Double_t b = beamY - m*xDiff;
              // printf("%lld %f %f %f %f %f %f %f\n", jentry, beamX_old, beamY_old, beamX, beamY, protonX, protonX_old, b);
              vertexPositionX = beamX;
              vertexPositionY = b;
              vertexPositionZ = 0.;
              foundVertex = true;
              break;
            }
            beamY_old = beamY;
            protonX_old = protonX;
            xDiff_old = xDiff;
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




      // Find the angle (use 2D angle for now)
      // Simply use the slope of the x-component and assume the beam is straight
      // The angle is then theta = atan(|m|) where m in the slope of the x-component (parsProton[1])
      Float_t m_xcomponent = fabs(parsProton[1]);
      angle = atan(m_xcomponent);

      TVector3 v1(0, 1, 0);
      TVector3 v2(parsProton[1], 1, parsProton[3]);
      // Float_t angle3d = v1.Angle(v2);
      // angle = v1.Angle(v2);

      // std::cout << jentry << '\t' << angle << '\t' << angle3d << std::endl;

      // printf("%lld %d %f\n", jentry, siDet, vertexPositionY);

      // if(houghAngleXY > 91) {

      //   Int_t nBinsX = 360;
      //   Int_t nBinsY = 360;
      //   Int_t minXY, maxXY, minYZ, maxYZ;
      //   // GetMinMaxD(protonTrack, minXY, maxXY, minYZ, maxYZ);
      //   GetMinMaxDRestricted(protonTrack, minXY, maxXY, minYZ, maxYZ, siDet);

      //   TH2I* houghXY = new TH2I(Form("houghXY_Event_%lld", jentry), Form("houghXY_Event_%lld", jentry), nBinsX, 0, 180, nBinsY, minXY, maxXY);
      //   TH2I* houghYZ = new TH2I(Form("houghYZ_Event_%lld", jentry), Form("houghYZ_Event_%lld", jentry), nBinsX, 0, 180, nBinsY, minYZ, maxYZ);

      //   // VisualizeHough(protonTrack, houghXY, houghYZ);
      //   VisualizeHoughRestricted(protonTrack, houghXY, houghYZ, siDet);
      //   houghXY->Write();

      //   std::vector<Double_t> houghAngleXY_, houghStdDevXY_;
      //   GetHoughStdDevXYRestricted(protonTrack, houghAngleXY_, houghStdDevXY_, siDet);
      //   TGraph *houghXYStdDevGraph = new TGraph(houghAngleXY_.size(), &houghAngleXY_[0], &houghStdDevXY_[0]);
      //   houghXYStdDevGraph->SetName(Form("Hough_StdDev_Event_%lld", jentry));
      //   houghXYStdDevGraph->Write();

      //   delete houghXY;
      //   delete houghYZ;

      //   TGraph *hTrackProtonRaw = new TGraph();
      //   for(UInt_t i = 0; i < protonTrack.size(); i++) {
      //     hTrackProtonRaw->SetPoint(i, protonTrack[i].xPosition, protonTrack[i].yPosition);
      //   }
      //   hTrackProtonRaw->SetMarkerStyle(20);
      //   hTrackProtonRaw->SetMarkerColor(2);
      //   hTrackProtonRaw->SetName(Form("Proton_Raw_Track_%lld", jentry));
      //   hTrackProtonRaw->Write();
      //   delete hTrackProtonRaw;

      //   Int_t n = 6000;
      //   Double_t t0 = -300.;
      //   Double_t dt = 300 - t0;
      //   TGraph *hTrackProtonFit = new TGraph();
      //   for(Int_t i = 0; i < n; i++) {
      //     Double_t t = t0 + dt*i/n;
      //     Double_t x, y, z;
      //     line(t, parsProton, x, y, z);
      //     // std::cout << i << '\t' << x << '\t' << y << '\t' << z << std::endl;
      //     hTrackProtonFit->SetPoint(i, x, y);
      //   }
      //   hTrackProtonFit->SetMarkerStyle(20);
      //   hTrackProtonFit->SetMarkerColor(3);
      //   hTrackProtonFit->SetName(Form("Proton_Fit_Track_%lld", jentry));
      //   hTrackProtonFit->Write();
      //   delete hTrackProtonFit;
      // }

      // TGraph *hTrackVertex = new TGraph();
      // hTrackVertex->SetPoint(0, vertexPositionX, vertexPositionY);
      // hTrackVertex->SetMarkerStyle(20);
      // hTrackVertex->SetMarkerColor(4);
      // hTrackVertex->SetName(Form("Vertex_%lld", jentry));
      // hTrackVertex->Write();
      // delete hTrackVertex;

      // TGraph *hTrackBeam = new TGraph();
      // for(UInt_t i = 0; i < mmCenterBeamTotal.size(); i++) {
      //   hTrackBeam->SetPoint(i, mmCenterBeamTotal[i].xPosition, mmCenterBeamTotal[i].yPosition);
      // }
      // hTrackBeam->SetMarkerStyle(20);
      // hTrackBeam->SetMarkerColor(1);
      // hTrackBeam->SetName(Form("Beam_Track_%lld", jentry));
      // hTrackBeam->Write();
      // delete hTrackBeam;
    }

    if(siDet < 10) {
      hVertexSiEForward[siDet]->Fill(siEnergy, vertexPositionY);
      hVertexSiEForwardCal[siDet]->Fill(siEnergyCal, vertexPositionY);

      hAngleEForward[siDet]->Fill(siEnergy, angle);
      hAngleEForwardCal[siDet]->Fill(siEnergyCal, angle);
      hAngleEForwardCalTotal[siDet]->Fill(totalEnergy, angle);

      hVertexAngleForward[siDet]->Fill(vertexPositionY, angle);
    }

    // Simply calculate the CS. Just using detectors 0, 1 and 9. Average x position is 124 mm.
    cmEnergy = 0.;
    if(siDet == 0 || siDet == 1 || siDet == 9) {
      // if(!dEEForwardCut[siDet]->IsInside(siEnergy, dE)) continue;
      // if(siDet == 0) {
      //   if(!angleTotEnergyCut_d0->IsInside(totalEnergy, angle)) continue;
      // }
      // else if(siDet == 1) {
      //   if(!angleTotEnergyCut_d1->IsInside(totalEnergy, angle)) continue;
      // }
      // else if(siDet == 9) {
      //   if(!angleTotEnergyCut_d9->IsInside(totalEnergy, angle)) continue;
      // }

      // if(!angleTotEnergyCut_d019->IsInside(totalEnergy, angle)) continue;

      Float_t protonPathLength = sqrt(124*124 + vertexPositionY*vertexPositionY);
      Float_t protonEnergy = protonMethane->AddBack(totalEnergy/1000., protonPathLength);
      Float_t cosAngle2 = cos(angle)*cos(angle);
      Float_t beamEnergy = protonEnergy*(m1 + m2)*(m1 + m2)/(4.*m1*m2*cosAngle2);
      cmEnergy = beamEnergy*m2/(m1 + m2);

      s1->Fill(cmEnergy);
    }

    //  ** End of event by event analysis ** //
    FillTree();
  }

  // TGraph boundaries
  TGraph *h_track_bound = new TGraph();
  h_track_bound->SetMarkerStyle(8);
  h_track_bound->SetMarkerColor(5);
  h_track_bound->SetPoint(0, -150, -300);
  h_track_bound->SetPoint(1, 150, 300);
  h_track_bound->SetName("Bound");
  h_track_bound->Write();
  delete h_track_bound;

  hIonizationChamberE->Write();
  hIonizationChamberT->Write();

  hMicroMegasCenterCumulative->Write();
  hMicroMegasCenterCumulativePosition->Write();

  // Forward Si Detectors
  for(UInt_t i = 0; i < 10 ; i++) {
    // hSiEForwardDet[i]->Write();
    // hSiTForwardDet[i]->Write();
    hSiEForwardDetCal[i]->Write();
    // for(int j=0; j<4; j++) {
      // hSiEForward[i][j]->Write();
      // hSiEForwardCal[i][j]->Write();
      // hSiTForward[i][j]->Write();
    // }
    // hCsIEForward[i]->Write();
    // hCsITForward[i]->Write();
    // hSiCsIForward[i]->Write();
  }

  // Left Si Detectors
  // for(UInt_t i = 0; i < 6 ; i++) {
    // hSiELeftDet[i]->Write();
    // hSiTLeftDet[i]->Write();
    // hSiELeftDetCal[i]->Write();
    // for(int j=0; j<4; j++) {
      // hSiELeft[i][j]->Write();
      // hSiTLeft[i][j]->Write();
    // }
    // hCsIELeft[i]->Write();
    // hCsITLeft[i]->Write();
    // hCsIETLeft[i]->Write();
  // }

  // hMicroMegasCenterCumulative->Write();
  // hMicroMegasCenterEnergyCumulative->Write();
  // hMicroMegasCenterEnergyAverage->Write();
  // hMicroMegasCenterEnergyAverageScaled->Write();
  // hMicroMegasCenterTimeAverage->Write();

  // Forward CsI Energy vs Time
  // for(UInt_t i = 0; i < 10; i++) {
  //   hCsIETForward[i]->Write();
  // }

  // Forward Si Energy vs CsI Energy
  // for(UInt_t i = 0; i < 10; i++) {
    // hSiCsIEForwardDet[i]->Write();
    // hSiCsIEForwardDetCal[i]->Write();
    // for(UInt_t j = 0; j < 4; j++) {
    //   hSiCsIEForward[i][j]->Write();
    //   hSiCsIEForwardCal[i][j]->Write();
    // }
  // }

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
  // for(UInt_t i = 0; i < 10; i++) {
    // hdEEForward[i]->Write();
    // hdEEForwardCal[i]->Write();
    // hdEEForwardCalTotal[i]->Write();
  // }

  // Forward Hough Angle
  for(UInt_t i = 0; i < 10; i++) {
    hHoughAngle[i]->Write();
  }

  // Forward Vertex vs Si Energy
  for(UInt_t i = 0; i < 10; i++) {
    hVertexSiEForward[i]->Write();
    hVertexSiEForwardCal[i]->Write();
  }

  // Forward Angle vs Si Energy
  for(UInt_t i = 0; i < 10; i++) {
    hAngleEForward[i]->Write();
    hAngleEForwardCal[i]->Write();
    hAngleEForwardCalTotal[i]->Write();
  }

  // Forward Vertex vs Angle
  // for(UInt_t i = 0; i < 10; i++) {
  //   hVertexAngleForward[i]->Write();
  // }

  DivideTargetThickness(s1);
  ReadSolidAngle();
  SolidAngle(s1);
  s1->Scale(1./numberB8);
  // s1->Write();
  WriteSpectrumToFile(s1, 3);

  WriteTree();

  file->Close();
}

void Spectra::StripChainMatch(std::vector<mmTrack> &stripChainMatched, std::vector<mmTrack> &stripChainRaw, std::vector<mmStripChain> chain_,
                              std::vector<mmStripChain> strip_, Bool_t leftSide, Double_t siTime) {
  stripChainMatched.clear();
  stripChainRaw.clear();
  std::vector<mmTrack> totalTime0;
  StripChainMatchingTime(totalTime0, chain_, strip_, leftSide, siTime, 0);
  StripChainMatchingTime(stripChainRaw, chain_, strip_, leftSide, siTime, 0);

  size_t numTimeBuckets = StripChainTime0TimeBuckets(totalTime0);
  // std::cout << "Time buckets: " << numTimeBuckets << std::endl;
  if(numTimeBuckets < 6) {
    StripChainMatchingBoxTime0(stripChainMatched, totalTime0);
    // std::cout << "BoxTime0: " << stripChainMatched.size() << std::endl;
  }
  // else if(chain_.size() > 12 && strip_.size() > 12) {
  //   StripChainMatchingTimeSlopeHough(stripChainMatched, chain_, strip_, leftSide, siTime, 10.);
  //   if(stripChainMatched.size() < 2) {
  //     stripChainMatched.clear();
  //     StripChainMatchingBoxTime0(stripChainMatched, totalTime0);
  //   }
  //   // std::cout << "SlopeHough0: " << stripChainMatched.size() << std::endl;
  // }
  else {
    StripChainMatchingTime(stripChainMatched, chain_, strip_, leftSide, siTime, 0);
    // std::cout << "MatchingTime: " << stripChainMatched.size() << std::endl;
  }
  // std::cout << chain_.size() << '\t' << strip_.size() << '\t' << stripChainMatched.size() << std::endl;
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
    if(strip[i].time < 4000 || strip[i].time > 7000) continue;
    stripReduced.push_back(strip[i]);
  }
  for(UInt_t i = 0; i < chain.size(); i++) {
    if(chain[i].time < 4000 || chain[i].time > 7000) continue;
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

void Spectra::DivideTargetThickness(TH1F *f) {
  Int_t i_size = f->GetSize();

  TAxis *x_axis = f->GetXaxis();
  for(Int_t i = 1; i <= i_size - 1; i++) {
    Float_t binLowEdge = x_axis->GetBinLowEdge(i);
    Float_t binUpEdge = x_axis->GetBinUpEdge(i);
    if(binLowEdge == 0.) binLowEdge += 0.001;
    binLowEdge *= (m1 + m2)/m2;
    binUpEdge *= (m1 + m2)/m2;
    Float_t binContent = f->GetBinContent(i);
    Float_t binError = f->GetBinError(i);
    Float_t delta_x = boronMethane->CalcRange(binUpEdge, binLowEdge);
    delta_x /= 10.;
    Float_t molarMassMethane = 0.01604;
    Float_t factor = 4.e-27*density*delta_x*TMath::Na()/molarMassMethane;
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
    Float_t binCenter = x_axis->GetBinCenter(i);
    Float_t binContent = f->GetBinContent(i);
    Float_t binError = f->GetBinError(i);

    Float_t simCS = reg3SA(binCenter);
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
    Float_t binCenter = x_axis->GetBinCenter(i);
    Float_t binContent = f->GetBinContent(i);
    Float_t binError = f->GetBinError(i);

    Float_t cmAngle = reg3CMAngle(binCenter);

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
    for(UInt_t i = 0; i < pointsXY.size(); i++) {
      Double_t d = pointsXY[i].x*cosj + pointsXY[i].y*sinj;
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