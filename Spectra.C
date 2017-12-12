#define Spectra_cxx
#include "Spectra.h"

void Spectra::Loop() {

  if (fChain == 0) return;

  InitChannelMap();
  InitHistograms();
  InitVariables();

  InitSiEForwardCalibration();

  InitCentralPadGainMatch();
  InitAverageBeamEnergy();

  bool individualMMHistograms = false;

  Long64_t nentries = fChain->GetEntriesFast();

  /////////////////
  // Set up cuts //
  /////////////////
  TFile* cutFile = TFile::Open("cuts.root");
  TCutG* det5_dEE_noPunchthroughCut = (TCutG*)cutFile->Get("det5_dEECut_noPunchthrough");
  TCutG* det56_dEECut = (TCutG*)cutFile->Get("det56_dEECut");
  cutFile->Close();

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

  // Make output file
  TFile* file = new TFile("spectra.root","recreate");

  // For average energy deposited
  std::vector<double> padCumulative[6][128];
  std::vector<double> padCumulativeScaled[6][128];
  std::vector<double> rowCumulativeScaled[128];
  std::vector<double> padCumulativeTime[6][128];

  printf("Starting Main Loop\n");

  Long64_t nbytes = 0, nb = 0;
  // for (Long64_t jentry=0; jentry<10; jentry++) {
  for (Long64_t jentry=0; jentry<nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;

    // if(jentry!=0 && jentry%200==0) printf("Processed %lld events\n",jentry);

    // Find hit in Si detectors
    std::vector<siDetect> siDetect_;

    // Find hit in CsI Detectors
    std::vector<csiDetect> csiDetect_;

    // Central Pads in Micromegas
    std::vector<mmCenter> mmCenter_;
    std::vector<mmCenter> mmCenterMatched_;

    // Beam Left Strips in Micromegas
    std::vector<mmStrip> mmLeftStrip_;

    // Beam Left Chains in Micromegas
    std::vector<mmChain> mmLeftChain_;

    // Beam Right Strips in Micromegas
    std::vector<mmStrip> mmRightStrip_;

    // Beam Right Chains in Micromegas
    std::vector<mmChain> mmRightChain_;

    // IC
    int icE = 0.;
    int icT = 0.;

    for(int i=0; i<mmMul; i++) {
      if(mmChan[i]==11 || mmChan[i]==22 || mmChan[i]==45 || mmChan[i]==56) continue; // Skip FPN Channels

      // Micromegas
      if(mmCobo[i]==0) {
        // Asad0 - All Center Pads
        if(mmAsad[i]==0) {
          // Aget0
          if(mmAget[i]==0) {
            std::pair<int,int> pad = MM_Map_Asad0_Aget0[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, mmEnergy[i], mmTime[i]};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, mmEnergy[i]*scale[pad.first][pad.second], mmTime[i]};
            mmCenterMatched_.push_back(mmHitMatched);
          }
          // Aget1
          else if(mmAget[i]==1) {
            std::pair<int,int> pad = MM_Map_Asad0_Aget1[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, mmEnergy[i], mmTime[i]};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, mmEnergy[i]*scale[pad.first][pad.second], mmTime[i]};
            mmCenterMatched_.push_back(mmHitMatched);
          }
          // Aget2
          else if(mmAget[i]==2) {
            std::pair<int,int> pad = MM_Map_Asad0_Aget2[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, mmEnergy[i], mmTime[i]};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, mmEnergy[i]*scale[pad.first][pad.second], mmTime[i]};
            mmCenterMatched_.push_back(mmHitMatched);
          }
          // Aget3
          else if(mmAget[i]==3) {
            std::pair<int,int> pad = MM_Map_Asad0_Aget3[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, mmEnergy[i], mmTime[i]};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, mmEnergy[i]*scale[pad.first][pad.second], mmTime[i]};
            mmCenterMatched_.push_back(mmHitMatched);
          }
        }
        // Asad1
        else if(mmAsad[i]==1) {
          // Aget0
          if(mmAget[i]==0) {
            std::pair<int,int> pad = MM_Map_Asad1_Aget0[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, mmEnergy[i], mmTime[i]};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, mmEnergy[i]*scale[pad.first][pad.second], mmTime[i]};
            mmCenterMatched_.push_back(mmHitMatched);
          }
          // Aget1
          else if(mmAget[i]==1) {
            std::pair<int,int> pad = MM_Map_Asad1_Aget1[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, mmEnergy[i], mmTime[i]};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, mmEnergy[i]*scale[pad.first][pad.second], mmTime[i]};
            mmCenterMatched_.push_back(mmHitMatched);
          }
          // Aget2
          else if(mmAget[i]==2) {
            std::pair<int,int> pad = MM_Map_Asad1_Aget2[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, mmEnergy[i], mmTime[i]};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, mmEnergy[i]*scale[pad.first][pad.second], mmTime[i]};
            mmCenterMatched_.push_back(mmHitMatched);
          }
          // Aget3
          else if(mmAget[i]==3) {
            std::pair<int,int> pad = MM_Map_Asad1_Aget3[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, mmEnergy[i], mmTime[i]};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, mmEnergy[i]*scale[pad.first][pad.second], mmTime[i]};
            mmCenterMatched_.push_back(mmHitMatched);
          }
        }
        // Asad2
        else if(mmAsad[i]==2) {
          // Aget0 - Strips and Chains
          if(mmAget[i]==0) {
            int bin = MM_Map_Asad2_Aget0[mmChan[i]];
          }
          // Aget1 - Strips and Chains
          else if(mmAget[i]==1) {
            int bin = MM_Map_Asad2_Aget1[mmChan[i]];
          }
          // Aget2 - Outside Central Pads
          else if(mmAget[i]==2) {
            std::pair<int,int> pad = MM_Map_Asad2_Aget2[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, mmEnergy[i], mmTime[i]};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, mmEnergy[i]*scale[pad.first][pad.second], mmTime[i]};
            mmCenterMatched_.push_back(mmHitMatched);
          }
          // Aget3 - Outside Central Pads
          else if(mmAget[i]==3) {
            std::pair<int,int> pad = MM_Map_Asad2_Aget3[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, mmEnergy[i], mmTime[i]};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, mmEnergy[i]*scale[pad.first][pad.second], mmTime[i]};
            mmCenterMatched_.push_back(mmHitMatched);
          }
        }
        // Asad3
        else if(mmAsad[i]==3) {
          // Aget0
          if(mmAget[i]==0) {
            int bin = MM_Map_Asad3_Aget0[mmChan[i]];
          }
          // Aget1
          else if(mmAget[i]==1) {
            int bin = MM_Map_Asad3_Aget1[mmChan[i]];
          }
          // Aget2
          else if(mmAget[i]==2) {
            std::pair<int,int> pad = MM_Map_Asad3_Aget2[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, mmEnergy[i], mmTime[i]};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, mmEnergy[i]*scale[pad.first][pad.second], mmTime[i]};
            mmCenterMatched_.push_back(mmHitMatched);
          }
          // Aget3
          else if(mmAget[i]==3) {
            std::pair<int,int> pad = MM_Map_Asad3_Aget3[mmChan[i]];
            mmCenter mmHit = {pad.first, pad.second, mmEnergy[i], mmTime[i]};
            mmCenter_.push_back(mmHit);
            mmCenter mmHitMatched = {pad.first, pad.second, mmEnergy[i]*scale[pad.first][pad.second], mmTime[i]};
            mmCenterMatched_.push_back(mmHitMatched);
          }
        }
      }

      else if(mmCobo[i]==1) {
        // Si Detectors
        if(mmAsad[i]==0) {
          // Forward Detectors
          if(mmAget[i]==0) {
            int detect = siForwardMap[mmChan[i]].first;
            int quad = siForwardMap[mmChan[i]].second;
            int channel = mmChan[i];
            if(detect>0 && quad>0){
              siDetect siHit = {detect, quad, channel, mmEnergy[i], mmTime[i]};
              siDetect_.push_back(siHit);
            }
          }
          // Beam Left Detectors
          else if(mmAget[i]==3) {
            int detect = siLeftMap[mmChan[i]].first;
            int quad = siLeftMap[mmChan[i]].second;
            int channel = mmChan[i];
            if(detect>0 && quad>0){
              siDetect siHit = {detect+10, quad, channel, mmEnergy[i], mmTime[i]};
              siDetect_.push_back(siHit);
            }
          }
        }

        // IC/CsI
        else if(mmAsad[i]==1) {
          // Ionization Chamber
          if(mmAget[i]==0 && mmChan[i]==10) {
            icE = mmEnergy[i];
            icT = mmTime[i];
          }
          // CsI Detectors
          else if(mmAget[i]==3) {
            if(mmChan[i]<43) {
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
    if(siDetect_.size()!=1) continue;

    hIonizationChamberE->Fill(icE);
    hIonizationChamberT->Fill(icT);

    // if(icE<1250) continue;
    if(icE < 1500 || icE > 1900) continue;
    // if(icT<5500 || icT>6500) continue;

    int siDet = siDetect_[0].detect;
    int siQuad = siDetect_[0].quad;
    int siChannel = siDetect_[0].channel;
    double siEnergy = siDetect_[0].energy;
    double siEnergyCal = 0.;
    if(siDet<11) {
      siEnergyCal = siDetect_[0].energy*siEForwardCalibration[siDet-1][siQuad-1].first + siEForwardCalibration[siDet-1][siQuad-1].second;
    }
    double siTime = siDetect_[0].time;

    if(siDet<11) hSiTForwardTotal[siDet-1]->Fill(siTime);

    // Find if CsI behind Si fired
    bool punchthrough = false;
    for(size_t i=0; i<csiDetect_.size(); i++) {
      if(csiDetect_[i].detect==siDet) punchthrough = true;
    }

    // Fill Silicon Histograms
    if(siDet<11) {
      hSiEForwardTotal[siDet-1]->Fill(siEnergy);
      hSiEForwardTotalCal[siDet-1]->Fill(siEnergyCal);
      hSiEForward[siDet-1][siQuad-1]->Fill(siEnergy);
      hSiEForwardCal[siDet-1][siQuad-1]->Fill(siEnergyCal);
    }
    else {
      hSiELeftTotal[siDet-11]->Fill(siEnergy);
      hSiELeftTotalCal[siDet-11]->Fill(siEnergyCal);
      hSiELeft[siDet-11][siQuad-1]->Fill(siEnergy);
      hSiELeftCal[siDet-1][siQuad-1]->Fill(siEnergyCal);
    }

    // Sum up rows in central pads
    float centerPadEnergy[128];
    float centerdEPadEnergy[6][4];
    for(size_t i=0; i<128; i++) {
      centerPadEnergy[i] = 0.;
    }

    // Fill MM Center Histograms
    float centerPadEnergydE = 0.;
    for(size_t i=0; i<mmCenterMatched_.size(); i++) {
      hMicroMegasCenterCumulative->Fill(mmCenterMatched_[i].column-3, mmCenterMatched_[i].row);
      hCenterPadTime->Fill(mmCenterMatched_[i].time);

      padCumulativeScaled[mmCenterMatched_[i].column][mmCenterMatched_[i].row].push_back(mmCenterMatched_[i].energy);
      padCumulativeTime[mmCenterMatched_[i].column][mmCenterMatched_[i].row].push_back(mmCenterMatched_[i].time);

      if(mmCenterMatched_[i].time>4000 && mmCenterMatched_[i].time<7000 && mmCenterMatched_[i].column!=0 && mmCenterMatched_[i].column!=5) centerPadEnergy[mmCenterMatched_[i].row] += mmCenterMatched_[i].energy;

      if(mmCenterMatched_[i].row>117 && mmCenterMatched_[i].row<124 && mmCenterMatched_[i].column!=0 && mmCenterMatched_[i].column!=5) {
        if(mmCenterMatched_[i].time>4000 && mmCenterMatched_[i].time<7000) centerPadEnergydE += mmCenterMatched_[i].energy;
      }
    }

    double padNumber[112];
    double centerPadEnergySubtracted[112];
    for(unsigned int i=0; i<112; i++) {
      padNumber[i] = i;
      if(centerPadEnergy[i] > 0) {
        centerPadEnergySubtracted[i] = centerPadEnergy[i] - averageBeamEnergy[i];
      }
      else {
        centerPadEnergySubtracted[i] = 0.;
      }
    }

    // Find point where centerPadEnergySubtracted goes over 200 (or else find the maximum)
    int vertexRow = -1;
    if(siDet == 5 && centerPadEnergydE<20000) {
      if(centerPadEnergy[0] == 0 || centerPadEnergy[1] == 0 || centerPadEnergy[2] == 0) continue;
      bool found = false;
      for(unsigned int i=0; i<112; i++) {
        if(centerPadEnergySubtracted[i] > 200 && !found) {
          vertexRow = i;
          found = true;
        }
      }

      if(vertexRow == -1) { // did not find a point above 200, not find the maximum
        double maxEnergy = -1000.;
        for(unsigned int i=0; i<112; i++) {
          if(centerPadEnergySubtracted[i] > maxEnergy) {
            maxEnergy = centerPadEnergySubtracted[i];
            vertexRow = i;
          }
        }
      }
      hVertexECal->Fill(vertexRow, siEnergyCal);
    }

    if(vertexRow == -1) continue;

    // Reduce MM Center to one entry per row
    std::map<int, double> centralPadPosition;
    std::map<int, double> centralPadTotalEnergy;
    std::map<int, double> centralPadTime;
    std::map<int, int> centralPadTotal;
    for(unsigned int i=0; i<128; i++) {
      centralPadTotal[i] = 0;
    }
    for(size_t i=0; i<mmCenterMatched_.size(); i++) {
      if(mmCenterMatched_[i].time>4000 && mmCenterMatched_[i].time<7000) {
        if(centralPadTotal[mmCenterMatched_[i].row] == 0) {
          centralPadPosition[mmCenterMatched_[i].row] = (mmCenterMatched_[i].column*3.0-7.5)*mmCenterMatched_[i].energy;
          centralPadTotalEnergy[mmCenterMatched_[i].row] = mmCenterMatched_[i].energy;
          centralPadTime[mmCenterMatched_[i].row] = mmCenterMatched_[i].time;
        }
        else {
          centralPadPosition[mmCenterMatched_[i].row] += (mmCenterMatched_[i].column*3.0-7.5)*mmCenterMatched_[i].energy;
          centralPadTotalEnergy[mmCenterMatched_[i].row] += mmCenterMatched_[i].energy;
          centralPadTime[mmCenterMatched_[i].row] += mmCenterMatched_[i].time;
        }
        centralPadTotal[mmCenterMatched_[i].row]++;
      }
    }

    std::vector<mmCenterTrack> mmCenterBeam;
    std::vector<mmCenterTrack> mmCenterProton;
    std::map<int, double>::iterator it;
    for(it = centralPadPosition.begin(); it != centralPadPosition.end(); it++) {
      int row = it->first;
      // std::cout << it->first << " : " << it->second << std::endl;
      if(row < 112) { // beam
        mmCenterTrack hit = {0., 0, 0., 0., 0., 0};
        hit.position = centralPadPosition[row]/centralPadTotalEnergy[row];
        hit.row = row;
        hit.time = centralPadTime[row]/centralPadTotal[row];
        hit.energy = centralPadTotalEnergy[row];
        hit.height = (zeroTime - hit.time)*driftVelocity;
        hit.total = centralPadTotal[row];
        mmCenterBeam.push_back(hit);
      }
      else if((row == 116) || (row > 117 && row<124)) { // proton
        mmCenterTrack hit = {0., 0, 0., 0., 0., 0};
        hit.position = centralPadPosition[row]/centralPadTotalEnergy[row];
        hit.row = row;
        hit.time = centralPadTime[row]/centralPadTotal[row];
        hit.energy = centralPadTotalEnergy[row];
        hit.height = (zeroTime - hit.time)*driftVelocity;
        hit.total = centralPadTotal[row];
        mmCenterProton.push_back(hit);
      }
    }

    // Find point on beamline corresponding to the vertex row
    for(size_t i=0; i<mmCenterBeam.size(); i++) {
      if(mmCenterBeam[i].row == vertexRow) {
        mmCenterProton.push_back(mmCenterBeam[i]);
      }
    }

    // // Plot event by event for vertex
    // if(det5_dEE_noPunchthroughCut->IsInside(siEnergyCal, centerPadEnergydE)) {
    //   TString name = Form("Vertex_Event_%lld", jentry);
    //   TGraph *gr = new TGraph(112, padNumber, centerPadEnergySubtracted);
    //   gr->SetName(name);
    //   gr->Write();
    //   delete gr;
    // }

    // for(unsigned int i=0; i<128; i++) {
    //   if(centerPadEnergy[i] > 5) {
    //     rowCumulativeScaled[i].push_back(centerPadEnergy[i]);
    //   }
    // }

    // // Make dE plots
    // if((siDet==5 || siDet==6) && !punchthrough && centerPadEnergydE>100) {
    //   hdEECenterPad->Fill(siEnergy, centerPadEnergydE);
    //   hdEECenterPadCal->Fill(siEnergyCal, centerPadEnergydE);
    // }
    // if((siDet==5 || siDet==6) && centerPadEnergydE>100) {
    //   hdEECenterPadPunchthrough->Fill(siEnergy, centerPadEnergydE);
    //   hdEECenterPadPunchthroughCal->Fill(siEnergyCal, centerPadEnergydE);
    // }

    // // Make sure the reaction happened over the micromegas and not before
    // float maxEnergy = 0;
    // int maxRow = 0;
    // if(centerPadEnergy[0]>0 && centerPadEnergy[1]>0 && centerPadEnergy[2]>0) {
    //   for(size_t i=0; i<110; i++) {
    //     if(centerPadEnergy[i]>maxEnergy) {
    //       maxEnergy = centerPadEnergy[i];
    //       maxRow = i;
    //     }
    //   }
    //   // Skip events with no maximum in micromegas
    //   if(maxRow==0) continue;

    //   if(siDet==5 && centerPadEnergydE<20000 && !punchthrough) {
    //     hVertexE->Fill(maxRow, siEnergy);
    //     // hVertexECal->Fill(maxRow, siEnergyCal);
    //   }
    //   if(siDet==5 && centerPadEnergydE<20000) {
    //     hVertexEPunchthrough->Fill(maxRow, siEnergy);
    //     hVertexEPunchthroughCal->Fill(maxRow, siEnergyCal);
    //   }
    // }

    // // Skip events with no maximum in micromegas
    // if(maxRow==0) continue;

    // // For events scattering over the micromegas, make 3d histogram (TH3F)
    // if(det5_dEE_noPunchthroughCut->IsInside(siEnergyCal, centerPadEnergydE)) {
    //   TString name = Form("Track_Event_%lld", jentry);
    //   TH3F *h_3d_track= new TH3F(name, name, 10, -5, 5, 128, 0, 128, 400, 2000, 10000);
    //   for(size_t i=0; i<mmCenterMatched_.size(); i++) {
    //     double xPosition = mmCenterMatched_
    //     h_3d_track->Fill(mmCenterMatched_[i].column-3, mmCenterMatched_[i].row, mmCenterMatched_[i].time);
    //   }
    //   h_3d_track->Write();
    //   delete h_3d_track;
    // }

    // For events scattering over the micromegas, make 3d graph (TGraph2D)
    if(det5_dEE_noPunchthroughCut->IsInside(siEnergyCal, centerPadEnergydE)) {
      // Loop over beam particles
      int N = 0;
      TString name = Form("Track_Event_%lld_Beam", jentry);
      TGraph2D *h_3d_track_beam= new TGraph2D();
      h_3d_track_beam->SetMarkerStyle(20);
      h_3d_track_beam->SetMarkerColor(2);
      for(size_t i=0; i<mmCenterBeam.size(); i++) {
        mmCenterTrack hit = mmCenterBeam[i];
        // h_3d_track_beam->SetPoint(i, hit.position, hit.row, hit.time);
        h_3d_track_beam->SetPoint(i, hit.position, hit.row, hit.height);
      }
      h_3d_track_beam->SetName(name);
      h_3d_track_beam->Write();
      delete h_3d_track_beam;

      // Loop over proton particles
      name = Form("Track_Event_%lld_Proton", jentry);
      TGraph2D *h_3d_track_proton= new TGraph2D();
      h_3d_track_proton->SetMarkerStyle(8);
      h_3d_track_proton->SetMarkerColor(4);
      for(size_t i=0; i<mmCenterProton.size(); i++) {
        mmCenterTrack hit = mmCenterProton[i];
        // h_3d_track_proton->SetPoint(i, hit.position, hit.row, hit.time);
        h_3d_track_proton->SetPoint(i, hit.position, hit.row, hit.height);
      }

      // Fit proton track
      double p0[4] = {10, 20, 1, 2};
      ROOT::Fit::Fitter fitter;
      SumDistance2 sdist(h_3d_track_proton);
      ROOT::Math::Functor fcn(sdist, 4);
      double pStart[4] = {1, 1, 1, 1};
      fitter.SetFCN(fcn, pStart);
      for(int i = 0; i < 4; i++) {
        fitter.Config().ParSettings(i).SetStepSize(0.01);
      }
      bool ok = fitter.FitFCN();
      if(!ok) {
        Error("line3Dfit", "Line3D Fit failed");
        continue;
      }
      const ROOT::Fit::FitResult &result = fitter.Result();
      const double *parFit = result.GetParams();
      h_3d_track_proton->SetName(name);
      h_3d_track_proton->Write();
      delete h_3d_track_proton;

      // Draw fitted line
      int n = 1000;
      double t0 = -100;
      double dt = 200;
      TGraph2D *l_proton = new TGraph2D(n);
      for(int i = 0; i < n; i++) {
        double t = t0 + dt*i/n;
        double x, y, z;
        line(t, parFit, x, y, z);
        std::cout << x << " " << y << " " << z << std::endl;
        l_proton->SetPoint(i, x, y, z);
      }
      l_proton->SetLineColor(kGreen);
      name = Form("Fit_Event_%lld_Proton", jentry);
      l_proton->SetName(name);
      l_proton->Write();
      delete l_proton;

      // TGraph2D boundaries
      TGraph2D *h_3d_track_bound = new TGraph2D();
      h_3d_track_bound->SetMarkerStyle(8);
      h_3d_track_bound->SetMarkerColor(5);
      h_3d_track_bound->SetPoint(0, -30, 0, -100);
      h_3d_track_bound->SetPoint(1, 30, 128, 100);
      h_3d_track_bound->SetName(Form("Track_Event_%lld_Bound", jentry));
      h_3d_track_bound->Write();
      delete h_3d_track_bound;

      printf("Entry: %lld Si Det: %d Si Quad: %d Si Chan: %d Vertex Row: %d\n", jentry, siDet, siQuad, siChannel, vertexRow);
    }

    // // Simple Cross Section calculation for detector 5
    // if(det5_dEE_noPunchthroughCut->IsInside(siEnergyCal, centerPadEnergydE)) {
    //   Double_t massFactor = m1*m2/((m1+m2)*(m1+m2));
    //   Double_t factor = 4.*massFactor; // missing cos^2 term but setting to 1 for now
    //   Double_t b8Energy = siEnergyCal/factor;
    //   Double_t cmEnergy = b8Energy*m2/(m1+m2);
    //   hCSDet5->Fill(cmEnergy/1000.);
    //   hCSDet5Counts->Fill(cmEnergy/1000.);
    //   if(maxRow>80 && siEnergyCal<1800) {
    //     printf("Event: %lld Row: %d Energy: %f SiEnergy: %f\n", jentry, maxRow, maxEnergy, siEnergyCal);
    //     TString name = Form("vertexEvent_%lld", jentry);
    //     TH1F* hVertexEvent = new TH1F(name, name, 115, 0, 115);
    //     for(size_t row=0; row<110; row++) {
    //       hVertexEvent->SetBinContent(row+1, centerPadEnergy[row]);
    //     }
    //     hVertexEvent->Write();
    //     delete hVertexEvent;
    //   }
    // }

    //  ** End of event by event analysis ** //
  }

  // For beam average
  // FILE* averageBeamEnergyFile = fopen("averageBeamEnergy.out", "w");
  // FILE* averageBeamTimeFile = fopen("averageBeamTime.out", "w");
  // Loop over rows
  // for(unsigned int j=0; j<128; j++) {
    // double averageEnergy = 0.;
    // double averageTime = 0.;
    // Loop over rows
    // for(size_t k=0; k<rowCumulativeScaled[j].size(); k++) {
      // averageEnergy += rowCumulativeScaled[j][k];
      // averageTime += padCumulativeTime[i][j][k];
    // }
    // if(rowCumulativeScaled[j].size() == 0) {
      // averageEnergy = 0.;
      // averageTime = 0.;
    // }
    // else {
      // averageEnergy /= rowCumulativeScaled[j].size();
      // averageTime /= padCumulativeScaled[i][j].size();
    // }
    // fprintf(averageBeamEnergyFile, "%d %f\n", j, averageEnergy);
    // fprintf(averageBeamTimeFile, "%d %f\n", j, averageTime);
  // }
  // fflush(averageBeamEnergyFile);
  // fflush(averageBeamTimeFile);
  // fclose(averageBeamEnergyFile);
  // fclose(averageBeamTimeFile);

  // // Simple Cross Section calculation for detector 5
  // SimpleCrossSection(hCSDet5);

  // FILE* cmCSFile;
  // cmCSFile = fopen("testCS.out", "w");
  // int i_size = hCSDet5->GetSize();
  // TAxis *xaxis = hCSDet5->GetXaxis();
  // for(int i=1; i<i_size; i++) {
  //   float binCenter = xaxis->GetBinCenter(i);
  //   float binContent = hCSDet5->GetBinContent(i);
  //   float binError = hCSDet5->GetBinError(i);
  //   fprintf(cmCSFile, "%f %f %f\n", binCenter, binContent, binError);
  // }
  // fflush(cmCSFile);
  // fclose(cmCSFile);

  // for(int i=0; i<10 ; i++) {
    // hSiEForwardTotal[i]->Write();
    // hSiTForwardTotal[i]->Write();
    // for(int j=0; j<4; j++) {
      // hSiEForward[i][j]->Write();
      // hSiEForwardCal[i][j]->Write();
      // hSiTForward[i][j]->Write();
    // }
    // hCsIEForward[i]->Write();
    // hCsITForward[i]->Write();
  // }

  // for(int i=0; i<6 ; i++) {
    // hSiELeftTotal[i]->Write();
    // for(int j=0; j<4; j++) {
      // hSiELeft[i][j]->Write();
      // hSiTLeft[i][j]->Write();
    // }
    // hCsIELeft[i]->Write();
    // hCsITLeft[i]->Write();
  // }

  // hIonizationChamberE->Write();
  // hIonizationChamberT->Write();

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

void Spectra::SimpleCrossSection(TH1F *f) {
  SimpleSolidAngleDet5(f);
  DivideTargetThickness(f);
  f->Scale(1./numberB8);
}

void Spectra::SimpleSolidAngleDet5(TH1F *f) {
  Int_t i_size = f->GetSize();

  TAxis *xaxis = f->GetXaxis();
  for(Int_t i=1; i<i_size; i++) {
    Double_t binCenter = xaxis->GetBinCenter(i);
    Double_t binContent = f->GetBinContent(i);
    Double_t binError = f->GetBinError(i);
    Double_t solidAngle = CalcSimpleSolidAngleDet5(binCenter);
    if(solidAngle==0.) {
      f->SetBinContent(i, 0.);
      f->SetBinError(i, 0.);
    }
    else {
      f->SetBinContent(i, binContent/solidAngle);
      f->SetBinError(i, binError/solidAngle);
    }

  }
}

Double_t Spectra::CalcSimpleSolidAngleDet5(Double_t cmEnergy) {
  Double_t b8Energy = cmEnergy*(m1+m2)/m2;
  Double_t dist = boronMethane->CalcRange(beamEnergy, b8Energy);

  Double_t distToSi = distanceHavarToSilicon - dist;

  TH1F* h = new TH1F("h", "h", 360, 0, 180);
  TH1F* h2 = new TH1F("h2", "h2", 360, 0, 180);

  Double_t angle;
  Double_t sum = 0;
  for(Double_t dx=0; dx<25.0; dx+=0.1) {
    for(Double_t dy=26.8; dy<76.8; dy+=0.1) {
      Double_t dr = sqrt(dx*dx+dy*dy+distToSi*distToSi);
      angle = acos(distToSi/dr)/3.14159*180.;
      Double_t protonEnergy = 4.*m1/(m1+m2)*distToSi*distToSi/(dr*dr)*cmEnergy;
      Double_t range = protonMethane->CalcRange(protonEnergy, 0.);
      if(range>=dr) {
        sum += 0.1*0.1*distToSi/(dr*dr*dr);
        h->Fill(angle);
        h2->Fill(180.-2.*angle);
      }
    }
  }
  sum *= 2.;

  if(sum==0) {
    delete h;
    delete h2;
    return 0.;
  }
  else {
    Double_t aFactor = 4.*cos(h->GetMean()*3.14159/180.);
    Double_t changeBinContent = sum*aFactor;
    delete h;
    delete h2;
    return changeBinContent;
  }
}

void Spectra::DivideTargetThickness(TH1F *f) {
  Int_t i_size = f->GetSize();

  TAxis *xaxis = f->GetXaxis();
  for(Int_t i=1; i<i_size; i++) {
    Double_t binLowEdge = xaxis->GetBinLowEdge(i);
    Double_t binUpEdge = xaxis->GetBinUpEdge(i);
    if(binLowEdge==0.) binLowEdge += 0.001;
    binLowEdge *= (m1+m2)/m2; // From C.M. to Lab frame
    binUpEdge *= (m1+m2)/m2; // From C.M. to Lab frame
    Double_t binContent = f->GetBinContent(i);
    Double_t binError = f->GetBinError(i);
    Double_t delta_x = boronMethane->CalcRange(binUpEdge, binLowEdge);
    delta_x /= 10.0;
    Double_t molarMassMethane = 0.01604;
    Double_t factor = 4.e-27*density*delta_x*TMath::Na()/molarMassMethane;
    binContent /= factor;
    binError /= factor;
    f->SetBinContent(i, binContent);
    f->SetBinError(i, binError);
  }
}
