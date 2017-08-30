#define Spectra_cxx
#include "Spectra.h"

void Spectra::Loop() {

  using namespace std;

  InitChannelMap();

  bool individualMMHistograms = true;

  // Read gainFile.dat
  ifstream inGainFile("gainFile.dat");
  assert(inGainFile.is_open());

  int varI, varJ;
  double varScale;
  double scale[6][128];
  while(inGainFile >> varI >> varJ >> varScale) {
    scale[varI][varJ] = varScale;
  }

  inGainFile.close();

  if (fChain == 0) return;

  TFile* file = new TFile("spectra.root","recreate");

  Long64_t nentries = fChain->GetEntriesFast();

  ////////////////
  // Histograms //
  ////////////////

  // Histograms for the Forward Si Detectors
  for(int i=0; i<10 ; i++) {
    for(int j=0; j<4; j++) {
      TString name = Form("SiEForward_d%d_q%d_ch%d",i+1,j+1,siForwardChannel[i][j]);
      hSiEForward[i][j] = new TH1F(name,name,1000,0,4000);
      hSiEForward[i][j]->GetXaxis()->SetTitle("Channels");
      hSiEForward[i][j]->GetYaxis()->SetTitle("Counts");

      name = Form("SiTForward_d%d_q%d_ch%d",i+1,j+1, siForwardChannel[i][j]);
      hSiTForward[i][j] = new TH1F(name,name,1000,0,20000);
      hSiTForward[i][j]->GetXaxis()->SetTitle("Channels");
      hSiTForward[i][j]->GetYaxis()->SetTitle("Counts");
    }
  }

  // Histograms for the Beam Left Si Detectors
  for(int i=0; i<8 ; i++) {
    for(int j=0; j<4; j++) {
      TString name = Form("SiELeft_d%d_q%d_ch%d",i+1,j+1,siLeftChannel[i][j]);
      hSiELeft[i][j] = new TH1F(name,name,1000,0,4000);
      hSiELeft[i][j]->GetXaxis()->SetTitle("Channels");
      hSiELeft[i][j]->GetYaxis()->SetTitle("Counts");

      name = Form("SiTLeft_d%d_q%d_ch%d",i+1,j+1,siLeftChannel[i][j]);
      hSiTLeft[i][j] = new TH1F(name,name,1000,0,20000);
      hSiTLeft[i][j]->GetXaxis()->SetTitle("Channels");
      hSiTLeft[i][j]->GetYaxis()->SetTitle("Counts");
    }
  }

  // Histograms for the Forward CsI Detectors
  for(int i=0; i<10 ; i++) {
    TString name = Form("CsIEForward_d%d_ch%d",i+1,CsIForwardChannel[i]);
    hCsIEForward[i] = new TH1F(name,name,1000,0,4000);
    hCsIEForward[i]->GetXaxis()->SetTitle("Channels");
    hCsIEForward[i]->GetYaxis()->SetTitle("Counts");

    name = Form("CsITForward_d%d_ch%d",i+1,CsIForwardChannel[i]);
    hCsITForward[i] = new TH1F(name,name,1000,0,20000);
    hCsITForward[i]->GetXaxis()->SetTitle("Channels");
    hCsITForward[i]->GetYaxis()->SetTitle("Counts");
  }

  // Histograms for the Beam Left CsI Detectors
  for(int i=0; i<8 ; i++) {
    TString name = Form("CsIELeft_d%d_ch%d",i+1,CsILeftChannel[i]);
    hCsIELeft[i] = new TH1F(name,name,1000,0,4000);
    hCsIELeft[i]->GetXaxis()->SetTitle("Channels");
    hCsIELeft[i]->GetYaxis()->SetTitle("Counts");

    name = Form("CsITLeft_d%d_ch%d",i+1,CsILeftChannel[i]);
    hCsITLeft[i] = new TH1F(name,name,1000,0,20000);
    hCsITLeft[i]->GetXaxis()->SetTitle("Channels");
    hCsITLeft[i]->GetYaxis()->SetTitle("Counts");
  }

  TH2I* hMicroMegasCenter[nentries];
  TH1I* hMicroMegasStrip[nentries];
  TH1I* hMicroMegasStripTime[nentries];
  TH1I* hMicroMegasChain[nentries];
  TH1I* hMicroMegasChainTime[nentries];
  TH2I* hMicroMegas[nentries];
  TH2I* hMicroMegasDetect10[nentries];
  TH2I* hMicroMegasDetect12[nentries];

  if(individualMMHistograms) {
    for(Long64_t i=0; i<nentries; i++) {
      TString name = Form("Event_%lld_MM_Center", i);
      hMicroMegasCenter[i] = new TH2I(name, name, 82, -41, 41, 128, 0, 128);
      name = Form("Event_%lld_MM_Strip", i);
      hMicroMegasStrip[i] = new TH1I(name, name, 128, 0, 128);
      name = Form("Event_%lld_MM_Chain", i);
      hMicroMegasChain[i] = new TH1I(name, name, 128, 0, 128);
      name = Form("Event_%lld_MM_Strip_Time", i);
      hMicroMegasStripTime[i] = new TH1I(name, name, 128, 0, 128);
      name = Form("Event_%lld_MM_Chain_Time", i);
      hMicroMegasChainTime[i] = new TH1I(name, name, 128, 0, 128);
      name = Form("Event_%lld_MM", i);
      hMicroMegas[i] = new TH2I(name, name, 140, -70, 70, 128, 0, 128);
      name = Form("Event_%lld_MM_Detect10", i);
      hMicroMegasDetect10[i] = new TH2I(name, name, 140, -70, 70, 128, 0, 128);
      name = Form("Event_%lld_MM_Detect12", i);
      hMicroMegasDetect12[i] = new TH2I(name, name, 140, -70, 70, 128, 0, 128);
    }
  }

  TH2I* hMicroMegasCenterCumulative = new TH2I("MM_Center_Cumulative", "MM_Center_Cumulative", 82, -41, 41, 128, 0, 128);
  TH2I* hMicroMegasCenterEnergyCumulative = new TH2I("MM_Center_Energy_Cumulative", "MM_Center_Energy_Cumulative", 20, -10, 10, 128, 0, 128);
  TH2I* hMicroMegasCenterEnergyAverage = new TH2I("MM_Center_Energy_Average", "MM_Center_Energy_Average", 20, -10, 10, 128, 0, 128);
  TH2I* hMicroMegasCenterEnergyAverageScaled = new TH2I("MM_Center_Energy_Average_Scaled", "MM_Center_Energy_Average_Scaled", 20, -10, 10, 128, 0, 128);
  TH2I* hMicroMegasCenterTimeAverage = new TH2I("MM_Center_Time_Average", "MM_Center_Time_Average", 20, -10, 10, 128, 0, 128);

  TH1I* hIonizationChamber = new TH1I("IC", "IC", 5000, 0, 5000);

  TH1I* hSiDetect10 = new TH1I("hSiDetect10", "hSiDetect10", 5000, 0, 5000);
  TH1I* hSiDetect12 = new TH1I("hSiDetect12", "hSiDetect12", 5000, 0, 5000);

  // For average energy deposited
  std::vector<double> padCumulative[6][128];
  std::vector<double> padCumulativeScaled[6][128];
  std::vector<double> padCumulativeTime[6][128];

  Long64_t nbytes = 0, nb = 0;
  // for (Long64_t jentry=0; jentry<10; jentry++) {
  for (Long64_t jentry=0; jentry<nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;

    if(jentry!=0 && jentry%100==0) printf("Processed %lld events\n",jentry);

    // Find hit in Si detectors
    vector<int> SiDetect_;
    vector<int> SiQuad_;
    vector<double> SiEnergy_;
    vector<double> SiTime_;

    // Find hit in CsI Detectors
    vector<int> CsIDetect_;
    vector<double> CsIEnergy_;
    vector<double> CsITime_;

    // IC
    double ICE;
    double ICT;

    for(int i=0; i<mmMul; i++) {
      if(mmChan[i]==11 || mmChan[i]==22 || mmChan[i]==45 || mmChan[i]==56) continue; // Skip FPN Channels

      if(mmCobo[i]==0) {

      }

      else if(mmCobo[i]==1) {
        // Si Detectors
        if(mmAsad[i]==0) {
          // Forward Detectors
          if(mmAget[i]==0) {
            int detect = siForwardMap[mmChan[i]].first;
            int quad = siForwardMap[mmChan[i]].second;
            if(detect>0 && quad>0){
              SiDetect_.push_back(detect);
              SiQuad_.push_back(quad);
              SiEnergy_.push_back(mmEnergy[i]);
              SiTime_.push_back(mmTime[i]);
            }
          }
          // Beam Left Detectors
          else if(mmAget[i]==3) {
            int detect = siLeftMap[mmChan[i]].first;
            int quad = siLeftMap[mmChan[i]].second;
            if(detect>0 && quad>0){
              SiDetect_.push_back(detect+10);
              SiQuad_.push_back(quad);
              SiEnergy_.push_back(mmEnergy[i]);
              SiTime_.push_back(mmTime[i]);
            }
          }
        }

        // IC/CsI
        else if(mmAsad[i]==1) {

        }
      }
      // Forward Si Detectors
      if(mmCobo[i]==1 && mmAsad[i]==0 && mmAget[i]==0) {
        int detect = siForwardMap[mmChan[i]].first;
        int quad = siForwardMap[mmChan[i]].second;
        if(detect>0 && quad>0){
          SiDetect_.push_back(detect);
          SiQuad_.push_back(quad);
          SiEnergy_.push_back(mmEnergy[i]);
          SiTime_.push_back(mmTime[i]);
        }
      }

      // Beam Left Si Detectors
      else if(mmCobo[i]==1 && mmAsad[i]==0 && mmAget[i]==3) {
        int detect = siLeftMap[mmChan[i]].first;
        int quad = siLeftMap[mmChan[i]].second;
        if(detect>0 && quad>0){
          SiDetect_.push_back(detect+10);
          SiQuad_.push_back(quad);
          SiEnergy_.push_back(mmEnergy[i]);
          SiTime_.push_back(mmTime[i]);
        }
      }

      else

      // Ionization Chamber
      if(mmCobo[i]==1 && mmAsad[i]==1 && mmAget[i]==0 && mmChan[i]==2) {
        hIonizationChamber->Fill(mmEnergy[i]);
      }
    }

    /*
    for(int i=0; i<mmMul; i++) {

      if(mmChan[i]==11 || mmChan[i]==22 || mmChan[i]==45 || mmChan[i]==56) continue; // Skip FPN Channels

      // Micromegas
      if(mmCobo[i]==0) {
        if(mmEnergy[i]<50) continue;

        // Asad0 - Central Pads
        if(mmAsad[i]==0) {
          pair<int, int> pad;
          if(mmAget[i]==0) {
            pad = MM_Map_Asad0_Aget0[mmChan[i]];
          }
          else if(mmAget[i]==1) {
            pad = MM_Map_Asad0_Aget1[mmChan[i]];
          }
          else if(mmAget[i]==2) {
            pad = MM_Map_Asad0_Aget2[mmChan[i]];
          }
          else if(mmAget[i]==3) {
            pad = MM_Map_Asad0_Aget3[mmChan[i]];
          }

          hMicroMegasCenterCumulative->Fill(pad.first-3, pad.second);
          double energyBinBefore = hMicroMegasCenterEnergyCumulative->GetBinContent(pad.first+8, pad.second+1);
          if(energyBinBefore>0) hMicroMegasCenterEnergyCumulative->SetBinContent(pad.first+8, pad.second+1, energyBinBefore+mmEnergy[i]*scale[pad.first][pad.second]);

          if(individualMMHistograms) {
            hMicroMegasCenter[jentry]->SetBinContent(pad.first+39, pad.second+1, mmEnergy[i]);
            hMicroMegas[jentry]->SetBinContent(pad.first*2+65, pad.second+1, mmEnergy[i]);
            hMicroMegas[jentry]->SetBinContent(pad.first*2+65+1, pad.second+1, mmEnergy[i]);
          }

          if(mmEnergy[i]>0) padCumulative[pad.first][pad.second].push_back(mmEnergy[i]);
          if(mmEnergy[i]>0) padCumulativeScaled[pad.first][pad.second].push_back(mmEnergy[i]*scale[pad.first][pad.second]);
          if(mmTime[i]>0) padCumulativeTime[pad.first][pad.second].push_back(mmTime[i]);
        }

        // Asad1 - Central Pads
        else if(mmAsad[i]==1) {
          pair<int, int> pad;
          if(mmAget[i]==0) {
            pad = MM_Map_Asad1_Aget0[mmChan[i]];
          }
          else if(mmAget[i]==1) {
            pad = MM_Map_Asad1_Aget1[mmChan[i]];
          }
          else if(mmAget[i]==2) {
            pad = MM_Map_Asad1_Aget2[mmChan[i]];
          }
          else if(mmAget[i]==3) {
            pad = MM_Map_Asad1_Aget3[mmChan[i]];
          }

          hMicroMegasCenterCumulative->Fill(pad.first-3, pad.second);
          double energyBinBefore = hMicroMegasCenterEnergyCumulative->GetBinContent(pad.first+8, pad.second+1);
          if(energyBinBefore>0) hMicroMegasCenterEnergyCumulative->SetBinContent(pad.first+8, pad.second+1, energyBinBefore+mmEnergy[i]*scale[pad.first][pad.second]);

          if(individualMMHistograms) {
            hMicroMegasCenter[jentry]->SetBinContent(pad.first+39, pad.second+1, mmEnergy[i]);
            hMicroMegas[jentry]->SetBinContent(pad.first*2+65, pad.second+1, mmEnergy[i]);
            hMicroMegas[jentry]->SetBinContent(pad.first*2+65+1, pad.second+1, mmEnergy[i]);
          }

          if(mmEnergy[i]>0) padCumulative[pad.first][pad.second].push_back(mmEnergy[i]);
          if(mmEnergy[i]>0) padCumulativeScaled[pad.first][pad.second].push_back(mmEnergy[i]*scale[pad.first][pad.second]);
          if(mmTime[i]>0) padCumulativeTime[pad.first][pad.second].push_back(mmTime[i]);
        }

        // Asad2 - Aget0&1 - Strips and Chains. Aget2&3 - Central Pads
        else if(mmAsad[i]==2) {
          if(mmAget[i]==0) {
            int bin = MM_Map_Asad2_Aget0[mmChan[i]];
            if(Aget_Map[mmChan[i]]<34) {
              hMicroMegasStrip[jentry]->SetBinContent(bin, mmEnergy[i]);
              hMicroMegasStripTime[jentry]->SetBinContent(bin, mmTime[i]);
            }
            else{
              hMicroMegasChain[jentry]->SetBinContent(bin, mmEnergy[i]);
              hMicroMegasChainTime[jentry]->SetBinContent(bin, mmTime[i]);
            }
          }
          else if(mmAget[i]==1) {
            int bin = MM_Map_Asad2_Aget1[mmChan[i]];
            if(Aget_Map[mmChan[i]]<34) {
              hMicroMegasStrip[jentry]->SetBinContent(bin, mmEnergy[i]);
              hMicroMegasStripTime[jentry]->SetBinContent(bin, mmTime[i]);
            }
            else{
              hMicroMegasChain[jentry]->SetBinContent(bin, mmEnergy[i]);
              hMicroMegasChainTime[jentry]->SetBinContent(bin, mmTime[i]);
            }
          }
          else if(mmAget[i]==2) {
            pair<int, int> pad = MM_Map_Asad2_Aget2[mmChan[i]];

            hMicroMegasCenterCumulative->Fill(pad.first-3, pad.second);
            double energyBinBefore = hMicroMegasCenterEnergyCumulative->GetBinContent(pad.first+8, pad.second+1);
            if(energyBinBefore>0) hMicroMegasCenterEnergyCumulative->SetBinContent(pad.first+8, pad.second+1, energyBinBefore+mmEnergy[i]*scale[pad.first][pad.second]);

            if(individualMMHistograms) {
              hMicroMegasCenter[jentry]->SetBinContent(pad.first+39, pad.second+1, mmEnergy[i]);
              hMicroMegas[jentry]->SetBinContent(pad.first*2+65, pad.second+1, mmEnergy[i]);
              hMicroMegas[jentry]->SetBinContent(pad.first*2+65+1, pad.second+1, mmEnergy[i]);
            }

            if(mmEnergy[i]>0) padCumulative[pad.first][pad.second].push_back(mmEnergy[i]);
            if(mmEnergy[i]>0) padCumulativeScaled[pad.first][pad.second].push_back(mmEnergy[i]*scale[pad.first][pad.second]);
            if(mmTime[i]>0) padCumulativeTime[pad.first][pad.second].push_back(mmTime[i]);
          }
          else if(mmAget[i]==3) {
            pair<int, int> pad = MM_Map_Asad2_Aget3[mmChan[i]];

            hMicroMegasCenterCumulative->Fill(pad.first-3, pad.second);
            double energyBinBefore = hMicroMegasCenterEnergyCumulative->GetBinContent(pad.first+8, pad.second+1);
            if(energyBinBefore>0) hMicroMegasCenterEnergyCumulative->SetBinContent(pad.first+8, pad.second+1, energyBinBefore+mmEnergy[i]*scale[pad.first][pad.second]);

            if(individualMMHistograms) {
              hMicroMegasCenter[jentry]->SetBinContent(pad.first+39, pad.second+1, mmEnergy[i]);
              hMicroMegas[jentry]->SetBinContent(pad.first*2+65, pad.second+1, mmEnergy[i]);
              hMicroMegas[jentry]->SetBinContent(pad.first*2+65+1, pad.second+1, mmEnergy[i]);
            }

            if(mmEnergy[i]>0) padCumulative[pad.first][pad.second].push_back(mmEnergy[i]);
            if(mmEnergy[i]>0) padCumulativeScaled[pad.first][pad.second].push_back(mmEnergy[i]*scale[pad.first][pad.second]);
            if(mmTime[i]>0) padCumulativeTime[pad.first][pad.second].push_back(mmTime[i]);
          }
        }

        // Asad3 - Aget0&1 - Strips and Chains. Aget2&3 - Central Pads
        else if(mmAsad[i]==3) {
          if(mmAget[i]==0) {
            int bin = MM_Map_Asad3_Aget0[mmChan[i]];
            if(Aget_Map[mmChan[i]]<34) {
              hMicroMegasStrip[jentry]->SetBinContent(bin+64, mmEnergy[i]);
              hMicroMegasStripTime[jentry]->SetBinContent(bin+64, mmTime[i]);
            }
            else{
              hMicroMegasChain[jentry]->SetBinContent(bin+64, mmEnergy[i]);
              hMicroMegasChainTime[jentry]->SetBinContent(bin+64, mmTime[i]);
            }
          }
          else if(mmAget[i]==1) {
            int bin = MM_Map_Asad3_Aget1[mmChan[i]];

            if(Aget_Map[mmChan[i]]<34) {

              hMicroMegasStrip[jentry]->SetBinContent(bin+64, mmEnergy[i]);
              hMicroMegasStripTime[jentry]->SetBinContent(bin+64, mmTime[i]);
            }
            else{
              hMicroMegasChain[jentry]->SetBinContent(bin+64, mmEnergy[i]);
              hMicroMegasChainTime[jentry]->SetBinContent(bin+64, mmTime[i]);
            }
          }
          else if(mmAget[i]==2) {
            pair<int, int> pad = MM_Map_Asad3_Aget2[mmChan[i]];

            hMicroMegasCenterCumulative->Fill(pad.first-3, pad.second);
            double energyBinBefore = hMicroMegasCenterEnergyCumulative->GetBinContent(pad.first+8, pad.second+1);
            if(energyBinBefore>0) hMicroMegasCenterEnergyCumulative->SetBinContent(pad.first+8, pad.second+1, energyBinBefore+mmEnergy[i]*scale[pad.first][pad.second]);

            if(individualMMHistograms) {
              hMicroMegasCenter[jentry]->SetBinContent(pad.first+39, pad.second+1, mmEnergy[i]);
              hMicroMegas[jentry]->SetBinContent(pad.first*2+65, pad.second+1, mmEnergy[i]);
              hMicroMegas[jentry]->SetBinContent(pad.first*2+65+1, pad.second+1, mmEnergy[i]);
            }

            if(mmEnergy[i]>0) padCumulative[pad.first][pad.second].push_back(mmEnergy[i]);
            if(mmEnergy[i]>0) padCumulativeScaled[pad.first][pad.second].push_back(mmEnergy[i]*scale[pad.first][pad.second]);
            if(mmTime[i]>0) padCumulativeTime[pad.first][pad.second].push_back(mmTime[i]);
          }
          else if(mmAget[i]==3) {
            pair<int, int> pad = MM_Map_Asad3_Aget3[mmChan[i]];

            hMicroMegasCenterCumulative->Fill(pad.first-3, pad.second);
            double energyBinBefore = hMicroMegasCenterEnergyCumulative->GetBinContent(pad.first+8, pad.second+1);
            if(energyBinBefore>0) hMicroMegasCenterEnergyCumulative->SetBinContent(pad.first+8, pad.second+1, energyBinBefore+mmEnergy[i]*scale[pad.first][pad.second]);

            if(individualMMHistograms) {
              hMicroMegasCenter[jentry]->SetBinContent(pad.first+39, pad.second+1, mmEnergy[i]);
              hMicroMegas[jentry]->SetBinContent(pad.first*2+65, pad.second+1, mmEnergy[i]);
              hMicroMegas[jentry]->SetBinContent(pad.first*2+65+1, pad.second+1, mmEnergy[i]);
            }

            if(mmEnergy[i]>0) padCumulative[pad.first][pad.second].push_back(mmEnergy[i]);
            if(mmEnergy[i]>0) padCumulativeScaled[pad.first][pad.second].push_back(mmEnergy[i]*scale[pad.first][pad.second]);
            if(mmTime[i]>0) padCumulativeTime[pad.first][pad.second].push_back(mmTime[i]);
          }
        }
      }

      // Forward Si Detectors
      if(mmCobo[i]==1 && mmAsad[i]==0 && mmAget[i]==0) {
        int detect = siForwardMap[mmChan[i]].first;
        int quad = siForwardMap[mmChan[i]].second;
        if(detect>0 && quad>0){
          hSiEForward[siForwardMap[mmChan[i]].first-1][siForwardMap[mmChan[i]].second-1]->Fill(mmEnergy[i]);
          hSiTForward[siForwardMap[mmChan[i]].first-1][siForwardMap[mmChan[i]].second-1]->Fill(mmTime[i]);
        }
      }

      // Beam Left Si Detectors
      if(mmCobo[i]==1 && mmAsad[i]==0 && mmAget[i]==3) {
        int detect = siForwardMap[mmChan[i]].first;
        int quad = siForwardMap[mmChan[i]].second;
        if(detect>0 && quad>0){
          hSiELeft[siForwardMap[mmChan[i]].first-1][siForwardMap[mmChan[i]].second-1]->Fill(mmEnergy[i]);
          hSiTLeft[siForwardMap[mmChan[i]].first-1][siForwardMap[mmChan[i]].second-1]->Fill(mmTime[i]);
        }
      }

      // CsI Detectors
      if(mmCobo[i]==1 && mmAsad[i]==1 && mmAget[i]==3) {
        int detect = CsIForwardMap[mmChan[i]];
        if(detect>0) {
          // Forward CsI Detectors
          if(mmChan[i]<43) {
            hCsIEForward[CsIForwardMap[mmChan[i]]-1]->Fill(mmEnergy[i]);
            hCsITForward[CsIForwardMap[mmChan[i]]-1]->Fill(mmTime[i]);
          }
          // Beam Left CsI Detectors
          else if (mmChan[i]>43){
            hCsIELeft[CsILeftMap[mmChan[i]]-1]->Fill(mmEnergy[i]);
            hCsITLeft[CsILeftMap[mmChan[i]]-1]->Fill(mmTime[i]);
          }
        }
      }

      // Ionization Chamber
      if(mmCobo[i]==1 && mmAsad[i]==1 && mmAget[i]==0 && mmChan[i]==2) {
        hIonizationChamber->Fill(mmEnergy[i]);
      }
    }
    */

    // Calculate average energy
    double padAverage[6][128];
    double padAverageScaled[6][128];
    double padAverageTime[6][128];
    for(int i=0; i<6; i++) {
      for(int j=0; j<128; j++) {
        if(i==4 && j==0) continue;
        padAverage[i][j] = 0.;
        padAverageScaled[i][j] = 0.;
        padAverageTime[i][j] = 0.;
        for(size_t k=0; k<padCumulative[i][j].size(); k++) {
          padAverage[i][j] += padCumulative[i][j][k];
        }
        for(size_t k=0; k<padCumulativeScaled[i][j].size(); k++) {
          padAverageScaled[i][j] += padCumulativeScaled[i][j][k];
        }
        for(size_t k=0; k<padCumulativeTime[i][j].size(); k++) {
          padAverageTime[i][j] += padCumulativeTime[i][j][k];
        }
        padAverage[i][j] = (padCumulative[i][j].size()>0) ? padAverage[i][j]/padCumulative[i][j].size() : 0.;
        padAverageScaled[i][j] = (padCumulativeScaled[i][j].size()>0) ? padAverageScaled[i][j]/padCumulativeScaled[i][j].size() : 0.;
        padAverageTime[i][j] = (padCumulativeTime[i][j].size()>0) ? padAverageTime[i][j]/padCumulativeTime[i][j].size() : 0.;
        if(padAverage[i][j]>0) hMicroMegasCenterEnergyAverage->SetBinContent(i+8, j+1, padAverage[i][j]);
        if(padAverageScaled[i][j]>0) hMicroMegasCenterEnergyAverageScaled->SetBinContent(i+8, j+1, padAverageScaled[i][j]);
        if(padAverageTime[i][j]>0) hMicroMegasCenterTimeAverage->SetBinContent(i+8, j+1, padAverageTime[i][j]);
      }
    }

    // // Gain match all pads
    // FILE* gainFile = fopen("gainFile.dat", "w");
    // double padBest = padAverage[2][64];
    // double rowEnergySlope = 0.01338555;
    // double rowEnergyIntercept = 7.7997787;
    // double padBestEnergy = 8.85;
    // double gainScale[6][128];
    // for(int i=0; i<6; i++) {
    //   for(int j=0; j<128; j++) {
    //     double energyRow = j*rowEnergySlope + rowEnergyIntercept;
    //     gainScale[i][j] = padBest*energyRow/(padAverage[i][j]*padBestEnergy);
    //     fprintf(gainFile, "%d %d %f\n", i, j, gainScale[i][j]);
    //   }
    // }
    // fflush(gainFile);
    // fclose(gainFile);

    int i_size_strip = hMicroMegasStrip[jentry]->GetSize();
    int lowBoundStrip = 1000;
    int highBoundStrip = 0;
    std::vector<pair<int, double> > stripEnergy;
    for(int i=1; i<i_size_strip; i++) {
      double binContent = hMicroMegasStrip[jentry]->GetBinContent(i);
      if(binContent<200) continue;
      if(i<lowBoundStrip) lowBoundStrip = i;
      if(i>highBoundStrip) highBoundStrip = i;
      stripEnergy.push_back(make_pair(i, binContent));
    }

    int i_size_chain = hMicroMegasChain[jentry]->GetSize();
    int lowBoundChain = 1000;
    int highBoundChain = 0;
    std::vector<pair<int, double> > chainEnergy;
    for(int i=1; i<i_size_chain; i++) {
      double binContent = hMicroMegasChain[jentry]->GetBinContent(i);
      if(binContent<200) continue;
      if(i<lowBoundChain) lowBoundChain = i;
      if(i>highBoundChain) highBoundChain = i;
      chainEnergy.push_back(make_pair(i, binContent));
    }

    for(size_t i=0; i<stripEnergy.size(); i++) {
      int stripBin = (stripEnergy[i].first<64) ? stripEnergy[i].first*2 : (stripEnergy[i].first-64)*2;
      for(size_t j=0; j<chainEnergy.size(); j++) {
        int chainBin = (chainEnergy[j].first<64) ? 64-chainEnergy[j].first : chainEnergy[j].first+12;
        hMicroMegas[jentry]->SetBinContent(chainBin+1, stripBin-1, stripEnergy[i].second);
        hMicroMegas[jentry]->SetBinContent(chainBin+1, stripBin, stripEnergy[i].second);
        if(highSiEDet==10) {
          hMicroMegasDetect10[jentry]->SetBinContent(chainBin+1, stripBin-1, stripEnergy[i].second);
          hMicroMegasDetect10[jentry]->SetBinContent(chainBin+1, stripBin, stripEnergy[i].second);
        }
        if(highSiEDet==12) {
          hMicroMegasDetect12[jentry]->SetBinContent(chainBin+1, stripBin-1, stripEnergy[i].second);
          hMicroMegasDetect12[jentry]->SetBinContent(chainBin+1, stripBin, stripEnergy[i].second);
        }
      }
    }
  }

  // Try to fix hMicroMegasCenterEnergyCumulative
  int i_size_cumulative_x = hMicroMegasCenterCumulative->GetXaxis()->GetNbins();
  int i_size_cumulative_y = hMicroMegasCenterCumulative->GetYaxis()->GetNbins();
  for(int i=0; i<i_size_cumulative_x; i++) {
    for(int j=0; j<i_size_cumulative_y; j++) {
      double binContent = hMicroMegasCenterCumulative->GetBinContent(i, j);
      if(binContent<0) hMicroMegasCenterCumulative->SetBinContent(i, j, 0.);
    }
  }

  for(Long64_t i=0; i<nentries; i++) {
    if(hMicroMegasCenter[i]->GetEntries()>1) {
      hMicroMegasCenter[i]->Write();
      hMicroMegas[i]->Write();
      if(hMicroMegasStrip[i]->GetEntries()>0) hMicroMegasStrip[i]->Write();
      if(hMicroMegasStripTime[i]->GetEntries()>0) hMicroMegasStripTime[i]->Write();
      if(hMicroMegasChain[i]->GetEntries()>0) hMicroMegasChain[i]->Write();
      if(hMicroMegasChainTime[i]->GetEntries()>0) hMicroMegasChainTime[i]->Write();
    }
  }

  for(int i=0; i<10 ; i++) {
    for(int j=0; j<4; j++) {
      hSiEForward[i][j]->Write();
      hSiTForward[i][j]->Write();
    }
    hCsIEForward[i]->Write();
    hCsITForward[i]->Write();
  }

  for(int i=0; i<6 ; i++) {
    for(int j=0; j<4; j++) {
      hSiELeft[i][j]->Write();
      hSiTLeft[i][j]->Write();
    }
    hCsIELeft[i]->Write();
    hCsITLeft[i]->Write();
  }

  hIonizationChamber->Write();

  hMicroMegasCenterCumulative->Write();
  hMicroMegasCenterEnergyCumulative->Write();
  hMicroMegasCenterEnergyAverage->Write();
  hMicroMegasCenterEnergyAverageScaled->Write();
  hMicroMegasCenterTimeAverage->Write();

  file->Close();
}
