#define Spectra_cxx
#include "Spectra.h"

void Spectra::Loop() {

  if (fChain == 0) return;

  InitChannelMap();

  Initialize();

  bool individualMMHistograms = false;

  // Read si calibrations
  std::pair<double, double> siEForwardCalibration[10][4] = {std::make_pair(0., 0.)};
  std::ifstream inSiCalFile("siCalibration.dat");
  assert(inSiCalFile.is_open());
  int var1, var2, var3;
  double slope, intercept;
  while(inSiCalFile >> var1 >> var2 >> var3 >> slope >> intercept) {
    siEForwardCalibration[var1-1][var2-1] = std::make_pair(slope, intercept);
  }
  inSiCalFile.close();

  // Read gainFile.dat
  std::ifstream inGainFile("gainFile.dat");
  assert(inGainFile.is_open());
  int varI, varJ;
  double varScale;
  double scale[6][128];
  while(inGainFile >> varI >> varJ >> varScale) {
    scale[varI][varJ] = varScale;
  }
  inGainFile.close();

  Long64_t nentries = fChain->GetEntriesFast();

  /////////////////
  // Set up cuts //
  /////////////////
  TFile* cutFile = TFile::Open("cuts.root");
  TCutG* det5_dEE_noPunchthroughCut = (TCutG*)cutFile->Get("det5_dEECut_noPunchthrough");
  cutFile->Close();

  ////////////////
  // Histograms //
  ////////////////

  // Histograms for the Forward Si Detectors
  for(int i=0; i<10; i++) {
    TString name = Form("SiEForward_d%d", i+1);
    hSiEForwardTotal[i] = new TH1F(name, name, 1000, 0, 4000);
    hSiEForwardTotal[i]->GetXaxis()->SetTitle("Channels"); hSiEForwardTotal[i]->GetXaxis()->CenterTitle();
    hSiEForwardTotal[i]->GetYaxis()->SetTitle("Counts"); hSiEForwardTotal[i]->GetXaxis()->CenterTitle(); hSiEForwardTotal[i]->GetYaxis()->SetTitleOffset(1.4);

    name = Form("SiTForward_d%d", i+1);
    hSiTForwardTotal[i] = new TH1F(name, name, 1000, 0, 20000);
    hSiTForwardTotal[i]->GetXaxis()->SetTitle("Timing [ns]"); hSiTForwardTotal[i]->GetXaxis()->CenterTitle();
    hSiTForwardTotal[i]->GetYaxis()->SetTitle("Counts"); hSiTForwardTotal[i]->GetXaxis()->CenterTitle(); hSiTForwardTotal[i]->GetYaxis()->SetTitleOffset(1.4);
    for(int j=0; j<4; j++) {
      TString name = Form("SiEForward_d%d_q%d_ch%d", i+1, j+1, siForwardChannel[i][j]);
      hSiEForward[i][j] = new TH1F(name, name, 1000, 0, 4000);
      hSiEForward[i][j]->GetXaxis()->SetTitle("Channels"); hSiEForward[i][j]->GetXaxis()->CenterTitle();
      hSiEForward[i][j]->GetYaxis()->SetTitle("Counts"); hSiEForward[i][j]->GetYaxis()->CenterTitle(); hSiEForward[i][j]->GetYaxis()->SetTitleOffset(1.4);

      name = Form("SiEForwardCal_d%d_q%d_ch%d", i+1, j+1, siForwardChannel[i][j]);
      hSiEForwardCal[i][j] = new TH1F(name, name, 1000, 0, 20000);
      hSiEForwardCal[i][j]->GetXaxis()->SetTitle("Energy [keV]"); hSiEForwardCal[i][j]->GetXaxis()->CenterTitle();
      hSiEForwardCal[i][j]->GetYaxis()->SetTitle("Counts"); hSiEForwardCal[i][j]->GetYaxis()->CenterTitle(); hSiEForwardCal[i][j]->GetYaxis()->SetTitleOffset(1.4);

      name = Form("SiTForward_d%d_q%d_ch%d", i+1, j+1, siForwardChannel[i][j]);
      hSiTForward[i][j] = new TH1F(name, name, 1000, 0, 20000);
      hSiTForward[i][j]->GetXaxis()->SetTitle("Channels"); hSiTForward[i][j]->GetXaxis()->CenterTitle();
      hSiTForward[i][j]->GetYaxis()->SetTitle("Counts"); hSiTForward[i][j]->GetYaxis()->CenterTitle(); hSiTForward[i][j]->GetYaxis()->SetTitleOffset(1.4);
    }
  }

  // Histograms for the Beam Left Si Detectors
  for(int i=0; i<8; i++) {
    TString name = Form("SiELeft_d%d", i+1);
    hSiELeftTotal[i] = new TH1F(name, name, 1000, 0, 4000);
    hSiELeftTotal[i]->GetXaxis()->SetTitle("Channels"); hSiELeftTotal[i]->GetXaxis()->CenterTitle();
    hSiELeftTotal[i]->GetYaxis()->SetTitle("Counts"); hSiELeftTotal[i]->GetYaxis()->CenterTitle(); hSiELeftTotal[i]->GetYaxis()->SetTitleOffset(1.4);
    for(int j=0; j<4; j++) {
      TString name = Form("SiELeft_d%d_q%d_ch%d", i+1, j+1, siLeftChannel[i][j]);
      hSiELeft[i][j] = new TH1F(name, name, 1000, 0, 4000);
      hSiELeft[i][j]->GetXaxis()->SetTitle("Channels"); hSiELeft[i][j]->GetXaxis()->CenterTitle();
      hSiELeft[i][j]->GetYaxis()->SetTitle("Counts"); hSiELeft[i][j]->GetYaxis()->CenterTitle(); hSiELeft[i][j]->GetYaxis()->SetTitleOffset(1.4);

      name = Form("SiELeftCal_d%d_q%d_ch%d", i+1, j+1, siLeftChannel[i][j]);
      hSiELeftCal[i][j] = new TH1F(name, name, 1000, 0, 4000);
      hSiELeftCal[i][j]->GetXaxis()->SetTitle("Channels"); hSiELeftCal[i][j]->GetXaxis()->CenterTitle();
      hSiELeftCal[i][j]->GetYaxis()->SetTitle("Counts"); hSiELeftCal[i][j]->GetYaxis()->CenterTitle(); hSiELeftCal[i][j]->GetYaxis()->SetTitleOffset(1.4);

      name = Form("SiTLeft_d%d_q%d_ch%d", i+1, j+1, siLeftChannel[i][j]);
      hSiTLeft[i][j] = new TH1F(name, name, 1000, 0, 20000);
      hSiTLeft[i][j]->GetXaxis()->SetTitle("Channels"); hSiTLeft[i][j]->GetXaxis()->CenterTitle();
      hSiTLeft[i][j]->GetYaxis()->SetTitle("Counts"); hSiTLeft[i][j]->GetYaxis()->CenterTitle(); hSiTLeft[i][j]->GetYaxis()->SetTitleOffset(1.4);
    }
  }

  // Histograms for the Forward CsI Detectors
  for(int i=0; i<10; i++) {
    TString name = Form("CsIEForward_d%d_ch%d", i+1, csiForwardChannel[i]);
    hCsIEForward[i] = new TH1F(name, name, 1000, 0, 4000);
    hCsIEForward[i]->GetXaxis()->SetTitle("Channels"); hCsIEForward[i]->GetXaxis()->CenterTitle();
    hCsIEForward[i]->GetYaxis()->SetTitle("Counts"); hCsIEForward[i]->GetYaxis()->CenterTitle(); hCsIEForward[i]->GetYaxis()->SetTitleOffset(1.4);

    name = Form("CsITForward_d%d_ch%d", i+1, csiForwardChannel[i]);
    hCsITForward[i] = new TH1F(name, name, 1000, 0, 20000);
    hCsITForward[i]->GetXaxis()->SetTitle("Channels"); hCsITForward[i]->GetXaxis()->CenterTitle();
    hCsITForward[i]->GetYaxis()->SetTitle("Counts"); hCsITForward[i]->GetYaxis()->CenterTitle(); hCsITForward[i]->GetYaxis()->SetTitleOffset(1.4);
  }

  // Histograms for the Beam Left CsI Detectors
  for(int i=0; i<8; i++) {
    TString name = Form("CsIELeft_d%d_ch%d", i+1, csiLeftChannel[i]);
    hCsIELeft[i] = new TH1F(name,name,1000,0,4000);
    hCsIELeft[i]->GetXaxis()->SetTitle("Channels"); hCsIELeft[i]->GetXaxis()->CenterTitle();
    hCsIELeft[i]->GetYaxis()->SetTitle("Counts"); hCsIELeft[i]->GetYaxis()->CenterTitle(); hCsIELeft[i]->GetYaxis()->SetTitleOffset(1.4);

    name = Form("CsITLeft_d%d_ch%d", i+1, csiLeftChannel[i]);
    hCsITLeft[i] = new TH1F(name,name,1000,0,20000);
    hCsITLeft[i]->GetXaxis()->SetTitle("Channels"); hCsITLeft[i]->GetXaxis()->CenterTitle();
    hCsITLeft[i]->GetYaxis()->SetTitle("Counts"); hCsITLeft[i]->GetYaxis()->CenterTitle(); hCsITLeft[i]->GetYaxis()->SetTitleOffset(1.4);
  }

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
  hVertexE->GetXaxis()->SetTitle("Energy [keV]"); hVertexE->GetXaxis()->SetTitleSize(0.05); hVertexE->GetXaxis()->CenterTitle();
  hVertexE->GetYaxis()->SetTitle("dE [Channels]"); hVertexE->GetYaxis()->SetTitleSize(0.05); hVertexE->GetYaxis()->CenterTitle(); hVertexE->GetYaxis()->SetTitleOffset(1.);
  TH2F* hVertexEPunchthrough = new TH2F("vertexEPunchthrough", "vertexEPunchthrough", 128, 0, 128, 500, 0, 4000);
  hVertexEPunchthrough->GetXaxis()->SetTitle("Energy [keV]"); hVertexEPunchthrough->GetXaxis()->SetTitleSize(0.05); hVertexEPunchthrough->GetXaxis()->CenterTitle();
  hVertexEPunchthrough->GetYaxis()->SetTitle("dE [Channels]"); hVertexEPunchthrough->GetYaxis()->SetTitleSize(0.05); hVertexEPunchthrough->GetYaxis()->CenterTitle(); hVertexEPunchthrough->GetYaxis()->SetTitleOffset(1.);
  TH2F* hVertexECal = new TH2F("vertexECal", "vertexECal", 128, 0, 128, 500, 0, 20000);
  hVertexECal->GetXaxis()->SetTitle("Energy [keV]"); hVertexECal->GetXaxis()->SetTitleSize(0.05); hVertexECal->GetXaxis()->CenterTitle();
  hVertexECal->GetYaxis()->SetTitle("dE [Channels]"); hVertexECal->GetYaxis()->SetTitleSize(0.05); hVertexECal->GetYaxis()->CenterTitle(); hVertexECal->GetYaxis()->SetTitleOffset(1.);
  TH2F* hVertexEPunchthroughCal = new TH2F("vertexEPunchthroughCal", "vertexEPunchthroughCal", 128, 0, 128, 500, 0, 20000);
  hVertexEPunchthroughCal->GetXaxis()->SetTitle("Energy [keV]"); hVertexEPunchthroughCal->GetXaxis()->SetTitleSize(0.05); hVertexEPunchthroughCal->GetXaxis()->CenterTitle();
  hVertexEPunchthroughCal->GetYaxis()->SetTitle("dE [Channels]"); hVertexEPunchthroughCal->GetYaxis()->SetTitleSize(0.05); hVertexEPunchthroughCal->GetYaxis()->CenterTitle(); hVertexEPunchthroughCal->GetYaxis()->SetTitleOffset(1.);

  hCSDet5 = new TH1F("CSDet5", "CSDDet5", 20, 0, 5); hCSDet5->Sumw2();
  hCSDet5->GetXaxis()->SetTitle("Center of Mass Energy [MeV]"); hCSDet5->GetXaxis()->SetTitleSize(0.05); hCSDet5->GetXaxis()->CenterTitle();
  hCSDet5->GetYaxis()->SetTitle("Cross Section [b/sr]"); hCSDet5->GetYaxis()->SetTitleSize(0.05); hCSDet5->GetYaxis()->CenterTitle(); hCSDet5->GetYaxis()->SetTitleOffset(1.);

  hCSDet5Counts = new TH1F("CSDet5Counts", "CSDDet5Counts", 20, 0, 5); hCSDet5Counts->Sumw2();
  hCSDet5Counts->GetXaxis()->SetTitle("Center of Mass Energy [MeV]"); hCSDet5Counts->GetXaxis()->SetTitleSize(0.05); hCSDet5Counts->GetXaxis()->CenterTitle();
  hCSDet5Counts->GetYaxis()->SetTitle("Counts"); hCSDet5Counts->GetYaxis()->SetTitleSize(0.05); hCSDet5Counts->GetYaxis()->CenterTitle(); hCSDet5Counts->GetYaxis()->SetTitleOffset(1.);

  // Make output file
  TFile* file = new TFile("spectra.root","recreate");

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

    if(jentry!=0 && jentry%5000==0) printf("Processed %lld events\n",jentry);

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
            if(detect>0 && quad>0){
              siDetect siHit = {detect, quad, mmEnergy[i], mmTime[i]};
              siDetect_.push_back(siHit);
            }
          }
          // Beam Left Detectors
          else if(mmAget[i]==3) {
            int detect = siLeftMap[mmChan[i]].first;
            int quad = siLeftMap[mmChan[i]].second;
            if(detect>0 && quad>0){
              siDetect siHit = {detect+10, quad, mmEnergy[i], mmTime[i]};
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

    if(icE<1250) continue;
    if(icT<5500 || icT>6500) continue;

    int siDet = siDetect_[0].detect;
    int siQuad = siDetect_[0].quad;
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
      hSiEForward[siDet-1][siQuad-1]->Fill(siEnergy);
      hSiEForwardCal[siDet-1][siQuad-1]->Fill(siEnergyCal);
    }
    else {
      hSiELeftTotal[siDet-11]->Fill(siEnergy);
      hSiELeft[siDet-11][siQuad-1]->Fill(siEnergy);
      hSiELeftCal[siDet-1][siQuad-1]->Fill(siEnergyCal);
    }

    // Sum up rows in central pads
    float centerPadEnergy[128];
    for(size_t i=0; i<128; i++) {
      centerPadEnergy[i] = 0.;
    }

    // Fill MM Center Histograms
    float centerPadEnergydE = 0.;
    for(size_t i=0; i<mmCenterMatched_.size(); i++) {
      hMicroMegasCenterCumulative->Fill(mmCenterMatched_[i].column-3, mmCenterMatched_[i].row);
      hCenterPadTime->Fill(mmCenterMatched_[i].time);

      if(mmCenterMatched_[i].time>4000 && mmCenterMatched_[i].time<7000) centerPadEnergy[mmCenterMatched_[i].row] += mmCenterMatched_[i].energy;

      if(mmCenterMatched_[i].row>117 && mmCenterMatched_[i].row<124 && mmCenterMatched_[i].column!=0 && mmCenterMatched_[i].column!=5) {
        if(mmCenterMatched_[i].time>4000 && mmCenterMatched_[i].time<7000) centerPadEnergydE += mmCenterMatched_[i].energy;
      }
    }

    // Make dE plots
    if(siDet==5 && !punchthrough && centerPadEnergydE>100) {
      hdEECenterPad->Fill(siEnergy, centerPadEnergydE);
      hdEECenterPadCal->Fill(siEnergyCal, centerPadEnergydE);
    }
    if(siDet==5 && centerPadEnergydE>100) {
      hdEECenterPadPunchthrough->Fill(siEnergy, centerPadEnergydE);
      hdEECenterPadPunchthroughCal->Fill(siEnergyCal, centerPadEnergydE);
    }

    // Make sure the reaction happened over the micromegas and not before
    if(centerPadEnergy[0]>0) {
      float maxEnergy = 0;
      int maxRow = 0;
      for(size_t i=0; i<110; i++) {
        if(centerPadEnergy[i]>maxEnergy) {
          maxEnergy = centerPadEnergy[i];
          maxRow = i;
        }
      }
      if(siDet==5 && centerPadEnergydE<20000 && !punchthrough) {
        hVertexE->Fill(maxRow, siEnergy);
        hVertexECal->Fill(maxRow, siEnergyCal);
      }
      if(siDet==5 && centerPadEnergydE<20000) {
        hVertexEPunchthrough->Fill(maxRow, siEnergy);
        hVertexEPunchthroughCal->Fill(maxRow, siEnergyCal);
      }
    }

    // Simple Cross Section calculation for detector 5
    if(det5_dEE_noPunchthroughCut->IsInside(siEnergyCal, centerPadEnergydE)) {
      Double_t massFactor = m1*m2/((m1+m2)*(m1+m2));
      Double_t factor = 4.*massFactor; // missing cos^2 term but setting to 1 for now
      Double_t b8Energy = siEnergyCal/factor;
      Double_t cmEnergy = b8Energy*m2/(m1+m2);
      hCSDet5->Fill(cmEnergy/1000.);
      hCSDet5Counts->Fill(cmEnergy/1000.);
    }
  }

  // Simple Cross Section calculation for detector 5
  SimpleCrossSection(hCSDet5);

  // // Try to fix hMicroMegasCenterEnergyCumulative
  // int i_size_cumulative_x = hMicroMegasCenterCumulative->GetXaxis()->GetNbins();
  // int i_size_cumulative_y = hMicroMegasCenterCumulative->GetYaxis()->GetNbins();
  // for(int i=0; i<i_size_cumulative_x; i++) {
  //   for(int j=0; j<i_size_cumulative_y; j++) {
  //     double binContent = hMicroMegasCenterCumulative->GetBinContent(i, j);
  //     if(binContent<0) hMicroMegasCenterCumulative->SetBinContent(i, j, 0.);
  //   }
  // }

  // for(Long64_t i=0; i<nentries; i++) {
  //   if(hMicroMegasCenter[i]->GetEntries()>1) {
  //     hMicroMegasCenter[i]->Write();
  //     hMicroMegas[i]->Write();
  //     if(hMicroMegasStrip[i]->GetEntries()>0) hMicroMegasStrip[i]->Write();
  //     if(hMicroMegasStripTime[i]->GetEntries()>0) hMicroMegasStripTime[i]->Write();
  //     if(hMicroMegasChain[i]->GetEntries()>0) hMicroMegasChain[i]->Write();
  //     if(hMicroMegasChainTime[i]->GetEntries()>0) hMicroMegasChainTime[i]->Write();
  //   }
  // }

  for(int i=0; i<10 ; i++) {
    // hSiEForwardTotal[i]->Write();
    hSiTForwardTotal[i]->Write();
    // for(int j=0; j<4; j++) {
      // hSiEForward[i][j]->Write();
      // hSiEForwardCal[i][j]->Write();
      // hSiTForward[i][j]->Write();
    // }
    // hCsIEForward[i]->Write();
    // hCsITForward[i]->Write();
  }

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
  hdEECenterPadCal->Write();
  hdEECenterPadPunchthroughCal->Write();
  // hCenterPadTime->Write();

  // hVertexE->Write();
  // hVertexEPunchthrough->Write();
  hVertexECal->Write();
  hVertexEPunchthroughCal->Write();

  hCSDet5->Write();
  hCSDet5Counts->Write();

  file->Close();
}

void Spectra::Initialize() {
  m1 = 8.; // AMU of projectile
  m2 = 1.; // AMU of target
  beamEnergy = 56.; // In MeV, after havar window
  density = 0.00038175; // in g/cm3, from LISE++ (Methane at 435 torr)
  numberB8 = 174809089.;

  distanceHavarToSilicon = 544.07; // Distance from Havar to Forward Silicon in mm

  // Initialize EnergyLoss
  boronMethane = new EnergyLoss("b8_methane.dat");
  protonMethane = new EnergyLoss("proton_methane.dat");
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
