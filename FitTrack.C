#include "FitTrack.h"

// Calculate distance line-point
Double_t distance2(Double_t x, Double_t y, Double_t z, const Double_t *p) {
  // distance line point is D = | (xp-x0) cross ux |
  // where ux is direction of line and x0 is a point in the line (like t = 0)
  ROOT::Math::XYZVector xp(x, y, z);
  ROOT::Math::XYZVector x0(p[0], 0., p[2]);
  ROOT::Math::XYZVector x1(p[0] + p[1], 1., p[2] + p[3]);
  ROOT::Math::XYZVector u = (x1 - x0).Unit();
  Double_t d2 = ((xp - x0).Cross(u)).Mag2();
  return d2;
}

FitTrack::FitTrack(): minDist(100000), fitLine(true), fitHough(true) {};

FitTrack::FitTrack(std::vector<mmTrack> track): minDist(100000), fitLine(true), fitHough(true) {
  tracks.push_back(track);
}

FitTrack::FitTrack(std::vector<std::vector<mmTrack> > listOfTracks): minDist(100000), fitLine(true), fitHough(true) {
  tracks = listOfTracks;
}

void FitTrack::AddTrack(std::vector<mmTrack> track) {
  tracks.push_back(track);
}

void FitTrack::DisableLineFit() {
  fitLine = false;
}

void FitTrack::DisableHoughFit() {
  fitHough = false;
}

Double_t FitTrack::MakeBestFit() {
  for(UInt_t i = 0; i < tracks.size(); i++) {
    std::vector<Double_t> parTrack;
    Double_t minDistTrack = Fit(tracks[i], parTrack);
    if(minDistTrack < minDist) {
      minDist = minDistTrack;
      std::copy(std::begin(parTrack), std::end(parTrack), std::begin(parFit));
    }
  }
  return minDist;
}

Double_t FitTrack::MakeFit() {
  parFit.clear();
  minDist = Fit(tracks[0], parFit);
  return minDist;
}

Double_t FitTrack::Fit(std::vector<mmTrack> track, std::vector<Double_t> &par) {
  Double_t minDistLine = 1.e15;
  std::vector<Double_t> parLine;
  Double_t minDistHough = 1.e15;
  std::vector<Double_t> parHough;

  if(fitLine) {
    minDistLine = FitLine(track, parLine);
  }
  if(fitHough) {
    minDistHough = FitHough(track, parHough);
  }
  if(!fitLine && !fitHough) {
    std::cout << "Error: Not using fitting method line or Hough. Why the hell are you fitting then?" << std::endl;
  }

  if(isnan(minDistLine)) minDistLine = 100000.;
  if(isnan(minDistHough)) minDistHough = 100000.;

  if(minDistLine < minDistHough) {
    par = parLine;
    return minDistLine;
  }
  else {
    par = parHough;
    return minDistHough;
  }
}

Double_t FitTrack::FitLine(std::vector<mmTrack> track, std::vector<Double_t> &par) {
  par.clear();
  TGraph2D *h_3d_track = new TGraph2D();
  for(UInt_t i = 0; i < track.size(); i++) {
    mmTrack hit = track[i];
    h_3d_track->SetPoint(i, hit.xPosition, hit.yPosition, hit.height);
  }

  // Fit track
  ROOT::Fit::Fitter fitter;
  SumDistance2 sdist(h_3d_track);
  ROOT::Math::Functor fcn(sdist, 4);
  Double_t pStart[4] = {1, 1, 1, 1};
  fitter.SetFCN(fcn, pStart);
  for(UInt_t i = 0; i < 4; i++) {
    fitter.Config().ParSettings(i).SetStepSize(0.01);
  }
  Bool_t ok = fitter.FitFCN();
  if(!ok) {
    par.push_back(0.);
    par.push_back(0.);
    par.push_back(0.);
    par.push_back(0.);
    return 100000.;
  }
  const ROOT::Fit::FitResult &result = fitter.Result();
  const Double_t *parFit = result.GetParams();
  delete h_3d_track;

  par.push_back(parFit[0]);
  par.push_back(parFit[1]);
  par.push_back(parFit[2]);
  par.push_back(parFit[3]);

  Double_t trackMinDistance = CalculateDistance(tracks[0], par);

  return trackMinDistance;
}

Double_t FitTrack::FitHough(std::vector<mmTrack> track, std::vector<Double_t> &par) {
  par.clear();
  Hough2D *houghEvent = new Hough2D(track);
  houghAngleXY = houghEvent->GetMaxThetaXY();
  houghDXY = houghEvent->GetMaxDXY();
  houghAngleYZ = houghEvent->GetMaxThetaYZ();
  houghDYZ = houghEvent->GetMaxDYZ();
  delete houghEvent;

  // std::cout << houghAngleXY << '\t' << houghDXY << std::endl;
  // std::cout << houghAngleYZ << '\t' << houghDYZ << std::endl;

  par.push_back(houghDXY/cos(houghAngleXY*M_PI/180.));
  par.push_back(-sin(houghAngleXY*M_PI/180.)/cos(houghAngleXY*M_PI/180.));
  par.push_back(houghDXY/sin(houghAngleYZ*M_PI/180.));
  par.push_back(-cos(houghAngleYZ*M_PI/180.)/sin(houghAngleYZ*M_PI/180.));

  Double_t trackMinDistance = CalculateDistance(tracks[0], par);

  return trackMinDistance;
}

Double_t FitTrack::CalculateDistance(std::vector<mmTrack> track, std::vector<Double_t> par) {
  Double_t dist = 0.;

  for(UInt_t i = 0; i < track.size(); i++) {
    ROOT::Math::XYZVector xp(track[i].xPosition, track[i].yPosition, track[i].height);
    ROOT::Math::XYZVector x0(par[0], 0., par[2]);
    ROOT::Math::XYZVector x1(par[0] + par[1], 1., par[2] + par[3]);
    ROOT::Math::XYZVector u = (x1 - x0).Unit();
    dist += ((xp - x0).Cross(u)).Mag2();
  }

  return dist/static_cast<Double_t>(track.size());
}