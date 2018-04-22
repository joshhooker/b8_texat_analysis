#include <TMath.h>
#include <TGraph2D.h>
#include <TROOT.h>

#include <Fit/Fitter.h>

#include <Math/Functor.h>
#include <Math/Vector3D.h>

#include "Hough2D.h"
#include "TypeDef.h"

#ifndef DISTANCE2FUNC
#define DISTANCE2FUNC
extern Double_t distance2(Double_t x, Double_t y, Double_t z, const Double_t *p);
#endif //DISTANCE2FUNC

#ifndef SUMDISTANCE2FUNC
#define SUMDISTANCE2FUNC
// function Object to be minimized
struct SumDistance2 {
  // the TGraph is a data memeber of the object
  TGraph2D *fGraph;
  SumDistance2(TGraph2D *g) : fGraph(g) {}

  // implementation of the function to be minimized
  Double_t operator() (const Double_t *par) {
    assert(fGraph != 0);
    Double_t *x = fGraph->GetX();
    Double_t *y = fGraph->GetY();
    Double_t *z = fGraph->GetZ();
    Int_t npoints = fGraph->GetN();
    Double_t sum = 0;
    for(Int_t i = 0; i < npoints; i++) {
      Double_t d = distance2(x[i], y[i], z[i], par);
      sum += d;
    }
    return sum;
  }
};
#endif //SUMDISTANCE2FUNC

#ifndef FITTRACK_H
#define FITTRACK_H

// Fit Track class is for fitting of one particle track (i.e. proton or beam before reaction).
// You can add separate definitions of the track (for instance, strip chain matched by time0 method or timeSlope method) and it fits them.
// Then there is a comparison of the fitted track to the first track in the variable 'tracks'.
// It is a good idea for tracks involving strips and chains to have the time0 as the first track in 'tracks'.
// Each track is fitted by a 2D Hough Transform and by minimization of points to line.
// The minimum distance from the line to each point in the track
// for each fit type and each track is calculated for the first track in the variable 'tracks'.

class FitTrack {
public:
  FitTrack();
  FitTrack(std::vector<mmTrack> track);
  FitTrack(std::vector<std::vector<mmTrack> > listOfTracks);
  void AddTrack(std::vector<mmTrack> track);
  Double_t MakeBestFit();
  Double_t MakeFit();
  std::vector<Double_t> GetPars();
  void DisableLineFit();
  void DisableHoughFit();

private:
  Double_t Fit(std::vector<mmTrack> track, std::vector<Double_t> &par);
  Double_t FitLine(std::vector<mmTrack> track, std::vector<Double_t> &par);
  Double_t FitHough(std::vector<mmTrack> track, std::vector<Double_t> &par);
  Double_t CalculateDistance(std::vector<mmTrack> track, std::vector<Double_t> par);
  std::vector<std::vector<mmTrack> > tracks;
  std::vector<Double_t> parFit;
  Double_t minDist;
  Bool_t fitLine;
  Bool_t fitHough;
};


#endif //FITTRACK_H
