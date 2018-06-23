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

#ifndef HOUGHTRACK_H
#define HOUGHTRACK_H

// Hough Track class is for fitting of one particle track (i.e. proton or beam before reaction).
// We restrict the Hough space in terms of angle depending on the detector that fired

class HoughTrack {
public:
  HoughTrack();
  HoughTrack(std::vector<mmTrack> inputTrack, Int_t det, Int_t quad);
  void AddTrack(std::vector<mmTrack> inputTrack, Int_t det, Int_t quad);

  Double_t Fit();
  Double_t FitRestricted();
  Double_t CalculateDistance();

  std::vector<Double_t> GetPars() { return par;};
  Double_t GetHoughAngleXY() {return houghAngleXY;};
  Double_t GetHoughDXY() {return houghDXY;};
  Double_t GetHoughAngleYZ() {return houghAngleYZ;};
  Double_t GetHoughDYZ() {return houghDYZ;};

private:
  std::vector<mmTrack> track;
  std::vector<Double_t> par;
  Double_t minDist;
  Int_t detFired;
  Int_t quadFired;

  Double_t houghAngleXY;
  Double_t houghDXY;
  Double_t houghAngleYZ;
  Double_t houghDYZ;
};


#endif //HOUGHTRACK_H