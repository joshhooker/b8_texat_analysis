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
extern double distance2(double x, double y, double z, const double *p);
#endif //DISTANCE2FUNC

#ifndef SUMDISTANCE2FUNC
#define SUMDISTANCE2FUNC
// function Object to be minimized
struct SumDistance2 {
  // the TGraph is a data memeber of the object
  TGraph2D *fGraph;
  SumDistance2(TGraph2D *g) : fGraph(g) {}

  // implementation of the function to be minimized
  double operator() (const double *par) {
    assert(fGraph != 0);
    double *x = fGraph->GetX();
    double *y = fGraph->GetY();
    double *z = fGraph->GetZ();
    int npoints = fGraph->GetN();
    double sum = 0;
    for(uint i = 0; i < npoints; i++) {
      double d = distance2(x[i], y[i], z[i], par);
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
  HoughTrack(std::vector<mmTrack> inputTrack, int det, int quad);
  void AddTrack(std::vector<mmTrack> inputTrack, int det, int quad);

  double Fit();
  double FitRestricted();
  double CalculateDistance();

  std::vector<double> GetPars() { return par;};
  double GetHoughAngleXY() {return houghAngleXY;};
  double GetHoughDXY() {return houghDXY;};
  double GetHoughAngleYZ() {return houghAngleYZ;};
  double GetHoughDYZ() {return houghDYZ;};

private:
  std::vector<mmTrack> track;
  std::vector<double> par;
  double minDist;
  int detFired;
  int quadFired;

  double houghAngleXY;
  double houghDXY;
  double houghAngleYZ;
  double houghDYZ;
};


#endif //HOUGHTRACK_H
