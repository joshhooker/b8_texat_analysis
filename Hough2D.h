#ifndef HOUGH2D_H
#define HOUGH2D_H

#include <cmath>
#include <iostream>
#include <map>
#include <vector>

#include <TH2.h>

#include "TypeDef.h"

class Hough2D {
public:
  Hough2D();
  ~Hough2D();
  Hough2D(std::vector<mmTrack>);
  Hough2D(std::vector<mmTrack>, Int_t, Int_t);

  void SetPoints(std::vector<mmTrack>);
  void CalculateHough();
  void CalculateHoughXY();
  void CalculateHoughYZ();

  Double_t GetMaxThetaXY();
  Double_t GetMaxDXY();
  TH2I GetHoughDiagramXY();

  Double_t GetMaxThetaYZ();
  Double_t GetMaxDYZ();
  TH2I GetHoughDiagramYZ();

  std::vector<xy> pointsXY;
  TH2I* hHoughXY = NULL;

  std::vector<yz> pointsYZ;
  TH2I* hHoughYZ = NULL;

private:
  Int_t nBinsX, nBinsY;
  Double_t maxThetaXY, maxDXY;
  Double_t maxThetaYZ, maxDYZ;
};

#endif //HOUGH2D_H
