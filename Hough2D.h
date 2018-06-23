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
  Hough2D(std::vector<mmTrack> initPoints);

  Hough2D(std::vector<mmTrack> initPoints, Int_t det, Int_t quad);
  Hough2D(std::vector<mmTrack> initPoints, Int_t det, Int_t quad, Int_t binsX, Int_t binsY);

  void SetPoints(std::vector<mmTrack> initPoints);
  void SetPoints(std::vector<mmTrack> initPoints, Int_t det);

  void CalculateHough();
  void CalculateHoughXY();
  void CalculateHoughYZ();

  void CalculateHoughRestricted();
  void CalculateHoughXYRestricted();
  void CalculateHoughYZRestricted();

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
  Int_t detFired;
  Double_t maxThetaXY, maxDXY;
  Double_t maxThetaYZ, maxDYZ;

  void InitHough2D();
  Float_t siXPosForward[10][4];
};

// This is a shitty way to do this since it's defined in Spectra class
inline void Hough2D::InitHough2D() {
  // X positions of Forward Si Detector Quadrants
  siXPosForward[0][0] = 124. + 12.5;
  siXPosForward[0][1] = 124. - 12.5;
  siXPosForward[0][2] = 124. - 12.5;
  siXPosForward[0][3] = 124. + 12.5;
  siXPosForward[1][0] = 124. + 12.5;
  siXPosForward[1][1] = 124. - 12.5;
  siXPosForward[1][2] = 124. - 12.5;
  siXPosForward[1][3] = 124. + 12.5;
  siXPosForward[2][0] = 62.2 + 12.5;
  siXPosForward[2][1] = 62.2 + 12.5;
  siXPosForward[2][2] = 62.2 - 12.5;
  siXPosForward[2][3] = 62.2 - 12.5;
  siXPosForward[3][0] = 62.2 + 12.5;
  siXPosForward[3][1] = 62.2 - 12.5;
  siXPosForward[3][2] = 62.2 - 12.5;
  siXPosForward[3][3] = 62.2 + 12.5;
  siXPosForward[4][0] = 0. + 12.5;
  siXPosForward[4][1] = 0. + 12.5;
  siXPosForward[4][2] = 0. + 12.5;
  siXPosForward[4][3] = 0. + 12.5;
  siXPosForward[5][0] = 0. + 12.5;
  siXPosForward[5][1] = 0. + 12.5;
  siXPosForward[5][2] = 0. + 12.5;
  siXPosForward[5][3] = 0. + 12.5;
  siXPosForward[6][0] = 62.2 - 12.5;
  siXPosForward[6][1] = 62.2 - 12.5;
  siXPosForward[6][2] = 62.2 + 12.5;
  siXPosForward[6][3] = 62.2 + 12.5;
  siXPosForward[7][0] = 62.2 - 12.5;
  siXPosForward[7][1] = 62.2 + 12.5;
  siXPosForward[7][2] = 62.2 + 12.5;
  siXPosForward[7][3] = 62.2 - 12.5;
  siXPosForward[8][0] = 124. - 12.5;
  siXPosForward[8][1] = 124. + 12.5;
  siXPosForward[8][2] = 124. + 12.5;
  siXPosForward[8][3] = 124. - 12.5;
  siXPosForward[9][0] = 124. - 12.5;
  siXPosForward[9][1] = 124. + 12.5;
  siXPosForward[9][2] = 124. + 12.5;
  siXPosForward[9][3] = 124. - 12.5;
}

#endif //HOUGH2D_H
