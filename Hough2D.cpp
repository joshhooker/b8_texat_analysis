#include "Hough2D.h"

Hough2D::Hough2D(): nBinsX(360), nBinsY(360) {}

Hough2D::~Hough2D() {
  delete hHoughXY;
  delete hHoughYZ;
}

Hough2D::Hough2D(std::vector<mmTrack> initPoints): nBinsX(360), nBinsY(1440) {
  for(auto point : initPoints) {
    xy initPointsXY = {point.xPosition, point.yPosition};
    yz initPointsYZ = {point.yPosition, point.height};
    pointsXY.push_back(initPointsXY);
    pointsYZ.push_back(initPointsYZ);
  }

  CalculateHough();
}

Hough2D::Hough2D(std::vector<mmTrack> initPoints, int det, int quad): nBinsX(360), nBinsY(1440) {
  for(auto point : initPoints) {
    xy initPointsXY = {point.xPosition, point.yPosition};
    yz initPointsYZ = {point.yPosition, point.height};
    pointsXY.push_back(initPointsXY);
    pointsYZ.push_back(initPointsYZ);
  }

  detFired = det;
  quadFired = quad;

  CalculateHoughRestricted();
}

Hough2D::Hough2D(std::vector<mmTrack> initPoints, int det, int quad, int binsX, int binsY): nBinsX(binsX),
  nBinsY(binsY) {
  for(auto point : initPoints) {
    xy initPointsXY = {point.xPosition, point.yPosition};
    yz initPointsYZ = {point.yPosition, point.height};
    pointsXY.push_back(initPointsXY);
    pointsYZ.push_back(initPointsYZ);
  }

  detFired = det;
  quadFired = quad;

  CalculateHoughRestricted();
}

Hough2D::Hough2D(std::vector<mmTrack> initPoints, int det, int quad, int help): nBinsX(360), nBinsY(1440) {
  for(auto point : initPoints) {
    xy initPointsXY = {point.xPosition, point.yPosition};
    yz initPointsYZ = {point.yPosition, point.height};
    pointsXY.push_back(initPointsXY);
    pointsYZ.push_back(initPointsYZ);
  }
  // Single Si Help
  if(help == 0) {
    xy initPointsXY = {siXPosForward[det][quad], 275.34};
    pointsXY.push_back(initPointsXY);
    pointsXY.push_back(initPointsXY);
  }
  // Full Si Help (3 points)
  else if(help == 1) {
    xy initPointsXY1 = {siXPosForward[det][quad] - 12.5, 275.34};
    xy initPointsXY2 = {siXPosForward[det][quad], 275.34};
    xy initPointsXY3 = {siXPosForward[det][quad] + 12.5, 275.34};
    pointsXY.push_back(initPointsXY1);
    pointsXY.push_back(initPointsXY2);
    pointsXY.push_back(initPointsXY3);
    pointsXY.push_back(initPointsXY1);
    pointsXY.push_back(initPointsXY2);
    pointsXY.push_back(initPointsXY3);
    pointsXY.push_back(initPointsXY1);
    pointsXY.push_back(initPointsXY2);
    pointsXY.push_back(initPointsXY3);
  }

  detFired = det;
  quadFired = quad;

  CalculateHoughRestricted();
}



void Hough2D::SetPoints(std::vector<mmTrack> initPoints) {
  for(auto point : initPoints) {
    xy initPointsXY = {point.xPosition, point.yPosition};
    yz initPointsYZ = {point.yPosition, point.height};
    pointsXY.push_back(initPointsXY);
    pointsYZ.push_back(initPointsYZ);
  }
}

void Hough2D::SetPoints(std::vector<mmTrack> initPoints, Int_t det) {
  for(auto point : initPoints) {
    xy initPointsXY = {point.xPosition, point.yPosition};
    yz initPointsYZ = {point.yPosition, point.height};
    pointsXY.push_back(initPointsXY);
    pointsYZ.push_back(initPointsYZ);
  }
  detFired = det;
}

void Hough2D::CalculateHough() {
  CalculateHoughXY();
  CalculateHoughYZ();
}

void Hough2D::CalculateHoughXY() {
  double minDHistogram = 1000;
  double maxDHistogram = -1000;
  for(uint j = 0; j < 180; j++) {
    double cosj = cos(j*M_PI/180.);
    double sinj = sin(j*M_PI/180.);
    for(uint i = 0; i < pointsXY.size(); i++) {
      double d = pointsXY[i].x*cosj + pointsXY[i].y*sinj;
      if(d < minDHistogram) minDHistogram = d;
      if(d > maxDHistogram) maxDHistogram = d;
    }
  }

  double thetaStep = 180./static_cast<double>(nBinsX);

  double diffD = maxDHistogram - minDHistogram;
  int nBinsY = static_cast<int>(diffD/0.005);

  // Fill Hough Matrix
  // Go through the angles and find the angle with the smallest standard deviation
  hHoughXY = new TH2I("houghXY", "houghXY", nBinsX, 0, 180, nBinsY, minDHistogram, maxDHistogram);
  hHoughXY->GetXaxis()->SetTitle("#Theta"); hHoughXY->GetXaxis()->CenterTitle();
  hHoughXY->GetYaxis()->SetTitle("d"); hHoughXY->GetYaxis()->CenterTitle();

  double smallestTheta;
  double smallestStdDev = 10000;
  for(double j = 0; j < 180; j += thetaStep) {
    if(j > 89 && j < 91) continue;
    double cosj = cos(j*M_PI/180.);
    double sinj = sin(j*M_PI/180.);
    std::vector<double> dVector;
    double mean = 0.;
    for(uint i = 0; i < pointsXY.size(); i++) {
      double d = pointsXY[i].x*cosj + pointsXY[i].y*sinj;
      mean += d;
      dVector.push_back(d);
      hHoughXY->Fill(j, d);
    }
    if(dVector.empty()) continue;
    mean /= static_cast<double>(pointsXY.size());
    double stdDev = 0.;
    for(uint i = 0; i < pointsXY.size(); i++) {
      double d = pointsXY[i].x*cosj + pointsXY[i].y*sinj;
      stdDev += (d - mean)*(d - mean);
    }
    stdDev /= static_cast<double>(pointsXY.size() - 1.);
    if(stdDev < smallestStdDev) {
      smallestStdDev = stdDev;
      smallestTheta = j;
    }
  }

  // For the theta with the smallest standard deviation, find the d with the max value
  int largestDValue = -1000;
  double largestDLocation;
  int binX = hHoughXY->GetXaxis()->FindBin(smallestTheta);
  for(int j = 0; j < hHoughXY->GetNbinsY() + 1; j++) {
    if(hHoughXY->GetBinContent(binX, j) > largestDValue) {
      largestDValue = hHoughXY->GetBinContent(binX, j);
      largestDLocation = hHoughXY->GetYaxis()->GetBinCenter(j);
    }
  }
  // Check if there are any more values with the maximum, then take the average
  double dSum = 0.;
  int dCount = 0;
  std::vector<double> dVec;
  for(int j = 0; j < hHoughXY->GetNbinsY() + 1; j++) {
    if(hHoughXY->GetBinContent(binX, j) == largestDValue) {
      dSum += hHoughXY->GetYaxis()->GetBinCenter(j);
      dCount++;
    }
  }

  largestDLocation = dSum/ static_cast<double>(dCount);

  maxThetaXY = smallestTheta;
  maxDXY = largestDLocation;

//  std::cout << "XY: Smallest std dev at " << maxThetaXY
//            << " with value " << smallestStdDev << " and max D at " << maxDXY << std::endl;
}

void Hough2D::CalculateHoughYZ() {
  double minDHistogram = 1000;
  double maxDHistogram = -1000;
  for(uint j = 0; j < 180; j++) {
    double cosj = cos(j*M_PI/180.);
    double sinj = sin(j*M_PI/180.);
    for(uint i = 0; i < pointsYZ.size(); i++) {
      double d = pointsYZ[i].y*cosj + pointsYZ[i].z*sinj;
      if(d < minDHistogram) minDHistogram = d;
      if(d > maxDHistogram) maxDHistogram = d;
    }
  }

  double thetaStep = 180./static_cast<double>(nBinsX);

  // Fill Hough Matrix
  // Go through the angles and find the angle with the smallest standard deviation
  hHoughYZ = new TH2I("houghYZ", "houghYZ", nBinsX, 0, 180, nBinsY, minDHistogram, maxDHistogram);
  hHoughYZ->GetXaxis()->SetTitle("#Theta"); hHoughYZ->GetXaxis()->CenterTitle();
  hHoughYZ->GetYaxis()->SetTitle("d"); hHoughYZ->GetYaxis()->CenterTitle();

  double smallestTheta;
  double smallestStdDev = 10000;
  for(double j = 0; j < 180; j += thetaStep) {
    if(j > 89 && j < 91) continue;
    double cosj = cos(j*M_PI/180.);
    double sinj = sin(j*M_PI/180.);
    std::vector<double> dVector;
    double mean = 0.;
    for(uint i = 0; i < pointsYZ.size(); i++) {
      double d = pointsYZ[i].y*cosj + pointsYZ[i].z*sinj;
      mean += d;
      dVector.push_back(d);
      hHoughYZ->Fill(j, d);
    }
    if(dVector.empty()) continue;
    mean /= static_cast<double>(pointsYZ.size());
    double stdDev = 0.;
    for(uint i = 0; i < pointsYZ.size(); i++) {
      double d = pointsYZ[i].y*cosj + pointsYZ[i].z*sinj;
      stdDev += (d - mean)*(d - mean);
    }
    stdDev /= static_cast<double>(pointsYZ.size() - 1.);
    if(stdDev < smallestStdDev) {
      smallestStdDev = stdDev;
      smallestTheta = j;
    }
  }

  // For the theta with the smallest standard deviation, find the d with the max value
  double largestDValue = -1000;
  double largestDLocation;
  int binX = hHoughYZ->GetXaxis()->FindBin(smallestTheta);
  for(int j = 0; j < hHoughYZ->GetNbinsY() + 1; j++) {
    if(hHoughYZ->GetBinContent(binX, j) > largestDValue) {
      largestDValue = hHoughYZ->GetBinContent(binX, j);
      largestDLocation = hHoughYZ->GetYaxis()->GetBinCenter(j);
    }
  }

  maxThetaYZ = smallestTheta;
  maxDYZ = largestDLocation;

//  std::cout << "XY: Smallest std dev at " << maxThetaYZ
//            << " with value " << smallestStdDev << " and max D at " << maxDYZ << std::endl;
}

void Hough2D::CalculateHoughRestricted() {
  CalculateHoughXYRestricted();
  CalculateHoughYZRestricted();
}

void Hough2D::CalculateHoughXYRestricted() {
  int minAngle = 0;
  int maxAngle = 0;
  if(detFired < 4) {
    minAngle = 91;
    maxAngle = 179;
  }
  else if(detFired > 5 && detFired < 10) {
    minAngle = 1;
    maxAngle = 89;
  }

  double minDHistogram = 1000;
  double maxDHistogram = -1000;
  for(int j = minAngle; j <= maxAngle; j++) {
    double cosj = cos(j*M_PI/180.);
    double sinj = sin(j*M_PI/180.);
    for(uint i = 0; i < pointsXY.size(); i++) {
      double d = pointsXY[i].x*cosj + pointsXY[i].y*sinj;
      if(d < minDHistogram) minDHistogram = d;
      if(d > maxDHistogram) maxDHistogram = d;
    }
  }

  double thetaRange = static_cast<double>(maxAngle) - static_cast<double>(minAngle);

  double thetaStep = thetaRange/static_cast<double>(nBinsX);

  // Fill Hough Matrix
  // Go through the angles and find the angle with the smallest standard deviation
  hHoughXY = new TH2I("houghXY", "houghXY", nBinsX, 0, 180, nBinsY, minDHistogram, maxDHistogram);
  hHoughXY->GetXaxis()->SetTitle("#Theta"); hHoughXY->GetXaxis()->CenterTitle();
  hHoughXY->GetYaxis()->SetTitle("d"); hHoughXY->GetYaxis()->CenterTitle();

  double smallestTheta;
  double smallestStdDev = 10000;
  for(double j = minAngle; j <= maxAngle; j += thetaStep) {
    double cosj = cos(j*M_PI/180.);
    double sinj = sin(j*M_PI/180.);
    std::vector<double> dVector;
    double mean = 0.;
    for(uint i = 0; i < pointsXY.size(); i++) {
      double d = pointsXY[i].x*cosj + pointsXY[i].y*sinj;
      mean += d;
      dVector.push_back(d);
      hHoughXY->Fill(j, d);
    }
    if(dVector.empty()) continue;
    mean /= static_cast<double>(pointsXY.size());
    double stdDev = 0.;
    for(uint i = 0; i < pointsXY.size(); i++) {
      double d = pointsXY[i].x*cosj + pointsXY[i].y*sinj;
      stdDev += (d - mean)*(d - mean);
    }
    stdDev /= static_cast<double>(pointsXY.size() - 1.);
    if(stdDev < smallestStdDev) {
      smallestStdDev = stdDev;
      smallestTheta = j;
    }
  }

  // For the theta with the smallest standard deviation, find the d with the max value
  int largestDValue = -1000;
  double largestDLocation;
  int binX = hHoughXY->GetXaxis()->FindBin(smallestTheta);
  for(int j = 0; j < hHoughXY->GetNbinsY() + 1; j++) {
    if(hHoughXY->GetBinContent(binX, j) > largestDValue) {
      largestDValue = hHoughXY->GetBinContent(binX, j);
      largestDLocation = hHoughXY->GetYaxis()->GetBinCenter(j);
    }
  }
  // Check if there are any more values with the maximum, then take the average
  double dSum = 0.;
  int dCount = 0;
  std::vector<double> dVec;
  for(int j = 0; j < hHoughXY->GetNbinsY() + 1; j++) {
    if(hHoughXY->GetBinContent(binX, j) == largestDValue) {
      dSum += hHoughXY->GetYaxis()->GetBinCenter(j);
      dCount++;
    }
  }

  largestDLocation = dSum/ static_cast<double>(dCount);

  maxThetaXY = smallestTheta;
  maxDXY = largestDLocation;

//  std::cout << "XY: Smallest std dev at " << maxThetaXY
//            << " with value " << smallestStdDev << " and max D at " << maxDXY << std::endl;
}

void Hough2D::CalculateHoughYZRestricted() {
  double minDHistogram = 1000;
  double maxDHistogram = -1000;
  for(int j = 0; j < 180; j++) {
    double cosj = cos(j*M_PI/180.);
    double sinj = sin(j*M_PI/180.);
    for(uint i = 0; i < pointsYZ.size(); i++) {
      double d = pointsYZ[i].y*cosj + pointsYZ[i].z*sinj;
      if(d < minDHistogram) minDHistogram = d;
      if(d > maxDHistogram) maxDHistogram = d;
    }
  }

  double thetaStep = 180./static_cast<double>(nBinsX);

  // Fill Hough Matrix
  // Go through the angles and find the angle with the smallest standard deviation
  hHoughYZ = new TH2I("houghYZ", "houghYZ", nBinsX, 0, 180, nBinsY, minDHistogram, maxDHistogram);
  hHoughYZ->GetXaxis()->SetTitle("#Theta"); hHoughYZ->GetXaxis()->CenterTitle();
  hHoughYZ->GetYaxis()->SetTitle("d"); hHoughYZ->GetYaxis()->CenterTitle();

  double smallestTheta;
  double smallestStdDev = 10000;
  for(double j = 0; j < 180; j += thetaStep) {
    if(j > 89 && j < 91) continue;
    double cosj = cos(j*M_PI/180.);
    double sinj = sin(j*M_PI/180.);
    std::vector<double> dVector;
    double mean = 0.;
    for(uint i = 0; i < pointsYZ.size(); i++) {
      double d = pointsYZ[i].y*cosj + pointsYZ[i].z*sinj;
      mean += d;
      dVector.push_back(d);
      hHoughYZ->Fill(j, d);
    }
    if(dVector.empty()) continue;
    mean /= static_cast<double>(pointsYZ.size());
    double stdDev = 0.;
    for(uint i = 0; i < pointsYZ.size(); i++) {
      double d = pointsYZ[i].y*cosj + pointsYZ[i].z*sinj;
      stdDev += (d - mean)*(d - mean);
    }
    stdDev /= static_cast<double>(pointsYZ.size() - 1.);
    if(stdDev < smallestStdDev) {
      smallestStdDev = stdDev;
      smallestTheta = j;
    }
  }

  // For the theta with the smallest standard deviation, find the d with the max value
  double largestDValue = -1000;
  double largestDLocation;
  int binX = hHoughYZ->GetXaxis()->FindBin(smallestTheta);
  for(int j = 0; j < hHoughYZ->GetNbinsY() + 1; j++) {
    if(hHoughYZ->GetBinContent(binX, j) > largestDValue) {
      largestDValue = hHoughYZ->GetBinContent(binX, j);
      largestDLocation = hHoughYZ->GetYaxis()->GetBinCenter(j);
    }
  }

  maxThetaYZ = smallestTheta;
  maxDYZ = largestDLocation;

//  std::cout << "XY: Smallest std dev at " << maxThetaYZ
//            << " with value " << smallestStdDev << " and max D at " << maxDYZ << std::endl;
}

double Hough2D::GetMaxThetaXY() {
  return maxThetaXY;
}

double Hough2D::GetMaxDXY() {
  return maxDXY;
}

TH2I Hough2D::GetHoughDiagramXY() {
  return *hHoughXY;
}

double Hough2D::GetMaxThetaYZ() {
  return maxThetaYZ;
}

double Hough2D::GetMaxDYZ() {
  return maxDYZ;
}

TH2I Hough2D::GetHoughDiagramYZ() {
  return *hHoughYZ;
}
