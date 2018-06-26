#include "Hough2D.h"

Hough2D::Hough2D(): nBinsX(360), nBinsY(360) {}

Hough2D::~Hough2D() {
  delete hHoughXY;
  delete hHoughYZ;
}

Hough2D::Hough2D(std::vector<mmTrack> initPoints): nBinsX(360), nBinsY(360) {
  for(UInt_t i = 0; i < initPoints.size(); i++) {
    xy initPointsXY = {initPoints[i].xPosition, initPoints[i].yPosition};
    yz initPointsYZ = {initPoints[i].yPosition, initPoints[i].height};
    pointsXY.push_back(initPointsXY);
    pointsYZ.push_back(initPointsYZ);
  }

  CalculateHough();
}

Hough2D::Hough2D(std::vector<mmTrack> initPoints, Int_t det, Int_t quad): nBinsX(360), nBinsY(360) {
  for(UInt_t i = 0; i < initPoints.size(); i++) {
    xy initPointsXY = {initPoints[i].xPosition, initPoints[i].yPosition};
    yz initPointsYZ = {initPoints[i].yPosition, initPoints[i].height};
    pointsXY.push_back(initPointsXY);
    if(i < initPoints.size() - 3) pointsYZ.push_back(initPointsYZ);
  }

  detFired = det;

  CalculateHoughRestricted();
}

Hough2D::Hough2D(std::vector<mmTrack> initPoints, Int_t det, Int_t quad, Int_t binsX, Int_t binsY): nBinsX(binsX),
  nBinsY(binsY) {
  for(UInt_t i = 0; i < initPoints.size(); i++) {
    xy initPointsXY = {initPoints[i].xPosition, initPoints[i].yPosition};
    yz initPointsYZ = {initPoints[i].yPosition, initPoints[i].height};
    pointsXY.push_back(initPointsXY);
    pointsYZ.push_back(initPointsYZ);
  }

  detFired = det;

  CalculateHoughRestricted();
}

void Hough2D::SetPoints(std::vector<mmTrack> initPoints) {
  for(UInt_t i = 0; i < initPoints.size(); i++) {
    xy initPointsXY = {initPoints[i].xPosition, initPoints[i].yPosition};
    yz initPointsYZ = {initPoints[i].yPosition, initPoints[i].height};
    pointsXY.push_back(initPointsXY);
    pointsYZ.push_back(initPointsYZ);
  }
}

void Hough2D::SetPoints(std::vector<mmTrack> initPoints, Int_t det) {
  for(UInt_t i = 0; i < initPoints.size(); i++) {
    xy initPointsXY = {initPoints[i].xPosition, initPoints[i].yPosition};
    yz initPointsYZ = {initPoints[i].yPosition, initPoints[i].height};
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
  Double_t minDHistogram = 1000;
  Double_t maxDHistogram = -1000;
  for(Int_t j = 0; j < 180; j++) {
    Double_t cosj = cos(j*M_PI/180.);
    Double_t sinj = sin(j*M_PI/180.);
    for(UInt_t i = 0; i < pointsXY.size(); i++) {
      Double_t d = pointsXY[i].x*cosj + pointsXY[i].y*sinj;
      if(d < minDHistogram) minDHistogram = d;
      if(d > maxDHistogram) maxDHistogram = d;
    }
  }

  Double_t thetaStep = 180./static_cast<Double_t>(nBinsX);

  // Fill Hough Matrix
  // Go through the angles and find the angle with the smallest standard deviation
  hHoughXY = new TH2I("houghXY", "houghXY", nBinsX, 0, 180, nBinsY, minDHistogram, maxDHistogram);
  hHoughXY->GetXaxis()->SetTitle("#Theta"); hHoughXY->GetXaxis()->CenterTitle();
  hHoughXY->GetYaxis()->SetTitle("d"); hHoughXY->GetYaxis()->CenterTitle();

  Double_t smallestTheta;
  Double_t smallestStdDev = 10000;
  for(Double_t j = 0; j < 180; j += thetaStep) {
    if(j > 89 && j < 91) continue;
    Double_t cosj = cos(j*M_PI/180.);
    Double_t sinj = sin(j*M_PI/180.);
    std::vector<Double_t> dVector;
    Double_t mean = 0.;
    for(UInt_t i = 0; i < pointsXY.size(); i++) {
      Double_t d = pointsXY[i].x*cosj + pointsXY[i].y*sinj;
      mean += d;
      dVector.push_back(d);
      hHoughXY->Fill(j, d);
    }
    if(dVector.empty()) continue;
    mean /= static_cast<Double_t>(pointsXY.size());
    Double_t stdDev = 0.;
    for(UInt_t i = 0; i < pointsXY.size(); i++) {
      Double_t d = pointsXY[i].x*cosj + pointsXY[i].y*sinj;
      stdDev += (d - mean)*(d - mean);
    }
    stdDev /= static_cast<Double_t>(pointsXY.size() - 1.);
    if(stdDev < smallestStdDev) {
      smallestStdDev = stdDev;
      smallestTheta = j;
    }
  }

  // For the theta with the smallest standard deviation, find the d with the max value
  Int_t largestDValue = -1000;
  Double_t largestDLocation;
  Int_t binX = hHoughXY->GetXaxis()->FindBin(smallestTheta);
  for(Int_t j = 0; j < hHoughXY->GetNbinsY() + 1; j++) {
    if(hHoughXY->GetBinContent(binX, j) > largestDValue) {
      largestDValue = hHoughXY->GetBinContent(binX, j);
      largestDLocation = hHoughXY->GetYaxis()->GetBinCenter(j);
    }
  }
  // Check if there are any more values with the maximum, then take the average
  Double_t dSum = 0.;
  Int_t dCount = 0;
  std::vector<Double_t> dVec;
  for(Int_t j = 0; j < hHoughXY->GetNbinsY() + 1; j++) {
    if(hHoughXY->GetBinContent(binX, j) == largestDValue) {
      dSum += hHoughXY->GetYaxis()->GetBinCenter(j);
      dCount++;
    }
  }

  largestDLocation = dSum/ static_cast<Double_t>(dCount);

  maxThetaXY = smallestTheta;
  maxDXY = largestDLocation;

//  std::cout << "XY: Smallest std dev at " << maxThetaXY
//            << " with value " << smallestStdDev << " and max D at " << maxDXY << std::endl;
}

void Hough2D::CalculateHoughYZ() {
  Double_t minDHistogram = 1000;
  Double_t maxDHistogram = -1000;
  for(Int_t j = 0; j < 180; j++) {
    Double_t cosj = cos(j*M_PI/180.);
    Double_t sinj = sin(j*M_PI/180.);
    for(UInt_t i = 0; i < pointsYZ.size(); i++) {
      Double_t d = pointsYZ[i].y*cosj + pointsYZ[i].z*sinj;
      if(d < minDHistogram) minDHistogram = d;
      if(d > maxDHistogram) maxDHistogram = d;
    }
  }

  Double_t thetaStep = 180./static_cast<Double_t>(nBinsX);

  // Fill Hough Matrix
  // Go through the angles and find the angle with the smallest standard deviation
  hHoughYZ = new TH2I("houghYZ", "houghYZ", nBinsX, 0, 180, nBinsY, minDHistogram, maxDHistogram);
  hHoughYZ->GetXaxis()->SetTitle("#Theta"); hHoughYZ->GetXaxis()->CenterTitle();
  hHoughYZ->GetYaxis()->SetTitle("d"); hHoughYZ->GetYaxis()->CenterTitle();

  Double_t smallestTheta;
  Double_t smallestStdDev = 10000;
  for(Double_t j = 0; j < 180; j += thetaStep) {
    if(j > 89 && j < 91) continue;
    Double_t cosj = cos(j*M_PI/180.);
    Double_t sinj = sin(j*M_PI/180.);
    std::vector<Double_t> dVector;
    Double_t mean = 0.;
    for(UInt_t i = 0; i < pointsYZ.size(); i++) {
      Double_t d = pointsYZ[i].y*cosj + pointsYZ[i].z*sinj;
      mean += d;
      dVector.push_back(d);
      hHoughYZ->Fill(j, d);
    }
    if(dVector.empty()) continue;
    mean /= static_cast<Double_t>(pointsYZ.size());
    Double_t stdDev = 0.;
    for(UInt_t i = 0; i < pointsYZ.size(); i++) {
      Double_t d = pointsYZ[i].y*cosj + pointsYZ[i].z*sinj;
      stdDev += (d - mean)*(d - mean);
    }
    stdDev /= static_cast<Double_t>(pointsYZ.size() - 1.);
    if(stdDev < smallestStdDev) {
      smallestStdDev = stdDev;
      smallestTheta = j;
    }
  }

  // For the theta with the smallest standard deviation, find the d with the max value
  Double_t largestDValue = -1000;
  Double_t largestDLocation;
  Int_t binX = hHoughYZ->GetXaxis()->FindBin(smallestTheta);
  for(Int_t j = 0; j < hHoughYZ->GetNbinsY() + 1; j++) {
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
  Int_t minAngle = 0;
  Int_t maxAngle = 0;
  if(detFired < 4) {
    minAngle = 91;
    maxAngle = 179;
  }
  else if(detFired > 5 && detFired < 10) {
    minAngle = 1;
    maxAngle = 89;
  }

  Double_t minDHistogram = 1000;
  Double_t maxDHistogram = -1000;
  for(Int_t j = minAngle; j <= maxAngle; j++) {
    Double_t cosj = cos(j*M_PI/180.);
    Double_t sinj = sin(j*M_PI/180.);
    for(UInt_t i = 0; i < pointsXY.size(); i++) {
      Double_t d = pointsXY[i].x*cosj + pointsXY[i].y*sinj;
      if(d < minDHistogram) minDHistogram = d;
      if(d > maxDHistogram) maxDHistogram = d;
    }
  }

  Double_t thetaRange = static_cast<Double_t>(maxAngle) - static_cast<Double_t>(minAngle);

  Double_t thetaStep = thetaRange/static_cast<Double_t>(nBinsX);

  // Fill Hough Matrix
  // Go through the angles and find the angle with the smallest standard deviation
  hHoughXY = new TH2I("houghXY", "houghXY", nBinsX, 0, 180, nBinsY, minDHistogram, maxDHistogram);
  hHoughXY->GetXaxis()->SetTitle("#Theta"); hHoughXY->GetXaxis()->CenterTitle();
  hHoughXY->GetYaxis()->SetTitle("d"); hHoughXY->GetYaxis()->CenterTitle();

  Double_t smallestTheta;
  Double_t smallestStdDev = 10000;
  for(Double_t j = minAngle; j <= maxAngle; j += thetaStep) {
    Double_t cosj = cos(j*M_PI/180.);
    Double_t sinj = sin(j*M_PI/180.);
    std::vector<Double_t> dVector;
    Double_t mean = 0.;
    for(UInt_t i = 0; i < pointsXY.size(); i++) {
      Double_t d = pointsXY[i].x*cosj + pointsXY[i].y*sinj;
      mean += d;
      dVector.push_back(d);
      hHoughXY->Fill(j, d);
    }
    if(dVector.empty()) continue;
    mean /= static_cast<Double_t>(pointsXY.size());
    Double_t stdDev = 0.;
    for(UInt_t i = 0; i < pointsXY.size(); i++) {
      Double_t d = pointsXY[i].x*cosj + pointsXY[i].y*sinj;
      stdDev += (d - mean)*(d - mean);
    }
    stdDev /= static_cast<Double_t>(pointsXY.size() - 1.);
    if(stdDev < smallestStdDev) {
      smallestStdDev = stdDev;
      smallestTheta = j;
    }
  }

  // For the theta with the smallest standard deviation, find the d with the max value
  Int_t largestDValue = -1000;
  Double_t largestDLocation;
  Int_t binX = hHoughXY->GetXaxis()->FindBin(smallestTheta);
  for(Int_t j = 0; j < hHoughXY->GetNbinsY() + 1; j++) {
    if(hHoughXY->GetBinContent(binX, j) > largestDValue) {
      largestDValue = hHoughXY->GetBinContent(binX, j);
      largestDLocation = hHoughXY->GetYaxis()->GetBinCenter(j);
    }
  }
  // Check if there are any more values with the maximum, then take the average
  Double_t dSum = 0.;
  Int_t dCount = 0;
  std::vector<Double_t> dVec;
  for(Int_t j = 0; j < hHoughXY->GetNbinsY() + 1; j++) {
    if(hHoughXY->GetBinContent(binX, j) == largestDValue) {
      dSum += hHoughXY->GetYaxis()->GetBinCenter(j);
      dCount++;
    }
  }

  largestDLocation = dSum/ static_cast<Double_t>(dCount);

  maxThetaXY = smallestTheta;
  maxDXY = largestDLocation;

//  std::cout << "XY: Smallest std dev at " << maxThetaXY
//            << " with value " << smallestStdDev << " and max D at " << maxDXY << std::endl;
}

void Hough2D::CalculateHoughYZRestricted() {
  Double_t minDHistogram = 1000;
  Double_t maxDHistogram = -1000;
  for(Int_t j = 0; j < 180; j++) {
    Double_t cosj = cos(j*M_PI/180.);
    Double_t sinj = sin(j*M_PI/180.);
    for(UInt_t i = 0; i < pointsYZ.size(); i++) {
      Double_t d = pointsYZ[i].y*cosj + pointsYZ[i].z*sinj;
      if(d < minDHistogram) minDHistogram = d;
      if(d > maxDHistogram) maxDHistogram = d;
    }
  }

  Double_t thetaStep = 180./static_cast<Double_t>(nBinsX);

  // Fill Hough Matrix
  // Go through the angles and find the angle with the smallest standard deviation
  hHoughYZ = new TH2I("houghYZ", "houghYZ", nBinsX, 0, 180, nBinsY, minDHistogram, maxDHistogram);
  hHoughYZ->GetXaxis()->SetTitle("#Theta"); hHoughYZ->GetXaxis()->CenterTitle();
  hHoughYZ->GetYaxis()->SetTitle("d"); hHoughYZ->GetYaxis()->CenterTitle();

  Double_t smallestTheta;
  Double_t smallestStdDev = 10000;
  for(Double_t j = 0; j < 180; j += thetaStep) {
    if(j > 89 && j < 91) continue;
    Double_t cosj = cos(j*M_PI/180.);
    Double_t sinj = sin(j*M_PI/180.);
    std::vector<Double_t> dVector;
    Double_t mean = 0.;
    for(UInt_t i = 0; i < pointsYZ.size(); i++) {
      Double_t d = pointsYZ[i].y*cosj + pointsYZ[i].z*sinj;
      mean += d;
      dVector.push_back(d);
      hHoughYZ->Fill(j, d);
    }
    if(dVector.empty()) continue;
    mean /= static_cast<Double_t>(pointsYZ.size());
    Double_t stdDev = 0.;
    for(UInt_t i = 0; i < pointsYZ.size(); i++) {
      Double_t d = pointsYZ[i].y*cosj + pointsYZ[i].z*sinj;
      stdDev += (d - mean)*(d - mean);
    }
    stdDev /= static_cast<Double_t>(pointsYZ.size() - 1.);
    if(stdDev < smallestStdDev) {
      smallestStdDev = stdDev;
      smallestTheta = j;
    }
  }

  // For the theta with the smallest standard deviation, find the d with the max value
  Double_t largestDValue = -1000;
  Double_t largestDLocation;
  Int_t binX = hHoughYZ->GetXaxis()->FindBin(smallestTheta);
  for(Int_t j = 0; j < hHoughYZ->GetNbinsY() + 1; j++) {
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

Double_t Hough2D::GetMaxThetaXY() {
  return maxThetaXY;
}

Double_t Hough2D::GetMaxDXY() {
  return maxDXY;
}

TH2I Hough2D::GetHoughDiagramXY() {
  return *hHoughXY;
}

Double_t Hough2D::GetMaxThetaYZ() {
  return maxThetaYZ;
}

Double_t Hough2D::GetMaxDYZ() {
  return maxDYZ;
}

TH2I Hough2D::GetHoughDiagramYZ() {
  return *hHoughYZ;
}