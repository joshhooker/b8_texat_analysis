#include "HoughTrack.h"

// Calculate distance line-point
double distance2(double x, double y, double z, const double *p) {
  // distance line point is D = | (xp-x0) cross ux |
  // where ux is direction of line and x0 is a point in the line (like t = 0)
  ROOT::Math::XYZVector xp(x, y, z);
  ROOT::Math::XYZVector x0(p[0], 0., p[2]);
  ROOT::Math::XYZVector x1(p[0] + p[1], 1., p[2] + p[3]);
  ROOT::Math::XYZVector u = (x1 - x0).Unit();
  double d2 = ((xp - x0).Cross(u)).Mag2();
  return d2;
}

HoughTrack::HoughTrack(): minDist(1e8) {};

HoughTrack::HoughTrack(std::vector<mmTrack> inputTrack, int det, int quad): minDist(1e8) {
  track.clear();
  track = inputTrack;
  detFired = det;
  quadFired = quad;
}

void HoughTrack::AddTrack(std::vector<mmTrack> inputTrack, int det, int quad) {
  track.clear();
  track = inputTrack;
  detFired = det;
  quadFired = quad;
}

double HoughTrack::Fit() {
  par.clear();
  Hough2D *houghEvent = new Hough2D(track);
  houghAngleXY = houghEvent->GetMaxThetaXY();
  houghDXY = houghEvent->GetMaxDXY();
  houghAngleYZ = houghEvent->GetMaxThetaYZ();
  houghDYZ = houghEvent->GetMaxDYZ();
  delete houghEvent;

  par.push_back(houghDXY/cos(houghAngleXY*M_PI/180.));
  par.push_back(-sin(houghAngleXY*M_PI/180.)/cos(houghAngleXY*M_PI/180.));
  par.push_back(houghDYZ/sin(houghAngleYZ*M_PI/180.));
  par.push_back(-cos(houghAngleYZ*M_PI/180.)/sin(houghAngleYZ*M_PI/180.));

  double trackMinDistance = CalculateDistance();

  return trackMinDistance;
}

double HoughTrack::FitRestricted() {
  par.clear();
  Hough2D *houghEvent = new Hough2D(track, detFired, quadFired);
  houghAngleXY = houghEvent->GetMaxThetaXY();
  houghDXY = houghEvent->GetMaxDXY();
  houghAngleYZ = houghEvent->GetMaxThetaYZ();
  houghDYZ = houghEvent->GetMaxDYZ();
  delete houghEvent;

  par.push_back(houghDXY/cos(houghAngleXY*M_PI/180.));
  par.push_back(-sin(houghAngleXY*M_PI/180.)/cos(houghAngleXY*M_PI/180.));
  par.push_back(houghDYZ/sin(houghAngleYZ*M_PI/180.));
  par.push_back(-cos(houghAngleYZ*M_PI/180.)/sin(houghAngleYZ*M_PI/180.));

  double trackMinDistance = CalculateDistance();

  return trackMinDistance;
}

double HoughTrack::FitRestrictedSingleHelp() {
  par.clear();
  Hough2D *houghEvent = new Hough2D(track, detFired, quadFired, 0);
  houghAngleXY = houghEvent->GetMaxThetaXY();
  houghDXY = houghEvent->GetMaxDXY();
  houghAngleYZ = houghEvent->GetMaxThetaYZ();
  houghDYZ = houghEvent->GetMaxDYZ();
  delete houghEvent;

  par.push_back(houghDXY/cos(houghAngleXY*M_PI/180.));
  par.push_back(-sin(houghAngleXY*M_PI/180.)/cos(houghAngleXY*M_PI/180.));
  par.push_back(houghDYZ/sin(houghAngleYZ*M_PI/180.));
  par.push_back(-cos(houghAngleYZ*M_PI/180.)/sin(houghAngleYZ*M_PI/180.));

  double trackMinDistance = CalculateDistance();

  return trackMinDistance;
}

double HoughTrack::FitRestrictedHelp() {
  par.clear();
  Hough2D *houghEvent = new Hough2D(track, detFired, quadFired, 1);
  houghAngleXY = houghEvent->GetMaxThetaXY();
  houghDXY = houghEvent->GetMaxDXY();
  houghAngleYZ = houghEvent->GetMaxThetaYZ();
  houghDYZ = houghEvent->GetMaxDYZ();
  delete houghEvent;

  par.push_back(houghDXY/cos(houghAngleXY*M_PI/180.));
  par.push_back(-sin(houghAngleXY*M_PI/180.)/cos(houghAngleXY*M_PI/180.));
  par.push_back(houghDYZ/sin(houghAngleYZ*M_PI/180.));
  par.push_back(-cos(houghAngleYZ*M_PI/180.)/sin(houghAngleYZ*M_PI/180.));

  double trackMinDistance = CalculateDistance();

  return trackMinDistance;
}

double HoughTrack::CalculateDistance() {
  double dist = 0.;

  for(uint i = 0; i < track.size(); i++) {
    ROOT::Math::XYZVector xp(track[i].xPosition, track[i].yPosition, track[i].height);
    ROOT::Math::XYZVector x0(par[0], 0., par[2]);
    ROOT::Math::XYZVector x1(par[0] + par[1], 1., par[2] + par[3]);
    ROOT::Math::XYZVector u = (x1 - x0).Unit();
    dist += ((xp - x0).Cross(u)).Mag2();
  }

  return dist/static_cast<double>(track.size());
}
