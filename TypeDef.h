#ifndef TYPEDEF_H
#define TYPEDEF_H

#include <TROOT.h>

typedef struct siDetect {
  int detect;
  int quad;
  int channel;
  double energy;
  double time;
  double cwt;
} siDetect;

typedef struct csiDetect {
  int detect;
  double energy;
  double time;
} csiDetect;

typedef struct protonDetect{
  int det;
  int quad;
  int siChannel;
  double siEnergy;
  double siEnergyCal;
  double siTime;
  double csiEnergy;
  double csiEnergyCal;
  double csiTime;
  double totalEnergy;
  bool punchthrough;
} protonDetect;

typedef struct mmCenter {
  int column;
  int row;
  double energy;
  double time;
  double cwt;
  int timeBucket;
} mmCenter;

typedef struct mmTrack {
  int row;
  double xPosition;
  double yPosition;
  double time;
  double energy;
  double height;
  int timeBucket;
  int total;
} mmTrack;

typedef struct mmChainStrip {
  int row;
  double energy;
  double time;
  int timeBucket;
} mmStripChain;

typedef struct xy {
  double x;
  double y;
} xy;

typedef struct yz {
  double y;
  double z;
} yz;

typedef struct centerDeriv {
  int row;
  double deriv;
} centerDeriv;

#endif //TYPEDEF_H
