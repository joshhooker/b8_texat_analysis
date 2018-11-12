#ifndef TYPEDEF_H
#define TYPEDEF_H

#include <TROOT.h>

typedef struct siDetect {
  Int_t detect;
  Int_t quad;
  Int_t channel;
  Double_t energy;
  Double_t time;
  Double_t cwt;
} siDetect;

typedef struct csiDetect {
  Int_t detect;
  Double_t energy;
  Double_t time;
} csiDetect;

typedef struct mmCenter {
  Int_t column;
  Int_t row;
  Double_t energy;
  Double_t time;
  Double_t cwt;
} mmCenter;

typedef struct mmTrack {
  Int_t row;
  Double_t xPosition;
  Double_t yPosition;
  Double_t time;
  Double_t energy;
  Double_t height;
  Int_t total;
} mmTrack;

typedef struct mmChainStrip {
  Int_t row;
  Double_t energy;
  Double_t time;
} mmStripChain;

typedef struct xy {
  Double_t x;
  Double_t y;
} xy;

typedef struct yz {
  Double_t y;
  Double_t z;
} yz;

typedef struct centerDeriv {
  Int_t row;
  Double_t deriv;
} centerDeriv;

#endif //TYPEDEF_H
