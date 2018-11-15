#ifndef ENERGY_LOSS_H
#define ENERGY_LOSS_H

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

#include "CubicSpline.h"

#define ASSERT_WITH_MESSAGE(condition, message)                                \
  do {                                                                         \
    if (!(condition)) {                                                        \
      printf((message));                                                       \
    }                                                                          \
    assert((condition));                                                       \
  } while (false)

class EnergyLoss {
public:
  EnergyLoss(const char *);
  EnergyLoss(const char *, double);

  void ReadInitParams();
  double CalcRemainder(double, double);
  double AddBack(double, double);
  double CalcRange(double, double);
  void CalcRemainderError(double);
  void AddBackError(double);
  void AddBackHigh(double);

  friend double GetdEdx(const EnergyLoss&, double);

  void UseGL16();
  void UseGL32();
  void UseGL64();
  void UseGL128();
  void UseGL256();
  void UseGL512();
  void UseGL1024();

private:
  std::vector<double> energy_;
  std::vector<double> dEdx_;
  CubicSpline energySpline;

  double x16[8], w16[8], x32[16], w32[16], x64[32], w64[32];
  double x128[64], w128[64], x256[128], w256[128];
  double x512[256], w512[256], x1024[512], w1024[512];

  double CalcRemainderErr;
  double AddBackErr;
  double AddBackHighPoint;

  bool GL16, GL32, GL64, GL128, GL256, GL512, GL1024;

  double CalcRangeGL16(double, double);
  double CalcRangeGL32(double, double);
  double CalcRangeGL64(double, double);
  double CalcRangeGL128(double, double);
  double CalcRangeGL256(double, double);
  double CalcRangeGL512(double, double);
  double CalcRangeGL1024(double, double);

  // Variables to read SRIM File
  const int MAX_CHARS_PER_LINE = 1024000;
  const char *const DELIMITER = " ";
  bool beforeMultiply = true;
  std::vector<char *> token;
  double stoppingConversionPower;
  std::vector<double> stoppingConversionPower_;
  std::vector<std::string> stoppingConversionUnit1_;
  std::vector<std::string> stoppingConversionUnit2_;

};

#endif
