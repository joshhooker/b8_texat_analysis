#include "EnergyLoss.h"

EnergyLoss::EnergyLoss(const char *srimFile) {
  std::ifstream inFile(srimFile);
  ASSERT_WITH_MESSAGE(inFile.is_open(), "Cannot find SRIM file!\n");

  while(!inFile.eof()) {
    char buf[MAX_CHARS_PER_LINE];
    inFile.getline(buf, MAX_CHARS_PER_LINE);
    char *segment = strtok(buf, DELIMITER);
    token.push_back(segment);

    if(token[0]) { //check if the line is good
      if(atof(token[0]) && beforeMultiply) { // check if number
        double energy = atof(token[0]);
        std::string energyUnit = strtok(0, DELIMITER);
        double dEElec = atof(strtok(0, DELIMITER));
        double dENuclear  = atof(strtok(0, DELIMITER));

        if(energyUnit == "eV") {
          energy_.push_back(energy*1.e-6);
        } else if(energyUnit == "keV") {
          energy_.push_back(energy*1.e-3);
        } else if(energyUnit == "MeV") {
          energy_.push_back(energy);
        } else if(energyUnit == "GeV") {
          energy_.push_back(energy*1.e3);
        }
        dEdx_.push_back(dEElec+dENuclear);
      }

      if(!strcmp(token[0], "Multiply")) {
        beforeMultiply = false;
      }

      if(atof(token[0]) && !beforeMultiply) {
        stoppingConversionPower_.push_back(atof(token[0]));
        stoppingConversionUnit1_.push_back(strtok(0, DELIMITER));
        strtok(0, DELIMITER);
        stoppingConversionUnit2_.push_back(strtok(0, DELIMITER));
      }
    }

    if(!token.empty()) {
      token.clear(); // clear token vector
    }
  }

  for(size_t i=0; i<stoppingConversionPower_.size(); i++) {
    if(stoppingConversionUnit1_[i]=="MeV" && stoppingConversionUnit2_[i]=="mm") {
      stoppingConversionPower = stoppingConversionPower_[i];
    }
  }
  std::transform(dEdx_.begin(), dEdx_.end(), dEdx_.begin(),
                 std::bind1st(std::multiplies<double>(), stoppingConversionPower));
  energySpline.SetPoints(energy_, dEdx_);

  GL16  = GL32 = GL64 = GL128 = GL512 = GL1024 = false;
  GL256 = true;

  CalcRemainderErr = 1.e-8;
  AddBackErr = 1.e-8;
  AddBackHighPoint = 250.;

  ReadInitParams();
}

EnergyLoss::EnergyLoss(const char *srimFile, double stoppingPower) {
  std::ifstream inFile(srimFile);
  ASSERT_WITH_MESSAGE(inFile.is_open(), "Cannot find SRIM file!\n");

  while(!inFile.eof()) {
    char buf[MAX_CHARS_PER_LINE];
    inFile.getline(buf, MAX_CHARS_PER_LINE);
    char *segment = strtok(buf, DELIMITER);
    token.push_back(segment);

    if(token[0]) { //check if the line is good
      if(atof(token[0]) && beforeMultiply) { // check if number
        double energy = atof(token[0]);
        std::string energyUnit = strtok(0, DELIMITER);
        double dEElec = atof(strtok(0, DELIMITER));
        double dENuclear = atof(strtok(0, DELIMITER));

        if(energyUnit == "eV") {
          energy_.push_back(energy*1.e-6);
        } else if(energyUnit == "keV") {
          energy_.push_back(energy*1.e-3);
        } else if(energyUnit == "MeV") {
          energy_.push_back(energy);
        } else if(energyUnit == "GeV") {
          energy_.push_back(energy*1.e3);
        }
        dEdx_.push_back(dEElec+dENuclear);
      }
    }

    if(!token.empty()) {
      token.clear(); // clear token vector
    }
  }

  std::transform(dEdx_.begin(), dEdx_.end(), dEdx_.begin(),
                 std::bind1st(std::multiplies<double>(), stoppingPower));
  energySpline.SetPoints(energy_, dEdx_);

  GL16  = GL32 = GL64 = GL128 = GL512 = GL1024 = false;
  GL256 = true;

  CalcRemainderErr = 1.e-8;
  AddBackErr = 1.e-8;
  AddBackHighPoint = 250.;

  ReadInitParams();
}

double EnergyLoss::CalcRemainder(double initialEnergy, double distance) {
  if(distance == 0.) return initialEnergy;
  if(initialEnergy < CalcRemainderErr) return 0.;

  double maxRange = CalcRange(initialEnergy, 0.);

  if(distance > maxRange) return 0.;

  double lowEnergy = 0.;
  double highEnergy = initialEnergy;
  double guessEnergy = (highEnergy+lowEnergy)/2.0;

  double range = CalcRange(initialEnergy, guessEnergy);
  while(fabs(range-distance) > CalcRemainderErr) {
    if(range > distance) {
      lowEnergy = guessEnergy;
      guessEnergy = (lowEnergy+highEnergy)/2.0;
    } else {
      highEnergy = guessEnergy;
      guessEnergy = (lowEnergy+highEnergy)/2.0;
    }
    range = CalcRange(initialEnergy, guessEnergy);
  }
  return guessEnergy;
}

double EnergyLoss::AddBack(double finalEnergy, double distance) {
  if(distance == 0.) return finalEnergy;

  double lowEnergy = finalEnergy;
  double highEnergy = AddBackHighPoint;
  double guessEnergy = (highEnergy+lowEnergy)/2.0;

  double range = CalcRange(guessEnergy, finalEnergy);
  while(fabs(range-distance) > AddBackErr) {
    if(range < distance) {
      lowEnergy = guessEnergy;
      guessEnergy = (lowEnergy+highEnergy)/2.0;
    } else {
      highEnergy = guessEnergy;
      guessEnergy = (lowEnergy+highEnergy)/2.0;
    }
    range = CalcRange(guessEnergy, finalEnergy);

    if(guessEnergy > (AddBackHighPoint-0.1)) {
      printf("Error: EnergyLoss::AddBack above starting high guess!\n");
      return guessEnergy;
    }
  }
  return guessEnergy;
}

double EnergyLoss::CalcRange(double initialEnergy, double remainder) {

  if(initialEnergy == remainder) return 0;
  double distance;

  if(GL16) distance = CalcRangeGL16(remainder, initialEnergy);
  else if(GL32) distance = CalcRangeGL32(remainder, initialEnergy);
  else if(GL64) distance = CalcRangeGL64(remainder, initialEnergy);
  else if(GL128) distance = CalcRangeGL128(remainder, initialEnergy);
  else if(GL256) distance = CalcRangeGL256(remainder, initialEnergy);
  else if(GL512) distance = CalcRangeGL512(remainder, initialEnergy);
  else if(GL1024) distance = CalcRangeGL1024(remainder, initialEnergy);
  else distance = CalcRangeGL256(remainder, initialEnergy);

  return distance;
}

void EnergyLoss::CalcRemainderError(double err) {
  CalcRemainderErr = err;
}

void EnergyLoss::AddBackError(double err) {
  AddBackErr = err;
}

void EnergyLoss::AddBackHigh(double high) {
  AddBackHighPoint = high;
}

void EnergyLoss::UseGL16() {
  GL32 = GL64 = GL128 = GL256 = GL512 = GL1024 = false;
  GL16 = true;
}

void EnergyLoss::UseGL32() {
  GL16 = GL64 = GL128 = GL256 = GL512 = GL1024 = false;
  GL32 = true;
}

void EnergyLoss::UseGL64() {
  GL16 = GL32 = GL128 = GL256 = GL512 = GL1024 = false;
  GL64 = true;
}

void EnergyLoss::UseGL128() {
  GL16 = GL32 = GL64 = GL256 = GL512 = GL1024 = false;
  GL128 = true;
}

void EnergyLoss::UseGL256() {
  GL16 = GL32 = GL64 = GL128 = GL512 = GL1024 = false;
  GL256 = true;
}

void EnergyLoss::UseGL512() {
  GL16 = GL32 = GL64 = GL128 = GL256 = GL1024 = false;
  GL512 = true;
}

void EnergyLoss::UseGL1024() {
  GL16 = GL32 = GL64 = GL128 = GL256 = GL512 = false;
  GL1024 = true;
}

double EnergyLoss::CalcRangeGL16(double a, double b) {
  double t0 = (a+b)/2.0;
  double dt = (b-a)/2.0;
  double result = 0.;
  for(int i=7; i>-1; i--) {
    result += w16[i]*(1./energySpline(t0 + dt*(-x16[i])));
  }
  for(int i=0; i<8; i++) {
    result += w16[i]*(1./energySpline(t0 + dt*x16[i]));
  }
  return result*dt;
}

double EnergyLoss::CalcRangeGL32(double a, double b) {
  double t0 = (a+b)/2.0;
  double dt = (b-a)/2.0;
  double result = 0.;
  for(int i=15; i>-1; i--) {
    result += w32[i]*(1./energySpline(t0 + dt*(-x32[i])));
  }
  for(int i=0; i<16; i++) {
    result += w32[i]*(1./energySpline(t0 + dt*x32[i]));
  }
  return result*dt;
}

double EnergyLoss::CalcRangeGL64(double a, double b) {
  double t0 = (a+b)/2.0;
  double dt = (b-a)/2.0;
  double result = 0.;
  for (int i=31; i>-1; i--) {
    result += w64[i]*(1./energySpline(t0 + dt*(-x64[i])));
  }
  for (int i=0; i<32; i++) {
    result += w64[i]*(1./energySpline(t0 + dt*x64[i]));
  }
  return result*dt;
}

double EnergyLoss::CalcRangeGL128(double a, double b) {
  double t0 = (a+b)/2.0;
  double dt = (b-a)/2.0;
  double result = 0.;
  for (int i=63; i>-1; i--) {
    result += w128[i]*(1./energySpline(t0 + dt*(-x128[i])));
  }
  for (int i=0; i<64; i++) {
    result += w128[i]*(1./energySpline(t0 + dt*x128[i]));
  }
  return result*dt;
}

double EnergyLoss::CalcRangeGL256(double a, double b) {
  double t0 = (a+b)/2.0;
  double dt = (b-a)/2.0;
  double result = 0.;
  for (int i=127; i>-1; i--) {
    result += w256[i]*(1./energySpline(t0 + dt*(-x256[i])));
  }
  for (int i=0; i<128; i++) {
    result += w256[i]*(1./energySpline(t0 + dt*x256[i]));
  }
  return result*dt;
}

double EnergyLoss::CalcRangeGL512(double a, double b) {
  double t0 = (a+b)/2.0;
  double dt = (b-a)/2.0;
  double result = 0.;
  for (int i=255; i>-1; i--) {
    result += w512[i]*(1./energySpline(t0 + dt*(-x512[i])));
  }
  for (int i=0; i<256; i++) {
    result += w512[i]*(1./energySpline(t0 + dt*x512[i]));
  }
  return result*dt;
}

double EnergyLoss::CalcRangeGL1024(double a, double b) {
  double t0 = (a+b)/2.0;
  double dt = (b-a)/2.0;
  double result = 0.;
  for (int i=511; i>-1; i--) {
    result += w1024[i]*(1./energySpline(t0 + dt*(-x1024[i])));
  }
  for (int i=0; i<512; i++) {
    result += w1024[i]*(1./energySpline(t0 + dt*x1024[i]));
  }
  return result * dt;
}

void EnergyLoss::ReadInitParams() {
  std::ifstream inFile("EnergyLoss.dat", std::ios::in | std::ifstream::binary);
  inFile.read(reinterpret_cast<char *> (x16), sizeof(x16));
  inFile.read(reinterpret_cast<char *> (w16), sizeof(w16));
  inFile.read(reinterpret_cast<char *> (x32), sizeof(x32));
  inFile.read(reinterpret_cast<char *> (w32), sizeof(w32));
  inFile.read(reinterpret_cast<char *> (x64), sizeof(x64));
  inFile.read(reinterpret_cast<char *> (w64), sizeof(w64));
  inFile.read(reinterpret_cast<char *> (x128), sizeof(x128));
  inFile.read(reinterpret_cast<char *> (w128), sizeof(w128));
  inFile.read(reinterpret_cast<char *> (x256), sizeof(x256));
  inFile.read(reinterpret_cast<char *> (w256), sizeof(w256));
  inFile.read(reinterpret_cast<char *> (x512), sizeof(x512));
  inFile.read(reinterpret_cast<char *> (w512), sizeof(w512));
  inFile.read(reinterpret_cast<char *> (x1024), sizeof(x1024));
  inFile.read(reinterpret_cast<char *> (w1024), sizeof(w1024));
  inFile.close();
}

double GetdEdx(const EnergyLoss& elClass, double energy) {
  return elClass.energySpline(energy);
}
