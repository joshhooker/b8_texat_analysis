#include <cassert>
#include <fstream>

void calcExcitationEnergy() {
  std::ifstream in("testCS.out");
  assert(in.is_open());

  std::vector<Double_t> binCenter, binContent, binError;
  double var1, var2, var3;
  while(in >> var1 >> var2 >> var3) {
    binCenter.push_back(var1+1.3);
    binContent.push_back(var2*1000.);
    binError.push_back(var3*1000.);
  }
  in.close();

  size_t n = binCenter.size();
  Double_t x[n], y[n], ex[n], ey[n];
  for(size_t i=0; i<n; i++) {
    x[i] = binCenter[i]+0.2;
    y[i] = binContent[i];
    ex[i] = 0.;
    ey[i] = binError[i];
  }

  TGraphErrors* gr = new TGraphErrors(n, x, y, ex, ey);
  gr->Draw("AP");
}