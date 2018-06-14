#include <fstream>
#include <sstream>

#include "CubicSpline.h"
#include "EnergyLoss.h"

void calcSolidAngle() {
    FILE* file = fopen("csReg3.out", "w");

    Float_t m1 = 8.;
    Float_t m2 = 1.;

    Float_t beamInitial = 56.; // In MeV, after havar window
    Float_t windowToSi = 544.;

    EnergyLoss *boronMethane = new EnergyLoss("b8_methane.dat");
    EnergyLoss *protonMethane = new EnergyLoss("proton_methane.dat");

    for(Float_t cmEnergy = 0.8; cmEnergy < 7.; cmEnergy += 0.05) {
        Float_t beamEnergy = cmEnergy*(m1 + m2)/m2;

        if(beamEnergy > beamInitial) continue;

        TH1F *h = new TH1F(Form("%f_lab_ik", cmEnergy), Form("%f_lab_ik", cmEnergy), 360, 0, 180);
        TH1F *h2 = new TH1F(Form("%f_cm_fk", cmEnergy), Form("%f_cm_fk", cmEnergy), 360, 0, 180);

        Float_t depth = boronMethane->CalcRange(beamInitial, beamEnergy);
        Float_t z = windowToSi - depth;

        std::cout << cmEnergy << '\t' << beamEnergy << '\t' << depth << '\t' << z << std::endl;

        Float_t angle;
        Float_t sum = 0;

        Float_t elementSize = 0.1;

        for(Float_t dx = 99; dx <= 149; dx += elementSize) {
            for(Float_t dy = 10.8; dy <= 60.8; dy += elementSize) {
                Float_t dr = sqrt(dx*dx + dy*dy + z*z);
                angle = acos(z/dr)/3.14159*180;
                Float_t protonEnergy = 4.*m1/(m1 + m2) * z*z/(dr*dr)*cmEnergy;
                Float_t range = protonMethane->CalcRange(protonEnergy, 0);
                if(range >= dr) {
                    sum += elementSize*elementSize*z/(dr*dr*dr);
                    h->Fill(angle);
                    h2->Fill(180. - 2.*angle);
                }
            }
        }
        sum *= 3.;

        if(sum == 0.) continue;

        Float_t aFactor = 4.*cos(h->GetMean()*3.14159/180.);
        Float_t change_bin_content = sum*aFactor;

        printf("%f %f %f %f\n", cmEnergy, change_bin_content, h->GetMean(), h2->GetMean());
        fprintf(file, "%f %f %f %f\n", cmEnergy, change_bin_content, h->GetMean(), h2->GetMean());

        delete h;
        delete h2;
    }

    fflush(file);
    fclose(file);
}