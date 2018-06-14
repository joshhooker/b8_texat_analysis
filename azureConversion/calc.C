#include <cassert>
#include <fstream>
#include <sstream>
#include <vector>

#include "CubicSpline.h"

void calc() {
    std::ifstream in("csReg3.out");
    assert(in.is_open());

    double var1, var2, var3, var4;

    std::vector<double> cmEnergy_, labAngle_, cmAngle_;
    while(in >> var1 >> var2 >> var3 >> var4) {
        cmEnergy_.push_back(var1);
        labAngle_.push_back(var3);
        cmAngle_.push_back(var4);
    }
    in.close();

    std::ifstream inCS("crossSectionReg3.out");
    assert(inCS.is_open());
    std::vector<double> reg3CMEnergy_, reg3CS_, reg3CSErr_;
    while(inCS >> var1 >> var2 >> var3 >> var4) {
        reg3CMEnergy_.push_back(var1);
        reg3CS_.push_back(var2);
        reg3CSErr_.push_back(var3);
    }
    inCS.close();

    CubicSpline splineLabAngle(cmEnergy_, labAngle_);
    CubicSpline splineCMAngle(cmEnergy_, cmAngle_);

    FILE* file = fopen("test.out","w");

    double m1 = 8.;
    double m2 = 1.;

    double x = m1/m2;

    double crossSection = 1.;
    double crossSectionErr = 0.1;

    for(Int_t i = 0; i < reg3CMEnergy_.size(); i++) {
        double cmE = reg3CMEnergy_[i];
        double cs = reg3CS_[i];
        double csErr = reg3CSErr_[i];

        double energyBeamReg = cmE*(m1 + m2)/m1;
        double energyBeamInv = cmE*(m1 + m2)/m2;

        double labAngleLightInv = splineLabAngle(cmE);
        double labAngleHeavyInv = atan(sin(2.*labAngleLightInv*M_PI/180.)/(x - cos(2.*labAngleLightInv*M_PI/180.)))*180./M_PI;

        double cmAngleLight = splineCMAngle(cmE);
        double cmAngleHeavy = 3.1415926 - cmAngleLight;

        double factor1 = m2*m2/((m1 + m2)*(m1 + m2));
        double factor2 = cos(labAngleLightInv*M_PI/180.) + sqrt((m1/m2)*(m1/m2) - sin(labAngleLightInv*M_PI/180.)*sin(labAngleLightInv*M_PI/180.));

        double energyLight = energyBeamReg*factor1*factor2*factor2;

        double labAngleReg = atan(sin(cmAngleLight*M_PI/180.)/(cos(cmAngleLight*M_PI/180.) + m2/m1))*360/(2 * 3.1415926) + 180;

        double crossSectionFactor = 4.*cos(labAngleLightInv*M_PI/180.);
        double crossSectionLab = crossSection/crossSectionFactor;
        double crossSectionLabErr = crossSectionErr/crossSectionFactor;

        fprintf(file,"%f %f %f %f\n", energyLight, labAngleReg, cs, csErr);
    }
    // for(double cmE = 1.; cmE <= 6.; cmE += 0.1) {
    //     double energyBeamReg = cmE*(m1 + m2)/m1;
    //     double energyBeamInv = cmE*(m1 + m2)/m2;

    //     double labAngleLightInv = splineLabAngle(cmE);
    //     double labAngleHeavyInv = atan(sin(2.*labAngleLightInv*M_PI/180.)/(x - cos(2.*labAngleLightInv*M_PI/180.)))*180./M_PI;

    //     double cmAngleLight = splineCMAngle(cmE);
    //     double cmAngleHeavy = 3.1415926 - cmAngleLight;

    //     double factor1 = m2*m2/((m1 + m2)*(m1 + m2));
    //     double factor2 = cos(labAngleLightInv*M_PI/180.) + sqrt((m1/m2)*(m1/m2) - sin(labAngleLightInv*M_PI/180.)*sin(labAngleLightInv*M_PI/180.));

    //     double energyLight = energyBeamReg*factor1*factor2*factor2;

    //     double labAngleReg = atan(sin(cmAngleLight*M_PI/180.)/(cos(cmAngleLight*M_PI/180.) + m2/m1))*360/(2 * 3.1415926) + 180;

    //     double crossSectionFactor = 4.*cos(labAngleLightInv*M_PI/180.);
    //     double crossSectionLab = crossSection/crossSectionFactor;
    //     double crossSectionLabErr = crossSectionErr/crossSectionFactor;

    //     fprintf(file,"%f %f %f %f\n", energyLight, labAngleReg, crossSectionLab, crossSectionLabErr);
    // }

    fflush(file);
    fclose(file);

}