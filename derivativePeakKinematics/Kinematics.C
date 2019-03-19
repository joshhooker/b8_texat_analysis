#include <cmath>

#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TMultiGraph.h>
#include <TRandom3.h>

#include "EnergyLoss.h"

double mmToRow = 1.75; // in mm

TRandom3 *r3 = new TRandom3();

double getDiffRowElastic(EnergyLoss* b8, double initB8Energy, double angleLight, double vertexMM);
int getDiffRowInelastic(EnergyLoss* b8, double initB8Energy, double angleLight);

double getB8EnergyPostVertexElastic(double initB8Energy, double angleLight);
double getB8EnergyPostVertexInelastic(double initB8Energy, double angleLight);

int main() {
    FILE* vertexFile = fopen("centralPeakVertex.out", "w");

    EnergyLoss* b8 = new EnergyLoss("b8_methane.dat");

    double b8StartEnergy = 56.; // in MeV
    double windowToMM = 268.73; // in mm

    TFile* file = new TFile("output.root", "recreate");

    TH2F* hVertexDiff = new TH2F("vertexDiff", "vertexDiff", 128, 0, 128, 128, 0, 128);
    hVertexDiff->GetXaxis()->SetTitle("Vertex Location [rows]"); hVertexDiff->GetXaxis()->CenterTitle();
    hVertexDiff->GetYaxis()->SetTitle("Difference in Vertex Location and Stopping Point [rows]"); hVertexDiff->GetYaxis()->CenterTitle();
    hVertexDiff->GetYaxis()->SetTitleOffset(1.4);
    hVertexDiff->SetStats(false);

    TH2F* hPeakDiff;
    hPeakDiff = new TH2F("peakDiff", "peakDiff", 128, 0, 128, 128, 0, 128);
    hPeakDiff->GetXaxis()->SetTitle("Peak Location [rows]"); hPeakDiff->GetXaxis()->CenterTitle();
    hPeakDiff->GetYaxis()->SetTitle("Difference in Vertex Location and Stopping Point [rows]"); hPeakDiff->GetYaxis()->CenterTitle();
    hPeakDiff->GetYaxis()->SetTitleOffset(1.4);
    hPeakDiff->SetStats(false);

    std::map<int, double> peakBeamEnergy;
    std::map<int, double> peakCMEnergy;
    std::map<int, double> peakVertexMM;
    std::map<int, int> peakTotal;

    for(int i = 0; i < 100; i++) {
        peakBeamEnergy[i] = 0.;
        peakCMEnergy[i] = 0.;
        peakVertexMM[i] = 0.;
        peakTotal[i] = 0;
    }

    for(double energy = b8StartEnergy; energy > 1.; energy -= 0.1) {
        // Calculate vertex location
        double vertexLocation = b8->CalcRange(b8StartEnergy, energy);

        double cmEnergy = energy/9.;

        // Get vertex location in relation to MM
        double vertexMM = vertexLocation - windowToMM;

        // Ignore events with vertex not over MM
        // if(vertexMM < 0) continue;

        // Convert to Row #
        int vertexRowMM = static_cast<int>(vertexMM/mmToRow);

        double peakMM = getDiffRowElastic(b8, energy, 0., vertexMM);
        int peakRowMM = static_cast<int>(peakMM/mmToRow);
        // hVertexDiff->Fill(rowMM, maxRows);
        // hPeakDiff->Fill(rowMM + maxRows, maxRows);

        if(peakRowMM < 1) continue;

        printf("energy: %f, cmEnergy; %f vertexMM: %f, peakMM: %f, peakRowMM: %d\n", energy, cmEnergy, vertexMM, peakMM, peakRowMM);
        // fprintf(vertexFile, "%f %f %f %f %d\n", energy, cmEnergy, vertexMM, peakMM, peakRowMM);

        peakBeamEnergy[peakRowMM] += energy;
        peakCMEnergy[peakRowMM] += cmEnergy;
        peakVertexMM[peakRowMM] += vertexMM;
        peakTotal[peakRowMM]++;
    }

    std::map<int, int>::iterator it;
    for(it = peakTotal.begin(); it != peakTotal.end(); it++) {
        int maxRow = it->first;
        if(peakTotal[maxRow] == 0) continue;
        double avgEnergy = peakBeamEnergy[maxRow]/static_cast<double>(peakTotal[maxRow]);
        double avgCMEnergy = peakCMEnergy[maxRow]/static_cast<double>(peakTotal[maxRow]);
        double avgVertex = peakVertexMM[maxRow]/static_cast<double>(peakTotal[maxRow]);
        fprintf(vertexFile, "%f %f %f %d\n", avgEnergy, avgCMEnergy, avgVertex, maxRow);
    }

    hVertexDiff->Write();
    hPeakDiff->Write();
    file->Close();

    fflush(vertexFile);
    fclose(vertexFile);
}

double getDiffRowElastic(EnergyLoss* b8, double energy, double angleLight, double vertexMM) {
    // Get 8B Energy with proton at zero degrees
    double recoilEnergy = getB8EnergyPostVertexElastic(energy, angleLight);

    double peakMM = vertexMM;
    double maxEnergyDep = -100.;
    while(recoilEnergy > 0.3) {
        double energyRemain = b8->CalcRemainder(recoilEnergy, mmToRow);
        double energyDep = recoilEnergy - energyRemain;
        if(energyDep > maxEnergyDep) {
            peakMM += mmToRow;
            maxEnergyDep = energyDep;
        }
        recoilEnergy = energyRemain;
    }

    return peakMM;
}

double getB8EnergyPostVertexElastic(double initB8Energy, double angleLight) {
    double m1 = 8.;
    double m2 = 1.;
    double kine = (4.*m1*m2/((m1+m2)*(m1+m2)))*cos(angleLight)*cos(angleLight)*initB8Energy;
    return initB8Energy - kine;
}

double getB8EnergyPostVertexInelastic(double initB8Energy, double angleLight) {
    double m1 = 8.;
    double m2 = 1.;
    double m3 = 1.;
    double m4 = 8.;
    double qValue = -0.770; // in MeV
    double E1 = initB8Energy;
    double ET = E1 + qValue;

    double A = m1*m4*(E1/ET)/((m1+m2)*(m3+m4));
    double B = m1*m3*(E1/ET)/((m1+m2)*(m3+m4));
    double C = (m2*m3/((m1+m2)*(m3+m4)))*(1. + m1*qValue/(m2*ET));
    double D = (m2*m4/((m1+m2)*(m3+m4)))*(1. + m1*qValue/(m2*ET));

    std::cout << r3->Uniform() << std::endl;

    double kine = (4.*m1*m2/((m1+m2)*(m1+m2)))*cos(angleLight)*cos(angleLight)*initB8Energy;
    return initB8Energy - kine;
}
