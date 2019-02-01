void plotReg3() {
    gROOT->SetStyle("Plain");

    double var1, var2, var3, var4, var5, var6, var7, var8, var9;

    std::vector<double> cm_05_25_15, exe_05_25_15;
    std::vector<double> fit_05_25_15;
    ifstream in_05_25_15("9C_Reg3_0.5-_2.5-_1.5-.out");
    while(in_05_25_15 >> var1 >> var2 >> var3 >> var4 >> var5 >> var6 >> var7 >> var8 >> var9) {
        cm_05_25_15.push_back(var1);
        exe_05_25_15.push_back(var2);
        fit_05_25_15.push_back(var4);
    }
    in_05_25_15.close();

    std::vector<double> cm_05_25_15_35, exe_05_25_15_35;
    std::vector<double> fit_05_25_15_35, cs_05_25_15_35, cserr_05_25_15_35;
    ifstream in_05_25_15_35("9C_Reg3_0.5-_2.5-_1.5+_3.5-.out");
    while(in_05_25_15_35 >> var1 >> var2 >> var3 >> var4 >> var5 >> var6 >> var7 >> var8 >> var9) {
        cm_05_25_15_35.push_back(var1);
        exe_05_25_15_35.push_back(var2);
        fit_05_25_15_35.push_back(var4);
        cs_05_25_15_35.push_back(var6);
        cserr_05_25_15_35.push_back(var7);
    }
    in_05_25_15_35.close();

    std::vector<double> cm_05_25_25_35, exe_05_25_25_35;
    std::vector<double> fit_05_25_25_35;
    ifstream in_05_25_25_35("9C_Reg3_0.5-_2.5-_2.5+_3.5-.out");
    while(in_05_25_25_35 >> var1 >> var2 >> var3 >> var4 >> var5 >> var6 >> var7 >> var8 >> var9) {
        cm_05_25_25_35.push_back(var1);
        exe_05_25_25_35.push_back(var2);
        fit_05_25_25_35.push_back(var4);
    }
    in_05_25_25_35.close();

    std::vector<double> cm_05_25_25, exe_05_25_25;
    std::vector<double> fit_05_25_25;
    ifstream in_05_25_25("9C_Reg3_0.5-_2.5-_2.5+.out");
    while(in_05_25_25 >> var1 >> var2 >> var3 >> var4 >> var5 >> var6 >> var7 >> var8 >> var9) {
        cm_05_25_25.push_back(var1);
        exe_05_25_25.push_back(var2);
        fit_05_25_25.push_back(var4);
    }
    in_05_25_25.close();

    std::vector<double> cm_05_25_35, exe_05_25_35;
    std::vector<double> fit_05_25_35;
    ifstream in_05_25_35("9C_Reg3_0.5-_2.5-_3.5-.out");
    while(in_05_25_35 >> var1 >> var2 >> var3 >> var4 >> var5 >> var6 >> var7 >> var8 >> var9) {
        cm_05_25_35.push_back(var1);
        exe_05_25_35.push_back(var2);
        fit_05_25_35.push_back(var4);
    }
    in_05_25_35.close();

    std::vector<double> cm_05_25, exe_05_25;
    std::vector<double> fit_05_25;
    ifstream in_05_25("9C_Reg3_0.5-_2.5-.out");
    while(in_05_25 >> var1 >> var2 >> var3 >> var4 >> var5 >> var6 >> var7 >> var8 >> var9) {
        cm_05_25.push_back(var1);
        exe_05_25.push_back(var2);
        fit_05_25.push_back(var4);
    }
    in_05_25.close();

    TCanvas* c = new TCanvas("c");
    int size = cm_05_25_15_35.size();

    TGraphErrors *region1 = new TGraphErrors();

    TGraph *gr_05_25_05_35 = new TGraph();
    TGraph *gr_05_25_15 = new TGraph();
    TGraph *gr_05_25_15_35 = new TGraph();
    TGraph *gr_05_25_25_35 = new TGraph();
    TGraph *gr_05_25_25 = new TGraph();
    TGraph *gr_05_25_35 = new TGraph();
    TGraph *gr_05_25 = new TGraph();

    for(int i = 0; i < size; i++) {
        region1->SetPoint(i, cm_05_25_15_35[i] + 1.3, cs_05_25_15_35[i]*1000.);
        region1->SetPointError(i, 0, cserr_05_25_15_35[i]*1000.);

        gr_05_25_15->SetPoint(i, cm_05_25_15[i] + 1.3, fit_05_25_15[i]*1000.);

        gr_05_25_15_35->SetPoint(i, cm_05_25_15_35[i] + 1.3, fit_05_25_15_35[i]*1000.);

        gr_05_25_25_35->SetPoint(i, cm_05_25_25_35[i] + 1.3, fit_05_25_25_35[i]*1000.);

        gr_05_25_25->SetPoint(i, cm_05_25_25[i] + 1.3, fit_05_25_25[i]*1000.);

        gr_05_25_35->SetPoint(i, cm_05_25_35[i] + 1.3, fit_05_25_35[i]*1000.);

        gr_05_25->SetPoint(i, cm_05_25[i] + 1.3, fit_05_25[i]*1000.);
    }

    region1->SetMarkerStyle(20);
    region1->SetTitle("; Excitation Energy [MeV]; Cross Section [mb/sr]");
    region1->GetXaxis()->CenterTitle();
    region1->GetYaxis()->CenterTitle();
    region1->SetMinimum(0);
    region1->SetMaximum(160);
    region1->Draw("APE");

    TLatex latex;
    latex.DrawLatex(2.5, 140., "105^{o}-145^{o}");

    // c->Print("9C_Region3_CS.pdf");

    //graph 1

    // 05_25_15
    gr_05_25_15->SetLineStyle(1);  gr_05_25_15->SetLineColor(1); gr_05_25_15->SetLineWidth(4);
    gr_05_25_15->Draw("PCsame");

    // 05_25_25
    gr_05_25_25->SetLineStyle(9); gr_05_25_25->SetLineColor(3); gr_05_25_25->SetLineWidth(4);
    gr_05_25_25->Draw("PCsame");

    // 05_25
    gr_05_25->SetLineStyle(10); gr_05_25->SetLineColor(4); gr_05_25->SetLineWidth(4);
    gr_05_25->Draw("PCsame");

    // 05_25_25_35
    gr_05_25_25_35->SetLineStyle(2); gr_05_25_25_35->SetLineColor(2); gr_05_25_25_35->SetLineWidth(5);
    gr_05_25_25_35->Draw("PCsame");

    c->Print("9C_Region3_CS_graph1.pdf");

    // graph 2

    // // 05_25_25_35
    // gr_05_25_25_35->SetLineStyle(2); gr_05_25_25_35->SetLineColor(2); gr_05_25_25_35->SetLineWidth(5);
    // gr_05_25_25_35->Draw("PCsame");

    // // 05_25_15_35
    // gr_05_25_15_35->SetLineStyle(9); gr_05_25_15_35->SetLineColor(3); gr_05_25_15_35->SetLineWidth(4);
    // gr_05_25_15_35->Draw("PCsame");

    // c->Print("9C_Region3_CS_graph2.pdf");
}