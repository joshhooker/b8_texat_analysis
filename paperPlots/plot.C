void plot() {
    gROOT->SetStyle("Plain");

    double var1, var2, var3, var4, var5, var6, var7, var8, var9;

    //////////////
    // Region 1 //
    //////////////
    std::vector<double> cm_reg1_05_25_15, exe_reg1_05_25_15;
    std::vector<double> fit_reg1_05_25_15;
    ifstream in_reg1_05_25_15("9C_Reg1_0.5-_2.5-_1.5-.out");
    while(in_reg1_05_25_15 >> var1 >> var2 >> var3 >> var4 >> var5 >> var6 >> var7 >> var8 >> var9) {
        cm_reg1_05_25_15.push_back(var1);
        exe_reg1_05_25_15.push_back(var2);
        fit_reg1_05_25_15.push_back(var4);
    }
    in_reg1_05_25_15.close();

    std::vector<double> cm_reg1_05_25_15_35, exe_reg1_05_25_15_35;
    std::vector<double> fit_reg1_05_25_15_35, cs_reg1_05_25_15_35, cserr_reg1_05_25_15_35;
    ifstream in_reg1_05_25_15_35("9C_Reg1_0.5-_2.5-_1.5+_3.5-.out");
    while(in_reg1_05_25_15_35 >> var1 >> var2 >> var3 >> var4 >> var5 >> var6 >> var7 >> var8 >> var9) {
        cm_reg1_05_25_15_35.push_back(var1);
        exe_reg1_05_25_15_35.push_back(var2);
        fit_reg1_05_25_15_35.push_back(var4);
        cs_reg1_05_25_15_35.push_back(var6);
        cserr_reg1_05_25_15_35.push_back(var7);
    }
    in_reg1_05_25_15_35.close();

    std::vector<double> cm_reg1_05_25_25_35, exe_reg1_05_25_25_35;
    std::vector<double> fit_reg1_05_25_25_35;
    ifstream in_reg1_05_25_25_35("9C_Reg1_0.5-_2.5-_2.5+_3.5-.out");
    while(in_reg1_05_25_25_35 >> var1 >> var2 >> var3 >> var4 >> var5 >> var6 >> var7 >> var8 >> var9) {
        cm_reg1_05_25_25_35.push_back(var1);
        exe_reg1_05_25_25_35.push_back(var2);
        fit_reg1_05_25_25_35.push_back(var4);
    }
    in_reg1_05_25_25_35.close();

    std::vector<double> cm_reg1_05_25_25, exe_reg1_05_25_25;
    std::vector<double> fit_reg1_05_25_25;
    ifstream in_reg1_05_25_25("9C_Reg1_0.5-_2.5-_2.5+.out");
    while(in_reg1_05_25_25 >> var1 >> var2 >> var3 >> var4 >> var5 >> var6 >> var7 >> var8 >> var9) {
        cm_reg1_05_25_25.push_back(var1);
        exe_reg1_05_25_25.push_back(var2);
        fit_reg1_05_25_25.push_back(var4);
    }
    in_reg1_05_25_25.close();

    std::vector<double> cm_reg1_05_25_35, exe_reg1_05_25_35;
    std::vector<double> fit_reg1_05_25_35;
    ifstream in_reg1_05_25_35("9C_Reg1_0.5-_2.5-_3.5-.out");
    while(in_reg1_05_25_35 >> var1 >> var2 >> var3 >> var4 >> var5 >> var6 >> var7 >> var8 >> var9) {
        cm_reg1_05_25_35.push_back(var1);
        exe_reg1_05_25_35.push_back(var2);
        fit_reg1_05_25_35.push_back(var4);
    }
    in_reg1_05_25_35.close();

    std::vector<double> cm_reg1_05_25, exe_reg1_05_25;
    std::vector<double> fit_reg1_05_25;
    ifstream in_reg1_05_25("9C_Reg1_0.5-_2.5-.out");
    while(in_reg1_05_25 >> var1 >> var2 >> var3 >> var4 >> var5 >> var6 >> var7 >> var8 >> var9) {
        cm_reg1_05_25.push_back(var1);
        exe_reg1_05_25.push_back(var2);
        fit_reg1_05_25.push_back(var4);
    }
    in_reg1_05_25.close();

    int size_reg1 = cm_reg1_05_25_15_35.size();

    TGraphErrors *region1 = new TGraphErrors();

    TGraph *gr_reg1_05_25_05_35 = new TGraph();
    TGraph *gr_reg1_05_25_15 = new TGraph();
    TGraph *gr_reg1_05_25_15_35 = new TGraph();
    TGraph *gr_reg1_05_25_25_35 = new TGraph();
    TGraph *gr_reg1_05_25_25 = new TGraph();
    TGraph *gr_reg1_05_25_35 = new TGraph();
    TGraph *gr_reg1_05_25 = new TGraph();

    for(int i = 0; i < size_reg1; i++) {
        region1->SetPoint(i, cm_reg1_05_25_15_35[i] + 1.3, cs_reg1_05_25_15_35[i]*1000.);
        region1->SetPointError(i, 0, cserr_reg1_05_25_15_35[i]*1000.);

        gr_reg1_05_25_15->SetPoint(i, cm_reg1_05_25_15[i] + 1.3, fit_reg1_05_25_15[i]*1000.);

        gr_reg1_05_25_15_35->SetPoint(i, cm_reg1_05_25_15_35[i] + 1.3, fit_reg1_05_25_15_35[i]*1000.);

        gr_reg1_05_25_25_35->SetPoint(i, cm_reg1_05_25_25_35[i] + 1.3, fit_reg1_05_25_25_35[i]*1000.);

        gr_reg1_05_25_25->SetPoint(i, cm_reg1_05_25_25[i] + 1.3, fit_reg1_05_25_25[i]*1000.);

        gr_reg1_05_25_35->SetPoint(i, cm_reg1_05_25_35[i] + 1.3, fit_reg1_05_25_35[i]*1000.);

        gr_reg1_05_25->SetPoint(i, cm_reg1_05_25[i] + 1.3, fit_reg1_05_25[i]*1000.);
    }

    region1->SetMarkerStyle(20);
    region1->SetTitle("; Excitation Energy [MeV]; Cross Section [mb/sr]");
    region1->GetXaxis()->CenterTitle();
    region1->GetYaxis()->CenterTitle();
    region1->SetMinimum(10);
    region1->SetMaximum(150);
    region1->GetXaxis()->SetLimits(1.5, 6.5);
    // region1->Draw("APE");

    //////////////
    // Region 3 //
    //////////////

    std::vector<double> cm_reg3_05_25_15, exe_reg3_05_25_15;
    std::vector<double> fit_reg3_05_25_15;
    ifstream in_reg3_05_25_15("9C_Reg3_0.5-_2.5-_1.5-.out");
    while(in_reg3_05_25_15 >> var1 >> var2 >> var3 >> var4 >> var5 >> var6 >> var7 >> var8 >> var9) {
        cm_reg3_05_25_15.push_back(var1);
        exe_reg3_05_25_15.push_back(var2);
        fit_reg3_05_25_15.push_back(var4);
    }
    in_reg3_05_25_15.close();

    std::vector<double> cm_reg3_05_25_15_35, exe_reg3_05_25_15_35;
    std::vector<double> fit_reg3_05_25_15_35, cs_reg3_05_25_15_35, cserr_reg3_05_25_15_35;
    ifstream in_reg3_05_25_15_35("9C_Reg3_0.5-_2.5-_1.5+_3.5-.out");
    while(in_reg3_05_25_15_35 >> var1 >> var2 >> var3 >> var4 >> var5 >> var6 >> var7 >> var8 >> var9) {
        cm_reg3_05_25_15_35.push_back(var1);
        exe_reg3_05_25_15_35.push_back(var2);
        fit_reg3_05_25_15_35.push_back(var4);
        cs_reg3_05_25_15_35.push_back(var6);
        cserr_reg3_05_25_15_35.push_back(var7);
    }
    in_reg3_05_25_15_35.close();

    std::vector<double> cm_reg3_05_25_25_35, exe_reg3_05_25_25_35;
    std::vector<double> fit_reg3_05_25_25_35;
    ifstream in_reg3_05_25_25_35("9C_Reg3_0.5-_2.5-_2.5+_3.5-.out");
    while(in_reg3_05_25_25_35 >> var1 >> var2 >> var3 >> var4 >> var5 >> var6 >> var7 >> var8 >> var9) {
        cm_reg3_05_25_25_35.push_back(var1);
        exe_reg3_05_25_25_35.push_back(var2);
        fit_reg3_05_25_25_35.push_back(var4);
    }
    in_reg3_05_25_25_35.close();

    std::vector<double> cm_reg3_05_25_25, exe_reg3_05_25_25;
    std::vector<double> fit_reg3_05_25_25;
    ifstream in_reg3_05_25_25("9C_Reg3_0.5-_2.5-_2.5+.out");
    while(in_reg3_05_25_25 >> var1 >> var2 >> var3 >> var4 >> var5 >> var6 >> var7 >> var8 >> var9) {
        cm_reg3_05_25_25.push_back(var1);
        exe_reg3_05_25_25.push_back(var2);
        fit_reg3_05_25_25.push_back(var4);
    }
    in_reg3_05_25_25.close();

    std::vector<double> cm_reg3_05_25_35, exe_reg3_05_25_35;
    std::vector<double> fit_reg3_05_25_35;
    ifstream in_reg3_05_25_35("9C_Reg3_0.5-_2.5-_3.5-.out");
    while(in_reg3_05_25_35 >> var1 >> var2 >> var3 >> var4 >> var5 >> var6 >> var7 >> var8 >> var9) {
        cm_reg3_05_25_35.push_back(var1);
        exe_reg3_05_25_35.push_back(var2);
        fit_reg3_05_25_35.push_back(var4);
    }
    in_reg3_05_25_35.close();

    std::vector<double> cm_reg3_05_25, exe_reg3_05_25;
    std::vector<double> fit_reg3_05_25;
    ifstream in_reg3_05_25("9C_Reg3_0.5-_2.5-.out");
    while(in_reg3_05_25 >> var1 >> var2 >> var3 >> var4 >> var5 >> var6 >> var7 >> var8 >> var9) {
        cm_reg3_05_25.push_back(var1);
        exe_reg3_05_25.push_back(var2);
        fit_reg3_05_25.push_back(var4);
    }
    in_reg3_05_25.close();

    int size_reg3 = cm_reg3_05_25_15_35.size();

    TGraphErrors *region3 = new TGraphErrors();

    TGraph *gr_reg3_05_25_05_35 = new TGraph();
    TGraph *gr_reg3_05_25_15 = new TGraph();
    TGraph *gr_reg3_05_25_15_35 = new TGraph();
    TGraph *gr_reg3_05_25_25_35 = new TGraph();
    TGraph *gr_reg3_05_25_25 = new TGraph();
    TGraph *gr_reg3_05_25_35 = new TGraph();
    TGraph *gr_reg3_05_25 = new TGraph();

    for(int i = 0; i < size_reg3; i++) {
        region3->SetPoint(i, cm_reg3_05_25_15_35[i] + 1.3, cs_reg3_05_25_15_35[i]*1000.);
        region3->SetPointError(i, 0, cserr_reg3_05_25_15_35[i]*1000.);

        gr_reg3_05_25_15->SetPoint(i, cm_reg3_05_25_15[i] + 1.3, fit_reg3_05_25_15[i]*1000.);

        gr_reg3_05_25_15_35->SetPoint(i, cm_reg3_05_25_15_35[i] + 1.3, fit_reg3_05_25_15_35[i]*1000.);

        gr_reg3_05_25_25_35->SetPoint(i, cm_reg3_05_25_25_35[i] + 1.3, fit_reg3_05_25_25_35[i]*1000.);

        gr_reg3_05_25_25->SetPoint(i, cm_reg3_05_25_25[i] + 1.3, fit_reg3_05_25_25[i]*1000.);

        gr_reg3_05_25_35->SetPoint(i, cm_reg3_05_25_35[i] + 1.3, fit_reg3_05_25_35[i]*1000.);

        gr_reg3_05_25->SetPoint(i, cm_reg3_05_25[i] + 1.3, fit_reg3_05_25[i]*1000.);
    }

    region3->SetMarkerStyle(20);
    region3->SetTitle("; Excitation Energy [MeV]; Cross Section [mb/sr]");
    region3->GetXaxis()->CenterTitle();
    region3->GetYaxis()->CenterTitle();
    region3->SetMinimum(10);
    region3->SetMaximum(150);
    region3->GetXaxis()->SetLimits(1.5, 6.5);
    // region3->Draw("APE");

    ///////////////////
    // Set up graphs //
    ///////////////////

    TCanvas* c = new TCanvas("c");

    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.5, 1, 1.0);
    pad1->SetBottomMargin(0);
    pad1->SetTopMargin(0.2);
    std::cout << pad1->GetTopMargin() << std::endl;
    pad1->Draw();
    pad1->cd();

    region1->Draw("APE");

    // 05_25_15
    gr_reg1_05_25_15->SetLineStyle(1);  gr_reg1_05_25_15->SetLineColor(1); gr_reg1_05_25_15->SetLineWidth(4);
    gr_reg1_05_25_15->Draw("PCsame");

    // 05_25_25
    gr_reg1_05_25_25->SetLineStyle(9); gr_reg1_05_25_25->SetLineColor(3); gr_reg1_05_25_25->SetLineWidth(4);
    gr_reg1_05_25_25->Draw("PCsame");

    // 05_25
    gr_reg1_05_25->SetLineStyle(10); gr_reg1_05_25->SetLineColor(4); gr_reg1_05_25->SetLineWidth(4);
    gr_reg1_05_25->Draw("PCsame");

    // 05_25_25_35
    gr_reg1_05_25_25_35->SetLineStyle(2); gr_reg1_05_25_25_35->SetLineColor(2); gr_reg1_05_25_25_35->SetLineWidth(5);
    gr_reg1_05_25_25_35->Draw("PCsame");

    TLatex latex_reg1;
    latex_reg1.DrawLatex(2, 120., "155^{o}-170^{o}");

    c->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.5);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->Draw();
    pad1->Draw();
    pad2->cd();

    region3->Draw("APE");

    // 05_25_15
    gr_reg3_05_25_15->SetLineStyle(1);  gr_reg3_05_25_15->SetLineColor(1); gr_reg3_05_25_15->SetLineWidth(4);
    gr_reg3_05_25_15->Draw("PCsame");

    // 05_25_25
    gr_reg3_05_25_25->SetLineStyle(9); gr_reg3_05_25_25->SetLineColor(3); gr_reg3_05_25_25->SetLineWidth(4);
    gr_reg3_05_25_25->Draw("PCsame");

    // 05_25
    gr_reg3_05_25->SetLineStyle(10); gr_reg3_05_25->SetLineColor(4); gr_reg3_05_25->SetLineWidth(4);
    gr_reg3_05_25->Draw("PCsame");

    // 05_25_25_35
    gr_reg3_05_25_25_35->SetLineStyle(2); gr_reg3_05_25_25_35->SetLineColor(2); gr_reg3_05_25_25_35->SetLineWidth(5);
    gr_reg3_05_25_25_35->Draw("PCsame");

    TLatex latex_reg3;
    latex_reg3.DrawLatex(2, 120., "105^{o}-145^{o}");

}