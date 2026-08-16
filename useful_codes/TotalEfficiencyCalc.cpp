#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <iostream>
#include <string>
#include <cmath>

void TotalEfficiencyCalc()
{
    int LS = 0;
    bool kornilov = true;
    bool cecil = false;
    const char* treeName_det = "LS_det"; 
    const char* treeName_tar = "Target"; 

    const char* branchName_tar = "Ekin";
    const char* particle_tarZ = "Z";
    const char* particle_tarA= "A";

    const char* branchName_det = "";  
    const char* particleZ = "";
    const char* particleA= "";

    double threshold = 0.;
    double threshold_max = 0.;
    double threshold_min = 0.;

    // Thresholds - Kornilov
    if(kornilov){
        if(LS==0){
            branchName_det = "totalEdep_0";  
            particleZ = "Z_0";
            particleA= "A_0";

            threshold = 0.240;
            threshold_max = 0.301;
            threshold_min = 0.162;}

        if(LS==1){
            branchName_det = "totalEdep_1";  
            particleZ = "Z_1";
            particleA= "A_1";

            threshold = 0.300;
            threshold_max = 0.346;
            threshold_min = 0.240;}

        if(LS==2){
            branchName_det = "totalEdep_2";  
            particleZ = "Z_2";
            particleA= "A_2";

            threshold = 0.240;
            threshold_max = 0.308;
            threshold_min = 0.136;}

        if(LS==3){
            branchName_det = "totalEdep_3";  
            particleZ = "Z_3";
            particleA= "A_3";

            threshold = 0.314;
            threshold_max = 0.346;
            threshold_min = 0.279;}
    }
    // Thresholds - Cecil
    if(cecil){
        if(LS==0){
            branchName_det = "totalEdep_0";  
            particleZ = "Z_0";
            particleA= "A_0";

            threshold = 0.213;
            threshold_max = 0.263;
            threshold_min = 0.153;}

        if(LS==1){
            branchName_det = "totalEdep_1";  
            particleZ = "Z_1";
            particleA= "A_1";

            threshold = 0.263;
            threshold_max = 0.303;
            threshold_min = 0.213;}

        if(LS==2){
            branchName_det = "totalEdep_2";  
            particleZ = "Z_2";
            particleA= "A_2";

            threshold = 0.213;
            threshold_max = 0.272;
            threshold_min = 0.136;}

        if(LS==3){
            branchName_det = "totalEdep_3";  
            particleZ = "Z_3";
            particleA= "A_3";

            threshold = 0.275;
            threshold_max = 0.303;
            threshold_min = 0.245;}
        }

    std::vector<double> energyVec;
    std::vector<double> efficiencyVec;
    std::vector<double> efficiencymaxVec;
    std::vector<double> efficiencyminVec;   

    // Energy ranges
    Int_t energy = 0;
    Int_t energyStep = 100; 
    Int_t maxEnergy = 10000; 
    Int_t minEnergy = 100;  

    char file_name[200], save_name[200];
    for(energy = minEnergy; energy <= maxEnergy; energy+=energyStep)
    {

        sprintf(file_name,"../bld/neutrons/4pi/n_%dkeV.root",energy);
        TFile* file = new TFile(file_name, "READ");

        TTree* tree_det = dynamic_cast<TTree*>(file->Get(treeName_det));
        TTree* tree_tar = dynamic_cast<TTree*>(file->Get(treeName_tar));

        // Variables to store the values from the branches
        Double_t totalEdep_det;
        Double_t Ekin_tar;
        Double_t Z_det;
        Double_t A_det;
        Double_t Z_tar;
        Double_t A_tar;

        // Set the branch addresses
        tree_det->SetBranchAddress(branchName_det, &totalEdep_det);
        tree_det->SetBranchAddress(particleZ, &Z_det);
        tree_det->SetBranchAddress(particleA, &A_det);

        tree_tar->SetBranchAddress(branchName_tar, &Ekin_tar);
        tree_tar->SetBranchAddress(particle_tarZ, &Z_tar);
        tree_tar->SetBranchAddress(particle_tarA, &A_tar);
        
        // Get the number of entries in the trees
        Long64_t numEntries_det = tree_det->GetEntries();
        Long64_t numEntries_hit = tree_tar->GetEntries();

        // Initialize a counter for entries above threshold
        Int_t numEntries_above_th = 0; 
        Int_t numEntries_above_th_max = 0;
        Int_t numEntries_above_th_min = 0; 
        Int_t numEntries_above0 = 0; 

        // Loop over the entries to count the number of entries above threshold
        for (Long64_t i = 0; i < numEntries_det; i++)
        {
            tree_det->GetEntry(i);

            // Check if totalEdep_det is above threshold
            if ((totalEdep_det > threshold && Z_det==1 && A_det==1)||(totalEdep_det > 3.0 && Z_det>=6 && A_det>=12)) numEntries_above_th++;
            if ((totalEdep_det > threshold_max && Z_det==1 && A_det==1)||(totalEdep_det > 3.0 && Z_det>=6 && A_det>=12)) numEntries_above_th_max++;
            if ((totalEdep_det > threshold_min && Z_det==1 && A_det==1)||(totalEdep_det > 3.0 && Z_det>=6 && A_det>=12)) numEntries_above_th_min++;
        }

        // Loop over the entries to count the number of entries above 0
        for (Long64_t i = 0; i < numEntries_hit; i++)
        {
            tree_tar->GetEntry(i);

            // Check if totalEdep_det is above 0
            if (Ekin_tar > 0)
            {
                numEntries_above0++;
            }
        }

        energyVec.push_back(energy);

        // Calculate the percentage value
        double efficiency = (numEntries_above_th / static_cast<double>(numEntries_above0)) * 100.0;
        double efficiency_max = (numEntries_above_th_max / static_cast<double>(numEntries_above0)) * 100.0;
        double efficiency_min = (numEntries_above_th_min / static_cast<double>(numEntries_above0)) * 100.0;

        efficiencyVec.push_back(efficiency);
        efficiencymaxVec.push_back(efficiency_max);
        efficiencyminVec.push_back(efficiency_min);

        // cout<< energy <<"\t"<<numEntries_above_th<<"\t"<<numEntries_above0<<endl;
        cout<<numEntries_above_th<<endl;

        file->Close();

        if(energy>=1500) energyStep=500;
    }

    // Create TGraph objects from the data arrays
    TGraph* graph = new TGraph(energyVec.size(), energyVec.data(), efficiencyVec.data());
    TGraph* graph_max = new TGraph(energyVec.size(), energyVec.data(), efficiencymaxVec.data());
    TGraph* graph_min = new TGraph(energyVec.size(), energyVec.data(), efficiencyminVec.data());

    // Create a TGraphErrors for the band
    TGraphErrors* graph_band = new TGraphErrors(energyVec.size(), energyVec.data(), efficiencymaxVec.data());
    for (int i = 0; i < energyVec.size(); i++) 
    {
        double x, y0, y1, y2;
        double y, y_err;

        graph->GetPoint(i, x, y0);
        graph_max->GetPoint(i, x, y1);
        graph_min->GetPoint(i, x, y2);
        y = (y1 + y2) / 2.0;
        y_err = std::abs(y1 - y2) / 2.0;

        graph_band->SetPoint(i, x, y);
        graph_band->SetPointError(i, 0.0, y_err);
    }

    // Set line and fill properties for the graphs
    graph->SetLineColor(kRed);
    graph_band->SetFillColorAlpha(kBlue, 0.3);
    graph_band->SetFillStyle(1001);

    // Create a TMultiGraph to hold the graphs
    TMultiGraph* mg = new TMultiGraph();
    mg->Add(graph_band);
    mg->Add(graph);

    // Set the font size for the axis titles and numbers
    // gStyle->SetLabelSize(0.05, "X");
    // gStyle->SetLabelSize(0.05, "Y");
    // gStyle->SetTitleSize(0.05, "X");
    // gStyle->SetTitleSize(0.05, "Y");

    // Create a canvas to draw the TMultiGraph
    TCanvas* canvas = new TCanvas("canvas", "Kornilov et al", 1400, 1400); 

    canvas->cd();
    mg->Draw("A3"); // Use "A3" option to include the band
    graph->Draw("L"); // Use "L" option to draw the graph as a line

    if(LS==3)
    {
        mg->SetTitle("Liquid Scintillator 3 (BR60)"); 
        mg->SetName("LS3");
    }
    if(LS==2)
    {
        mg->SetTitle("Liquid Scintillator 2 (BR30)"); 
        mg->SetName("LS2");
    }
    if(LS==1)
    {
        mg->SetTitle("Liquid Scintillator 1 (BL30)"); 
        mg->SetName("LS1");
    }
    if(LS==0)
    {
        mg->SetTitle("Liquid Scintillator 0 (BL60)"); 
        mg->SetName("LS0");
    }
    mg->GetXaxis()->SetTitle("Emitted Neutron Energy (keV)");
    mg->GetYaxis()->SetTitle("Efficiency (%)");

    canvas->Update();
    sprintf(save_name,"bld/neutrons/4pi/totalEff_LS%d_Kornilov.root",LS);
    TFile *fout = new TFile(save_name, "RECREATE");
    
    mg->Write();
    graph->Write();

    sprintf(save_name,"bld/neutrons/4pi/totalEff_LS%d_Kornilov.png",LS);
    canvas->SaveAs(save_name);
    
    return canvas;
    // Clean up
    delete graph;
    delete graph_band;
    delete mg;

    fout->Close();
}
