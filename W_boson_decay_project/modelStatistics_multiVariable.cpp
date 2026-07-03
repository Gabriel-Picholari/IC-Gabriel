#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMVA/Reader.h>
#include <TGraph.h>
#include <vector>
#include <TLine.h>

void modelStatistics_multiVariable(const char* inputFileName, std::string switch_string, std::string contaminatingGluonMode, const float threshold = 0.7) 
{
    std::string sfx;
    std::string long_sfx;
    std::string short_sfx;
    std::string uppercase_switch_string;
    std::string datasetName;

    if (switch_string == "strange")
    {
        sfx = "_s";
        long_sfx = "_strange";
        short_sfx = "s";
        uppercase_switch_string = "Strange";
    } 
    else if (switch_string == "charm")
    {
        sfx = "_c";
        long_sfx = "_charm";
        short_sfx = "c";
        uppercase_switch_string = "Charm";
    }

    if (switch_string == "charm")  
    {
        if (contaminatingGluonMode == "include") datasetName = "dataset_c_multiVariable_gluonJetsIncluded";
        if (contaminatingGluonMode == "exclude") datasetName = "dataset_c_multiVariable_gluonJetsExcluded";
    }
    if (switch_string == "strange")
    {
        if (contaminatingGluonMode == "include") datasetName = "dataset_s_multiVariable_gluonJetsIncluded";
        if (contaminatingGluonMode == "exclude") datasetName = "dataset_s_multiVariable_gluonJetsExcluded";
    }

    //---------------------------------------------------------------------------------------------------------
    // Criação do objeto Reader para leitura de resultados 
    //---------------------------------------------------------------------------------------------------------

    TMVA::Reader* reader = new TMVA::Reader("!Color:!Silent");

    Float_t pT, nConst, eta, phi, mass, label, eventID, score, nRho = 0;

    reader->AddVariable(("pT" + sfx).c_str(), &pT);
    reader->AddVariable(("nRho" + sfx).c_str(), &nRho);
    reader->AddVariable(("nConst" + sfx).c_str(), &nConst);

    reader->AddSpectator(("eta" + sfx).c_str(), &eta);
    reader->AddSpectator(("phi" + sfx).c_str(), &phi);
    reader->AddSpectator(("label" + sfx).c_str(), &label);
    reader->AddSpectator(("eventID" + sfx).c_str(), &eventID);
    reader->BookMVA("GradBoost", (datasetName + "/weights/TMVAClassification_GradBoost.weights.xml").c_str());

    //---------------------------------------------------------------------------------------------------------
    // Recuperação de TTrees de entrada 
    //---------------------------------------------------------------------------------------------------------

    TFile* inputFile = TFile::Open(inputFileName, "READ");
    TTree* signalTree = (TTree*)inputFile->Get(("SignalTree" + sfx).c_str());
    TTree* backgroundTree = (TTree*)inputFile->Get(("BackgroundTree" + sfx).c_str());

    signalTree->SetBranchAddress(("pT" + sfx).c_str(), &pT);
    signalTree->SetBranchAddress(("eta" + sfx).c_str(), &eta);
    signalTree->SetBranchAddress(("phi" + sfx).c_str(), &phi);
    signalTree->SetBranchAddress(("nConst" + sfx).c_str(), &nConst);
    signalTree->SetBranchAddress(("nRho" + sfx).c_str(), &nRho);
    signalTree->SetBranchAddress(("label" + sfx).c_str(), &label);
    signalTree->SetBranchAddress(("eventID" + sfx).c_str(), &eventID);

    backgroundTree->SetBranchAddress(("pT" + sfx).c_str(), &pT);
    backgroundTree->SetBranchAddress(("eta" + sfx).c_str(), &eta);
    backgroundTree->SetBranchAddress(("phi" + sfx).c_str(), &phi);
    backgroundTree->SetBranchAddress(("nConst" + sfx).c_str(), &nConst);
    backgroundTree->SetBranchAddress(("nRho" + sfx).c_str(), &nRho);
    backgroundTree->SetBranchAddress(("label" + sfx).c_str(), &label);
    backgroundTree->SetBranchAddress(("eventID" + sfx).c_str(), &eventID);

    //---------------------------------------------------------------------------------------------------------
    // Analise de desempenho (threshold único)
    //---------------------------------------------------------------------------------------------------------
    Int_t VP = 0, FN = 0, FP = 0, VN = 0;

    TH1F* h_signal = new TH1F("h_signal", "TMVA response for classifier: GradBoost;GradBoost response;Events", 100, -1, 1);
    TH1F* h_background = new TH1F("h_background", "TMVA response for classifier: GradBoost;GradBoost response;Events", 100, -1, 1);

    // Vamos guardar todos os scores e labels
    std::vector<float> allScores;
    std::vector<int>   allLabels;

    for (Long64_t i = 0; i < signalTree->GetEntries(); ++i) 
    {
        signalTree->GetEntry(i);
        score = reader->EvaluateMVA("GradBoost");
        h_signal->Fill(score);
        allScores.push_back(score);
        allLabels.push_back((int)label);

        if (label == 1) 
        {
            if (score >= threshold) VP++; else FN++;
        } 
        else 
        {
            if (score >= threshold) FP++; else VN++;
        }
    }

    for (Long64_t i = 0; i < backgroundTree->GetEntries(); ++i) 
    {
        backgroundTree->GetEntry(i);
        score = reader->EvaluateMVA("GradBoost");
        h_background->Fill(score);
        allScores.push_back(score);
        allLabels.push_back((int)label);

        if (label == 1) 
        {
            if (score >= threshold) VP++; else FN++;
        } 
        else 
        {
            if (score >= threshold) FP++; else VN++;
        }
    }

    Float_t eficiencia;
    if (VP + FN > 0) 
    {
        eficiencia = (float)VP / (VP + FN);
    } 
    else 
    {
        eficiencia = 0;
    }

    Float_t pureza;
    if (VP + FP > 0) 
    {
        pureza = (float)VP / (VP + FP);
    } 
    else 
    {
        pureza = 0;
    }

    std::cout << "--------------------------------------" << std::endl;
    std::cout << "\nMatriz de Confusão (threshold = " << threshold << "):" << std::endl;
    std::cout << "--------------------------------------" << std::endl;
    std::cout << "             | Pred: S | Pred: B     " << std::endl;
    std::cout << "--------------------------------------" << std::endl;
    std::cout << "Real: S (1)  |   " << VP << "     |   " << FN << std::endl;
    std::cout << "Real: B (0)  |   " << FP << "     |   " << VN << std::endl;
    std::cout << "--------------------------------------" << std::endl;
    std::cout << "Eficiência (Recall) = " << eficiencia << std::endl;
    std::cout << "Pureza (Precision)  = " << pureza << std::endl;

//---------------------------------------------------------------------------------------------------------
// Plot original histograms
//---------------------------------------------------------------------------------------------------------

    TCanvas* c1 = new TCanvas("c1", ("GradBoost Score Distribution (" + uppercase_switch_string + ")").c_str(), 900, 700);
    c1->SetGrid();

    h_background->SetLineColor(kRed);
    h_signal->SetLineColor(kGreen);

    h_signal->SetTitle(("GradBoost Score Distribution for " + uppercase_switch_string + " Jets;Score;Entries").c_str());
    h_background->SetTitle(("GradBoost Score Distribution for " + uppercase_switch_string + " Jets;Score;Entries").c_str());

    h_background->DrawCopy();
    h_signal->DrawCopy("same");

    //c1->SaveAs(("GradBoost_Score_Distribution" + long_sfx + ".png").c_str());

    //---------------------------------------------------------------------------------------------------------
    // Scan em thresholds
    //---------------------------------------------------------------------------------------------------------

    int nSteps = 1000;
    float tmin = -1.0, tmax = 1.0;

    std::vector<double> vx, vEff, vPur;
    vx.reserve(nSteps); vEff.reserve(nSteps); vPur.reserve(nSteps);

    for (int k=0; k<nSteps; ++k) 
    {
        Float_t thr = tmin + (tmax-tmin)*k/(nSteps-1);
        int vp=0, fn=0, fp=0, vn=0;

        for (size_t j=0;j<allScores.size();++j) 
        {
            bool predS = (allScores[j] >= thr);
            if (allLabels[j]==1) 
            { 
                if(predS) vp++; 
                else fn++; 
            }
            else 
            { 
                if(predS) fp++; 
                else vn++; 
            }
        }

        Float_t eff, pur;
        if (vp + fn > 0) 
        {
            eff = (float)vp / (vp + fn);
        } 
        else 
        {
            eff = 0.0;
        }

        if (vp + fp > 0) 
        {
            pur = (float)vp / (vp + fp);
        } 
        else 
        {
            pur = 0.0;
        }

        vx.push_back(thr);
        vEff.push_back(eff);
        vPur.push_back(pur);
    }

    TH1F* hEff = new TH1F("hEff",("Efficiency vs Threshold (" + uppercase_switch_string + ");Threshold;Efficiency").c_str(), nSteps, tmin, tmax);
    TH1F* hPur = new TH1F("hPur",("Purity vs Threshold (" + uppercase_switch_string + ");Threshold;Purity").c_str(), nSteps, tmin, tmax);

    for (int k=0; k<nSteps; ++k) {
        float thr = tmin + (tmax-tmin)*k/(nSteps-1);
        hEff->SetBinContent(k+1, vEff[k]);
        hPur->SetBinContent(k+1, vPur[k]);
    }

    TCanvas* c2 = new TCanvas("c2", ("Efficiency and Purity vs Threshold (" + uppercase_switch_string + ")").c_str(), 2500, 2500);
    c2->SetGrid();

    TGraph* gEff = new TGraph(nSteps, vx.data(), vEff.data());
    gEff->SetTitle(("Efficiency and Purity vs Threshold for " + uppercase_switch_string + " Jets;Threshold;Value").c_str());
    gEff->SetLineColor(kBlue);
    gEff->SetLineWidth(2);
    gEff->Draw("AL");

    TGraph* gPur = new TGraph(nSteps, vx.data(), vPur.data());
    gPur->SetLineColor(kMagenta);
    gPur->SetLineWidth(2);
    gPur->Draw("L same");

    TLine* l1 = new TLine(threshold,0,threshold,1);
    l1->SetLineStyle(2);
    l1->Draw("same");

    //c2->SaveAs(("Efficiency_Purity_vs_Threshold_" + long_sfx + ".png").c_str());

    inputFile->Close();
    delete reader;
}
