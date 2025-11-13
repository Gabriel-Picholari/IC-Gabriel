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

void applyAndAnalyzeModel_3var_charm(const char* inputFileName, float threshold = 0.7) {

    //---------------------------------------------------------------------------------------------------------
    // Criação do objeto Reader para leitura de resultados 
    //---------------------------------------------------------------------------------------------------------

    TMVA::Reader* reader = new TMVA::Reader("!Color:!Silent");

    Float_t pT_c, nConst_c, eta_c, phi_c, mass_c, label_c, eventID_c, score, nRho_c = 0;

    reader->AddVariable("pT_c", &pT_c);
    reader->AddVariable("nRho_c", &nRho_c);
    reader->AddVariable("nConst_c", &nConst_c);

    reader->AddSpectator("eta_c", &eta_c);
    reader->AddSpectator("phi_c", &phi_c);
    reader->AddSpectator("mass_c", &mass_c);
    reader->AddSpectator("label_c", &label_c);
    reader->AddSpectator("eventID_c", &eventID_c);
    reader->BookMVA("GradBoost", "dataset_c_3var/weights/TMVAClassification_GradBoost.weights.xml");

    //---------------------------------------------------------------------------------------------------------
    // Recuperação de TTrees de entrada 
    //---------------------------------------------------------------------------------------------------------

    TFile* inputFile = TFile::Open(inputFileName, "READ");
    TTree* signalTree = (TTree*)inputFile->Get("SignalTree_c");
    TTree* backgroundTree = (TTree*)inputFile->Get("BackgroundTree_c");

    signalTree->SetBranchAddress("pT_c", &pT_c);
    signalTree->SetBranchAddress("eta_c", &eta_c);
    signalTree->SetBranchAddress("phi_c", &phi_c);
    signalTree->SetBranchAddress("mass_c", &mass_c);
    signalTree->SetBranchAddress("nConst_c", &nConst_c);
    signalTree->SetBranchAddress("nRho_c", &nRho_c);
    signalTree->SetBranchAddress("label_c", &label_c);
    signalTree->SetBranchAddress("eventID_c", &eventID_c);

    backgroundTree->SetBranchAddress("pT_c", &pT_c);
    backgroundTree->SetBranchAddress("eta_c", &eta_c);
    backgroundTree->SetBranchAddress("phi_c", &phi_c);
    backgroundTree->SetBranchAddress("mass_c", &mass_c);
    backgroundTree->SetBranchAddress("nConst_c", &nConst_c);
    backgroundTree->SetBranchAddress("nRho_c", &nRho_c);
    backgroundTree->SetBranchAddress("label_c", &label_c);
    backgroundTree->SetBranchAddress("eventID_c", &eventID_c);

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
        allLabels.push_back((int)label_c);

        if (label_c == 1) 
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
        allLabels.push_back((int)label_c);

        if (label_c == 1) 
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

    std::cout << "--------------Charm--------------" << std::endl;
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

    TCanvas* c1 = new TCanvas("c1", "GradBoost Score Distribution (Charm)", 900, 700);
    c1->SetGrid();

    h_background->SetLineColor(kRed);
    h_signal->SetLineColor(kGreen);

    h_signal->SetTitle("GradBoost Score Distribution for Charm Jets;Score;Number of Events");
    h_background->SetTitle("GradBoost Score Distribution for Charm Jets;Score;Number of Events");

    h_background->DrawCopy();
    h_signal->DrawCopy("same");

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

    TH1F* hEff = new TH1F("hEff","Efficiency vs Threshold (Charm);Threshold;Efficiency", nSteps, tmin, tmax);
    TH1F* hPur = new TH1F("hPur","Purity vs Threshold (Charm);Threshold;Purity", nSteps, tmin, tmax);

    for (int k=0; k<nSteps; ++k) {
        float thr = tmin + (tmax-tmin)*k/(nSteps-1);
        hEff->SetBinContent(k+1, vEff[k]);
        hPur->SetBinContent(k+1, vPur[k]);
    }

    TCanvas* c2 = new TCanvas("c2", "Efficiency and Purity vs Threshold (Charm)", 900, 700);
    c2->SetGrid();  

    TGraph* gEff = new TGraph(nSteps, vx.data(), vEff.data());
    gEff->SetTitle("Efficiency and Purity vs Threshold for Charm Jets;Threshold;Value");
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

    inputFile->Close();
    delete reader;
}
