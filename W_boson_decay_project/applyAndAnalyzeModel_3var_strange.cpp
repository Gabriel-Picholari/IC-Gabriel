#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMVA/Reader.h>

// >>> ADIÇÃO
#include <TGraph.h>
#include <TLine.h>
#include <vector>

void applyAndAnalyzeModel_3var_strange(const char* inputFileName, float threshold = 0.3) {

    //---------------------------------------------------------------------------------------------------------
    // Criação do objeto Reader para leitura de resultados 
    //---------------------------------------------------------------------------------------------------------
    
    TMVA::Reader* reader = new TMVA::Reader("!Color:!Silent");

    Float_t pT_s, nConst_s, eta_s, phi_s, mass_s, label_s, eventID_s, score, maxRho_s, nRho_s = 0;

    reader->AddVariable("pT_s", &pT_s);
    reader->AddVariable("nRho_s", &nRho_s);
    reader->AddVariable("nConst_s", &nConst_s);
    //reader->AddVariable("maxRho_s", &maxRho_s); Discontinued
    reader->AddSpectator("eta_s", &eta_s);
    reader->AddSpectator("phi_s", &phi_s);
    reader->AddSpectator("mass_s", &mass_s);
    reader->AddSpectator("label_s", &label_s);
    reader->AddSpectator("eventID_s", &eventID_s);
    reader->BookMVA("GradBoost", "dataset_s_3var/weights/TMVAClassification_GradBoost.weights.xml");

    //---------------------------------------------------------------------------------------------------------
    // Recuperação de TTrees de entrada 
    //---------------------------------------------------------------------------------------------------------

    TFile* inputFile = TFile::Open(inputFileName, "READ");
    TTree* signalTree = (TTree*)inputFile->Get("SignalTree_s");
    TTree* backgroundTree = (TTree*)inputFile->Get("BackgroundTree_s");

    signalTree->SetBranchAddress("pT_s", &pT_s);
    signalTree->SetBranchAddress("eta_s", &eta_s);
    signalTree->SetBranchAddress("phi_s", &phi_s);
    signalTree->SetBranchAddress("mass_s", &mass_s);
    signalTree->SetBranchAddress("nConst_s", &nConst_s);
    signalTree->SetBranchAddress("nRho_s", &nRho_s);
    //signalTree->SetBranchAddress("maxRho_s", &maxRho_s); Discontinued
    signalTree->SetBranchAddress("label_s", &label_s);
    signalTree->SetBranchAddress("eventID_s", &eventID_s);

    backgroundTree->SetBranchAddress("pT_s", &pT_s);
    backgroundTree->SetBranchAddress("eta_s", &eta_s);
    backgroundTree->SetBranchAddress("phi_s", &phi_s);
    backgroundTree->SetBranchAddress("mass_s", &mass_s);
    backgroundTree->SetBranchAddress("nConst_s", &nConst_s);
    backgroundTree->SetBranchAddress("nRho_s", &nRho_s);
    //backgroundTree->SetBranchAddress("maxRho_s", &maxRho_s); Discontinued
    backgroundTree->SetBranchAddress("label_s", &label_s);
    backgroundTree->SetBranchAddress("eventID_s", &eventID_s);

    //---------------------------------------------------------------------------------------------------------
    // Análise de desempenho
    //---------------------------------------------------------------------------------------------------------

    Int_t VP = 0, FN = 0, FP = 0, VN = 0;

    TH1F* h_signal = new TH1F("h_signal", "TMVA response for classifier: GradBoost;GradBoost response;Events", 100, -1, 1);
    TH1F* h_background = new TH1F("h_background", "TMVA response for classifier: GradBoost;GradBoost response;Events", 100, -1, 1);

    // Armazenando todos os scores e labels
    std::vector<float> allScores;
    std::vector<int>   allLabels;

    for (Long64_t i = 0; i < signalTree->GetEntries(); ++i) {
        signalTree->GetEntry(i);
        score = reader->EvaluateMVA("GradBoost");
        h_signal->Fill(score);

        allScores.push_back(score);
        allLabels.push_back((int)label_s);

        if (label_s == 1) 
        {
            if (score >= threshold) VP++;
            else FN++;
        } 
        else 
        {
            if (score >= threshold) FP++;
            else VN++;
        }
    }

    for (Long64_t i = 0; i < backgroundTree->GetEntries(); ++i) 
    {
        backgroundTree->GetEntry(i);
        score = reader->EvaluateMVA("GradBoost");
        h_background->Fill(score);

        allScores.push_back(score);
        allLabels.push_back((int)label_s);

        if (label_s == 1) 
        {
            if (score >= threshold) VP++;
            else FN++;
        } 
        else 
        {
            if (score >= threshold) FP++;
            else VN++;
        }
    }

    //---------------------------------------------------------------------------------------------------------
    // Normalização e métricas
    //---------------------------------------------------------------------------------------------------------

    // (substitui ternários por if/else como você pediu)
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

    std::cout << "-------------Strange-------------" << std::endl;
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
    // Plotar histogramas
    //---------------------------------------------------------------------------------------------------------

    TCanvas* c1 = new TCanvas("c1", "Strange score distribution", 900, 700);
    c1->SetGrid();

    h_background->SetLineColor(kRed);
    h_signal->SetLineColor(kGreen);

    h_background->DrawCopy();
    h_signal->DrawCopy("same");

    //---------------------------------------------------------------------------------------------------------
    // >>> ADIÇÃO: scan de thresholds e gráficos sobrepostos (Eficiência & Pureza)
    //---------------------------------------------------------------------------------------------------------
    int   nSteps = 100;
    float tmin   = -1.0f;
    float tmax   =  1.0f;

    std::vector<double> vx, vEff, vPur;
    vx.reserve(nSteps);
    vEff.reserve(nSteps);
    vPur.reserve(nSteps);

    for (int k = 0; k < nSteps; ++k) 
    {
        float thr = tmin + (tmax - tmin) * k / (float)(nSteps - 1);

        int vp = 0, fn = 0, fp = 0, vn = 0;
        for (size_t j = 0; j < allScores.size(); ++j) 
        {
            bool predS = (allScores[j] >= thr);
            if (allLabels[j] == 1) 
            {
                if (predS) vp++;
                else       fn++;
            } 
            else 
            {
                if (predS) fp++;
                else       vn++;
            }
        }

        float eff;
        if (vp + fn > 0) 
        {
            eff = (float)vp / (vp + fn);
        } 
        else 
        {
            eff = 0.0f;
        }

        float pur;
        if (vp + fp > 0) 
        {
            pur = (float)vp / (vp + fp);
        } 
        else 
        {
            pur = 0.0f;
        }

        vx.push_back(thr);
        vEff.push_back(eff);
        vPur.push_back(pur);
    }

    TCanvas* c2 = new TCanvas("c2","Eficiência e Pureza vs Threshold",900,700);
    c2->SetGrid();

    TGraph* gEff = new TGraph(nSteps, vx.data(), vEff.data());
    gEff->SetTitle("Eficiência e Pureza vs Threshold;Threshold;Valor");
    gEff->SetLineColor(kBlue);
    gEff->SetLineWidth(2);
    gEff->Draw("AL");

    TGraph* gPur = new TGraph(nSteps, vx.data(), vPur.data());
    gPur->SetLineColor(kMagenta);
    gPur->SetLineWidth(2);
    gPur->Draw("L same");

    TLine* lth = new TLine(threshold, 0, threshold, 1);
    lth->SetLineStyle(2);
    lth->Draw("same");

    inputFile->Close();
    delete reader;
}
