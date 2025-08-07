#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMVA/Reader.h>

void applyAndAnalyzeModel_strange(const char* inputFileName, float threshold = 0.5) {

    //---------------------------------------------------------------------------------------------------------
    // Criação do objeto Reader para leitura de resultados 
    //---------------------------------------------------------------------------------------------------------
    
    TMVA::Reader* reader = new TMVA::Reader("!Color:!Silent");

    Float_t pT_s, nConst_s, eta_s, phi_s, mass_s, label_s, eventID_s, score,maxRho_s = 0;

    reader->AddVariable("pT_s", &pT_s);
    reader->AddVariable("nConst_s", &nConst_s);
    //reader->AddVariable("maxRho_s", &maxRho_s);
    reader->AddSpectator("eta_s", &eta_s);
    reader->AddSpectator("phi_s", &phi_s);
    reader->AddSpectator("mass_s", &mass_s);
    reader->AddSpectator("label_s", &label_s);
    reader->AddSpectator("eventID_s", &eventID_s);
    reader->BookMVA("GradBoost", "dataset_s/weights/TMVAClassification_GradBoost.weights.xml");

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
    //signalTree->SetBranchAddress("maxRho_s", &maxRho_s);
    signalTree->SetBranchAddress("label_s", &label_s);
    signalTree->SetBranchAddress("eventID_s", &eventID_s);

    backgroundTree->SetBranchAddress("pT_s", &pT_s);
    backgroundTree->SetBranchAddress("eta_s", &eta_s);
    backgroundTree->SetBranchAddress("phi_s", &phi_s);
    backgroundTree->SetBranchAddress("mass_s", &mass_s);
    backgroundTree->SetBranchAddress("nConst_s", &nConst_s);
    //backgroundTree->SetBranchAddress("maxRho_s", &maxRho_s);
    backgroundTree->SetBranchAddress("label_s", &label_s);
    backgroundTree->SetBranchAddress("eventID_s", &eventID_s);

    //---------------------------------------------------------------------------------------------------------
    // Análise de desempenho
    //---------------------------------------------------------------------------------------------------------

    Int_t VP = 0, FN = 0, FP = 0, VN = 0;

    TH1F* h_signal = new TH1F("h_signal", "TMVA response for classifier: GradBoost;GradBoost response;(1/N) dN / dx", 100, -1, 1);
    TH1F* h_background = new TH1F("h_background", "TMVA response for classifier: GradBoost;GradBoost response;(1/N) dN / dx", 100, -1, 1);

    for (Long64_t i = 0; i < signalTree->GetEntries(); ++i) {
        signalTree->GetEntry(i);
        score = reader->EvaluateMVA("GradBoost");
        h_signal->Fill(score);

        if (label_s == 1) {
            if (score >= threshold) VP++;
            else FN++;
        } else {
            if (score >= threshold) FP++;
            else VN++;
        }
    }

    for (Long64_t i = 0; i < backgroundTree->GetEntries(); ++i) {
        backgroundTree->GetEntry(i);
        score = reader->EvaluateMVA("GradBoost");
        h_background->Fill(score);

        if (label_s == 1) {
            if (score >= threshold) VP++;
            else FN++;
        } else {
            if (score >= threshold) FP++;
            else VN++;
        }
    }

    //---------------------------------------------------------------------------------------------------------
    // Normalização e métricas
    //---------------------------------------------------------------------------------------------------------

    Float_t eficiencia = (VP + FN > 0) ? (float)VP / (VP + FN) : 0;
    Float_t pureza = (VP + FP > 0) ? (float)VP / (VP + FP) : 0;

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

    TCanvas* c1 = new TCanvas("c1", "Distribuição dos scores", 900, 700);
    c1->SetGrid();

    h_signal->SetLineColor(kBlue + 1);
    h_signal->SetFillColorAlpha(kBlue, 0.35);
    h_signal->SetLineWidth(2);

    h_background->SetLineColor(kRed + 1);
    h_background->SetFillColorAlpha(kRed, 0.35);
    h_background->SetLineWidth(2);

    h_signal->SetStats(0);
    h_background->SetStats(0);

    //h_signal->Scale(1.0 / h_signal->Integral());
    //h_background->Scale(1.0 / h_background->Integral());

    h_background->DrawCopy("HIST");
    h_signal->DrawCopy("HIST SAME");

    TLegend* leg = new TLegend(0.65, 0.75, 0.88, 0.90);
    leg->SetBorderSize(1);
    leg->SetFillColorAlpha(0, 0.3);
    leg->SetTextSize(0.03);
    leg->AddEntry(h_signal, "Sinal (label = 1)", "f");
    leg->AddEntry(h_background, "Fundo (label = 0)", "f");
    leg->Draw();

    // Finalizar
    inputFile->Close();
    delete reader;
}
