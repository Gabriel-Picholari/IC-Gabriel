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

    Float_t pT_s, nConst_s, eta_s, phi_s, mass_s, label_s, eventID_s, score = 0;

    reader->AddVariable("pT_s", &pT_s);
    reader->AddVariable("nConst_s", &nConst_s);
    reader->AddSpectator("eta_s", &eta_s);
    reader->AddSpectator("phi_s", &phi_s);
    reader->AddSpectator("mass_s", &mass_s);
    reader->AddSpectator("label_s", &label_s);
    reader->AddSpectator("eventID_s", &eventID_s);
    reader->BookMVA("Likelihood", "dataset_s/weights/TMVAClassification_Likelihood.weights.xml");

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
    signalTree->SetBranchAddress("label_s", &label_s);
    signalTree->SetBranchAddress("eventID_s", &eventID_s);

    backgroundTree->SetBranchAddress("pT_s", &pT_s);
    backgroundTree->SetBranchAddress("eta_s", &eta_s);
    backgroundTree->SetBranchAddress("phi_s", &phi_s);
    backgroundTree->SetBranchAddress("mass_s", &mass_s);
    backgroundTree->SetBranchAddress("nConst_s", &nConst_s);
    backgroundTree->SetBranchAddress("label_s", &label_s);
    backgroundTree->SetBranchAddress("eventID_s", &eventID_s);

    //---------------------------------------------------------------------------------------------------------
    // Criação de TTrees de saída
    //---------------------------------------------------------------------------------------------------------

    TFile* outputFile = TFile::Open("scoredOutput_2var_strange.root", "RECREATE");

    TTree* outputSignal = new TTree("ScoredSignalTree_s", "Signal tree with ML score");
    outputSignal->Branch("pT_s", &pT_s);
    outputSignal->Branch("eta_s", &eta_s);
    outputSignal->Branch("phi_s", &phi_s);
    outputSignal->Branch("mass_s", &mass_s);
    outputSignal->Branch("nConst_s", &nConst_s);
    outputSignal->Branch("label_s", &label_s);
    outputSignal->Branch("score", &score);
    outputSignal->Branch("eventID_s", &eventID_s);


    TTree* outputBackground = new TTree("ScoredBackgroundTree_s", "Background tree with ML score");
    outputBackground->Branch("pT_s", &pT_s);
    outputBackground->Branch("eta_s", &eta_s);
    outputBackground->Branch("phi_s", &phi_s);
    outputBackground->Branch("mass_s", &mass_s);
    outputBackground->Branch("nConst_s", &nConst_s);
    outputBackground->Branch("label_s", &label_s);
    outputBackground->Branch("score", &score);
    outputBackground->Branch("eventID_s", &eventID_s);


    //---------------------------------------------------------------------------------------------------------
    // Aplicação do modelo aos dados externos
    //---------------------------------------------------------------------------------------------------------

    for (Long64_t i = 0; i < signalTree->GetEntries(); ++i) {
        signalTree->GetEntry(i);
        score = reader->EvaluateMVA("Likelihood");
        outputSignal->Fill();
    }

    for (Long64_t i = 0; i < backgroundTree->GetEntries(); ++i) {
        backgroundTree->GetEntry(i);
        score = reader->EvaluateMVA("Likelihood");
        outputBackground->Fill();
    }

    outputFile->cd();
    outputSignal->Write();
    outputBackground->Write();

    //---------------------------------------------------------------------------------------------------------
    // Analize de desempenho
    //---------------------------------------------------------------------------------------------------------

    Int_t VP = 0, FN = 0, FP = 0, VN = 0;

    TH1F* h_signal = new TH1F("h_signal", "Scores para eventos de sinal;Score;Eventos", 100, 0, 1);
    TH1F* h_background = new TH1F("h_background", "Scores para eventos de fundo;Score;Eventos", 100, 0, 1);

    for (Long64_t i = 0; i < outputSignal->GetEntries(); ++i) {
        outputSignal->GetEntry(i);
        h_signal->Fill(score);

        if (label_s == 1) {
            if (score >= threshold) {
                VP++;
            } else {
                FN++;
            }
        } else {
            if (score >= threshold) {
                FP++;
            } else {
                VN++;
            }
        }
    }

    for (Long64_t i = 0; i < outputBackground->GetEntries(); ++i) {
        outputBackground->GetEntry(i);
        h_background->Fill(score);

        if (label_s == 1) {
            if (score >= threshold) {
                VP++;
            } else {
                FN++;
            }
        } else {
            if (score >= threshold) {
                FP++;
            } else {
                VN++;
            }
        }
    }

    Float_t eficiencia;
    if (VP + FN > 0) {
        eficiencia = (float)VP / (VP + FN);
    } else {
        eficiencia = 0;
    }

    Float_t pureza;
    if (VP + FP > 0) {
        pureza = (float)VP / (VP + FP);
    } else {
        pureza = 0;
    }

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

    TCanvas* c1 = new TCanvas("c1", "Distribuição dos scores", 800, 600);
    h_signal->SetLineColor(kRed);
    h_background->SetLineColor(kBlue);
    h_background->DrawCopy();
    h_signal->DrawCopy("same");

    TLegend* leg = new TLegend(0.65, 0.75, 0.88, 0.88);
    leg->AddEntry(h_signal, "Sinal (label = 1)", "l");
    leg->AddEntry(h_background, "Fundo (label = 0)", "l");
    leg->Draw();

    // Finalizar
    outputFile->Close();
    inputFile->Close();
    delete reader;
}
