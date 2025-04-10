#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMVA/Reader.h>

void applyAndAnalyzeModel_charm(const char* inputFileName, float threshold = 0.3) {

    //---------------------------------------------------------------------------------------------------------
    // Criação do objeto Reader para leitura de resultados 
    //---------------------------------------------------------------------------------------------------------
    
    TMVA::Reader* reader = new TMVA::Reader("!Color:!Silent");

    Float_t pT_c, nConst_c, eta_c, phi_c, mass_c, label_c, eventID_c, score = 0;

    reader->AddVariable("pT_c", &pT_c);
    reader->AddVariable("nConst_c", &nConst_c);
    reader->AddSpectator("eta_c", &eta_c);
    reader->AddSpectator("phi_c", &phi_c);
    reader->AddSpectator("mass_c", &mass_c);
    reader->AddSpectator("label_c", &label_c);
    reader->AddSpectator("eventID_c", &eventID_c);
    reader->BookMVA("LogisticRegression", "dataset_c/weights/TMVARegression_LogisticRegression.weights.xml");


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
    signalTree->SetBranchAddress("label_c", &label_c);
    signalTree->SetBranchAddress("eventID_c", &eventID_c);

    backgroundTree->SetBranchAddress("pT_c", &pT_c);
    backgroundTree->SetBranchAddress("eta_c", &eta_c);
    backgroundTree->SetBranchAddress("phi_c", &phi_c);
    backgroundTree->SetBranchAddress("mass_c", &mass_c);
    backgroundTree->SetBranchAddress("nConst_c", &nConst_c);
    backgroundTree->SetBranchAddress("label_c", &label_c);
    backgroundTree->SetBranchAddress("eventID_c", &eventID_c);

    //---------------------------------------------------------------------------------------------------------
    // Criação de TTrees de saída
    //---------------------------------------------------------------------------------------------------------

    TFile* outputFile = TFile::Open("scoredOutput_2var_charm.root", "RECREATE");

    TTree* outputSignal = new TTree("ScoredSignalTree_c", "Signal tree with ML score");
    outputSignal->Branch("pT_c", &pT_c);
    outputSignal->Branch("eta_c", &eta_c);
    outputSignal->Branch("phi_c", &phi_c);
    outputSignal->Branch("mass_c", &mass_c);
    outputSignal->Branch("nConst_c", &nConst_c);
    outputSignal->Branch("label_c", &label_c);
    outputSignal->Branch("score", &score);
    outputSignal->Branch("eventID_c", &eventID_c);


    TTree* outputBackground = new TTree("ScoredBackgroundTree_c", "Background tree with ML score");
    outputBackground->Branch("pT_c", &pT_c);
    outputBackground->Branch("eta_c", &eta_c);
    outputBackground->Branch("phi_c", &phi_c);
    outputBackground->Branch("mass_c", &mass_c);
    outputBackground->Branch("nConst_c", &nConst_c);
    outputBackground->Branch("label_c", &label_c);
    outputBackground->Branch("score", &score);
    outputBackground->Branch("eventID_c", &eventID_c);
    //---------------------------------------------------------------------------------------------------------
    // Aplicação do modelo aos dados externos
    //---------------------------------------------------------------------------------------------------------

    for (Long64_t i = 0; i < signalTree->GetEntries(); ++i) {
        signalTree->GetEntry(i);
        score = reader->EvaluateMVA("LogisticRegression");
        outputSignal->Fill();
    }

    for (Long64_t i = 0; i < backgroundTree->GetEntries(); ++i) {
        backgroundTree->GetEntry(i);
        score = reader->EvaluateMVA("LogisticRegression");
        outputBackground->Fill();
    }

    outputFile->cd();
    outputSignal->Write();
    outputBackground->Write();

    //---------------------------------------------------------------------------------------------------------
    // Analize de desempenho
    //---------------------------------------------------------------------------------------------------------

    Int_t VP = 0, FN = 0, FP = 0, VN = 0;

    TH1F* h_signal = new TH1F("h_signal", "Scores para eventos de sinal; Score ;Eventos", 100, 0, 1);
    TH1F* h_background = new TH1F("h_background", "Scores para eventos de fundo; Score; Eventos", 100, 0, 1);

    for (Long64_t i = 0; i < outputSignal->GetEntries(); ++i) {
        outputSignal->GetEntry(i);
        h_signal->Fill(score);

        if (label_c == 1) {
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

        if (label_c == 1) {
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

    TLegend* leg = new TLegend(0.65, 0.75, 0.85, 0.90);
    leg->AddEntry(h_signal, "Sinal (label = 1)", "l");
    leg->AddEntry(h_background, "Fundo (label = 0)", "l");
    leg->Draw();

    // Finalizar
    outputFile->Close();
    inputFile->Close();
    delete reader;
}
