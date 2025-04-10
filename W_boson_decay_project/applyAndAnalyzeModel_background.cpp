#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMVA/Reader.h>

void applyAndAnalyzeModel_background(const char* inputFileName, float threshold = 0.5) {

    //---------------------------------------------------------------------------------------------------------
    // Criação do objeto Reader para leitura de resultados 
    //---------------------------------------------------------------------------------------------------------
    
    TMVA::Reader* reader = new TMVA::Reader("!Color:!Silent");

    Float_t pT_bkg, nConst_bkg, eta_bkg, phi_bkg, mass_bkg, label_bkg, eventID_bkg, score = 0;

    reader->AddVariable("pT_bkg", &pT_bkg);
    reader->AddVariable("nConst_bkg", &nConst_bkg);
    reader->AddSpectator("eta_bkg", &eta_bkg);
    reader->AddSpectator("phi_bkg", &phi_bkg);
    reader->AddSpectator("mass_bkg", &mass_bkg);
    reader->AddSpectator("label_bkg", &label_bkg);
    reader->AddSpectator("eventID_bkg", &eventID_bkg);
    reader->BookMVA("LogisticRegression", "dataset_bkg/weights/TMVARegression_LogisticRegression.weights.xml");

    //---------------------------------------------------------------------------------------------------------
    // Recuperação de TTrees de entrada 
    //---------------------------------------------------------------------------------------------------------

    TFile* inputFile = TFile::Open(inputFileName, "READ");
    TTree* signalTree = (TTree*)inputFile->Get("SignalTree_bkg");
    TTree* backgroundTree = (TTree*)inputFile->Get("BackgroundTree_bkg");

    signalTree->SetBranchAddress("pT_bkg", &pT_bkg);
    signalTree->SetBranchAddress("eta_bkg", &eta_bkg);
    signalTree->SetBranchAddress("phi_bkg", &phi_bkg);
    signalTree->SetBranchAddress("mass_bkg", &mass_bkg);
    signalTree->SetBranchAddress("nConst_bkg", &nConst_bkg);
    signalTree->SetBranchAddress("label_bkg", &label_bkg);
    signalTree->SetBranchAddress("eventID_bkg", &eventID_bkg);

    backgroundTree->SetBranchAddress("pT_bkg", &pT_bkg);
    backgroundTree->SetBranchAddress("eta_bkg", &eta_bkg);
    backgroundTree->SetBranchAddress("phi_bkg", &phi_bkg);
    backgroundTree->SetBranchAddress("mass_bkg", &mass_bkg);
    backgroundTree->SetBranchAddress("nConst_bkg", &nConst_bkg);
    backgroundTree->SetBranchAddress("label_bkg", &label_bkg);
    backgroundTree->SetBranchAddress("eventID_bkg", &eventID_bkg);

    //---------------------------------------------------------------------------------------------------------
    // Criação de TTrees de saída
    //---------------------------------------------------------------------------------------------------------

    TFile* outputFile = TFile::Open("scoredOutput_2var_background.root", "RECREATE");

    TTree* outputSignal = new TTree("ScoredSignalTree_bkg", "Signal tree with ML score");
    outputSignal->Branch("pT_bkg", &pT_bkg);
    outputSignal->Branch("eta_bkg", &eta_bkg);
    outputSignal->Branch("phi_bkg", &phi_bkg);
    outputSignal->Branch("mass_bkg", &mass_bkg);
    outputSignal->Branch("nConst_bkg", &nConst_bkg);
    outputSignal->Branch("label_bkg", &label_bkg);
    outputSignal->Branch("score", &score);
    outputSignal->Branch("eventID_bkg", &eventID_bkg);

    TTree* outputBackground = new TTree("ScoredBackgroundTree_bkg", "Background tree with ML score");
    outputBackground->Branch("pT_bkg", &pT_bkg);
    outputBackground->Branch("eta_bkg", &eta_bkg);
    outputBackground->Branch("phi_bkg", &phi_bkg);
    outputBackground->Branch("mass_bkg", &mass_bkg);
    outputBackground->Branch("nConst_bkg", &nConst_bkg);
    outputBackground->Branch("label_bkg", &label_bkg);
    outputBackground->Branch("score", &score);
    outputBackground->Branch("eventID_bkg", &eventID_bkg);
    // Até aqui foi debugado e está tudo certo!
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

        if (label_bkg == 1) {
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

        if (label_bkg == 1) {
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
    h_signal->DrawCopy();
    h_background->DrawCopy("same");

    TLegend* leg = new TLegend(0.65, 0.75, 0.85, 0.90);
    leg->AddEntry(h_signal, "Sinal (label = 1)", "l");
    leg->AddEntry(h_background, "Fundo (label = 0)", "l");
    leg->Draw();

    // Finalizar
    outputFile->Close();
    inputFile->Close();
    delete reader;
}
