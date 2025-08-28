#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMVA/Reader.h>

void applyAndAnalyzeModel_3var_charm(const char* inputFileName, float threshold = 0.0) {

    //---------------------------------------------------------------------------------------------------------
    // Criação do objeto Reader para leitura de resultados 
    //---------------------------------------------------------------------------------------------------------
    
    TMVA::Reader* reader = new TMVA::Reader("!Color:!Silent");

    Float_t pT_c, nConst_c, eta_c, phi_c, mass_c, label_c, eventID_c, score, maxRho_c, nRho_c = 0;

    reader->AddVariable("pT_c", &pT_c);
    reader->AddVariable("nRho_c", &nRho_c);
    reader->AddVariable("nConst_c", &nConst_c);
    //reader->AddVariable("maxRho_c", &maxRho_c); Discontinued

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
    //signalTree->SetBranchAddress("maxRho_c", &maxRho_c); Discontinued
    signalTree->SetBranchAddress("label_c", &label_c);
    signalTree->SetBranchAddress("eventID_c", &eventID_c);

    backgroundTree->SetBranchAddress("pT_c", &pT_c);
    backgroundTree->SetBranchAddress("eta_c", &eta_c);
    backgroundTree->SetBranchAddress("phi_c", &phi_c);
    backgroundTree->SetBranchAddress("mass_c", &mass_c);
    backgroundTree->SetBranchAddress("nConst_c", &nConst_c);
    backgroundTree->SetBranchAddress("nRho_c", &nRho_c);
    //backgroundTree->SetBranchAddress("maxRho_c", &maxRho_c); Discontinued
    backgroundTree->SetBranchAddress("label_c", &label_c);
    backgroundTree->SetBranchAddress("eventID_c", &eventID_c);

    //---------------------------------------------------------------------------------------------------------
    // Analize de desempenho
    //---------------------------------------------------------------------------------------------------------

    Int_t VP = 0, FN = 0, FP = 0, VN = 0;

    TH1F* h_signal = new TH1F("h_signal", "TMVA response for classifier: GradBoost;GradBoost response;Events", 100, -1, 1);
    TH1F* h_background = new TH1F("h_background", "TMVA response for classifier: GradBoost;GradBoost response;Events", 100, -1, 1);


    for (Long64_t i = 0; i < signalTree->GetEntries(); ++i) {
        signalTree->GetEntry(i);
        score = reader->EvaluateMVA("GradBoost");
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

    for (Long64_t i = 0; i < backgroundTree->GetEntries(); ++i) {
        backgroundTree->GetEntry(i);
        score = reader->EvaluateMVA("GradBoost");
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
    // Plotar histogramas
    //---------------------------------------------------------------------------------------------------------

    TCanvas* c1 = new TCanvas("c1", "Charmed score distribution", 900, 700);
    c1->SetGrid();
    h_signal->SetLineColor(kGreen );
    h_background->SetLineColor(kRed );

    //h_signal->Scale(1.0 / h_signal->Integral());
    //h_background->Scale(1.0 / h_background->Integral());


    h_signal->SetTitle("Distribuição do Score do GradBoost");
    h_signal->GetXaxis()->SetTitle("Score");
    h_signal->GetYaxis()->SetTitle("Eventos");

    // Primeiro fundo, depois sinal por cima
    h_background->DrawCopy();
    h_signal->DrawCopy("same");

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
