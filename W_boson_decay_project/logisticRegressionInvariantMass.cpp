#include <map>
#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMVA/Reader.h>
#include <TLorentzVector.h>

void logisticRegressionInvariantMass(const char* inputFileName_c, const char* inputFileName_s, float threshold = 0.3) {

    //---------------------------------------------------------------------------------------------------------
    // Criação dos objetos Readers para leitura de resultados referentes tanto ao c como ao s
    //---------------------------------------------------------------------------------------------------------

    TMVA::Reader* reader_c = new TMVA::Reader("!Color:!Silent");

    Float_t eventID_c, pT_c, eta_c, phi_c, mass_c, nConst_c = 0;
    Float_t label_c, score_c = 0;

    reader_c->AddVariable("pT_c", &pT_c);
    reader_c->AddVariable("nConst_c", &nConst_c);

    reader_c->AddSpectator("eta_c", &eta_c);
    reader_c->AddSpectator("phi_c", &phi_c);
    reader_c->AddSpectator("mass_c", &mass_c);
    reader_c->AddSpectator("label_c", &label_c);
    reader_c->AddSpectator("eventID_c", &eventID_c);
    reader_c->BookMVA("LogisticRegression", "dataset_c/weights/TMVARegression_LogisticRegression.weights.xml");

    TMVA::Reader* reader_s = new TMVA::Reader("!Color:!Silent");

    Float_t eventID_s, pT_s, eta_s, phi_s, mass_s, nConst_s = 0;
    Float_t label_s, score_s = 0;

    reader_s->AddVariable("pT_s", &pT_s);
    reader_s->AddVariable("nConst_s", &nConst_s);

    reader_s->AddSpectator("eta_s", &eta_s);
    reader_s->AddSpectator("phi_s", &phi_s);
    reader_s->AddSpectator("mass_s", &mass_s);
    reader_s->AddSpectator("label_s", &label_s);
    reader_s->AddSpectator("eventID_s", &eventID_s);
    reader_s->BookMVA("LogisticRegression", "dataset_s/weights/TMVARegression_LogisticRegression.weights.xml");

    //---------------------------------------------------------------------------------------------------------
    // Recuperação de TTrees de entrada 
    //---------------------------------------------------------------------------------------------------------

    TFile* inputFile_c = TFile::Open(inputFileName_c, "READ");
    TTree* signalTree_c = (TTree*)inputFile_c->Get("SignalTree_c");
    TTree* backgroundTree_c = (TTree*)inputFile_c->Get("BackgroundTree_c");

    signalTree_c->SetBranchAddress("pT_c", &pT_c);
    signalTree_c->SetBranchAddress("eta_c", &eta_c);
    signalTree_c->SetBranchAddress("phi_c", &phi_c);
    signalTree_c->SetBranchAddress("mass_c", &mass_c);
    signalTree_c->SetBranchAddress("nConst_c", &nConst_c);
    signalTree_c->SetBranchAddress("label_c", &label_c);
    signalTree_c->SetBranchAddress("eventID_c", &eventID_c);

    backgroundTree_c->SetBranchAddress("pT_c", &pT_c);
    backgroundTree_c->SetBranchAddress("eta_c", &eta_c);
    backgroundTree_c->SetBranchAddress("phi_c", &phi_c);
    backgroundTree_c->SetBranchAddress("mass_c", &mass_c);
    backgroundTree_c->SetBranchAddress("nConst_c", &nConst_c);
    backgroundTree_c->SetBranchAddress("label_c", &label_c);
    backgroundTree_c->SetBranchAddress("eventID_c", &eventID_c);

    TFile* inputFile_s = TFile::Open(inputFileName_s, "READ");
    TTree* signalTree_s = (TTree*)inputFile_s->Get("SignalTree_s");
    TTree* backgroundTree_s = (TTree*)inputFile_s->Get("BackgroundTree_s");

    signalTree_s->SetBranchAddress("pT_s", &pT_s);
    signalTree_s->SetBranchAddress("eta_s", &eta_s);
    signalTree_s->SetBranchAddress("phi_s", &phi_s);
    signalTree_s->SetBranchAddress("mass_s", &mass_s);
    signalTree_s->SetBranchAddress("nConst_s", &nConst_s);
    signalTree_s->SetBranchAddress("label_s", &label_s);
    signalTree_s->SetBranchAddress("eventID_s", &eventID_s);

    backgroundTree_s->SetBranchAddress("pT_s", &pT_s);
    backgroundTree_s->SetBranchAddress("eta_s", &eta_s);
    backgroundTree_s->SetBranchAddress("phi_s", &phi_s);
    backgroundTree_s->SetBranchAddress("mass_s", &mass_s);
    backgroundTree_s->SetBranchAddress("nConst_s", &nConst_s);
    backgroundTree_s->SetBranchAddress("label_s", &label_s);
    backgroundTree_s->SetBranchAddress("eventID_s", &eventID_s);

    //---------------------------------------------------------------------------------------------------------
    // Aplicação do modelo aos dados externos
    //---------------------------------------------------------------------------------------------------------


    // Charm block

    std::multimap<Int_t, TLorentzVector> jatos_c;

    for (Long64_t i = 0; i < signalTree_c->GetEntries(); ++i) // By construction, both signal and background TTrees for either c or s have the same lenght, that is, they refer to the same amount of events
    {
        signalTree_c->GetEntry(i);
        score_c = reader_c->EvaluateMVA("LogisticRegression");

        if (score_c > threshold) 
        {
            TLorentzVector jet;
            jet.SetPtEtaPhiM(pT_c, eta_c, phi_c, mass_c);
            jatos_c.insert({(Int_t)eventID_c, jet});
        }
    }
    
    // Obviously, ideally, no jet from this TTree should contribute to "jatos_c". 
    // However, since our model is, of course, not perfect, we shall consider those aberrations.

    for (Long64_t i = 0; i < backgroundTree_c->GetEntries(); ++i) 
    {
        backgroundTree_c->GetEntry(i);
        score_c = reader_c->EvaluateMVA("LogisticRegression");

        if (score_c > threshold) 
        {
            TLorentzVector jet;
            jet.SetPtEtaPhiM(pT_c, eta_c, phi_c, mass_c);
            jatos_c.insert({(Int_t)eventID_c, jet});
        }
    }

    // Strange block

        std::multimap<Int_t, TLorentzVector> jatos_s;

        for (Long64_t i = 0; i < signalTree_s->GetEntries(); ++i) 
        {
            signalTree_s->GetEntry(i);
            score_s = reader_s->EvaluateMVA("LogisticRegression");

            if (score_s > threshold) 
            {
                TLorentzVector jet;
                jet.SetPtEtaPhiM(pT_s, eta_s, phi_s, mass_s);
                jatos_s.insert({(Int_t)eventID_s, jet});
            }
        }

        for (Long64_t i = 0; i < backgroundTree_s->GetEntries(); ++i) 
        {
            backgroundTree_s->GetEntry(i);
            score_s = reader_s->EvaluateMVA("LogisticRegression");

            if (score_s > threshold) 
            {
                TLorentzVector jet;
                jet.SetPtEtaPhiM(pT_s, eta_s, phi_s, mass_s);
                jatos_s.insert({(Int_t)eventID_s, jet});
            }
        }

    //---------------------------------------------------------------------------------------------------------
    // Reconstrução combinatória da massa do W
    //---------------------------------------------------------------------------------------------------------
    
    TH1F* h_massW = new TH1F("h_massW", "Invariant Mass of W;M_{cs} [GeV];Events", 100, 0, 100);

    for (auto const& [eventID, _] : jatos_c) 
    {
        // Recuperar todos os jatos c e s com esse eventID
        auto range_c = jatos_c.equal_range(eventID);
        auto range_s = jatos_s.equal_range(eventID);

        // Se não houver jato s para esse evento, pula
        if (range_s.first == range_s.second) continue;

        for (auto it_c = range_c.first; it_c != range_c.second; ++it_c) // Combinatorial loop
        {
            for (auto it_s = range_s.first; it_s != range_s.second; ++it_s) 
            {
                TLorentzVector W = it_c->second + it_s->second;
                h_massW->Fill(W.M());
            }
        }
    }

    //---------------------------------------------------------------------------------------------------------
    // Plotar histogramas
    //---------------------------------------------------------------------------------------------------------

    TCanvas *c1 = new TCanvas("c1", "Reconstructed invariant mass distribution", 2500, 2500);
    c1->Divide(1, 1);

    c1->cd(1);
    h_massW->SetTitle("Reconstructed jet's invariant mass spectrum");
    h_massW->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    h_massW->GetYaxis()->SetTitle("Frequency");
    h_massW->DrawCopy();

    inputFile_c->Close();
    inputFile_s->Close();
    delete reader_c;
    delete reader_s;
}
