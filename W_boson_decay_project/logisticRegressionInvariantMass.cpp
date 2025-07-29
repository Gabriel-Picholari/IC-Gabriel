#include <map>
#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <algorithm>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMVA/Reader.h>
#include <TLorentzVector.h>

void logisticRegressionInvariantMass(const char* inputFileName_c, const char* inputFileName_s, const char* inputFileName_bkg, float threshold = 0.3) {

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
    reader_c->BookMVA("Likelihood", "dataset_c/weights/TMVAClassification_Likelihood.weights.xml");

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
    reader_s->BookMVA("Likelihood", "dataset_s/weights/TMVAClassification_Likelihood.weights.xml");
    
    TMVA::Reader* reader_bkg = new TMVA::Reader("!Color:!Silent");

    Float_t eventID_bkg, pT_bkg, eta_bkg, phi_bkg, mass_bkg, nConst_bkg = 0;
    Float_t label_bkg, score_bkg = 0;

    reader_bkg->AddVariable("pT_bkg", &pT_bkg);
    reader_bkg->AddVariable("nConst_bkg", &nConst_bkg);

    reader_bkg->AddSpectator("eta_bkg", &eta_bkg);
    reader_bkg->AddSpectator("phi_bkg", &phi_bkg);
    reader_bkg->AddSpectator("mass_bkg", &mass_bkg);
    reader_bkg->AddSpectator("label_bkg", &label_bkg);
    reader_bkg->AddSpectator("eventID_bkg", &eventID_bkg);
    reader_bkg->BookMVA("Likelihood", "dataset_bkg/weights/TMVAClassification_Likelihood.weights.xml");


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

    TFile* inputFile_bkg = TFile::Open(inputFileName_bkg, "READ");
    TTree* signalTree_bkg = (TTree*)inputFile_bkg->Get("SignalTree_bkg");
    TTree* backgroundTree_bkg = (TTree*)inputFile_bkg->Get("BackgroundTree_bkg");

    signalTree_bkg->SetBranchAddress("pT_bkg", &pT_bkg);
    signalTree_bkg->SetBranchAddress("eta_bkg", &eta_bkg);
    signalTree_bkg->SetBranchAddress("phi_bkg", &phi_bkg);
    signalTree_bkg->SetBranchAddress("mass_bkg", &mass_bkg);
    signalTree_bkg->SetBranchAddress("nConst_bkg", &nConst_bkg);
    signalTree_bkg->SetBranchAddress("label_bkg", &label_bkg);
    signalTree_bkg->SetBranchAddress("eventID_bkg", &eventID_bkg);

    backgroundTree_bkg->SetBranchAddress("pT_bkg", &pT_bkg);
    backgroundTree_bkg->SetBranchAddress("eta_bkg", &eta_bkg);
    backgroundTree_bkg->SetBranchAddress("phi_bkg", &phi_bkg);
    backgroundTree_bkg->SetBranchAddress("mass_bkg", &mass_bkg);
    backgroundTree_bkg->SetBranchAddress("nConst_bkg", &nConst_bkg);
    backgroundTree_bkg->SetBranchAddress("label_bkg", &label_bkg);
    backgroundTree_bkg->SetBranchAddress("eventID_bkg", &eventID_bkg);

    //---------------------------------------------------------------------------------------------------------
    // Histogramas
    //---------------------------------------------------------------------------------------------------------

    TH2F* signalS_set_correlation = new TH2F("signalS_corr", "Correlation between charmed and strange scores over the known strange signal data", 100, 0, 1, 100, 0, 1);
    TH2F* backgroundS_set_correlation = new TH2F("backgroundS_corr", "Correlation between charmed and strange scores over the known strange background data", 100, 0, 1, 100, 0, 1);
    TH2F* signalC_set_correlation = new TH2F("signalC_corr", "Correlation between charmed and strange scores over the known charmed signal data", 100, 0, 1, 100, 0, 1);
    TH2F* backgroundC_set_correlation = new TH2F("backgroundC_corr", "Correlation between charmed and strange scores over the known charmed background data", 100, 0, 1, 100, 0, 1);

    TH1F* h_massW = new TH1F("h_massW", "Invariant Mass of W;M_{cs} [GeV];Events", 100, 0, 100);

    //---------------------------------------------------------------------------------------------------------
    // Aplicação do modelo aos dados externos
    //---------------------------------------------------------------------------------------------------------

    // Charm block

    std::multimap<Int_t, TLorentzVector> jatos_c;
    std::multimap<Int_t, TLorentzVector> jatos_s;

    for (Long64_t i = 0; i < signalTree_c->GetEntries(); ++i) // By construction, both signal and background TTrees for either c or s have the same lenght, that is, they refer to the same amount of events
    {
        signalTree_c->GetEntry(i);
        score_c = reader_c->EvaluateMVA("Likelihood");
        score_s = reader_s->EvaluateMVA("Likelihood");
        signalC_set_correlation->Fill(score_c, score_s);

        score_bkg = reader_bkg->EvaluateMVA("Likelihood");

        if ( score_c == std::max({score_c, score_s, score_bkg}) ) 
        {
            TLorentzVector jet;
            jet.SetPtEtaPhiM(pT_c, eta_c, phi_c, mass_c);
            jatos_c.insert({(Int_t)eventID_c, jet});
        }
        else if ( score_s == std::max({score_c, score_s, score_bkg}) )
        {
            TLorentzVector jet;
            jet.SetPtEtaPhiM(pT_c, eta_c, phi_c, mass_c);
            jatos_s.insert({(Int_t)eventID_c, jet});
        }
    }
    
    // Obviously, ideally, no jet from this TTree should contribute to "jatos_c". 
    // However, since our model is, of course, not perfect, we shall consider those aberrations.

    for (Long64_t i = 0; i < backgroundTree_c->GetEntries(); ++i) 
    {
        backgroundTree_c->GetEntry(i);
        score_c = reader_c->EvaluateMVA("Likelihood");
        score_s = reader_s->EvaluateMVA("Likelihood");
        backgroundC_set_correlation->Fill(score_c, score_s);

        score_bkg = reader_bkg->EvaluateMVA("Likelihood");

        if ( score_c == std::max({score_c, score_s, score_bkg}) ) 
        {
            TLorentzVector jet;
            jet.SetPtEtaPhiM(pT_c, eta_c, phi_c, mass_c);
            jatos_c.insert({(Int_t)eventID_c, jet});
        }
        else if ( score_s == std::max({score_c, score_s, score_bkg}) )
        {
            TLorentzVector jet;
            jet.SetPtEtaPhiM(pT_c, eta_c, phi_c, mass_c);
            jatos_s.insert({(Int_t)eventID_c, jet});
        }
    }

    // Strange block

    for (Long64_t i = 0; i < signalTree_s->GetEntries(); ++i) 
    {
        signalTree_s->GetEntry(i);
        score_c = reader_c->EvaluateMVA("Likelihood");
        score_s = reader_s->EvaluateMVA("Likelihood");
        signalS_set_correlation->Fill(score_c, score_s);

        score_bkg = reader_bkg->EvaluateMVA("Likelihood");

        if ( score_c == std::max({score_c, score_s, score_bkg}) ) 
        {
            TLorentzVector jet;
            jet.SetPtEtaPhiM(pT_s, eta_s, phi_s, mass_s);
            jatos_c.insert({(Int_t)eventID_s, jet});
        }
        else if ( score_s == std::max({score_c, score_s, score_bkg}) )
        {
            TLorentzVector jet;
            jet.SetPtEtaPhiM(pT_s, eta_s, phi_s, mass_s);
            jatos_s.insert({(Int_t)eventID_s, jet});
        }
    }

    for (Long64_t i = 0; i < backgroundTree_s->GetEntries(); ++i) 
    {
        backgroundTree_s->GetEntry(i);
        score_c = reader_c->EvaluateMVA("Likelihood");
        score_s = reader_s->EvaluateMVA("Likelihood");
        backgroundS_set_correlation->Fill(score_c, score_s);

        score_bkg = reader_bkg->EvaluateMVA("Likelihood");

        if ( score_c == std::max({score_c, score_s, score_bkg}) ) 
        {
            TLorentzVector jet;
            jet.SetPtEtaPhiM(pT_s, eta_s, phi_s, mass_s);
            jatos_c.insert({(Int_t)eventID_s, jet});
        }
        else if ( score_s == std::max({score_c, score_s, score_bkg}) )
        {
            TLorentzVector jet;
            jet.SetPtEtaPhiM(pT_s, eta_s, phi_s, mass_s);
            jatos_s.insert({(Int_t)eventID_s, jet});
        }
    }
    
    // Background block ( I suppose ? ) I'd say it's necessary, once our model can fail while classifying background as signal as well
    // The result of such a failiure would, then, be a jet classified as charm or strange   
    /*
    for (Long64_t i = 0; i < signalTree_s->GetEntries(); ++i) 
    {
        signalTree_bkg->GetEntry(i);
        score_c = reader_c->EvaluateMVA("Likelihood");
        score_s = reader_s->EvaluateMVA("Likelihood");
        score_bkg = reader_bkg->EvaluateMVA("Likelihood");

        if ( score_c == std::max({score_c, score_s, score_bkg}) ) 
        {
            TLorentzVector jet;
            jet.SetPtEtaPhiM(pT_bkg, eta_bkg, phi_bkg, mass_bkg);
            jatos_c.insert({(Int_t)eventID_bkg, jet});
        }
        else if ( score_s == std::max({score_c, score_s, score_bkg}) )
        {
            TLorentzVector jet;
            jet.SetPtEtaPhiM(pT_bkg, eta_bkg, phi_bkg, mass_bkg);
            jatos_s.insert({(Int_t)eventID_bkg, jet});
        }
    }

    for (Long64_t i = 0; i < backgroundTree_s->GetEntries(); ++i) 
    {
        backgroundTree_bkg->GetEntry(i);
        score_c = reader_c->EvaluateMVA("Likelihood");
        score_s = reader_s->EvaluateMVA("Likelihood");
        score_bkg = reader_bkg->EvaluateMVA("Likelihood");

        if ( score_c == std::max({score_c, score_s, score_bkg}) ) 
        {
            TLorentzVector jet;
            jet.SetPtEtaPhiM(pT_bkg, eta_bkg, phi_bkg, mass_bkg);
            jatos_c.insert({(Int_t)eventID_bkg, jet});
        }
        else if ( score_s == std::max({score_c, score_s, score_bkg}) )
        {
            TLorentzVector jet;
            jet.SetPtEtaPhiM(pT_bkg, eta_bkg, phi_bkg, mass_bkg);
            jatos_s.insert({(Int_t)eventID_bkg, jet});
        }
    }
    */

    //---------------------------------------------------------------------------------------------------------
    // Reconstrução combinatória da massa do W
    //---------------------------------------------------------------------------------------------------------

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

    TCanvas *c2 = new TCanvas("c2", "Correlation distributions", 2500, 2500);
    c2->Divide(2, 2);

    c2->cd(1);
    signalS_set_correlation->SetTitle("Charm x Strange scores over strange signal data");
    signalS_set_correlation->GetXaxis()->SetTitle("Score c");
    signalS_set_correlation->GetYaxis()->SetTitle("Score s");
    signalS_set_correlation->DrawCopy();

    c2->cd(2);
    backgroundS_set_correlation->SetTitle("Charm x Strange scores over strange background data");
    backgroundS_set_correlation->GetXaxis()->SetTitle("Score c");
    backgroundS_set_correlation->GetYaxis()->SetTitle("Score s");
    backgroundS_set_correlation->DrawCopy();
    
    c2->cd(3);
    signalC_set_correlation->SetTitle("Charm x Strange scores over charm signal data");
    signalC_set_correlation->GetXaxis()->SetTitle("Score c");
    signalC_set_correlation->GetYaxis()->SetTitle("Score s");
    signalC_set_correlation->DrawCopy();

    c2->cd(4);
    backgroundC_set_correlation->SetTitle("Charm x Strange scores over charm background data");
    backgroundC_set_correlation->GetXaxis()->SetTitle("Score c");
    backgroundC_set_correlation->GetYaxis()->SetTitle("Score s");
    backgroundC_set_correlation->DrawCopy();


    inputFile_c->Close();
    inputFile_s->Close();
    delete reader_c;
    delete reader_s;
}
