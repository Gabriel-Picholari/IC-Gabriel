#include <map>
#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <iomanip> 
#include <iostream>
#include <algorithm>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMVA/Reader.h>
#include <TLorentzVector.h>

void printJets(const std::multimap<Int_t, TLorentzVector>& jatos, const std::string& nome) {
    std::cout << "-------------------------------------------\n";
    std::cout << "     Informações dos jatos: " << nome << "\n";
    std::cout << "-------------------------------------------\n";

    Int_t last_event = -1;
    for (const auto& [eventID, jet] : jatos) {
        if (eventID != last_event) {
            std::cout << "\nEventID: " << eventID << "\n";
            last_event = eventID;
        }

        std::cout << std::fixed << std::setprecision(4);
        std::cout << "Jet: pT = " << jet.Pt() << ", Mass = " << jet.M() << std::endl;
    }
}

void invariantMassDistribution_3var_BTD(const char* inputFileName_c, const char* inputFileName_s, float t_c = 0.7, float t_s = 0.1) 
{

    //---------------------------------------------------------------------------------------------------------
    // Criação dos objetos Readers para leitura de resultados referentes tanto ao c como ao s
    //---------------------------------------------------------------------------------------------------------

    Float_t eventID, pT, eta, phi, mass, nConst, maxRho, nRho = 0;
    Float_t score_c, score_s, score_bkg = 0;
    Float_t label, score = 0;

    TMVA::Reader* reader_c = new TMVA::Reader("!Color:!Silent");

    reader_c->AddVariable("pT_c", &pT);
    reader_c->AddVariable("nRho_c", &nRho);
    reader_c->AddVariable("nConst_c", &nConst);
    //reader_c->AddVariable("maxRho_c", &maxRho); Discontinued

    reader_c->AddSpectator("eta_c", &eta);
    reader_c->AddSpectator("phi_c", &phi);
    reader_c->AddSpectator("mass_c", &mass);
    reader_c->AddSpectator("label_c", &label);
    reader_c->AddSpectator("eventID_c", &eventID);
    reader_c->BookMVA("GradBoost", "dataset_c_3var/weights/TMVAClassification_GradBoost.weights.xml");

    TMVA::Reader* reader_s = new TMVA::Reader("!Color:!Silent");

    reader_s->AddVariable("pT_s", &pT);
    reader_s->AddVariable("nRho_s", &nRho);
    reader_s->AddVariable("nConst_s", &nConst);
    //reader_s->AddVariable("maxRho_s", &maxRho); Discontinued

    reader_s->AddSpectator("eta_s", &eta);
    reader_s->AddSpectator("phi_s", &phi);
    reader_s->AddSpectator("mass_s", &mass);
    reader_s->AddSpectator("label_s", &label);
    reader_s->AddSpectator("eventID_s", &eventID);
    reader_s->BookMVA("GradBoost", "dataset_s_3var/weights/TMVAClassification_GradBoost.weights.xml");

    //---------------------------------------------------------------------------------------------------------
    // Recuperação de TTrees de entrada 
    //---------------------------------------------------------------------------------------------------------

    TFile* inputFile_c = TFile::Open(inputFileName_c, "READ");
    TTree* signalTree_c = (TTree*)inputFile_c->Get("SignalTree_c");
    TTree* backgroundTree_c = (TTree*)inputFile_c->Get("BackgroundTree_c");

    signalTree_c->SetBranchAddress("pT_c", &pT);
    signalTree_c->SetBranchAddress("eta_c", &eta);
    signalTree_c->SetBranchAddress("phi_c", &phi);
    signalTree_c->SetBranchAddress("mass_c", &mass);
    signalTree_c->SetBranchAddress("nConst_c", &nConst);
    signalTree_c->SetBranchAddress("nRho_c", &nRho);
    //signalTree_c->SetBranchAddress("maxRho_c", &maxRho); Discontinued
    signalTree_c->SetBranchAddress("label_c", &label);
    signalTree_c->SetBranchAddress("eventID_c", &eventID);

    backgroundTree_c->SetBranchAddress("pT_c", &pT);
    backgroundTree_c->SetBranchAddress("eta_c", &eta);
    backgroundTree_c->SetBranchAddress("phi_c", &phi);
    backgroundTree_c->SetBranchAddress("mass_c", &mass);
    backgroundTree_c->SetBranchAddress("nConst_c", &nConst);
    backgroundTree_c->SetBranchAddress("nRho_c", &nRho);
    //backgroundTree_c->SetBranchAddress("maxRho_c", &maxRho); Discontinued
    backgroundTree_c->SetBranchAddress("label_c", &label);
    backgroundTree_c->SetBranchAddress("eventID_c", &eventID);

    TFile* inputFile_s = TFile::Open(inputFileName_s, "READ");
    TTree* signalTree_s = (TTree*)inputFile_s->Get("SignalTree_s");
    TTree* backgroundTree_s = (TTree*)inputFile_s->Get("BackgroundTree_s");

    signalTree_s->SetBranchAddress("pT_s", &pT);
    signalTree_s->SetBranchAddress("eta_s", &eta);
    signalTree_s->SetBranchAddress("phi_s", &phi);
    signalTree_s->SetBranchAddress("mass_s", &mass);
    signalTree_s->SetBranchAddress("nConst_s", &nConst);
    signalTree_s->SetBranchAddress("nRho_s", &nRho);
    //signalTree_s->SetBranchAddress("maxRho_s", &maxRho); Discontinued
    signalTree_s->SetBranchAddress("label_s", &label);
    signalTree_s->SetBranchAddress("eventID_s", &eventID);

    backgroundTree_s->SetBranchAddress("pT_s", &pT);
    backgroundTree_s->SetBranchAddress("eta_s", &eta);
    backgroundTree_s->SetBranchAddress("phi_s", &phi);
    backgroundTree_s->SetBranchAddress("mass_s", &mass);
    backgroundTree_s->SetBranchAddress("nConst_s", &nConst);
    backgroundTree_s->SetBranchAddress("nRho_s", &nRho);
    //backgroundTree_s->SetBranchAddress("maxRho_s", &maxRho); Discontinued
    backgroundTree_s->SetBranchAddress("label_s", &label);
    backgroundTree_s->SetBranchAddress("eventID_s", &eventID);

    //---------------------------------------------------------------------------------------------------------
    // Histogramas
    //---------------------------------------------------------------------------------------------------------

    TH2F* signalS_set_correlation = new TH2F("signalS_corr", "Correlation between charmed and strange scores over the known strange signal data", 100, -1, 1, 100, -1, 1);
    TH2F* backgroundS_set_correlation = new TH2F("backgroundS_corr", "Correlation between charmed and strange scores over the known strange background data", 100, -1, 1, 100, -1, 1);
    TH2F* signalC_set_correlation = new TH2F("signalC_corr", "Correlation between charmed and strange scores over the known charmed signal data", 100, -1, 1, 100, -1, 1);
    TH2F* backgroundC_set_correlation = new TH2F("backgroundC_corr", "Correlation between charmed and strange scores over the known charmed background data", 100, -1, 1, 100, -1, 1);

    TH1F* h_massW = new TH1F("h_massW", "Invariant Mass of W;M_{cs} [GeV];Events", 100, 0, 100);

    //---------------------------------------------------------------------------------------------------------
    // Aplicação do modelo aos dados externos
    //---------------------------------------------------------------------------------------------------------

    // Charm block

    std::multimap<Int_t, TLorentzVector> jatos_c;
    std::multimap<Int_t, TLorentzVector> jatos_s;

    auto iguais = [](const TLorentzVector& a, const TLorentzVector& b) 
    {
        return (std::fabs(a.Pt() - b.Pt())  < 1e-4 && std::fabs(a.M() - b.M())   < 1e-6);
    };

    for (Long64_t i = 0; i < signalTree_c->GetEntries(); ++i) // By construction, both signal and background TTrees for either c or s have the same lenght, that is, they refer to the same amount of events
    {
        signalTree_c->GetEntry(i);
        score_c = reader_c->EvaluateMVA("GradBoost");
        score_s = reader_s->EvaluateMVA("GradBoost");
        if (score_c < t_c && score_s < t_s) continue;

        // Caso 1: passou apenas no corte do charm
        if (score_c >= t_c && score_s < t_s) 
        {
            TLorentzVector jet;
            jet.SetPtEtaPhiM(pT, eta, phi, mass);

            auto range = jatos_c.equal_range((Int_t)eventID);

            bool existe = false;
            for (auto it = range.first; it != range.second; ++it) 
            {
                if (iguais(it->second, jet)) 
                {
                    existe = true;
                    break;
                }
            }

            if (!existe) 
            {
                jatos_c.insert({(Int_t)eventID, jet});
            }
        }

        // Caso 2: passou apenas no corte do strange
        else if (score_s >= t_s && score_c < t_c) 
        {
            TLorentzVector jet;
            jet.SetPtEtaPhiM(pT, eta, phi, mass);

            auto range = jatos_s.equal_range((Int_t)eventID);

            bool existe = false;
            for (auto it = range.first; it != range.second; ++it) 
            {
                if (iguais(it->second, jet)) 
                {
                    existe = true;
                    break;
                }
            }
            if (!existe) 
            {
                jatos_s.insert({(Int_t)eventID, jet});
            }
        }

        // Caso 3: ambos passam os cortes
        else if (score_c >= t_c && score_s >= t_s) 
        {
            if (score_c > score_s)
            {
                TLorentzVector jet;
                jet.SetPtEtaPhiM(pT, eta, phi, mass);

                auto range = jatos_c.equal_range((Int_t)eventID);

                bool existe = false;
                for (auto it = range.first; it != range.second; ++it) 
                {
                    if (iguais(it->second, jet)) 
                    {
                        existe = true;
                        break;
                    }
                }
                if (!existe) 
                {
                    jatos_c.insert({(Int_t)eventID, jet});
                }
            }
            else
            {
                TLorentzVector jet;
                jet.SetPtEtaPhiM(pT, eta, phi, mass);

                auto range = jatos_s.equal_range((Int_t)eventID);

                bool existe = false;
                for (auto it = range.first; it != range.second; ++it) 
                {
                    if (iguais(it->second, jet)) 
                    {
                        existe = true;
                        break;
                    }
                }
                if (!existe) 
                {
                    jatos_s.insert({(Int_t)eventID, jet});
                }
            }
        }
        //std::cout << "Signal for C - Score charm: " << score_c << ", Score S: " << score_s << std::endl;
        signalC_set_correlation->Fill(score_c, score_s);

    }
    
    // Obviously, ideally, no jet from this TTree should contribute to "jatos_c". 
    // However, since our model is, of course, not perfect, we shall consider those aberrations.

    for (Long64_t i = 0; i < backgroundTree_c->GetEntries(); ++i) 
    {
        backgroundTree_c->GetEntry(i);
        score_c = reader_c->EvaluateMVA("GradBoost");
        score_s = reader_s->EvaluateMVA("GradBoost");
        if (score_c < t_c && score_s < t_s) continue;
        
        // Caso 1: passou apenas no corte do charm
        if (score_c >= t_c && score_s < t_s) 
        {
            TLorentzVector jet;
            jet.SetPtEtaPhiM(pT, eta, phi, mass);

            auto range = jatos_c.equal_range((Int_t)eventID);

            bool existe = false;
            for (auto it = range.first; it != range.second; ++it) 
            {
                if (iguais(it->second, jet)) 
                {
                    existe = true;
                    break;
                }
            }

            if (!existe) 
            {
                jatos_c.insert({(Int_t)eventID, jet});
            }
        }

        // Caso 2: passou apenas no corte do strange
        else if (score_s >= t_s && score_c < t_c) 
        {
            TLorentzVector jet;
            jet.SetPtEtaPhiM(pT, eta, phi, mass);

            auto range = jatos_s.equal_range((Int_t)eventID);

            bool existe = false;
            for (auto it = range.first; it != range.second; ++it) 
            {
                if (iguais(it->second, jet)) 
                {
                    existe = true;
                    break;
                }
            }
            if (!existe) 
            {
                jatos_s.insert({(Int_t)eventID, jet});
            }
        }

        // Caso 3: ambos passam os cortes
        else if (score_c >= t_c && score_s >= t_s) 
        {
            if (score_c > score_s)
            {
                TLorentzVector jet;
                jet.SetPtEtaPhiM(pT, eta, phi, mass);

                auto range = jatos_c.equal_range((Int_t)eventID);

                bool existe = false;
                for (auto it = range.first; it != range.second; ++it) 
                {
                    if (iguais(it->second, jet)) 
                    {
                        existe = true;
                        break;
                    }
                }
                if (!existe) 
                {
                    jatos_c.insert({(Int_t)eventID, jet});
                }
            }
            else
            {
                TLorentzVector jet;
                jet.SetPtEtaPhiM(pT, eta, phi, mass);

                auto range = jatos_s.equal_range((Int_t)eventID);

                bool existe = false;
                for (auto it = range.first; it != range.second; ++it) 
                {
                    if (iguais(it->second, jet)) 
                    {
                        existe = true;
                        break;
                    }
                }
                if (!existe) 
                {
                    jatos_s.insert({(Int_t)eventID, jet});
                }
            }
        }
        //std::cout << "Backgound for C - Score charm: " << score_c << ", Score S: " << score_s << std::endl;
        backgroundC_set_correlation->Fill(score_c, score_s);
    }

    // Strange block

    for (Long64_t i = 0; i < signalTree_s->GetEntries(); ++i) 
    {
        signalTree_s->GetEntry(i);
        score_c = reader_c->EvaluateMVA("GradBoost");
        score_s = reader_s->EvaluateMVA("GradBoost");
        if (score_c < t_c && score_s < t_s) continue;
        
        // Caso 1: passou apenas no corte do charm
        if (score_c >= t_c && score_s < t_s) 
        {
            TLorentzVector jet;
            jet.SetPtEtaPhiM(pT, eta, phi, mass);

            auto range = jatos_c.equal_range((Int_t)eventID);

            bool existe = false;
            for (auto it = range.first; it != range.second; ++it) 
            {
                if (iguais(it->second, jet)) 
                {
                    existe = true;
                    break;
                }
            }

            if (!existe) 
            {
                jatos_c.insert({(Int_t)eventID, jet});
            }
        }

        // Caso 2: passou apenas no corte do strange
        else if (score_s >= t_s && score_c < t_c) 
        {
            TLorentzVector jet;
            jet.SetPtEtaPhiM(pT, eta, phi, mass);

            auto range = jatos_s.equal_range((Int_t)eventID);

            bool existe = false;
            for (auto it = range.first; it != range.second; ++it) 
            {
                if (iguais(it->second, jet)) 
                {
                    existe = true;
                    break;
                }
            }
            if (!existe) 
            {
                jatos_s.insert({(Int_t)eventID, jet});
            }
        }

        // Caso 3: ambos passam os cortes
        else if (score_c >= t_c && score_s >= t_s) 
        {
            if (score_c > score_s)
            {
                TLorentzVector jet;
                jet.SetPtEtaPhiM(pT, eta, phi, mass);

                auto range = jatos_c.equal_range((Int_t)eventID);

                bool existe = false;
                for (auto it = range.first; it != range.second; ++it) 
                {
                    if (iguais(it->second, jet)) 
                    {
                        existe = true;
                        break;
                    }
                }
                if (!existe) 
                {
                    jatos_c.insert({(Int_t)eventID, jet});
                }
            }
            else
            {
                TLorentzVector jet;
                jet.SetPtEtaPhiM(pT, eta, phi, mass);

                auto range = jatos_s.equal_range((Int_t)eventID);

                bool existe = false;
                for (auto it = range.first; it != range.second; ++it) 
                {
                    if (iguais(it->second, jet)) 
                    {
                        existe = true;
                        break;
                    }
                }
                if (!existe) 
                {
                    jatos_s.insert({(Int_t)eventID, jet});
                }
            }
        }
        //std::cout << "Signal for S - Score charm: " << score_c << ", Score S: " << score_s << std::endl;
        signalS_set_correlation->Fill(score_c, score_s);
    }
    for (Long64_t i = 0; i < backgroundTree_s->GetEntries(); ++i) 
    {
        backgroundTree_s->GetEntry(i);
        score_c = reader_c->EvaluateMVA("GradBoost");
        score_s = reader_s->EvaluateMVA("GradBoost");
        if (score_c < t_c && score_s < t_s) continue;

        // Caso 1: passou apenas no corte do charm
        if (score_c >= t_c && score_s < t_s) 
        {
            TLorentzVector jet;
            jet.SetPtEtaPhiM(pT, eta, phi, mass);

            auto range = jatos_c.equal_range((Int_t)eventID);

            bool existe = false;
            for (auto it = range.first; it != range.second; ++it) 
            {
                if (iguais(it->second, jet)) 
                {
                    existe = true;
                    break;
                }
            }

            if (!existe) 
            {
                jatos_c.insert({(Int_t)eventID, jet});
            }
        }

        // Caso 2: passou apenas no corte do strange
        else if (score_s >= t_s && score_c < t_c) 
        {
            TLorentzVector jet;
            jet.SetPtEtaPhiM(pT, eta, phi, mass);

            auto range = jatos_s.equal_range((Int_t)eventID);

            bool existe = false;
            for (auto it = range.first; it != range.second; ++it) 
            {
                if (iguais(it->second, jet)) 
                {
                    existe = true;
                    break;
                }
            }
            if (!existe) 
            {
                jatos_s.insert({(Int_t)eventID, jet});
            }
        }

        // Caso 3: ambos passam os cortes
        else if (score_c >= t_c && score_s >= t_s) 
        {
            if (score_c > score_s)
            {
                TLorentzVector jet;
                jet.SetPtEtaPhiM(pT, eta, phi, mass);

                auto range = jatos_c.equal_range((Int_t)eventID);

                bool existe = false;
                for (auto it = range.first; it != range.second; ++it) 
                {
                    if (iguais(it->second, jet)) 
                    {
                        existe = true;
                        break;
                    }
                }
                if (!existe) 
                {
                    jatos_c.insert({(Int_t)eventID, jet});
                }
            }
            else
            {
                TLorentzVector jet;
                jet.SetPtEtaPhiM(pT, eta, phi, mass);

                auto range = jatos_s.equal_range((Int_t)eventID);

                bool existe = false;
                for (auto it = range.first; it != range.second; ++it) 
                {
                    if (iguais(it->second, jet)) 
                    {
                        existe = true;
                        break;
                    }
                }
                if (!existe) 
                {
                    jatos_s.insert({(Int_t)eventID, jet});
                }
            }
        }
        //std::cout << "Background for S - Score charm: " << score_c << ", Score S: " << score_s << std::endl;
        backgroundS_set_correlation->Fill(score_c, score_s);
    }


    //printJets(jatos_c, "Charm");
    //printJets(jatos_s, "Strange");
    
    //---------------------------------------------------------------------------------------------------------
    // Reconstrução combinatória da massa do W
    //---------------------------------------------------------------------------------------------------------
    
    for (auto it = jatos_c.begin(); it != jatos_c.end(); it = jatos_c.upper_bound(it->first)) 
    {
        Int_t eventID = it->first;

        auto range_c = jatos_c.equal_range(eventID);
        auto range_s = jatos_s.equal_range(eventID);

        if (range_c.first == range_c.second) continue; // não deveria acontecer aqui, mas ok
        if (range_s.first == range_s.second) continue; // sem s: não há par

        for (auto it_c = range_c.first; it_c != range_c.second; ++it_c) 
        {
            for (auto it_s = range_s.first; it_s != range_s.second; ++it_s) {
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
