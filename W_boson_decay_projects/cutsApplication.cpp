#include <cmath>
#include <string>
#include "TH1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TMath.h"
#include "MyJet.h"
#include "TFile.h"
#include <iostream>
#include "TSystem.h"
#include "MyQuark.h"
#include "TCanvas.h"
#include "TPythia8.h"
#include "TParticle.h"
#include "TClonesArray.h"
#include <TLorentzVector.h>
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

void cutsApplication(const char* fileName)
{
    gSystem->Load("libEG");
    gSystem->Load("libEGPythia8");

    //---------------------------------------------------------------------------------------------------------
    //Inicializacao das variaveis
    //---------------------------------------------------------------------------------------------------------
    
    Float_t fpPt, fpEta, fpPhi, fpE, fpPx, fpPy, fpPz, fpMass = 0;
    Float_t jetPt, jetEta, jetPhi, jetE, jetPx, jetPy, jetPz, jetMass, jetNConst, pT_LeadConst = 0;

    //---------------------------------------------------------------------------------------------------------
    // Inicializacao dos histogramas
    //---------------------------------------------------------------------------------------------------------

    TH1F *noCut_IMS = new TH1F("h1", "Unfiltered Jet's Invariant Mass Spectrum", 100, 0, 100);
    TH1F *pTCut_IMS = new TH1F("h2", "Jet's Invariant Mass Sprectrum with cut in p_{T}", 100, 0, 100);
    TH1F *nConstCut_IMS = new TH1F("h3", "Jet's Invariant Mass Sprectrum with a cut in the number of constituents", 100, 0, 100);
    TH1F *second_pTCut_IMS = new TH1F("h4", "Jet's Invariant Mass Sprectrum with an arbitrary cut in p_{T}", 100, 0, 100);
    TH1F *pT_and_nConstCut_IMS = new TH1F("h5", "Jet's Invariant Mass Sprectrum with cuts both in the number of constituents and p_{T}", 100, 0, 100);
    TH1F *second_pTCut_and_nConst_IMS = new TH1F("h6", "Jet's Invariant Mass Sprectrum with secondary cuts both in the number of constituents and p_{T}", 100, 0, 100);

    TH1F *hist_jetPt = new TH1F("h7","Jet p_{T}", 100, 0, 100);
    TH1F *hist_jetNConst = new TH1F("h8","Jet number of constituents", 100, 0, 100);

    //---------------------------------------------------------------------------------------------------------
    // Inicializacoes e configuracoes do FastJet:
    //---------------------------------------------------------------------------------------------------------

    Float_t jetR = 0.5;

    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, jetR);

    std::vector<fastjet::PseudoJet> particles_fastjet;
    std::vector<fastjet::PseudoJet> jets;

    //---------------------------------------------------------------------------------------------------------
    // Inicializacao do arquivo.root e das TTrees
    //---------------------------------------------------------------------------------------------------------

    TFile *file = TFile::Open(fileName, "READ");
    TTree *ttree = dynamic_cast<TTree *>(file->Get("W decay TTree"));

    TClonesArray *jets_array = new TClonesArray("MyJet");
    TClonesArray *quarks = new TClonesArray("MyQuark");

    ttree->SetBranchAddress("jets_array", &jets_array);
    ttree->SetBranchAddress("quarks", &quarks);

    Long64_t ne = ttree->GetEntries();
    
    //---------------------------------------------------------------------------------------------------------
    // Loop equivalente ao loop de eventos
    //---------------------------------------------------------------------------------------------------------

    for ( Long64_t ni = 0; ni < ne; ni++)
    {
        ttree->GetEntry(ni);

        //---------------------------------------------------------------------------------------------------------
        // Loop equivalente ao loop de eventos
        //---------------------------------------------------------------------------------------------------------

        for (Int_t nj = 0; nj < jets_array->GetEntries(); nj++)
        {
        MyJet *fp = static_cast<MyJet *>(jets_array->At(nj));

        fpPx = fp->fPx;
        fpPy = fp->fPy;
        fpPz = fp->fPz;
        fpE = fp->fE;

        fastjet::PseudoJet particle(fpPx, fpPy, fpPz, fpE);
        particles_fastjet.push_back(particle);
        }

        fastjet::ClusterSequence clusterSeq(particles_fastjet, jet_def);
        jets = clusterSeq.inclusive_jets();
    

        for (Int_t m = 0; m < jets.size() - 1; m++)
        {

            for (Int_t k = m+1; k < jets.size(); k++)
            {

                fastjet::PseudoJet jet_m = jets[m];
                fastjet::PseudoJet jet_k = jets[k];

                TLorentzVector vec_m(0,0,0,0);
                vec_m.SetPxPyPzE(jet_m.px(), jet_m.py(), jet_m.pz(), jet_m.e());
                
                Float_t mpT = vec_m.Pt();
                Float_t mNConst = jet_m.constituents().size();

                hist_jetPt->Fill(mpT);
                hist_jetNConst->Fill(mNConst);

                TLorentzVector vec_k(0,0,0,0);
                vec_k.SetPxPyPzE(jet_k.px(), jet_k.py(), jet_k.pz(), jet_k.e());

                TLorentzVector vec_mom = vec_m + vec_k;

                jetPt = vec_mom.Pt();
                jetEta = vec_mom.Eta();
                jetPhi = vec_mom.Phi();
                jetMass = vec_mom.M();
                jetPx = vec_mom.Px();
                jetPy = vec_mom.Py();
                jetPz = vec_mom.Pz();
                jetE = vec_mom.E();
                
                jetNConst = jet_m.constituents().size() + jet_k.constituents().size();
                
                /* 
                pT_LeadConst = 0.0;
                for (const fastjet::PseudoJet &constituent : jet_m.constituents())
                {
                    if (constituent.pt() > pT_LeadConst)
                    {
                        pT_LeadConst = constituent.pt();
                    }
                }
                */
                
                noCut_IMS->Fill( vec_mom.M() );

                if ( jetNConst > 3 ) nConstCut_IMS->Fill( vec_mom.M() );

                if ( jetPt < 1.5 ) continue;
                pTCut_IMS->Fill( vec_mom.M() ); 
                
                if ( jetNConst > 3 ) pT_and_nConstCut_IMS->Fill( vec_mom.M() ); 

                if (jetPt < 10 ) continue;
                second_pTCut_IMS->Fill( vec_mom.M() );

                if ( jetNConst > 3 ) second_pTCut_and_nConst_IMS->Fill( vec_mom.M() );
                
            }
        }
 
        particles_fastjet.clear();
        jets.clear();
        jets_array->Clear();
        quarks->Clear();
    }
    /* 
    TCanvas *c1 = new TCanvas("c1", "Invariant Mass Spectrum", 2500, 2500);
    c1->Divide(1, 2);

    c1->cd(1);
    noCut_IMS->SetTitle("Invariant Mass Spectrum of Unfiltered Jets");
    noCut_IMS->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    noCut_IMS->GetYaxis()->SetTitle("Frequency");
    noCut_IMS->Draw();
    */

    TCanvas *c2 = new TCanvas("c2", "Variants of Invariant Mass Spectrum", 2500, 2500);
    c2->Divide(2, 2);

    c2->cd(1);
    noCut_IMS->SetTitle("Invariant Mass Spectrum of Unfiltered Jets");
    noCut_IMS->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    noCut_IMS->GetYaxis()->SetTitle("Frequency");
    noCut_IMS->Draw();

    c2->cd(2);
    pTCut_IMS->SetTitle("Invariant Mass Spectrum of Jets with p_{T} Cut");
    pTCut_IMS->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    pTCut_IMS->GetYaxis()->SetTitle("Frequency");
    pTCut_IMS->Draw();

    c2->cd(3);
    nConstCut_IMS->SetTitle("Invariant Mass Spectrum of Jets with Cut on Number of Constituents");
    nConstCut_IMS->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    nConstCut_IMS->GetYaxis()->SetTitle("Frequency");
    nConstCut_IMS->Draw();

    /* 
    c2->cd(3);
    second_pTCut_IMS->SetTitle("Invariant Mass Spectrum of Jets with Arbitrary p_{T} Cut");
    second_pTCut_IMS->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    second_pTCut_IMS->GetYaxis()->SetTitle("Frequency");
    second_pTCut_IMS->Draw();
    */

    c2->cd(4);
    pT_and_nConstCut_IMS->SetTitle("Invariant Mass Spectrum of Jets with Cuts on Number of Constituents and p_{T}");
    pT_and_nConstCut_IMS->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    pT_and_nConstCut_IMS->GetYaxis()->SetTitle("Frequency");
    pT_and_nConstCut_IMS->Draw();

    TCanvas *c3 = new TCanvas("c3", "Simple histograms", 2500, 2500);
    c3->Divide(1, 2);

    c3->cd(1);
    hist_jetPt->SetTitle("Jet p_{T}");
    hist_jetPt->GetXaxis()->SetTitle("Pt [GeV/c]");
    hist_jetPt->GetYaxis()->SetTitle("Frequency");
    hist_jetPt->Draw();

    c3->cd(2);
    hist_jetNConst->SetTitle("Jet Number of Constituents");
    hist_jetNConst->GetXaxis()->SetTitle("Number of Constituents");
    hist_jetNConst->GetYaxis()->SetTitle("Frequency");
    hist_jetNConst->Draw();

}

