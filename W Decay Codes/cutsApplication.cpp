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
    Float_t match_R = 0.1;

    // Criacao dos TLorentzVector

    TLorentzVector vec_unfiltered(0,0,0,0);
    
    TLorentzVector vec_pTCut(0,0,0,0);
    TLorentzVector vec_nConstCut(0,0,0,0);
    TLorentzVector vec_second_pTCut(0,0,0,0);
    TLorentzVector vec_pT_and_nConstCut(0,0,0,0);

    //---------------------------------------------------------------------------------------------------------
    // Inicializacao dos histogramas
    //---------------------------------------------------------------------------------------------------------

    TH1F *noCut_IMS = new TH1F("h1", "Unfiltered Jet's Invariant Mass Spectrum", 100, 0, 100);
    TH1F *pTCut_IMS = new TH1F("h2", "Jet's Invariant Mass Sprectrum with cut in p_{T}", 100, 0, 100);
    TH1F *nConstCut_IMS = new TH1F("h3", "Jet's Invariant Mass Sprectrum with a cut in the number of constituents", 100, 0, 100);
    TH1F *second_pTCut_IMS = new TH1F("h4", "Jet's Invariant Mass Sprectrum with an arbitrary cut in p_{T}", 100, 0, 100);
    TH1F *pT_and_nConstCut_IMS = new TH1F("h5", "Jet's Invariant Mass Sprectrum with cuts both in the number of constituents and p_{T}", 100, 0, 100);

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

        vec_unfiltered.SetPtEtaPhiM(10,0,0,(3.141592));

        vec_pTCut.SetPtEtaPhiM(10,0,0,(3.141592));
        vec_nConstCut.SetPtEtaPhiM(10,0,0,(3.141592));
        vec_second_pTCut.SetPtEtaPhiM(10,0,0,(3.141592));
        vec_pT_and_nConstCut.SetPtEtaPhiM(10,0,0,(3.141592));
        

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

        for (const fastjet::PseudoJet& jet : jets)
        {
        jetPt = jet.pt();
        jetEta = jet.eta();
        jetPhi = jet.phi();
        jetMass = jet.m();
        jetPx = jet.px();
        jetPy = jet.py();
        jetPz = jet.pz();
        jetE = jet.E();
        jetNConst = jet.constituents().size();
        pT_LeadConst = 0.0;
        for (const fastjet::PseudoJet &constituent : jet.constituents())
        {
            if (constituent.pt() > pT_LeadConst)
            {
            pT_LeadConst = constituent.pt();
            }
        }

            //---------------------------------------------------------------------------------------------------------
            // Preenchimento dos TLorentzVectors dos jatos do evento
            //---------------------------------------------------------------------------------------------------------
            
            vec_unfiltered = TLorentzVector(jetPx, jetPy, jetPz, jetE);
            
            if ( jetNConst > 2.5 )
            {
                vec_nConstCut = TLorentzVector(jetPx, jetPy, jetPz, jetE);
            }

            if ( jetPt < 1.5 ) continue;
            vec_pTCut = TLorentzVector(jetPx, jetPy, jetPz, jetE);
            
            if ( jetNConst > 2.5 )
            {
                vec_pT_and_nConstCut = TLorentzVector(jetPx, jetPy, jetPz, jetE);
            }

            if (jetPt < 10 ) continue;
            vec_second_pTCut = TLorentzVector(jetPx, jetPy, jetPz, jetE);


        } // End of individual jets creation
        
        const Float_t tolerance = 1e-5; // IMS stands for Invariant Mass Spectrum

        if ( fabs( vec_unfiltered.M() - 3.141592) > tolerance ) 
        {
        noCut_IMS->Fill( vec_unfiltered.M() ); 
        }

        if ( fabs( vec_nConstCut.M() - 3.141592) > tolerance ) 
        {
        nConstCut_IMS->Fill( vec_nConstCut.M() ); 
        }

        if ( fabs( vec_pTCut.M() - 3.141592) > tolerance ) 
        {
        pTCut_IMS->Fill( vec_pTCut.M() ); 
        }

        if ( fabs( vec_pT_and_nConstCut.M() - 3.141592) > tolerance ) 
        {
        pT_and_nConstCut_IMS->Fill( vec_pT_and_nConstCut.M() ); 
        }

        if ( fabs( vec_second_pTCut.M() - 3.141592) > tolerance ) 
        {
        second_pTCut_IMS->Fill( vec_second_pTCut.M() ); 
        }


        particles_fastjet.clear();
        jets.clear();
        jets_array->Clear();
        quarks->Clear();
    }

    TCanvas *c1 = new TCanvas("c1", "Invariant Mass Spectrum", 2500, 2500);
    c1->Divide(1, 1);

    c1->cd(1);
    noCut_IMS->SetTitle("Unfiltered Jet's Invariant Mass Spectrum");
    noCut_IMS->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    noCut_IMS->GetYaxis()->SetTitle("Frequency");
    noCut_IMS->Draw();

    TCanvas *c2 = new TCanvas("c2", "Invariant Mass Spectrum variants", 2500, 2500);
    c2->Divide(2, 2);

    c2->cd(1);
    pTCut_IMS->SetTitle("Jet's Invariant Mass Sprectrum with cut in p_{T}");
    pTCut_IMS->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    pTCut_IMS->GetYaxis()->SetTitle("Frequency");
    pTCut_IMS->Draw();

    c2->cd(2);
    nConstCut_IMS->SetTitle("Jet's Invariant Mass Sprectrum with a cut in the number of constituents");
    nConstCut_IMS->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    nConstCut_IMS->GetYaxis()->SetTitle("Frequency");
    nConstCut_IMS->Draw();

    c2->cd(3);
    second_pTCut_IMS->SetTitle("Jet's Invariant Mass Sprectrum with an arbitrary cut in p_{T}");
    second_pTCut_IMS->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    second_pTCut_IMS->GetYaxis()->SetTitle("Frequency");
    second_pTCut_IMS->Draw();

    c2->cd(4);
    pT_and_nConstCut_IMS->SetTitle("Jet's Invariant Mass Sprectrum with cuts both in the number of constituents and p_{T}");
    pT_and_nConstCut_IMS->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    pT_and_nConstCut_IMS->GetYaxis()->SetTitle("Frequency");
    pT_and_nConstCut_IMS->Draw();
}

