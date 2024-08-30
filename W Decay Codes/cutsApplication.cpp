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

    //---------------------------------------------------------------------------------------------------------
    // Inicializacao dos histogramas
    //---------------------------------------------------------------------------------------------------------

    TH1F *unfiltered_Jets_InvariantMass = new TH1F("h1", "Massa invariante de jatos não filtrados em [GeV/c^{2}]", 100, 0, 100);
    TH1F *cutPt_Jets_InvariantMass = new TH1F("h2", "Massa invariante de jatos com corte de pT de 1.5 [GeV/c] em [GeV/c^{2}]", 100, 0, 100);
    TH1F *arbitraryPtCut_Jets_InvariantMass = new TH1F("h2", "Massa invariante de jatos com corte de pT [GeV/c] arbitrário em [GeV/c^{2}]", 100, 0, 100);


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

        vec_unfiltered.SetPtEtaPhiM(0,0,0,0);
        vec_pTCut.SetPtEtaPhiM(0,0,0,0);
        

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
            if ( jetPt < 1.5 /* && jetNConst < 2.5 && pT_LeadConst < 0.5 */ ) continue;
            vec_pTCut = TLorentzVector(jetPx, jetPy, jetPz, jetE);

        } // End of individual jets creation
        
        unfiltered_Jets_InvariantMass->Fill(vec_unfiltered.M());
        cutPt_Jets_InvariantMass->Fill(vec_pTCut.M());

        particles_fastjet.clear();
        jets.clear();
        jets_array->Clear();
        quarks->Clear();
    }

    TCanvas *c1 = new TCanvas("c1", "Histograms and distributions", 2500, 2500);
    c1->Divide(1, 2);

    c1->cd(1);
    unfiltered_Jets_InvariantMass->SetTitle("Unfiltered jets invariant mass [GeV/c^{2}]");
    unfiltered_Jets_InvariantMass->GetXaxis()->SetTitle("Transverse Momentum [GeV/c]");
    unfiltered_Jets_InvariantMass->GetYaxis()->SetTitle("Frequency");
    unfiltered_Jets_InvariantMass->Draw();

    c1->cd(2);
    cutPt_Jets_InvariantMass->SetTitle("Jets invariant mass with cut in p_{T} [GeV/c^{2}]");
    cutPt_Jets_InvariantMass->GetXaxis()->SetTitle("Transverse Momentum [GeV/c]");
    cutPt_Jets_InvariantMass->GetYaxis()->SetTitle("Frequency");
    cutPt_Jets_InvariantMass->Draw();

}