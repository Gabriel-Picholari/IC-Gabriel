
//===========================================================================================================================================================================
// W boson analysis:
// This macro sole purpose is to analyse the first W boson we are forcing to decay into a charm and strange quark
//===========================================================================================================================================================================

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
#include <algorithm>
#include "TPythia8.h"
#include "TRandom3.h"
#include "TParticle.h"
#include <unordered_set>
#include "TClonesArray.h"
#include <TLorentzVector.h>
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

void wBoson_analysis(const char* fileName)
{
    gSystem->Load("libEG");
    gSystem->Load("libEGPythia8");

    //---------------------------------------------------------------------------------------------------------
    // Initialization of histograms
    //---------------------------------------------------------------------------------------------------------

    TH1F *invariantMass = new TH1F("boson_invariantMass_hist", "W^{+-} invariant mass spectrum [GeV/c^{2}]", 100, 0, 120);
    TH1F *boson_energy_distribution = new TH1F("boson_energy_hist", "Initial W^{+-} boson energy distribution", 100, 0, 1000);
    TH1F *boson_eta_distribution = new TH1F("boson_eta_hist", "Initial W^{+-} boson pseudorapidity distribution", 100, -10, 10);
    TH1F *boson_pT_distribution = new TH1F("boson_pT_hist", "Initial W^{+-} boson transverse momentum distribution", 100, 0, 100);

    //---------------------------------------------------------------------------------------------------------
    // Initialization of the .root file and TTrees
    //---------------------------------------------------------------------------------------------------------

    Float_t boson_pT, boson_eta, boson_phi, boson_energy, boson_eventID;

    TFile *file = TFile::Open(fileName, "READ");

    TTree* boson_ttree = (TTree *) file->Get("boson_TTree");
    boson_ttree->SetBranchAddress("boson_pT", &boson_pT);
    boson_ttree->SetBranchAddress("boson_eta", &boson_eta);
    boson_ttree->SetBranchAddress("boson_phi", &boson_phi);
    boson_ttree->SetBranchAddress("boson_energy", &boson_energy);
    boson_ttree->SetBranchAddress("boson_eventID", &boson_eventID);

    //---------------------------------------------------------------------------------------------------------
    // Boson W analysis
    //---------------------------------------------------------------------------------------------------------

    for (Long64_t i = 0; i < boson_ttree->GetEntries(); i++)
    {
        boson_ttree->GetEntry(i);
        TLorentzVector bosonVec;
        bosonVec.SetPtEtaPhiE(boson_pT, boson_eta, boson_phi, boson_energy);
        invariantMass->Fill(bosonVec.M());
        if (abs(boson_eta) > 2) continue;
        boson_energy_distribution->Fill(boson_energy);
        boson_eta_distribution->Fill(boson_eta);
        boson_pT_distribution->Fill(boson_pT);
    }

    //---------------------------------------------------------------------------------------------------------
    // Plotting histograms
    //---------------------------------------------------------------------------------------------------------

    TCanvas *c1 = new TCanvas("c1", "Boson W analysis", 2500, 2500);
    c1->Divide(2, 2);

    c1->cd(1);
    invariantMass->SetTitle("Boson W^{#pm} invariant mass spectrum");
    invariantMass->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    invariantMass->GetYaxis()->SetTitle("Frequency");
    invariantMass->Draw();

    c1->cd(2);
    boson_energy_distribution->SetTitle("Initial W^{#pm} boson energy distribution");
    boson_energy_distribution->GetXaxis()->SetTitle("Energy [GeV]");
    boson_energy_distribution->GetYaxis()->SetTitle("Frequency");
    boson_energy_distribution->Draw(); 

    c1->cd(3);
    boson_eta_distribution->SetTitle("Initial W^{#pm} boson pseudorapidity distribution");
    boson_eta_distribution->GetXaxis()->SetTitle("Pseudorapidity");
    boson_eta_distribution->GetYaxis()->SetTitle("Frequency");
    boson_eta_distribution->Draw(); 

    c1->cd(4);
    boson_pT_distribution->SetTitle("Initial W^{#pm} boson transverse momentum distribution");
    boson_pT_distribution->GetXaxis()->SetTitle("Transverse momentum [GeV/c]");
    boson_pT_distribution->GetYaxis()->SetTitle("Frequency");
    boson_pT_distribution->Draw();

    c1->SaveAs("wBoson_analysis_results.png");
}
