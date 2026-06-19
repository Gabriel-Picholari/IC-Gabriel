//===========================================================================================================================================================================
// Cuts in jets distributions analysis:
// This macro has the purpose of analysing the distributions of the discriminatory variables chosen to train the classifiers. The cuts determined here through analysis of
// the signal and background distributions will be applied in a classical analysis of the W boson decay. Thus, we will be able to compare the classical and the machine
// learning approaches.
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

void cuts_analysis(const char* charm_file, const char* strange_file)
{
    gSystem->Load("libEG");
    gSystem->Load("libEGPythia8");

    //---------------------------------------------------------------------------------------------------------
    // Initialization of the .root file and TTrees
    //---------------------------------------------------------------------------------------------------------

    Float_t pT_s, nConst_s, nRho_s, eventID_s;
    Float_t pT_c, nConst_c, nRho_c, eventID_c;
    std::vector<Float_t> *jetVerticesInvariantMasses_s = nullptr; 
    std::vector<Float_t> *jetVerticesInvariantMasses_c = nullptr;

    TFile *s_file = TFile::Open(strange_file, "READ");

    TTree *signalTree_s = dynamic_cast<TTree *>(s_file->Get("SignalTree_s"));
    signalTree_s->SetBranchAddress("pT_s", &pT_s);
    signalTree_s->SetBranchAddress("nConst_s", &nConst_s);
    signalTree_s->SetBranchAddress("nRho_s", &nRho_s);
    signalTree_s->SetBranchAddress("eventID_s", &eventID_s);
    signalTree_s->SetBranchAddress("jetVerticesInvariantMasses_s", &jetVerticesInvariantMasses_s);

    TTree *backgroundTree_s = dynamic_cast<TTree *>(s_file->Get("BackgroundTree_s")); 
    backgroundTree_s->SetBranchAddress("pT_s", &pT_s);
    backgroundTree_s->SetBranchAddress("nConst_s", &nConst_s);
    backgroundTree_s->SetBranchAddress("nRho_s", &nRho_s);
    backgroundTree_s->SetBranchAddress("eventID_s", &eventID_s);
    backgroundTree_s->SetBranchAddress("jetVerticesInvariantMasses_s", &jetVerticesInvariantMasses_s);

    Long64_t ne_signal_s = signalTree_s->GetEntries();
    Long64_t ne_background_s = backgroundTree_s->GetEntries();

    //---------------------------------------------------------------------------------------------------------

    TFile *c_file = TFile::Open(charm_file, "READ");

    TTree *signalTree_c = dynamic_cast<TTree *>(c_file->Get("SignalTree_c"));
    signalTree_c->SetBranchAddress("pT_c", &pT_c);
    signalTree_c->SetBranchAddress("nConst_c", &nConst_c);
    signalTree_c->SetBranchAddress("nRho_c", &nRho_c);
    signalTree_c->SetBranchAddress("eventID_c", &eventID_c);
    signalTree_c->SetBranchAddress("jetVerticesInvariantMasses_c", &jetVerticesInvariantMasses_c);

    TTree *backgroundTree_c = dynamic_cast<TTree *>(c_file->Get("BackgroundTree_c"));
    backgroundTree_c->SetBranchAddress("pT_c", &pT_c);
    backgroundTree_c->SetBranchAddress("nConst_c", &nConst_c);
    backgroundTree_c->SetBranchAddress("nRho_c", &nRho_c);
    backgroundTree_c->SetBranchAddress("eventID_c", &eventID_c);
    backgroundTree_c->SetBranchAddress("jetVerticesInvariantMasses_c", &jetVerticesInvariantMasses_c);

    Long64_t ne_signal_c = signalTree_c->GetEntries();
    Long64_t ne_background_c = backgroundTree_c->GetEntries();

    //---------------------------------------------------------------------------------------------------------
    // Analysis
    //---------------------------------------------------------------------------------------------------------

    // From my perspective the optimal approach to dealing with this pre analysis is to go through each of the trees separately, build the histograms and then plot
    // them in conjoint. Thats what will be done from now on!


    Float_t pTCut = 20;
    Int_t nConstCut = 15;
    Int_t nRhoCut = 2;

    TH1F *strange_signal_pT = new TH1F("strange_signal_pT", "p_{T} distribution of strange jets; p_{T} [GeV/c]; Events", 100, 0, 100);
    TH1F *strange_signal_nConst = new TH1F("strange_signal_nConst", "Number of constituents distribution of strange jets; nConst; Events", 50, 0, 50);
    TH1F *strange_signal_nRho = new TH1F("strange_signal_nRho", "N_{#rho} distribution of strange jets; N_{#rho}; Events", 10, 0, 10);
    TH1F *strange_signal_jetVerticesInvariantMasses = new TH1F("strange_signal_jetVerticesInvariantMasses", "Jet vertices invariant masses distribution of strange jets; Jet vertices invariant mass (GeV/c^{2}); Events", 100, 0, 5);

    TH1F *strange_signal_pT_pTCut = new TH1F("strange_signal_pT_pTCut", ("p_{T} distribution of strange jets (p_{T} > " + std::to_string(pTCut) + " GeV/c); p_{T} [GeV/c]; Events").c_str(), 100, 0, 100);
    TH1F *strange_signal_nConst_pTCut = new TH1F("strange_signal_nConst_pTCut", ("Number of constituents distribution of strange jets (p_{T} > " + std::to_string(pTCut) + " GeV/c); nConst; Events").c_str(), 50, 0, 50);
    TH1F *strange_signal_nRho_pTCut = new TH1F("strange_signal_nRho_pTCut", ("N_{#rho} distribution of strange jets (p_{T} > " + std::to_string(pTCut) + " GeV/c); N_{#rho}; Events").c_str(), 10, 0, 10);
    TH1F *strange_signal_jetVerticesInvariantMasses_pTCut = new TH1F("strange_signal_jetVerticesInvariantMasses_pTCut", ("Jet vertices invariant masses distribution of strange jets (p_{T} > " + std::to_string(pTCut) + " GeV/c); Jet vertices invariant mass (GeV/c^{2}); Events").c_str(), 100, 0, 5);

    TH1F *strange_signal_pT_nConstCut = new TH1F("strange_signal_pT_nConstCut", ("p_{T} distribution of strange jets (N_{const} > " + std::to_string(nConstCut) + "); p_{T} [GeV/c]; Events").c_str(), 100, 0, 100);
    TH1F *strange_signal_nConst_nConstCut = new TH1F("strange_signal_nConst_nConstCut", ("Number of constituents distribution of strange jets (N_{const} > " + std::to_string(nConstCut) + "); nConst; Events").c_str(), 50, 0, 50);
    TH1F *strange_signal_nRho_nConstCut = new TH1F("strange_signal_nRho_nConstCut", ("N_{#rho} distribution of strange jets (N_{const} > " + std::to_string(nConstCut) + "); N_{#rho}; Events").c_str(), 10, 0, 10);
    TH1F *strange_signal_jetVerticesInvariantMasses_nConstCut = new TH1F("strange_signal_jetVerticesInvariantMasses_nConstCut", ("Jet vertices invariant masses distribution of strange jets (N_{const} > " + std::to_string(nConstCut) + "); Jet vertices invariant mass (GeV/c^{2}); Events").c_str(), 100, 0, 5);

    TH1F *strange_signal_pT_nRhoCut = new TH1F("strange_signal_pT_nRhoCut", ("p_{T} distribution of strange jets (N_{#rho} > " + std::to_string(nRhoCut) + "); p_{T} [GeV/c]; Events").c_str(), 100, 0, 100);
    TH1F *strange_signal_nConst_nRhoCut = new TH1F("strange_signal_nConst_nRhoCut", ("Number of constituents distribution of strange jets (N_{#rho} > " + std::to_string(nRhoCut) + "); nConst; Events").c_str(), 50, 0, 50);
    TH1F *strange_signal_nRho_nRhoCut = new TH1F("strange_signal_nRho_nRhoCut", ("N_{#rho} distribution of strange jets (N_{#rho} > " + std::to_string(nRhoCut) + "); N_{#rho}; Events").c_str(), 10, 0, 10);
    TH1F *strange_signal_jetVerticesInvariantMasses_nRhoCut = new TH1F("strange_signal_jetVerticesInvariantMasses_nRhoCut", ("Jet vertices invariant masses distribution of strange jets (N_{#rho} > " + std::to_string(nRhoCut) + "); Jet vertices invariant mass (GeV/c^{2}); Events").c_str(), 100, 0, 5);

    TH1F *strange_background_pT = new TH1F("strange_background_pT", "p_{T} distribution of non-strange jets; p_{T} [GeV/c]; Events", 100, 0, 100);
    TH1F *strange_background_nConst = new TH1F("strange_background_nConst", "Number of constituents distribution of non-strange jets; nConst; Events", 50, 0, 50);
    TH1F *strange_background_nRho = new TH1F("strange_background_nRho", "N_{#rho} distribution of non-strange jets; N_{#rho}; Events", 10, 0, 10);
    TH1F *strange_background_jetVerticesInvariantMasses = new TH1F("strange_background_jetVerticesInvariantMasses", "Jet vertices invariant masses distribution of non-strange jets; Jet vertices invariant mass (GeV/c^{2}); Events", 100, 0, 5);

    TH1F *strange_background_pT_pTCut = new TH1F("strange_background_pT_pTCut", ("p_{T} distribution of non-strange jets (p_{T} > " + std::to_string(pTCut) + " GeV/c); p_{T} [GeV/c]; Events").c_str(), 100, 0, 100);
    TH1F *strange_background_nConst_pTCut = new TH1F("strange_background_nConst_pTCut", ("Number of constituents distribution of non-strange jets (p_{T} > " + std::to_string(pTCut) + " GeV/c); nConst; Events").c_str(), 50, 0, 50);
    TH1F *strange_background_nRho_pTCut = new TH1F("strange_background_nRho_pTCut", ("N_{#rho} distribution of non-strange jets (p_{T} > " + std::to_string(pTCut) + " GeV/c); N_{#rho}; Events").c_str(), 10, 0, 10);
    TH1F *strange_background_jetVerticesInvariantMasses_pTCut = new TH1F("strange_background_jetVerticesInvariantMasses_pTCut", ("Jet vertices invariant masses distribution of non-strange jets (p_{T} > " + std::to_string(pTCut) + " GeV/c); Jet vertices invariant mass (GeV/c^{2}); Events").c_str(), 100, 0, 5);

    TH1F *strange_background_pT_nConstCut = new TH1F("strange_background_pT_nConstCut", ("p_{T} distribution of non-strange jets (N_{const} > " + std::to_string(nConstCut) + "); p_{T} [GeV/c]; Events").c_str(), 100, 0, 100);
    TH1F *strange_background_nConst_nConstCut = new TH1F("strange_background_nConst_nConstCut", ("Number of constituents distribution of non-strange jets (N_{const} > " + std::to_string(nConstCut) + "); nConst; Events").c_str(), 50, 0, 50);
    TH1F *strange_background_nRho_nConstCut = new TH1F("strange_background_nRho_nConstCut", ("N_{#rho} distribution of non-strange jets (N_{const} > " + std::to_string(nConstCut) + "); N_{#rho}; Events").c_str(), 10, 0, 10);
    TH1F *strange_background_jetVerticesInvariantMasses_nConstCut = new TH1F("strange_background_jetVerticesInvariantMasses_nConstCut", ("Jet vertices invariant masses distribution of non-strange jets (N_{const} > " + std::to_string(nConstCut) + "); Jet vertices invariant mass (GeV/c^{2}); Events").c_str(), 100, 0, 5);

    TH1F *strange_background_pT_nRhoCut = new TH1F("strange_background_pT_nRhoCut", ("p_{T} distribution of non-strange jets (N_{#rho} > " + std::to_string(nRhoCut) + "); p_{T} [GeV/c]; Events").c_str(), 100, 0, 100);
    TH1F *strange_background_nConst_nRhoCut = new TH1F("strange_background_nConst_nRhoCut", ("Number of constituents distribution of non-strange jets (N_{#rho} > " + std::to_string(nRhoCut) + "); nConst; Events").c_str(), 50, 0, 50);
    TH1F *strange_background_nRho_nRhoCut = new TH1F("strange_background_nRho_nRhoCut", ("N_{#rho} distribution of non-strange jets (N_{#rho} > " + std::to_string(nRhoCut) + "); N_{#rho}; Events").c_str(), 10, 0, 10);
    TH1F *strange_background_jetVerticesInvariantMasses_nRhoCut = new TH1F("strange_background_jetVerticesInvariantMasses_nRhoCut", ("Jet vertices invariant masses distribution of non-strange jets (N_{#rho} > " + std::to_string(nRhoCut) + "); Jet vertices invariant mass (GeV/c^{2}); Events").c_str(), 100, 0, 5);
    
    for (Long64_t i = 0; i < ne_signal_s; i++)
    {
        signalTree_s->GetEntry(i);

        strange_signal_pT->Fill(pT_s);
        strange_signal_nConst->Fill(nConst_s);
        strange_signal_nRho->Fill(nRho_s);

        for (const Float_t &mass : *jetVerticesInvariantMasses_s)
        {
            strange_signal_jetVerticesInvariantMasses->Fill(mass);
        }

        if (pT_s > pTCut)
        {
            strange_signal_pT_pTCut->Fill(pT_s);
            strange_signal_nConst_pTCut->Fill(nConst_s);
            strange_signal_nRho_pTCut->Fill(nRho_s);

            for (const Float_t &mass : *jetVerticesInvariantMasses_s)
            {
                strange_signal_jetVerticesInvariantMasses_pTCut->Fill(mass);
            }
        }

        if (nConst_s > nConstCut)
        {
            strange_signal_pT_nConstCut->Fill(pT_s);
            strange_signal_nConst_nConstCut->Fill(nConst_s);
            strange_signal_nRho_nConstCut->Fill(nRho_s);

            for (const Float_t &mass : *jetVerticesInvariantMasses_s)
            {
                strange_signal_jetVerticesInvariantMasses_nConstCut->Fill(mass);
            }
        }

        if (nRho_s > nRhoCut)
        {
            strange_signal_pT_nRhoCut->Fill(pT_s);
            strange_signal_nConst_nRhoCut->Fill(nConst_s);
            strange_signal_nRho_nRhoCut->Fill(nRho_s);

            for (const Float_t &mass : *jetVerticesInvariantMasses_s)
            {
                strange_signal_jetVerticesInvariantMasses_nRhoCut->Fill(mass);
            }
        }
    }

    for (Long64_t i = 0; i < ne_background_s; i++)
    {
        backgroundTree_s->GetEntry(i);

        strange_background_pT->Fill(pT_s);
        strange_background_nConst->Fill(nConst_s);
        strange_background_nRho->Fill(nRho_s);

        for (const Float_t &mass : *jetVerticesInvariantMasses_s)
        {
            strange_background_jetVerticesInvariantMasses->Fill(mass);
        }

        if (pT_s > pTCut)
        {
            strange_background_pT_pTCut->Fill(pT_s);
            strange_background_nConst_pTCut->Fill(nConst_s);
            strange_background_nRho_pTCut->Fill(nRho_s);
        }

        if (nConst_s > nConstCut)
        {
            strange_background_pT_nConstCut->Fill(pT_s);
            strange_background_nConst_nConstCut->Fill(nConst_s);
            strange_background_nRho_nConstCut->Fill(nRho_s);
        }

        if (nRho_s > nRhoCut)
        {
            strange_background_pT_nRhoCut->Fill(pT_s);
            strange_background_nConst_nRhoCut->Fill(nConst_s);
            strange_background_nRho_nRhoCut->Fill(nRho_s);
        }
    }

    //---------------------------------------------------------------------------------------------------------

    TH1F *charm_signal_pT = new TH1F("charm_signal_pT", "p_{T} distribution of charm jets; p_{T} [GeV/c]; Events", 100, 0, 100);
    TH1F *charm_signal_nConst = new TH1F("charm_signal_nConst", "Number of constituents distribution of charm jets; nConst; Events", 50, 0, 50);
    TH1F *charm_signal_nRho = new TH1F("charm_signal_nRho", "N_{#rho} distribution of charm jets; N_{#rho}; Events", 10, 0, 10);
    TH1F *charm_signal_jetVerticesInvariantMasses = new TH1F("charm_signal_jetVerticesInvariantMasses", "Jet vertices invariant masses distribution of charm jets; Jet vertices invariant mass (GeV/c^{2}); Events", 100, 0, 5);

    TH1F *charm_signal_pT_pTCut = new TH1F("charm_signal_pT_pTCut", ("p_{T} distribution of charm jets (p_{T} > " + std::to_string(pTCut) + " GeV/c); p_{T} [GeV/c]; Events").c_str(), 100, 0, 100);
    TH1F *charm_signal_nConst_pTCut = new TH1F("charm_signal_nConst_pTCut", ("Number of constituents distribution of charm jets (p_{T} > " + std::to_string(pTCut) + " GeV/c); nConst; Events").c_str(), 50, 0, 50);
    TH1F *charm_signal_nRho_pTCut = new TH1F("charm_signal_nRho_pTCut", ("N_{#rho} distribution of charm jets (p_{T} > " + std::to_string(pTCut) + " GeV/c); N_{#rho}; Events").c_str(), 10, 0, 10);
    TH1F *charm_signal_jetVerticesInvariantMasses_pTCut = new TH1F("charm_signal_jetVerticesInvariantMasses_pTCut", ("Jet vertices invariant masses distribution of charm jets (p_{T} > " + std::to_string(pTCut) + " GeV/c); Jet vertices invariant mass (GeV/c^{2}); Events").c_str(), 100, 0, 5);

    TH1F *charm_signal_pT_nConstCut = new TH1F("charm_signal_pT_nConstCut", ("p_{T} distribution of charm jets (N_{const} > " + std::to_string(nConstCut) + "); p_{T} [GeV/c]; Events").c_str(), 100, 0, 100);
    TH1F *charm_signal_nConst_nConstCut = new TH1F("charm_signal_nConst_nConstCut", ("Number of constituents distribution of charm jets (N_{const} > " + std::to_string(nConstCut) + "); nConst; Events").c_str(), 50, 0, 50);
    TH1F *charm_signal_nRho_nConstCut = new TH1F("charm_signal_nRho_nConstCut", ("N_{#rho} distribution of charm jets (N_{const} > " + std::to_string(nConstCut) + "); N_{#rho}; Events").c_str(), 10, 0, 10);
    TH1F *charm_signal_jetVerticesInvariantMasses_nConstCut = new TH1F("charm_signal_jetVerticesInvariantMasses_nConstCut", ("Jet vertices invariant masses distribution of charm jets (N_{const} > " + std::to_string(nConstCut) + "); Jet vertices invariant mass (GeV/c^{2}); Events").c_str(), 100, 0, 5);

    TH1F *charm_signal_pT_nRhoCut = new TH1F("charm_signal_pT_nRhoCut", ("p_{T} distribution of charm jets (N_{#rho} > " + std::to_string(nRhoCut) + "); p_{T} [GeV/c]; Events").c_str(), 100, 0, 100);
    TH1F *charm_signal_nConst_nRhoCut = new TH1F("charm_signal_nConst_nRhoCut", ("Number of constituents distribution of charm jets (N_{#rho} > " + std::to_string(nRhoCut) + "); nConst; Events").c_str(), 50, 0, 50);
    TH1F *charm_signal_nRho_nRhoCut = new TH1F("charm_signal_nRho_nRhoCut", ("N_{#rho} distribution of charm jets (N_{#rho} > " + std::to_string(nRhoCut) + "); N_{#rho}; Events").c_str(), 10, 0, 10);
    TH1F *charm_signal_jetVerticesInvariantMasses_nRhoCut = new TH1F("charm_signal_jetVerticesInvariantMasses_nRhoCut", ("Jet vertices invariant masses distribution of charm jets (N_{#rho} > " + std::to_string(nRhoCut) + "); Jet vertices invariant mass (GeV/c^{2}); Events").c_str(), 100, 0, 5);

    TH1F *charm_background_pT = new TH1F("charm_background_pT", "p_{T} distribution of non-charm jets; p_{T} [GeV/c]; Events", 100, 0, 100);
    TH1F *charm_background_nConst = new TH1F("charm_background_nConst", "Number of constituents distribution of non-charm jets; nConst; Events", 50, 0, 50);
    TH1F *charm_background_nRho = new TH1F("charm_background_nRho", "N_{#rho} distribution of non-charm jets; N_{#rho}; Events", 10, 0, 10);
    TH1F *charm_background_jetVerticesInvariantMasses = new TH1F("charm_background_jetVerticesInvariantMasses", "Jet vertices invariant masses distribution of charm jets; Jet vertices invariant mass (GeV/c^{2}); Events", 100, 0, 5);

    TH1F *charm_background_pT_pTCut = new TH1F("charm_background_pT_pTCut", ("p_{T} distribution of non-charm jets (p_{T} > " + std::to_string(pTCut) + " GeV/c); p_{T} [GeV/c]; Events").c_str(), 100, 0, 100);
    TH1F *charm_background_nConst_pTCut = new TH1F("charm_background_nConst_pTCut", ("Number of constituents distribution of non-charm jets (p_{T} > " + std::to_string(pTCut) + " GeV/c); nConst; Events").c_str(), 50, 0, 50);
    TH1F *charm_background_nRho_pTCut = new TH1F("charm_background_nRho_pTCut", ("N_{#rho} distribution of non-charm jets (p_{T} > " + std::to_string(pTCut) + " GeV/c); N_{#rho}; Events").c_str(), 10, 0, 10);
    TH1F *charm_background_jetVerticesInvariantMasses_pTCut = new TH1F("charm_background_jetVerticesInvariantMasses_pTCut", ("Jet vertices invariant masses distribution of non-charm jets (p_{T} > " + std::to_string(pTCut) + " GeV/c); Jet vertices invariant mass (GeV/c^{2}); Events").c_str(), 100, 0, 5);

    TH1F *charm_background_pT_nConstCut = new TH1F("charm_background_pT_nConstCut", ("p_{T} distribution of non-charm jets (N_{const} > " + std::to_string(nConstCut) + "); p_{T} [GeV/c]; Events").c_str(), 100, 0, 100);
    TH1F *charm_background_nConst_nConstCut = new TH1F("charm_background_nConst_nConstCut", ("Number of constituents distribution of non-charm jets (N_{const} > " + std::to_string(nConstCut) + "); nConst; Events").c_str(), 50, 0, 50);
    TH1F *charm_background_nRho_nConstCut = new TH1F("charm_background_nRho_nConstCut", ("N_{#rho} distribution of non-charm jets (N_{const} > " + std::to_string(nConstCut) + "); N_{#rho}; Events").c_str(), 10, 0, 10);
    TH1F *charm_background_jetVerticesInvariantMasses_nConstCut = new TH1F("charm_background_jetVerticesInvariantMasses_nConstCut", ("Jet vertices invariant masses distribution of non-charm jets (N_{const} > " + std::to_string(nConstCut) + "); Jet vertices invariant mass (GeV/c^{2}); Events").c_str(), 100, 0, 5);

    TH1F *charm_background_pT_nRhoCut = new TH1F("charm_background_pT_nRhoCut", ("p_{T} distribution of non-charm jets (N_{#rho} > " + std::to_string(nRhoCut) + "); p_{T} [GeV/c]; Events").c_str(), 100, 0, 100);
    TH1F *charm_background_nConst_nRhoCut = new TH1F("charm_background_nConst_nRhoCut", ("Number of constituents distribution of non-charm jets (N_{#rho} > " + std::to_string(nRhoCut) + "); nConst; Events").c_str(), 50, 0, 50);
    TH1F *charm_background_nRho_nRhoCut = new TH1F("charm_background_nRho_nRhoCut", ("N_{#rho} distribution of non-charm jets (N_{#rho} > " + std::to_string(nRhoCut) + "); N_{#rho}; Events").c_str(), 10, 0, 10);
    TH1F *charm_background_jetVerticesInvariantMasses_nRhoCut = new TH1F("charm_background_jetVerticesInvariantMasses_nRhoCut", ("Jet vertices invariant masses distribution of non-charm jets (N_{#rho} > " + std::to_string(nRhoCut) + "); Jet vertices invariant mass (GeV/c^{2}); Events").c_str(), 100, 0, 5);
    
    for (Long64_t i = 0; i < ne_signal_c; i++)
    {
        signalTree_c->GetEntry(i);
        charm_signal_pT->Fill(pT_c);
        charm_signal_nConst->Fill(nConst_c);
        charm_signal_nRho->Fill(nRho_c);

        for (const Float_t &mass : *jetVerticesInvariantMasses_c)
        {
            charm_signal_jetVerticesInvariantMasses->Fill(mass);
        }

        if (pT_c > pTCut)
        {
            charm_signal_pT_pTCut->Fill(pT_c);
            charm_signal_nConst_pTCut->Fill(nConst_c);
            charm_signal_nRho_pTCut->Fill(nRho_c); 
        }

        if (nConst_c > nConstCut)
        {
            charm_signal_pT_nConstCut->Fill(pT_c);
            charm_signal_nConst_nConstCut->Fill(nConst_c);
            charm_signal_nRho_nConstCut->Fill(nRho_c);
        }

        if (nRho_c > nRhoCut)
        {
            charm_signal_pT_nRhoCut->Fill(pT_c);
            charm_signal_nConst_nRhoCut->Fill(nConst_c);
            charm_signal_nRho_nRhoCut->Fill(nRho_c);  
        }
    }

    for (Long64_t i = 0; i < ne_background_c; i++)
    {
        backgroundTree_c->GetEntry(i);

        charm_background_pT->Fill(pT_c);
        charm_background_nConst->Fill(nConst_c);
        charm_background_nRho->Fill(nRho_c);

        for (const Float_t &mass : *jetVerticesInvariantMasses_c)
        {
            charm_background_jetVerticesInvariantMasses->Fill(mass);
        }

        if (pT_c > pTCut)
        {
            charm_background_pT_pTCut->Fill(pT_c);
            charm_background_nConst_pTCut->Fill(nConst_c);
            charm_background_nRho_pTCut->Fill(nRho_c);
        }

        if (nConst_c > nConstCut)
        {
            charm_background_pT_nConstCut->Fill(pT_c);
            charm_background_nConst_nConstCut->Fill(nConst_c);
            charm_background_nRho_nConstCut->Fill(nRho_c);
        }

        if (nRho_c > nRhoCut)
        {
            charm_background_pT_nRhoCut->Fill(pT_c);
            charm_background_nConst_nRhoCut->Fill(nConst_c);
            charm_background_nRho_nRhoCut->Fill(nRho_c);
        }
    }

    // The background distributions are much more populated than the signal ones. Then, to make the signal and background distributions comparable, we will rescale
    // the background histograms accordingly.

    Float_t scaleS = (Float_t)(ne_signal_s) / (Float_t)(ne_background_s);
    
    strange_background_pT->Scale(scaleS);

    strange_background_nConst->Scale(scaleS);
    strange_background_nRho->Scale(scaleS);
    strange_background_jetVerticesInvariantMasses->Scale(scaleS);

    strange_background_pT_pTCut->Scale(scaleS);
    strange_background_nConst_pTCut->Scale(scaleS);
    strange_background_nRho_pTCut->Scale(scaleS);

    strange_background_pT_nConstCut->Scale(scaleS);
    strange_background_nConst_nConstCut->Scale(scaleS);
    strange_background_nRho_nConstCut->Scale(scaleS);

    strange_background_pT_nRhoCut->Scale(scaleS);
    strange_background_nConst_nRhoCut->Scale(scaleS);
    strange_background_nRho_nRhoCut->Scale(scaleS);

    Float_t scaleC = (Float_t)(ne_signal_c) / (Float_t)(ne_background_c);
    charm_background_pT->Scale(scaleC);
    charm_background_nConst->Scale(scaleC);
    charm_background_nRho->Scale(scaleC);
    charm_background_jetVerticesInvariantMasses->Scale(scaleC);

    charm_background_pT_pTCut->Scale(scaleC);
    charm_background_nConst_pTCut->Scale(scaleC);
    charm_background_nRho_pTCut->Scale(scaleC);

    charm_background_pT_nConstCut->Scale(scaleC);
    charm_background_nConst_nConstCut->Scale(scaleC);
    charm_background_nRho_nConstCut->Scale(scaleC);

    charm_background_pT_nRhoCut->Scale(scaleC);
    charm_background_nConst_nRhoCut->Scale(scaleC);
    charm_background_nRho_nRhoCut->Scale(scaleC);

    //=====================================================================
    // Strange jets canvases
    //=====================================================================

    TCanvas *c1 = new TCanvas("c1", "Transverse momentum strange variables distributions", 2500, 2500);
    c1->Divide(2, 2);

    c1->cd(1);
    strange_background_pT->SetTitle("Strange jets: p_{T} distribution;p_{T} [GeV/c];Events");
    strange_background_pT->GetYaxis()->SetRangeUser(std::min(strange_background_pT->GetMinimum(), strange_signal_pT->GetMinimum()) * 1.5, std::max(strange_background_pT->GetMaximum(), strange_signal_pT->GetMaximum()) * 1.5);
    strange_background_pT->SetLineColor(kRed);
    strange_background_pT->SetLineStyle(1);
    strange_background_pT->DrawCopy(); // Desenha o fundo primeiro

    strange_signal_pT->SetLineColor(kBlue);
    strange_signal_pT->SetLineStyle(1);
    strange_signal_pT->DrawCopy("same");

    TLegend *legend1 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend1->AddEntry(strange_signal_pT, "Strange signal jets", "l");
    legend1->AddEntry(strange_background_pT, "Strange background jets", "l");
    legend1->Draw();

    std::ostringstream oss;
    oss << std::fixed << std::setprecision(1) << pTCut;

    std::string title = "Strange jets: p_{T} distribution (p_{T} > " + oss.str() + " GeV/c);p_{T} [GeV/c];Events";

    c1->cd(2);
    strange_background_pT_pTCut->SetTitle(title.c_str());
    strange_background_pT_pTCut->GetYaxis()->SetRangeUser(std::min(strange_background_pT_pTCut->GetMinimum(), strange_signal_pT_pTCut->GetMinimum()) * 1.5, std::max(strange_background_pT_pTCut->GetMaximum(), strange_signal_pT_pTCut->GetMaximum()) * 1.5);
    strange_background_pT_pTCut->SetLineColor(kRed);
    strange_background_pT_pTCut->SetLineStyle(1);
    strange_background_pT_pTCut->DrawCopy();

    strange_signal_pT_pTCut->SetLineColor(kBlue);
    strange_signal_pT_pTCut->SetLineStyle(1);
    strange_signal_pT_pTCut->DrawCopy("same");

    TLegend *legend2 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend1->AddEntry(strange_signal_pT_pTCut, "Strange signal jets", "l");
    legend1->AddEntry(strange_background_pT_pTCut, "Strange background jets", "l");
    legend1->Draw();

    c1->cd(3);
    strange_background_pT_nConstCut->SetTitle(("Strange jets: p_{T} distribution (N_{const} > " + std::to_string(nConstCut) + ");p_{T} [GeV/c];Events").c_str());
    strange_background_pT_nConstCut->GetYaxis()->SetRangeUser(std::min(strange_background_pT_nConstCut->GetMinimum(), strange_signal_pT_nConstCut->GetMinimum()) * 1.5, std::max(strange_background_pT_nConstCut->GetMaximum(), strange_signal_pT_nConstCut->GetMaximum()) * 1.5);
    strange_background_pT_nConstCut->SetLineColor(kRed);
    strange_background_pT_nConstCut->SetLineStyle(1);
    strange_background_pT_nConstCut->DrawCopy();
    
    strange_signal_pT_nConstCut->SetLineColor(kBlue);
    strange_signal_pT_nConstCut->SetLineStyle(1);
    strange_signal_pT_nConstCut->DrawCopy("same");

    TLegend *legend3 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend3->AddEntry(strange_signal_pT_nConstCut, "Strange signal jets", "l");
    legend3->AddEntry(strange_background_pT_nConstCut, "Strange background jets", "l");
    legend3->Draw();

    c1->cd(4);
    strange_background_pT_nRhoCut->SetTitle(("Strange jets: p_{T} distribution (N_{#rho} > " + std::to_string(nRhoCut) + ");p_{T} [GeV/c];Events").c_str());
    strange_background_pT_nRhoCut->GetYaxis()->SetRangeUser(std::min(strange_background_pT_nRhoCut->GetMinimum(), strange_signal_pT_nRhoCut->GetMinimum()) * 1.5, std::max(strange_background_pT_nRhoCut->GetMaximum(), strange_signal_pT_nRhoCut->GetMaximum()) * 1.5);
    strange_background_pT_nRhoCut->SetLineColor(kRed);
    strange_background_pT_nRhoCut->SetLineStyle(1);
    strange_background_pT_nRhoCut->DrawCopy();

    strange_signal_pT_nRhoCut->SetLineColor(kBlue);
    strange_signal_pT_nRhoCut->SetLineStyle(1);
    strange_signal_pT_nRhoCut->DrawCopy("same");

    TLegend *legend4 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend4->AddEntry(strange_signal_pT_nRhoCut, "Strange signal jets", "l");
    legend4->AddEntry(strange_background_pT_nRhoCut, "Strange background jets", "l");
    legend4->Draw();

    TCanvas *c2 = new TCanvas("c2", "Number of constituents strange variables distributions", 2500, 2500);
    c2->Divide(2,2);

    c2->cd(1);
    strange_background_nConst->SetTitle("Strange jets: number of constituents distribution;N_{const};Events");
    strange_background_nConst->GetYaxis()->SetRangeUser(std::min(strange_background_nConst->GetMinimum(), strange_signal_nConst->GetMinimum()) * 1.5, std::max(strange_background_nConst->GetMaximum(), strange_signal_nConst->GetMaximum()) * 1.5);
    strange_background_nConst->SetLineColor(kRed);
    strange_background_nConst->SetLineStyle(1);
    strange_background_nConst->DrawCopy();

    strange_signal_nConst->SetLineColor(kBlue);
    strange_signal_nConst->SetLineStyle(1);
    strange_signal_nConst->DrawCopy("same");

    TLegend *legend5 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend2->AddEntry(strange_signal_nConst, "Strange signal jets", "l");
    legend2->AddEntry(strange_background_nConst, "Strange background jets", "l");
    legend2->Draw();

    std::string title2 ="Strange jets: number of constituents distribution (p_{T} > " + oss.str() + " GeV/c);N_{const};Events";

    c2->cd(2);
    strange_background_nConst_pTCut->SetTitle(title2.c_str());
    strange_background_nConst_pTCut->GetYaxis()->SetRangeUser(std::min(strange_background_nConst_pTCut->GetMinimum(), strange_signal_nConst_pTCut->GetMinimum()) * 1.5, std::max(strange_background_nConst_pTCut->GetMaximum(), strange_signal_nConst_pTCut->GetMaximum()) * 1.5);
    strange_background_nConst_pTCut->SetLineColor(kRed);
    strange_background_nConst_pTCut->SetLineStyle(1);
    strange_background_nConst_pTCut->DrawCopy();

    strange_signal_nConst_pTCut->SetLineColor(kBlue);
    strange_signal_nConst_pTCut->SetLineStyle(1);
    strange_signal_nConst_pTCut->DrawCopy("same");

    TLegend *legend6 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend6->AddEntry(strange_signal_nConst_pTCut, "Strange signal jets", "l");
    legend6->AddEntry(strange_background_nConst_pTCut, "Strange background jets", "l");
    legend6->Draw();

    c2->cd(3);
    strange_background_nConst_nConstCut->SetTitle(("Strange jets: number of constituents distribution (N_{const} > " + std::to_string((Int_t)nConstCut) + ");N_{const};Events").c_str());
    strange_background_nConst_nConstCut->GetYaxis()->SetRangeUser(std::min(strange_background_nConst_nConstCut->GetMinimum(), strange_signal_nConst_nConstCut->GetMinimum()) * 1.5, std::max(strange_background_nConst_nConstCut->GetMaximum(), strange_signal_nConst_nConstCut->GetMaximum()) * 1.5);
    strange_background_nConst_nConstCut->SetLineColor(kRed);
    strange_background_nConst_nConstCut->SetLineStyle(1);
    strange_background_nConst_nConstCut->DrawCopy();

    strange_signal_nConst_nConstCut->SetLineColor(kBlue);
    strange_signal_nConst_nConstCut->SetLineStyle(1);
    strange_signal_nConst_nConstCut->DrawCopy("same");

    TLegend *legend7 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend7->AddEntry(strange_signal_nConst_nConstCut, "Strange signal jets", "l");
    legend7->AddEntry(strange_background_nConst_nConstCut, "Strange background jets", "l");
    legend7->Draw();

    c2->cd(4);
    strange_background_nConst_nRhoCut->SetTitle(("Strange jets: number of constituents distribution (N_{#rho} > " + std::to_string(nRhoCut) + ");N_{const};Events").c_str());
    strange_background_nConst_nRhoCut->GetYaxis()->SetRangeUser(std::min(strange_background_nConst_nRhoCut->GetMinimum(), strange_signal_nConst_nRhoCut->GetMinimum()) * 1.5, std::max(strange_background_nConst_nRhoCut->GetMaximum(), strange_signal_nConst_nRhoCut->GetMaximum()) * 1.5);
    strange_background_nConst_nRhoCut->SetLineColor(kRed);
    strange_background_nConst_nRhoCut->SetLineStyle(1);
    strange_background_nConst_nRhoCut->DrawCopy();

    strange_signal_nConst_nRhoCut->SetLineColor(kBlue);
    strange_signal_nConst_nRhoCut->SetLineStyle(1);
    strange_signal_nConst_nRhoCut->DrawCopy("same");

    TLegend *legend8 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend8->AddEntry(strange_signal_nConst_nRhoCut, "Strange signal jets", "l");
    legend8->AddEntry(strange_background_nConst_nRhoCut, "Strange background jets", "l");
    legend8->Draw();
    
    TCanvas *c3 = new TCanvas("c3", "N_{#rho} strange variables distributions", 2500, 2500);
    c3->Divide(2,2);

    c3->cd(1);
    strange_background_nRho->SetTitle("Strange jets: N_{#rho} distribution;N_{#rho};Events");
    strange_background_nRho->GetYaxis()->SetRangeUser(std::min(strange_background_nRho->GetMinimum(), strange_signal_nRho->GetMinimum()) * 1.5, std::max(strange_background_nRho->GetMaximum(), strange_signal_nRho->GetMaximum()) * 1.5);
    strange_background_nRho->SetLineColor(kRed);
    strange_background_nRho->SetLineStyle(1);
    strange_background_nRho->DrawCopy();

    strange_signal_nRho->SetLineColor(kBlue);
    strange_signal_nRho->SetLineStyle(1);
    strange_signal_nRho->DrawCopy("same");

    TLegend *legend9 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend9->AddEntry(strange_signal_nRho, "Strange signal jets", "l");
    legend9->AddEntry(strange_background_nRho, "Strange background jets", "l");
    legend9->Draw();

    std::string title3 ="Strange jets: N_{#rho} distribution (p_{T} > " + oss.str() + " GeV/c);N_{#rho};Events";

    c3->cd(2);
    strange_background_nRho_pTCut->SetTitle(title3.c_str());
    strange_background_nRho_pTCut->GetYaxis()->SetRangeUser(std::min(strange_background_nRho_pTCut->GetMinimum(), strange_signal_nRho_pTCut->GetMinimum()) * 1.5, std::max(strange_background_nRho_pTCut->GetMaximum(), strange_signal_nRho_pTCut->GetMaximum()) * 1.5);
    strange_background_nRho_pTCut->SetLineColor(kRed);
    strange_background_nRho_pTCut->SetLineStyle(1);
    strange_background_nRho_pTCut->DrawCopy();

    strange_signal_nRho_pTCut->SetLineColor(kBlue);
    strange_signal_nRho_pTCut->SetLineStyle(1);
    strange_signal_nRho_pTCut->DrawCopy("same");

    TLegend *legend10 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend10->AddEntry(strange_signal_nRho_pTCut, "Strange signal jets", "l");
    legend10->AddEntry(strange_background_nRho_pTCut, "Strange background jets", "l");
    legend10->Draw();

    c3->cd(3);
    strange_background_nRho_nConstCut->SetTitle(("Strange jets: N_{#rho} distribution (N_{const} > " + std::to_string(nConstCut) + ");N_{#rho};Events").c_str());
    strange_background_nRho_nConstCut->GetYaxis()->SetRangeUser(std::min(strange_background_nRho_nConstCut->GetMinimum(), strange_signal_nRho_nConstCut->GetMinimum()) * 1.5, std::max(strange_background_nRho_nConstCut->GetMaximum(), strange_signal_nRho_nConstCut->GetMaximum()) * 1.5);
    strange_background_nRho_nConstCut->SetLineColor(kRed);
    strange_background_nRho_nConstCut->SetLineStyle(1);
    strange_background_nRho_nConstCut->DrawCopy();

    strange_signal_nRho_nConstCut->SetLineColor(kBlue);
    strange_signal_nRho_nConstCut->SetLineStyle(1);
    strange_signal_nRho_nConstCut->DrawCopy("same");

    TLegend *legend11 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend11->AddEntry(strange_signal_nRho_nConstCut, "Strange signal jets", "l");
    legend11->AddEntry(strange_background_nRho_nConstCut, "Strange background jets", "l");
    legend11->Draw();

    c3->cd(4);
    strange_background_nRho_nRhoCut->SetTitle(("Strange jets: N_{#rho} distribution (N_{#rho} > " + std::to_string(nRhoCut) + " GeV/c);N_{#rho};Events").c_str());
    strange_background_nRho_nRhoCut->GetYaxis()->SetRangeUser(std::min(strange_background_nRho_nRhoCut->GetMinimum(), strange_signal_nRho_nRhoCut->GetMinimum()) * 1.5, std::max(strange_background_nRho_nRhoCut->GetMaximum(), strange_signal_nRho_nRhoCut->GetMaximum()) * 1.5);
    strange_background_nRho_nRhoCut->SetLineColor(kRed);
    strange_background_nRho_nRhoCut->SetLineStyle(1);
    strange_background_nRho_nRhoCut->DrawCopy();

    strange_signal_nRho_nRhoCut->SetLineColor(kBlue);
    strange_signal_nRho_nRhoCut->SetLineStyle(1);
    strange_signal_nRho_nRhoCut->DrawCopy("same");

    TLegend *legend12 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend12->AddEntry(strange_signal_nRho_nRhoCut, "Strange signal jets", "l");
    legend12->AddEntry(strange_background_nRho_nRhoCut, "Strange background jets", "l");
    legend12->Draw();

    //=====================================================================
    // Charm jets canvases
    //=====================================================================

    TCanvas *c4 = new TCanvas("c4", "Transverse momentum charm variables distributions", 2500, 2500);
    c4->Divide(2, 2);

    c4->cd(1);
    charm_background_pT->SetTitle("Charm jets: p_{T} distribution;p_{T} [GeV/c];Events");
    charm_background_pT->GetYaxis()->SetRangeUser(std::min(charm_background_pT->GetMinimum(), charm_signal_pT->GetMinimum()) * 1.5, std::max(charm_background_pT->GetMaximum(), charm_signal_pT->GetMaximum()) * 1.5);
    charm_background_pT->SetLineColor(kRed);
    charm_background_pT->DrawCopy();

    charm_signal_pT->SetLineColor(kBlue);
    charm_signal_pT->DrawCopy("same");

    TLegend *legend13 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend13->AddEntry(charm_signal_pT, "Charm signal jets", "l");
    legend13->AddEntry(charm_background_pT, "Charm background jets", "l");
    legend13->Draw();

    std::string title4 ="Charm jets: p_{T} distribution (p_{T} > " + oss.str() + " GeV/c);p_{T} [GeV/c];Events";

    c4->cd(2);
    charm_background_pT_pTCut->SetTitle(title4.c_str());
    charm_background_pT_pTCut->GetYaxis()->SetRangeUser(std::min(charm_background_pT_pTCut->GetMinimum(), charm_signal_pT_pTCut->GetMinimum()) * 1.5, std::max(charm_background_pT_pTCut->GetMaximum(), charm_signal_pT_pTCut->GetMaximum()) * 1.5);
    charm_background_pT_pTCut->SetLineColor(kRed);
    charm_background_pT_pTCut->DrawCopy();

    charm_signal_pT_pTCut->SetLineColor(kBlue);
    charm_signal_pT_pTCut->DrawCopy("same");

    TLegend *legend14 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend14->AddEntry(charm_signal_pT_pTCut, "Charm signal jets", "l");
    legend14->AddEntry(charm_background_pT_pTCut, "Charm background jets", "l");
    legend14->Draw();

    c4->cd(3);
    charm_background_pT_nConstCut->SetTitle(("Charm jets: p_{T} distribution (N_{const} > " + std::to_string(nConstCut) + ");p_{T} [GeV/c];Events").c_str());
    charm_background_pT_nConstCut->GetYaxis()->SetRangeUser(std::min(charm_background_pT_nConstCut->GetMinimum(), charm_signal_pT_nConstCut->GetMinimum()) * 1.5, std::max(charm_background_pT_nConstCut->GetMaximum(), charm_signal_pT_nConstCut->GetMaximum()) * 1.5);
    charm_background_pT_nConstCut->SetLineColor(kRed);
    charm_background_pT_nConstCut->DrawCopy();

    charm_signal_pT_nConstCut->SetLineColor(kBlue);
    charm_signal_pT_nConstCut->DrawCopy("same");

    TLegend *legend15 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend15->AddEntry(charm_signal_pT_nConstCut, "Charm signal jets", "l");
    legend15->AddEntry(charm_background_pT_nConstCut, "Charm background jets", "l");
    legend15->Draw();

    c4->cd(4);
    charm_background_pT_nRhoCut->SetTitle(("Charm jets: p_{T} distribution (N_{#rho} > " + std::to_string(nRhoCut) + ");p_{T} [GeV/c];Events").c_str());
    charm_background_pT_nRhoCut->GetYaxis()->SetRangeUser(std::min(charm_background_pT_nRhoCut->GetMinimum(), charm_signal_pT_nRhoCut->GetMinimum()) * 1.5, std::max(charm_background_pT_nRhoCut->GetMaximum(), charm_signal_pT_nRhoCut->GetMaximum()) * 1.5);
    charm_background_pT_nRhoCut->SetLineColor(kRed);
    charm_background_pT_nRhoCut->DrawCopy();

    charm_signal_pT_nRhoCut->SetLineColor(kBlue);
    charm_signal_pT_nRhoCut->DrawCopy("same");

    TLegend *legend16 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend16->AddEntry(charm_signal_pT_nRhoCut, "Charm signal jets", "l");
    legend16->AddEntry(charm_background_pT_nRhoCut, "Charm background jets", "l");
    legend16->Draw();

    //=====================================================================

    TCanvas *c5 = new TCanvas("c5", "Number of constituents charm variables distributions", 2500, 2500);
    c5->Divide(2, 2);

    c5->cd(1);
    charm_background_nConst->SetTitle("Charm jets: number of constituents distribution;N_{const};Events");
    charm_background_nConst->GetYaxis()->SetRangeUser(std::min(charm_background_nConst->GetMinimum(), charm_signal_nConst->GetMinimum()) * 1.5, std::max(charm_background_nConst->GetMaximum(), charm_signal_nConst->GetMaximum()) * 1.5);
    charm_background_nConst->SetLineColor(kRed);
    charm_background_nConst->DrawCopy();

    charm_signal_nConst->SetLineColor(kBlue);
    charm_signal_nConst->DrawCopy("same");

    TLegend *legend17 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend17->AddEntry(charm_signal_nConst, "Charm signal jets", "l");
    legend17->AddEntry(charm_background_nConst, "Charm background jets", "l");
    legend17->Draw();

    std::string title5 ="Charm jets: number of constituents distribution (p_{T} > " + oss.str() + " GeV/c);N_{const};Events";

    c5->cd(2);
    charm_background_nConst_pTCut->SetTitle(title5.c_str());
    charm_background_nConst_pTCut->GetYaxis()->SetRangeUser(std::min(charm_background_nConst_pTCut->GetMinimum(), charm_signal_nConst_pTCut->GetMinimum()) * 1.5, std::max(charm_background_nConst_pTCut->GetMaximum(), charm_signal_nConst_pTCut->GetMaximum()) * 1.5);
    charm_background_nConst_pTCut->SetLineColor(kRed);
    charm_background_nConst_pTCut->DrawCopy();

    charm_signal_nConst_pTCut->SetLineColor(kBlue);
    charm_signal_nConst_pTCut->DrawCopy("same");

    TLegend *legend18 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend18->AddEntry(charm_signal_nConst_pTCut, "Charm signal jets", "l");
    legend18->AddEntry(charm_background_nConst_pTCut, "Charm background jets", "l");
    legend18->Draw();

    c5->cd(3);
    charm_background_nConst_nConstCut->SetTitle(("Charm jets: number of constituents distribution (N_{const} > " + std::to_string(nConstCut) + ");N_{const};Events").c_str());
    charm_background_nConst_nConstCut->GetYaxis()->SetRangeUser(std::min(charm_background_nConst_nConstCut->GetMinimum(), charm_signal_nConst_nConstCut->GetMinimum()) * 1.5, std::max(charm_background_nConst_nConstCut->GetMaximum(), charm_signal_nConst_nConstCut->GetMaximum()) * 1.5);
    charm_background_nConst_nConstCut->SetLineColor(kRed);
    charm_background_nConst_nConstCut->DrawCopy();

    charm_signal_nConst_nConstCut->SetLineColor(kBlue);
    charm_signal_nConst_nConstCut->DrawCopy("same");

    TLegend *legend19 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend19->AddEntry(charm_signal_nConst_nConstCut, "Charm signal jets", "l");
    legend19->AddEntry(charm_background_nConst_nConstCut, "Charm background jets", "l");
    legend19->Draw();

    c5->cd(4);
    charm_background_nConst_nRhoCut->SetTitle(("Charm jets: number of constituents distribution (N_{#rho} > " + std::to_string(nRhoCut) + ");N_{const};Events").c_str());
    charm_background_nConst_nRhoCut->GetYaxis()->SetRangeUser(std::min(charm_background_nConst_nRhoCut->GetMinimum(), charm_signal_nConst_nRhoCut->GetMinimum()) * 1.5, std::max(charm_background_nConst_nRhoCut->GetMaximum(), charm_signal_nConst_nRhoCut->GetMaximum()) * 1.5);
    charm_background_nConst_nRhoCut->SetLineColor(kRed);
    charm_background_nConst_nRhoCut->DrawCopy();

    charm_signal_nConst_nRhoCut->SetLineColor(kBlue);
    charm_signal_nConst_nRhoCut->DrawCopy("same");

    TLegend *legend20 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend20->AddEntry(charm_signal_nConst_nRhoCut, "Charm signal jets", "l");
    legend20->AddEntry(charm_background_nConst_nRhoCut, "Charm background jets", "l");
    legend20->Draw();

    //=====================================================================

    TCanvas *c6 = new TCanvas("c6", "N_{#rho} charm variables distributions", 2500, 2500);
    c6->Divide(2, 2);

    c6->cd(1);
    charm_background_nRho->SetTitle("Charm jets: N_{#rho} distribution;N_{#rho};Events");
    charm_background_nRho->SetLineColor(kRed);
    charm_background_nRho->GetYaxis()->SetRangeUser(std::min(charm_background_nRho->GetMinimum(), charm_signal_nRho->GetMinimum()) * 1.5, std::max(charm_background_nRho->GetMaximum(), charm_signal_nRho->GetMaximum()) * 1.5);
    charm_background_nRho->DrawCopy();

    charm_signal_nRho->SetLineColor(kBlue);
    charm_signal_nRho->DrawCopy("same");

    TLegend *legend21 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend21->AddEntry(charm_signal_nRho, "Charm signal jets", "l");
    legend21->AddEntry(charm_background_nRho, "Charm background jets", "l");
    legend21->Draw();

    std::string title6 = "Charm jets: N_{#rho} distribution (p_{T} > " + std::to_string(pTCut) + " GeV/c);N_{#rho};Events";

    c6->cd(2);
    charm_background_nRho_pTCut->SetTitle(title6.c_str());
    charm_background_nRho_pTCut->SetLineColor(kRed);
    charm_background_nRho_pTCut->GetYaxis()->SetRangeUser(std::min(charm_background_nRho_pTCut->GetMinimum(), charm_signal_nRho_pTCut->GetMinimum()) * 1.5, std::max(charm_background_nRho_pTCut->GetMaximum(), charm_signal_nRho_pTCut->GetMaximum()) * 1.5);
    charm_background_nRho_pTCut->DrawCopy();

    charm_signal_nRho_pTCut->SetLineColor(kBlue);
    charm_signal_nRho_pTCut->DrawCopy("same");

    TLegend *legend22 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend22->AddEntry(charm_signal_nRho_pTCut, "Charm signal jets", "l");
    legend22->AddEntry(charm_background_nRho_pTCut, "Charm background jets", "l");
    legend22->Draw();

    c6->cd(3);
    charm_background_nRho_nConstCut->SetTitle(("Charm jets: N_{#rho} distribution (N_{const} > " + std::to_string(nConstCut) + ");N_{#rho};Events").c_str());
    charm_background_nRho_nConstCut->GetYaxis()->SetRangeUser(std::min(charm_background_nRho_nConstCut->GetMinimum(), charm_signal_nRho_nConstCut->GetMinimum()) * 1.5, std::max(charm_background_nRho_nConstCut->GetMaximum(), charm_signal_nRho_nConstCut->GetMaximum()) * 1.5);
    charm_background_nRho_nConstCut->DrawCopy();

    charm_signal_nRho_nConstCut->SetLineColor(kBlue);
    charm_signal_nRho_nConstCut->DrawCopy("same");

    TLegend *legend23 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend23->AddEntry(charm_signal_nRho_nConstCut, "Charm signal jets", "l");
    legend23->AddEntry(charm_background_nRho_nConstCut, "Charm background jets", "l");
    legend23->Draw();

    c6->cd(4);
    charm_background_nRho_nRhoCut->SetTitle(("Charm jets: N_{#rho} distribution (N_{#rho} > " + std::to_string(nRhoCut) + ");N_{#rho};Events").c_str());
    charm_background_nRho_nRhoCut->GetYaxis()->SetRangeUser(std::min(charm_background_nRho_nRhoCut->GetMinimum(), charm_signal_nRho_nRhoCut->GetMinimum()) * 1.5, std::max(charm_background_nRho_nRhoCut->GetMaximum(), charm_signal_nRho_nRhoCut->GetMaximum()) * 1.5);
    charm_background_nRho_nRhoCut->DrawCopy();

    charm_signal_nRho_nRhoCut->SetLineColor(kBlue);
    charm_signal_nRho_nRhoCut->DrawCopy("same");

    TLegend *legend24 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend24->AddEntry(charm_signal_nRho_nRhoCut, "Charm signal jets", "l");
    legend24->AddEntry(charm_background_nRho_nRhoCut, "Charm background jets", "l");
    legend24->Draw();

    //=====================================================================

    TCanvas *c7 = new TCanvas("c7", "Jet vertices invariant masses charm distributions", 2500, 2500);
    c7->Divide(1, 1);

    c7->cd(1);
    charm_background_jetVerticesInvariantMasses->SetTitle("Charm jets: jet vertices invariant masses distribution;Jet vertices invariant mass (GeV/c^{2});Events");
    charm_background_jetVerticesInvariantMasses->GetYaxis()->SetRangeUser(std::min(charm_background_jetVerticesInvariantMasses->GetMinimum(), charm_signal_jetVerticesInvariantMasses->GetMinimum()) * 1.5, std::max(charm_background_jetVerticesInvariantMasses->GetMaximum(), charm_signal_jetVerticesInvariantMasses->GetMaximum()) * 1.5);
    charm_background_jetVerticesInvariantMasses->DrawCopy();

    charm_signal_jetVerticesInvariantMasses->SetLineColor(kBlue);
    charm_signal_jetVerticesInvariantMasses->DrawCopy("same");

    TLegend *legend25 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend25->AddEntry(charm_signal_jetVerticesInvariantMasses, "Charm signal jets", "l");
    legend25->AddEntry(charm_background_jetVerticesInvariantMasses, "Charm background jets", "l");
    legend25->Draw();

    /*
    c7->cd(2);
    charm_background_jetVerticesInvariantMasses_pTCut->SetTitle(("Charm jets: jet vertices invariant masses distribution (p_{T} > " + std::to_string(pTCut) + " GeV/c);Jet vertices invariant mass (GeV/c^{2});Events").c_str());
    
    charm_background_jetVerticesInvariantMasses_pTCut->SetLineColor(kRed);
    charm_background_jetVerticesInvariantMasses_pTCut->DrawCopy();

    charm_signal_jetVerticesInvariantMasses_pTCut->SetLineColor(kBlue);
    charm_signal_jetVerticesInvariantMasses_pTCut->DrawCopy("same");

    TLegend *legend26 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend26->AddEntry(charm_signal_jetVerticesInvariantMasses_pTCut, "Charm signal jets", "l");
    legend26->AddEntry(charm_background_jetVerticesInvariantMasses_pTCut, "Charm background jets", "l");
    legend26->Draw();

    c7->cd(3);
    charm_background_jetVerticesInvariantMasses_nConstCut->SetTitle(("Charm jets: jet vertices invariant masses distribution (N_{const} > " + std::to_string(nConstCut) + ");Jet vertices invariant mass (GeV/c^{2});Events").c_str());
    charm_background_jetVerticesInvariantMasses_nConstCut->SetLineColor(kRed);
    charm_background_jetVerticesInvariantMasses_nConstCut->DrawCopy();

    charm_signal_jetVerticesInvariantMasses_nConstCut->SetLineColor(kBlue);
    charm_signal_jetVerticesInvariantMasses_nConstCut->DrawCopy("same");

    TLegend *legend27 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend27->AddEntry(charm_signal_jetVerticesInvariantMasses_nConstCut, "Charm signal jets", "l");
    legend27->AddEntry(charm_background_jetVerticesInvariantMasses_nConstCut, "Charm background jets", "l");
    legend27->Draw();

    c7->cd(4);
    charm_background_jetVerticesInvariantMasses_nRhoCut->SetTitle(("Charm jets: jet vertices invariant masses distribution (N_{#rho} > " + std::to_string(nRhoCut) + ");Jet vertices invariant mass (GeV/c^{2});Events").c_str());
    charm_background_jetVerticesInvariantMasses_nRhoCut->SetLineColor(kRed);
    charm_background_jetVerticesInvariantMasses_nRhoCut->DrawCopy();

    charm_signal_jetVerticesInvariantMasses_nRhoCut->SetLineColor(kBlue);
    charm_signal_jetVerticesInvariantMasses_nRhoCut->DrawCopy("same");

    TLegend *legend28 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend28->AddEntry(charm_signal_jetVerticesInvariantMasses_nRhoCut, "Charm signal jets", "l");
    legend28->AddEntry(charm_background_jetVerticesInvariantMasses_nRhoCut, "Charm background jets", "l");
    legend28->Draw();
    */ // Not a discriminatory variable, cuts not applied
}
