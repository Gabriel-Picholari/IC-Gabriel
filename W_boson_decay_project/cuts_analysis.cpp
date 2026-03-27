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

    TFile *s_file = TFile::Open(strange_file, "READ");

    TTree *signalTree_s = dynamic_cast<TTree *>(s_file->Get("SignalTree_s"));
    signalTree_s->SetBranchAddress("pT_s", &pT_s);
    signalTree_s->SetBranchAddress("nConst_s", &nConst_s);
    signalTree_s->SetBranchAddress("nRho_s", &nRho_s);
    signalTree_s->SetBranchAddress("eventID_s", &eventID_s);

    TTree *backgroundTree_s = dynamic_cast<TTree *>(s_file->Get("BackgroundTree_s")); 
    backgroundTree_s->SetBranchAddress("pT_s", &pT_s);
    backgroundTree_s->SetBranchAddress("nConst_s", &nConst_s);
    backgroundTree_s->SetBranchAddress("nRho_s", &nRho_s);
    backgroundTree_s->SetBranchAddress("eventID_s", &eventID_s);

    Long64_t ne_signal_s = signalTree_s->GetEntries();
    Long64_t ne_background_s = backgroundTree_s->GetEntries();

    //---------------------------------------------------------------------------------------------------------

    TFile *c_file = TFile::Open(charm_file, "READ");

    TTree *signalTree_c = dynamic_cast<TTree *>(c_file->Get("SignalTree_c"));
    signalTree_c->SetBranchAddress("pT_c", &pT_c);
    signalTree_c->SetBranchAddress("nConst_c", &nConst_c);
    signalTree_c->SetBranchAddress("nRho_c", &nRho_c);
    signalTree_c->SetBranchAddress("eventID_c", &eventID_c);

    TTree *backgroundTree_c = dynamic_cast<TTree *>(c_file->Get("BackgroundTree_c"));
    backgroundTree_c->SetBranchAddress("pT_c", &pT_c);
    backgroundTree_c->SetBranchAddress("nConst_c", &nConst_c);
    backgroundTree_c->SetBranchAddress("nRho_c", &nRho_c);
    backgroundTree_c->SetBranchAddress("eventID_c", &eventID_c);

    Long64_t ne_signal_c = signalTree_c->GetEntries();
    Long64_t ne_background_c = backgroundTree_c->GetEntries();

    //---------------------------------------------------------------------------------------------------------
    // Analysis
    //---------------------------------------------------------------------------------------------------------

    // From my perspective the optimal approach to dealing with this pre analysis is to go through each of the trees separately, build the histograms and then plot
    // them in conjoint. Thats what will be done from now on!

    TH1F *strange_signal_pT = new TH1F("strange_signal_pT", "p_{T} distribution of strange jets; p_{T} [GeV/c]; Events", 100, 0, 100);
    TH1F *strange_signal_nConst = new TH1F("strange_signal_nConst", "Number of constituents distribution of strange jets; nConst; Events", 50, 0, 50);
    TH1F *strange_signal_nRho = new TH1F("strange_signal_nRho", "N_{#rho} distribution of strange jets; N_{#rho}; Events", 10, 0, 10);

    TH1F *strange_background_pT = new TH1F("strange_background_pT", "p_{T} distribution of non-strange jets; p_{T} [GeV/c]; Events", 100, 0, 100);
    TH1F *strange_background_nConst = new TH1F("strange_background_nConst", "Number of constituents distribution of non-strange jets; nConst; Events", 50, 0, 50);
    TH1F *strange_background_nRho = new TH1F("strange_background_nRho", "N_{#rho} distribution of non-strange jets; N_{#rho}; Events", 10, 0, 10);

    for (Long64_t i = 0; i < ne_signal_s; i++)
    {
        signalTree_s->GetEntry(i);
        strange_signal_pT->Fill(pT_s);
        strange_signal_nConst->Fill(nConst_s);
        strange_signal_nRho->Fill(nRho_s);
    }

    for (Long64_t i = 0; i < ne_background_s; i++)
    {
        backgroundTree_s->GetEntry(i);
        strange_background_pT->Fill(pT_s);
        strange_background_nConst->Fill(nConst_s);
        strange_background_nRho->Fill(nRho_s);
    }

    //---------------------------------------------------------------------------------------------------------

    TH1F *charm_signal_pT = new TH1F("charm_signal_pT", "p_{T} distribution of charm jets; p_{T} [GeV/c]; Events", 100, 0, 100);
    TH1F *charm_signal_nConst = new TH1F("charm_signal_nConst", "Number of constituents distribution of charm jets; nConst; Events", 50, 0, 50);
    TH1F *charm_signal_nRho = new TH1F("charm_signal_nRho", "N_{#rho} distribution of charm jets; N_{#rho}; Events", 10, 0, 10);

    TH1F *charm_background_pT = new TH1F("charm_background_pT", "p_{T} distribution of non-charm jets; p_{T} [GeV/c]; Events", 100, 0, 100);
    TH1F *charm_background_nConst = new TH1F("charm_background_nConst", "Number of constituents distribution of non-charm jets; nConst; Events", 50, 0, 50);
    TH1F *charm_background_nRho = new TH1F("charm_background_nRho", "N_{#rho} distribution of non-charm jets; N_{#rho}; Events", 10, 0, 10);

    for (Long64_t i = 0; i < ne_signal_c; i++)
    {
        signalTree_c->GetEntry(i);
        charm_signal_pT->Fill(pT_c);
        charm_signal_nConst->Fill(nConst_c);
        charm_signal_nRho->Fill(nRho_c);
    }

    for (Long64_t i = 0; i < ne_background_c; i++)
    {
        backgroundTree_c->GetEntry(i);
        charm_background_pT->Fill(pT_c);
        charm_background_nConst->Fill(nConst_c);
        charm_background_nRho->Fill(nRho_c);

    }

    // The background distributions are much more populated than the signal ones. Then, to make the signal and background distributions comparable, we will rescale
    // the background histograms accordingly.

    Float_t scaleS = (Float_t)(ne_signal_s) / (Float_t)(ne_background_s);
    strange_background_pT->Scale(scaleS);
    strange_background_nConst->Scale(scaleS);
    strange_background_nRho->Scale(scaleS);

    Float_t scaleC = (Float_t)(ne_signal_c) / (Float_t)(ne_background_c);
    charm_background_pT->Scale(scaleC);
    charm_background_nConst->Scale(scaleC);
    charm_background_nRho->Scale(scaleC);

TCanvas *c1 = new TCanvas("c1", "Discriminatory variables distributions", 2500, 2500);
    c1->Divide(3, 2);

    c1->cd(1);
    strange_background_pT->SetTitle("Strange jets: p_{T} distribution;p_{T} [GeV/c];Events");
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

    c1->cd(2);
    strange_background_nConst->SetTitle("Strange jets: number of constituents distribution;N_{const};Events");
    strange_background_nConst->SetLineColor(kRed);
    strange_background_nConst->SetLineStyle(1);
    strange_background_nConst->DrawCopy();

    strange_signal_nConst->SetLineColor(kBlue);
    strange_signal_nConst->SetLineStyle(1);
    strange_signal_nConst->DrawCopy("same");

    TLegend *legend2 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend2->AddEntry(strange_signal_nConst, "Strange signal jets", "l");
    legend2->AddEntry(strange_background_nConst, "Strange background jets", "l");
    legend2->Draw();

    c1->cd(3);
    strange_background_nRho->SetTitle("Strange jets: N_{#rho} distribution;N_{#rho};Events");
    strange_background_nRho->SetLineColor(kRed);
    strange_background_nRho->SetLineStyle(1);
    strange_background_nRho->DrawCopy();

    strange_signal_nRho->SetLineColor(kBlue);
    strange_signal_nRho->SetLineStyle(1);
    strange_signal_nRho->DrawCopy("same");

    TLegend *legend3 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend3->AddEntry(strange_signal_nRho, "Strange signal jets", "l");
    legend3->AddEntry(strange_background_nRho, "Strange background jets", "l");
    legend3->Draw();

    c1->cd(4);
    charm_background_pT->SetTitle("Charm jets: p_{T} distribution;p_{T} [GeV/c];Events");
    charm_background_pT->SetLineColor(kRed);
    charm_background_pT->SetLineStyle(1);
    charm_background_pT->DrawCopy();

    charm_signal_pT->SetLineColor(kBlue);
    charm_signal_pT->SetLineStyle(1);
    charm_signal_pT->DrawCopy("same");

    TLegend *legend4 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend4->AddEntry(charm_signal_pT, "Charm signal jets", "l");
    legend4->AddEntry(charm_background_pT, "Charm background jets", "l");
    legend4->Draw();

    c1->cd(5);
    charm_background_nConst->SetTitle("Charm jets: number of constituents distribution;N_{const};Events");
    charm_background_nConst->SetLineColor(kRed);
    charm_background_nConst->SetLineStyle(1);
    charm_background_nConst->DrawCopy();

    charm_signal_nConst->SetLineColor(kBlue);
    charm_signal_nConst->SetLineStyle(1);
    charm_signal_nConst->DrawCopy("same");

    TLegend *legend5 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend5->AddEntry(charm_signal_nConst, "Charm signal jets", "l");
    legend5->AddEntry(charm_background_nConst, "Charm background jets", "l");
    legend5->Draw();

    c1->cd(6);
    charm_background_nRho->SetTitle("Charm jets: N_{#rho} distribution;N_{#rho};Events");
    charm_background_nRho->SetLineColor(kRed);
    charm_background_nRho->SetLineStyle(1);
    charm_background_nRho->DrawCopy();

    charm_signal_nRho->SetLineColor(kBlue);
    charm_signal_nRho->SetLineStyle(1);
    charm_signal_nRho->DrawCopy("same");

    TLegend *legend6 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend6->AddEntry(charm_signal_nRho, "Charm signal jets", "l");
    legend6->AddEntry(charm_background_nRho, "Charm background jets", "l");
    legend6->Draw();
}
