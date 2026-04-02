//===========================================================================================================================================================================
// Classical purity and efficiency analysis:
// This macro intends to study the performance of varios combinations of cuts on the discriminatory variables through the analysis of classical purity and efficiency
// distributions.
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

void purity_efficiency_analysis(std::string switch_string, const char* strange_file)
{
    gSystem->Load("libEG");
    gSystem->Load("libEGPythia8");

    //---------------------------------------------------------------------------------------------------------
    // Initialization of the .root file and TTrees
    //---------------------------------------------------------------------------------------------------------

    Float_t pT, nConst, nRho, eventID;
    TTree *signalTree;
    TTree *backgroundTree;

    std::string sfx;

    if (switch_string == "strange") 
    {
        sfx = "_s";
    } 
    else if (switch_string == "charm")
    {
        sfx = "_c";
    }

    TFile *file = TFile::Open(strange_file, "READ");

    
    signalTree = dynamic_cast<TTree *>(file->Get(("SignalTree" + sfx).c_str()));
    signalTree->SetBranchAddress(("pT" + sfx).c_str(), &pT);
    signalTree->SetBranchAddress(("nConst" + sfx).c_str(), &nConst);
    signalTree->SetBranchAddress(("nRho" + sfx).c_str(), &nRho);
    signalTree->SetBranchAddress(("eventID" + sfx).c_str(), &eventID);

    backgroundTree = dynamic_cast<TTree *>(file->Get(("BackgroundTree" + sfx).c_str()));
    backgroundTree->SetBranchAddress(("pT" + sfx).c_str(), &pT);
    backgroundTree->SetBranchAddress(("nConst" + sfx).c_str(), &nConst);
    backgroundTree->SetBranchAddress(("nRho" + sfx).c_str(), &nRho);
    backgroundTree->SetBranchAddress(("eventID" + sfx).c_str(), &eventID);

    Long64_t ne_signal = signalTree->GetEntries();
    Long64_t ne_background = backgroundTree->GetEntries();

    //---------------------------------------------------------------------------------------------------------
    // Histograms
    //---------------------------------------------------------------------------------------------------------

    TH1F* pT_cut_purity = new TH1F("pT_cut_purity", "Purity for pT cuts; pT cut (GeV/c); Purity", 60, 0, 60);
    TH1F* pT_cut_efficiency = new TH1F("pT_cut_efficiency", "Efficiency for pT cuts; pT cut (GeV/c); Efficiency", 40, 0, 40);

    //---------------------------------------------------------------------------------------------------------
    // Analysis
    //---------------------------------------------------------------------------------------------------------

    std::vector<Float_t> signal_pT, signal_nConst, signal_nRho;

    for (Long64_t i = 0; i < ne_signal; i++)
    {
        signalTree->GetEntry(i);
        signal_pT.push_back(pT);
        signal_nRho.push_back(nRho);
        signal_nConst.push_back(nConst);
    }

    std::vector<Float_t> background_pT, background_nConst, background_nRho;
    for (Long64_t i = 0; i < ne_background; i++)
    {
        backgroundTree->GetEntry(i);
        background_pT.push_back(pT);
        background_nRho.push_back(nRho);
        background_nConst.push_back(nConst);
    }
    
    Float_t pT_maxThreshold = 40;
    Float_t nConstCut = 20;
    Float_t nRhoCut = 1;

    for (Float_t pT_thr = 0; pT_thr <= pT_maxThreshold; pT_thr += 0.01)
    {
        Int_t TP = 0, FN = 0, TN = 0, FP = 0;
        for (Long64_t i = 0; i < signal_pT.size(); i++)
        {
            if (signal_pT[i] >= pT_thr)
            {
                TP++;
            }
            else
            {
                FN++;
            }
        }
        for (Long64_t i = 0; i < background_pT.size(); i++)
        {
            if (background_pT[i] >= pT_thr)
            {
                FP++;
            }
            else
            {
                TN++;
            }
        }

        Float_t purity = (Float_t)TP / (TP + FP);
        Float_t efficiency = (Float_t)TP / (TP + FN);

        pT_cut_purity->Fill(pT_thr, purity);
        pT_cut_efficiency->Fill(pT_thr, efficiency);
    }

    TCanvas *c1 = new TCanvas("c1", "Purity and efficiency for pT cuts only", 2500, 2500);
    c1->Divide(1,1);

    c1->cd(1);
    pT_cut_efficiency->SetTitle(("Purity and Efficiency for pT cuts only (" + switch_string + " jets)").c_str());
    pT_cut_efficiency->GetXaxis()->SetTitle("pT cut (GeV/c)");
    pT_cut_efficiency->GetYaxis()->SetTitle("Purity and Efficiency");
    pT_cut_efficiency->SetLineColor(kRed);
    pT_cut_efficiency->DrawCopy("HIST");

    pT_cut_purity->SetLineColor(kBlue);
    pT_cut_purity->DrawCopy("HIST SAME");

}
