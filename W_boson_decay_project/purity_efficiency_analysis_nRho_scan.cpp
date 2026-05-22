//===========================================================================================================================================================================
// Classical purity and efficiency analysis:
// This macro intends to study the performance of varios combinations of cuts on the discriminatory variables through the analysis of classical purity and efficiency
// distributions. The scan is performed on the N_rho variable.
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

void purity_efficiency_analysis_nRho_scan(const char* strange_file, std::string switch_string)
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

    TH1F* nRho_cut_purity = new TH1F("nRho_cut_purity", "Purity for nRho scan; nRho cut (GeV/c); Purity", 7, 0, 7);
    TH1F* nRho_cut_efficiency = new TH1F("nRho_cut_efficiency", "Efficiency for nRho scan; nRho cut (GeV/c); Efficiency", 7, 0, 7);


    TH1F* nRho_nConst_cut_purity1 = new TH1F("nRho_nConst_cut_purity1", "Purity for nRho scan with nConst cut (first cut); nRho cut (GeV/c); Purity", 7, 0, 7);
    TH1F* nRho_nConst_cut_efficiency1 = new TH1F("nRho_nConst_cut_efficiency1", "Efficiency for nRho scan with nConst cut (first cut); nRho cut (GeV/c); Efficiency", 7, 0, 7);

    TH1F* nRho_nConst_cut_purity2 = new TH1F("nRho_nConst_cut_purity2", "Purity for nRho scan with nConst cut (second cut); nRho cut (GeV/c); Purity", 7, 0, 7);
    TH1F* nRho_nConst_cut_efficiency2 = new TH1F("nRho_nConst_cut_efficiency2", "Efficiency for nRho scan with nConst cut (second cut); nRho cut (GeV/c); Efficiency", 7, 0, 7);


    TH1F* nRho_pT_cut_purity1 = new TH1F("nRho_pT_cut_purity1", "Purity for nRho scan with pT cut (first cut); nRho cut (GeV/c); Purity", 7, 0, 7);
    TH1F* nRho_pT_cut_efficiency1 = new TH1F("nRho_pT_cut_efficiency1", "Efficiency for nRho scan with pT cut (first cut); nRho cut (GeV/c); Efficiency", 5, 0, 5);

    TH1F* nRho_pT_cut_purity2 = new TH1F("nRho_pT_cut_purity2", "Purity for nRho scan with pT cut (second cut); nRho cut (GeV/c); Purity", 7, 0, 7);
    TH1F* nRho_pT_cut_efficiency2 = new TH1F("nRho_pT_cut_efficiency2", "Efficiency for nRho scan with pT cut (second cut); nRho cut (GeV/c); Efficiency", 7, 0, 7);


    TH1F* nRho_nConst_pT_cut_purity1 = new TH1F("nRho_nConst_pT_cut_purity1", "Purity for nRho scan with nConst and pT cut (first cuts); nRho cut (GeV/c); Purity", 7, 0, 7);
    TH1F* nRho_nConst_pT_cut_efficiency1 = new TH1F("nRho_nConst_pT_cut_efficiency1", "Efficiency for nRho scan with nConst and pT cut (first cuts); nRho cut (GeV/c); Efficiency", 7, 0, 7);

    TH1F* nRho_nConst_pT_cut_purity2 = new TH1F("nRho_nConst_pT_cut_purity2", "Purity for nRho scan with nConst and pT cut (second cuts); nRho cut (GeV/c); Purity", 7, 0, 7);
    TH1F* nRho_nConst_pT_cut_efficiency2 = new TH1F("nRho_nConst_pT_cut_efficiency2", "Efficiency for nRho scan with nConst and pT cut (second cuts); nRho cut (GeV/c); Efficiency", 7, 0, 7);

    TH1F* nRho_nConst_pT_cut_purity3 = new TH1F("nRho_nConst_pT_cut_purity3", "Purity for nRho scan with nConst and pT cut (third cuts); nRho cut (GeV/c); Purity", 7, 0, 7);
    TH1F* nRho_nConst_pT_cut_efficiency3 = new TH1F("nRho_nConst_pT_cut_efficiency3", "Efficiency for nRho scan with nConst and pT cut (third cuts); nRho cut (GeV/c); Efficiency", 7, 0, 7);

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
    
    Float_t maxThreshold = 5;
    Float_t threshold_increment = 1;
    
    // Unidimensional nRho scan (no cuts included)
    for (Float_t nRho_thr = 0; nRho_thr <= maxThreshold; nRho_thr += threshold_increment)
    {
        Int_t TP = 0, FN = 0, TN = 0, FP = 0;
        for (Long64_t i = 0; i < signal_nRho.size(); i++)
        {
            if (signal_nRho[i] >= nRho_thr)
            {
                TP++;
            }
            else
            {
                FN++;
            }
        }
        for (Long64_t i = 0; i < background_nRho.size(); i++)
        {
            if (background_nRho[i] >= nRho_thr)
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

        nRho_cut_purity->Fill(nRho_thr, purity);
        nRho_cut_efficiency->Fill(nRho_thr, efficiency);
    }

    Float_t first_nConstCut = 0, second_nConstCut = 0;

    if(switch_string == "strange")
    {
        first_nConstCut = 10;
        second_nConstCut = 20;
    }
    else if(switch_string == "charm")
    {
        first_nConstCut = 10;
        second_nConstCut = 20;
    }

    // Unidimensional nRho scan with nConst cuts
    for (Float_t nRho_thr = 0; nRho_thr <= maxThreshold; nRho_thr += threshold_increment)
    {
        Int_t TP1 = 0, FN1 = 0, TN1 = 0, FP1 = 0;
        Int_t TP2 = 0, FN2 = 0, TN2 = 0, FP2 = 0;

        for (Long64_t i = 0; i < signal_nRho.size(); i++)
        {
            if (signal_nRho[i] >= nRho_thr && signal_nConst[i] >= first_nConstCut)
            {
                TP1++;
            }
            else
            {
                FN1++;
            }

            if (signal_nRho[i] >= nRho_thr && signal_nConst[i] >= second_nConstCut)
            {
                TP2++;
            }
            else
            {
                FN2++;
            }

        }

        for (Long64_t i = 0; i < background_nRho.size(); i++)
        {
            if (background_nRho[i] >= nRho_thr && background_nConst[i] >= first_nConstCut)
            {
                FP1++;
            }
            else
            {
                TN1++;
            }

            if (background_nRho[i] >= nRho_thr && background_nConst[i] >= second_nConstCut)
            {
                FP2++;
            }
            else
            {
                TN2++;
            }

        }

        Float_t purity1 = (Float_t)TP1 / (TP1 + FP1);
        Float_t efficiency1 = (Float_t)TP1 / (TP1 + FN1);

        Float_t purity2 = (Float_t)TP2 / (TP2 + FP2);
        Float_t efficiency2 = (Float_t)TP2 / (TP2 + FN2);

        nRho_nConst_cut_purity1->Fill(nRho_thr, purity1);
        nRho_nConst_cut_efficiency1->Fill(nRho_thr, efficiency1);

        nRho_nConst_cut_purity2->Fill(nRho_thr, purity2);
        nRho_nConst_cut_efficiency2->Fill(nRho_thr, efficiency2);

    }

    Float_t first_pTCut = 0, second_pTCut = 0;

    if (switch_string == "strange") 
    {
        first_pTCut = 15;
        second_pTCut = 20;
    } 
    else if (switch_string == "charm")
    {
        first_pTCut = 15;
        second_pTCut = 20;
    }


    // Unidimensional nRho scan with pT cuts
    for (Float_t nRho_thr = 0; nRho_thr <= maxThreshold; nRho_thr += threshold_increment)
    {
        Int_t TP1 = 0, FN1 = 0, TN1 = 0, FP1 = 0;
        Int_t TP2 = 0, FN2 = 0, TN2 = 0, FP2 = 0;

        for (Long64_t i = 0; i < signal_nRho.size(); i++)
        {
            if (signal_nRho[i] >= nRho_thr && signal_pT[i] >= first_pTCut)
            {
                TP1++;
            }
            else
            {
                FN1++;
            }

            if (signal_nRho[i] >= nRho_thr && signal_pT[i] >= second_pTCut)
            {
                TP2++;
            }
            else
            {
                FN2++;
            }

        }

        for (Long64_t i = 0; i < background_nRho.size(); i++)
        {
            if (background_nRho[i] >= nRho_thr && background_pT[i] >= first_pTCut)
            {
                FP1++;
            }
            else
            {
                TN1++;
            }

            if (background_nRho[i] >= nRho_thr && background_pT[i] >= second_pTCut)
            {
                FP2++;
            }
            else
            {
                TN2++;
            }

        }

        Float_t purity1 = (Float_t)TP1 / (TP1 + FP1);
        Float_t efficiency1 = (Float_t)TP1 / (TP1 + FN1);

        Float_t purity2 = (Float_t)TP2 / (TP2 + FP2);
        Float_t efficiency2 = (Float_t)TP2 / (TP2 + FN2);

        nRho_pT_cut_purity1->Fill(nRho_thr, purity1);
        nRho_pT_cut_efficiency1->Fill(nRho_thr, efficiency1);

        nRho_pT_cut_purity2->Fill(nRho_thr, purity2);
        nRho_pT_cut_efficiency2->Fill(nRho_thr, efficiency2);
    }

    // Unidimensional pT scan with nConst and nRho cuts
    for (Float_t nRho_thr = 0; nRho_thr <= maxThreshold; nRho_thr += threshold_increment)
    {
        Int_t TP1 = 0, FN1 = 0, TN1 = 0, FP1 = 0;
        Int_t TP2 = 0, FN2 = 0, TN2 = 0, FP2 = 0;
        Int_t TP3 = 0, FN3 = 0, TN3 = 0, FP3 = 0;

        for (Long64_t i = 0; i < signal_pT.size(); i++)
        {
            if (signal_pT[i] >= nRho_thr && signal_nConst[i] >= first_nConstCut && signal_pT[i] >= first_pTCut)
            {
                TP1++;
            }
            else
            {
                FN1++;
            }

            if (signal_pT[i] >= nRho_thr && signal_nConst[i] >= first_nConstCut && signal_pT[i] >= second_pTCut)
            {
                TP2++;
            }
            else
            {
                FN2++;
            }

            if (signal_pT[i] >= nRho_thr && signal_nConst[i] >= second_nConstCut && signal_pT[i] >= second_pTCut)
            {
                TP3++;
            }
            else
            {
                FN3++;
            }

        }

        for (Long64_t i = 0; i < background_pT.size(); i++)
        {
            if (background_pT[i] >= nRho_thr && background_nConst[i] >= first_nConstCut && background_pT[i] >= first_pTCut)
            {
                FP1++;
            }
            else
            {
                TN1++;
            }

            if (background_pT[i] >= nRho_thr && background_nConst[i] >= first_nConstCut && background_pT[i] >= second_pTCut)
            {
                FP2++;
            }
            else
            {
                TN2++;
            }

            if (background_pT[i] >= nRho_thr && background_nConst[i] >= second_nConstCut && background_pT[i] >= second_pTCut)
            {
                FP3++;
            }
            else
            {
                TN3++;
            }

        }

        Float_t purity1 = (Float_t)TP1 / (TP1 + FP1);
        Float_t efficiency1 = (Float_t)TP1 / (TP1 + FN1);

        Float_t purity2 = (Float_t)TP2 / (TP2 + FP2);
        Float_t efficiency2 = (Float_t)TP2 / (TP2 + FN2);

        Float_t purity3 = (Float_t)TP3 / (TP3 + FP3);
        Float_t efficiency3 = (Float_t)TP3 / (TP3 + FN3);

        nRho_nConst_pT_cut_purity1->Fill(nRho_thr, purity1);
        nRho_nConst_pT_cut_efficiency1->Fill(nRho_thr, efficiency1);

        nRho_nConst_pT_cut_purity2->Fill(nRho_thr, purity2);
        nRho_nConst_pT_cut_efficiency2->Fill(nRho_thr, efficiency2);

        nRho_nConst_pT_cut_purity3->Fill(nRho_thr, purity3);
        nRho_nConst_pT_cut_efficiency3->Fill(nRho_thr, efficiency3);

    }

    TCanvas *c1 = new TCanvas("c1", "Purity and efficiency: N_{#rho} scan", 2500, 2500);
    c1->Divide(1,1);

    c1->cd(1);
    nRho_cut_efficiency->SetTitle(("Purity and Efficiency for N_{#rho} scan (" + switch_string + " jets)").c_str());
    nRho_cut_efficiency->GetXaxis()->SetTitle("N_{#rho} threshold (GeV/c)");
    nRho_cut_efficiency->GetYaxis()->SetTitle("Purity and Efficiency");
    nRho_cut_efficiency->SetLineColor(kRed);
    nRho_cut_efficiency->DrawCopy("HIST");

    nRho_cut_purity->SetLineColor(kBlue);
    nRho_cut_purity->DrawCopy("HIST SAME");

    TLegend *legend1 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend1->AddEntry(nRho_cut_purity, "Purity", "l");
    legend1->AddEntry(nRho_cut_efficiency, "Efficiency", "l");
    legend1->Draw();

    TCanvas *c2 = new TCanvas("c2", "Purity and efficiency: N_{#rho} scan with N_{const} cuts", 2500, 2500);
    c2->Divide(1,1);

    c2->cd(2);
    nRho_nConst_cut_efficiency1->SetTitle(("Purity and Efficiency for N_{#rho} scan with N_{const} cuts (" + switch_string + " jets)").c_str());
    nRho_nConst_cut_efficiency1->GetXaxis()->SetTitle("N_{#rho} threshold (GeV/c)");
    nRho_nConst_cut_efficiency1->GetYaxis()->SetTitle("Purity and Efficiency");
    nRho_nConst_cut_efficiency1->SetLineColor(kRed);
    nRho_nConst_cut_efficiency1->DrawCopy("HIST");

    nRho_nConst_cut_purity1->SetLineColor(kBlue);
    nRho_nConst_cut_purity1->DrawCopy("HIST SAME");

    nRho_nConst_cut_efficiency2->SetMarkerStyle(25);
    nRho_nConst_cut_efficiency2->SetMarkerColor(kGreen);
    nRho_nConst_cut_efficiency2->DrawCopy("HIST P SAME");

    nRho_nConst_cut_purity2->SetMarkerStyle(25);
    nRho_nConst_cut_purity2->SetMarkerColor(kMagenta);
    nRho_nConst_cut_purity2->DrawCopy("HIST P SAME");

    TLegend *legend2 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend2->AddEntry(nRho_nConst_cut_purity1, Form("Purity (N_{const} cut = %.1f)", first_nConstCut), "l");
    legend2->AddEntry(nRho_nConst_cut_efficiency1, Form("Efficiency (N_{const} cut = %.1f)", first_nConstCut), "l");
    legend2->AddEntry(nRho_nConst_cut_purity2, Form("Purity (N_{const} cut = %.1f)", second_nConstCut), "p");
    legend2->AddEntry(nRho_nConst_cut_efficiency2, Form("Efficiency (N_{const} cut = %.1f)", second_nConstCut), "p");
    legend2->Draw();

    TCanvas *c3 = new TCanvas("c3", "Purity and efficiency: N_{#rho} scan with p_{T} cuts", 2500, 2500);
    c3->Divide(1,1);

    c3->cd(1);
    nRho_pT_cut_efficiency1->SetTitle(("Purity and Efficiency for N_{#rho} scan with p_{T} cuts (" + switch_string + " jets)").c_str());
    nRho_pT_cut_efficiency1->GetXaxis()->SetTitle("N_{#rho} threshold (GeV/c)");
    nRho_pT_cut_efficiency1->GetYaxis()->SetTitle("Purity and Efficiency");
    nRho_pT_cut_efficiency1->SetLineColor(kMagenta);
    nRho_pT_cut_efficiency1->DrawCopy("HIST"); 

    nRho_pT_cut_purity1->SetLineColor(kBlue);
    nRho_pT_cut_purity1->DrawCopy("HIST SAME");

    nRho_pT_cut_efficiency2->SetMarkerStyle(25);
    nRho_pT_cut_efficiency2->SetMarkerColor(kRed);
    nRho_pT_cut_efficiency2->DrawCopy("HIST P SAME");

    nRho_pT_cut_purity2->SetMarkerStyle(25);
    nRho_pT_cut_purity2->SetMarkerColor(kSpring);
    nRho_pT_cut_purity2->DrawCopy("HIST P SAME");

    TLegend *legend3 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend3->AddEntry(nRho_pT_cut_purity1, Form("Purity (p_{T} cut = %.1f)", first_pTCut), "l");
    legend3->AddEntry(nRho_pT_cut_efficiency1, Form("Efficiency (p_{T} cut = %.1f)", first_pTCut), "l");
    legend3->AddEntry(nRho_pT_cut_purity2, Form("Purity (p_{T} cut = %.1f)", second_pTCut), "p");
    legend3->AddEntry(nRho_pT_cut_efficiency2, Form("Efficiency (p_{T} cut = %.1f)", second_pTCut), "p");
    legend3->Draw();

    TCanvas *c4 = new TCanvas("c4", "Purity and efficiency: N_{#rho} scan with N_{const} and p_{T} cuts", 2500, 2500);
    c4->Divide(1,1);

    c4->cd(1);
    nRho_nConst_pT_cut_efficiency1->SetTitle(("Purity and Efficiency for N_{#rho} scan with N_{const} and p_{T} cuts (" + switch_string + " jets)").c_str());
    nRho_nConst_pT_cut_efficiency1->GetXaxis()->SetTitle("N_{#rho} threshold (GeV/c)");
    nRho_nConst_pT_cut_efficiency1->GetYaxis()->SetTitle("Purity and Efficiency");
    nRho_nConst_pT_cut_efficiency1->SetLineColor(kRed);
    nRho_nConst_pT_cut_efficiency1->DrawCopy("HIST");

    nRho_nConst_pT_cut_purity1->SetLineColor(kBlue);
    nRho_nConst_pT_cut_purity1->DrawCopy("HIST SAME");

    nRho_nConst_pT_cut_efficiency2->SetMarkerStyle(25);
    nRho_nConst_pT_cut_efficiency2->SetMarkerColor(kGreen);
    nRho_nConst_pT_cut_efficiency2->DrawCopy("HIST P SAME");

    nRho_nConst_pT_cut_purity2->SetMarkerStyle(25);
    nRho_nConst_pT_cut_purity2->SetMarkerColor(kMagenta);
    nRho_nConst_pT_cut_purity2->DrawCopy("HIST P SAME");

    nRho_nConst_pT_cut_efficiency3->SetMarkerStyle(26);
    nRho_nConst_pT_cut_efficiency3->SetMarkerColor(kCyan);
    nRho_nConst_pT_cut_efficiency3->DrawCopy("HIST P SAME");

    nRho_nConst_pT_cut_purity3->SetMarkerStyle(26);
    nRho_nConst_pT_cut_purity3->SetMarkerColor(kBlack);
    nRho_nConst_pT_cut_purity3->DrawCopy("HIST P SAME");

    TLegend *legend4 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend4->AddEntry(nRho_nConst_pT_cut_purity1, Form("Purity (N_{const} cut = %.1f, p_{T} cut = %.1f)", first_nConstCut, first_pTCut), "l");
    legend4->AddEntry(nRho_nConst_pT_cut_efficiency1, Form("Efficiency (N_{const} cut = %.1f, p_{T} cut = %.1f)", first_nConstCut, first_pTCut), "l");
    legend4->AddEntry(nRho_nConst_pT_cut_purity2, Form("Purity (N_{const} cut = %.1f, p_{T} cut = %.1f)", first_nConstCut, second_pTCut), "p");
    legend4->AddEntry(nRho_nConst_pT_cut_efficiency2, Form("Efficiency (N_{const} cut = %.1f, p_{T} cut = %.1f)", first_nConstCut, second_pTCut), "p");
    legend4->AddEntry(nRho_nConst_pT_cut_purity3, Form("Purity (N_{const} cut = %.1f, p_{T} cut = %.1f)", second_nConstCut, second_pTCut), "p");
    legend4->AddEntry(nRho_nConst_pT_cut_efficiency3, Form("Efficiency (N_{const} cut = %.1f, p_{T} cut = %.1f)", second_nConstCut, second_pTCut), "p");
    legend4->Draw();

}
