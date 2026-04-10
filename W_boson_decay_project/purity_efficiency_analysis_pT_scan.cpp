//===========================================================================================================================================================================
// Classical purity and efficiency analysis:
// This macro intends to study the performance of varios combinations of cuts on the discriminatory variables through the analysis of classical purity and efficiency
// distributions. The scan is performed on the pT variable.
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

void purity_efficiency_analysis_pT_scan(std::string switch_string, const char* strange_file)
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

    TH1F* pT_cut_purity = new TH1F("pT_cut_purity", "Purity for pT scan; pT cut (GeV/c); Purity", 120, 0, 60);
    TH1F* pT_cut_efficiency = new TH1F("pT_cut_efficiency", "Efficiency for pT scan; pT cut (GeV/c); Efficiency", 120, 0, 60);


    TH1F* pT_nConst_cut_purity1 = new TH1F("pT_nConst_cut_purity1", "Purity for pT scan with nConst cut (first cut); pT cut (GeV/c); Purity", 120, 0, 60);
    TH1F* pT_nConst_cut_efficiency1 = new TH1F("pT_nConst_cut_efficiency1", "Efficiency for pT scan with nConst cut (first cut); pT cut (GeV/c); Efficiency", 120, 0, 60);

    TH1F* pT_nConst_cut_purity2 = new TH1F("pT_nConst_cut_purity2", "Purity for pT scan with nConst cut (second cut); pT cut (GeV/c); Purity", 120, 0, 60);
    TH1F* pT_nConst_cut_efficiency2 = new TH1F("pT_nConst_cut_efficiency2", "Efficiency for pT scan with nConst cut (second cut); pT cut (GeV/c); Efficiency", 120, 0, 60);


    TH1F* pT_nRho_cut_purity1 = new TH1F("pT_nRho_cut_purity1", "Purity for pT scan with nRho cut (first cut); pT cut (GeV/c); Purity", 120, 0, 60);
    TH1F* pT_nRho_cut_efficiency1 = new TH1F("pT_nRho_cut_efficiency1", "Efficiency for pT scan with nRho cut (first cut); pT cut (GeV/c); Efficiency", 120, 0, 60);

    TH1F* pT_nRho_cut_purity2 = new TH1F("pT_nRho_cut_purity2", "Purity for pT scan with nRho cut (second cut); pT cut (GeV/c); Purity", 120, 0, 60);
    TH1F* pT_nRho_cut_efficiency2 = new TH1F("pT_nRho_cut_efficiency2", "Efficiency for pT scan with nRho cut (second cut); pT cut (GeV/c); Efficiency", 120, 0, 60);


    TH1F* pT_nConst_nRho_cut_purity1 = new TH1F("pT_nConst_nRho_cut_purity1", "Purity for pT scan with nConst and nRho cut (first cuts); pT cut (GeV/c); Purity", 120, 0, 60);
    TH1F* pT_nConst_nRho_cut_efficiency1 = new TH1F("pT_nConst_nRho_cut_efficiency1", "Efficiency for pT scan with nConst and nRho cut (first cuts); pT cut (GeV/c); Efficiency", 120, 0, 60);

    TH1F* pT_nConst_nRho_cut_purity2 = new TH1F("pT_nConst_nRho_cut_purity2", "Purity for pT scan with nConst and nRho cut (second cuts); pT cut (GeV/c); Purity", 120, 0, 60);
    TH1F* pT_nConst_nRho_cut_efficiency2 = new TH1F("pT_nConst_nRho_cut_efficiency2", "Efficiency for pT scan with nConst and nRho cut (second cuts); pT cut (GeV/c); Efficiency", 120, 0, 60);

    TH1F* pT_nConst_nRho_cut_purity3 = new TH1F("pT_nConst_nRho_cut_purity3", "Purity for pT scan with nConst and nRho cut (third cuts); pT cut (GeV/c); Purity", 120, 0, 60);
    TH1F* pT_nConst_nRho_cut_efficiency3 = new TH1F("pT_nConst_nRho_cut_efficiency3", "Efficiency for pT scan with nConst and nRho cut (third cuts); pT cut (GeV/c); Efficiency", 120, 0, 60);

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
    
    Float_t pT_maxThreshold = 100;
    
    // Unidimensional pT scan (no cuts included)
    for (Float_t pT_thr = 0; pT_thr <= pT_maxThreshold; pT_thr += 0.1)
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

    // Unidimensional pT scan with nConst cuts
    for (Float_t pT_thr = 0; pT_thr <= pT_maxThreshold; pT_thr += 0.1)
    {
        Int_t TP1 = 0, FN1 = 0, TN1 = 0, FP1 = 0;
        Int_t TP2 = 0, FN2 = 0, TN2 = 0, FP2 = 0;

        for (Long64_t i = 0; i < signal_pT.size(); i++)
        {
            if (signal_pT[i] >= pT_thr && signal_nConst[i] >= first_nConstCut)
            {
                TP1++;
            }
            else
            {
                FN1++;
            }

            if (signal_pT[i] >= pT_thr && signal_nConst[i] >= second_nConstCut)
            {
                TP2++;
            }
            else
            {
                FN2++;
            }

        }

        for (Long64_t i = 0; i < background_pT.size(); i++)
        {
            if (background_pT[i] >= pT_thr && background_nConst[i] >= first_nConstCut)
            {
                FP1++;
            }
            else
            {
                TN1++;
            }

            if (background_pT[i] >= pT_thr && background_nConst[i] >= second_nConstCut)
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

        pT_nConst_cut_purity1->Fill(pT_thr, purity1);
        pT_nConst_cut_efficiency1->Fill(pT_thr, efficiency1);

        pT_nConst_cut_purity2->Fill(pT_thr, purity2);
        pT_nConst_cut_efficiency2->Fill(pT_thr, efficiency2);

    }

    Float_t first_nRhoCut = 0, second_nRhoCut = 0;

    if (switch_string == "strange") 
    {
        first_nRhoCut = 0.5;
        second_nRhoCut = 1;
    } 
    else if (switch_string == "charm")
    {
        first_nRhoCut = 1;
        second_nRhoCut = 2;
    }


    // Unidimensional pT scan with nConst cuts
    for (Float_t pT_thr = 0; pT_thr <= pT_maxThreshold; pT_thr += 0.1)
    {
        Int_t TP1 = 0, FN1 = 0, TN1 = 0, FP1 = 0;
        Int_t TP2 = 0, FN2 = 0, TN2 = 0, FP2 = 0;

        for (Long64_t i = 0; i < signal_pT.size(); i++)
        {
            if (signal_pT[i] >= pT_thr && signal_nRho[i] >= first_nRhoCut)
            {
                TP1++;
            }
            else
            {
                FN1++;
            }

            if (signal_pT[i] >= pT_thr && signal_nRho[i] >= second_nRhoCut)
            {
                TP2++;
            }
            else
            {
                FN2++;
            }

        }

        for (Long64_t i = 0; i < background_pT.size(); i++)
        {
            if (background_pT[i] >= pT_thr && background_nRho[i] >= first_nRhoCut)
            {
                FP1++;
            }
            else
            {
                TN1++;
            }

            if (background_pT[i] >= pT_thr && background_nRho[i] >= second_nRhoCut)
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

        pT_nRho_cut_purity1->Fill(pT_thr, purity1);
        pT_nRho_cut_efficiency1->Fill(pT_thr, efficiency1);

        pT_nRho_cut_purity2->Fill(pT_thr, purity2);
        pT_nRho_cut_efficiency2->Fill(pT_thr, efficiency2);
    }

    // Unidimensional pT scan with nConst and nRho cuts
    for (Float_t pT_thr = 0; pT_thr <= pT_maxThreshold; pT_thr += 0.1)
    {
        Int_t TP1 = 0, FN1 = 0, TN1 = 0, FP1 = 0;
        Int_t TP2 = 0, FN2 = 0, TN2 = 0, FP2 = 0;
        Int_t TP3 = 0, FN3 = 0, TN3 = 0, FP3 = 0;

        for (Long64_t i = 0; i < signal_pT.size(); i++)
        {
            if (signal_pT[i] >= pT_thr && signal_nConst[i] >= first_nConstCut && signal_nRho[i] >= first_nRhoCut)
            {
                TP1++;
            }
            else
            {
                FN1++;
            }

            if (signal_pT[i] >= pT_thr && signal_nConst[i] >= first_nConstCut && signal_nRho[i] >= second_nRhoCut)
            {
                TP2++;
            }
            else
            {
                FN2++;
            }

            if (signal_pT[i] >= pT_thr && signal_nConst[i] >= second_nConstCut && signal_nRho[i] >= second_nRhoCut)
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
            if (background_pT[i] >= pT_thr && background_nConst[i] >= first_nConstCut && background_nRho[i] >= first_nRhoCut)
            {
                FP1++;
            }
            else
            {
                TN1++;
            }

            if (background_pT[i] >= pT_thr && background_nConst[i] >= first_nConstCut && background_nRho[i] >= second_nRhoCut)
            {
                FP2++;
            }
            else
            {
                TN2++;
            }

            if (background_pT[i] >= pT_thr && background_nConst[i] >= second_nConstCut && background_nRho[i] >= second_nRhoCut)
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

        pT_nConst_nRho_cut_purity1->Fill(pT_thr, purity1);
        pT_nConst_nRho_cut_efficiency1->Fill(pT_thr, efficiency1);

        pT_nConst_nRho_cut_purity2->Fill(pT_thr, purity2);
        pT_nConst_nRho_cut_efficiency2->Fill(pT_thr, efficiency2);

        pT_nConst_nRho_cut_purity3->Fill(pT_thr, purity3);
        pT_nConst_nRho_cut_efficiency3->Fill(pT_thr, efficiency3);

    }

    TCanvas *c1 = new TCanvas("c1", "Purity and efficiency: p_{T} scan", 2500, 2500);
    c1->Divide(1,1);

    c1->cd(1);
    pT_cut_efficiency->SetTitle(("Purity and Efficiency for p_{T} scan (" + switch_string + " jets)").c_str());
    pT_cut_efficiency->GetXaxis()->SetTitle("p_{T} threshold (GeV/c)");
    pT_cut_efficiency->GetYaxis()->SetTitle("Purity and Efficiency");
    pT_cut_efficiency->SetLineColor(kRed);
    pT_cut_efficiency->DrawCopy("HIST");

    pT_cut_purity->SetLineColor(kBlue);
    pT_cut_purity->DrawCopy("HIST SAME");

    TLegend *legend1 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend1->AddEntry(pT_cut_purity, "Purity", "l");
    legend1->AddEntry(pT_cut_efficiency, "Efficiency", "l");
    legend1->Draw();

    TCanvas *c2 = new TCanvas("c2", "Purity and efficiency: p_{T} scan with nConst cuts", 2500, 2500);
    c2->Divide(1,1);

    c2->cd(2);
    pT_nConst_cut_efficiency1->SetTitle(("Purity and Efficiency for p_{T} scan with N_{const} cuts (" + switch_string + " jets)").c_str());
    pT_nConst_cut_efficiency1->GetXaxis()->SetTitle("p_{T} threshold (GeV/c)");
    pT_nConst_cut_efficiency1->GetYaxis()->SetTitle("Purity and Efficiency");
    pT_nConst_cut_efficiency1->SetLineColor(kRed);
    pT_nConst_cut_efficiency1->DrawCopy("HIST");

    pT_nConst_cut_purity1->SetLineColor(kBlue);
    pT_nConst_cut_purity1->DrawCopy("HIST SAME");

    pT_nConst_cut_efficiency2->SetMarkerStyle(25);
    pT_nConst_cut_efficiency2->SetMarkerColor(kGreen);
    pT_nConst_cut_efficiency2->DrawCopy("HIST P SAME");

    pT_nConst_cut_purity2->SetMarkerStyle(25);
    pT_nConst_cut_purity2->SetMarkerColor(kMagenta);
    pT_nConst_cut_purity2->DrawCopy("HIST P SAME");

    TLegend *legend2 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend2->AddEntry(pT_nConst_cut_purity1, Form("Purity (N_{const} cut = %.1f)", first_nConstCut), "l");
    legend2->AddEntry(pT_nConst_cut_efficiency1, Form("Efficiency (N_{const} cut = %.1f)", first_nConstCut), "l");
    legend2->AddEntry(pT_nConst_cut_purity2, Form("Purity (N_{const} cut = %.1f)", second_nConstCut), "p");
    legend2->AddEntry(pT_nConst_cut_efficiency2, Form("Efficiency (N_{const} cut = %.1f)", second_nConstCut), "p");
    legend2->Draw();

    TCanvas *c3 = new TCanvas("c3", "Purity and efficiency: p_{T} scan with N_{#rho} cuts", 2500, 2500);
    c3->Divide(1,1);

    c3->cd(1);
    pT_nRho_cut_efficiency1->SetTitle(("Purity and Efficiency for p_{T} scan with N_{#rho} cuts (" + switch_string + " jets)").c_str());
    pT_nRho_cut_efficiency1->GetXaxis()->SetTitle("p_{T} threshold (GeV/c)");
    pT_nRho_cut_efficiency1->GetYaxis()->SetTitle("Purity and Efficiency");
    pT_nRho_cut_efficiency1->SetLineColor(kMagenta);
    pT_nRho_cut_efficiency1->DrawCopy("HIST"); 

    pT_nRho_cut_purity1->SetLineColor(kBlue);
    pT_nRho_cut_purity1->DrawCopy("HIST SAME");

    pT_nRho_cut_efficiency2->SetMarkerStyle(25);
    pT_nRho_cut_efficiency2->SetMarkerColor(kRed);
    pT_nRho_cut_efficiency2->DrawCopy("HIST P SAME");

    pT_nRho_cut_purity2->SetMarkerStyle(25);
    pT_nRho_cut_purity2->SetMarkerColor(kSpring);
    pT_nRho_cut_purity2->DrawCopy("HIST P SAME");

    TLegend *legend3 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend3->AddEntry(pT_nRho_cut_purity1, Form("Purity (N_{#rho} cut = %.1f)", first_nRhoCut), "l");
    legend3->AddEntry(pT_nRho_cut_efficiency1, Form("Efficiency (N_{#rho} cut = %.1f)", first_nRhoCut), "l");
    legend3->AddEntry(pT_nRho_cut_purity2, Form("Purity (N_{#rho} cut = %.1f)", second_nRhoCut), "p");
    legend3->AddEntry(pT_nRho_cut_efficiency2, Form("Efficiency (N_{#rho} cut = %.1f)", second_nRhoCut), "p");
    legend3->Draw();

    TCanvas *c4 = new TCanvas("c4", "Purity and efficiency: pT scan with N_{const} and N_{#rho} cuts", 2500, 2500);
    c4->Divide(1,1);

    c4->cd(1);
    pT_nConst_nRho_cut_efficiency1->SetTitle(("Purity and Efficiency for p_{T} scan with N_{const} and N_{#rho} cuts (" + switch_string + " jets)").c_str());
    pT_nConst_nRho_cut_efficiency1->GetXaxis()->SetTitle("p_{T} threshold (GeV/c)");
    pT_nConst_nRho_cut_efficiency1->GetYaxis()->SetTitle("Purity and Efficiency");
    pT_nConst_nRho_cut_efficiency1->SetLineColor(kRed);
    pT_nConst_nRho_cut_efficiency1->DrawCopy("HIST");

    pT_nConst_nRho_cut_purity1->SetLineColor(kBlue);
    pT_nConst_nRho_cut_purity1->DrawCopy("HIST SAME");

    pT_nConst_nRho_cut_efficiency2->SetMarkerStyle(25);
    pT_nConst_nRho_cut_efficiency2->SetMarkerColor(kGreen);
    pT_nConst_nRho_cut_efficiency2->DrawCopy("HIST P SAME");

    pT_nConst_nRho_cut_purity2->SetMarkerStyle(25);
    pT_nConst_nRho_cut_purity2->SetMarkerColor(kMagenta);
    pT_nConst_nRho_cut_purity2->DrawCopy("HIST P SAME");

    pT_nConst_nRho_cut_efficiency3->SetMarkerStyle(26);
    pT_nConst_nRho_cut_efficiency3->SetMarkerColor(kCyan);
    pT_nConst_nRho_cut_efficiency3->DrawCopy("HIST P SAME");

    pT_nConst_nRho_cut_purity3->SetMarkerStyle(26);
    pT_nConst_nRho_cut_purity3->SetMarkerColor(kBlack);
    pT_nConst_nRho_cut_purity3->DrawCopy("HIST P SAME");

    TLegend *legend4 = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend4->AddEntry(pT_nConst_nRho_cut_purity1, Form("Purity (N_{const} cut = %.1f, N_{#rho} cut = %.1f)", first_nConstCut, first_nRhoCut), "l");
    legend4->AddEntry(pT_nConst_nRho_cut_efficiency1, Form("Efficiency (N_{const} cut = %.1f, N_{#rho} cut = %.1f)", first_nConstCut, first_nRhoCut), "l");
    legend4->AddEntry(pT_nConst_nRho_cut_purity2, Form("Purity (N_{const} cut = %.1f, N_{#rho} cut = %.1f)", first_nConstCut, second_nRhoCut), "p");
    legend4->AddEntry(pT_nConst_nRho_cut_efficiency2, Form("Efficiency (N_{const} cut = %.1f, N_{#rho} cut = %.1f)", first_nConstCut, second_nRhoCut), "p");
    legend4->AddEntry(pT_nConst_nRho_cut_purity3, Form("Purity (N_{const} cut = %.1f, N_{#rho} cut = %.1f)", second_nConstCut, second_nRhoCut), "p");
    legend4->AddEntry(pT_nConst_nRho_cut_efficiency3, Form("Efficiency (N_{const} cut = %.1f, N_{#rho} cut = %.1f)", second_nConstCut, second_nRhoCut), "p");
    legend4->Draw();

}
