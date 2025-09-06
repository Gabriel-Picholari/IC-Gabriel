//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
/*

  This code determines the cuts for multiple potential discriminating variables through cumulative histograms. The differences from the previous version include the addition of a new discriminating 
variable and the creation of cumulative histograms for each of them. This macro is intended to be used in conjunction with wdecayTTree2.cpp.

*/
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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

Float_t cut_binCalculation(TH1F* signalHist, TH1F* backgroundHist)
{
    Float_t minRatio = std::numeric_limits<Float_t>::max();
    Int_t minRatioBin = 0;

    Int_t nBins /*Same for both cumulative histograms*/ = signalHist->GetNbinsX();

    for (Int_t i = 0; i <= nBins; i++)
    {
        Double_t signalValue = signalHist->GetBinContent(i);
        Double_t backgroundValue = backgroundHist->GetBinContent(i);
        Double_t ratio = signalValue / backgroundValue;

        if (ratio < minRatio)
        {
            minRatio = ratio;
            minRatioBin = i;
        }
    }

    Float_t chosenBinValue = signalHist->GetBinCenter(minRatioBin);

    return chosenBinValue;
}


void read_wdecayTTree2(const char* fileName)
{

    gSystem->Load("libEG");
    gSystem->Load("libEGPythia8");

    //---------------------------------------------------------------------------------------------------------
    //Inicializacao das variaveis
    //---------------------------------------------------------------------------------------------------------

    Float_t fpPt, fpEta, fpPhi, fpE, fpPx, fpPy, fpPz, fpMass = 0;
    Float_t jetPt, jetEta, jetPhi, jetE, jetPx, jetPy, jetPz, jetMass, jetNConst, pT_LeadConst = 0;
    Float_t s_pT, s_Eta, s_Phi, c_pT, c_Eta, c_Phi = 0;
    Float_t sbar_pT, sbar_Eta, sbar_Phi = 0;
    Float_t cbar_pT, cbar_Eta, cbar_Phi = 0;
    Float_t maxRho, nVert;
    Float_t match_R = 0.5;

    // Criacao dos TLorentzVector

    TLorentzVector vec_s(0,0,0,0);
    TLorentzVector vec_c(0,0,0,0);
    TLorentzVector vec_sbar(0,0,0,0);
    TLorentzVector vec_cbar(0,0,0,0);

    //---------------------------------------------------------------------------------------------------------
    // Inicializacao dos histogramas
    //---------------------------------------------------------------------------------------------------------

    TH1F *cJet_pT = new TH1F("cJet_pT", "Jet p_{T} from c quark [GeV/c] (signal)", 100, 0, 100);
    TH1F *cbarJet_pT = new TH1F("cbarJet_pT", "p_{T} of jet from anti-c quark [GeV/c]", 100, 0, 100);

    TH1F *sJet_pT = new TH1F("sJet_pT", "p_{T} of jet from s quark [GeV/c]", 100, 0, 100);
    TH1F *sbarJet_pT = new TH1F("sbarJet_pT", "Jet p_{T} from anti-s quark [GeV/c] (signal)", 100, 0, 100);

    TH1F *invariantMass_cbar_s = new TH1F("invariantMass_cbar_s", "Invariant mass spectrum of jets from anti-c and s quarks [GeV/c^{2}]", 100, 0, 100);
    TH1F *invariantMass_c_sbar = new TH1F("invariantMass_c_sbar", "Invariant mass spectrum of jets from c and anti-s quarks [GeV/c^{2}]", 100, 0, 100);
    TH1F *invariantMass = new TH1F("invariantMass", "Invariant mass spectrum [GeV/c^{2}]", 98, 2, 100);



    //---------------------------------------------------------------------------------------------------------
    // Inicializacoes e configuracoes do FastJet:
    //---------------------------------------------------------------------------------------------------------

    Float_t jetR = 0.4;

    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, jetR);

    std::vector<fastjet::PseudoJet> particles_fastjet;
    std::vector<fastjet::PseudoJet> jets;

    //---------------------------------------------------------------------------------------------------------
    // Inicializacao do arquivo.root e das TTrees
    //---------------------------------------------------------------------------------------------------------

    TFile *file = TFile::Open(fileName, "READ");
    TTree *ttree = dynamic_cast<TTree *>(file->Get("W decay TTree 2"));

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

        vec_s.SetPtEtaPhiM(10,0,0,(3.141592));
        vec_c.SetPtEtaPhiM(10,0,0,(3.141592));
        vec_sbar.SetPtEtaPhiM(10,0,0,(3.141592));
        vec_cbar.SetPtEtaPhiM(10,0,0,(3.141592));


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
            if (jetPt < 5) continue;

            jetEta = jet.eta();
            jetPhi = jet.phi();
            jetMass = jet.m(); // Invariant mass
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
            // Comparacao dos jatos com os quarks -> Match entre ambos
            //---------------------------------------------------------------------------------------------------------

            for (Int_t nj2 = 0; nj2 < quarks->GetEntries(); nj2++)
            {

                MyQuark *quarkPart = static_cast<MyQuark *>(quarks->At(nj2));
                Int_t nj2Pdg = quarkPart->qPdg;

                
                if (nj2Pdg == 3)
                {
                    MyQuark *sq = static_cast<MyQuark *>(quarks->At(nj2));

                    s_pT = sq->qPdg;
                    s_Eta = sq->qEta;
                    s_Phi = sq->qPhi;

                    Float_t distancia_s = TMath::Sqrt(TMath::Power(jetPhi - s_Phi, 2) + TMath::Power(jetEta - s_Eta, 2));

                    if (distancia_s <= match_R && jetPt)
                    {
                        sJet_pT->Fill(jetPt);
                        //cout << "sJet: " << jetPt << ", " << jetMass << endl;
                        vec_s = TLorentzVector(jetPx, jetPy, jetPz, jetE);
                    }
                }

                //---------------------------------------------------------------------------------------------------------
                // QUARK ANTI-S
                //---------------------------------------------------------------------------------------------------------

                if (nj2Pdg == -3)
                {
                    MyQuark *sbarq = static_cast<MyQuark *>(quarks->At(nj2));

                    sbar_pT = sbarq->qpT;
                    sbar_Eta = sbarq->qEta;
                    sbar_Phi = sbarq->qPhi;

                    Float_t distancia_sbar = TMath::Sqrt(TMath::Power(jetPhi - sbar_Phi, 2) + TMath::Power(jetEta - sbar_Eta, 2));

                    if (distancia_sbar <= match_R)
                    {
                        sbarJet_pT->Fill(jetPt);
                        vec_sbar = TLorentzVector(jetPx, jetPy, jetPz, jetE);
                    }
                }
                //---------------------------------------------------------------------------------------------------------
                // QUARK C
                //---------------------------------------------------------------------------------------------------------

                if (nj2Pdg == 4)
                {
                    MyQuark *cq = static_cast<MyQuark *>(quarks->At(nj2));

                    c_pT = cq->qpT;
                    c_Eta = cq->qEta;
                    c_Phi = cq->qPhi;

                    Float_t distancia_c = TMath::Sqrt(TMath::Power(jetPhi - c_Phi, 2) + TMath::Power(jetEta - c_Eta, 2));

                    if (distancia_c <= match_R)
                    {
                        cJet_pT->Fill(jetPt);
                        vec_c = TLorentzVector(jetPx, jetPy, jetPz, jetE);
                    }
                }

                //---------------------------------------------------------------------------------------------------------
                
                if (nj2Pdg == -4)
                {
                    MyQuark *cbarq = static_cast<MyQuark *>(quarks->At(nj2));

                    cbar_pT = cbarq->qpT;
                    cbar_Eta = cbarq->qEta;
                    cbar_Phi = cbarq->qPhi;

                    Float_t distancia_cbar = TMath::Sqrt(TMath::Power(jetPhi - cbar_Phi, 2) + TMath::Power(jetEta - cbar_Eta, 2));

                    if (distancia_cbar <= match_R)
                    {
                        cbarJet_pT->Fill(jetPt);
                        //cout << "cbarJet " << jetPt << ", " << jetMass << endl;
                        vec_cbar = TLorentzVector(jetPx, jetPy, jetPz, jetE);
                    }
                }
                

            } // End of particles and jets particular matches

        } // End of individual jet creation

        const Float_t tolerance = 1e-5;
        TLorentzVector totalVector;

        if ( fabs(vec_cbar.M() - 3.141592) > tolerance && fabs(vec_s.M() - 3.141592) > tolerance ) 
        {
            TLorentzVector vec_cbar_s_event = vec_cbar + vec_s;
            invariantMass_cbar_s->Fill( vec_cbar_s_event.M() );
            totalVector = totalVector + vec_cbar_s_event;
        }
        if ( fabs(vec_sbar.M() - 3.141592) > tolerance && fabs(vec_c.M() - 3.141592) > tolerance )
        {
            TLorentzVector vec_c_sbar_event = vec_c + vec_sbar;
            invariantMass_c_sbar->Fill( vec_c_sbar_event.M() );
            totalVector = totalVector + vec_c_sbar_event;
        }

        invariantMass->Fill(totalVector.M());

        particles_fastjet.clear();
        jets.clear();
        jets_array->Clear();
        quarks->Clear();

    } // End of event loop equivalent

    //---------------------------------------------------------------------------------------------------------
    // Canvas creation and plotting
    //---------------------------------------------------------------------------------------------------------

    TCanvas *c1 = new TCanvas("c1", "Quark's pT histograms", 2500, 2500);
    c1->Divide(2, 2);

    c1->cd(1);
    sJet_pT->SetTitle("s Quark Jet's Transverse Momentum");
    sJet_pT->GetXaxis()->SetTitle("Transverse Momentum [GeV/c]");
    sJet_pT->GetYaxis()->SetTitle("Frequency");
    sJet_pT->Draw();

    c1->cd(2);
    sbarJet_pT->SetTitle("Anti-s Quark Jet's Transverse Momentum");
    sbarJet_pT->GetXaxis()->SetTitle("Transverse Momentum [GeV/c]");
    sbarJet_pT->GetYaxis()->SetTitle("Frequency");
    sbarJet_pT->Draw();

    c1->cd(3);
    cJet_pT->SetTitle("c Quark Jet's Transverse Momentum");
    cJet_pT->GetXaxis()->SetTitle("Transverse Momentum [GeV/c]");
    cJet_pT->GetYaxis()->SetTitle("Frequency");
    cJet_pT->Draw();

    c1->cd(4);
    cbarJet_pT->SetTitle("Anti-c Quark Jet's Transverse Momentum");
    cbarJet_pT->GetXaxis()->SetTitle("Transverse Momentum [GeV/c]");
    cbarJet_pT->GetYaxis()->SetTitle("Frequency");
    cbarJet_pT->Draw();

    //---------------------------------------------------------------------------------------------------------

    TCanvas *c2 = new TCanvas("c2", "Invariant mass distributions", 2500, 2500);
    c2->Divide(1, 3);

    c2->cd(1);
    invariantMass_cbar_s->SetTitle("Jet's invariant mass spectrum - anti-c and s");
    invariantMass_cbar_s->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    invariantMass_cbar_s->GetYaxis()->SetTitle("Frequency");
    invariantMass_cbar_s->Draw();

    c2->cd(2);
    invariantMass_c_sbar->SetTitle("Jet's invariant mass spectrum - c and anti-s");
    invariantMass_c_sbar->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    invariantMass_c_sbar->GetYaxis()->SetTitle("Frequency");
    invariantMass_c_sbar->Draw();

    c2->cd(3);
    invariantMass->SetTitle("Jet's invariant mass spectrum");
    invariantMass->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    invariantMass->GetYaxis()->SetTitle("Frequency");
    invariantMass->Draw();


    //---------------------------------------------------------------------------------------------------------
    // Saving histograms
    //---------------------------------------------------------------------------------------------------------

    /* 
    TFile *outPlots = new TFile("Histogramas_read_wdecayTTree2.root", "RECREATE");
    
    for (int t = 1; t <= 40; t++) 
    {
        TH1F *histo = (TH1F*) gDirectory->Get(Form("h%d", t));
        if (histo) {
            histo->Write();
        } 
        else 
        {
            std::cout << "Histograma h" << t << " não encontrado ou nulo!" << std::endl;
        }
    }

    outPlots->Close();
    */

    file->Close();
}