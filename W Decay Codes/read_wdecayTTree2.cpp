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

    TH1F *cJet_pT = new TH1F("h1", "Jet p_{T} from c quark [GeV/c] (signal)", 100, 0, 100);
    TH1F *cJet_nConst = new TH1F("h2", "Number of constituents of jet from c quark (signal)", 100, 0, 100);
    TH1F *cJet_leadConst_pT = new TH1F("h3", "p_{T} of lead constituent of jet from c quark [GeV/c] (signal)", 100, 0, 100);
    TH1F *cJet_eta = new TH1F("h4", "Jet #eta from c quark (signal)", 100, -2, 2);
    TH1F *cJet_phi = new TH1F("h5", "Jet #phi from c quark (signal)", 100, -3.14, 3.14);
    TH1F *cJet_energy = new TH1F("h6", "Jet Energy from c quark [GeV] (signal)", 100, 0, 100);
    TH1F *cJet_sigmaKT = new TH1F("h7", "Jet sigma_{KT} from c quark (signal)", 100, 0, 1);
    TH1F *cJet_MaxRho = new TH1F("h8", "Jet Max Rho from c quark (signal)", 100, 0, 100);
    TH1F *cJet_vertex = new TH1F("h9", "Jet Vertex from c quark (signal)", 100, 0, 100);
    TH1F *cJet_angularity = new TH1F("h10", "Jet Angularity from c quark (signal)", 100, 0, 1);

    TH1F *cJet_pT_bkg = new TH1F("h11", "Jet p_{T} from c quark [GeV/c] (background)", 100, 0, 100);
    TH1F *cJet_nConst_bkg = new TH1F("h12", "Number of constituents of jet from c quark (background)", 100, 0, 100);
    TH1F *cJet_leadConst_pT_bkg = new TH1F("h13", "p_{T} of lead constituent of jet from c quark [GeV/c] (background)", 100, 0, 100);
    TH1F *cJet_eta_bkg = new TH1F("h14", "Jet #eta from c quark (background)", 100, -2, 2);
    TH1F *cJet_phi_bkg = new TH1F("h15", "Jet #phi from c quark (background)", 100, -3.14, 3.14);
    TH1F *cJet_energy_bkg = new TH1F("h16", "Jet Energy from c quark [GeV] (background)", 100, 0, 200);
    TH1F *cJet_sigmaKT_bkg = new TH1F("h17", "Jet sigma_{KT} from c quark (background)", 100, 0, 1);
    TH1F *cJet_MaxRho_bkg = new TH1F("h18", "Jet Max Rho from c quark (background)", 100, 0, 100);
    TH1F *cJet_vertex_bkg = new TH1F("h19", "Jet Vertex from c quark (background)", 100, 0, 100);
    TH1F *cJet_angularity_bkg = new TH1F("h20", "Jet Angularity from c quark (background)", 100, 0, 1);

    TH1F *sbarJet_pT = new TH1F("h21", "Jet p_{T} from anti-s quark [GeV/c] (signal)", 100, 0, 100);
    TH1F *sbarJet_nConst = new TH1F("h22", "Number of constituents of jet from anti-s quark (signal)", 100, 0, 100);
    TH1F *sbarJet_leadConst_pT = new TH1F("h23", "p_{T} of lead constituent of jet from anti-s quark [GeV/c] (signal)", 100, 0, 100);
    TH1F *sbarJet_eta = new TH1F("h24", "Jet #eta from anti-s quark (signal)", 100, -2, 2);
    TH1F *sbarJet_phi = new TH1F("h25", "Jet #phi from anti-s quark (signal)", 100, -3.14, 3.14);
    TH1F *sbarJet_energy = new TH1F("h26", "Jet Energy from anti-s quark [GeV] (signal)", 100, 0, 200);
    TH1F *sbarJet_sigmaKT = new TH1F("h27", "Jet sigma_{KT} from anti-s quark (signal)", 100, 0, 1);
    TH1F *sbarJet_MaxRho = new TH1F("h28", "Jet Max Rho from anti-s quark (signal)", 100, 0, 100);
    TH1F *sbarJet_vertex = new TH1F("h29", "Jet Vertex from anti-s quark (signal)", 100, 0, 100);
    TH1F *sbarJet_angularity = new TH1F("h30", "Jet Angularity from anti-s quark (signal)", 100, 0, 1);

    TH1F *sbarJet_pT_bkg = new TH1F("h31", "Jet p_{T} from anti-s quark [GeV/c] (background)", 100, 0, 100);
    TH1F *sbarJet_nConst_bkg = new TH1F("h32", "Number of constituents of jet from anti-s quark (background)", 100, 0, 100);
    TH1F *sbarJet_leadConst_pT_bkg = new TH1F("h33", "p_{T} of lead constituent of jet from anti-s quark [GeV/c] (background)", 100, 0, 100);
    TH1F *sbarJet_eta_bkg = new TH1F("h34", "Jet #eta from anti-s quark (background)", 100, -2, 2);
    TH1F *sbarJet_phi_bkg = new TH1F("h35", "Jet #phi from anti-s quark (background)", 100, -3.14, 3.14);
    TH1F *sbarJet_energy_bkg = new TH1F("h36", "Jet Energy from anti-s quark [GeV] (background)", 100, 0, 200);
    TH1F *sbarJet_sigmaKT_bkg = new TH1F("h37", "Jet sigma_{KT} from anti-s quark (background)", 100, 0, 1);
    TH1F *sbarJet_MaxRho_bkg = new TH1F("h38", "Jet Max Rho from anti-s quark (background)", 100, 0, 100);
    TH1F *sbarJet_vertex_bkg = new TH1F("h39", "Jet Vertex from anti-s quark (background)", 100, 0, 100);
    TH1F *sbarJet_angularity_bkg = new TH1F("h40", "Jet Angularity from anti-s quark (background)", 100, 0, 1);

    TH1F *sJet_pT = new TH1F("h41", "p_{T} of jet from s quark [GeV/c]", 100, 0, 100);
    TH1F *cbarJet_pT = new TH1F("h42", "p_{T} of jet from anti-c quark [GeV/c]", 100, 0, 100);

    TH1F *invariantMass_cbar_s = new TH1F("h43", "Invariant mass spectrum of jets from anti-c and s quarks [GeV/c^{2}]", 100, 0, 100);
    TH1F *invariantMass_c_sbar = new TH1F("h44", "Invariant mass spectrum of jets from c and anti-s quarks [GeV/c^{2}]", 100, 0, 100);


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

            Float_t angAve, sigmaKT = 0;

            for (Int_t i = 0; i < jetNConst; ++i) 
            {
                Double_t pt_constituentes = jet.constituents()[i].pt();
                Double_t eta_constituentes = jet.constituents()[i].eta();
                Double_t phi_constituentes = jet.constituents()[i].phi();

                Double_t angPart = TMath::Sqrt(TMath::Power(TMath::Abs(phi_constituentes) - TMath::Abs(jetPhi), 2) + TMath::Power(eta_constituentes - jetEta, 2));

                angAve = angAve + angPart;
                sigmaKT = sigmaKT + (TMath::Power(pt_constituentes - (jetPt / jetNConst), 2));
            }

            angAve = angAve / jetNConst;
            sigmaKT = sigmaKT / jetNConst;

            //---------------------------------------------------------------------------------------------------------
            // Comparacao dos jatos com os quarks -> Match entre ambos
            //---------------------------------------------------------------------------------------------------------

            for (Int_t nj2 = 0; nj2 < quarks->GetEntries(); nj2++)
            {

                MyQuark *quarkPdg = static_cast<MyQuark *>(quarks->At(nj2));
                Int_t nj2Pdg = quarkPdg->qPdg;

                
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
                        sbarJet_eta->Fill(jetEta);
                        sbarJet_phi->Fill(jetPhi);
                        sbarJet_energy->Fill(jetE);
                        sbarJet_vertex->Fill(nVert);
                        sbarJet_MaxRho->Fill(maxRho);
                        sbarJet_sigmaKT->Fill(sigmaKT);
                        sbarJet_nConst->Fill(jetNConst);
                        sbarJet_angularity->Fill(angAve);
                        sbarJet_leadConst_pT->Fill(pT_LeadConst);

                        vec_sbar = TLorentzVector(jetPx, jetPy, jetPz, jetE);
                    }
                    else
                    {
                        sbarJet_pT_bkg->Fill(jetPt);
                        sbarJet_eta_bkg->Fill(jetEta);
                        sbarJet_phi_bkg->Fill(jetPhi);
                        sbarJet_energy_bkg->Fill(jetE);
                        sbarJet_vertex_bkg->Fill(nVert);
                        sbarJet_MaxRho_bkg->Fill(maxRho);
                        sbarJet_sigmaKT_bkg->Fill(sigmaKT);
                        sbarJet_nConst_bkg->Fill(jetNConst);
                        sbarJet_angularity_bkg->Fill(angAve);
                        sbarJet_leadConst_pT_bkg->Fill(pT_LeadConst);
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
                        cJet_eta->Fill(jetEta);
                        cJet_phi->Fill(jetPhi);
                        cJet_energy->Fill(jetE);
                        cJet_vertex->Fill(nVert);
                        cJet_MaxRho->Fill(maxRho);
                        cJet_sigmaKT->Fill(sigmaKT);
                        cJet_nConst->Fill(jetNConst);
                        cJet_angularity->Fill(angAve);
                        cJet_leadConst_pT->Fill(pT_LeadConst);

                        vec_c = TLorentzVector(jetPx, jetPy, jetPz, jetE);
                    }
                    else
                    {
                        cJet_pT_bkg->Fill(jetPt);
                        cJet_eta_bkg->Fill(jetEta);
                        cJet_phi_bkg->Fill(jetPhi);
                        cJet_energy_bkg->Fill(jetE);
                        cJet_vertex_bkg->Fill(nVert);
                        cJet_MaxRho_bkg->Fill(maxRho);
                        cJet_sigmaKT_bkg->Fill(sigmaKT);
                        cJet_nConst_bkg->Fill(jetNConst);
                        cJet_angularity_bkg->Fill(angAve);
                        cJet_leadConst_pT_bkg->Fill(pT_LeadConst);
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

        if ( fabs(vec_cbar.M() - 3.141592) > tolerance && fabs(vec_s.M() - 3.141592) > tolerance ) 
        {
            TLorentzVector vec_cbar_s_event = vec_cbar + vec_s;
            invariantMass_cbar_s->Fill( vec_cbar_s_event.M() ); 
        }
        if ( fabs(vec_sbar.M() - 3.141592) > tolerance && fabs(vec_c.M() - 3.141592) > tolerance )
        {
            TLorentzVector vec_c_sbar_event = vec_c + vec_sbar;
            invariantMass_c_sbar->Fill( vec_c_sbar_event.M() );
        }

        particles_fastjet.clear();
        jets.clear();
        jets_array->Clear();
        quarks->Clear();

    } // End of event loop equivalent

    //---------------------------------------------------------------------------------------------------------
    // Creation of the cumulative histograms regarding the classical analysis target object
    //---------------------------------------------------------------------------------------------------------

    TH1F *cJet_pT_cumulative = (TH1F *)cJet_pT->GetCumulative();
    cJet_pT_cumulative->Scale(1.0/cJet_pT_cumulative->GetMaximum());

    TH1F *cJet_nConst_cumulative = (TH1F *)cJet_nConst->GetCumulative();
    cJet_nConst_cumulative->Scale(1.0/cJet_nConst_cumulative->GetMaximum());

    TH1F *cJet_leadConst_pT_cumulative = (TH1F *)cJet_leadConst_pT->GetCumulative();
    cJet_leadConst_pT_cumulative->Scale(1.0/cJet_leadConst_pT_cumulative->GetMaximum());

    TH1F *cJet_eta_cumulative = (TH1F *)cJet_eta->GetCumulative();
    cJet_eta_cumulative->Scale(1.0/cJet_eta_cumulative->GetMaximum());

    TH1F *cJet_phi_cumulative = (TH1F *)cJet_phi->GetCumulative();
    cJet_phi_cumulative->Scale(1.0/cJet_phi_cumulative->GetMaximum());

    TH1F *cJet_energy_cumulative = (TH1F *)cJet_energy->GetCumulative();
    cJet_energy_cumulative->Scale(1.0/cJet_energy_cumulative->GetMaximum());

    TH1F *cJet_sigmaKT_cumulative = (TH1F *)cJet_sigmaKT->GetCumulative();
    cJet_sigmaKT_cumulative->Scale(1.0/cJet_sigmaKT_cumulative->GetMaximum());

    TH1F *cJet_MaxRho_cumulative = (TH1F *)cJet_MaxRho->GetCumulative();
    cJet_MaxRho_cumulative->Scale(1.0/cJet_MaxRho_cumulative->GetMaximum());

    TH1F *cJet_vertex_cumulative = (TH1F *)cJet_vertex->GetCumulative();
    cJet_vertex_cumulative->Scale(1.0/cJet_vertex_cumulative->GetMaximum());

    TH1F *cJet_angularity_cumulative = (TH1F *)cJet_angularity->GetCumulative();
    cJet_angularity_cumulative->Scale(1.0/cJet_angularity_cumulative->GetMaximum());

    //---------------------------------------------------------------------------------------------------------

    TH1F *cJet_pT_cumulative_bkg = (TH1F *)cJet_pT_bkg->GetCumulative();
    cJet_pT_cumulative_bkg->Scale(1.0/cJet_pT_cumulative_bkg->GetMaximum());

    TH1F *cJet_nConst_cumulative_bkg = (TH1F *)cJet_nConst_bkg->GetCumulative();
    cJet_nConst_cumulative_bkg->Scale(1.0/cJet_nConst_cumulative_bkg->GetMaximum());

    TH1F *cJet_leadConst_pT_cumulative_bkg = (TH1F *)cJet_leadConst_pT_bkg->GetCumulative();
    cJet_leadConst_pT_cumulative_bkg->Scale(1.0/cJet_leadConst_pT_cumulative_bkg->GetMaximum());

    TH1F *cJet_eta_cumulative_bkg = (TH1F *)cJet_eta_bkg->GetCumulative();
    cJet_eta_cumulative_bkg->Scale(1.0/cJet_eta_cumulative_bkg->GetMaximum());

    TH1F *cJet_phi_cumulative_bkg = (TH1F *)cJet_phi_bkg->GetCumulative();
    cJet_phi_cumulative_bkg->Scale(1.0/cJet_phi_cumulative_bkg->GetMaximum());

    TH1F *cJet_energy_cumulative_bkg = (TH1F *)cJet_energy_bkg->GetCumulative();
    cJet_energy_cumulative_bkg->Scale(1.0/cJet_energy_cumulative_bkg->GetMaximum());

    TH1F *cJet_sigmaKT_cumulative_bkg = (TH1F *)cJet_sigmaKT_bkg->GetCumulative();
    cJet_sigmaKT_cumulative_bkg->Scale(1.0/cJet_sigmaKT_cumulative_bkg->GetMaximum());

    TH1F *cJet_MaxRho_cumulative_bkg = (TH1F *)cJet_MaxRho_bkg->GetCumulative();
    cJet_MaxRho_cumulative_bkg->Scale(1.0/cJet_MaxRho_cumulative_bkg->GetMaximum());

    TH1F *cJet_vertex_cumulative_bkg = (TH1F *)cJet_vertex_bkg->GetCumulative();
    cJet_vertex_cumulative_bkg->Scale(1.0/cJet_vertex_cumulative_bkg->GetMaximum());

    TH1F *cJet_angularity_cumulative_bkg = (TH1F *)cJet_angularity_bkg->GetCumulative();
    cJet_angularity_cumulative_bkg->Scale(1.0/cJet_angularity_cumulative_bkg->GetMaximum());

    //---------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------

    TH1F *sbarJet_pT_cumulative = (TH1F *)sbarJet_pT->GetCumulative();
    sbarJet_pT_cumulative->Scale(1.0/sbarJet_pT_cumulative->GetMaximum());

    TH1F *sbarJet_nConst_cumulative = (TH1F *)sbarJet_nConst->GetCumulative();
    sbarJet_nConst_cumulative->Scale(1.0/sbarJet_nConst_cumulative->GetMaximum());

    TH1F *sbarJet_leadConst_pT_cumulative = (TH1F *)sbarJet_leadConst_pT->GetCumulative();
    sbarJet_leadConst_pT_cumulative->Scale(1.0/sbarJet_leadConst_pT_cumulative->GetMaximum());

    TH1F *sbarJet_eta_cumulative = (TH1F *)sbarJet_eta->GetCumulative();
    sbarJet_eta_cumulative->Scale(1.0/sbarJet_eta_cumulative->GetMaximum());

    TH1F *sbarJet_phi_cumulative = (TH1F *)sbarJet_phi->GetCumulative();
    sbarJet_phi_cumulative->Scale(1.0/sbarJet_phi_cumulative->GetMaximum());

    TH1F *sbarJet_energy_cumulative = (TH1F *)sbarJet_energy->GetCumulative();
    sbarJet_energy_cumulative->Scale(1.0/sbarJet_energy_cumulative->GetMaximum());

    TH1F *sbarJet_sigmaKT_cumulative = (TH1F *)sbarJet_sigmaKT->GetCumulative();
    sbarJet_sigmaKT_cumulative->Scale(1.0/sbarJet_sigmaKT_cumulative->GetMaximum());

    TH1F *sbarJet_MaxRho_cumulative = (TH1F *)sbarJet_MaxRho->GetCumulative();
    sbarJet_MaxRho_cumulative->Scale(1.0/sbarJet_MaxRho_cumulative->GetMaximum());

    TH1F *sbarJet_vertex_cumulative = (TH1F *)sbarJet_vertex->GetCumulative();
    sbarJet_vertex_cumulative->Scale(1.0/sbarJet_vertex_cumulative->GetMaximum());

    TH1F *sbarJet_angularity_cumulative = (TH1F *)sbarJet_angularity->GetCumulative();
    sbarJet_angularity_cumulative->Scale(1.0/sbarJet_angularity_cumulative->GetMaximum());

    //---------------------------------------------------------------------------------------------------------

    TH1F *sbarJet_pT_cumulative_bkg = (TH1F *)sbarJet_pT_bkg->GetCumulative();
    sbarJet_pT_cumulative_bkg->Scale(1.0/sbarJet_pT_cumulative_bkg->GetMaximum());

    TH1F *sbarJet_nConst_cumulative_bkg = (TH1F *)sbarJet_nConst_bkg->GetCumulative();
    sbarJet_nConst_cumulative_bkg->Scale(1.0/sbarJet_nConst_cumulative_bkg->GetMaximum());

    TH1F *sbarJet_leadConst_pT_cumulative_bkg = (TH1F *)sbarJet_leadConst_pT_bkg->GetCumulative();
    sbarJet_leadConst_pT_cumulative_bkg->Scale(1.0/sbarJet_leadConst_pT_cumulative_bkg->GetMaximum());

    TH1F *sbarJet_eta_cumulative_bkg = (TH1F *)sbarJet_eta_bkg->GetCumulative();
    sbarJet_eta_cumulative_bkg->Scale(1.0/sbarJet_eta_cumulative_bkg->GetMaximum());

    TH1F *sbarJet_phi_cumulative_bkg = (TH1F *)sbarJet_phi_bkg->GetCumulative();
    sbarJet_phi_cumulative_bkg->Scale(1.0/sbarJet_phi_cumulative_bkg->GetMaximum());

    TH1F *sbarJet_energy_cumulative_bkg = (TH1F *)sbarJet_energy_bkg->GetCumulative();
    sbarJet_energy_cumulative_bkg->Scale(1.0/sbarJet_energy_cumulative_bkg->GetMaximum());

    TH1F *sbarJet_sigmaKT_cumulative_bkg = (TH1F *)sbarJet_sigmaKT_bkg->GetCumulative();
    sbarJet_sigmaKT_cumulative_bkg->Scale(1.0/sbarJet_sigmaKT_cumulative_bkg->GetMaximum());

    TH1F *sbarJet_MaxRho_cumulative_bkg = (TH1F *)sbarJet_MaxRho_bkg->GetCumulative();
    sbarJet_MaxRho_cumulative_bkg->Scale(1.0/sbarJet_MaxRho_cumulative_bkg->GetMaximum());

    TH1F *sbarJet_vertex_cumulative_bkg = (TH1F *)sbarJet_vertex_bkg->GetCumulative();
    sbarJet_vertex_cumulative_bkg->Scale(1.0/sbarJet_vertex_cumulative_bkg->GetMaximum());

    TH1F *sbarJet_angularity_cumulative_bkg = (TH1F *)sbarJet_angularity_bkg->GetCumulative();
    sbarJet_angularity_cumulative_bkg->Scale(1.0/sbarJet_angularity_cumulative_bkg->GetMaximum());

    //---------------------------------------------------------------------------------------------------------
    // Determination of cuts
    //---------------------------------------------------------------------------------------------------------

   std::cout << "Jets from c quark" << std::endl;

    Float_t cut_pT_c = cut_binCalculation(cJet_pT_cumulative, cJet_pT_cumulative_bkg);
    std::cout << "Transverse momentum cut: " << cut_pT_c << "." << std::endl;

    Float_t cut_nConst_c = cut_binCalculation(cJet_nConst_cumulative, cJet_nConst_cumulative_bkg);
    std::cout << "Number of constituents cut: " << cut_nConst_c << "." << std::endl;

    Float_t cut_pTLead_c = cut_binCalculation(cJet_leadConst_pT_cumulative, cJet_leadConst_pT_cumulative_bkg);
    std::cout << "Lead constituent transverse momentum cut: " << cut_pTLead_c << "." << std::endl;

    Float_t cut_eta_c = cut_binCalculation(cJet_eta_cumulative, cJet_eta_cumulative_bkg);
    std::cout << "Jet #eta cut: " << cut_eta_c << "." << std::endl;

    Float_t cut_phi_c = cut_binCalculation(cJet_phi_cumulative, cJet_phi_cumulative_bkg);
    std::cout << "Jet #phi cut: " << cut_phi_c << "." << std::endl;

    Float_t cut_energy_c = cut_binCalculation(cJet_energy_cumulative, cJet_energy_cumulative_bkg);
    std::cout << "Jet Energy cut: " << cut_energy_c << "." << std::endl;

    Float_t cut_sigmaKT_c = cut_binCalculation(cJet_sigmaKT_cumulative, cJet_sigmaKT_cumulative_bkg);
    std::cout << "Jet sigma_{KT} cut: " << cut_sigmaKT_c << "." << std::endl;

    Float_t cut_MaxRho_c = cut_binCalculation(cJet_MaxRho_cumulative, cJet_MaxRho_cumulative_bkg);
    std::cout << "Jet Max Rho cut: " << cut_MaxRho_c << "." << std::endl;

    Float_t cut_vertex_c = cut_binCalculation(cJet_vertex_cumulative, cJet_vertex_cumulative_bkg);
    std::cout << "Jet Vertex cut: " << cut_vertex_c << "." << std::endl;

    Float_t cut_angularity_c = cut_binCalculation(cJet_angularity_cumulative, cJet_angularity_cumulative_bkg);
    std::cout << "Jet Angularity cut: " << cut_angularity_c << "." << std::endl;

    std::cout << "========================================================================" << std::endl;

    std::cout << "Jets from sbar quark" << std::endl;

    Float_t cut_pT_sbar = cut_binCalculation(sbarJet_pT_cumulative, sbarJet_pT_cumulative_bkg);
    std::cout << "Transverse momentum cut: " << cut_pT_sbar << "." << std::endl;

    Float_t cut_nConst_sbar = cut_binCalculation(sbarJet_nConst_cumulative, sbarJet_nConst_cumulative_bkg);
    std::cout << "Number of constituents cut: " << cut_nConst_sbar << "." << std::endl;

    Float_t cut_pTLead_sbar = cut_binCalculation(sbarJet_leadConst_pT_cumulative, sbarJet_leadConst_pT_cumulative_bkg);
    std::cout << "Lead constituent transverse momentum cut: " << cut_pTLead_sbar << "." << std::endl;

    Float_t cut_eta_sbar = cut_binCalculation(sbarJet_eta_cumulative, sbarJet_eta_cumulative_bkg);
    std::cout << "Jet #eta cut: " << cut_eta_sbar << "." << std::endl;

    Float_t cut_phi_sbar = cut_binCalculation(sbarJet_phi_cumulative, sbarJet_phi_cumulative_bkg);
    std::cout << "Jet #phi cut: " << cut_phi_sbar << "." << std::endl;

    Float_t cut_energy_sbar = cut_binCalculation(sbarJet_energy_cumulative, sbarJet_energy_cumulative_bkg);
    std::cout << "Jet Energy cut: " << cut_energy_sbar << "." << std::endl;

    Float_t cut_sigmaKT_sbar = cut_binCalculation(sbarJet_sigmaKT_cumulative, sbarJet_sigmaKT_cumulative_bkg);
    std::cout << "Jet sigma_{KT} cut: " << cut_sigmaKT_sbar << "." << std::endl;

    Float_t cut_MaxRho_sbar = cut_binCalculation(sbarJet_MaxRho_cumulative, sbarJet_MaxRho_cumulative_bkg);
    std::cout << "Jet Max Rho cut: " << cut_MaxRho_sbar << "." << std::endl;

    Float_t cut_vertex_sbar = cut_binCalculation(sbarJet_vertex_cumulative, sbarJet_vertex_cumulative_bkg);
    std::cout << "Jet Vertex cut: " << cut_vertex_sbar << "." << std::endl;

    Float_t cut_angularity_sbar = cut_binCalculation(sbarJet_angularity_cumulative, sbarJet_angularity_cumulative_bkg);
    std::cout << "Jet Angularity cut: " << cut_angularity_sbar << "." << std::endl;

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
    c2->Divide(1, 2);

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


    //---------------------------------------------------------------------------------------------------------

    TCanvas *c3 = new TCanvas("c3", "Cumulative Histograms for c quark - Part I", 2500, 2500);
    c3->Divide(3, 2);

    c3->cd(1);
    cJet_pT_cumulative->SetTitle("Cumulative p_{T} of c-quark jets (signal)");
    cJet_pT_cumulative->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    cJet_pT_cumulative->GetYaxis()->SetTitle("Cumulative Frequency");
    cJet_pT_cumulative->SetEntries(100);
    cJet_pT_cumulative->DrawCopy();

    c3->cd(2);
    cJet_nConst_cumulative->SetTitle("Cumulative number of constituents of c-quark jets (signal)");
    cJet_nConst_cumulative->GetXaxis()->SetTitle("Number of Constituents");
    cJet_nConst_cumulative->GetYaxis()->SetTitle("Cumulative Frequency");
    cJet_nConst_cumulative->SetEntries(100);
    cJet_nConst_cumulative->DrawCopy();

    c3->cd(3);
    cJet_leadConst_pT_cumulative->SetTitle("Cumulative p_{T} of leading constituent of c-quark jets (signal)");
    cJet_leadConst_pT_cumulative->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    cJet_leadConst_pT_cumulative->GetYaxis()->SetTitle("Cumulative Frequency");
    cJet_leadConst_pT_cumulative->SetEntries(100);
    cJet_leadConst_pT_cumulative->DrawCopy();

    c3->cd(4);
    cJet_pT_cumulative_bkg->SetTitle("Cumulative p_{T} of c-quark jets (background)");
    cJet_pT_cumulative_bkg->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    cJet_pT_cumulative_bkg->GetYaxis()->SetTitle("Cumulative Frequency");
    cJet_pT_cumulative_bkg->SetEntries(100);
    cJet_pT_cumulative_bkg->DrawCopy();

    c3->cd(5);
    cJet_nConst_cumulative_bkg->SetTitle("Cumulative number of constituents of c-quark jets (background)");
    cJet_nConst_cumulative_bkg->GetXaxis()->SetTitle("Number of Constituents");
    cJet_nConst_cumulative_bkg->GetYaxis()->SetTitle("Cumulative Frequency");
    cJet_nConst_cumulative_bkg->SetEntries(100);
    cJet_nConst_cumulative_bkg->DrawCopy();

    c3->cd(6);
    cJet_leadConst_pT_cumulative_bkg->SetTitle("Cumulative p_{T} of leading constituent of c-quark jets (background)");
    cJet_leadConst_pT_cumulative_bkg->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    cJet_leadConst_pT_cumulative_bkg->GetYaxis()->SetTitle("Cumulative Frequency");
    cJet_leadConst_pT_cumulative_bkg->SetEntries(100);
    cJet_leadConst_pT_cumulative_bkg->DrawCopy();

    TCanvas *c4 = new TCanvas("c4", "Cumulative Histograms for c quark - Part II", 2500, 2500);
    c4->Divide(2, 2);

    c4->cd(1);
    cJet_eta_cumulative->SetTitle("Cumulative #eta of c-quark jets (signal)");
    cJet_eta_cumulative->GetXaxis()->SetTitle("#eta");
    cJet_eta_cumulative->GetYaxis()->SetTitle("Cumulative Frequency");
    cJet_eta_cumulative->SetEntries(100);
    cJet_eta_cumulative->DrawCopy();

    c4->cd(2);
    cJet_phi_cumulative->SetTitle("Cumulative #phi of c-quark jets (signal)");
    cJet_phi_cumulative->GetXaxis()->SetTitle("#phi");
    cJet_phi_cumulative->GetYaxis()->SetTitle("Cumulative Frequency");
    cJet_phi_cumulative->SetEntries(100);
    cJet_phi_cumulative->DrawCopy();

    c4->cd(3);
    cJet_eta_cumulative_bkg->SetTitle("Cumulative #eta of c-quark jets (background)");
    cJet_eta_cumulative_bkg->GetXaxis()->SetTitle("#eta");
    cJet_eta_cumulative_bkg->GetYaxis()->SetTitle("Cumulative Frequency");
    cJet_eta_cumulative_bkg->SetEntries(100);
    cJet_eta_cumulative_bkg->DrawCopy();

    c4->cd(4);
    cJet_phi_cumulative_bkg->SetTitle("Cumulative #phi of c-quark jets (background)");
    cJet_phi_cumulative_bkg->GetXaxis()->SetTitle("#phi");
    cJet_phi_cumulative_bkg->GetYaxis()->SetTitle("Cumulative Frequency");
    cJet_phi_cumulative_bkg->SetEntries(100);
    cJet_phi_cumulative_bkg->DrawCopy();

    //---------------------------------------------------------------------------------------------------------

    TCanvas *c5 = new TCanvas("c5", "Cumulative Histograms for c quark - Part III", 2500, 2500);
    c5->Divide(3, 2);

    c5->cd(1);
    cJet_energy_cumulative->SetTitle("Cumulative Jet Energy from c quark [GeV] (signal)");
    cJet_energy_cumulative->GetXaxis()->SetTitle("Energy [GeV]");
    cJet_energy_cumulative->GetYaxis()->SetTitle("Cumulative Frequency");
    cJet_energy_cumulative->SetEntries(100);
    cJet_energy_cumulative->DrawCopy();

    c5->cd(2);
    cJet_sigmaKT_cumulative->SetTitle("Cumulative Jet sigma_{KT} from c quark (signal)");
    cJet_sigmaKT_cumulative->GetXaxis()->SetTitle("sigma_{KT}");
    cJet_sigmaKT_cumulative->GetYaxis()->SetTitle("Cumulative Frequency");
    cJet_sigmaKT_cumulative->SetEntries(100);
    cJet_sigmaKT_cumulative->DrawCopy();

    c5->cd(3);
    cJet_angularity_cumulative->SetTitle("Cumulative Jet Angularity from c quark (signal)");
    cJet_angularity_cumulative->GetXaxis()->SetTitle("Angularity");
    cJet_angularity_cumulative->GetYaxis()->SetTitle("Cumulative Frequency");
    cJet_angularity_cumulative->SetEntries(100);
    cJet_angularity_cumulative->DrawCopy();

    c5->cd(4);
    cJet_energy_cumulative_bkg->SetTitle("Cumulative Jet Energy from c quark [GeV] (background)");
    cJet_energy_cumulative_bkg->GetXaxis()->SetTitle("Energy [GeV]");
    cJet_energy_cumulative_bkg->GetYaxis()->SetTitle("Cumulative Frequency");
    cJet_energy_cumulative_bkg->SetEntries(100);
    cJet_energy_cumulative_bkg->DrawCopy();

    c5->cd(5);
    cJet_sigmaKT_cumulative_bkg->SetTitle("Cumulative Jet sigma_{KT} from c quark (background)");
    cJet_sigmaKT_cumulative_bkg->GetXaxis()->SetTitle("sigma_{KT}");
    cJet_sigmaKT_cumulative_bkg->GetYaxis()->SetTitle("Cumulative Frequency");
    cJet_sigmaKT_cumulative_bkg->SetEntries(100);
    cJet_sigmaKT_cumulative_bkg->DrawCopy();

    c5->cd(6);
    cJet_angularity_cumulative_bkg->SetTitle("Cumulative Jet Angularity from c quark (background)");
    cJet_angularity_cumulative_bkg->GetXaxis()->SetTitle("Angularity");
    cJet_angularity_cumulative_bkg->GetYaxis()->SetTitle("Cumulative Frequency");
    cJet_angularity_cumulative_bkg->SetEntries(100);
    cJet_angularity_cumulative_bkg->DrawCopy();

    TCanvas *c6 = new TCanvas("c6", "Cumulative Histograms for c quark - Part IV", 2500, 2500);
    c6->Divide(2, 2);

    c6->cd(1);
    cJet_MaxRho_cumulative->SetTitle("Cumulative Jet Max Rho from c quark (signal)");
    cJet_MaxRho_cumulative->GetXaxis()->SetTitle("Max Rho");
    cJet_MaxRho_cumulative->GetYaxis()->SetTitle("Cumulative Frequency");
    cJet_MaxRho_cumulative->SetEntries(100);
    cJet_MaxRho_cumulative->DrawCopy();

    c6->cd(2);
    cJet_vertex_cumulative->SetTitle("Cumulative Jet Vertex from c quark (signal)");
    cJet_vertex_cumulative->GetXaxis()->SetTitle("Vertex");
    cJet_vertex_cumulative->GetYaxis()->SetTitle("Cumulative Frequency");
    cJet_vertex_cumulative->SetEntries(100);
    cJet_vertex_cumulative->DrawCopy();

    c6->cd(3);
    cJet_MaxRho_cumulative_bkg->SetTitle("Cumulative Jet Max Rho from c quark (background)");
    cJet_MaxRho_cumulative_bkg->GetXaxis()->SetTitle("Max Rho");
    cJet_MaxRho_cumulative_bkg->GetYaxis()->SetTitle("Cumulative Frequency");
    cJet_MaxRho_cumulative_bkg->SetEntries(100);
    cJet_MaxRho_cumulative_bkg->DrawCopy();

    c6->cd(4);
    cJet_vertex_cumulative_bkg->SetTitle("Cumulative Jet Vertex from c quark (background)");
    cJet_vertex_cumulative_bkg->GetXaxis()->SetTitle("Vertex");
    cJet_vertex_cumulative_bkg->GetYaxis()->SetTitle("Cumulative Frequency");
    cJet_vertex_cumulative_bkg->SetEntries(100);
    cJet_vertex_cumulative_bkg->DrawCopy();

    //---------------------------------------------------------------------------------------------------------

    TCanvas *c7 = new TCanvas("c7", "Cumulative Histograms for s-bar quark - Part I", 2500, 2500);
    c7->Divide(3, 2);

    c7->cd(1);
    sbarJet_pT_cumulative->SetTitle("Cumulative p_{T} of s-bar-quark jets (signal)");
    sbarJet_pT_cumulative->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    sbarJet_pT_cumulative->GetYaxis()->SetTitle("Cumulative Frequency");
    sbarJet_pT_cumulative->SetEntries(100);
    sbarJet_pT_cumulative->DrawCopy();

    c7->cd(2);
    sbarJet_nConst_cumulative->SetTitle("Cumulative number of constituents of s-bar-quark jets (signal)");
    sbarJet_nConst_cumulative->GetXaxis()->SetTitle("Number of Constituents");
    sbarJet_nConst_cumulative->GetYaxis()->SetTitle("Cumulative Frequency");
    sbarJet_nConst_cumulative->SetEntries(100);
    sbarJet_nConst_cumulative->DrawCopy();

    c7->cd(3);
    sbarJet_leadConst_pT_cumulative->SetTitle("Cumulative p_{T} of leading constituent of s-bar-quark jets (signal)");
    sbarJet_leadConst_pT_cumulative->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    sbarJet_leadConst_pT_cumulative->GetYaxis()->SetTitle("Cumulative Frequency");
    sbarJet_leadConst_pT_cumulative->SetEntries(100);
    sbarJet_leadConst_pT_cumulative->DrawCopy();

    c7->cd(4);
    sbarJet_pT_cumulative_bkg->SetTitle("Cumulative p_{T} of s-bar-quark jets (background)");
    sbarJet_pT_cumulative_bkg->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    sbarJet_pT_cumulative_bkg->GetYaxis()->SetTitle("Cumulative Frequency");
    sbarJet_pT_cumulative_bkg->SetEntries(100);
    sbarJet_pT_cumulative_bkg->DrawCopy();

    c7->cd(5);
    sbarJet_nConst_cumulative_bkg->SetTitle("Cumulative number of constituents of s-bar-quark jets (background)");
    sbarJet_nConst_cumulative_bkg->GetXaxis()->SetTitle("Number of Constituents");
    sbarJet_nConst_cumulative_bkg->GetYaxis()->SetTitle("Cumulative Frequency");
    sbarJet_nConst_cumulative_bkg->SetEntries(100);
    sbarJet_nConst_cumulative_bkg->DrawCopy();

    c7->cd(6);
    sbarJet_leadConst_pT_cumulative_bkg->SetTitle("Cumulative p_{T} of leading constituent of s-bar-quark jets (background)");
    sbarJet_leadConst_pT_cumulative_bkg->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    sbarJet_leadConst_pT_cumulative_bkg->GetYaxis()->SetTitle("Cumulative Frequency");
    sbarJet_leadConst_pT_cumulative_bkg->SetEntries(100);
    sbarJet_leadConst_pT_cumulative_bkg->DrawCopy();

    TCanvas *c8 = new TCanvas("c8", "Cumulative Histograms for s-bar quark - Part II", 2500, 2500);
    c8->Divide(2, 2);

    c8->cd(1);
    sbarJet_eta_cumulative->SetTitle("Cumulative #eta of s-bar-quark jets (signal)");
    sbarJet_eta_cumulative->GetXaxis()->SetTitle("#eta");
    sbarJet_eta_cumulative->GetYaxis()->SetTitle("Cumulative Frequency");
    sbarJet_eta_cumulative->SetEntries(100);
    sbarJet_eta_cumulative->DrawCopy();

    c8->cd(2);
    sbarJet_phi_cumulative->SetTitle("Cumulative #phi of s-bar-quark jets (signal)");
    sbarJet_phi_cumulative->GetXaxis()->SetTitle("#phi");
    sbarJet_phi_cumulative->GetYaxis()->SetTitle("Cumulative Frequency");
    sbarJet_phi_cumulative->SetEntries(100);
    sbarJet_phi_cumulative->DrawCopy();

    c8->cd(3);
    sbarJet_eta_cumulative_bkg->SetTitle("Cumulative #eta of s-bar-quark jets (background)");
    sbarJet_eta_cumulative_bkg->GetXaxis()->SetTitle("#eta");
    sbarJet_eta_cumulative_bkg->GetYaxis()->SetTitle("Cumulative Frequency");
    sbarJet_eta_cumulative_bkg->SetEntries(100);
    sbarJet_eta_cumulative_bkg->DrawCopy();

    c8->cd(4);
    sbarJet_phi_cumulative_bkg->SetTitle("Cumulative #phi of s-bar-quark jets (background)");
    sbarJet_phi_cumulative_bkg->GetXaxis()->SetTitle("#phi");
    sbarJet_phi_cumulative_bkg->GetYaxis()->SetTitle("Cumulative Frequency");
    sbarJet_phi_cumulative_bkg->SetEntries(100);
    sbarJet_phi_cumulative_bkg->DrawCopy();

    //---------------------------------------------------------------------------------------------------------

    TCanvas *c9 = new TCanvas("c9", "Cumulative Histograms for s-bar quark - Part III", 2500, 2500);
    c9->Divide(3, 2);

    c9->cd(1);
    sbarJet_energy_cumulative->SetTitle("Cumulative Jet Energy from s-bar quark [GeV] (signal)");
    sbarJet_energy_cumulative->GetXaxis()->SetTitle("Energy [GeV]");
    sbarJet_energy_cumulative->GetYaxis()->SetTitle("Cumulative Frequency");
    sbarJet_energy_cumulative->SetEntries(100);
    sbarJet_energy_cumulative->DrawCopy();

    c9->cd(2);
    sbarJet_sigmaKT_cumulative->SetTitle("Cumulative Jet sigma_{KT} from s-bar quark (signal)");
    sbarJet_sigmaKT_cumulative->GetXaxis()->SetTitle("sigma_{KT}");
    sbarJet_sigmaKT_cumulative->GetYaxis()->SetTitle("Cumulative Frequency");
    sbarJet_sigmaKT_cumulative->SetEntries(100);
    sbarJet_sigmaKT_cumulative->DrawCopy();

    c9->cd(3);
    sbarJet_angularity_cumulative->SetTitle("Cumulative Jet Angularity from s-bar quark (signal)");
    sbarJet_angularity_cumulative->GetXaxis()->SetTitle("Angularity");
    sbarJet_angularity_cumulative->GetYaxis()->SetTitle("Cumulative Frequency");
    sbarJet_angularity_cumulative->SetEntries(100);
    sbarJet_angularity_cumulative->DrawCopy();

    c9->cd(4);
    sbarJet_energy_cumulative_bkg->SetTitle("Cumulative Jet Energy from s-bar quark [GeV] (background)");
    sbarJet_energy_cumulative_bkg->GetXaxis()->SetTitle("Energy [GeV]");
    sbarJet_energy_cumulative_bkg->GetYaxis()->SetTitle("Cumulative Frequency");
    sbarJet_energy_cumulative_bkg->SetEntries(100);
    sbarJet_energy_cumulative_bkg->DrawCopy();

    c9->cd(5);
    sbarJet_sigmaKT_cumulative_bkg->SetTitle("Cumulative Jet sigma_{KT} from s-bar quark (background)");
    sbarJet_sigmaKT_cumulative_bkg->GetXaxis()->SetTitle("sigma_{KT}");
    sbarJet_sigmaKT_cumulative_bkg->GetYaxis()->SetTitle("Cumulative Frequency");
    sbarJet_sigmaKT_cumulative_bkg->SetEntries(100);
    sbarJet_sigmaKT_cumulative_bkg->DrawCopy();

    c9->cd(6);
    sbarJet_angularity_cumulative_bkg->SetTitle("Cumulative Jet Angularity from s-bar quark (background)");
    sbarJet_angularity_cumulative_bkg->GetXaxis()->SetTitle("Angularity");
    sbarJet_angularity_cumulative_bkg->GetYaxis()->SetTitle("Cumulative Frequency");
    sbarJet_angularity_cumulative_bkg->SetEntries(100);
    sbarJet_angularity_cumulative_bkg->DrawCopy();


    TCanvas *c10 = new TCanvas("c10", "Cumulative Histograms for s-bar quark - Part IV", 2500, 2500);
    c10->Divide(2, 2);

    c10->cd(1);
    sbarJet_MaxRho_cumulative->SetTitle("Cumulative Jet Max Rho from s-bar quark (signal)");
    sbarJet_MaxRho_cumulative->GetXaxis()->SetTitle("Max Rho");
    sbarJet_MaxRho_cumulative->GetYaxis()->SetTitle("Cumulative Frequency");
    sbarJet_MaxRho_cumulative->SetEntries(100);
    sbarJet_MaxRho_cumulative->DrawCopy();

    c10->cd(2);
    sbarJet_vertex_cumulative->SetTitle("Cumulative Jet Vertex from s-bar quark (signal)");
    sbarJet_vertex_cumulative->GetXaxis()->SetTitle("Vertex");
    sbarJet_vertex_cumulative->GetYaxis()->SetTitle("Cumulative Frequency");
    sbarJet_vertex_cumulative->SetEntries(100);
    sbarJet_vertex_cumulative->DrawCopy();

    c10->cd(3);
    sbarJet_MaxRho_cumulative_bkg->SetTitle("Cumulative Jet Max Rho from s-bar quark (background)");
    sbarJet_MaxRho_cumulative_bkg->GetXaxis()->SetTitle("Max Rho");
    sbarJet_MaxRho_cumulative_bkg->GetYaxis()->SetTitle("Cumulative Frequency");
    sbarJet_MaxRho_cumulative_bkg->SetEntries(100);
    sbarJet_MaxRho_cumulative_bkg->DrawCopy();

    c10->cd(4);
    sbarJet_vertex_cumulative_bkg->SetTitle("Cumulative Jet Vertex from s-bar quark (background)");
    sbarJet_vertex_cumulative_bkg->GetXaxis()->SetTitle("Vertex");
    sbarJet_vertex_cumulative_bkg->GetYaxis()->SetTitle("Cumulative Frequency");
    sbarJet_vertex_cumulative_bkg->SetEntries(100);
    sbarJet_vertex_cumulative_bkg->DrawCopy();
    
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