//---------------------------------------------------------------------------------------------------------//
//                                                                                                         //
// This code reads the original data from the pp collision and is currently being developed to perform an  //
// introductory classical analysis of the data, so it can be compared with the TMVA analysis when          //
// available. In summary, it determines successive and cumulative cuts in the appropriate data histograms  //
// so we can better separate signal and background information.                                            //
//                                                                                                         //
//---------------------------------------------------------------------------------------------------------//

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

Float_t cut_binCalculation(TH1F* signalHist, TH1F* backgroundHist, TH1F* ratioHist)
{
  Float_t minRatio = std::numeric_limits<Float_t>::infinity();
  Int_t minRatioBin = -1;

  Int_t nBins /*Same for both cumulative histograms*/ = signalHist->GetNbinsX();

  for (Int_t i = 0; i <= nBins; i++)
  {
    Double_t signalValue = signalHist->GetBinContent(i);
    Double_t backgroundValue = backgroundHist->GetBinContent(i);
    Double_t ratio = signalValue / backgroundValue;


    if (ratio =! 0)
    {
      ratioHist->Fill(ratio);
    }
    if (ratio < minRatio)
    {
      minRatio = ratio;
      minRatioBin = i;
    }
  }

  Float_t chosenBinValue = signalHist->GetBinCenter(minRatioBin);
  
  return chosenBinValue;
}

void read_wdecayTTree(const char* fileName)
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
  Float_t match_R = 0.1;

  // Criacao dos TLorentzVector

  TLorentzVector vec_s(0,0,0,0);
  TLorentzVector vec_c(0,0,0,0);
  TLorentzVector vec_sbar(0,0,0,0);
  TLorentzVector vec_cbar(0,0,0,0);

  //---------------------------------------------------------------------------------------------------------
  // Inicializacao dos histogramas
  //---------------------------------------------------------------------------------------------------------

  TH1F *ratioHist_pT_c = new TH1F("ratioHist_pT", "Razão cumulativa de p_{T} do sinal pelo fundo", 100, 0, 100);
  TH1F *ratioHist_nConst_c = new TH1F("ratioHist_nConst", "Razão cumulativa de número de constituintes do sinal pelo fundo", 100, 0, 100);
  TH1F *ratioHist_pTLead_c = new TH1F("ratioHist_pTLead", "Razão cumulativa de pT do lead do sinal pelo fundo", 100, 0, 100);

  TH1F *cJet_pT = new TH1F("h1", "p_{T} do jato proveniente do quark c [GeV/c]", 100, 0, 100);
  TH1F *cJet_nConst = new TH1F("h7", "Numero de constituintes do jato proveniente do quark c", 100, 0, 100);
  TH1F *cJet_leadConst_pT = new TH1F("h8", "p_{T} do lead constituent do jato proveniente do quark c [GeV/c]", 100, 0, 100);

  TH1F *cJet_pT_bkg = new TH1F("h12", "p_{T} do jato proveniente do quark c [GeV/c^{2}]", 100, 0, 100);
  TH1F *cJet_nConst_bkg = new TH1F("h72", "Numero de constituintes do jato proveniente do quark c", 100, 0, 100);
  TH1F *cJet_leadConst_pT_bkg = new TH1F("h82", "p_{T} do lead constituent do jato proveniente do quark c [GeV/c]", 100, 0, 100);

  TH1F *sJet_pT = new TH1F("h2", "p_{T} do jato proveniente do quark s [GeV/c]", 100, 0, 100);
  TH1F *cbarJet_pT = new TH1F("h3", "p_{T} do jato proveniente do quark anti-c [GeV/c]", 100, 0, 100);
  TH1F *sbarJet_pT = new TH1F("h4", "p_{T} do jato proveniente do quark anti-s [GeV/c]", 100, 0, 100);

  TH1F *invariantMass_cbar_s = new TH1F("h5", "Espectro de massa invariante dos jatos provenientes dos quarks anti-c e s [GeV/c^{2}]", 100, 0, 100);
  TH1F *invariantMass_c_sbar = new TH1F("h6", "Espectro de massa invariante dos jatos provenientes dos quarks c e anti-s [GeV/c^{2}]", 100, 0, 100);
  


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
            //cout << "sbarJet: " << jetPt << ", " << jetMass << endl;
            vec_sbar = TLorentzVector(jetPx, jetPy, jetPz, jetE);
          }
        }
        //---------------------------------------------------------------------------------------------------------
        // QUARK C ---> Classical analysis interest object
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
            //std::cout << jetPt << ", " << jetNConst << ", " << pT_LeadConst << std::endl;
            cJet_pT->Fill(jetPt);
            cJet_nConst->Fill(jetNConst);
            cJet_leadConst_pT->Fill(pT_LeadConst);

            //cout << "cJet " << jetPt << ", " << jetMass << endl;
            vec_c = TLorentzVector(jetPx, jetPy, jetPz, jetE);
          }
          else
          {
            cJet_pT_bkg->Fill(jetPt);
            cJet_nConst_bkg->Fill(jetNConst);
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

  TH1F *cJet_pT_cumulative_bkg = (TH1F *)cJet_pT_bkg->GetCumulative();
  cJet_pT_cumulative_bkg->Scale(1.0/cJet_pT_cumulative_bkg->GetMaximum());

  //---------------------------------------------------------------------------------------------------------

  TH1F *cJet_nConst_cumulative = (TH1F *)cJet_nConst->GetCumulative();
  cJet_nConst_cumulative->Scale(1.0/cJet_nConst_cumulative->GetMaximum());

  TH1F *cJet_nConst_cumulative_bkg = (TH1F *)cJet_nConst_bkg->GetCumulative();
  cJet_nConst_cumulative_bkg->Scale(1.0/cJet_nConst_cumulative_bkg->GetMaximum());

  //---------------------------------------------------------------------------------------------------------

  TH1F *cJet_leadConst_pT_cumulative = (TH1F *)cJet_leadConst_pT->GetCumulative();
  cJet_leadConst_pT_cumulative->Scale(1.0/cJet_leadConst_pT_cumulative->GetMaximum());

  TH1F *cJet_leadConst_pT_cumulative_bkg = (TH1F *)cJet_leadConst_pT_bkg->GetCumulative();
  cJet_leadConst_pT_cumulative_bkg->Scale(1.0/cJet_leadConst_pT_cumulative_bkg->GetMaximum());

  std::string cumulativeHistName = cJet_pT_cumulative->GetName();

  //---------------------------------------------------------------------------------------------------------
  // Determination of cuts
  //---------------------------------------------------------------------------------------------------------

  Float_t cut_pT_c = cut_binCalculation(cJet_pT_cumulative, cJet_pT_cumulative_bkg, ratioHist_pT_c);
  std::cout << "Transverse momentum cut: " << cut_pT_c << "." << std::endl;

  Float_t cut_nConst_c = cut_binCalculation(cJet_nConst_cumulative, cJet_nConst_cumulative_bkg, ratioHist_nConst_c);
  std::cout << "Number of constituents cut: " << cut_nConst_c << "." << std::endl;

  Float_t cut_pTLead_c = cut_binCalculation(cJet_leadConst_pT_cumulative, cJet_leadConst_pT_cumulative_bkg, ratioHist_pTLead_c);
  std::cout << "Lead constituent transverse momentum cut: " << cut_pTLead_c << "." << std::endl;

  //---------------------------------------------------------------------------------------------------------
  // Canvas creation and plotting
  //---------------------------------------------------------------------------------------------------------
  
  /* pT and invariant mass histograms

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

  TCanvas *c2 = new TCanvas("c2", "Invariant mass distributions for W", 2500, 2500);
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
  
  */

  //---------------------------------------------------------------------------------------------------------
  
  TCanvas *c3 = new TCanvas("c3", "Cumulative histograms", 2500, 2500);
  c3->Divide(3, 2);

  c3->cd(1);
  cJet_pT_cumulative->SetTitle("c Quark Jet's Cumulative Transverse Momentum (Signal)");
  cJet_pT_cumulative->GetXaxis()->SetTitle("Transverse Momentum [GeV/c]");
  cJet_pT_cumulative->GetYaxis()->SetTitle("Cumulated Frequency");
  cJet_pT_cumulative->Draw();

  c3->cd(2);
  cJet_nConst_cumulative->SetTitle("c Quark Jet's Cumulative Numeber of Constituents (Signal)");
  cJet_nConst_cumulative->GetXaxis()->SetTitle("Number of Constituents");
  cJet_nConst_cumulative->GetYaxis()->SetTitle("Cumulated Frequency");
  cJet_nConst_cumulative->Draw();
  
  c3->cd(3);
  cJet_leadConst_pT_cumulative->SetTitle("c Quark Jet's Cumulative pT of the Lead Constituent (Signal)");
  cJet_leadConst_pT_cumulative->GetXaxis()->SetTitle("Transverse Momentum [GeV/c]");
  cJet_leadConst_pT_cumulative->GetYaxis()->SetTitle("Cumulated Frequency");
  cJet_leadConst_pT_cumulative->Draw();

  c3->cd(4);
  cJet_pT_cumulative_bkg->SetTitle("c Quark Jet's Cumulative Transverse Momentum (Background)");
  cJet_pT_cumulative_bkg->GetXaxis()->SetTitle("Transverse Momentum [GeV/c]");
  cJet_pT_cumulative_bkg->GetYaxis()->SetTitle("Cumulated Frequency");
  cJet_pT_cumulative_bkg->Draw();

  c3->cd(5);
  cJet_nConst_cumulative_bkg->SetTitle("c Quark Jet's Cumulative Numeber of Constituents (Background)");
  cJet_nConst_cumulative_bkg->GetXaxis()->SetTitle("Number of Constituents");
  cJet_nConst_cumulative_bkg->GetYaxis()->SetTitle("Cumulated Frequency");
  cJet_nConst_cumulative_bkg->Draw();
  
  c3->cd(6);
  cJet_leadConst_pT_cumulative_bkg->SetTitle("c Quark Jet's Cumulative pT of the Lead Constituent (Background)");
  cJet_leadConst_pT_cumulative_bkg->GetXaxis()->SetTitle("Transverse Momentum [GeV/c]");
  cJet_leadConst_pT_cumulative_bkg->GetYaxis()->SetTitle("Cumulated Frequency");
  cJet_leadConst_pT_cumulative_bkg->Draw();

  //---------------------------------------------------------------------------------------------------------

  TCanvas *c4 = new TCanvas("c4", "Cumulative ratio histograms", 2500, 2500); // "Ratio" refers to the division of signal data by bkg data
  c4->Divide(1,3);

  c4->cd(1);
  ratioHist_pT_c->SetTitle("c Quark Jet's p_{T} Cumulative Ratio");
  ratioHist_pT_c->SetXaxis()->SetTitle("Transverse Momentum [Gev/c]");
  ratioHist_pT_c->SetYaxis()->SetTitle("Frequency");
  ratioHist_pT_c->Draw();

  c4->cd(2);
  ratioHist_nConst_c->SetTitle("c Quark Jet's Number of Constituents Cumulative Ratio");
  ratioHist_nConst_c->SetXaxis()->SetTitle("Number of Constituents");
  ratioHist_nConst_c->SetYaxis()->SetTitle("Frequency");
  ratioHist_nConst_c->Draw();

  c4->cd(3);
  ratioHist_pT_c->SetTitle("c Quark Jet's Lead Constituent p_{T} Cumulative Ratio");
  ratioHist_pT_c->SetXaxis()->SetTitle("Transverse Momentum [Gev/c]");
  ratioHist_pT_c->SetYaxis()->SetTitle("Frequency");
  ratioHist_pT_c->Draw();

  //---------------------------------------------------------------------------------------------------------

  file->Close();
}
