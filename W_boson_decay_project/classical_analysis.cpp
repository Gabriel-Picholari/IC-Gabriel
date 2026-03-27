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

void MakeInvariantMassHist(const std::vector<fastjet::PseudoJet> &jet, TH1F* invariantMass)
{

    for (size_t ic = 0; ic < jet.size(); ++ic)
    {
        TLorentzVector vc(jet[ic].px(), jet[ic].py(), jet[ic].pz(), jet[ic].E());

        for (size_t is = ic + 1; is < jet.size(); ++is)
        {
            TLorentzVector vs(jet[is].px(), jet[is].py(), jet[is].pz(), jet[is].E());

            TLorentzVector vW = vc + vs;
            invariantMass->Fill(vW.M());
        }
    }
}

class JetInfo : public fastjet::PseudoJet::UserInfoBase
{
    public:
        JetInfo(const Float_t& vx, const Float_t& vy, const Float_t& vz) : fpVx(vx), fpVy(vy), fpVz(vz) {}

        void setVx(Float_t vx) { fpVx = vx; }
        Float_t getVx() const { return fpVx; }

        void setVy(Float_t vy) { fpVy = vy; }
        Float_t getVy() const { return fpVy; }

        void setVz(Float_t vz) { fpVz = vz; }
        Float_t getVz() const { return fpVz; }

    private:
        Float_t fpVx;
        Float_t fpVy;
        Float_t fpVz;
};

void classical_analysis(const char* fileName)
{

    gSystem->Load("libEG");
    gSystem->Load("libEGPythia8");

    //---------------------------------------------------------------------------------------------------------
    // Initialization of variables
    //---------------------------------------------------------------------------------------------------------

    Float_t fpPt, fpEta, fpPhi, fpE, fpPx, fpPy, fpPz, fpMass;
    Float_t jetPt, jetEta, jetPhi, jetE, jetPx, jetPy, jetPz, jetMass, jetNConst, pT_LeadConst;
    Float_t vertex_x, vertex_y, vertex_z;
    
    TLorentzVector vec_s(0,0,0,0);
    TLorentzVector vec_c(0,0,0,0);

    //---------------------------------------------------------------------------------------------------------
    // Initialization of histograms
    //---------------------------------------------------------------------------------------------------------

    TH1F *no_cuts_invariantMass = new TH1F("invariant_mass", "W^{+-} invariant mass spectrum [GeV/c^{2}]", 100, 0, 120);
    TH1F *pT_cut_invariantMass = new TH1F("invariant_mass_pTcut", "W^{+-} invariant mass spectrum with p_{T} cut [GeV/c^{2}]", 100, 0, 120);
    TH1F *nConst_cut_invariantMass = new TH1F("invariant_mass_nConstCut", "W^{+-} invariant mass spectrum with number of constituents cut [GeV/c^{2}]", 100, 0, 120);
    TH1F *nRho_cut_invariantMass = new TH1F("invariant_mass_nRhoCut", "W^{+-} invariant mass spectrum with nRho cut [GeV/c^{2}]", 100, 0, 120);
    
    TH1F *pT_and_nConst_cuts_invariantMass = new TH1F("invariant_mass_pT_and_nConstCuts", "W^{+-} invariant mass spectrum with p_{T} and number of constituents cuts [GeV/c^{2}]", 100, 0, 120);
    TH1F *combined_cuts_invariantMass = new TH1F("invariant_mass_combinedCuts", "W^{+-} invariant mass spectrum with combined cuts [GeV/c^{2}]", 100, 0, 120);

    TH1F *nRho_distribution = new TH1F("nRho_distribution", "N_{rho} jet distribution", 10, 0, 10);

    //---------------------------------------------------------------------------------------------------------
    // Initializations and FastJet configurations:
    //---------------------------------------------------------------------------------------------------------

    Float_t jetR = 0.7;

    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, jetR);

    std::vector<fastjet::PseudoJet> particles_fastjet;
    std::vector<fastjet::PseudoJet> jets;

    //---------------------------------------------------------------------------------------------------------
    // Initialization of the .root file and TTrees
    //---------------------------------------------------------------------------------------------------------

    TFile *file = TFile::Open(fileName, "READ");
    TTree *ttree = dynamic_cast<TTree *>(file->Get("W decay TTree 2"));

    TClonesArray *jets_array = new TClonesArray("MyJet");
    TClonesArray *quarks = new TClonesArray("MyQuark");

    ttree->SetBranchAddress("jets_array", &jets_array);
    ttree->SetBranchAddress("quarks", &quarks);

    Long64_t ne = ttree->GetEntries();

    //---------------------------------------------------------------------------------------------------------
    // Event loop equivalent
    //---------------------------------------------------------------------------------------------------------

    for ( Long64_t ni = 0; ni < ne; ni++)
    {
        ttree->GetEntry(ni);

        vec_s.SetPtEtaPhiM(10,0,0,(3.141592));
        vec_c.SetPtEtaPhiM(10,0,0,(3.141592));

        std::vector<fastjet::PseudoJet> tagged_c_jets;
        std::vector<fastjet::PseudoJet>tagged_s_jets;

        Int_t count = 0;

        //---------------------------------------------------------------------------------------------------------
        // Quark information ( it's unique per event - a single pair of cbar(c) - s(sbar) quarks )
        //---------------------------------------------------------------------------------------------------------

        Float_t charmPt = 0;
        Float_t strangePt = 0;

        //std::cout << quarks->GetEntries() << std::endl;
        for (Int_t nk = 0; nk < quarks->GetEntries(); nk++)
        {
            MyQuark *mq = static_cast<MyQuark *>(quarks->At(nk));
            Int_t quarkPdg = mq->qPdg;
            
            if ( abs(quarkPdg) == 4) charmPt = mq->qpT;
            //std::cout << charmPt << std::endl;
            if ( abs(quarkPdg) == 3) strangePt = mq->qpT;
            //std::cout << strangePt << std::endl;
        }

        std::vector<fastjet::PseudoJet> pT_cut_selectedJets;
        std::vector<fastjet::PseudoJet> nConst_cut_selectedJets;
        std::vector<fastjet::PseudoJet> nRho_cut_selectedJets;
        std::vector<fastjet::PseudoJet> pT_and_nConst_cuts_selectedJets;
        std::vector<fastjet::PseudoJet> combined_cuts_selectedJets;


        //---------------------------------------------------------------------------------------------------------
        // Particle loop equivalent
        //---------------------------------------------------------------------------------------------------------

        for (Int_t nj = 0; nj < jets_array->GetEntries(); nj++)
        {
            MyJet *fp = static_cast<MyJet *>(jets_array->At(nj));

            
            fpPx = fp->fPx;
            fpPy = fp->fPy;
            fpPz = fp->fPz;
            fpE = fp->fE;
            vertex_x = fp->fVx;
            vertex_y = fp->fVy;
            vertex_z = fp->fVz;

            fastjet::PseudoJet particle(fpPx, fpPy, fpPz, fpE);
            
            JetInfo* jetInfo = new JetInfo(vertex_x, vertex_y, vertex_z);
            particle.set_user_info(jetInfo);

            particles_fastjet.push_back(particle);
        }

        fastjet::ClusterSequence clusterSeq(particles_fastjet, jet_def);
        jets = clusterSeq.inclusive_jets();

        for (const fastjet::PseudoJet& jet : jets)
        {
            //---------------------------------------------------------------------------------------------------------
            // Collecting jet's general information
            //---------------------------------------------------------------------------------------------------------

            //std::cout << "Jet number: " << count << std::endl;
            count++;

            jetPt = jet.pt();
            if (jetPt < 5) continue; // Standard cut on jet pT

            jetEta = jet.eta();

            Float_t absEta = TMath::Abs(jetEta);
            //if (absEta > 2) continue; // Standard cut on jet eta
        
            // It seems that the cut on eta causes the number of jet entries with null pT to increase drastically

            jetPhi = jet.phi();
            jetMass = jet.m();
            jetPx = jet.px();
            jetPy = jet.py();
            jetPz = jet.pz();
            jetE = jet.E();
            jetNConst = jet.constituents().size();
            
            pT_LeadConst = 0;
            for (const fastjet::PseudoJet &constituent : jet.constituents())
            {
                if (constituent.pt() > pT_LeadConst)
                {
                    pT_LeadConst = constituent.pt();
                }
            }

            Int_t nRho = 0;
            Int_t maxRho = 2;
            pT_LeadConst = 0;
            Float_t rhoLowerBound = 0.01;
            Float_t rhoUpperBound = 2;

            for (const fastjet::PseudoJet &constituent : jet.constituents())
            {
                Float_t vx = constituent.user_info<JetInfo>().getVx();
                Float_t vy = constituent.user_info<JetInfo>().getVy();
                
                Double_t Rho = TMath::Sqrt(pow(vx, 2) + pow(vy, 2));

                if (Rho >= rhoLowerBound && Rho < rhoUpperBound)
                {
                    nRho++;
                }
            }
            
            nRho_distribution->Fill(nRho);

            Float_t pT_cut_pT = jetPt;
            if (pT_cut_pT > 15)
            {
                pT_cut_selectedJets.push_back(jet);
            }

            Float_t nConst_cut_nConst = jetNConst;
            if (nConst_cut_nConst > 15)
            {
                nConst_cut_selectedJets.push_back(jet);
            }

            Float_t nRho_cut_nRho = nRho;
            if (nRho_cut_nRho > 1)
            {
                nRho_cut_selectedJets.push_back(jet);
            }

            if (pT_cut_pT > 15 && nConst_cut_nConst > 15)
            {
                pT_and_nConst_cuts_selectedJets.push_back(jet);
            }

            if (pT_cut_pT > 15 && nConst_cut_nConst > 15 && nRho_cut_nRho > 1)
            {
                combined_cuts_selectedJets.push_back(jet);
            }
        }

        MakeInvariantMassHist(jets, no_cuts_invariantMass);
        MakeInvariantMassHist(pT_cut_selectedJets, pT_cut_invariantMass);
        MakeInvariantMassHist(nConst_cut_selectedJets, nConst_cut_invariantMass);
        MakeInvariantMassHist(nRho_cut_selectedJets, nRho_cut_invariantMass);
        MakeInvariantMassHist(pT_and_nConst_cuts_selectedJets, pT_and_nConst_cuts_invariantMass);
        MakeInvariantMassHist(combined_cuts_selectedJets, combined_cuts_invariantMass);

        particles_fastjet.clear();
        jets.clear();
        jets_array->Clear();
        quarks->Clear();
    
    } // End of event loop equivalent

    
    //---------------------------------------------------------------------------------------------------------
    // Plotting histograms
    //---------------------------------------------------------------------------------------------------------
    
    TCanvas *c1 = new TCanvas("c1", "Invariant mass distribution with individual cuts", 2500, 2500);
    c1->Divide(2, 2);

    c1->cd(1);
    no_cuts_invariantMass->SetTitle("Boson W^{#pm} no cuts invariant mass spectrum");
    no_cuts_invariantMass->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    no_cuts_invariantMass->GetYaxis()->SetTitle("Frequency");
    no_cuts_invariantMass->DrawCopy();

    c1->cd(2);
    pT_cut_invariantMass->SetTitle("Boson W^{#pm} invariant mass spectrum with p_{T} cut");
    pT_cut_invariantMass->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    pT_cut_invariantMass->GetYaxis()->SetTitle("Frequency");
    pT_cut_invariantMass->DrawCopy();

    c1->cd(3);
    nConst_cut_invariantMass->SetTitle("Boson W^{#pm} invariant mass spectrum with number of constituents cut");
    nConst_cut_invariantMass->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    nConst_cut_invariantMass->GetYaxis()->SetTitle("Frequency");
    nConst_cut_invariantMass->DrawCopy();

    c1->cd(4);
    nRho_cut_invariantMass->SetTitle("Boson W^{#pm} invariant mass spectrum with nRho cut");
    nRho_cut_invariantMass->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    nRho_cut_invariantMass->GetYaxis()->SetTitle("Frequency");
    nRho_cut_invariantMass->DrawCopy();

    TCanvas *c2 = new TCanvas("c2", "Invariant mass distribution with combined cuts", 2500, 2500);
    c2->Divide(2, 1);

    c2->cd(1);
    pT_and_nConst_cuts_invariantMass->SetTitle("Boson W^{#pm} invariant mass spectrum with p_{T} and number of constituents cuts");
    pT_and_nConst_cuts_invariantMass->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    pT_and_nConst_cuts_invariantMass->GetYaxis()->SetTitle("Frequency");
    pT_and_nConst_cuts_invariantMass->DrawCopy();

    c2->cd(2);
    combined_cuts_invariantMass->SetTitle("Boson W^{#pm} invariant mass spectrum with combined cuts");
    combined_cuts_invariantMass->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    combined_cuts_invariantMass->GetYaxis()->SetTitle("Frequency");
    combined_cuts_invariantMass->DrawCopy();

    TCanvas *c3 = new TCanvas("c3", "N_{rho} distribution", 2500, 2500);
    c3->Divide(1, 1);

    c3->cd(1);
    nRho_distribution->SetTitle("N_{rho} distribution");
    nRho_distribution->GetXaxis()->SetTitle("N_{rho}");
    nRho_distribution->GetYaxis()->SetTitle("Frequency");
    nRho_distribution->DrawCopy();

    TFile *second_outputFile = new TFile("modified2_histogramas_jetR_07_fullList.root", "RECREATE");
    no_cuts_invariantMass->Write();
    pT_cut_invariantMass->Write();
    nConst_cut_invariantMass->Write();
    nRho_cut_invariantMass->Write();
    pT_and_nConst_cuts_invariantMass->Write();
    combined_cuts_invariantMass->Write();
    second_outputFile->Close();
    file->Close();
}

