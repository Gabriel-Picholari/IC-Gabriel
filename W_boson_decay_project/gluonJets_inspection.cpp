#include <cmath>
#include "TH1.h"
#include "TH1F.h"
#include "TH1D.h"
#include <string>
#include <fstream>
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

class JetInfo : public fastjet::PseudoJet::UserInfoBase
{
    public:
        JetInfo(const TString& type = "", const Int_t& pdg = 0, const Int_t& motherPdg = 0, const Int_t& secondMotherPdg = 0, const Int_t& thirdMotherPdg = 0) : signalType(type), finalParticlePdg(pdg), finalParticleMotherPdg(motherPdg), finalParticleSecondMotherPdg(secondMotherPdg), finalParticleThirdMotherPdg(thirdMotherPdg){}

        void setSignalType(const TString& type) { signalType = type; }
        TString getSignalType() const { return signalType; }

        void setFinalParticlePdg(Int_t pdg) { finalParticlePdg = pdg; }
        Int_t getFinalParticlePdg() const { return finalParticlePdg; }

        void setFinalParticleMotherPdg(Int_t motherPdg) { finalParticleMotherPdg = motherPdg; }
        Int_t getFinalParticleMotherPdg() const { return finalParticleMotherPdg; }

        void setFinalParticleSecondMotherPdg(Int_t secondMotherPdg) { finalParticleSecondMotherPdg = secondMotherPdg; }
        Int_t getFinalParticleSecondMotherPdg() const { return finalParticleSecondMotherPdg; }

        void setFinalParticleThirdMotherPdg(Int_t thirdMotherPdg) { finalParticleThirdMotherPdg = thirdMotherPdg; }
        Int_t getFinalParticleThirdMotherPdg() const { return finalParticleThirdMotherPdg; }

    private:
        TString signalType;
        Int_t finalParticlePdg;
        Int_t finalParticleMotherPdg;
        Int_t finalParticleSecondMotherPdg;
        Int_t finalParticleThirdMotherPdg;
};

struct JetTagResult
{
    Bool_t isCharmTagged = false;
    Bool_t isStrangeTagged = false;
    Int_t  nCharmConst = 0;
    Int_t  nStrangeConst = 0;
};

JetTagResult classifyJet(
    const fastjet::PseudoJet& jet,
    std::vector<fastjet::PseudoJet>& tagged_c_jets,
    std::vector<fastjet::PseudoJet>& tagged_s_jets,
    TLorentzVector& vec_c,
    TLorentzVector& vec_s
)
{
    JetTagResult result;

    TLorentzVector currentJet(jet.px(), jet.py(), jet.pz(), jet.E());
    if (currentJet.M() < 0) return result;

    // Step one: count the number of constituents of each type
    for (const fastjet::PseudoJet &constituent : jet.constituents())
    {
        TString signalType = constituent.user_info<JetInfo>().getSignalType();

        if (signalType == "charm")   result.nCharmConst++;
        else if (signalType == "strange") result.nStrangeConst++;
    }

    // Step two: classify the jet based on the counts
    if (result.nCharmConst > result.nStrangeConst && result.nCharmConst > 0)
    {
        result.isCharmTagged = true;
        tagged_c_jets.push_back(jet);
        vec_c = currentJet;
    }
    else if (result.nStrangeConst > result.nCharmConst && result.nStrangeConst > 0)
    {
        result.isStrangeTagged = true;
        tagged_s_jets.push_back(jet);
        vec_s = currentJet;
    }
    // In case of a tie or if both counts are zero, the jet is not tagged as either type (ends up as background)

    return result;
}

void gluonJets_inspection(const char* fileName)
{

    gSystem->Load("libEG");
    gSystem->Load("libEGPythia8");

    //---------------------------------------------------------------------------------------------------------
    // Initialization of variables
    //---------------------------------------------------------------------------------------------------------

    Float_t fpPt, fpEta, fpPhi, fpE, fpPx, fpPy, fpPz, fpMass = 0;
    Float_t jetPt, jetEta, jetPhi, jetE, jetPx, jetPy, jetPz, jetMass, jetNConst, pT_LeadConst = 0;
    Float_t maxRho, nVert;
    TString signalType = "";
    Int_t finalParticlePdg = 0;
    Int_t finalParticleMotherPdg = 0;
    Int_t finalParticleSecondMotherPdg = 0;
    Int_t finalParticleThirdMotherPdg = 0;
    
    TLorentzVector vec_s(0,0,0,0);
    TLorentzVector vec_c(0,0,0,0);

    //---------------------------------------------------------------------------------------------------------
    // Initialization of histograms
    //---------------------------------------------------------------------------------------------------------

    TH1F *strange_nConstDistribution_lowerRatio = new TH1F("strange_nConstDistribution_lowerRatio", "Strange jets: N_{const} distribution for jets such that 0.1 < p_{T}^{jet}/p_{T}^{quark} < 0.4", 60, 0, 60);
    TH1F *strange_nConstDistribution_higherRatio = new TH1F("strange_nConstDistribution_higherRatio", "Strange jets: N_{const} distribution for jets such that 0.8 < p_{T}^{jet}/p_{T}^{quark} < 1.2", 60, 0, 60);

    TH1F *strange_deltaEtaDistribution_lowerRatio = new TH1F("strange_deltaEtaDistribution_lowerRatio", "Strange jets: #Delta#eta distribution for jets such that 0.1 < p_{T}^{jet}/p_{T}^{quark} < 0.4", 32, -4, 4);
    TH1F *strange_deltaEtaDistribution_higherRatio = new TH1F("strange_deltaEtaDistribution_higherRatio", "Strange jets: #Delta#eta distribution for jets such that 0.8 < p_{T}^{jet}/p_{T}^{quark} < 1.2", 32, -4, 4);

    TH1F *strange_deltaPhiDistribution_lowerRatio = new TH1F("strange_deltaPhiDistribution_lowerRatio", "Strange jets: #Delta#phi distribution for jets such that 0.1 < p_{T}^{jet}/p_{T}^{quark} < 0.4", 32, -4, 4);
    TH1F *strange_deltaPhiDistribution_higherRatio = new TH1F("strange_deltaPhiDistribution_higherRatio", "Strange jets: #Delta#phi distribution for jets such that 0.8 < p_{T}^{jet}/p_{T}^{quark} < 1.2", 32, -4, 4);

    TH1F *strange_deltaRDistribution_lowerRatio = new TH1F("strange_deltaRDistribution_lowerRatio", "Strange jets: #DeltaR distribution for jets such that 0.1 < p_{T}^{jet}/p_{T}^{quark} < 0.4", 32, 0, 4);
    TH1F *strange_deltaRDistribution_higherRatio = new TH1F("strange_deltaRDistribution_higherRatio", "Strange jets: #DeltaR distribution for jets such that 0.8 < p_{T}^{jet}/p_{T}^{quark} < 1.2", 32, 0, 4);

    TH1F *charm_nConstDistribution_lowerRatio = new TH1F("charm_nConstDistribution_lowerRatio", "Charm jets: N_{const} distribution for jets such that 0.1 < p_{T}^{jet}/p_{T}^{quark} < 0.4", 60, 0, 60);
    TH1F *charm_nConstDistribution_higherRatio = new TH1F("charm_nConstDistribution_higherRatio", "Charm jets: N_{const} distribution for jets such that 0.8 < p_{T}^{jet}/p_{T}^{quark} < 1.2", 60, 0, 60);

    TH1F *charm_deltaEtaDistribution_lowerRatio = new TH1F("charm_deltaEtaDistribution_lowerRatio", "Charm jets: #Delta#eta distribution for jets such that 0.1 < p_{T}^{jet}/p_{T}^{quark} < 0.4", 32, -4, 4);
    TH1F *charm_deltaEtaDistribution_higherRatio = new TH1F("charm_deltaEtaDistribution_higherRatio", "Charm jets: #Delta#eta distribution for jets such that 0.8 < p_{T}^{jet}/p_{T}^{quark} < 1.2", 32, -4, 4);

    TH1F *charm_deltaPhiDistribution_lowerRatio = new TH1F("charm_deltaPhiDistribution_lowerRatio", "Charm jets: #Delta#phi distribution for jets such that 0.1 < p_{T}^{jet}/p_{T}^{quark} < 0.4", 32, -4, 4);
    TH1F *charm_deltaPhiDistribution_higherRatio = new TH1F("charm_deltaPhiDistribution_higherRatio", "Charm jets: #Delta#phi distribution for jets such that 0.8 < p_{T}^{jet}/p_{T}^{quark} < 1.2", 32, -4, 4);

    TH1F *charm_deltaRDistribution_lowerRatio = new TH1F("charm_deltaRDistribution_lowerRatio", "Charm jets: #DeltaR distribution for jets such that 0.1 < p_{T}^{jet}/p_{T}^{quark} < 0.4", 32, 0, 4);
    TH1F *charm_deltaRDistribution_higherRatio = new TH1F("charm_deltaRDistribution_higherRatio", "Charm jets: #DeltaR distribution for jets such that 0.8 < p_{T}^{jet}/p_{T}^{quark} < 1.2", 32, 0, 4);

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

        Float_t charmPt = 0, charmEta = 0, charmPhi = 0;
        Float_t strangePt = 0, strangeEta = 0, strangePhi = 0;

        //std::cout << quarks->GetEntries() << std::endl;
        for (Int_t nk = 0; nk < quarks->GetEntries(); nk++)
        {
            MyQuark *mq = static_cast<MyQuark *>(quarks->At(nk));
            Int_t quarkPdg = mq->qPdg;
            
            if ( abs(quarkPdg) == 4)
            {
                charmPt = mq->qpT;
                charmEta = mq->qEta;
                charmPhi = mq->qPhi;
            }
            if ( abs(quarkPdg) == 3)
            {
                strangePt = mq->qpT;
                strangeEta = mq->qEta;
                strangePhi = mq->qPhi;
            }
        }

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
            signalType = fp->signalType;
            finalParticlePdg = fp->finalParticlePdg;
            finalParticleMotherPdg = fp->finalParticleMotherPdg;
            finalParticleSecondMotherPdg = fp->finalParticleSecondMotherPdg;
            finalParticleThirdMotherPdg = fp->finalParticleThirdMotherPdg;

            fastjet::PseudoJet particle(fpPx, fpPy, fpPz, fpE);
            
            JetInfo* jetInfo = new JetInfo(signalType, finalParticlePdg, finalParticleMotherPdg, finalParticleSecondMotherPdg, finalParticleThirdMotherPdg);
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
            if (jetPt < 10) continue; // Basic cut on jet pT

            jetEta = jet.eta();

            Float_t absEta = TMath::Abs(jetEta);
            if (absEta > 1.3) continue; // Basic cut on jet eta
            
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

            //---------------------------------------------------------------------------------------------------------
            // Jet classification block (based on constituents info)
            //---------------------------------------------------------------------------------------------------------

            JetTagResult tagResult = classifyJet(jet, tagged_c_jets, tagged_s_jets, vec_c, vec_s);

        } // End of individual jet creation

        const std::unordered_set<int> charmPdgSet = {411, 421, 413, 423, 415, 425, 431, 433, 435, 10411, 10421, 413, 423, 10413, 10423, 20413, 20423, 415, 425, 431, 10431, 433, 10433, 20433, 435, 4122, 4222, 4212, 4112, 4224, 4214, 4114, 4232, 4132, 4322, 4312, 4324, 4314, 4332, 4334, 4412, 4422, 4414, 4424, 4432, 4434, 4444};
        const std::unordered_set<int> strangePdgSet = {130, 310, 311, 321, 313, 323, 315, 325, 317, 327, 319, 329, 9000311, 9000321, 10311, 10321, 100311, 100321, 9010311, 9010321, 9020311, 9020321, 313, 323, 10313, 10323, 20313, 20323, 100313, 100323, 9000313, 9000323, 30313, 30323, 315, 325, 9000315, 9000325, 10315, 10325, 20315, 20325, 9010315, 9010325, 9020315, 9020325, 317, 327, 9010317, 9010327, 319, 329, 3122, 3222, 3212, 3112, 3224, 3214, 3114, 3322, 3312, 3324, 3314, 3334};
        
        //---------------------------------------------------------------------------------------------------------
        // --- c-tagged Jets Filling ---
        //---------------------------------------------------------------------------------------------------------
        for (const fastjet::PseudoJet &jet : tagged_c_jets)
        {
            TLorentzVector cJet(jet.px(), jet.py(), jet.pz(), jet.E());
            Float_t charmRatio = jet.pt() / charmPt;

            Bool_t hasCharmConstituent = false;

            for (const fastjet::PseudoJet &constituent : jet.constituents())
            {
                Int_t constituentPdg                = constituent.user_info<JetInfo>().getFinalParticlePdg();
                Int_t constituentMotherPdg          = constituent.user_info<JetInfo>().getFinalParticleMotherPdg();
                Int_t constituentSecondMotherPdg    = constituent.user_info<JetInfo>().getFinalParticleSecondMotherPdg();
                Int_t constituentThirdMotherPdg     = constituent.user_info<JetInfo>().getFinalParticleThirdMotherPdg();

                Int_t abs_constituentPdg                = abs(constituentPdg);
                Int_t abs_constituentMotherPdg          = abs(constituentMotherPdg);
                Int_t abs_constituentSecondMotherPdg    = abs(constituentSecondMotherPdg);
                Int_t abs_constituentThirdMotherPdg     = abs(constituentThirdMotherPdg);

                if (charmPdgSet.count(abs_constituentPdg) || charmPdgSet.count(abs_constituentMotherPdg) || charmPdgSet.count(abs_constituentSecondMotherPdg) || charmPdgSet.count(abs_constituentThirdMotherPdg))
                {
                    hasCharmConstituent = true;
                    break;
                }
            }

            Int_t nConst = jet.constituents().size();
            Float_t deltaEta = (jet.eta() - charmEta);
            Float_t deltaPhi = (jet.phi() - charmPhi);
            Float_t deltaR = TMath::Sqrt(deltaEta * deltaEta + deltaPhi * deltaPhi);

            if (hasCharmConstituent)
            {
                if(charmRatio > 0.1 && charmRatio < 0.4)
                {
                    charm_nConstDistribution_lowerRatio->Fill(nConst);
                    charm_deltaEtaDistribution_lowerRatio->Fill(deltaEta);
                    charm_deltaPhiDistribution_lowerRatio->Fill(deltaPhi);
                    charm_deltaRDistribution_lowerRatio->Fill(deltaR);
                }
                else if(charmRatio > 0.8 && charmRatio < 1.2)
                {
                    charm_nConstDistribution_higherRatio->Fill(nConst);
                    charm_deltaEtaDistribution_higherRatio->Fill(deltaEta);
                    charm_deltaPhiDistribution_higherRatio->Fill(deltaPhi);
                    charm_deltaRDistribution_higherRatio->Fill(deltaR);
                }
            }
        }

        //std::cout << "\n--- s-tagged Jets ---" << std::endl;
        for (const fastjet::PseudoJet &jet : tagged_s_jets)
        {
            TLorentzVector sJet(jet.px(), jet.py(), jet.pz(), jet.E());
            Float_t strangeRatio = jet.pt() / strangePt;

           Bool_t hasStrangeConstituent = false;
           Int_t missingStrangePdg = 0;

            for (const fastjet::PseudoJet &constituent : jet.constituents())
            {
                Int_t constituentPdg                = constituent.user_info<JetInfo>().getFinalParticlePdg();
                Int_t constituentMotherPdg          = constituent.user_info<JetInfo>().getFinalParticleMotherPdg();
                Int_t constituentSecondMotherPdg    = constituent.user_info<JetInfo>().getFinalParticleSecondMotherPdg();
                Int_t constituentThirdMotherPdg     = constituent.user_info<JetInfo>().getFinalParticleThirdMotherPdg();

                Int_t abs_constituentPdg                = abs(constituentPdg);
                Int_t abs_constituentMotherPdg          = abs(constituentMotherPdg);
                Int_t abs_constituentSecondMotherPdg    = abs(constituentSecondMotherPdg);
                Int_t abs_constituentThirdMotherPdg     = abs(constituentThirdMotherPdg);

                if (strangePdgSet.count(abs_constituentPdg) || strangePdgSet.count(abs_constituentMotherPdg) || strangePdgSet.count(abs_constituentSecondMotherPdg) || strangePdgSet.count(abs_constituentThirdMotherPdg))
                {
                    hasStrangeConstituent = true;
                    break;
                }
            }

            Int_t nConst = jet.constituents().size();
            Float_t deltaEta = (jet.eta() - strangeEta);
            Float_t deltaPhi = (jet.phi() - strangePhi);
            Float_t deltaR = TMath::Sqrt(deltaEta * deltaEta + deltaPhi * deltaPhi);

            

            if (hasStrangeConstituent)
            {
                if(strangeRatio > 0.1 && strangeRatio < 0.4)
                {
                    strange_nConstDistribution_lowerRatio->Fill(nConst);
                    strange_deltaEtaDistribution_lowerRatio->Fill(deltaEta);
                    strange_deltaPhiDistribution_lowerRatio->Fill(deltaPhi);
                    strange_deltaRDistribution_lowerRatio->Fill(deltaR);
                }
                else if(strangeRatio > 0.8 && strangeRatio < 1.2)
                {
                    strange_nConstDistribution_higherRatio->Fill(nConst);
                    strange_deltaEtaDistribution_higherRatio->Fill(deltaEta);
                    strange_deltaPhiDistribution_higherRatio->Fill(deltaPhi);
                    strange_deltaRDistribution_higherRatio->Fill(deltaR);
                }
            }
        }

        particles_fastjet.clear();
        jets.clear();
        jets_array->Clear();
        quarks->Clear();
    
    } // End of event loop equivalent
    
    //---------------------------------------------------------------------------------------------------------
    // Plotting histograms
    //---------------------------------------------------------------------------------------------------------

    
    // Configurações globais de estilo (coloque antes de renderizar os Canvas)
    //gStyle->SetOptStat(0);           // Desativa caixas de estatísticas para limpar o visual
    gStyle->SetTextFont(42);          // Fonte Helvetica estável
    gStyle->SetPadLeftMargin(0.13);   // Espaço para os títulos do eixo Y
    gStyle->SetPadRightMargin(0.05);  // Margem direita limpa
    gStyle->SetPadBottomMargin(0.12); // Espaço para títulos do eixo X
    gStyle->SetPadTopMargin(0.10);    // Espaço para os títulos principais dos plots
    gStyle->SetOptStat(0);

    TCanvas *c1 = new TCanvas("c1", "Strange jets histograms", 1800, 600);
    c1->Divide(2, 1);

    c1->cd(1);
    //strange_nConstDistribution_lowerRatio->SetTitle("N_{const} distribution of strange-tagged jet candidates;N_{const};Entries");
    strange_nConstDistribution_lowerRatio->SetTitle("");
    strange_nConstDistribution_lowerRatio->GetXaxis()->SetTitle("N_{const}");
    strange_nConstDistribution_lowerRatio->GetYaxis()->SetTitle("Entries");
    strange_nConstDistribution_lowerRatio->GetXaxis()->SetTitleSize(0.05);  // Eixo X
    strange_nConstDistribution_lowerRatio->GetYaxis()->SetTitleSize(0.05);  // Eixo Y

    strange_nConstDistribution_lowerRatio->GetYaxis()->SetRangeUser(std::min(strange_nConstDistribution_lowerRatio->GetMinimum(), strange_nConstDistribution_higherRatio->GetMinimum()) * 0.5, std::max(strange_nConstDistribution_lowerRatio->GetMaximum(), strange_nConstDistribution_higherRatio->GetMaximum()) * 1.05);
    strange_nConstDistribution_lowerRatio->SetLineColor(kBlue);
    strange_nConstDistribution_lowerRatio->SetLineWidth(2);
    strange_nConstDistribution_lowerRatio->Draw();
    strange_nConstDistribution_higherRatio->SetLineColor(kRed);
    strange_nConstDistribution_higherRatio->SetLineWidth(2);
    strange_nConstDistribution_higherRatio->Draw("same");
    
    /*
    TLegend *legend = new TLegend(0.6, 0.7, 0.89, 0.89);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.05);
    legend->AddEntry(strange_nConstDistribution_lowerRatio, "0.1 < p_{T}^{jet}/p_{T}^{quark} < 0.4", "l");
    legend->AddEntry(strange_nConstDistribution_higherRatio, "0.8 < p_{T}^{jet}/p_{T}^{quark} < 1.2", "l");
    legend->Draw();
    */

    /*
    c1->cd(2);
    strange_deltaEtaDistribution_lowerRatio->SetTitle("'Strange jets' #Delta#eta distribution for 0.1 < p_{T}^{jet}/p_{T}^{quark} < 0.4 (blue) and 0.8 < p_{T}^{jet}/p_{T}^{quark} < 1.2 (red);#Delta#eta;Entries");
    strange_deltaEtaDistribution_lowerRatio->GetYaxis()->SetRangeUser(std::min(strange_deltaEtaDistribution_lowerRatio->GetMinimum(), strange_deltaEtaDistribution_higherRatio->GetMinimum()) * 0.5, std::max(strange_deltaEtaDistribution_lowerRatio->GetMaximum(), strange_deltaEtaDistribution_higherRatio->GetMaximum()) * 1.05);
    strange_deltaEtaDistribution_lowerRatio->SetLineColor(kBlue);
    strange_deltaEtaDistribution_lowerRatio->Draw();
    strange_deltaEtaDistribution_higherRatio->SetLineColor(kRed);
    strange_deltaEtaDistribution_higherRatio->Draw("same");

    legend = new TLegend(0.6, 0.7, 0.89, 0.89);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.04);
    legend->AddEntry(strange_deltaEtaDistribution_lowerRatio, "0.1 < p_{T}^{jet}/p_{T}^{quark} < 0.4", "l");
    legend->AddEntry(strange_deltaEtaDistribution_higherRatio, "0.8 < p_{T}^{jet}/p_{T}^{quark} < 1.2", "l");
    legend->Draw();

    c1->cd(3);
    strange_deltaPhiDistribution_lowerRatio->SetTitle("'Strange jets' #Delta#phi distribution for 0.1 < p_{T}^{jet}/p_{T}^{quark} < 0.4 (blue) and 0.8 < p_{T}^{jet}/p_{T}^{quark} < 1.2 (red);#Delta#phi;Entries");
    strange_deltaPhiDistribution_lowerRatio->GetYaxis()->SetRangeUser(std::min(strange_deltaPhiDistribution_lowerRatio->GetMinimum(), strange_deltaPhiDistribution_higherRatio->GetMinimum()) * 0.5, std::max(strange_deltaPhiDistribution_lowerRatio->GetMaximum(), strange_deltaPhiDistribution_higherRatio->GetMaximum()) * 1.05);
    strange_deltaPhiDistribution_lowerRatio->SetLineColor(kBlue);
    strange_deltaPhiDistribution_lowerRatio->Draw();
    strange_deltaPhiDistribution_higherRatio->SetLineColor(kRed);
    strange_deltaPhiDistribution_higherRatio->Draw("same");

    legend = new TLegend(0.6, 0.7, 0.89, 0.89);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.04);
    legend->AddEntry(strange_deltaPhiDistribution_lowerRatio, "0.1 < p_{T}^{jet}/p_{T}^{quark} < 0.4", "l");
    legend->AddEntry(strange_deltaPhiDistribution_higherRatio, "0.8 < p_{T}^{jet}/p_{T}^{quark} < 1.2", "l");
    legend->Draw();
    */

    c1->cd(2);
    //strange_deltaRDistribution_lowerRatio->SetTitle("Strange-tagged jet candidates #DeltaR relative to the initial strange quark;#DeltaR;Entries");
    strange_deltaRDistribution_lowerRatio->SetTitle("");
    strange_deltaRDistribution_lowerRatio->GetYaxis()->SetTitle("");
    strange_deltaRDistribution_lowerRatio->GetXaxis()->SetTitle("#DeltaR");
    strange_deltaRDistribution_lowerRatio->GetXaxis()->SetTitleSize(0.05);  // Eixo X
    //strange_deltaRDistribution_lowerRatio->GetYaxis()->SetTitleSize(0.05);  // Eixo Y

    strange_deltaRDistribution_lowerRatio->GetYaxis()->SetRangeUser(std::min(strange_deltaRDistribution_lowerRatio->GetMinimum(), strange_deltaRDistribution_higherRatio->GetMinimum()) * 0.5, std::max(strange_deltaRDistribution_lowerRatio->GetMaximum(), strange_deltaRDistribution_higherRatio->GetMaximum()) * 1.05);
    strange_deltaRDistribution_lowerRatio->SetLineColor(kBlue);
    strange_deltaRDistribution_lowerRatio->SetLineWidth(2);
    strange_deltaRDistribution_lowerRatio->Draw();
    strange_deltaRDistribution_higherRatio->SetLineColor(kRed);
    strange_deltaRDistribution_higherRatio->SetLineWidth(2);
    strange_deltaRDistribution_higherRatio->Draw("same");


    TLegend *legend = new TLegend(0.6, 0.7, 0.89, 0.89);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.06);
    legend->AddEntry(strange_deltaRDistribution_lowerRatio, "0.1 < p_{T}^{jet}/p_{T}^{quark} < 0.4", "l");
    legend->AddEntry(strange_deltaRDistribution_higherRatio, "0.8 < p_{T}^{jet}/p_{T}^{quark} < 1.2", "l");
    legend->Draw();
    

    /* 
    TCanvas *c2 = new TCanvas("c2", "Charm jets histograms", 1200, 800);
    c2->Divide(2, 2);

    c2->cd(1);
    charm_nConstDistribution_lowerRatio->SetTitle("'Charm jets' N_{const} distribution for 0.1 < p_{T}^{jet}/p_{T}^{quark} < 0.4 (blue) and 0.8 < p_{T}^{jet}/p_{T}^{quark} < 1.2 (red);N_{const};Entries");
    charm_nConstDistribution_lowerRatio->GetYaxis()->SetRangeUser(std::min(charm_nConstDistribution_lowerRatio->GetMinimum(), charm_nConstDistribution_higherRatio->GetMinimum()) * 0.5, std::max(charm_nConstDistribution_lowerRatio->GetMaximum(), charm_nConstDistribution_higherRatio->GetMaximum()) * 1.05);
    charm_nConstDistribution_lowerRatio->SetLineColor(kBlue);
    charm_nConstDistribution_lowerRatio->Draw();
    charm_nConstDistribution_higherRatio->SetLineColor(kRed);
    charm_nConstDistribution_higherRatio->Draw("same");

    TLegend *legendC = new TLegend(0.6, 0.7, 0.89, 0.89);
    legendC->SetBorderSize(0);
    legendC->SetFillStyle(0);
    legendC->SetTextSize(0.04);
    legendC->AddEntry(charm_nConstDistribution_lowerRatio, "0.1 < p_{T}^{jet}/p_{T}^{quark} < 0.4", "l");
    legendC->AddEntry(charm_nConstDistribution_higherRatio, "0.8 < p_{T}^{jet}/p_{T}^{quark} < 1.2", "l");
    legendC->Draw();

    c2->cd(2);
    charm_deltaEtaDistribution_lowerRatio->SetTitle("'Charm jets' #Delta#eta distribution for 0.1 < p_{T}^{jet}/p_{T}^{quark} < 0.4 (blue) and 0.8 < p_{T}^{jet}/p_{T}^{quark} < 1.2 (red);#Delta#eta;Entries");
    charm_deltaEtaDistribution_lowerRatio->GetYaxis()->SetRangeUser(std::min(charm_deltaEtaDistribution_lowerRatio->GetMinimum(), charm_deltaEtaDistribution_higherRatio->GetMinimum()) * 0.5, std::max(charm_deltaEtaDistribution_lowerRatio->GetMaximum(), charm_deltaEtaDistribution_higherRatio->GetMaximum()) * 1.05);
    charm_deltaEtaDistribution_lowerRatio->SetLineColor(kBlue);
    charm_deltaEtaDistribution_lowerRatio->Draw();
    charm_deltaEtaDistribution_higherRatio->SetLineColor(kRed);
    charm_deltaEtaDistribution_higherRatio->Draw("same");

    legendC = new TLegend(0.6, 0.7, 0.89, 0.89);
    legendC->SetBorderSize(0);
    legendC->SetFillStyle(0);
    legendC->SetTextSize(0.04);
    legendC->AddEntry(charm_deltaEtaDistribution_lowerRatio, "0.1 < p_{T}^{jet}/p_{T}^{quark} < 0.4", "l");
    legendC->AddEntry(charm_deltaEtaDistribution_higherRatio, "0.8 < p_{T}^{jet}/p_{T}^{quark} < 1.2", "l");
    legendC->Draw();

    c2->cd(3);
    charm_deltaPhiDistribution_lowerRatio->SetTitle("'Charm jets' #Delta#phi distribution for 0.1 < p_{T}^{jet}/p_{T}^{quark} < 0.4 (blue) and 0.8 < p_{T}^{jet}/p_{T}^{quark} < 1.2 (red);#Delta#phi;Entries");
    charm_deltaPhiDistribution_lowerRatio->GetYaxis()->SetRangeUser(std::min(charm_deltaPhiDistribution_lowerRatio->GetMinimum(), charm_deltaPhiDistribution_higherRatio->GetMinimum()) * 0.5, std::max(charm_deltaPhiDistribution_lowerRatio->GetMaximum(), charm_deltaPhiDistribution_higherRatio->GetMaximum()) * 1.05);
    charm_deltaPhiDistribution_lowerRatio->SetLineColor(kBlue);
    charm_deltaPhiDistribution_lowerRatio->Draw();
    charm_deltaPhiDistribution_higherRatio->SetLineColor(kRed);
    charm_deltaPhiDistribution_higherRatio->Draw("same");

    legendC = new TLegend(0.6, 0.7, 0.89, 0.89);
    legendC->SetBorderSize(0);
    legendC->SetFillStyle(0);
    legendC->SetTextSize(0.04);
    legendC->AddEntry(charm_deltaPhiDistribution_lowerRatio, "0.1 < p_{T}^{jet}/p_{T}^{quark} < 0.4", "l");
    legendC->AddEntry(charm_deltaPhiDistribution_higherRatio, "0.8 < p_{T}^{jet}/p_{T}^{quark} < 1.2", "l");
    legendC->Draw();

    c2->cd(4);
    charm_deltaRDistribution_lowerRatio->SetTitle("'Charm jets' #DeltaR distribution for 0.1 < p_{T}^{jet}/p_{T}^{quark} < 0.4 (blue) and 0.8 < p_{T}^{jet}/p_{T}^{quark} < 1.2 (red);#DeltaR;Entries");
    charm_deltaRDistribution_lowerRatio->GetYaxis()->SetRangeUser(std::min(charm_deltaRDistribution_lowerRatio->GetMinimum(), charm_deltaRDistribution_higherRatio->GetMinimum()) * 0.5, std::max(charm_deltaRDistribution_lowerRatio->GetMaximum(), charm_deltaRDistribution_higherRatio->GetMaximum()) * 1.05);
    charm_deltaRDistribution_lowerRatio->SetLineColor(kBlue);
    charm_deltaRDistribution_lowerRatio->Draw();
    charm_deltaRDistribution_higherRatio->SetLineColor(kRed);
    charm_deltaRDistribution_higherRatio->Draw("same");

    legendC = new TLegend(0.6, 0.7, 0.89, 0.89);
    legendC->SetBorderSize(0);
    legendC->SetFillStyle(0);
    legendC->SetTextSize(0.04);
    legendC->AddEntry(charm_deltaRDistribution_lowerRatio, "0.1 < p_{T}^{jet}/p_{T}^{quark} < 0.4", "l");
    legendC->AddEntry(charm_deltaRDistribution_higherRatio, "0.8 < p_{T}^{jet}/p_{T}^{quark} < 1.2", "l");
    legendC->Draw();

    */


    TFile *second_outputFile = new TFile("modified2_histogramas_jetR_07_fullList.root", "RECREATE");
    strange_nConstDistribution_lowerRatio->Write();
    strange_nConstDistribution_higherRatio->Write();
    strange_deltaEtaDistribution_lowerRatio->Write();
    strange_deltaEtaDistribution_higherRatio->Write();
    strange_deltaPhiDistribution_lowerRatio->Write();
    strange_deltaPhiDistribution_higherRatio->Write();
    strange_deltaRDistribution_lowerRatio->Write();
    strange_deltaRDistribution_higherRatio->Write();
    second_outputFile->Close();

    file->Close();
}

