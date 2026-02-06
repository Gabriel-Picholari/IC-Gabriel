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


fastjet::PseudoJet quenchedJet(const fastjet::PseudoJet &jet, double deltaE)
{
    TLorentzVector jetVec(jet.px(), jet.py(), jet.pz(), jet.E());

    Double_t dPt = deltaE / TMath::CosH(jet.eta());
    TLorentzVector dEVec;
    dEVec.SetPtEtaPhiE(dPt, jet.eta(), jet.phi(), deltaE);
    TLorentzVector quenchedVec = jetVec - dEVec;
    fastjet::PseudoJet q(quenchedVec.Px(), quenchedVec.Py(), quenchedVec.Pz(), quenchedVec.E());

    return q;
}

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

void jetClassification1(const char* fileName)
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
    Int_t wPtFlag;
    Float_t PbNucleusRadius = 7; // fm
    
    TLorentzVector vec_s(0,0,0,0);
    TLorentzVector vec_c(0,0,0,0);

    TRandom3 generator(0);

    //---------------------------------------------------------------------------------------------------------
    // Initialization of histograms
    //---------------------------------------------------------------------------------------------------------

    TH1F *invariantMass = new TH1F("h1", "W^{+-} invariant mass spectrum [GeV/c^{2}]", 100, 0, 100);
    TH1F *primary_CharmRatioHist = new TH1F("primaryCharmRatio", "Charm: Ratio between tagged jets (with correct hadron) and quark p_{T}", 50, 0, 2);
    TH1F *secondary_CharmRatioHist = new TH1F("secondaryCharmRatio", "Charm: Ratio between tagged jets (whitout correct hadron) and quark p_{T}", 50, 0, 2);

    TH1F *primary_StrangeRatioHist = new TH1F("primaryStrangeRatio", "Strange: Ratio between tagged jets (with correct hadron) and quark p_{T}", 50, 0, 2);
    TH1F *secondary_StrangeRatioHist = new TH1F("secondaryStrangeRatio", "Strange: Ratio between tagged jets (whitout correct hadron) and quark p_{T}", 50, 0, 2);

    TH1F *missingStrangeConstituentsPdgMap = new TH1F("strangeMissingMap","Potential particles' PDGs missing in the strange list", 5000, -2500, 2500);
    TH1F *missingCharmConstituentsPdgMap = new TH1F("charmMissingMap","Potential particles' PDGs missing in the charm list", 5000, -2500, 2500);

    TH1F *observable_F_sc_Distribution_wBosonPt_greaterThan10 = new TH1F("observable_F_sc_Distribution_wBosonPt_greaterThan10", "Strange jet by charm jet p_{T} - F_{sc} distribution for events with W^{+-} boson p_{T} > 10", 100, 0, 10);    
    TH1F *observable_F_sc_Distribution_wBosonPt_smallerThan10 = new TH1F("observable_F_sc_Distribution_wBosonPt_smallerThan10", "Strange jet by charm jet p_{T} - F_{sc} distribution for events with W^{+-} boson p_{T} <= 10", 100, 0, 10);

    TH1F *observable_F_sc_Distribution_wBosonPt_greaterThan10_diffEnergyLosses = new TH1F("observable_F_sc_Distribution_wBosonPt_greaterThan10_diffEnergyLosses", "Strange jet by charm jet p_{T} - F_{sc} distribution for events with W^{+-} boson p_{T} > 10 for different energy loss funtions and null impact parameter", 100, 0, 10);    
    TH1F *observable_F_sc_Distribution_wBosonPt_smallerThan10_diffEnergyLosses = new TH1F("observable_F_sc_Distribution_wBosonPt_smallerThan10_diffEnergyLosses", "Strange jet by charm jet p_{T} - F_{sc} distribution for events with W^{+-} boson p_{T} <= 10 for different energy loss funtions and null impact parameter", 100, 0, 10);

    TH1F *observable_F_sc_Distribution_wBosonPt_greaterThan10_sameEnergyLosses = new TH1F("observable_F_sc_Distribution_wBosonPt_greaterThan10_sameEnergyLosses", "Strange jet by charm jet p_{T} - F_{sc} distribution for events with W^{+-} boson p_{T} > 10 for the same energy loss function and null impact parameter", 100, 0, 10);    
    TH1F *observable_F_sc_Distribution_wBosonPt_smallerThan10_sameEnergyLosses = new TH1F("observable_F_sc_Distribution_wBosonPt_smallerThan10_sameEnergyLosses", "Strange jet by charm jet p_{T} - F_{sc} distribution for events with W^{+-} boson p_{T} <= 10 for the same energy loss function and null impact parameter", 100, 0, 10);

    TH1F *second_observable_F_sc_Distribution_wBosonPt_greaterThan10_diffEnergyLosses = new TH1F("second_observable_F_sc_Distribution_wBosonPt_greaterThan10_diffEnergyLosses", "Strange jet by charm jet p_{T} - F_{sc} distribution for events with W^{+-} boson p_{T} > 10 for different energy loss funtions and non null impact parameter", 100, 0, 10);    
    TH1F *second_observable_F_sc_Distribution_wBosonPt_smallerThan10_diffEnergyLosses = new TH1F("second_observable_F_sc_Distribution_wBosonPt_smallerThan10_diffEnergyLosses", "Strange jet by charm jet p_{T} - F_{sc} distribution for events with W^{+-} boson p_{T} <= 10 for different energy loss funtions and non null impact parameter", 100, 0, 10);

    TH1F *second_observable_F_sc_Distribution_wBosonPt_greaterThan10_sameEnergyLosses = new TH1F("second_observable_F_sc_Distribution_wBosonPt_greaterThan10_sameEnergyLosses", "Strange jet by charm jet p_{T} - F_{sc} distribution for events with W^{+-} boson p_{T} > 10 for the same energy loss function and non null impact parameter", 100, 0, 10);    
    TH1F *second_observable_F_sc_Distribution_wBosonPt_smallerThan10_sameEnergyLosses = new TH1F("second_observable_F_sc_Distribution_wBosonPt_smallerThan10_sameEnergyLosses", "Strange jet by charm jet p_{T} - F_{sc} distribution for events with W^{+-} boson p_{T} <= 10 for the same energy loss function and non null impact parameter", 100, 0, 10);

    TH1F* h_rho_b0  = new TH1F("h_rho_b0",  "#rho (b=0);#rho [fm];Counts", 100, 0, 10);
    TH1F* h_phi_b0  = new TH1F("h_phi_b0",  "#phi (b=0);#phi;Counts",     100, -4, 4);

    TH1F* h_rho_bn0 = new TH1F("h_rho_bn0", "#rho (b#neq0);#rho [fm];Counts", 100, 0, 10);
    TH1F* h_phi_bn0 = new TH1F("h_phi_bn0", "#phi (b#neq0);#phi;Counts",      100, -4, 4);

    TH2F* xy_b0_jetSpawningCoordinates = new TH2F("xy_b0_jetSpawningCoordinates", "Jet Spawning Coordinates (b=0)", 100, -10, 10, 100, -10, 10);
    TH2F* xy_b_jetSpawningCoordinates = new TH2F("xy_b_jetSpawningCoordinates", "Jet Spawning Coordinates (b#neq0)", 100, -10, 10, 100, -10, 10);

    TH1F* pT_charmDistribution_beforeQuenching = new TH1F("pT_charmDistribution_beforeQuenching", "Charm Jet p_{T} Distribution before Quenching; p_{T} [GeV/c]; Frequency", 100, 0, 100);
    TH1F* pT_strangeDistribution_beforeQuenching = new TH1F("pT_strangeDistribution_beforeQuenching", "Strange Jet p_{T} Distribution before Quenching; p_{T} [GeV/c]; Frequency", 100, 0, 100);
    

    TH1F* null_b_pT_charmDistribution_afterQuenching_diff = new TH1F("null_b_pT_charmDistribution_afterQuenching_diff", "Charm Jet p_{T} Distribution after Quenching with Different Energy Losses (b=0); p_{T} [GeV/c]; Frequency", 100, 0, 100);
    TH1F* null_b_pT_strangeDistribution_afterQuenching_diff = new TH1F("null_b_pT_strangeDistribution_afterQuenching_diff", "Strange Jet p_{T} Distribution after Quenching with Different Energy Losses (b=0); p_{T} [GeV/c]; Frequency", 100, 0, 100);

    TH1F* non_null_b_pT_charmDistribution_afterQuenching_diff = new TH1F("non_null_b_pT_charmDistribution_afterQuenching_diff", "Charm Jet p_{T} Distribution after Quenching Different Energy Losses (b#neq0); p_{T} [GeV/c]; Frequency", 100, 0, 100);
    TH1F* non_null_b_pT_strangeDistribution_afterQuenching_diff = new TH1F("non_null_b_pT_strangeDistribution_afterQuenching_diff", "Strange Jet p_{T} Distribution after Quenching Different Energy Losses (b#neq0); p_{T} [GeV/c]; Frequency", 100, 0, 100);

    auto normalize = [](TH1F* h) {
        Double_t integral = h->Integral();
        if (integral > 0) h->Scale(1.0 / integral);
    };

    auto path_length_in_circle = [&](Double_t xc, Double_t yc, Double_t x0, Double_t y0, Double_t phi, Double_t R) -> Double_t
    {
        Double_t dx = x0 - xc;
        Double_t dy = y0 - yc;

        Double_t A = dx * std::cos(phi) + dy * std::sin(phi);
        Double_t B = dx*dx + dy*dy - R*R;

        Double_t Delta = A*A - B;
        if (Delta < 0.0) 
        {
            return 0.0;
        }

        Double_t t_exit = -A + std::sqrt(Delta);
        return t_exit;
    };


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
            //std::cout << strangePt << std::endl;                  // Check over that later (not a serious problem I suppose)

        }

        //---------------------------------------------------------------------------------------------------------
        // Particle loop equivalent
        //---------------------------------------------------------------------------------------------------------

        for (Int_t nj = 0; nj < jets_array->GetEntries(); nj++)
        {
            MyJet *fp = static_cast<MyJet *>(jets_array->At(nj));

            
            wPtFlag = fp->wPtFlag; // Again, even though it is rewritten for each (final) particle kept in the TTree, it is the same for all particles within the event nj
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
            if (jetPt < 5) continue; // Basic cut on jet pT

            jetEta = jet.eta();

            Float_t absEta = TMath::Abs(jetEta);
            if (absEta < 1) continue; // Basic cut on jet eta
            
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

            TLorentzVector currentJet(jetPx, jetPy, jetPz, jetE);
            if (currentJet.M() < 0) continue;
            
            Bool_t isCharmTagged = false;
            Bool_t isStrangeTagged = false;

            for (const fastjet::PseudoJet &constituent : jet.constituents())
            {
                TString signalType_jet = constituent.user_info<JetInfo>().getSignalType();

                if (signalType_jet == "charm" && !isCharmTagged) // Then one of the final particles that is part of this jet has a distant mother in a c quark
                {   
                    tagged_c_jets.push_back(jet);
                    vec_c = TLorentzVector(jetPx, jetPy, jetPz, jetE); // Just to check consistency in the number of entries (will be reconsidered)
                    //std::cout << "Added charm jet with mass: " << currentJet.M() << std::endl;
                    isCharmTagged = true;
                }
                else if (signalType_jet == "strange" && !isStrangeTagged) // Then it's a strange jet
                {
                    tagged_s_jets.push_back(jet);
                    vec_s = TLorentzVector(jetPx, jetPy, jetPz, jetE); // Just to check consistency in the number of entries (will be reconsidered)
                    //std::cout << "Added strange jet with mass: " << currentJet.M() << std::endl;
                    isStrangeTagged = true;
                }
            }

        } // End of individual jet creation

        const std::unordered_set<int> charmPdgSet = {411, 421, 413, 423, 415, 425, 431, 433, 435, 10411, 10421, 413, 423, 10413, 10423, 20413, 20423, 415, 425, 431, 10431, 433, 10433, 20433, 435, 4122, 4222, 4212, 4112, 4224, 4214, 4114, 4232, 4132, 4322, 4312, 4324, 4314, 4332, 4334, 4412, 4422, 4414, 4424, 4432, 4434, 4444};
        const std::unordered_set<int> strangePdgSet = {130, 310, 311, 321, 313, 323, 315, 325, 317, 327, 319, 329, 9000311, 9000321, 10311, 10321, 100311, 100321, 9010311, 9010321, 9020311, 9020321, 313, 323, 10313, 10323, 20313, 20323, 100313, 100323, 9000313, 9000323, 30313, 30323, 315, 325, 9000315, 9000325, 10315, 10325, 20315, 20325, 9010315, 9010325, 9020315, 9020325, 317, 327, 9010317, 9010327, 319, 329, 3122, 3222, 3212, 3112, 3224, 3214, 3114, 3322, 3312, 3324, 3314, 3334};
        
        std::vector<fastjet::PseudoJet> good_c_jets;
        std::vector<fastjet::PseudoJet> good_s_jets;

        //std::cout << "\n--- c-tagged Jets ---" << std::endl;
        for (const fastjet::PseudoJet  &jet : tagged_c_jets) // Opening the jet vector: the analysis object is a jet
        {
            TLorentzVector cJet(jet.px(), jet.py(), jet.pz(), jet.E());
            Float_t charmRatio = jet.pt() / charmPt;

            //std::cout << "----- New Jet -----" << std::endl;
            //std::cout << "Jet Mass: " << cJet.M() << std::endl;
            //std::cout << "Jet pT: " << cJet.Pt() << std::endl;
            //std::cout << "Charm ratio: " << charmRatio << std::endl;
            //std::cout << std::endl;
            

            Bool_t hasCharmConstituent = false;

            for (const fastjet::PseudoJet &constituent : jet.constituents()) // Opening the jet itself: the analysis object is a constituent of the jet
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

            if (hasCharmConstituent && charmRatio > 0.60) // Please, refer to 5th september research log for the reasoning behind the 0.6 lower limit
            {
                good_c_jets.push_back(jet);
            }

            if (hasCharmConstituent) 
            {
                primary_CharmRatioHist->Fill(charmRatio);
            }
            else
            {
                secondary_CharmRatioHist->Fill(charmRatio);

                if (charmRatio < 0.85) continue;
                for (const fastjet::PseudoJet &constituent : jet.constituents())
                {
                    // Since we're only interested in the jets that "went missing" with more than, let' say, 75% of the quark's pT, we can place our attention to jets that go above this number
                    Int_t constituentPdg = constituent.user_info<JetInfo>().getFinalParticlePdg();
                    //std::cout << "Missing jet constituents: " << constituentPdg << std::endl;
                    missingCharmConstituentsPdgMap->Fill(constituentPdg);
                }
            }
        }

        //std::cout << "\n--- s-tagged Jets ---" << std::endl;
        for (const fastjet::PseudoJet &jet : tagged_s_jets)
        {
            TLorentzVector sJet(jet.px(), jet.py(), jet.pz(), jet.E());
            Float_t strangeRatio = jet.pt() / strangePt;

            //std::cout << "----- New Jet -----" << std::endl;
            //std::cout << "Jet Mass: " << sJet.M() << std::endl;
            //std::cout << "Jet pT: " << sJet.Pt() << std::endl;
            //std::cout << "Stange ratio: " << strangeRatio << std::endl;
            //std::cout << std::endl;

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

            if (hasStrangeConstituent && strangeRatio > 0.60) // Same as before
            {
                good_s_jets.push_back(jet);
            }

            if (hasStrangeConstituent)
            {
                primary_StrangeRatioHist->Fill(strangeRatio);
            }
            else
            {
                secondary_StrangeRatioHist->Fill(strangeRatio);

                // If the event reaches this point, it means we have a jet that does not have any hadron within our list. We are, then, interested in the constituents related to this jet in particular
                // Since we're only interested in the jets that "went missing" with more than, let' say, 75% of the quark's pT, we can place our attention to jets that go above this number
                // Thus, I guess we can just reopen the jet into it's constituents and get them directly from there

                if (strangeRatio < 0.85) continue;
                for (const fastjet::PseudoJet &constituent : jet.constituents())
                {
                    // Since we're only interested in the jets that "went missing" with more than, let' say, 75% of the quark's pT, we can place our attention to jets that go above this number
                    Int_t constituentPdg = constituent.user_info<JetInfo>().getFinalParticlePdg();
                    //std::cout << "Missing jet constituents: " << constituentPdg << std::endl;
                    missingStrangeConstituentsPdgMap->Fill(constituentPdg);
                }
            }
        }

        // Combinatoral loop to buid the W boson invariant mass scpectrum

        for (size_t ic = 0; ic < good_c_jets.size(); ++ic) 
        {
            TLorentzVector vc(good_c_jets[ic].px(), good_c_jets[ic].py(), good_c_jets[ic].pz(), good_c_jets[ic].E());

            for (size_t is = 0; is < good_s_jets.size(); ++is) 
            {
                TLorentzVector vs(good_s_jets[is].px(), good_s_jets[is].py(), good_s_jets[is].pz(), good_s_jets[is].E());
                TLorentzVector vW = vc + vs;
                invariantMass->Fill(vW.M());
            }
        } // Notice how it was a very straight forward implementation: all jets of charm were combined with all jets of strange in the same event

        //---------------------------------------------------------------------------------------------------------
        // F_sc observable calculation
        //---------------------------------------------------------------------------------------------------------

        Float_t backToback_lowerLimit = 7*TMath::Pi()/8;
        Float_t backToback_upperLimit = 11*TMath::Pi()/8;

        fastjet::PseudoJet event_strange_jet;
        fastjet::PseudoJet event_charmed_jet;
        Double_t max_pt_c = -1.0; // We use a negative value so the fist jet in the vector will always be the maximum pT jet
        Double_t max_pt_s = -1.0;

        for (const fastjet::PseudoJet &jet : good_c_jets)
        {
            if (jet.pt() > max_pt_c)
            {
                max_pt_c = jet.pt();
                event_charmed_jet = jet;
            }
        }
        
        for (const fastjet::PseudoJet &jet : good_s_jets)
        {
            if (jet.pt() > max_pt_s)
            {
                max_pt_s = jet.pt();
                event_strange_jet = jet;
            }
        }

        // At this point, we have the leading jets for both flavors well defined

        Float_t deltaPhi = TMath::Abs(event_charmed_jet.phi() - event_strange_jet.phi());

        if (deltaPhi >= backToback_lowerLimit && deltaPhi < backToback_upperLimit) // Ensuring the back-to-back condition is fulfilled for all jets from now on
        {
            //---------------------------------------------------------------------------------------------------------
            // No energy loss block
            //---------------------------------------------------------------------------------------------------------
            
            Float_t F_sc = event_strange_jet.pt() / event_charmed_jet.pt(); // New observable - for more information, check November 6th, 2025 research log

            pT_charmDistribution_beforeQuenching->Fill(event_charmed_jet.pt());
            pT_strangeDistribution_beforeQuenching->Fill(event_strange_jet.pt());

            if (wPtFlag == 1) // Then the W boson pT is greater than 10 GeV/c
            {
                observable_F_sc_Distribution_wBosonPt_greaterThan10->Fill(F_sc);
            }
            else // Then the W boson pT is smaller than or equal to 10 GeV/c
            {
                observable_F_sc_Distribution_wBosonPt_smallerThan10->Fill(F_sc);
            }
            
            // The plot above is the one without the QGP energy loss consideration
            // Now, we shall procede to calculate the entries of a plot which contemplates the energy loss due to the medium created in an artificial Pb-Pb collision

            //---------------------------------------------------------------------------------------------------------
            // Ignorable block - outdated rho and phi generation for (b=0)
            //---------------------------------------------------------------------------------------------------------

            // For more information on the equations implemented bellow, check on January 13th, 2026 research log

            //Float_t U = generator.Uniform(0, 1);
            //Float_t V = generator.Uniform(0, 1);
            //Float_t randomPhi = 2*TMath::Pi()*V;
            //Float_t randomRho = PbNucleusRadius * TMath::Sqrt(1 - TMath::Power(1 - U, 2/3));

            //---------------------------------------------------------------------------------------------------------
            // General rho and phi generation method breakdown:
            //---------------------------------------------------------------------------------------------------------

            // Since the method above is analytical and limited when the interest raises to non null impact parameter energy loss regions, I will implement a more general one whose steps are explained below:

            // First step:  we define a region D of energy loss given by the overlapping area of two circles of radius R (Pb nucleus radius) separated by a distance b (impact parameter). Notice that
            // this region will explicitly depend on the impact parameter, thus making it usable for both central and peripheral collisions

            // Second step: we uniformly generate points (x,y) within a box going from -R to R and check if they belong to the region D. If positive, we keep the point, if negative, 
            // we discard it and generate a new one 

            // Third step: we apply a weight function, w(x,y), to the points kept. The functional form of this function is yet to be defined. In any case, we find the maximum of w(r), w_max, which
            // will be used as reference point. In this step we generate another uniform random number u_4 and impose the condition u_4 < w(r) / w_max to finally accept or reject the point.

            // Fourth step: with a valid point in hands, we project it to the xy plane and calculate the rho and phi as usual to be used in further calculations

            // Note: In first intance, we will only work with x and y coordinates

            //---------------------------------------------------------------------------------------------------------

            // We take the jets we own and place them into a new polar coordinate system. The jets are then placed in (randomRho, randomPhi).
            // Since by hypotesis we'll take them as back to back, this point will be the starting point for both of them in this new system.
            // Our intentions are to calculate the path lenght each jet would have within a Pb nucleus if they were created at this random point and traveled in opposite directions.
            // We can account for two situations: either the impact parameter is zero (central collision) or it has a finite value (peripheral collision). For the time being, we'll take b = 0.

            //---------------------------------------------------------------------------------------------------------
            // Energy loss block (b=0)
            //---------------------------------------------------------------------------------------------------------

            Float_t randomRho, randomPhi;
            Float_t xMain, yMain;
            
            while(true)
            {
                Float_t u1 = generator.Uniform(0, 1);
                Float_t u2 = generator.Uniform(0, 1);

                Float_t x, y;

                x = - PbNucleusRadius + (2 * PbNucleusRadius) * u1;
                y = - PbNucleusRadius + (2 * PbNucleusRadius) * u2;
                // This way, we're mapping U in X or U in [c, d] to X in [a, b] by X = a + (b-a)*U, making sure X is uniform as well as U

                if ( x*x + y*y > PbNucleusRadius*PbNucleusRadius) continue; // In the b=0 case, only one inequality will make it

                Float_t w = 1;
                Float_t w_max = 1;
                Float_t u4 = generator.Uniform(0, 1);

                if (u4 > w / w_max) continue; // Trivially satisfied and shown here only for consistency with the non null impact parameter case we shall get in detail next

                randomRho = std::sqrt(x*x + y*y);
                randomPhi   = std::atan2(y, x);
                xMain = x;
                yMain = y;

                break;
            }
            
            h_rho_b0->Fill(randomRho);
            h_phi_b0->Fill(randomPhi);
            xy_b0_jetSpawningCoordinates->Fill(xMain, yMain);

            // Charm jet path length calculation
            Float_t phi_tot_charm = randomPhi - event_charmed_jet.phi();
            Float_t path_length_charm;
            Float_t path_length_charm_1 = -randomRho * TMath::Cos(phi_tot_charm) + TMath::Sqrt( TMath::Power(PbNucleusRadius, 2) - TMath::Power(randomRho * TMath::Sin(phi_tot_charm), 2) );
            Float_t path_length_charm_2 = -randomRho * TMath::Cos(phi_tot_charm) - TMath::Sqrt( TMath::Power(PbNucleusRadius, 2) - TMath::Power(randomRho * TMath::Sin(phi_tot_charm), 2) );

            if (path_length_charm_1 > 0)
            {
                path_length_charm = path_length_charm_1;
            }
            else
            {
                path_length_charm = path_length_charm_2;
            }

            // Strange jet path length calculation
            Float_t phi_tot_strange = randomPhi - event_strange_jet.phi();
            Float_t path_length_strange;
            Float_t path_length_strange_1 = -randomRho * TMath::Cos(phi_tot_strange) + TMath::Sqrt( TMath::Power(PbNucleusRadius, 2) - TMath::Power(randomRho * TMath::Sin(phi_tot_strange), 2) );
            Float_t path_length_strange_2 = -randomRho * TMath::Cos(phi_tot_strange) - TMath::Sqrt( TMath::Power(PbNucleusRadius, 2) - TMath::Power(randomRho * TMath::Sin(phi_tot_strange), 2) );

            if (path_length_strange_1 > 0)
            {
                path_length_strange = path_length_strange_1;
            }
            else
            {
                path_length_strange = path_length_strange_2;
            }

            // Applying a simplified energy loss model to check consistency
            Float_t dE_dx_charm = 2; // GeV/fm
            Float_t dE_dx_strange = 8; // GeV/fm
            Float_t dE_dx_both = 4; // GeV/fm

            Float_t deltaE_charm = dE_dx_charm * path_length_charm;
            Float_t deltaE_strange = dE_dx_strange * path_length_strange;

            Float_t deltaE_charm_both = dE_dx_both * path_length_charm;
            Float_t deltaE_strange_both = dE_dx_both * path_length_strange;

            //With the energy shift calculated, we can now redefine the jet 4-momenta accordingly using the quenchedJet function defined at the beginning of this code
            fastjet::PseudoJet quenched_charm_jet = quenchedJet(event_charmed_jet, deltaE_charm);
            fastjet::PseudoJet quenched_strange_jet = quenchedJet(event_strange_jet, deltaE_strange);
            fastjet::PseudoJet quenched_charm_jet_both = quenchedJet(event_charmed_jet, deltaE_charm_both);
            fastjet::PseudoJet quenched_strange_jet_both = quenchedJet(event_strange_jet, deltaE_strange_both);

            null_b_pT_charmDistribution_afterQuenching_diff->Fill(quenched_charm_jet.pt());
            null_b_pT_strangeDistribution_afterQuenching_diff->Fill(quenched_strange_jet.pt());

            // Now we can recalculate F_sc with the quenched jets
            Float_t F_sc_quenched = quenched_strange_jet.pt() / quenched_charm_jet.pt();
            Float_t F_sc_quenched_both = quenched_strange_jet_both.pt() / quenched_charm_jet_both.pt();

            if (wPtFlag == 1) // Then the W boson pT is greater than 10 GeV/c
            {
                observable_F_sc_Distribution_wBosonPt_greaterThan10_diffEnergyLosses->Fill(F_sc_quenched);
                observable_F_sc_Distribution_wBosonPt_greaterThan10_sameEnergyLosses->Fill(F_sc_quenched_both);
            }
            else // Then the W boson pT is smaller than or equal to 10 GeV/c
            {
                observable_F_sc_Distribution_wBosonPt_smallerThan10_diffEnergyLosses->Fill(F_sc_quenched);
                observable_F_sc_Distribution_wBosonPt_smallerThan10_sameEnergyLosses->Fill(F_sc_quenched_both);
            }

            //---------------------------------------------------------------------------------------------------------
            // Energy loss block (b different from 0)
            //---------------------------------------------------------------------------------------------------------

            Float_t b = 9; // Since the radius is approximately 7 fm, b ranges from 0 to 14 fm
            Float_t x0, y0;

            while(true)
            {
                Float_t u1 = generator.Uniform(0, 1);
                Float_t u2 = generator.Uniform(0, 1);

                Float_t x, y;

                x = - PbNucleusRadius + (2 * PbNucleusRadius) * u1;
                y = - PbNucleusRadius + (2 * PbNucleusRadius) * u2;
                // This way, we're mapping U in X or U in [c, d] to X in [a, b] by X = a + (b-a)*U, making sure X is uniform as well as U

                // Inside the sphere A with center in x = -b/2
                Float_t dxA = x + 0.5f*b;
                bool insideA = (dxA*dxA + y*y <= PbNucleusRadius*PbNucleusRadius);

                // Inside the sphere B with center in x = +b/2
                Float_t dxB = x - 0.5f*b;
                bool insideB = (dxB*dxB + y*y <= PbNucleusRadius*PbNucleusRadius);

                if (!(insideA && insideB)) continue;

                Float_t w = 1;
                Float_t w_max = 1;
                Float_t u4 = generator.Uniform(0, 1);

                if (u4 > w / w_max) continue; // Trivially satisfied and shown here only for consistency with the non null impact parameter case we shall get in detail next

                x0 = x;
                y0 = y;

                break;
            }

            Float_t rho_bn0 = std::sqrt(x0*x0 + y0*y0);
            Float_t phi_bn0 = std::atan2(y0, x0);

            h_rho_bn0->Fill(rho_bn0);
            h_phi_bn0->Fill(phi_bn0);
            xy_b_jetSpawningCoordinates->Fill(x0, y0);

            /* Previous method for (x0, y0) generation - outdated

            Double_t bmax = 2.0 * PbNucleusRadius;
            Double_t u = generator.Uniform(0.0, 1.0);
            Double_t b = bmax * std::sqrt(u);
            
            Double_t x0, y0;

            while (true) 
            {
                // 1) Initial guess for (x0, y0)
                x0 = generator.Uniform(-PbNucleusRadius, PbNucleusRadius);
                y0 = generator.Uniform(-PbNucleusRadius, PbNucleusRadius);

                // 2) Check if the point is inside both nuclei (overlap region)
                Double_t dxA = x0 + 0.5 * b;
                Double_t dyA = y0;
                Double_t dxB = x0 - 0.5 * b;
                Double_t dyB = y0;

                bool insideA = (dxA*dxA + dyA*dyA <= PbNucleusRadius*PbNucleusRadius);
                bool insideB = (dxB*dxB + dyB*dyB <= PbNucleusRadius*PbNucleusRadius);

                if (insideA && insideB) break; // Valid point encontered
            }
            */

            Double_t xcA = -0.5 * b;
            Double_t ycA = 0.0;
            Double_t xcB = +0.5 * b;
            Double_t ycB = 0.0;

            Double_t L_A_charm = path_length_in_circle(xcA, ycA, x0, y0, event_charmed_jet.phi(), PbNucleusRadius);
            Double_t L_B_charm = path_length_in_circle(xcB, ycB, x0, y0, event_charmed_jet.phi(), PbNucleusRadius);
            Double_t L_charm = std::min(L_A_charm, L_B_charm);

            Double_t L_A_strange = path_length_in_circle(xcA, ycA, x0, y0, event_strange_jet.phi(), PbNucleusRadius);
            Double_t L_B_strange = path_length_in_circle(xcB, ycB, x0, y0, event_strange_jet.phi(), PbNucleusRadius);
            Double_t L_strange = std::min(L_A_strange, L_B_strange);

            Float_t second_deltaE_charm = dE_dx_charm * L_charm;
            Float_t second_deltaE_strange = dE_dx_strange * L_strange;

            Float_t second_deltaE_charm_both = dE_dx_both * L_charm;
            Float_t second_deltaE_strange_both = dE_dx_both * L_strange;

            fastjet::PseudoJet second_quenched_charm_jet = quenchedJet(event_charmed_jet, second_deltaE_charm);
            fastjet::PseudoJet second_quenched_strange_jet = quenchedJet(event_strange_jet, second_deltaE_strange);
            fastjet::PseudoJet second_quenched_charm_jet_both = quenchedJet(event_charmed_jet, second_deltaE_charm_both);
            fastjet::PseudoJet second_quenched_strange_jet_both = quenchedJet(event_strange_jet, second_deltaE_strange_both);

            non_null_b_pT_charmDistribution_afterQuenching_diff->Fill(second_quenched_charm_jet.pt());
            non_null_b_pT_strangeDistribution_afterQuenching_diff->Fill(second_quenched_strange_jet.pt());

            // Now we can recalculate F_sc with the quenched jets
            Float_t second_F_sc_quenched = second_quenched_strange_jet.pt() / second_quenched_charm_jet.pt();
            Float_t second_F_sc_quenched_both = second_quenched_strange_jet_both.pt() / second_quenched_charm_jet_both.pt();

            if (wPtFlag == 1) // Then the W boson pT is greater than 10 GeV/c
            {
                second_observable_F_sc_Distribution_wBosonPt_greaterThan10_diffEnergyLosses->Fill(second_F_sc_quenched);
                second_observable_F_sc_Distribution_wBosonPt_greaterThan10_sameEnergyLosses->Fill(second_F_sc_quenched_both);
            }
            else // Then the W boson pT is smaller than or equal to 10 GeV/c
            {
                second_observable_F_sc_Distribution_wBosonPt_smallerThan10_diffEnergyLosses->Fill(second_F_sc_quenched);
                second_observable_F_sc_Distribution_wBosonPt_smallerThan10_sameEnergyLosses->Fill(second_F_sc_quenched_both);
            }

        }

        particles_fastjet.clear();
        jets.clear();
        jets_array->Clear();
        quarks->Clear();
    
    } // End of event loop equivalent


    TH1F* observable_F_sc_g10_norm = (TH1F*) observable_F_sc_Distribution_wBosonPt_greaterThan10->Clone("observable_F_sc_g10_norm");
    TH1F* observable_F_sc_l10_norm = (TH1F*) observable_F_sc_Distribution_wBosonPt_smallerThan10->Clone("observable_F_sc_l10_norm");
    TH1F* observable_F_sc_g10_diff_norm = (TH1F*) observable_F_sc_Distribution_wBosonPt_greaterThan10_diffEnergyLosses->Clone("observable_F_sc_g10_diff_norm");
    TH1F* observable_F_sc_l10_diff_norm = (TH1F*) observable_F_sc_Distribution_wBosonPt_smallerThan10_diffEnergyLosses->Clone("observable_F_sc_l10_diff_norm");

    TH1F* observable_F_sc_g10_same_norm = (TH1F*) observable_F_sc_Distribution_wBosonPt_greaterThan10_sameEnergyLosses->Clone("observable_F_sc_g10_same_norm");
    TH1F* observable_F_sc_l10_same_norm = (TH1F*) observable_F_sc_Distribution_wBosonPt_smallerThan10_sameEnergyLosses->Clone("observable_F_sc_l10_same_norm");

    // Normalize the copies
    normalize(observable_F_sc_g10_norm);
    normalize(observable_F_sc_l10_norm);
    normalize(observable_F_sc_g10_diff_norm);
    normalize(observable_F_sc_l10_diff_norm);
    normalize(observable_F_sc_g10_same_norm);
    normalize(observable_F_sc_l10_same_norm);

    TH1F* ratio_greater10_diff = (TH1F*) observable_F_sc_g10_diff_norm->Clone("ratio_greater10_diff");
    ratio_greater10_diff->Divide(observable_F_sc_g10_norm);

    TH1F* ratio_smaller10_diff = (TH1F*) observable_F_sc_l10_diff_norm->Clone("ratio_smaller10_diff");
    ratio_smaller10_diff->Divide(observable_F_sc_l10_norm);

    TH1F* ratio_greater10_same = (TH1F*) observable_F_sc_g10_same_norm->Clone("ratio_greater10_same");
    ratio_greater10_same->Divide(observable_F_sc_g10_norm);

    TH1F* ratio_smaller10_same = (TH1F*) observable_F_sc_l10_same_norm->Clone("ratio_smaller10_same");
    ratio_smaller10_same->Divide(observable_F_sc_l10_norm);

    // ---------------------------------------------------------------------------------------------------------

    TH1F* second_F_sc_g10_diff_norm = (TH1F*) second_observable_F_sc_Distribution_wBosonPt_greaterThan10_diffEnergyLosses->Clone("second_F_sc_g10_diff_norm");
    TH1F* second_F_sc_l10_diff_norm = (TH1F*) second_observable_F_sc_Distribution_wBosonPt_smallerThan10_diffEnergyLosses->Clone("second_F_sc_l10_diff_norm");
    TH1F* second_F_sc_g10_same_norm = (TH1F*) second_observable_F_sc_Distribution_wBosonPt_greaterThan10_sameEnergyLosses->Clone("second_F_sc_g10_same_norm");
    TH1F* second_F_sc_l10_same_norm = (TH1F*) second_observable_F_sc_Distribution_wBosonPt_smallerThan10_sameEnergyLosses->Clone("second_F_sc_l10_same_norm");

    normalize(second_F_sc_g10_diff_norm);
    normalize(second_F_sc_l10_diff_norm);
    normalize(second_F_sc_g10_same_norm);
    normalize(second_F_sc_l10_same_norm);

    TH1F* second_ratio_greater10_diff = (TH1F*) second_F_sc_g10_diff_norm->Clone("second_ratio_greater10_diff");
    second_ratio_greater10_diff->Divide(observable_F_sc_g10_norm);

    TH1F* second_ratio_smaller10_diff = (TH1F*) second_F_sc_l10_diff_norm->Clone("second_ratio_smaller10_diff");
    second_ratio_smaller10_diff->Divide(observable_F_sc_l10_norm);

    TH1F* second_ratio_greater10_same = (TH1F*) second_F_sc_g10_same_norm->Clone("second_ratio_greater10_same");
    second_ratio_greater10_same->Divide(observable_F_sc_g10_norm);

    TH1F* second_ratio_smaller10_same = (TH1F*) second_F_sc_l10_same_norm->Clone("second_ratio_smaller10_same");
    second_ratio_smaller10_same->Divide(observable_F_sc_l10_norm);
    
    //---------------------------------------------------------------------------------------------------------
    // Plotting histograms
    //---------------------------------------------------------------------------------------------------------

    TCanvas *c1 = new TCanvas("c1", "Invariant mass distribution", 2500, 2500);
    c1->Divide(1, 1);

    c1->cd(1);
    invariantMass->SetTitle("Boson W^{#pm} invariant mass spectrum");
    invariantMass->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    invariantMass->GetYaxis()->SetTitle("Frequency");
    invariantMass->DrawCopy();

    TCanvas *c2 = new TCanvas("c2", "Charm", 2500, 2500);
    c2->Divide(1, 1);

    c2->cd(1);
    //secondary_CharmRatioHist->SetTitle("Distribution of charm jet-to-quark p_{T} ratio with no p_{T} cut applied"); // (Title is correct when no cut is wanted!)
    primary_CharmRatioHist->SetTitle("Distribution of charm jet-to-quark p_{T} ratio with a 5 GeV/c jet p_{T} cut"); // (Title is correct when a cut is wanted!)

    primary_CharmRatioHist->GetXaxis()->SetTitle("Ratio");
    primary_CharmRatioHist->GetYaxis()->SetTitle("Frequency");
    primary_CharmRatioHist->SetLineColor(kGreen);
    primary_CharmRatioHist->DrawCopy();

    primary_CharmRatioHist->SetLineColor(kRed);
    secondary_CharmRatioHist->DrawCopy("same");

    TCanvas *c3 = new TCanvas("c3", "Strange", 2500, 2500);
    c3->Divide(1, 1);

    c3->cd(1);
    //secondary_StrangeRatioHist->SetTitle("Distribution of strange jet-to-quark p_{T} ratio with no p_{T} cut applied"); // (Title is correct when no cut is wanted!)

    primary_StrangeRatioHist->SetTitle("Distribution of strange jet-to-quark p_{T} ratio with a 5 GeV/c jet p_{T} cut"); // (Title is correct when a cut is wanted!)
    
    primary_StrangeRatioHist->GetXaxis()->SetTitle("Ratio");
    primary_StrangeRatioHist->GetYaxis()->SetTitle("Frequency");
    primary_StrangeRatioHist->SetLineColor(kGreen);
    primary_StrangeRatioHist->DrawCopy();
    
    secondary_StrangeRatioHist->SetLineColor(kRed);
    secondary_StrangeRatioHist->DrawCopy("same");

    TCanvas *c4 = new TCanvas("c4", "Missing particles PDG", 2500, 2500);
    c4->Divide(1, 2);

    c4->cd(1);
    missingStrangeConstituentsPdgMap->SetTitle("Potentially missing particles PDGs for strange jets list");
    missingStrangeConstituentsPdgMap->GetXaxis()->SetTitle("PDG");
    missingStrangeConstituentsPdgMap->GetYaxis()->SetTitle("Frequency");
    missingStrangeConstituentsPdgMap->DrawCopy();

    c4->cd(2);
    missingCharmConstituentsPdgMap->SetTitle("Potentially missing particles PDGs for charm jets list");
    missingCharmConstituentsPdgMap->GetXaxis()->SetTitle("PDG");
    missingCharmConstituentsPdgMap->GetYaxis()->SetTitle("Frequency");
    missingCharmConstituentsPdgMap->DrawCopy();

    TCanvas *c5 = new TCanvas("c5", "Observable F_{sc} distributions no energy losses", 2500, 2500);
    c5->Divide(1, 2);

    c5->cd(1);
    observable_F_sc_Distribution_wBosonPt_greaterThan10->SetTitle("Observable F_{sc} distribution for events with W^{+-} boson p_{T} > 10");
    observable_F_sc_Distribution_wBosonPt_greaterThan10->GetXaxis()->SetTitle("p_{T s}^{jet} / p_{T c}^{jet}");
    observable_F_sc_Distribution_wBosonPt_greaterThan10->GetYaxis()->SetTitle("Frequency");
    observable_F_sc_Distribution_wBosonPt_greaterThan10->DrawCopy();

    c5->cd(2);
    observable_F_sc_Distribution_wBosonPt_smallerThan10->SetTitle("Observable F_{sc} distribution for events with W^{+-} boson p_{T} <= 10");
    observable_F_sc_Distribution_wBosonPt_smallerThan10->GetXaxis()->SetTitle("p_{T s}^{jet} / p_{T c}^{jet}");
    observable_F_sc_Distribution_wBosonPt_smallerThan10->GetYaxis()->SetTitle("Frequency");
    observable_F_sc_Distribution_wBosonPt_smallerThan10->DrawCopy();

    // The "modified" plots are those with the updated randomization method in rho and phi for the jet production point in the energy loss region

    TCanvas *c6 = new TCanvas("c6", "Observable F_{sc} distributions different energy losses for null impact parameter (modified)", 2500, 2500);
    c6->Divide(1, 2);

    c6->cd(1);
    observable_F_sc_Distribution_wBosonPt_greaterThan10_diffEnergyLosses->SetTitle("Observable F_{sc} distribution for events with W^{+-} boson p_{T} > 10 for different energy loss functions (modified)");
    observable_F_sc_Distribution_wBosonPt_greaterThan10_diffEnergyLosses->GetXaxis()->SetTitle("p_{T s}^{jet} / p_{T c}^{jet}");
    observable_F_sc_Distribution_wBosonPt_greaterThan10_diffEnergyLosses->GetYaxis()->SetTitle("Frequency");
    observable_F_sc_Distribution_wBosonPt_greaterThan10_diffEnergyLosses->DrawCopy();

    c6->cd(2);
    observable_F_sc_Distribution_wBosonPt_smallerThan10_diffEnergyLosses->SetTitle("Observable F_{sc} distribution for events with W^{+-} boson p_{T} <= 10 for different energy loss functions (modified)");
    observable_F_sc_Distribution_wBosonPt_smallerThan10_diffEnergyLosses->GetXaxis()->SetTitle("p_{T s}^{jet} / p_{T c}^{jet}");
    observable_F_sc_Distribution_wBosonPt_smallerThan10_diffEnergyLosses->GetYaxis()->SetTitle("Frequency");
    observable_F_sc_Distribution_wBosonPt_smallerThan10_diffEnergyLosses->DrawCopy();

    TCanvas *c7 = new TCanvas("c7", "Observable F_{sc} distributions same energy losses for null impact parameter (modified)", 2500, 2500);
    c7->Divide(1, 2);

    c7->cd(1);
    observable_F_sc_Distribution_wBosonPt_greaterThan10_sameEnergyLosses->SetTitle("Observable F_{sc} distribution for events with W^{+-} boson p_{T} > 10 for the same energy loss function (modified)");
    observable_F_sc_Distribution_wBosonPt_greaterThan10_sameEnergyLosses->GetXaxis()->SetTitle("p_{T s}^{jet} / p_{T c}^{jet}");
    observable_F_sc_Distribution_wBosonPt_greaterThan10_sameEnergyLosses->GetYaxis()->SetTitle("Frequency");
    observable_F_sc_Distribution_wBosonPt_greaterThan10_sameEnergyLosses->DrawCopy();

    c7->cd(2);
    observable_F_sc_Distribution_wBosonPt_smallerThan10_sameEnergyLosses->SetTitle("Observable F_{sc} distribution for events with W^{+-} boson p_{T} <= 10 for the same energy loss function (modified)");
    observable_F_sc_Distribution_wBosonPt_smallerThan10_sameEnergyLosses->GetXaxis()->SetTitle("p_{T s}^{jet} / p_{T c}^{jet}");
    observable_F_sc_Distribution_wBosonPt_smallerThan10_sameEnergyLosses->GetYaxis()->SetTitle("Frequency");
    observable_F_sc_Distribution_wBosonPt_smallerThan10_sameEnergyLosses->DrawCopy();

    TCanvas *c8 = new TCanvas("c8", "Observable F_{sc} distributions of normalized ratio for same energy losses and null impact parameter (modified)", 2500, 2500);
    c8->Divide(1, 2);

    c8->cd(1);
    ratio_greater10_same->SetTitle("Normalized ratio of F_{sc} distributions for events with W^{+-} boson p_{T} > 10 for same energy loss functions (modified)");
    ratio_greater10_same->GetXaxis()->SetTitle("p_{T s}^{jet} / p_{T c}^{jet}");
    ratio_greater10_same->GetYaxis()->SetTitle("Frequency");
    ratio_greater10_same->DrawCopy();

    c8->cd(2);
    ratio_smaller10_same->SetTitle("Normalized ratio of F_{sc} distributions for events with W^{+-} boson p_{T} <= 10 for same energy loss functions (modified)");
    ratio_smaller10_same->GetXaxis()->SetTitle("p_{T s}^{jet} / p_{T c}^{jet}");
    ratio_smaller10_same->GetYaxis()->SetTitle("Frequency");
    ratio_smaller10_same->DrawCopy();

    TCanvas *c9 = new TCanvas("c9", "Observable F_{sc} distributions of normalized ratio for distinct energy losses and null impact parameter (modified)", 2500, 2500);
    c9->Divide(1, 2);

    c9->cd(1);
    ratio_greater10_diff->SetTitle("Normalized ratio of F_{sc} distributions for events with W^{+-} boson p_{T} > 10 for different energy loss functions (modified)");
    ratio_greater10_diff->GetXaxis()->SetTitle("p_{T s}^{jet} / p_{T c}^{jet}");
    ratio_greater10_diff->GetYaxis()->SetTitle("Frequency");
    ratio_greater10_diff->DrawCopy();

    c9->cd(2);
    ratio_smaller10_diff->SetTitle("Normalized ratio of F_{sc} distributions for events with W^{+-} boson p_{T} <= 10 for different energy loss functions (modified)");
    ratio_smaller10_diff->GetXaxis()->SetTitle("p_{T s}^{jet} / p_{T c}^{jet}");
    ratio_smaller10_diff->GetYaxis()->SetTitle("Frequency");
    ratio_smaller10_diff->DrawCopy();

    TCanvas *c10 = new TCanvas("c10", "Observable F_{sc} distributions different energy losses for non null impact parameter", 2500, 2500);
    c10->Divide(1, 2);

    c10->cd(1);
    second_observable_F_sc_Distribution_wBosonPt_greaterThan10_diffEnergyLosses->SetTitle("Observable F_{sc} distribution for events with W^{+-} boson p_{T} > 10 for different energy loss functions and non null impact parameter");
    second_observable_F_sc_Distribution_wBosonPt_greaterThan10_diffEnergyLosses->GetXaxis()->SetTitle("p_{T s}^{jet} / p_{T c}^{jet}");
    second_observable_F_sc_Distribution_wBosonPt_greaterThan10_diffEnergyLosses->GetYaxis()->SetTitle("Frequency");
    second_observable_F_sc_Distribution_wBosonPt_greaterThan10_diffEnergyLosses->DrawCopy();

    c10->cd(2);
    second_observable_F_sc_Distribution_wBosonPt_smallerThan10_diffEnergyLosses->SetTitle("Observable F_{sc} distribution for events with W^{+-} boson p_{T} <= 10 for different energy loss functions and non null impact parameter");
    second_observable_F_sc_Distribution_wBosonPt_smallerThan10_diffEnergyLosses->GetXaxis()->SetTitle("p_{T s}^{jet} / p_{T c}^{jet}");
    second_observable_F_sc_Distribution_wBosonPt_smallerThan10_diffEnergyLosses->GetYaxis()->SetTitle("Frequency");
    second_observable_F_sc_Distribution_wBosonPt_smallerThan10_diffEnergyLosses->DrawCopy();

    TCanvas *c11 = new TCanvas("c11", "Observable F_{sc} distributions same energy losses for non null impact parameter", 2500, 2500);
    c11->Divide(1, 2);

    c11->cd(1);
    second_observable_F_sc_Distribution_wBosonPt_greaterThan10_sameEnergyLosses->SetTitle("Observable F_{sc} distribution for events with W^{+-} boson p_{T} > 10 for the same energy loss function and non null impact parameter");
    second_observable_F_sc_Distribution_wBosonPt_greaterThan10_sameEnergyLosses->GetXaxis()->SetTitle("p_{T s}^{jet} / p_{T c}^{jet}");
    second_observable_F_sc_Distribution_wBosonPt_greaterThan10_sameEnergyLosses->GetYaxis()->SetTitle("Frequency");
    second_observable_F_sc_Distribution_wBosonPt_greaterThan10_sameEnergyLosses->DrawCopy();

    c11->cd(2);
    second_observable_F_sc_Distribution_wBosonPt_smallerThan10_sameEnergyLosses->SetTitle("Observable F_{sc} distribution for events with W^{+-} boson p_{T} <= 10 for the same energy loss function and non null impact parameter");
    second_observable_F_sc_Distribution_wBosonPt_smallerThan10_sameEnergyLosses->GetXaxis()->SetTitle("p_{T s}^{jet} / p_{T c}^{jet}");
    second_observable_F_sc_Distribution_wBosonPt_smallerThan10_sameEnergyLosses->GetYaxis()->SetTitle("Frequency");
    second_observable_F_sc_Distribution_wBosonPt_smallerThan10_sameEnergyLosses->DrawCopy();

    TCanvas *c12 = new TCanvas("c12", "Observable F_{sc} distributions of normalized ratio for same energy losses and non null impact parameter (modified)", 2500, 2500);
    c12->Divide(1, 2);

    c12->cd(1);
    second_ratio_greater10_same->SetTitle("Normalized ratio of F_{sc} distributions for events with W^{+-} boson p_{T} > 10 for same energy loss functions (modified)");
    second_ratio_greater10_same->GetXaxis()->SetTitle("p_{T s}^{jet} / p_{T c}^{jet}");
    second_ratio_greater10_same->GetYaxis()->SetTitle("Frequency");
    second_ratio_greater10_same->DrawCopy();

    c12->cd(2);
    second_ratio_smaller10_same->SetTitle("Normalized ratio of F_{sc} distributions for events with W^{+-} boson p_{T} <= 10 for same energy loss functions (modified)");
    second_ratio_smaller10_same->GetXaxis()->SetTitle("p_{T s}^{jet} / p_{T c}^{jet}");
    second_ratio_smaller10_same->GetYaxis()->SetTitle("Frequency");
    second_ratio_smaller10_same->DrawCopy();

    TCanvas *c13 = new TCanvas("c13", "Observable F_{sc} distributions of normalized ratio for distinct energy losses and non null impact parameter (modified)", 2500, 2500);
    c13->Divide(1, 2);

    c13->cd(1);
    second_ratio_greater10_diff->SetTitle("Normalized ratio of F_{sc} distributions for events with W^{+-} boson p_{T} > 10 for different energy loss functions (modified)");
    second_ratio_greater10_diff->GetXaxis()->SetTitle("p_{T s}^{jet} / p_{T c}^{jet}");
    second_ratio_greater10_diff->GetYaxis()->SetTitle("Frequency");
    second_ratio_greater10_diff->DrawCopy();

    c13->cd(2);
    second_ratio_smaller10_diff->SetTitle("Normalized ratio of F_{sc} distributions for events with W^{+-} boson p_{T} <= 10 for different energy loss functions (modified)");
    second_ratio_smaller10_diff->GetXaxis()->SetTitle("p_{T s}^{jet} / p_{T c}^{jet}");
    second_ratio_smaller10_diff->GetYaxis()->SetTitle("Frequency");
    second_ratio_smaller10_diff->DrawCopy();

    TCanvas *c14 = new TCanvas("c14", "Distributions of rho and phi (b=0 and b#neq0)", 2500, 2500);

    c14->Divide(2,2);

    c14->cd(1);
    h_rho_b0->SetTitle("#rho distribution (b = 0)");
    h_rho_b0->GetXaxis()->SetTitle("#rho [fm]");
    h_rho_b0->GetYaxis()->SetTitle("Counts");
    h_rho_b0->DrawCopy();

    c14->cd(2);
    h_phi_b0->SetTitle("#phi distribution (b = 0)");
    h_phi_b0->GetXaxis()->SetTitle("#phi");
    h_phi_b0->GetYaxis()->SetTitle("Counts");
    h_phi_b0->DrawCopy();

    c14->cd(3);
    h_rho_bn0->SetTitle("#rho distribution (b #neq 0)");
    h_rho_bn0->GetXaxis()->SetTitle("#rho [fm]");
    h_rho_bn0->GetYaxis()->SetTitle("Counts");
    h_rho_bn0->DrawCopy();

    c14->cd(4);
    h_phi_bn0->SetTitle("#phi distribution (b #neq 0)");
    h_phi_bn0->GetXaxis()->SetTitle("#phi");
    h_phi_bn0->GetYaxis()->SetTitle("Counts");
    h_phi_bn0->DrawCopy();

    TCanvas *c15 = new TCanvas("c15", "Distributions of (x,y) jet spawning coordinates", 2500, 2500);

    c15->Divide(2,1);

    c15->cd(1);
    xy_b0_jetSpawningCoordinates->SetTitle("(x,y) jet spawning coordinates distribution (b = 0)");
    xy_b0_jetSpawningCoordinates->GetXaxis()->SetTitle("x [fm]");
    xy_b0_jetSpawningCoordinates->GetYaxis()->SetTitle("y [fm]");
    xy_b0_jetSpawningCoordinates->DrawCopy("colz");

    c15->cd(2);
    xy_b_jetSpawningCoordinates->SetTitle("(x,y) jet spawning coordinates distribution (b #neq 0)");
    xy_b_jetSpawningCoordinates->GetXaxis()->SetTitle("x [fm]");
    xy_b_jetSpawningCoordinates->GetYaxis()->SetTitle("y [fm]");
    xy_b_jetSpawningCoordinates->DrawCopy("colz");

    TCanvas *c16 = new TCanvas("c16", "pT distributions before/after quenching", 2500, 2500);

    c16->Divide(2,3);

    c16->cd(1);
    pT_charmDistribution_beforeQuenching->SetTitle("Charm jet p_{T} distribution prior to the quenching");
    pT_charmDistribution_beforeQuenching->DrawCopy();

    c16->cd(2);
    pT_strangeDistribution_beforeQuenching->SetTitle("Strange jet p_{T} distribution prior to the quenching");
    pT_strangeDistribution_beforeQuenching->DrawCopy();

    c16->cd(3);
    null_b_pT_charmDistribution_afterQuenching_diff->SetTitle("Charm jet p_{T} distribution after quenching (b = 0, different energy losses)");
    null_b_pT_charmDistribution_afterQuenching_diff->DrawCopy();

    c16->cd(4);
    null_b_pT_strangeDistribution_afterQuenching_diff->SetTitle("Strange jet p_{T} distribution after quenching (b = 0, different energy losses)");
    null_b_pT_strangeDistribution_afterQuenching_diff->DrawCopy();

    c16->cd(5);
    non_null_b_pT_charmDistribution_afterQuenching_diff->SetTitle("Charm jet p_{T} distribution after quenching (b #neq 0, different energy losses)");
    non_null_b_pT_charmDistribution_afterQuenching_diff->DrawCopy();

    c16->cd(6);
    non_null_b_pT_strangeDistribution_afterQuenching_diff->SetTitle("Strange jet p_{T} distribution after quenching (b #neq 0, different energy losses)");
    non_null_b_pT_strangeDistribution_afterQuenching_diff->DrawCopy();

    TFile *second_outputFile = new TFile("modified2_histogramas_jetR_07_fullList.root", "RECREATE");
    observable_F_sc_Distribution_wBosonPt_greaterThan10->Write();
    observable_F_sc_Distribution_wBosonPt_smallerThan10->Write();
    observable_F_sc_Distribution_wBosonPt_greaterThan10_diffEnergyLosses->Write();
    observable_F_sc_Distribution_wBosonPt_smallerThan10_diffEnergyLosses->Write();
    observable_F_sc_Distribution_wBosonPt_greaterThan10_sameEnergyLosses->Write();
    observable_F_sc_Distribution_wBosonPt_smallerThan10_sameEnergyLosses->Write();
    ratio_greater10_same->Write();
    ratio_smaller10_same->Write();
    ratio_greater10_diff->Write();
    ratio_smaller10_diff->Write();
    second_observable_F_sc_Distribution_wBosonPt_greaterThan10_diffEnergyLosses->Write();
    second_observable_F_sc_Distribution_wBosonPt_smallerThan10_diffEnergyLosses->Write();
    second_observable_F_sc_Distribution_wBosonPt_greaterThan10_sameEnergyLosses->Write();
    second_observable_F_sc_Distribution_wBosonPt_smallerThan10_sameEnergyLosses->Write();
    second_ratio_greater10_same->Write();
    second_ratio_smaller10_same->Write();
    second_ratio_greater10_diff->Write();
    second_ratio_smaller10_diff->Write();
    h_rho_b0->Write();
    h_phi_b0->Write();
    h_rho_bn0->Write();
    h_phi_bn0->Write();
    xy_b0_jetSpawningCoordinates->Write();
    xy_b_jetSpawningCoordinates->Write();
    pT_charmDistribution_beforeQuenching->Write();
    pT_strangeDistribution_beforeQuenching->Write();
    null_b_pT_charmDistribution_afterQuenching_diff->Write();
    null_b_pT_strangeDistribution_afterQuenching_diff->Write();
    non_null_b_pT_charmDistribution_afterQuenching_diff->Write();
    non_null_b_pT_strangeDistribution_afterQuenching_diff->Write();
    second_outputFile->Close();

    file->Close();
}

