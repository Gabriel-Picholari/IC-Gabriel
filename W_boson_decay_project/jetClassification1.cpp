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
#include "TRandom3.h"
#include "TParticle.h"
#include <unordered_set>
#include "TClonesArray.h"
#include <TLorentzVector.h>
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"


fastjet::PseudoJet quenchedJet(const fastjet::PseudoJet &jet, double deltaE)
{
    Double_t E_old = jet.E();
    Double_t E_new = E_old - deltaE;

    Double_t m  = jet.m();
    Double_t p2 = jet.px()*jet.px() + jet.py()*jet.py() + jet.pz()*jet.pz();
    Double_t p  = std::sqrt(p2);

    Double_t p_new2 = E_new*E_new - m*m;

    // We suppose that the energy loss does not change the direction of the new 3-momentum vector that shows compatibility (within the Minkowski metric) with the new energy.
    // Thus, we calculate the factor that keeps the proportionality between the old and new momentum vectors and then apply it to each component of the old momentum vector.
    Double_t p_new = std::sqrt(p_new2);
    Double_t scale = p_new / p;

    Double_t px_new = scale * jet.px();
    Double_t py_new = scale * jet.py();
    Double_t pz_new = scale * jet.pz();

    fastjet::PseudoJet q(px_new, py_new, pz_new, E_new); // This composes the new 4-momentum of the quenched jet with redefined energy and momentum
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

    TH1F *observable_F_sc_Distribution_wBosonPt_greaterThan10_diffEnergyLosses = new TH1F("observable_F_sc_Distribution_wBosonPt_greaterThan10_diffEnergyLosses", "Strange jet by charm jet p_{T} - F_{sc} distribution for events with W^{+-} boson p_{T} > 10 for different energy loss funtions", 100, 0, 10);    
    TH1F *observable_F_sc_Distribution_wBosonPt_smallerThan10_diffEnergyLosses = new TH1F("observable_F_sc_Distribution_wBosonPt_smallerThan10_diffEnergyLosses", "Strange jet by charm jet p_{T} - F_{sc} distribution for events with W^{+-} boson p_{T} <= 10 for different energy loss funtions", 100, 0, 10);

    TH1F *observable_F_sc_Distribution_wBosonPt_greaterThan10_sameEnergyLosses = new TH1F("observable_F_sc_Distribution_wBosonPt_greaterThan10_sameEnergyLosses", "Strange jet by charm jet p_{T} - F_{sc} distribution for events with W^{+-} boson p_{T} > 10 for the same energy loss function", 100, 0, 10);    
    TH1F *observable_F_sc_Distribution_wBosonPt_smallerThan10_sameEnergyLosses = new TH1F("observable_F_sc_Distribution_wBosonPt_smallerThan10_sameEnergyLosses", "Strange jet by charm jet p_{T} - F_{sc} distribution for events with W^{+-} boson p_{T} <= 10 for the same energy loss function", 100, 0, 10);

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

            if (hasCharmConstituent && charmRatio > 0.6) // Please, refer to 5th september research log for the reasoning behind the 0.6 lower limit
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

            if (hasStrangeConstituent && strangeRatio > 0.6) // Same as before
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
                // Thus, I guess we can just reopen the jet into it's constituents and get them from directly from there

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
        Double_t max_pt_c = -1.0; // We use a negative value so the fist jet in the vector will always be the macimum pT jet
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

        if (deltaPhi >= backToback_lowerLimit || deltaPhi < backToback_upperLimit) // Ensuring the back-to-back condition is fulfilled for all jets from now on
        { 

            //---------------------------------------------------------------------------------------------------------
            // No energy loss block
            //---------------------------------------------------------------------------------------------------------
            
            Float_t F_sc = event_strange_jet.pt() / event_charmed_jet.pt(); // New observable - for more information, check November 6th, 2025 research log

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
            // Energy loss block (b=0)
            //---------------------------------------------------------------------------------------------------------

            Float_t randomRho = generator.Uniform(0, PbNucleusRadius);
            Float_t randomPhi = generator.Uniform(0, 2*TMath::Pi());

            // We take the jets we own and place them into a new polar coordinate system. The jets are then placed in (randomRho, randomPhi).
            // Since by hypotesis we'll take them as back to back, this point will be the starting point for both of them in this new system.
            // Our intentions are to calculate the path lenght each jet would have within a Pb nucleus if they were created at this random point and traveled in opposite directions.
            // We can account for two situations: either the impact parameter is zero (central collision) or it has a finite value (peripheral collision). For the time being, we'll take b = 0.

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
        }

        particles_fastjet.clear();
        jets.clear();
        jets_array->Clear();
        quarks->Clear();
    
    } // End of event loop equivalent
    
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
    primary_CharmRatioHist->Draw();

    primary_CharmRatioHist->SetLineColor(kRed);
    secondary_CharmRatioHist->Draw("same");

    TCanvas *c3 = new TCanvas("c3", "Strange", 2500, 2500);
    c3->Divide(1, 1);

    c3->cd(1);
    //secondary_StrangeRatioHist->SetTitle("Distribution of strange jet-to-quark p_{T} ratio with no p_{T} cut applied"); // (Title is correct when no cut is wanted!)

    primary_StrangeRatioHist->SetTitle("Distribution of strange jet-to-quark p_{T} ratio with a 5 GeV/c jet p_{T} cut"); // (Title is correct when a cut is wanted!)
    
    primary_StrangeRatioHist->GetXaxis()->SetTitle("Ratio");
    primary_StrangeRatioHist->GetYaxis()->SetTitle("Frequency");
    primary_StrangeRatioHist->SetLineColor(kGreen);
    primary_StrangeRatioHist->Draw();
    
    secondary_StrangeRatioHist->SetLineColor(kRed);
    secondary_StrangeRatioHist->Draw("same");

    TCanvas *c4 = new TCanvas("c4", "Missing particles PDG", 2500, 2500);
    c4->Divide(1, 2);

    c4->cd(1);
    missingStrangeConstituentsPdgMap->SetTitle("Potentially missing particles PDGs for strange jets list");
    missingStrangeConstituentsPdgMap->GetXaxis()->SetTitle("PDG");
    missingStrangeConstituentsPdgMap->GetYaxis()->SetTitle("Frequency");
    missingStrangeConstituentsPdgMap->Draw();

    c4->cd(2);
    missingCharmConstituentsPdgMap->SetTitle("Potentially missing particles PDGs for charm jets list");
    missingCharmConstituentsPdgMap->GetXaxis()->SetTitle("PDG");
    missingCharmConstituentsPdgMap->GetYaxis()->SetTitle("Frequency");
    missingCharmConstituentsPdgMap->Draw();

    TFile *outputFile = new TFile("histogramas_jetR_07_fullList.root", "RECREATE");
    invariantMass->Write();
    primary_CharmRatioHist->Write();
    secondary_CharmRatioHist->Write();
    primary_StrangeRatioHist->Write();
    secondary_StrangeRatioHist->Write();
    outputFile->Close();

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

    TCanvas *c6 = new TCanvas("c6", "Observable F_{sc} distributions different energy losses", 2500, 2500);
    c6->Divide(1, 2);

    c6->cd(1);
    observable_F_sc_Distribution_wBosonPt_greaterThan10_diffEnergyLosses->SetTitle("Observable F_{sc} distribution for events with W^{+-} boson p_{T} > 10 for different energy loss functions");
    observable_F_sc_Distribution_wBosonPt_greaterThan10_diffEnergyLosses->GetXaxis()->SetTitle("p_{T s}^{jet} / p_{T c}^{jet}");
    observable_F_sc_Distribution_wBosonPt_greaterThan10_diffEnergyLosses->GetYaxis()->SetTitle("Frequency");
    observable_F_sc_Distribution_wBosonPt_greaterThan10_diffEnergyLosses->DrawCopy();

    c6->cd(2);
    observable_F_sc_Distribution_wBosonPt_smallerThan10_diffEnergyLosses->SetTitle("Observable F_{sc} distribution for events with W^{+-} boson p_{T} <= 10 for different energy loss functions");
    observable_F_sc_Distribution_wBosonPt_smallerThan10_diffEnergyLosses->GetXaxis()->SetTitle("p_{T s}^{jet} / p_{T c}^{jet}");
    observable_F_sc_Distribution_wBosonPt_smallerThan10_diffEnergyLosses->GetYaxis()->SetTitle("Frequency");
    observable_F_sc_Distribution_wBosonPt_smallerThan10_diffEnergyLosses->DrawCopy();

    TCanvas *c7 = new TCanvas("c7", "Observable F_{sc} distributions same energy losses", 2500, 2500);
    c7->Divide(1, 2);

    c7->cd(1);
    observable_F_sc_Distribution_wBosonPt_greaterThan10_sameEnergyLosses->SetTitle("Observable F_{sc} distribution for events with W^{+-} boson p_{T} > 10 for the same energy loss function");
    observable_F_sc_Distribution_wBosonPt_greaterThan10_sameEnergyLosses->GetXaxis()->SetTitle("p_{T s}^{jet} / p_{T c}^{jet}");
    observable_F_sc_Distribution_wBosonPt_greaterThan10_sameEnergyLosses->GetYaxis()->SetTitle("Frequency");
    observable_F_sc_Distribution_wBosonPt_greaterThan10_sameEnergyLosses->DrawCopy();

    c7->cd(2);
    observable_F_sc_Distribution_wBosonPt_smallerThan10_sameEnergyLosses->SetTitle("Observable F_{sc} distribution for events with W^{+-} boson p_{T} <= 10 for the same energy loss function");
    observable_F_sc_Distribution_wBosonPt_smallerThan10_sameEnergyLosses->GetXaxis()->SetTitle("p_{T s}^{jet} / p_{T c}^{jet}");
    observable_F_sc_Distribution_wBosonPt_smallerThan10_sameEnergyLosses->GetYaxis()->SetTitle("Frequency");
    observable_F_sc_Distribution_wBosonPt_smallerThan10_sameEnergyLosses->DrawCopy();

    file->Close();
}

