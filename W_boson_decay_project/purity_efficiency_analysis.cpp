//===========================================================================================================================================================================
// Cuts in jets distributions analysis:
// This macro has the purpose of analysing the distributions of the discriminatory variables chosen to train the classifiers. The cuts determined here through analysis of
// the signal and background distributions will be applied in a classical analysis of the W boson decay. Thus, we will be able to compare the classical and the machine
// learning approaches.
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

void cuts_analysis(std::string switch_string, const char* strange_file)
{
    gSystem->Load("libEG");
    gSystem->Load("libEGPythia8");

    //---------------------------------------------------------------------------------------------------------
    // Initialization of the .root file and TTrees
    //---------------------------------------------------------------------------------------------------------

    Float_t pT_s, nConst_s, nRho_s, eventID_s;
    Float_t pT_c, nConst_c, nRho_c, eventID_c;

    TFile *s_file = TFile::Open(strange_file, "READ");

    TTree *signalTree_s = dynamic_cast<TTree *>(s_file->Get("SignalTree_s"));
    signalTree_s->SetBranchAddress("pT_s", &pT_s);
    signalTree_s->SetBranchAddress("nConst_s", &nConst_s);
    signalTree_s->SetBranchAddress("nRho_s", &nRho_s);
    signalTree_s->SetBranchAddress("eventID_s", &eventID_s);

    TTree *backgroundTree_s = dynamic_cast<TTree *>(s_file->Get("BackgroundTree_s")); 
    backgroundTree_s->SetBranchAddress("pT_s", &pT_s);
    backgroundTree_s->SetBranchAddress("nConst_s", &nConst_s);
    backgroundTree_s->SetBranchAddress("nRho_s", &nRho_s);
    backgroundTree_s->SetBranchAddress("eventID_s", &eventID_s);

    Long64_t ne_signal_s = signalTree_s->GetEntries();
    Long64_t ne_background_s = backgroundTree_s->GetEntries();

    //---------------------------------------------------------------------------------------------------------

    TFile *c_file = TFile::Open(charm_file, "READ");

    TTree *signalTree_c = dynamic_cast<TTree *>(c_file->Get("SignalTree_c"));
    signalTree_c->SetBranchAddress("pT_c", &pT_c);
    signalTree_c->SetBranchAddress("nConst_c", &nConst_c);
    signalTree_c->SetBranchAddress("nRho_c", &nRho_c);
    signalTree_c->SetBranchAddress("eventID_c", &eventID_c);

    TTree *backgroundTree_c = dynamic_cast<TTree *>(c_file->Get("BackgroundTree_c"));
    backgroundTree_c->SetBranchAddress("pT_c", &pT_c);
    backgroundTree_c->SetBranchAddress("nConst_c", &nConst_c);
    backgroundTree_c->SetBranchAddress("nRho_c", &nRho_c);
    backgroundTree_c->SetBranchAddress("eventID_c", &eventID_c);

    Long64_t ne_signal_c = signalTree_c->GetEntries();
    Long64_t ne_background_c = backgroundTree_c->GetEntries();

    //---------------------------------------------------------------------------------------------------------
    // Analysis
    //---------------------------------------------------------------------------------------------------------

    for(Float_t threshold = 0; threshold < 100; threshold += 0.01)
    {

        for (Long64_t i = 0; i < ne_signal_s; i++)
        {
            signalTree_s->GetEntry(i);
            
        }

        for (Long64_t i = 0; i < ne_background_s; i++)
        {
            backgroundTree_s->GetEntry(i);
            
        }

        //---------------------------------------------------------------------------------------------------------

        for (Long64_t i = 0; i < ne_signal_c; i++)
        {
            signalTree_c->GetEntry(i);
            
        }

        for (Long64_t i = 0; i < ne_background_c; i++)
        {
            backgroundTree_c->GetEntry(i);

        }
    }

}
