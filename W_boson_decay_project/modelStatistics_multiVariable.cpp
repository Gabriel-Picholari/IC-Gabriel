#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMVA/Reader.h>
#include <TGraph.h>
#include <vector>
#include <TLine.h>

void modelStatistics_multiVariable(const char* inputFileName, std::string switch_string, std::string contaminatingGluonMode, const float threshold = 0.60) 
{
    std::string sfx;
    std::string long_sfx;
    std::string short_sfx;
    std::string uppercase_switch_string;
    std::string datasetName;

    if (switch_string == "strange")
    {
        sfx = "_s";
        long_sfx = "_strange";
        short_sfx = "s";
        uppercase_switch_string = "Strange";
    } 
    else if (switch_string == "charm")
    {
        sfx = "_c";
        long_sfx = "_charm";
        short_sfx = "c";
        uppercase_switch_string = "Charm";
    }

    if (switch_string == "charm")  
    {
        if (contaminatingGluonMode == "include") datasetName = "dataset_c_multiVariable_gluonJetsIncluded";
        if (contaminatingGluonMode == "exclude") datasetName = "dataset_c_multiVariable_gluonJetsExcluded";
    }
    if (switch_string == "strange")
    {
        if (contaminatingGluonMode == "include") datasetName = "dataset_s_multiVariable_gluonJetsIncluded";
        if (contaminatingGluonMode == "exclude") datasetName = "dataset_s_multiVariable_gluonJetsExcluded";
    }

    //---------------------------------------------------------------------------------------------------------
    // Criação do objeto Reader para leitura de resultados 
    //---------------------------------------------------------------------------------------------------------

    TMVA::Reader* reader = new TMVA::Reader("!Color:!Silent");

    Float_t pT, nConst, eta, phi, mass, label, eventID, score, nRho;
    Int_t flavor;

    reader->AddVariable(("pT" + sfx).c_str(), &pT);
    reader->AddVariable(("nRho" + sfx).c_str(), &nRho);
    reader->AddVariable(("nConst" + sfx).c_str(), &nConst);

    reader->AddSpectator(("eta" + sfx).c_str(), &eta);
    reader->AddSpectator(("phi" + sfx).c_str(), &phi);
    reader->AddSpectator(("label" + sfx).c_str(), &label);
    reader->AddSpectator(("eventID" + sfx).c_str(), &eventID);
    reader->AddSpectator(("flavor" + sfx).c_str(), &flavor);
    reader->BookMVA("GradBoost", (datasetName + "/weights/TMVAClassification_GradBoost.weights.xml").c_str());

    //---------------------------------------------------------------------------------------------------------
    // Recuperação de TTrees de entrada 
    //---------------------------------------------------------------------------------------------------------

    TFile* inputFile = TFile::Open(inputFileName, "READ");
    TTree* signalTree = (TTree*)inputFile->Get(("SignalTree" + sfx).c_str());
    TTree* backgroundTree = (TTree*)inputFile->Get(("BackgroundTree" + sfx).c_str());

    signalTree->SetBranchAddress(("pT" + sfx).c_str(), &pT);
    signalTree->SetBranchAddress(("eta" + sfx).c_str(), &eta);
    signalTree->SetBranchAddress(("phi" + sfx).c_str(), &phi);
    signalTree->SetBranchAddress(("nConst" + sfx).c_str(), &nConst);
    signalTree->SetBranchAddress(("nRho" + sfx).c_str(), &nRho);
    signalTree->SetBranchAddress(("label" + sfx).c_str(), &label);
    signalTree->SetBranchAddress(("eventID" + sfx).c_str(), &eventID);
    signalTree->SetBranchAddress(("flavor" + sfx).c_str(), &flavor);

    backgroundTree->SetBranchAddress(("pT" + sfx).c_str(), &pT);
    backgroundTree->SetBranchAddress(("eta" + sfx).c_str(), &eta);
    backgroundTree->SetBranchAddress(("phi" + sfx).c_str(), &phi);
    backgroundTree->SetBranchAddress(("nConst" + sfx).c_str(), &nConst);
    backgroundTree->SetBranchAddress(("nRho" + sfx).c_str(), &nRho);
    backgroundTree->SetBranchAddress(("label" + sfx).c_str(), &label);
    backgroundTree->SetBranchAddress(("eventID" + sfx).c_str(), &eventID);
    backgroundTree->SetBranchAddress(("flavor" + sfx).c_str(), &flavor);

    //---------------------------------------------------------------------------------------------------------
    // Analise de desempenho (threshold único)
    //---------------------------------------------------------------------------------------------------------
    Int_t VP = 0, FN = 0, FP = 0, VN = 0;

    TH1F* h_signal = new TH1F("h_signal", "TMVA response for classifier: GradBoost;GradBoost response;Events", 100, -1, 1);
    TH1F* h_background = new TH1F("h_background", "TMVA response for classifier: GradBoost;GradBoost response;Events", 100, -1, 1);

    // Vamos guardar todos os scores e labels
    std::vector<float> allScores;
    std::vector<int>   allLabels;

    for (Long64_t i = 0; i < signalTree->GetEntries(); ++i) 
    {
        signalTree->GetEntry(i);
        score = reader->EvaluateMVA("GradBoost");
        h_signal->Fill(score);
        allScores.push_back(score);
        allLabels.push_back((int)label);

        if (label == 1) 
        {
            if (score >= threshold) VP++; else FN++;
        } 
        else 
        {
            if (score >= threshold) FP++; else VN++;
        }
    }

    for (Long64_t i = 0; i < backgroundTree->GetEntries(); ++i) 
    {
        backgroundTree->GetEntry(i);
        score = reader->EvaluateMVA("GradBoost");
        h_background->Fill(score);
        allScores.push_back(score);
        allLabels.push_back((int)label);

        if (label == 1) 
        {
            if (score >= threshold) VP++; else FN++;
        } 
        else 
        {
            if (score >= threshold) FP++; else VN++;
        }
    }

    Float_t eficiencia;
    if (VP + FN > 0) 
    {
        eficiencia = (float)VP / (VP + FN);
    } 
    else 
    {
        eficiencia = 0;
    }

    Float_t pureza;
    if (VP + FP > 0) 
    {
        pureza = (float)VP / (VP + FP);
    } 
    else 
    {
        pureza = 0;
    }

    std::cout << "--------------------------------------" << std::endl;
    std::cout << "\nMatriz de Confusão (threshold = " << threshold << "):" << std::endl;
    std::cout << "--------------------------------------" << std::endl;
    std::cout << "             | Pred: S | Pred: B     " << std::endl;
    std::cout << "--------------------------------------" << std::endl;
    std::cout << "Real: S (1)  |   " << VP << "     |   " << FN << std::endl;
    std::cout << "Real: B (0)  |   " << FP << "     |   " << VN << std::endl;
    std::cout << "--------------------------------------" << std::endl;
    std::cout << "Eficiência (Recall) = " << eficiencia << std::endl;
    std::cout << "Pureza (Precision)  = " << pureza << std::endl;

//---------------------------------------------------------------------------------------------------------
// Plot original histograms
//---------------------------------------------------------------------------------------------------------

    h_background->Scale(h_signal->Integral() / h_background->Integral());

    TCanvas* c1 = new TCanvas("c1", ("GradBoost Score Distribution (" + uppercase_switch_string + ")").c_str(), 900, 700);
    c1->SetGrid();

    h_background->SetLineColor(kRed);
    h_background->SetLineWidth(3);
    h_background->SetMarkerColor(kRed);
    h_background->SetMarkerStyle(4);   
    h_background->SetMarkerSize(1);



    h_signal->SetLineColor(kBlue);
    h_signal->SetLineWidth(2);

    h_signal->SetTitle(("GradBoost Score Distribution for " + uppercase_switch_string + " Jets;Score;Entries").c_str());
    h_background->GetXaxis()->SetTitleSize(0.05);
    h_background->GetYaxis()->SetTitleSize(0.05);
    h_background->SetTitle(("GradBoost Score Distribution for " + uppercase_switch_string + " Jets;Score;Entries").c_str());

    h_background->DrawCopy("E1P");
    h_signal->DrawCopy("same");

    
    TLegend* leg1 = new TLegend(0.15, 0.75, 0.38, 0.88);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->SetTextSize(0.06);
    leg1->AddEntry(h_signal, "Signal", "l");
    leg1->AddEntry(h_background, "Background", "l");
    leg1->Draw();
    

    //c1->SaveAs(("GradBoost_Score_Distribution" + long_sfx + ".png").c_str());

    //---------------------------------------------------------------------------------------------------------
    // Scan em thresholds
    //---------------------------------------------------------------------------------------------------------

    int nSteps = 1000;
    float tmin = -1.0, tmax = 1.0;

    std::vector<double> vx, vEff, vPur;
    vx.reserve(nSteps); vEff.reserve(nSteps); vPur.reserve(nSteps);

    for (int k=0; k<nSteps; ++k) 
    {
        Float_t thr = tmin + (tmax-tmin)*k/(nSteps-1);
        int vp=0, fn=0, fp=0, vn=0;

        for (size_t j=0;j<allScores.size();++j) 
        {
            bool predS = (allScores[j] >= thr);
            if (allLabels[j]==1) 
            { 
                if(predS) vp++; 
                else fn++; 
            }
            else 
            { 
                if(predS) fp++; 
                else vn++; 
            }
        }

        Float_t eff, pur;
        if (vp + fn > 0) 
        {
            eff = (float)vp / (vp + fn);
        } 
        else 
        {
            eff = 0.0;
        }

        if (vp + fp > 0) 
        {
            pur = (float)vp / (vp + fp);
        } 
        else 
        {
            pur = 0.0;
        }

        vx.push_back(thr);
        vEff.push_back(eff);
        vPur.push_back(pur);
    }

    gStyle->SetOptStat(0);

    TCanvas* c2 = new TCanvas("c2", ("Efficiency and Purity vs Threshold (" + uppercase_switch_string + ")").c_str(), 1800, 600);

    c2->SetGrid();
    c2->SetLeftMargin(0.12);
    c2->SetBottomMargin(0.12);
    c2->SetRightMargin(0.05);
    c2->SetTopMargin(0.08);

    TGraph* gEff = new TGraph(nSteps, vx.data(), vEff.data()); gEff->SetTitle(("Efficiency and Purity vs Threshold for " + uppercase_switch_string + " Jets").c_str());

    gEff->SetLineColor(kBlue+1);
    gEff->SetMarkerColor(kBlue+1);
    gEff->SetLineWidth(3);

    gEff->GetXaxis()->SetTitle("Threshold");
    gEff->GetYaxis()->SetTitle("Efficiency / Purity");
    gEff->GetYaxis()->SetRangeUser(0.0,1.05);
    gEff->GetXaxis()->SetRangeUser(-1,1);


    gEff->GetXaxis()->SetTitleSize(0.06);
    gEff->GetYaxis()->SetTitleSize(0.06);
    gEff->GetXaxis()->SetLabelSize(0.04);
    gEff->GetYaxis()->SetLabelSize(0.04);

    gEff->Draw("ALP");

    TGraph* gPur = new TGraph(nSteps, vx.data(), vPur.data());
    gPur->SetLineColor(kMagenta+2);
    gPur->SetLineWidth(3);
    gPur->Draw("LP SAME");

    // Threshold
    TLine* l1 = new TLine(threshold,0,threshold,1.05);
    l1->SetLineStyle(2);
    l1->SetLineWidth(3);
    l1->SetLineColor(kGray+2);
    l1->Draw();

    // Legenda
    TLegend* leg = new TLegend(0.15,0.72,0.38,0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(gEff,"Efficiency","lp");
    leg->AddEntry(gPur,"Purity","lp");
    leg->AddEntry(l1,Form("Threshold = %.3f",threshold),"l");
    leg->SetTextSize(0.06);
    leg->Draw();

    c2->Update();

    //c2->SaveAs(("Efficiency_Purity_vs_Threshold_" + long_sfx + ".png").c_str());

    inputFile->Close();
    delete reader;
}
