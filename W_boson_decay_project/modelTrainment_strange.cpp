//---------------------------------------------------------------------------------------------------------//
//                                                                                                         //
// Macro that utilizes the training and testing trees from both signal and background data to finally      //
// train, test and evaluate the logistic regression method. Note that the library used is TMVA. The macro  //
// also reads the results and print the confusion matrix.                                                  //
//                                                                                                         //
//---------------------------------------------------------------------------------------------------------//

#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TMVA/Reader.h>
#include <TMVA/Factory.h>
#include <TMVA/DataLoader.h>

void modelTrainment_strange(const char *fileName) {
    
    TFile* outputFile = TFile::Open("TMVAOutput_2var_strange.root", "RECREATE");

    TFile* inputFile = TFile::Open(fileName, "READ");

    TTree* signalTree_s = dynamic_cast<TTree *>(inputFile->Get("SignalTree_s"));
    TTree* backgroundTree_s = dynamic_cast<TTree *>(inputFile->Get("BackgroundTree_s"));

    TMVA::Factory factory("TMVAClassification", outputFile, "AnalysisType=Classification");
    TMVA::DataLoader loader("dataset_s");

    loader.AddVariable("pT_s", 'F');
    loader.AddVariable("nConst_s", 'F');
    
    
    loader.AddSpectator("eta_s", "F");
    loader.AddSpectator("phi_s", "F");
    loader.AddSpectator("mass_s", "F");
    loader.AddSpectator("label_s", "F");
    loader.AddSpectator("eventID_s", "F");

    loader.AddSignalTree(signalTree_s, 1.0);
    loader.AddBackgroundTree(backgroundTree_s, 1.0);

    //---------------------------------------------------------------------------------------------------------
    // Divisao das Trees de sinal e fundo em treinamento e teste
    //---------------------------------------------------------------------------------------------------------

    TCut mycut = "";
    Int_t nTest_Signal = signalTree_s->GetEntries() * 0.2;
    Int_t nTrain_Signal = signalTree_s->GetEntries() * 0.8;

    Int_t nTest_Background = backgroundTree_s->GetEntries() * 0.2;
    Int_t nTrain_Background = backgroundTree_s->GetEntries() * 0.8;
    
    TString options = TString::Format("SplitMode=Random:SplitSeed=0:NormMode=None:nTrain_Signal=%d:nTrain_Background=%d:nTest_Signal=%d:nTest_Background=%d:!V", nTrain_Signal, nTrain_Background, nTest_Signal, nTest_Background);
    loader.PrepareTrainingAndTestTree(mycut, options.Data());

    //---------------------------------------------------------------------------------------------------------
    // Escolha do metodo de treinamento e concretizacao de teste, treino e evaluacao do modelo
    //---------------------------------------------------------------------------------------------------------

    factory.BookMethod(&loader, TMVA::Types::kLikelihood, "Likelihood", "H:!V");

    factory.TrainAllMethods();
    factory.TestAllMethods();
    factory.EvaluateAllMethods();

    inputFile->Close();
    outputFile->Close();    
    }
