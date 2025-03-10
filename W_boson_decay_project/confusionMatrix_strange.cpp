//---------------------------------------------------------------------------------------------------------//
//                                                                                                         //
// This code retrieves the trees of testing and verifies their compatibility with the training trees       //
// in order to determine true and false positives and negatives. Then, a confusion matrix is built.        //
// The code utilizes a loop counter strategy. Important to mention that the file read by this macro is     //
// the one that creates both the signal and bgr Trees and NOT the modelTraining.cpp output.                //
//                                                                                                         //
//---------------------------------------------------------------------------------------------------------//

#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <TROOT.h>
#include <TMVA/Reader.h>
#include <TMVA/Results.h>
#include <TMVA/DataLoader.h>

void confusionMatrix_strange(const char *fileName)
{
    //---------------------------------------------------------------------------------------------------------
    // Recuperacao da Tree de teste
    //---------------------------------------------------------------------------------------------------------

    TFile* inputFile = TFile::Open(fileName, "READ");
    TTree* testTree = dynamic_cast<TTree*>(inputFile->Get("dataset_s/TestTree"));

    //---------------------------------------------------------------------------------------------------------
    // Inicializacao do objeto reader recuperando o arquivo xml gerado pela macro modelTrainment.cpp
    //---------------------------------------------------------------------------------------------------------

    TMVA::Reader reader;

    Float_t pT_s, label_s, nConst_s;

    reader.AddVariable("pT_s", &pT_s);
    reader.AddVariable("nConst_s", &nConst_s);
    reader.AddSpectator("label_s", &label_s);

    reader.BookMVA("LogisticRegression", "dataset_s/weights/TMVARegression_LogisticRegression.weights.xml");

    //---------------------------------------------------------------------------------------------------------
    // Loop que contabiliza os verdadeiros positivos (TP), verdadeiros negativos (TN), falsos positivos (TP) e
    // falsos negativos (FN) usando o objeto reader
    //---------------------------------------------------------------------------------------------------------

    int TP = 0, TN = 0, FP = 0, FN = 0;

    testTree->SetBranchAddress("label_s", &label_s);
    testTree->SetBranchAddress("pT_s", &pT_s);
    testTree->SetBranchAddress("nConst_s", &nConst_s);


    Long64_t nEntries = testTree->GetEntries();

    for (Long64_t i = 0; i < nEntries; ++i) 
    {
        testTree->GetEntry(i);
        double response = reader.EvaluateMVA("LogisticRegression");

        int predictedLabel;
        if ( response >= 0.5 )
        {
            predictedLabel = 1;
        }
        else
        {
            predictedLabel = 0;
        }

        if (predictedLabel == 1 && label_s == 1) 
        {
            TP++;
        } 
        else if (predictedLabel == 0 && label_s == 0) 
        {
            TN++;
        } 
        else if (predictedLabel == 1 && label_s == 0) 
        {
            FP++;
        } 
        else if (predictedLabel == 0 && label_s == 1) 
        {
            FN++;
        }
        
    }

    //---------------------------------------------------------------------------------------------------------
    // Matriz de confusao
    //---------------------------------------------------------------------------------------------------------
    
    std::cout << "Confusion Matrix:" << std::endl;
    std::cout << "TP: " << TP << "  FN: " << FN << std::endl;
    std::cout << "FP: " << FP << "  TN: " << TN << std::endl;

    inputFile->Close();
}