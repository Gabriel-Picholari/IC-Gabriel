// Corrigido

#include "TSystem.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TMath.h"

void decaimento_JPsi(){

  // Randomizacao do angulo e do momento
  TRandom3 *rphi = new TRandom3(0);
  TRandom3 *rpJPsi = new TRandom3(0); // No laboratorio

Int_t nev = 100000;

// Histogramas para o decaimento via randomizacao por dois metodos: uniforme e exponencial E graficos de dispersao
TH1F *hist_ele = new TH1F("hist_ele", "Momentos dos Eletrons (Uni)", 100, 0, 0);
TH1F *hist_pos = new TH1F("hist_pos","Momentos dos Positrons (Uni)", 100, 0, 0);
TH1F *hist_JPsi = new TH1F("hist_JPsi", "Momentos dos J/Psi (Uni)", 100, 0, 0);

TH1F *hist_ele2 = new TH1F("hist_ele2", "Momentos dos Eletrons (Expo)", 100, 0, 0);
TH1F *hist_pos2 = new TH1F("hist_pos2","Momentos dos Positrons (Expo)", 100, 0, 0);
TH1F *hist_JPsi2 = new TH1F("hist_JPsi2", "Momentos dos J/Psi (Expo)", 100, 0, 0);

TH2F *dis_ele_pos_uni = new TH2F("dis1","Distribuicao (Metodo Uniforme), Momento Eletron (GeV/c), Momento Positron (GeV/c)", 100, 0, 10, 100, 0, 10);
TH2F *dis_ele_pos_exp = new TH2F("dis2","Distribuicao (Metodo Exponencial), Momento Eletron (GeV/c), Momento Positron (GeV/c)", 100, 0, 10, 100, 0, 10);

  for (Int_t iev = 0; iev < nev; iev++ ){

  // Dados (ambos os metodos)
  Double_t ephi = rphi->Uniform(0, TMath::TwoPi());
  Double_t epJPsi = rpJPsi->Uniform(0, 5);
  Double_t epJPsi2 = rpJPsi->Exp(1);

  Double_t m_JPsi = 3.0969; // GeV
  Double_t m_ele = 0.00051099895; // GeV c

  // Valores necessarios para a transformacao de Lorentz (ambos os metodos)
  Double_t beta_lor = epJPsi / TMath::Sqrt(TMath::Power(epJPsi, 2) + TMath::Power(m_JPsi, 2));
  Double_t gamma_lor = 1 / TMath::Sqrt(1 - TMath::Power(beta_lor, 2));

  Double_t beta_lor2 = epJPsi2 / TMath::Sqrt(TMath::Power(epJPsi2, 2) + TMath::Power(m_JPsi, 2));
  Double_t gamma_lor2 = 1 / TMath::Sqrt(1 - TMath::Power(beta_lor2, 2));

  //========================================================================================

  // Momentos do positron e eletron no referencial do centro de massa para o metodo uniforme

  // Momento total do eletron e do positron
  Double_t ptot_ele = TMath::Sqrt( (TMath::Power(m_JPsi, 2)/4) - TMath::Power(m_ele, 2) );
  Double_t ptot_pos = ptot_ele;

  // Componentes dos momentos do eletron
  Double_t px_ele = ptot_ele * cos(ephi);
  Double_t py_ele = ptot_ele * sin(ephi);
  Double_t pz_ele = 0;

  // Componentes dos momentos do positron
  Double_t px_pos = ptot_ele * cos(ephi + TMath::Pi());
  Double_t py_pos = ptot_ele * sin(ephi + TMath::Pi());
  Double_t pz_pos = 0;

  //---------------------------------------------------------------------------------------

  // Momentos do positron e eletron no referencial do laboratorio

  // Energia
  Double_t energ_lab_ele = (m_JPsi/2)*(gamma_lor - beta_lor * gamma_lor);
  Double_t energ_lab_pos = energ_lab_ele;

  // Componentes do momentos eletron no lab e momento total
  Double_t px_lab_ele = gamma_lor * (px_ele - beta_lor * m_JPsi);
  Double_t py_lab_ele = py_ele;
  Double_t pz_lab_ele = pz_ele;
  Double_t ptot_lab_ele = TMath::Sqrt(px_lab_ele * px_lab_ele + py_lab_ele * py_lab_ele);

  // Componentes do momento positron no laboratorio e momento total
  Double_t px_lab_pos = gamma_lor * (px_pos - beta_lor * m_JPsi);
  Double_t py_lab_pos = py_pos;
  Double_t pz_lab_pos = pz_pos;
  Double_t ptot_lab_pos = TMath::Sqrt(px_lab_pos * px_lab_pos + py_lab_pos * py_lab_pos);

  //=======================================================================================

  // Momentos do positron e eletron no referencial do centro de massa para o metodo exponencial

  // Momento total do eletron e do positron
  Double_t ptot_ele2 = TMath::Sqrt( (TMath::Power(m_JPsi, 2)/4) - TMath::Power(m_ele, 2) );
  Double_t ptot_pos2 = ptot_ele2;

  // Componentes dos momentos do eletron
  Double_t px_ele2 = ptot_ele2 * cos(ephi);
  Double_t py_ele2 = ptot_ele2 * sin(ephi);
  Double_t pz_ele2 = 0;

  // Componentes dos momentos do positron
  Double_t px_pos2 = ptot_ele2 * cos(ephi + TMath::Pi());
  Double_t py_pos2 = ptot_ele2 * sin(ephi + TMath::Pi());
  Double_t pz_pos2 = 0;

  //---------------------------------------------------------------------------------------

  // Momentos do positron e eletron no referencial do laboratorio

  // Energia
  Double_t energ_lab_ele2 = (m_JPsi/2)*(gamma_lor2 - beta_lor2 * gamma_lor2);
  Double_t energ_lab_pos2 = energ_lab_ele2;

  // Componentes do momentos eletron no lab e momento total
  Double_t px_lab_ele2 = gamma_lor2 * (px_ele2 - beta_lor2 * m_JPsi);
  Double_t py_lab_ele2 = py_ele2;
  Double_t pz_lab_ele2 = pz_ele2;
  Double_t ptot_lab_ele2 = TMath::Sqrt(px_lab_ele2 * px_lab_ele2 + py_lab_ele2 * py_lab_ele2);



  // Componentes do momento positron no laboratorio e momento total
  Double_t px_lab_pos2 = gamma_lor2 * (px_pos2 - beta_lor2 * m_JPsi);
  Double_t py_lab_pos2 = py_pos2;
  Double_t pz_lab_pos2 = pz_pos2;
  Double_t ptot_lab_pos2 = TMath::Sqrt(px_lab_pos2 * px_lab_pos2 + py_lab_pos2 * py_lab_pos2);

//=======================================================================================

  // Preenchendo os histogramas
  hist_ele->Fill(ptot_lab_ele);
  hist_pos->Fill(ptot_lab_pos);
  hist_JPsi->Fill(epJPsi);

  hist_ele2->Fill(ptot_lab_ele2);
  hist_pos2->Fill(ptot_lab_pos2);
  hist_JPsi2->Fill(epJPsi2);

  dis_ele_pos_uni->Fill(ptot_lab_ele, ptot_lab_pos);
  dis_ele_pos_exp->Fill(ptot_lab_ele2, ptot_lab_pos2);
}

  TCanvas* canvas = new TCanvas("canvas", "Histogramas", 2500, 2500);

  canvas->Divide(2, 4);

  canvas->cd(1);
  hist_JPsi->SetTitle("Momentos do J/Psi (Metodo Uniforme)");
  hist_JPsi->GetXaxis()->SetTitle("Momento (GeV/c)");
  hist_JPsi->GetYaxis()->SetTitle("Frequencia");
  hist_JPsi->Draw();

  canvas->cd(2);
  hist_JPsi2->SetTitle("Momentos do J/Psi (Metodo Exponencial)");
  hist_JPsi2->GetXaxis()->SetTitle("Momento (GeV/c)");
  hist_JPsi2->GetYaxis()->SetTitle("Frequencia");
  hist_JPsi2->Draw();

  canvas->cd(3);
  hist_ele->SetTitle("Momentos do eletron (Metodo Uniforme)");
  hist_ele->GetXaxis()->SetTitle("Momento (GeV/c)");
  hist_ele->GetYaxis()->SetTitle("Frequencia");
  hist_ele->Draw();

  canvas->cd(4);
  hist_ele2->SetTitle("Momentos do eletron (Metodo Exponencial)");
  hist_ele2->GetXaxis()->SetTitle("Momento (GeV/c)");
  hist_ele2->GetYaxis()->SetTitle("Frequencia");
  hist_ele2->Draw();

  canvas->cd(5);
  hist_pos->SetTitle("Momentos do positron (Metodo Uniforme)");
  hist_pos->GetXaxis()->SetTitle("Momento (GeV/c)");
  hist_pos->GetYaxis()->SetTitle("Frequencia");
  hist_pos->Draw();

  canvas->cd(6);
  hist_pos2->SetTitle("Momentos do positron (Metodo Exponencial)");
  hist_pos2->GetXaxis()->SetTitle("Momento (GeV/c)");
  hist_pos2->GetYaxis()->SetTitle("Frequencia");
  hist_pos2->Draw();

  canvas->cd(7);
  dis_ele_pos_uni->SetTitle("Dispersao do decaimento no lab (Metodo Uniforme)");
  dis_ele_pos_uni->GetXaxis()->SetTitle("Momento Eletron (GeV/c)");
  dis_ele_pos_uni->GetYaxis()->SetTitle("Momento Positron (GeV/c)");
  dis_ele_pos_uni->Draw("colz");

  canvas->cd(8);
  dis_ele_pos_exp->SetTitle("Dispersao do decaimento no lab (Metodo Exponencial)");
  dis_ele_pos_exp->GetXaxis()->SetTitle("Momento Eletron (GeV/c)");
  dis_ele_pos_exp->GetYaxis()->SetTitle("Momento Positron (GeV/c)");
  dis_ele_pos_exp->Draw("colz");
}

// Simulacoes com sementes diferentes (duas) e para cada evento somar as cargas de todas as particulas de estado final e fazer um histograma com elas ---> valor esperado e +2
// (dois protons iniciais), i.e., esperamos uma barra vertical no 2