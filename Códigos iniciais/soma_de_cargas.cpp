#include "TSystem.h"
#include "TH1F.h"
#include "TPythia8.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"

void soma_de_cargas(Int_t nev = 10, Int_t ndeb = 1)
{
  // Importamos bibliotecas para a melhor visualizacao do canvas e integracao com o pythia8
  gSystem->Load("libEG");
  gSystem->Load("libEGPythia8");

  // Criamos um objeto pythia8 que estara sujeito a receber metodos da biblioteca referente ao pythia
  TPythia8 *pythia8_1 = new TPythia8();
  TPythia8 *pythia8_2 = new TPythia8();

  // Habilitamos interacoes entre forcas fortes e uma semente para um evento de modo que possa ser reproduzido
  pythia8_1->ReadString("HardQCD:all = on");
  pythia8_1->ReadString("Random:setSeed = on");
  pythia8_1->ReadString("Random:seed = 42");

  pythia8_2->ReadString("HardQCD:all = on");
  pythia8_2->ReadString("Random:setSeed = on");
  pythia8_2->ReadString("Random:seed = 84");

  // Inicializamos as particulas que comporao os eventos, definindo-as, bem como sua energia em GeV
  pythia8_1->Initialize(2212, 2212, 14000);
  pythia8_2->Initialize(2212, 2212, 14000);

  // Criamos um array de armazenamento para as particulas resultantes dos eventos
  TClonesArray *particles = new TClonesArray("TParticle", 10000);
  TClonesArray *particles_2 = new TClonesArray("TParticle", 10000);

  // Criamos uma variavel de interesse
  Float_t sum_of_charges = 0;
  Float_t sum_of_charges_2 = 0;

  // Criamos um histograma de interesse
  TH1F *hist_sum_of_charges = new TH1F("hist_sum_of_charges", "Soma das cargas das particulas do evento (Seed 42)", 100, 0, 0);
  TH1F *hist_charges = new TH1F("hist_charges", "Cargas das particulas finais do evento (Seed 42)", 100, 0, 0);

  TH1F *hist_sum_of_charges_2 = new TH1F("hist_sum_of_charges_2", "Soma das cargas das particulas do evento (Seed 84)", 100, 0, 0);
  TH1F *hist_charges_2 = new TH1F("hist_charges_2", "Cargas das particulas finais do evento (Seed 84)", 100, 0, 0);


  // Loop gerador de eventos (comum a ambas as seeds)
  for (Int_t iev = 0; iev < nev; iev++)
  {
    pythia8_1->GenerateEvent(); // Geramos o evento de indice iev, resultado possivel de uma colisao p-p
    pythia8_2->GenerateEvent();
    sum_of_charges_2 = 0;
    sum_of_charges = 0;

    if (iev < ndeb) // Listaremos detalhadamente apenas o primeiro evento
    {
      pythia8_1->EventListing();
      pythia8_2->EventListing();
    }

    pythia8_1->ImportParticles(particles, "All"); // Importamos todas as particulas criadas no evento para o array criado
    pythia8_2->ImportParticles(particles_2, "All");

    Int_t np = particles->GetEntriesFast(); // Obtemos o numero efetivo de particulas geradas no evento - 10000 e apenas um espaco predefinido e flexivel
    Int_t np_2 = particles_2->GetEntriesFast();

    // Loop leitura de particulas para a seed 42
    for (Int_t ip = 0; ip < np; ip++)
    {
      TParticle *part = (TParticle*) particles->At(ip);
      Int_t status = part->GetStatusCode(); // Tomamos o codigo da particula

      Int_t pdg = part->GetPdgCode(); // Usando o ponteiro criado, obtemos o codigo PDG da particula de indice ip
      Float_t charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge(); //Dividir por 3

      if (pdg >= 1 && pdg < 6) printf ("PDG=%d, charge=%1.1f\n",pdg, charge);
      // Codigos positivos indicam particulas finais, que sao as que nos interessam
      if (status > 0)
      {
        sum_of_charges = sum_of_charges + charge;
        hist_charges->Fill(charge);
        // printf ("PDG=%d, charge=%1.1f\n",pdg, charge);
      }
    }

    // Loop leitura de particulas para a seed 84
    for (Int_t ip_2 = 0; ip_2 < np_2; ip_2++)
    {
      TParticle *part_2 = (TParticle*) particles_2->At(ip_2);
      Int_t stats_2 = part_2->GetStatusCode();

      if (stats_2 > 0)
      {
        Int_t pdg_2 = part_2->GetPdgCode();
        Float_t charge_2 = TDatabasePDG::Instance()->GetParticle(pdg_2)->Charge();
        sum_of_charges_2 = sum_of_charges_2 + charge_2;
        hist_charges_2->Fill(charge_2);
      }
    }
  }
  hist_sum_of_charges->Fill(sum_of_charges);
  hist_sum_of_charges_2->Fill(sum_of_charges_2);


  TCanvas *c1 = new TCanvas("c1", "Soma de cargas de particulas do evento", 2500, 2500);
  c1->Divide(2, 2);

  c1->cd(1);
  hist_sum_of_charges->SetTitle("Soma das cargas das particulas finais da colisao (Seed 42)");
  hist_sum_of_charges->GetXaxis()->SetTitle("Soma (Seed 42)");
  hist_sum_of_charges->GetYaxis()->SetTitle("Frequencia");
  hist_sum_of_charges->Draw();

  c1->cd(2);
  hist_charges->SetTitle("Cargas das particulas finais da colisao (Seed 42)");
  hist_charges->GetXaxis()->SetTitle("Cargas (Seed 42)");
  hist_charges->GetYaxis()->SetTitle("Frequencia");
  hist_charges->Draw();

  c1->cd(3);
  hist_sum_of_charges_2->SetTitle("Soma das cargas das particulas finais da colisao (Seed 84)");
  hist_sum_of_charges_2->GetXaxis()->SetTitle("Soma (Seed 84)");
  hist_sum_of_charges_2->GetYaxis()->SetTitle("Frequencia");
  hist_sum_of_charges_2->Draw();

  c1->cd(4);
  hist_charges_2->SetTitle("Cargas das particulas finais da colisao (Seed 84)");
  hist_charges_2->GetXaxis()->SetTitle("Cargas (Seed 84)");
  hist_charges_2->GetYaxis()->SetTitle("Frequencia");
  hist_charges_2->Draw();
}