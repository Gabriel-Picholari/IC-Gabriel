#include "TSystem.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TPythia8.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"


void fastJet(Int_t nev = 10, Int_t ndeb = 1 /* Listagem */ )
{
  //---------------------------------------------------------------------------------------------------------
  // Inicializacoes e configuracoes do Pythia:
  //---------------------------------------------------------------------------------------------------------

  TPythia8 *pythia8 = new TPythia8();
  pythia8->ReadString("HardQCD:all = on");
  pythia8->ReadString("Random:setSeed = on");
  pythia8->ReadString("Random:seed = 43");


  pythia8->Initialize(2212 /* Proton */, 2212 /* Proton */, 14000 /* TeV */); /* 14000 TeV = 14000000 GeV */

  TClonesArray *particles = new TClonesArray("TParticle", 1000);

  //---------------------------------------------------------------------------------------------------------
  // Inicializacoes e configuracoes do FastJet:
  //---------------------------------------------------------------------------------------------------------

  Double_t R = 0.4;

  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);

  std::vector<fastjet::PseudoJet> particles_fastjet;
  std::vector<fastjet::PseudoJet> jets;

  //---------------------------------------------------------------------------------------------------------

  // Loop de eventos
  for (Int_t iev = 0; iev < nev; iev++)
  {
    pythia8->GenerateEvent();
    if (iev < ndeb) pythia8->EventListing();
    pythia8->ImportParticles(particles, "All");
    Int_t np = particles->GetEntriesFast();

    particles_fastjet.clear();

    //Loop de particulas
    for (Int_t ip = 0; ip < np; ip++)
    {
      TParticle *part = (TParticle*) particles->At(ip);

      Int_t situacao = part->GetStatusCode();

      if (situacao > 0)
      {
        double px = part->Px();
        double py = part->Py();
        double pz = part->Pz();
        double e = part->Energy();

        fastjet::PseudoJet particle(px, py, pz, e);
        particles_fastjet.push_back(particle);
      }
    }
    // Agrupamento de jatos segundo  jet_def definido anteriormente e os jatos obtidos e armazenados em particles_fastjet
    fastjet::ClusterSequence clusterSeq(particles_fastjet, jet_def);
    jets = clusterSeq.inclusive_jets();

    // Loop sobre os jatos e imprima informações
    for (const fastjet::PseudoJet& jet : jets) {
    std::cout << "Jet pT: " << jet.pt() << " GeV/c" << std::endl;
    std::cout << "Jet rapidity: " << jet.rap() << std::endl;
    std::cout << "Jet phi: " << jet.phi() << std::endl;
  }
}
}