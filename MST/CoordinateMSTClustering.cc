#include "CoordinateMSTClustering.hh"
#include "Repulsion.hh"
#include "ExcitationEnergy.hh"
#include <limits>
#include <stack>

std::unique_ptr<cola::EventData> CoordinateMSTClustering::get_clusters(std::unique_ptr<cola::EventData>&& edata, const Node& rootnode) {
  cola::EventParticles nucleons;
  cola::EventParticles nucleons_B;

  cola::EventIniState inistate = (*edata).iniState;
  double pZA = inistate.pZA;
  double pZB = inistate.pZB;
  double Ea;
  double Eb;
  cola::EventParticles particles = (*edata).particles;
  pls_ = particles;

  sourceA = cola::pdgToAZ(inistate.pdgCodeA).first;
  sourceAb = cola::pdgToAZ(inistate.pdgCodeB).first;
  A = 0;
  Ab = 0;
  Z = 0;
  Zb = 0;

  cola::LorentzVector pA = {0.0, 0.0, 0.0, 0.0};
  cola::LorentzVector pB = {0.0, 0.0, 0.0, 0.0};
  for (auto& particle : particles) {
    auto az = particle.getAZ();
    if (particle.pClass == cola::ParticleClass::spectatorA) {
      pA += particle.momentum;
      Ea += particle.momentum.e;
      ++A;
      if (az.second == 1) {
        ++Z;
      }
    }
    if (particle.pClass == cola::ParticleClass::spectatorB) {
      pB += particle.momentum;
      Eb += particle.momentum.e;
      ++Ab;
      if (az.second == 1) {
        ++Zb;
      }
    }
  }
  // std::cout << pA.z / A << " " << pZA << "\n";
  auto boostA = this->get_boost(CLHEP::Hep3Vector(pA.x, pA.y, pA.z), Ea, static_cast<double>(A));
  auto boostB = this->get_boost(CLHEP::Hep3Vector(pB.x, pB.y, pB.z), Eb, static_cast<double>(Ab));
  // std::cout << boostA.mag() << "\n";
  for (auto& particle : particles) {
    if (particle.pClass == cola::ParticleClass::spectatorA) {
      particle.momentum = particle.momentum.boost(-boostA.x(), -boostA.y(), -boostA.z());
      nucleons.push_back(particle);
    }
    if (particle.pClass == cola::ParticleClass::spectatorB) {
      particle.momentum = particle.momentum.boost(-boostB.x(), -boostB.y(), -boostB.z());
      nucleons_B.push_back(particle);
    }
  }

  auto ExEns = this->get_exens();
  double ExA = ExEns.first;
  double ExB = ExEns.second;

  std::vector<std::vector<cola::Particle*>> outClusters;
  std::vector<cola::Particle*> output_vector_A;
  std::vector<cola::Particle*> output_vector_B;
  std::vector<int> rmapsA;
  std::vector<int> rmapsB;
  int rcountA = 0;
  int rcountB = 0;
  cola::EventParticles rnucsA;
  cola::EventParticles rnucsB;

  std::vector<std::vector<uint>> clusters = this->get_comps(this->get_cd(ExA, A), cola::ParticleClass::spectatorA);
  std::vector<std::vector<uint>> clusters_B = this->get_comps(this->get_cd(ExB, Ab), cola::ParticleClass::spectatorB);

  output_vector_A = this->fragments_from_clusters(clusters);
  output_vector_B = this->fragments_from_clusters(clusters_B);

  outClusters.push_back(output_vector_A);
  outClusters.push_back(output_vector_B);


  for (uint i = 0; i < clusters.size(); ++i) {
    for (uint j = 0; j < clusters[i].size(); ++j) {
      cola::Particle* nucleon = &(pls_.at((clusters[i])[j]));
      if (nucleon->getAZ().second == 1) {
        rnucsA.push_back(*nucleon);
        rmapsA.push_back(rcountA);
      }
    }
    ++rcountA;
  }

  for (uint i = 0; i < clusters_B.size(); ++i) {
    for (uint j = 0; j < clusters_B[i].size(); ++j) {
      cola::Particle* nucleon = &(pls_.at((clusters_B[i])[j]));
      if (nucleon->getAZ().second == 1) {
        rnucsB.push_back(*nucleon);
        rmapsB.push_back(rcountB);
      }
    }
    ++rcountB;
  }

  cola::EventParticles cFragments = this->calculate_momentum(outClusters, ExA, ExB, boostA, boostB, rnucsA, rnucsB, rmapsA, rmapsB);
  return std::make_unique<cola::EventData>(cola::EventData{cola::EventIniState{inistate.pdgCodeA, inistate.pdgCodeB, inistate.pZA, inistate.pZB, inistate.energy, inistate.sectNN, inistate.b, inistate.nColl, inistate.nCollPP, inistate.nCollPN, inistate.nCollNN, inistate.nPart, inistate.nPartA, inistate.nPartB, inistate.phiRotA, inistate.thetaRotA, inistate.phiRotB, inistate.thetaRotB, cFragments}, cFragments});
}

std::vector<MSTClustering::edge> CoordinateMSTClustering::get_edges(const cola::EventData& edata) {
  cola::EventParticles particles = edata.particles;
  std::vector<edge> edges;
  for (size_t i = 0; i < particles.size(); i++) {
    for (size_t j = i + 1; j < particles.size(); j++) {
      double temp_dist = std::numeric_limits<double>::infinity();
      if (particles[i].pClass == particles[j].pClass) {
        temp_dist = std::pow(particles[i].position.x - particles[j].position.x, 2);
        temp_dist += std::pow(particles[i].position.y - particles[j].position.y, 2);
        temp_dist += std::pow(particles[i].position.z - particles[j].position.z, 2);
        temp_dist = std::sqrt(temp_dist);
      }
      edges.emplace_back(temp_dist, std::make_pair(i, j));
      edges.emplace_back(temp_dist, std::make_pair(j, i));
    }
  }
  return edges;
}

cola::EventParticles CoordinateMSTClustering::calculate_momentum(std::vector<std::vector<cola::Particle*>> noMomClusters, double ExEnA, double ExEnB, CLHEP::Hep3Vector boostA, CLHEP::Hep3Vector boostB, cola::EventParticles rnucsA, cola::EventParticles rnucsB, std::vector<int> rmapsA, std::vector<int> rmapsB) {
  cola::EventParticles particles;

  auto momClusters = noMomClusters;
  std::vector<double> MstMassVector_A;
  MstMassVector_A.reserve(noMomClusters.at(0).size());

  std::vector<double> MstMassVector_B;
  MstMassVector_B.reserve(noMomClusters.at(1).size());

  double SumMassMst = 0;
  double SumMassMstEx = 0;

  if (!noMomClusters.at(0).empty()) {
    for (uint i = 0; i < noMomClusters.at(0).size(); ++i) {
      uint clfrag_A = (noMomClusters.at(0)[i])->getAZ().first;
      uint clfrag_Z = (noMomClusters.at(0)[i])->getAZ().second;

      double energy = 0;
      if (!((clfrag_Z == 0 && (clfrag_A == 1)) || (clfrag_Z == 1 && (clfrag_A == 1)))) {
        energy = ExEnA * double(clfrag_A) / double(SpecAa);
      }

      double NuclearMass = G4NucleiProperties::GetNuclearMass(static_cast<G4int>(clfrag_A), static_cast<G4int>(clfrag_Z)) + energy;

      SumMassMst += G4NucleiProperties::GetNuclearMass(static_cast<G4int>(clfrag_A), static_cast<G4int>(clfrag_Z));

      SumMassMstEx += NuclearMass;

      MstMassVector_A.push_back(NuclearMass);
    }

    double PrefragmentMass_A = SumMassMst + ExEnA;

    std::vector<G4LorentzVector*>* momentumVectorA;

    if (PrefragmentMass_A < (SumMassMstEx + 1e-5 * MeV)) {
      PrefragmentMass_A += 1e-5 * MeV;
    }
    momentumVectorA = phaseSpaceDecay.Decay(PrefragmentMass_A, MstMassVector_A);

    for (int I = 0; I < momClusters.at(0).size(); ++I) {
      momClusters.at(0).at(I)->momentum = ToColaLorentzVector(*(momentumVectorA->at(I)));
    }

    momentumVectorA->clear();

    if (consider_rep_) {
      momClusters.at(0) = RepulsionStage::CalculateRepulsion(momClusters.at(0), rnucsA, rmapsA);
    }
    for (int I = 0; I < momClusters.at(0).size(); ++I) {
      auto fragment = momClusters.at(0).at(I);
      cola::LorentzVector momentum = fragment->momentum;
      momentum.boost(boostA.x(), boostA.y(), boostA.z());
      fragment->momentum = momentum;
      fragment->pClass = cola::ParticleClass::spectatorA;

      particles.push_back(*fragment);
    }
  }

  SumMassMst = 0;
  SumMassMstEx = 0;
  if (!noMomClusters.at(1).empty()) {
    for (uint i = 0; i < noMomClusters.at(1).size(); ++i) {
      uint clfrag_A = (noMomClusters.at(1)[i])->getAZ().first;
      uint clfrag_Z = (noMomClusters.at(1)[i])->getAZ().second;

      double energy = 0;
      if (!((clfrag_Z == 0 && (clfrag_A == 1)) || (clfrag_Z == 1 && (clfrag_A == 1)))) {
        energy = ExEnB * double(clfrag_A) / double(SpecAb);
      }

      double NuclearMass = G4NucleiProperties::GetNuclearMass(static_cast<G4int>(clfrag_A), static_cast<G4int>(clfrag_Z)) + energy;

      SumMassMst += G4NucleiProperties::GetNuclearMass(static_cast<G4int>(clfrag_A), static_cast<G4int>(clfrag_Z));

      SumMassMstEx += NuclearMass;

      MstMassVector_B.push_back(NuclearMass);
    }

    double PrefragmentMass_B = SumMassMst + ExEnB;
    std::vector<G4LorentzVector*> *momentumVectorB;

    if (PrefragmentMass_B < (SumMassMstEx + 1e-5 * MeV)) {
      PrefragmentMass_B += 1e-5 * MeV;
    }
    momentumVectorB = phaseSpaceDecay.Decay(PrefragmentMass_B, MstMassVector_B);


    for (int I = 0; I < momClusters.at(1).size(); ++I) {
      momClusters.at(1).at(I)->momentum = ToColaLorentzVector(*momentumVectorB->at(I));
    }

    momentumVectorB->clear();

    if (consider_rep_) {
      momClusters.at(1) = RepulsionStage::CalculateRepulsion(momClusters.at(1), rnucsB, rmapsB);
    }

    for (int I = 0; I < momClusters.at(1).size(); ++I) {
      auto fragment = momClusters.at(1).at(I);
      cola::LorentzVector momentum = fragment->momentum;
      momentum.boost(boostB.x(), boostB.y(), boostB.z());
      fragment->momentum = momentum;
      fragment->pClass = cola::ParticleClass::spectatorB;

      particles.push_back(*fragment);
    }
  }

  MstMassVector_A.clear();
  MstMassVector_B.clear();
  noMomClusters.clear();
  momClusters.clear();

  return particles;
}

double CoordinateMSTClustering::get_cd(double Ex, uint A) {
  if ((Ex / double(A)) < 2.17 * MeV) {
    return d0;
  }
  double ex = Ex / double(A);
  double dep = std::exp(-std::pow((ex / b_opt), a_opt)) * c_opt + d_opt;
  return d0 * std::pow(dep, 1. / 3.);
}

std::vector<cola::Particle*> CoordinateMSTClustering::fragments_from_clusters(const std::vector<std::vector<uint>>& clusters) {
  std::vector<cola::Particle*> fragments;

  for (uint i = 0; i < clusters.size(); ++i) {
    uint Z_clust = 0;
    uint A_clust = 0;
    cola::LorentzVector position = {0.0, 0.0, 0.0, 0.0};
    for (uint j = 0; j < clusters[i].size(); ++j) {
      auto nucleon = &(pls_.at(clusters[i][j]));
      position += nucleon->position;
      if (nucleon->getAZ().second == 1) {
        Z_clust += 1;
      }
      A_clust += 1;
    }

    cola::Particle* frag = new cola::Particle();
    frag->pdgCode = cola::AZToPdg({A_clust, Z_clust});
    frag->position = (position / double(A_clust));
    fragments.push_back(frag);
  }
  return fragments;
}

std::vector<std::vector<uint>> CoordinateMSTClustering::get_comps(double cd, cola::ParticleClass pClass) {
  std::vector<std::vector<uint>> components;
  std::vector<bool> visited(tree.size(), false);

  for (size_t i = 0; i < tree.size(); ++i) {
    if (!visited[i] && (pls_[i].pClass == pClass)) {
      std::vector<uint> component;
      dfs(tree[i], visited, component, cd, pClass);
      components.push_back(component);
    }
  }
  return components;
}

void CoordinateMSTClustering::dfs(std::shared_ptr<Node> node, std::vector<bool>& visited, std::vector<uint>& component, double cd, cola::ParticleClass pClass) {
  std::stack<std::shared_ptr<Node>> stack;
  stack.push(node);

  while (!stack.empty()) {
    auto current = stack.top();
    stack.pop();

    for (uint v : current->get_vertices()) {
      if (!visited[v] && (pls_[v].pClass == pClass)) {
        visited[v] = true;
        component.push_back(v);
      }
    }

    auto children = current->get_children();
    if (children.first && children.first->get_height() <= cd) {
      stack.push(children.first);
    }
    if (children.second && children.second->get_height() <= cd) {
      stack.push(children.second);
    }
  }
}

std::pair<double, double> CoordinateMSTClustering::get_exens() {
  auto ExEnA = std::make_unique<ExcitationEnergy>(stat_exen_type_, sourceA);
  auto ExEnB = std::make_unique<ExcitationEnergy>(stat_exen_type_, sourceAb);

  if (stat_exen_type_ > 2) {
    double e_0 = 11.5 * MeV;
    double sigma0 = 0.005;
    double b0 = 2;
    double sigmaE0 = 1 * MeV;
    double c0 = 0.1;
    double Pe = 24 * MeV;
    double Pm = 0.2;

    ExEnA->SetParametersALADIN(e_0, sigma0, b0);
    ExEnB->SetParametersALADIN(e_0, sigma0, b0);
    ExEnA->SetParametersParabolicApproximation(Pe, Pm, sigma0, c0, 0.01);
    ExEnB->SetParametersParabolicApproximation(Pe, Pm, sigma0, c0, 0.01);

    ExEnA->SetParametersHybridFit(11.46648905 * MeV, -1.84830078 * MeV, -58.53674677 * MeV, 284.66431513 * MeV, -637.51406293 * MeV, 652.80324427 * MeV, -251.28205381 * MeV, 0.4 * MeV, 0.5, 0.2);
    ExEnB->SetParametersHybridFit(11.46648905 * MeV, -1.84830078 * MeV, -58.53674677 * MeV, 284.66431513 * MeV, -637.51406293 * MeV, 652.80324427 * MeV, -251.28205381 * MeV, 0.4 * MeV, 0.5, 0.2);
  }

  double energy_A = ExEnA->GetEnergy(A);
  double energy_B = ExEnB->GetEnergy(Ab);

  return std::make_pair(energy_A, energy_B);
}

CLHEP::Hep3Vector CoordinateMSTClustering::get_boost(CLHEP::Hep3Vector p, double E, double A) {
  if(A != 0.0) {
    // std::cout << p.x() / A << "\n";
    CLHEP::HepLorentzVector futureBoost(p.x(), p.y(), p.z(), E);
    if(!futureBoost.isTimelike()) {
      std::cout << "Not TimeLike Hep3LorentsVector at " << E << " nucl = " << A << std::endl;
    }
    // std::cout << p.x() / A << "\n";
    // std::cout << futureBoost.boostVector().mag() << " " << p.x() / A << " " << p.z() / A << "\n";
    // std::cout << E << "\n";
    return futureBoost.boostVector();
  }
  else {
    return CLHEP::Hep3Vector(0.0, 0.0, 0.0);
  }
}

cola::LorentzVector CoordinateMSTClustering::ToColaLorentzVector(G4LorentzVector lv) {
  cola::LorentzVector vec;
  vec.e = lv.e();
  vec.x = lv.x();
  vec.y = lv.y();
  vec.z = lv.z();
  return vec;
}
