//
// Created by _amp_ on 10/21/24.
//

#include "MSTClustering.h"

void MSTClustering::construct_tree(std::vector<edge>&& verticeData, size_t size) {
  tree.clear();
  for (int i = 0; i < size; i++) tree.push_back(std::make_shared<Node>(i));
  std::sort(verticeData.begin(), verticeData.end());
  for (auto it : verticeData) {
    uint v1 = it.second.first;
    uint v2 = it.second.second;
    if (tree[v1].get() != tree[v2].get()) {
      auto newNode = std::make_shared<Node>(std::move(tree[v1]), std::move(tree[v2]), it.first);
      for (const auto vertice : newNode->get_vertices()) tree[vertice] = newNode;
    }
  }
}

std::vector<edge> GMSTClustering::get_vertices(const cola::EventData& edata) {
  cola::EventParticles particles = edata.particles;
  std::vector<edge> edges;
  for (size_t i = 0; i < particles.size(); i++) {
    for (size_t j = i + 1; j < particles.size(); i++) {
      double temp_dist = 0.0;
      temp_dist += std::pow(particles[i].position.x - particles[j].position.x, 2);
      temp_dist += std::pow(particles[i].position.y - particles[j].position.y, 2);
      temp_dist += std::pow(particles[i].position.z - particles[j].position.z, 2);
      temp_dist = std::sqrt(temp_dist) edges.emplace_back({temp_dist, {i, j}});
      edges.emplace_back({temp_dist, {j, i}});
    }
  }
  return edges;
}

std::unique_ptr<cola::EventData> GMSTClustering::get_clusters(std::unique_ptr<cola::EventData>&& edata, const Node& rootnode) {
  cola::EventParticles nucleons;
  cola::EventParticles nucleons_B;

  cola::EventIniState inidata = edata.iniState;
  cola::EventParticles particles = edata.particles;

  sourceA = pdgToAZ(inidata.pdgCodeA).first;
  sourceAb = pdgToAZ(inidata.pdgCodeB).first;
  A = 0;
  Ab = 0;
  Z = 0;
  Zb = 0;

  for (auto& particle : particles) {
    auto az = particle.getAZ();
    if (particle.pClass.spectatorA) {
      nucleons.push_back(particle);
      ++A;
      if (az.second == 1) {
        ++Z;
      }
    }
    if (particle.pClass.spectatorB) {
      nucleons_B.push_back(particle);
      ++Ab;
      if (az.second == 1) {
        ++Zb;
      }
    }
  }

  auto ExEns = this->get_exens(4);  // Setting default stat type 4
  double ExA = ExEns.first;
  double ExB = ExEns.second;

  std::vector<G4FragmentVector> outClusters;
  std::vector<std::vector<LorentzVector>> positions;

  std::vector<std::vector<uint>> clusters = this->get_connected_components(this->get_cd(ExA, A));
  std::vector<std::vector<uint>> clusters_B = this->get_connected_components(this->get_cd(ExB, Ab));

  auto [output_vector_A, positionsA] = this->fragments_from_clusters(clusters, nucleons);
  auto [output_vector_B, positionsB] = this->fragments_from_clusters(clusters_B, nucleons_B);

  outClusters.push_back(output_vector_A);
  outClusters.push_back(output_vector_B);

  positions.push_back(positionsA);
  positions.push_back(positionsB);

  auto [boostA, boostB] = this->get_boosts();

  cola::EventParticles cFragments = this->calculate_momentum(outClusters, ExA, ExB, boostA, boostB, positions);

  delete nucleons;
  delete nucleons_B;

  return std::make_unique<cola::EventData>(cola::EventData{cola::EventIniState{inistate.pdgCodeA, inistate.pdgCodeB, inistate.pZA, inistate.pZB, inistate.energy, inistate.sectNN, inistate.b, inistate.nColl, initstate.nCollP, inistate.nCollPN, inistate.nCollNN, inistate.nPart, inistate.nPartA, inistate.nPartB, inistate.phiRotA, inistate.thetaRotA, inistate.phiRotB, inistate.thetaRotB, cFragments}, cFragments});
}

EventParticles GMSTClustering::calculate_momentum(std::vector<G4FragmentVector> noMomClusters, double ExEnA, double ExEnB, CLHEP::Hep3Vector boostA, CLHEP::Hep3Vector boostB, std::vector<std::vector<CLHEP::Hep3Vector>> positions) {
  cola::EventParticles particles;

  std::vector<double> MstMassVector_A;
  MstMassVector_A.reserve(noMomClusters.at(0).size());

  std::vector<double> MstMassVector_B;
  MstMassVector_B.reserve(noMomClusters.at(1).size());

  double SumMassMst = 0;
  double SumMassMstEx = 0;

  if (!noMomClusters.at(0).empty()) {
    for (uint i = 0; i < noMomClusters.at(0).size(); ++i) {
      uint clfrag_A = (noMomClusters.at(0)[i])->GetA();
      uint clfrag_Z = (noMomClusters.at(0)[i])->GetZ();

      double energy = 0;
      if (!((clfrag_Z == 0 && (clfrag_A == 1)) || (clfrag_Z == 1 && (clfrag_A == 1)))) {
        energy = ExEnA * double(clfrag_A) / double(SpecAa);
      }

      double NuclearMass = G4NucleiProperties::GetNuclearMass(clfrag_A, clfrag_Z) + energy;

      SumMassMst += G4NucleiProperties::GetNuclearMass(clfrag_A, clfrag_Z);

      SumMassMstEx += NuclearMass;

      MstMassVector_A.push_back(NuclearMass);
    }

    double PrefragmentMass_A = SumMassMst + ExEnA;

    std::vector<cola::LorentzVector*>* momentumVectorA;

    if (PrefragmentMass_A < (SumMassMstEx + 1e-5 * MeV)) {
      PrefragmentMass_A += 1e-5 * MeV;
    }
    momentumVectorA = phaseSpaceDecay.Decay(PrefragmentMass_A, MstMassVector_A);

    for (int I = 0; I < noMomClusters.at(0).size(); ++I) {
      auto fragment = noMomClusters.at(0).at(I);
      fragment->SetMomentum((*momentumVectorA->at(I)).boost(boostA.x(), boostA.y(), boostA.z()));

      cola::Particle particle = this->fragment_to_particle(fragment);
      particle.push_back(particle);
    }
    momentumVectorA->clear();
  }

  // side B

  SumMassMst = 0;
  SumMassMstEx = 0;
  if (!noMomClusters.at(1).empty()) {
    for (uint i = 0; i < noMomClusters.at(1).size(); ++i) {
      uint clfrag_A = (noMomClusters.at(1)[i])->GetA();
      uint clfrag_Z = (noMomClusters.at(1)[i])->GetZ();
      double energy = 0;
      if (!((clfrag_Z == 0 && (clfrag_A == 1)) || (clfrag_Z == 1 && (clfrag_A == 1)))) {
        energy = ExEnB * double(clfrag_A) / double(SpecAb);
      }

      double NuclearMass = G4NucleiProperties::GetNuclearMass(clfrag_A, clfrag_Z) + energy;

      SumMassMst += G4NucleiProperties::GetNuclearMass(clfrag_A, clfrag_Z);

      SumMassMstEx += NuclearMass;

      MstMassVector_B.push_back(NuclearMass);
    }

    double PrefragmentMass_B = SumMassMst + ExEnB;
    std::vector<cola::LorentzVector*>* momentumVectorB;

    if (PrefragmentMass_B < (SumMassMstEx + 1e-5 * MeV)) {
      PrefragmentMass_B += 1e-5 * MeV;
    }
    momentumVectorB = phaseSpaceDecay.Decay(PrefragmentMass_B, MstMassVector_B);

    for (int I = 0; I < noMomClusters.at(1).size(); ++I) {
      auto fragment = noMomClusters.at(1).at(I);
      fragment->SetMomentum((*momentumVectorA->at(I)).boost(boostB.x(), boostB.y(), boostB.z()));
      z
      cola::Particle particle = this->fragment_to_particle(fragment);
      particle.push_back(particle);
    }

    momentumVectorB->clear();
  }

  MstMassVector_A.clear();
  MstMassVector_B.clear();
  noMomClusters.clear();

  return particles;
}

// Get Critical Distance for 1 nuclei
double GMSTClustering::get_cd(double Ex, uint A) {
  if (Ex / double(A) < 2.17 * MeV) {
    return d0;
  }
  double ex = Ex / double(A);
  double dep = std::exp(-std::pow((ex / b_opt), a_opt)) * c_opt + d_opt;
  return d0 * std::pow(dep, 1. / 3.);
}

// Get fragments and their positions from clusters after Krusal Algorithm for 1 nuclei
std::pair<G4FragmentVector, std::vector<std::vector<cola::LorentzVector>>> GMSTClustering::fragments_from_clusters(std::vector<std::vector<uint>> clusters, cola::EventParticles nucleons) {
  G4FragmentVector fragments;
  std::vector<cola::LorentzVector>> positions;

  for (uint i = 0; i < clusters.size(); ++i) {
    cola::LorentzVector position;
    uint Z_clust = 0;
    uint A_clust = 0;
    for (uint j = 0; j < clusters[i].size(); ++j) {
      cola::Particle* nucleon = &(nucleons.at(clusters[i][j]));
      position += nucleon->position;

      if (nucleon.getAZ().second == 1) {
        Z_clust += 1;
      }
      A_clust += 1;
    }

    position /= double(A_clust);

    G4Fragment* frag = new G4Fragment();
    frag->SetA(A_clust);
    frag->SetZ(Z_clust);
    fragments.push_back(frag);
    positions.push_back(position);
  }
  return std::make_pair<fragments, positions>;
}

// Get boost vectors for 1 nuclei
std::pair<CLHEP::Hep3Vector, CLHEP::Hep3Vector> GMSTClustering::get_boosts() {
  // Need arguments
  
  // std::string model_type = "G";
  // if (histoManager.GetReaderID()) model_type = "VAL";
  // FermiMomentum FermiMom(ain, model_type);

  // FermiMom.SetPzPerNucleon(histoManager.GetInitialContidions().GetPzA() / sourceA, histoManager.GetInitialContidions().GetPzB() / sourceAb);

  // CLHEP::HepLorentzVector Fermi4MomA = FermiMom.GetLorentzVector("A");
  // CLHEP::Hep3Vector boostA = Fermi4MomA.boostVector();
  // event.FermiMomA_x = Fermi4MomA.px();
  // event.FermiMomA_y = Fermi4MomA.py();
  // event.FermiMomA_z = Fermi4MomA.pz();
  // CLHEP::HepLorentzVector Fermi4MomB = FermiMom.GetLorentzVector("B");
  // CLHEP::Hep3Vector boostB = Fermi4MomB.boostVector();
  // event.FermiMomB_x = Fermi4MomB.px();
  // event.FermiMomB_y = Fermi4MomB.py();
  // event.FermiMomB_z = Fermi4MomB.pz();
  // return std::make_pair(boostA, boostB);
}

// DFS for Kruskal Algorithm
void GMSTClustering::dfs(std::shared_ptr<Node> node, std::unordered_set<uint>& visited, std::vector<uint>& component) {
  if (!node || visited.count(node->get_vertices().front())) return;

  for (const auto& vertice : node->get_vertices()) {
    visited.insert(vertice);
    component.push_back(vertice);
  }
  auto [left, right] = node->get_children();
  dfs(left, visited, component);
  dfs(right, visited, component);
}

std::vector<std::vector<uint>> GMSTClustering::get_connected_components() {
  std::vector<std::vector<uint>> components;
  std::unordered_set<uint> visited;

  for (const auto& node : tree) {
    if (!node || visited.count(node->get_vertices().front())) continue;

    std::vector<uint> component;
    dfs(node, visited, component);
    components.push_back(component);
  }

  return components;
}

// Get Excitation energies for projectile and target
std::pair<double, double> GMSTClustering::get_exens(uint statType) {
  auto ExEnA = std::make_unique<ExcitationEnergy>(statType, sourceA);
  auto ExEnB = std::make_unique<ExcitationEnergy>(statType, sourceAb);

  if (statType > 2) {
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
  double ExEnA = energy_A / double(A);
  double ExEnB = energy_B / double(Ab);

  delete ExEnA;
  delete ExEnB;

  return std::make_pair(energy_A, energy_B);
}

cola::Particle GMSTClustering::fragment_to_particle(G4Fragment fragment) {
  cola::Particle particle;
  particle.position = *(positions.at(0).at(I));
  particle.momentum = fragment->GetMomentum();
  particle.pdgCode = AZToPdg(std::make_pair(fragment->GetA_asInt(), fragment->GetZ_asInt()));
  particle.pClass = cola::ParticleClass::produced;
  particles.push_back(particle);
  return particle;
}