/**
 * Copyright (c) 2026 Savva Savenkov, Artemii Novikov
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
 * associated documentation files (the "Software"), to deal in the Software without restriction,
 * including without limitation the rights to use, copy, modify, merge, publish, distribute,
 * sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies or
 * substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
 * NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "CoordinateMSTClustering.hh"

#include "ExcitationEnergy.hh"
#include "Repulsion.hh"

#include <G4SystemOfUnits.hh>

#include <limits>
#include <queue>

// constants for implementation
constexpr double eps0 = 2.17 * MeV;
constexpr double alphaPow = -1.02;
constexpr double d0 = 2.7;
constexpr double aColRel = 5.315;
constexpr double aSurfRel = 0.017;
constexpr double aVol = 0.054;

constexpr double SpecAa = 0;
constexpr double SpecAb = 0;

constexpr double a_opt = 2.243;
constexpr double b_opt = 3.183 * MeV;
constexpr double c_opt = 0.99;
constexpr double d_opt = 0.29041;

using namespace cola;

std::unique_ptr<cola::EventData> CoordinateMSTClustering::get_clusters(std::unique_ptr<cola::EventData>&& data) {
  cola::EventParticles clustersA;
  cola::EventParticles clustersB;

  // get clusters
  if (rootA != nullptr) {
    clustersA = _process_side(*data, cola::ParticleClass::kSpectatorA);
  }
  if (rootB != nullptr) {
    clustersB = _process_side(*data, cola::ParticleClass::kSpectatorB);
  }

  // erase kSpectator nucleons
  data->particles.erase(spectIterA != endIter ? spectIterA : spectIterB, endIter);

  // append clusters
  if (rootA != nullptr) {
    data->particles.insert(data->particles.end(), std::make_move_iterator(clustersA.begin()),
                           std::make_move_iterator(clustersA.end()));
  }
  if (rootB != nullptr) {
    data->particles.insert(data->particles.end(), std::make_move_iterator(clustersB.begin()),
                           std::make_move_iterator(clustersB.end()));
  }

  return std::move(data);
}

std::vector<MSTClustering::Edge> CoordinateMSTClustering::get_edges(const cola::EventData&)
// Notice that we don't use EventData in this implementation since we have the needed iterators. The possibility is
// still there though
{
  std::vector<Edge> edges;
  edges.reserve(std::pow(std::distance(spectIterA, spectIterB), 2) + std::pow(std::distance(spectIterB, endIter), 2));

  // particle vector is sorted, process kSpectatorA nucleons (check for no kSpectatorA nucleons)
  if (spectIterA != endIter) {
    for (auto iter = spectIterA; iter != spectIterB; ++iter) {
      for (auto jter = iter + 1; jter != spectIterB; ++jter) {
        auto delta = iter->position - jter->position;
        edges.emplace_back(std::make_pair(&(*iter), &(*jter)),
                           std::sqrt(delta.x * delta.x + delta.y * delta.y + delta.z * delta.z),
                           cola::ParticleClass::kSpectatorA);
      }
    }
  }
  // repeat for kSpectatorB nucleons
  for (auto iter = spectIterB; iter != endIter; ++iter) {
    for (auto jter = iter + 1; jter != endIter; ++jter) {
      auto delta = iter->position - jter->position;
      edges.emplace_back(std::make_pair(&(*iter), &(*jter)),
                         std::sqrt(delta.x * delta.x + delta.y * delta.y + delta.z * delta.z),
                         cola::ParticleClass::kSpectatorB);
    }
  }
  return edges;
}

double get_cd(double Ex, uint32_t A) {
  const auto Afloat = static_cast<double>(A);
  if (Ex / Afloat < eps0) {
    return d0;
  }
  double ex = Ex / Afloat;
  double dep = std::exp(-std::pow((ex / b_opt), a_opt)) * c_opt + d_opt;
  return d0 * std::pow(dep, 1. / 3.);
}

cola::LorentzVector ToColaLorentzVector(const G4LorentzVector& lv) {
  return {
      lv.e(),
      lv.x(),
      lv.y(),
      lv.z(),
  };
}

cola::EventParticles CoordinateMSTClustering::_process_side(const cola::EventData& data, cola::ParticleClass side) {
  cola::EventParticles clusters;

  auto& root = side == cola::ParticleClass::kSpectatorA ? rootA : rootB;
  auto bIter = side == cola::ParticleClass::kSpectatorA ? spectIterA : spectIterB;
  auto eIter = side == cola::ParticleClass::kSpectatorA ? spectIterB : endIter;

  uint32_t sourceA =
      cola::PdgToAZ(side == cola::ParticleClass::kSpectatorA ? data.ini_state.pdg_code_a : data.ini_state.pdg_code_b)
          .first;
  // Boost to rest frame for each set of kSpectators
  cola::LorentzVector pNucleus = {0.0, 0.0, 0.0, 0.0};
  for (auto particleIt = bIter; particleIt != eIter; ++particleIt) {
    pNucleus += particleIt->momentum;
  }
  for (auto particleIt = bIter; particleIt != eIter; ++particleIt) {
    particleIt->momentum.Boost(-pNucleus);
  }
  const auto count = std::distance(bIter, eIter);

  // get excitation energy
  double exEn = ExcitationEnergy(_stat_exen_type, sourceA).GetEnergy(count);

  // at this point the construct_tree() method has already built up MST trees for both sides

  double cd = get_cd(exEn, count);
  auto unprocessed = std::queue<Node*>();
  unprocessed.push(root);

  while (!unprocessed.empty()) {
    auto* topView = unprocessed.front();
    if (topView->height <= cd) {
      cola::AZ clusterAZ = {0, 0};
      cola::LorentzVector position{0, 0, 0, 0}, momentum{0, 0, 0, 0};

      for (const auto* nucleon : topView->vertices) {
        cola::AZ componentAZ = nucleon->GetAZ();
        clusterAZ.first += componentAZ.first;
        clusterAZ.second += componentAZ.second;
        position += nucleon->position;
        momentum += nucleon->momentum;
      }
      cola::Particle cluster;
      cluster.position = position / topView->vertices.size();
      cluster.momentum = momentum;
      cluster.pdg_code = cola::AZToPdg(clusterAZ);
      cluster.p_class = side;
      clusters.push_back(cluster);
    } else if (topView->children.has_value()) {
      unprocessed.push(topView->children.value().first);
      unprocessed.push(topView->children.value().second);
    }
    unprocessed.pop();
  }

  // now that we have defined clusters, we need to recalculate mass (and add momentum)

  std::vector<double> masses;
  masses.reserve(clusters.size());
  double totalMass = .0;
  double totalMassEx = 0.;

  for (const auto& cluster : clusters) {
    double mass = G4NucleiProperties::GetNuclearMass(static_cast<G4int>(cluster.GetAZ().first),
                                                     static_cast<G4int>(cluster.GetAZ().second));
    totalMass += mass;
    if (cluster.pdg_code != 2212 && cluster.pdg_code != 2112) mass += exEn * cluster.GetAZ().first / sourceA;
    masses.push_back(mass);
    totalMassEx += mass;
  }

  if (_extra_momentum && clusters.size() > 1) {
    double preFragmentMass = totalMass + exEn;
    if (preFragmentMass < totalMassEx + 1e-5 * MeV) preFragmentMass += 1e-5 * MeV;

    const auto generatedMomentum = phaseSpaceDecay.CalculateDecay(G4LorentzVector(preFragmentMass, {}), masses);

    for (size_t i = 0; i < clusters.size(); ++i) {
      clusters[i].momentum = ToColaLorentzVector(generatedMomentum.at(i));
    }
  } else {
    for (size_t i = 0; i < clusters.size(); ++i) {
      clusters[i].momentum.e = std::sqrt(std::pow(masses[i], 2) + clusters[i].momentum.SpatialPart().Mag2());
    }
  }

  // repulsion calculation
  if (_consider_rep) {
    RepulsionStage::CalculateRepulsion(std::move(clusters));
  }

  // Boost back from rest frame
  for (auto& cluster : clusters) {
    cluster.momentum.Boost(pNucleus);
  }

  return clusters;
}
