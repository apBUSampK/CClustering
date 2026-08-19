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
#include "MSTClustering.hh"
#include "Repulsion.hh"

#include <EventData.hh>
#include <G4LorentzVector.hh>
#include <G4NucleiProperties.hh>  //NOLINT
#include <G4SystemOfUnits.hh>
#include <G4Types.hh>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <memory>
#include <queue>
#include <utility>
#include <vector>

// constants for implementation
constexpr double kEps0 = 2.17 * MeV;
constexpr double kAlphaPow = -1.02;
constexpr double kD0 = 2.7;
constexpr double kAColRel = 5.315;
constexpr double kASurfRel = 0.017;
constexpr double kAVol = 0.054;

constexpr double kSpecAa = 0;
constexpr double kSpecAb = 0;

constexpr double kAOpt = 2.243;
constexpr double kBOpt = 3.183 * MeV;
constexpr double kCOpt = 0.99;
constexpr double kDOpt = 0.29041;

using namespace cola;

std::unique_ptr<cola::EventData> CoordinateMSTClustering::GetClusters(std::unique_ptr<cola::EventData>&& data) {
  cola::EventParticles clusters_a;
  cola::EventParticles clusters_b;

  // get clusters
  if (root_a_ != nullptr) {
    clusters_a = ProcessSide(*data, cola::ParticleClass::kSpectatorA);
  }
  if (root_b_ != nullptr) {
    clusters_b = ProcessSide(*data, cola::ParticleClass::kSpectatorB);
  }

  // erase kSpectator nucleons
  data->particles.erase(spect_iter_a_ != end_iter_ ? spect_iter_a_ : spect_iter_b_, end_iter_);

  // append clusters
  if (root_a_ != nullptr) {
    data->particles.insert(data->particles.end(), std::make_move_iterator(clusters_a.begin()),
                           std::make_move_iterator(clusters_a.end()));
  }
  if (root_b_ != nullptr) {
    data->particles.insert(data->particles.end(), std::make_move_iterator(clusters_b.begin()),
                           std::make_move_iterator(clusters_b.end()));
  }

  return std::move(data);
}

std::vector<MSTClustering::Edge> CoordinateMSTClustering::GetEdges(const cola::EventData& /*unused*/)
// Notice that we don't use EventData in this implementation since we have the needed iterators. The possibility is
// still there though
{
  std::vector<Edge> edges;
  edges.reserve(std::pow(std::distance(spect_iter_a_, spect_iter_b_), 2) +
                std::pow(std::distance(spect_iter_b_, end_iter_), 2));

  // particle vector is sorted, process kSpectatorA nucleons (check for no kSpectatorA nucleons)
  if (spect_iter_a_ != end_iter_) {
    for (auto iter = spect_iter_a_; iter != spect_iter_b_; ++iter) {
      for (auto jter = iter + 1; jter != spect_iter_b_; ++jter) {
        auto delta = iter->position - jter->position;
        edges.emplace_back(std::make_pair(&(*iter), &(*jter)),
                           std::sqrt(delta.x * delta.x + delta.y * delta.y + delta.z * delta.z),
                           cola::ParticleClass::kSpectatorA);
      }
    }
  }
  // repeat for kSpectatorB nucleons
  for (auto iter = spect_iter_b_; iter != end_iter_; ++iter) {
    for (auto jter = iter + 1; jter != end_iter_; ++jter) {
      auto delta = iter->position - jter->position;
      edges.emplace_back(std::make_pair(&(*iter), &(*jter)),
                         std::sqrt(delta.x * delta.x + delta.y * delta.y + delta.z * delta.z),
                         cola::ParticleClass::kSpectatorB);
    }
  }
  return edges;
}

namespace {

  double GetCd(double ex, uint32_t a) {
    const auto afloat = static_cast<double>(a);
    if (ex / afloat < kEps0) {
      return kD0;
    }
    ex = ex / afloat;
    double dep = std::exp(-std::pow((ex / kBOpt), kAOpt)) * kCOpt + kDOpt;
    return kD0 * std::pow(dep, 1. / 3.);
  }

  cola::LorentzVector ToColaLorentzVector(const G4LorentzVector& lv) {
    return {
        lv.e(),
        lv.x(),
        lv.y(),
        lv.z(),
    };
  }
}  // namespace

cola::EventParticles CoordinateMSTClustering::ProcessSide(const cola::EventData& data, cola::ParticleClass side) {
  cola::EventParticles clusters;

  auto& root = side == cola::ParticleClass::kSpectatorA ? root_a_ : root_b_;
  auto b_iter = side == cola::ParticleClass::kSpectatorA ? spect_iter_a_ : spect_iter_b_;
  auto e_iter = side == cola::ParticleClass::kSpectatorA ? spect_iter_b_ : end_iter_;

  uint32_t source_a =
      cola::PdgToAZ(side == cola::ParticleClass::kSpectatorA ? data.ini_state.pdg_code_a : data.ini_state.pdg_code_b)
          .first;
  // Boost to rest frame for each set of kSpectators
  cola::LorentzVector p_nucleus = {0.0, 0.0, 0.0, 0.0};
  for (auto particle_it = b_iter; particle_it != e_iter; ++particle_it) {
    p_nucleus += particle_it->momentum;
  }
  for (auto particle_it = b_iter; particle_it != e_iter; ++particle_it) {
    particle_it->momentum.Boost(-p_nucleus);
  }
  const auto count = std::distance(b_iter, e_iter);

  // get excitation energy
  double ex_en = ExcitationEnergy(stat_exen_type_, source_a).GetEnergy(count);

  // at this point the construct_tree() method has already built up MST trees for both sides

  double cd = GetCd(ex_en, count);
  auto unprocessed = std::queue<Node*>();
  unprocessed.push(root);

  while (!unprocessed.empty()) {
    auto* top_view = unprocessed.front();
    if (top_view->height <= cd) {
      cola::AZ cluster_az = {0, 0};
      cola::LorentzVector position{0, 0, 0, 0};
      cola::LorentzVector momentum{0, 0, 0, 0};

      for (const auto* nucleon : top_view->vertices) {
        cola::AZ component_az = nucleon->GetAZ();
        cluster_az.first += component_az.first;
        cluster_az.second += component_az.second;
        position += nucleon->position;
        momentum += nucleon->momentum;
      }
      cola::Particle cluster;
      cluster.position = position / top_view->vertices.size();
      cluster.momentum = momentum;
      cluster.pdg_code = cola::AZToPdg(cluster_az);
      cluster.p_class = side;
      clusters.push_back(cluster);
    } else if (top_view->children.has_value()) {
      unprocessed.push(top_view->children.value().first);
      unprocessed.push(top_view->children.value().second);
    }
    unprocessed.pop();
  }

  // now that we have defined clusters, we need to recalculate mass (and add momentum)

  std::vector<double> masses;
  masses.reserve(clusters.size());
  double total_mass = .0;
  double total_mass_ex = 0.;

  for (const auto& cluster : clusters) {
    double mass = G4NucleiProperties::GetNuclearMass(static_cast<G4int>(cluster.GetAZ().first),
                                                     static_cast<G4int>(cluster.GetAZ().second));
    total_mass += mass;
    if (cluster.pdg_code != 2212 && cluster.pdg_code != 2112) {
      mass += ex_en * cluster.GetAZ().first / source_a;
    }
    masses.push_back(mass);
    total_mass_ex += mass;
  }

  if (extra_momentum_ && clusters.size() > 1) {
    double pre_fragment_mass = total_mass + ex_en;
    if (pre_fragment_mass < total_mass_ex + 1e-5 * MeV) {
      pre_fragment_mass += 1e-5 * MeV;
    }

    const auto generated_momentum = phase_space_decay_.CalculateDecay(G4LorentzVector(pre_fragment_mass, {}), masses);

    for (size_t i = 0; i < clusters.size(); ++i) {
      clusters[i].momentum = ToColaLorentzVector(generated_momentum.at(i));
    }
  } else {
    for (size_t i = 0; i < clusters.size(); ++i) {
      clusters[i].momentum.e = std::sqrt(std::pow(masses[i], 2) + clusters[i].momentum.SpatialPart().Mag2());
    }
  }

  // repulsion calculation
  if (consider_rep_) {
    clusters = repulsion_stage::CalculateRepulsion(std::move(clusters));
  }

  // Boost back from rest frame
  for (auto& cluster : clusters) {
    cluster.momentum.Boost(p_nucleus);
  }

  return clusters;
}
