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

#include "ExcitationEnergy.hh"

#include "G4SystemOfUnits.hh"

#include <CLHEP/Random/RandGauss.h>
#include <CLHEP/Random/RandGeneral.h>
#include <G4Exception.hh>
#include <G4ExceptionSeverity.hh>
#include <G4Types.hh>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

// default values for excitation energy models parameters
namespace exc_en_default_vals {
  constexpr double kE0 = 11.5 * MeV;
  constexpr double kSigma0 = 0.005;
  constexpr double kB0 = 2;
  constexpr double kSigmaE0 = 1 * MeV;
  constexpr double kC0 = 0.1;
  constexpr double kPe = 24 * MeV;
  constexpr double kPm = 0.2;
}  // namespace exc_en_default_vals

namespace {
  G4double GaimardSchmidt(G4double energy, G4double evaporation_energy, G4int a_final, G4int a_initial) {
    G4double g0 = 16.;
    G4double g1 = 0.7;
    G4int removed_nucleons = a_initial - a_final;
    G4double res = 0.;

    if (energy > 0.5 * evaporation_energy * removed_nucleons) {  // Energy restriction for easier caclulations
      return res;
    }

    // First term calculation at H-S formula
    G4double term = g0;

    for (G4int k = 2; k < removed_nucleons - 1; k++) {
      term *= (g0 * energy) / static_cast<G4double>(k * (k - 1));
    }

    res = term;

    // Rest sum calculating with recursive form of H-S formula
    for (G4int i = 1; i < removed_nucleons; i++) {
      term *= (g1 * energy / g0) * ((-1) * static_cast<G4double>(removed_nucleons - i + 1) /
                                    (static_cast<G4double>(i) * static_cast<G4double>(removed_nucleons + i - 1)));
      // G4cout<<testRes<<G4endl; //NaN is result of inf-inf becouse term becomes inf at some step.
      res += term;
    }

    if (res != res) {
      res = 0;
    }  //- костыль!!!
    res = std::max<G4double>(res, 0);

    return res;
  }

  G4double Ericson(G4double energy, G4double evaporation_energy, G4int a_final, G4int a_initial) {
    G4double g0 =
        16;  // 16 was in Shidenberger - not influence calculations at all cause its freeze out in distr normalization.
    G4int removed_nucleons = a_initial - a_final;

    if (energy > evaporation_energy * (a_initial - a_final)) {
      G4double s_var = 0;
      return s_var;
    }

    G4double s_var = g0;

    for (G4int k = 2; k <= removed_nucleons; k++) {
      s_var *= (g0 * energy) / static_cast<G4double>(k * (k - 1));
    }

    return s_var;
  }
}  // namespace

ExcitationEnergy::ExcitationEnergy(G4int ex_en_label_in, G4int init_a_in) {
  ex_en_label_ = ex_en_label_in;
  init_a_ = init_a_in;
  low_ex_en_ = 0 * init_a_;
  up_ex_en_ = 100 * init_a_;
  ebound_ = 40;

#ifdef DATA_INSTALL
  std::string filepath(DATA_INSTALL);
#else
  std::string filepath("./ExcitationEnergy/");
#endif

  filepath += "CorrectedALADINParameters.dat";
  param_file_.open(filepath.c_str());

  SetParametersCorrectedALADINFromFile();
  SetParametersALADIN(exc_en_default_vals::kE0, exc_en_default_vals::kSigma0, exc_en_default_vals::kB0);
  SetParametersParabolicApproximation(exc_en_default_vals::kPe, exc_en_default_vals::kPm, exc_en_default_vals::kSigma0,
                                      exc_en_default_vals::kC0, 0.01);
  SetParametersHybridFit(11.46648905 * MeV, -1.84830078 * MeV, -58.53674677 * MeV, 284.66431513 * MeV,
                         -637.51406293 * MeV, 652.80324427 * MeV, -251.28205381 * MeV, 0.4 * MeV, 0.5, 0.2);
}

ExcitationEnergy::~ExcitationEnergy() = default;

void ExcitationEnergy::SetInitNuclMass(G4int init_a_in) { init_a_ = init_a_in; }

void ExcitationEnergy::SetParametersALADIN(G4double e0_in, G4double sigma0_in, G4double b0_in) {
  e0_ = e0_in;
  sigma0_ = sigma0_in;
  b0_ = b0_in;

  alpha_switch_ = (e0_ / ebound_) * (e0_ / ebound_);
  for (int k = 0; k < 100; k++) {
    alpha_switch_ = pow((e0_ / ebound_) * (alpha_switch_ * (1 - 1 / (static_cast<G4double>(init_a_))) +
                                           1 / static_cast<G4double>(init_a_) - alpha_switch_ * alpha_switch_),
                        0.66666666666666);
  }
  // std::cout<<"alphaSwitch is set equal to "<<alphaSwitch<<"\n";
}

void ExcitationEnergy::SetParametersEricson(G4double g0_in) { g0_ = g0_in; }

void ExcitationEnergy::SetParametersGaimardSchmidt(G4double g0_in, G4double g1_in) {
  g0_ = g0_in;
  g1_ = g1_in;
}

void ExcitationEnergy::SetParametersCorrectedALADIN(G4double e0_in, G4double c0_in, G4double sigma0_in, G4double b0_in,
                                                    G4double b1_in) {
  e0_ = e0_in;
  c0_ = c0_in;
  sigma0_ = sigma0_in;
  b0_ = b0_in;
  b1_ = b1_in;
}

void ExcitationEnergy::SetParametersCorrectedALADINFromFile() {
  std::vector<G4double> param_vect;
  G4double param = 0;
  G4int iter = 0;
  while (true) {
    param_file_ >> param;
    param_vect.push_back(param);
    if (!param_file_.good()) {
      break;
    }
    iter++;
  }
  SetParametersCorrectedALADIN(param_vect.at(0), param_vect.at(1), param_vect.at(2), param_vect.at(3),
                               param_vect.at(4));
}

void ExcitationEnergy::SetParametersParabolicApproximation(G4double pe_in, G4double pm_in, G4double sigma_p_in,
                                                           G4double b_p0_in, G4double b_p1_in) {
  pe_ = pe_in;
  pm_ = pm_in;
  sigma_p_ = sigma_p_in;
  b_p0_ = b_p0_in;
  b_p1_ = b_p1_in;
}

void ExcitationEnergy::SetParametersHybridFit(G4double a0_in, G4double a1_in, G4double a2_in, G4double a3_in,
                                              G4double a4_in, G4double a5_in, G4double a6_in, G4double sigma1_in,
                                              G4double sigma2_in, G4double sigma3_in) {
  a0_ = a0_in;
  a1_ = a1_in;
  a2_ = a2_in;
  a3_ = a3_in;
  a4_ = a4_in;
  a5_ = a5_in;
  a6_ = a6_in;
  sigma1_ = sigma1_in;
  sigma2_ = sigma2_in;
  sigma3_ = sigma3_in;
}

G4double ExcitationEnergy::GetEnergyALADIN(G4int a) const {
  CLHEP::RandGauss rand_gauss(nullptr, 1);
  G4double energy = 0;
  G4double alpha = static_cast<G4double>(a) / static_cast<G4double>(init_a_);
  G4double sigma_a = CLHEP::RandGauss::shoot() * sigma0_ * (1 + b0_ * (1 - alpha));
  G4double sigma_e = CLHEP::RandGauss::shoot() * sigma0_ * (1 + b0_ * (1 - alpha));
  energy = e0_ * a * pow(1 - alpha - sigma_a, 0.5) + 0 * a * sigma_e;

  while (energy < 0 || energy != energy) {
    sigma_e = CLHEP::RandGauss::shoot() * sigma0_ * (1 + b0_ * (1 - alpha));
    sigma_a = CLHEP::RandGauss::shoot() * sigma0_ * (1 + b0_ * (1 - alpha));
    energy = e0_ * a * pow(1 - alpha - sigma_a, 0.5) + 0 * a * sigma_e;
  }
  return energy;
}

G4double ExcitationEnergy::GetEnergyCorrectedALADIN(G4int a) const {
  CLHEP::RandGauss rand_gauss(nullptr, 1);
  G4double energy = 0;
  G4double alpha = static_cast<G4double>(a) / static_cast<G4double>(init_a_);
  G4double sigma_e = CLHEP::RandGauss::shoot() * sigma0_ * (1 + b0_ * (1 - alpha) + b1_ * (1 - alpha) * (1 - alpha));
  energy = a * std::sqrt(std::sqrt(e0_ * e0_ + c0_ * (1 - alpha)) - e0_) + a * sigma_e;
  while (energy < 0 || energy != energy) {
    sigma_e = CLHEP::RandGauss::shoot() * sigma0_ * (1 + b0_ * (1 - alpha) + b1_ * (1 - alpha) * (1 - alpha));
    energy = a * std::sqrt(std::sqrt(e0_ * e0_ - c0_ * (1 - alpha)) - e0_) + a * sigma_e;
  }
  return energy;
}

G4double ExcitationEnergy::GetEnergyEricson(G4int a) const {
  G4double excitation_energy_distribution[10000];  // NOLINT

  for (G4int n = 0; n < 10000; n++) {
    G4double sum = Ericson(static_cast<G4double>(n) * ((up_ex_en_ - low_ex_en_) / static_cast<G4double>(10000)),
                           ebound_, a, init_a_) *
                   (up_ex_en_ - low_ex_en_) / static_cast<G4double>(10000);
    excitation_energy_distribution[n] = sum;
  }
  CLHEP::RandGeneral rand_general(excitation_energy_distribution, 10000);

  G4double energy = rand_general.shoot() * (up_ex_en_ - low_ex_en_) + low_ex_en_;

  return energy;
}

G4double ExcitationEnergy::GetEnergyGaimardSchmidt(G4int a) const {
  G4double excitation_energy_distribution[10000];  // NOLINT

  for (G4int n = 0; n < 10000; n++) {
    G4double sum = GaimardSchmidt(static_cast<G4double>(n) * ((up_ex_en_ - low_ex_en_) / static_cast<G4double>(10000)),
                                  ebound_, a, init_a_) *
                   (up_ex_en_ - low_ex_en_) / static_cast<G4double>(10000);
    excitation_energy_distribution[n] = sum;
  }
  CLHEP::RandGeneral rand_general(excitation_energy_distribution, 10000);

  G4double energy = rand_general.shoot() * (up_ex_en_ - low_ex_en_) + low_ex_en_;

  return energy;
}
G4double ExcitationEnergy::GetEnergyParabolicApproximation(G4int a) const {
  CLHEP::RandGauss rand_gauss(nullptr, 1);
  G4double energy;
  G4double alpha = static_cast<G4double>(a) / static_cast<G4double>(init_a_);
  G4double sigma_e =
      CLHEP::RandGauss::shoot() * sigma_p_ * (1 + b_p0_ * (1 - alpha) + b_p1_ * (1 - alpha) * (1 - alpha));
  energy = pe_ * static_cast<G4double>(a) * (1 - alpha) * (alpha + pm_) + sigma_e;
  return energy;
}

G4double ExcitationEnergy::GetEnergyDampEricson(G4int a) const {
  G4double energy = 0;
  if (init_a_ - a < alpha_switch_ * init_a_) {
    energy = GetEnergyEricson(a);
  } else {
    energy = GetEnergyALADIN(a);
  }
  return energy;
}

G4double ExcitationEnergy::GetEnergyHybridFit(G4int a) const {
  CLHEP::RandGauss rand_gauss(nullptr, 1);
  G4double alpha = static_cast<G4double>(a) / static_cast<G4double>(init_a_);
  G4double energy = -1;
  G4double sigma_e;
  while (energy < 0) {
    G4double rand = CLHEP::RandGauss::shoot();
    if (alpha > alpha_switch_) {
      sigma_e = rand *
                ((sigma3_ - 1) / (1 - alpha_switch_) * alpha + (1 - sigma3_ * alpha_switch_) / (1 - alpha_switch_)) *
                sigma1_ * (1 + sigma2_ * (1 - alpha));
    } else {
      sigma_e = rand * sigma1_ * (1 + sigma2_ * (1 - alpha));
    }
    energy = a * (a0_ + a1_ * alpha + a2_ * pow(alpha, 2) + a3_ * pow(alpha, 3) + a4_ * pow(alpha, 4) +
                  a5_ * pow(alpha, 5) + a6_ * pow(alpha, 6) + sigma_e);
  }
  return energy;
}

G4double ExcitationEnergy::GetEnergy(G4int a) const {
  G4double energy;
  switch (ex_en_label_) {
    case 1: {
      return GetEnergyEricson(a);
      break;
    }
    case 2: {
      return GetEnergyGaimardSchmidt(a);
      break;
    }
    case 3: {
      return GetEnergyALADIN(a);
      break;
    }
    case 4: {
      return GetEnergyDampEricson(a);
      break;
    }
    case 5: {
      return GetEnergyParabolicApproximation(a);
      break;
    }
    case 6: {
      return GetEnergyCorrectedALADIN(a);
      break;
    }
    case 7: {
      return GetEnergyHybridFit(a);
      break;
    }
    default: {
      G4Exception("Statistics label", "GRATE-1", FatalException, "Statistics label is invalid");
      return .0;
      break;
    }
  }
}