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

#include <G4Types.hh>
#ifndef ExcitationEnergy_h
#define ExcitationEnergy_h 1

#include <cmath>
#include <fstream>

class ExcitationEnergy {
 public:
  ExcitationEnergy(G4int ex_en_label_in, G4int init_a_in);

  ~ExcitationEnergy();

  G4double GetEnergy(G4int a);

  G4double GetEnergyEricson(G4int a) const;

  G4double GetEnergyGaimardSchmidt(G4int a) const;

  G4double GetEnergyALADIN(G4int a) const;

  G4double GetEnergyCorrectedALADIN(G4int a) const;

  G4double GetEnergyParabolicApproximation(G4int a) const;

  G4double GetEnergyDampEricson(G4int a) const;

  G4double GetEnergyHybridFit(G4int a) const;

  void SetInitNuclMass(G4int init_a_in);

  void SetParametersEricson(G4double g0_in);

  void SetParametersGaimardSchmidt(G4double g0_in, G4double g1_in);

  void SetParametersALADIN(G4double e0_in, G4double sigma0_in, G4double b0_in);

  void SetParametersCorrectedALADIN(G4double e0_in, G4double c0_in, G4double sigma0_in, G4double b0_in, G4double b1_in);

  void SetParametersCorrectedALADINFromFile();

  void SetParametersParabolicApproximation(G4double pe_in, G4double pm_in, G4double sigma_p_in, G4double b_p0_in,
                                           G4double b_p1_in);

  void SetParametersHybridFit(G4double a0_in, G4double a1_in, G4double a2_in, G4double a3_in, G4double a4_in,
                              G4double a5_in, G4double a6_in, G4double sigma1_in, G4double sigma2_in,
                              G4double sigma3_in);

 private:
  G4double g0_;
  G4double g1_;
  G4double e0_;
  G4double sigma0_;
  G4double b0_;
  G4double c0_;
  G4double b1_;
  G4double pe_;
  G4double pm_;
  G4double sigma_p_;
  G4double b_p0_;
  G4double b_p1_;
  G4double alpha_switch_;

  G4int init_a_;
  G4int ex_en_label_{3};

  G4double low_ex_en_;
  G4double up_ex_en_;
  G4double ebound_;

  G4double a0_;
  G4double a1_;
  G4double a2_;
  G4double a3_;
  G4double a4_;
  G4double a5_;
  G4double a6_;
  G4double sigma1_;
  G4double sigma2_;
  G4double sigma3_;

  std::ifstream param_file_;
};
#endif
