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

#ifndef CCLUSTERING_GMSTCLUSTERING_H
#define CCLUSTERING_GMSTCLUSTERING_H

#include "MSTClustering.hh"

#include <CLHEP/Units/SystemOfUnits.h>
#include <EventData.hh>
#include <G4FermiPhaseDecay.hh>

#include <cstdint>
#include <memory>
#include <vector>

namespace cola {

  class CoordinateMSTClustering : public MSTClustering {
   public:
    CoordinateMSTClustering() = delete;
    CoordinateMSTClustering(bool consider_rep, bool extra_momentum, int stat_exen_type)
        : consider_rep_(consider_rep), extra_momentum_(extra_momentum), stat_exen_type_(stat_exen_type) {};

   private:
    static constexpr double kNucleonAverMass = 0.93891875434 * CLHEP::GeV;
    bool consider_rep_;
    bool extra_momentum_;
    uint32_t stat_exen_type_;

    std::vector<Edge> GetEdges(const EventData& /*unused*/) final;
    std::unique_ptr<EventData> GetClusters(std::unique_ptr<EventData>&& /*data*/) final;
    EventParticles ProcessSide(const EventData& /*data*/, ParticleClass /*side*/);

    G4FermiPhaseDecay phase_space_decay_;
  };

}  // namespace cola
#endif  // CCLUSTERING_GMSTCLUSTERING_H