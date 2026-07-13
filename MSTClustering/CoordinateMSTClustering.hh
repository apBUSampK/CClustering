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

#include <G4ExcitationHandler.hh>
#include <G4FermiPhaseDecay.hh>
#include <G4NucleiProperties.hh>
#include <G4ReactionProductVector.hh>

#include <algorithm>
#include <memory>

namespace cola {

  class CoordinateMSTClustering : public MSTClustering {
   public:
    CoordinateMSTClustering() = delete;
    CoordinateMSTClustering(bool consider_rep, bool extra_momentum, int stat_exen_type)
        : _consider_rep(consider_rep), _extra_momentum(extra_momentum), _stat_exen_type(stat_exen_type) {};

   private:
    static constexpr double nucleonAverMass = 0.93891875434 * CLHEP::GeV;
    bool _consider_rep;
    bool _extra_momentum;
    uint32_t _stat_exen_type;

    std::vector<Edge> get_edges(const EventData&) final;
    std::unique_ptr<EventData> get_clusters(std::unique_ptr<EventData>&&) final;
    EventParticles _process_side(const EventData&, ParticleClass);

    G4FermiPhaseDecay phaseSpaceDecay;
  };

}  // namespace cola
#endif  // CCLUSTERING_GMSTCLUSTERING_H