#include "CoordinateMSTClustering.hh"

#include "Repulsion.hh"
#include "ExcitationEnergy.hh"

#include <limits>
#include <queue>

#include "G4SystemOfUnits.hh"

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

std::unique_ptr<cola::EventData> CoordinateMSTClustering::get_clusters(std::unique_ptr<cola::EventData> &&data)
{
    // get clusters
    auto clustersA = _process_side(*data, cola::ParticleClass::spectatorA);
    auto clustersB = _process_side(*data, cola::ParticleClass::spectatorB);
    // erase spectator nucleons
    data->particles.erase(spectIterA, endIter);
    // append clusters
    data->particles.insert(data->particles.end(), clustersA.begin(), clustersA.end());
    data->particles.insert(data->particles.end(), clustersB.begin(), clustersB.end());

    return std::move(data);
}

std::vector<MSTClustering::Edge> CoordinateMSTClustering::get_edges(const cola::EventData&)
// Notice that we don't use EventData in this implementation since we have the needed iterators. The possibility is still there though
{
    std::vector<Edge> edges;
    
    // particle vector is sorted, process spectatorA nucleons
    for (auto iter = spectIterA; iter != spectIterB; iter++)
    {
        for (auto jter = iter + 1; jter != spectIterB; jter++)
        {
            auto delta = iter->position - jter->position;
            edges.emplace_back(std::make_pair(&(*iter), &(*jter)), std::sqrt(delta.x * delta.x + delta.y * delta.y + delta.z * delta.z), cola::ParticleClass::spectatorA);
        }
    }
    // repeat for spectatorB nucleons
    for (auto iter = spectIterB; iter != endIter; iter++)
    {
        for (auto jter = iter + 1; jter != endIter; jter++)
        {
            auto delta = iter->position - jter->position;
            edges.emplace_back(std::make_pair(&(*iter), &(*jter)), std::sqrt(delta.x * delta.x + delta.y * delta.y + delta.z * delta.z), cola::ParticleClass::spectatorB);
        }
    }
    return edges;
}

double get_cd(double Ex, uint A)
{
    if ((Ex / double(A)) < eps0)
    {
        return d0;
    }
    double ex = Ex / double(A);
    double dep = std::exp(-std::pow((ex / b_opt), a_opt)) * c_opt + d_opt;
    return d0 * std::pow(dep, 1. / 3.);
}

cola::LorentzVector ToColaLorentzVector(G4LorentzVector& lv)
{
    cola::LorentzVector vec;
    vec.e = lv.e();
    vec.x = lv.x();
    vec.y = lv.y();
    vec.z = lv.z();
    return vec;
}

cola::EventParticles CoordinateMSTClustering::_process_side(const cola::EventData& data, cola::ParticleClass side) 
{
    cola::EventParticles clusters;
    uint count = 0;

    auto& root = side == cola::ParticleClass::spectatorA ? rootA : rootB;
    uint sourceA = cola::pdgToAZ(side == cola::ParticleClass::spectatorA ? data.iniState.pdgCodeA : data.iniState.pdgCodeB).first;
    if(root.get() == nullptr)
    {
        return clusters;
    }
    // boost to rest frame for each set of spectators
    cola::LorentzVector pNucleus = {0.0, 0.0, 0.0, 0.0};
    for (auto particle = spectIterA; particle != spectIterB; particle++)
        pNucleus += particle->momentum;
    for (auto particle = spectIterA; particle != endIter; particle++)
    {
        particle->momentum.boost(-pNucleus);
        count++;
    }

    // get excitation energy
    double exEn = ExcitationEnergy(_stat_exen_type, sourceA).GetEnergy(count); 

    // at this point the construct_tree() method has already built up MST trees for both sides
    
    double cd = get_cd(exEn, count);
    auto unprocessed = std::queue<std::shared_ptr<Node>>();
    unprocessed.push(root);

    while (!unprocessed.empty()) {
        Node* topView = unprocessed.front().get();
        if (topView->height <= cd)
        {
            cola::AZ clusterAZ = {0, 0};
            cola::LorentzVector position, momentum;

            for (auto nucleon : topView->vertices)
            {
                cola::AZ componentAZ = nucleon->getAZ();
                clusterAZ.first += componentAZ.first;
                clusterAZ.second += componentAZ.second;
                position += nucleon->position;
                momentum += nucleon->momentum;
            }

            cola::Particle cluster;
            cluster.position = position / topView->vertices.size();
            cluster.momentum = momentum;
            cluster.pdgCode = cola::AZToPdg(clusterAZ);
            cluster.pClass = side;
            clusters.push_back(cluster);
        } else if (topView->children.has_value()) 
        {
            unprocessed.push(topView->children.value().first);
            unprocessed.push(topView->children.value().second);
        }
        unprocessed.pop();
    }
    
    // now that we have defined clusters, we need to calculate additional momentum

    if (_extra_momentum) {
        std::vector<double> mass;
        double totalMass = .0;

        for (auto& cluster: clusters)
        {
            mass.push_back(G4NucleiProperties::GetNuclearMass(static_cast<G4int>(cluster.getAZ().first), static_cast<G4int>(cluster.getAZ().second)) + exEn * cluster.getAZ().first / sourceA);
            totalMass += mass.back();
        }

        totalMass += 1e-5*MeV; // fix for segfault
        auto extra_momentum = phaseSpaceDecay.Decay(totalMass, mass);

        for (size_t i = 0; i < clusters.size(); i++)
        {
            clusters.at(i).momentum += ToColaLorentzVector(*extra_momentum->at(i));
        }

        // clear generated vector
        for (auto pointer: *extra_momentum)
        {
            delete pointer;
        }
        delete extra_momentum;
    }

    if (_consider_rep) 
    {
        // do something
    }

    // make mass consistent with G4 tables and boost back from rest frame
    for (auto& cluster: clusters)
    {
        cluster.momentum.e = std::sqrt(
            std::pow(G4NucleiProperties::GetNuclearMass(static_cast<G4int>(cluster.getAZ().first), static_cast<G4int>(cluster.getAZ().second)), 2) +
            cluster.momentum.x*cluster.momentum.x + cluster.momentum.y*cluster.momentum.y + cluster.momentum.z*cluster.momentum.z);
        cluster.momentum.boost(pNucleus);
    }

    return clusters;
}