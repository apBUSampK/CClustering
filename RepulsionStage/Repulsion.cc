#include "Repulsion.hh"

namespace RepulsionStage {

cola::EventParticles CalculateRepulsion(cola::EventParticles&& frags) {
  if (frags.empty()) {
    return frags;
  }

  double time = 0.0;
  double delta_time = totalTime / static_cast<double>(iterations);

  while (time < totalTime) {
    BHTree bhtree(frags);

    double temp_timedelta = std::min(delta_time, bhtree.GetAdaptiveTimeDelta());

    auto r_delta = bhtree.Iterate(temp_timedelta);

    /*for (auto& frag : frags) {
      nucleons[i].position.x += r_delta[maps[i]].x();
      nucleons[i].position.y += r_delta[maps[i]].y();
      nucleons[i].position.z += r_delta[maps[i]].z();
    }*/
    time += temp_timedelta;
  }
  return frags;
}

BHTree::BHTree(cola::EventParticles& frags) : frags_(frags) {
  BuildBHTree(frags);
  GetForces(rootnode_.get());
}

std::unique_ptr<BHNode> BHTree::InitializeRoot(const cola::EventParticles& frags) {
  cola::Vector3 minVector;
  cola::Vector3 maxVector;
  minVector.x = std::min_element(frags.begin(), frags.end(), [](const cola::Particle& l, const cola::Particle& r) {return l.position.x < r.position.x;})->position.x;
  minVector.y = std::min_element(frags.begin(), frags.end(), [](const cola::Particle& l, const cola::Particle& r) {return l.position.y < r.position.y;})->position.y;
  minVector.z = std::min_element(frags.begin(), frags.end(), [](const cola::Particle& l, const cola::Particle& r) {return l.position.z < r.position.z;})->position.z;
  maxVector.x = std::max_element(frags.begin(), frags.end(), [](const cola::Particle& l, const cola::Particle& r) {return l.position.x < r.position.x;})->position.x;
  maxVector.y = std::max_element(frags.begin(), frags.end(), [](const cola::Particle& l, const cola::Particle& r) {return l.position.y < r.position.y;})->position.y;
  maxVector.z = std::max_element(frags.begin(), frags.end(), [](const cola::Particle& l, const cola::Particle& r) {return l.position.z < r.position.z;})->position.z;

  auto range = maxVector - minVector;
  double maxRange = std::max(range.x, range.y, range.z);
  return std::make_unique<BHNode>(maxRange, (minVector + maxVector) / 2);
}

void BHTree::BuildBHTree(const cola::EventParticles& frags) {
  rootnode_ = InitializeRoot(frags);

  for (size_t i = 0; i < frags.size(); i++) {
    InsertNucleon(rootnode_, frags.at(i).position.spatialPart(), i);
  }
}

void BHTree::InsertNucleon(const std::unique_ptr<BHNode>& node, const cola::Vector3<double>& cords, int pIndex) {
  if (node->Z == 0) {
    node->Z = 1;
    node->cr = cords;
    node->index = pIndex;
    return;
  }
  if (node->Z == 1) {
    node->Divide();

    int index = 0;
    if (node->cr.x > node->ctr.x) index |= 1;
    if (node->cr.y > node->ctr.y) index |= 2;
    if (node->cr.z > node->ctr.z) index |= 4;
    InsertNucleon(node->children[index], node->cr, node->index);

    index = 0;
    if (cords.x > node->ctr.x) index |= 1;
    if (cords.y > node->ctr.y) index |= 2;
    if (cords.z > node->ctr.z) index |= 4;
    InsertNucleon(node->children[index], cords, pIndex);

    node->cr = (node->cr + cords) / 2;
    node->Z += 1;
    node->index = -1;
    return;
  }
  node->cr = (node->cr * node->Z + cords) / (node->Z + 1);
  node->Z += 1;

  int index = 0;
  if (cords.x > node->ctr.x) index |= 1;
  if (cords.y > node->ctr.y) index |= 2;
  if (cords.z > node->ctr.z) index |= 4;
  InsertNucleon(node->children[index], cords, pIndex);
}

double BHTree::GetAdaptiveTimeDelta() const {
  double min_time = max_adaptive_delta;
  for (size_t i = 0; i < frags_.size(); i++) {
    auto mnt = frags_.at(i).momentum;
    double mval = std::sqrt(mnt.x * mnt.x + mnt.y * mnt.y + mnt.z * mnt.z);
    if (!frags_.at(i).getAZ().second || !fs_[i].mag() || !mval) {
      continue;
    }
    min_time = std::min(min_time, 0.05 *  mval / fs_[i].mag());
  }
  return min_time;
}

std::vector<cola::Vector3<double>> BHTree::Iterate(double time_delta) {
  std::vector<cola::Vector3<double>> r_delta(frags_.size());

  for (size_t i = 0; i < frags_.size(); ++i) {
    if (frags_.at(i).getAZ().second == 0) {
      continue;
    }
    cola::LorentzVector p = frags_.at(i).momentum;
    cola::Vector3<double> half_dp = fs_[i] * time_delta * 0.5;

    cola::Vector3<double> mid_v = (p.spatialPart() + half_dp) / std::sqrt((p.spatialPart() + half_dp).mag2() + p.mag2());
    r_delta[i] = cola::Vector3<double>{time_delta * mid_v.x, time_delta * mid_v.y, time_delta * mid_v.z};

    cola::Vector3<double> pvec = p.spatialPart() + 2 * half_dp;
    p.x = pvec.x;
    p.y = pvec.y;
    p.z = pvec.z;
    p.e = std::sqrt((p.spatialPart() + 2 * half_dp).mag2() + p.mag2());

    frags_.at(i).momentum = p;
  }
  return r_delta;
}

void BHTree::GetForces(const BHNode* node) {
  if (node->Z == 0) {
    return;
  }
  if (node->Z == 1) {
    fs_[maps_node->index)] += Force(rootnode_.get(), node);
    return;
  }
  for (const auto& child : node->children) {
    GetForces(child.get());
  }
}

cola::Vector3<double> BHTree::Force(const BHNode* rootnode, const BHNode* node) const {
  if ((rootnode->Z == 0) || ((rootnode->index != -1) && maps_->at(rootnode->index) == maps_->at(node->index))) {
    return {0.0, 0.0, 0.0};
  }
  if (rootnode->Z == 1) {
    return DuoForce(node->cr - rootnode->cr, rootnode->Z);
  }
  if ((rootnode->size / (node->cr - rootnode->cr).mag()) < theta) {
    return DuoForce(node->cr - rootnode->cr, rootnode->Z);
  }
  cola::Vector3<double> totalForce = {0.0, 0.0, 0.0};
  for (const auto& child : rootnode->children) {
    totalForce += Force(child.get(), node);
  }
  return totalForce;
}

cola::Vector3<double> BHTree::DuoForce(const cola::Vector3<double> vec, const double& from_Z) const {
  cola::Vector3<double> fos = vec * CLHEP::elm_coupling * from_Z / std::pow(vec.mag(), 3);
  return fos;
}

void BHNode::Divide() {
  children.resize(8);
  for (size_t i = 0; i < 8; i++) {
    cola::Vector3<double> offset(
      (i & 1 ? 1 : -1) * size / 4.0,
      (i & 2 ? 1 : -1) * size / 4.0,
      (i & 4 ? 1 : -1) * size / 4.0
    );
    children[i] = std::make_unique<BHNode>(size / 2.0, ctr + offset);
  }
}
}