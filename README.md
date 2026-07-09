# CPreEquilibrium
## Preequilibrium de-excitation module for COLA

Configurable parameters:
* `clustering` - clustering type, only GMST is available now.
* `stat_exen_type` - model used to calculate additional energy from missing nucleons in the pre-fragment. Available options: 1 to 7, defaults to 4. See `ExcitationEnergy` files for more info.
* `simulate_mometum` - whether to simulate additional fragments momentum based on excitation energy, 0 or 1, defaults to 1 
* `consider_coulomb` - whether to consider coulomb repulsion after clusterization, 0 or 1, defaults to 0. Warning: even with Barnes-Hut algorithm implementation, this is very compute heavy.