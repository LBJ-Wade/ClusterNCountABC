#!/usr/bin/python3

import threading
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm
from numpy import pi
import os
import os.path

Ncm.cfg_init ()

z_max = 1.32
z_min = 0.3
area = 2500

cosmo = Nc.HICosmo.new_from_name (Nc.HICosmo, "NcHICosmoDEXcdm{'H0':<71.15>, 'Omegac':<0.218>, 'Omegab':<0.044>, 'ns':<0.97>, 'sigma8':<0.807>, 'w':<-1.01>, 'Omegac-fit':<1>, 'sigma8-fit':<1>, 'w-fit':<1>}")
cosmo.omega_x2omega_k ()
cosmo.param_set_by_name ("Omegak", 0.0);

(Found, w_i) = cosmo.param_index_from_name ("w")
cosmo.param_set_lower_bound (w_i, -3.0)
cosmo.param_set_upper_bound (w_i,  0.0)

(Found, Omegac_i) = cosmo.param_index_from_name ("Omegac")
cosmo.param_set_lower_bound (Omegac_i, 0.0)
cosmo.param_set_upper_bound (Omegac_i, 1.0)

(Found, sigma8_i) = cosmo.param_index_from_name ("sigma8")
cosmo.param_set_lower_bound (sigma8_i, 0.5)
cosmo.param_set_upper_bound (sigma8_i, 1.1)

dist = Nc.Distance.new (z_max * 2.0)
wp =  Nc.Window.new_from_name ("NcWindowTophat")
tf = Nc.TransferFunc.new_from_name ("NcTransferFuncEH")
vp = Nc.MatterVar.new (Nc.MatterVarStrategy.FFT, wp, tf)
gf = Nc.GrowthFunc.new ()
mulf = Nc.MultiplicityFunc.new_from_name ("NcMultiplicityFuncTinkerCrit{'Delta':<500.0>}")
mf = Nc.MassFunction.new (dist, vp, gf, mulf)

del dist
del wp
del tf
del vp
del mulf

cluster_m = Nc.ClusterMass.new_from_name ("NcClusterMassBenson{'M0':<3e14>, 'z0':<0.6>, 'signif-obs-min':<5.0>, 'Asz':<6.24>, 'Bsz':<1.33>, 'Csz':<0.83>, 'Dsz':<0.24>}")
cluster_z = Nc.ClusterRedshift.new_from_name ("NcClusterPhotozGaussGlobal{'pz-min':<%f>, 'pz-max':<%f>, 'z-bias':<0.0>, 'sigma0':<0.05>}" % (z_min, z_max))
cad = Nc.ClusterAbundance.new (mf, None, cluster_z, cluster_m)

ncount = Nc.DataClusterNCount.new (cad)

del cad

mset = Ncm.MSet ()
mset.set (cosmo)
mset.set (cluster_m)
mset.set (cluster_z)

rng = Ncm.RNG.pool_get ("example_ca_sampling");
rng.set_seed (123)

ncount.init_from_sampling (mset, cluster_z, cluster_m, area * (pi / 180.0)**2, rng)
ncount.true_data (False)
dset = Ncm.Dataset.new ()
dset.append_data (ncount)

del cosmo
del cluster_m
del cluster_z
del ncount

mset.prepare_fparam_map ()

prior = Ncm.MSetTransKernFlat.new ()
prior.set_mset (mset)
prior.set_prior_from_mset ()

abc = Nc.ABCClusterNCount.new (mset, prior, dset);
abc.props.epsilon = 1.2
abc.set_scale_cov (False)
abc.set_bin_uniform (200, 200)
#q = Ncm.Vector.new_array ([0.05 * (i+1) for i in range(19)])
#abc.set_bin_quantile (q)
#abc.props.epsilon_update_type = Nc.ABCClusterNCountEpsilonUpdate.UNIFORM
abc.set_epsilon_update (0.75)

#del q
del mset
del prior
del dset

abc.set_nthreads (10000);

abc.set_mtype (Ncm.FitRunMsgs.FULL) #FULL, SIMPLE
abc.set_rng (rng)
del rng

Ncm.func_eval_set_max_threads (4)
Ncm.func_eval_log_pool_stats ()

abc.start_run ()
abc.run (100)
abc.end_run ()

while True:
  abc.start_update ()
  abc.update ()
  abc.end_update ()
  if abc.get_accept_rate () < 0.002:
    i = 0
    filename = "nc_cluster_run_n.fits"
    while os.path.exists (filename):
      filename = "nc_cluster_run_n_%d.fits" % i
      i += 1
    abc.set_data_file (filename)
    break    

del abc

