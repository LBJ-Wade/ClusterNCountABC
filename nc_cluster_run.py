#!/usr/bin/python3

try:
  import gi
  gi.require_version('NumCosmo', '1.0')
  gi.require_version('NumCosmoMath', '1.0')
except:
  pass

import math
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm

#import threading
#from numpy import pi
import os
import os.path
import math

Ncm.cfg_init ()

z_max = 1.32
z_min = 0.3
area  = 2500

#
#  New homogeneous and isotropic cosmological model NcHICosmoDEXcdm 
#
cosmo = Nc.HICosmo.new_from_name (Nc.HICosmo, "NcHICosmoDEXcdm{'H0':<71.15>, 'Omegac':<0.218>, 'Omegab':<0.044>, 'w':<-1.01>, 'Omegac-fit':<1>, 'w-fit':<1>}")
cosmo.omega_x2omega_k ()
cosmo.param_set_by_name ("Omegak", 0.0);

(Found, w_i) = cosmo.param_index_from_name ("w")
cosmo.param_set_lower_bound (w_i, -3.0)
cosmo.param_set_upper_bound (w_i,  0.0)

(Found, Omegac_i) = cosmo.param_index_from_name ("Omegac")
cosmo.param_set_lower_bound (Omegac_i, 0.15)
cosmo.param_set_upper_bound (Omegac_i, 0.40)

#
#  New homogeneous and isotropic reionization object.
#
reion = Nc.HIReionCamb.new () 

#
#  New homogeneous and isotropic primordial object.
#
prim = Nc.HIPrimPowerLaw.new () 
prim.props.n_SA = 0.97

#
# Adding submodels to the main cosmological model.
#
cosmo.add_submodel (reion)
cosmo.add_submodel (prim)

#
# Calc objects
#
dist = Nc.Distance.new (z_max * 2.0)

# Uncomment below to use CLASS power spectrum
ps    = Nc.PowspecMLTransfer.new (Nc.TransferFuncEH.new ())
#ps   = Nc.PowspecMLCBE.new ()

ps.require_kmin (1.0e-6)
ps.require_kmax (1.0e2)

#
# Apply a tophat filter to the psml object, set best output interval.
#
psf  = Ncm.PowspecFilter.new (ps, Ncm.PowspecFilterType.TOPHAT)
psf.set_best_lnr0 ()

mulf = Nc.MultiplicityFunc.new_from_name ("NcMultiplicityFuncTinkerCrit{'Delta':<500.0>}")
mf   = Nc.HaloMassFunction.new (dist, psf, mulf)
mf.props.prec = 1.0e-3

cluster_m = Nc.ClusterMass.new_from_name ("NcClusterMassBenson{'M0':<3e14>, 'z0':<0.6>, 'signif-obs-min':<5.0>, 'Asz':<6.24>, 'Bsz':<1.33>, 'Csz':<0.83>, 'Dsz':<0.24>}")
cluster_z = Nc.ClusterRedshift.new_from_name ("NcClusterPhotozGaussGlobal{'pz-min':<%f>, 'pz-max':<%f>, 'z-bias':<0.0>, 'sigma0':<0.05>}" % (z_min, z_max))
cad       = Nc.ClusterAbundance.new (mf, None)

ncount = Nc.DataClusterNCount.new (cad)

mset = Ncm.MSet ()
mset.set (cosmo)
mset.set (cluster_m)
mset.set (cluster_z)

rng = Ncm.RNG.pool_get ("example_ca_sampling");
rng.set_seed (123)

ncount.init_from_sampling (mset, area * (math.pi / 180.0)**2, rng)
ncount.true_data (False)

dset = Ncm.Dataset.new ()
dset.append_data (ncount)

mset.prepare_fparam_map ()

prior = Ncm.MSetTransKernFlat.new ()
prior.set_mset (mset)
prior.set_prior_from_mset ()

#print "# sigma8 % 22.15g" % cosmo.sigma8 (psf)
print (ncount.get_desc ())
#print (ncount.get_lnM_obs()) 
#print (ncount.get_lnM_obs().ncols ()) 
#print (ncount.get_lnM_obs().nrows ()) 
#print (ncount.get_lnM_obs().get (0, 0))
#print (ncount.get_lnM_obs().get (1, 0))
#print (ncount.get_z_obs().ncols ()) 
#print (ncount.get_z_obs().nrows ()) 
#print (ncount.get_z_obs().get (0, 0))
#print (ncount.get_z_obs().get (1, 0))
#print ("Exitting...")
#exit()

abc = Nc.ABCClusterNCount.new (mset, prior, dset);
abc.props.epsilon = 1.5
abc.props.summary_type = Nc.ABCClusterNCountSummary.GAUSS_RBF # BIN_UNIFORM, GAUSS_RBF
abc.set_scale_cov (False)

#abc.set_bin_uniform (200, 200)
#q = Ncm.Vector.new_array ([0.05 * (i+1) for i in range(19)])
#abc.set_bin_quantile (q)
#abc.props.epsilon_update_type = Nc.ABCClusterNCountEpsilonUpdate.UNIFORM

abc.set_epsilon_update (0.75)
abc.set_nthreads (2)

abc.set_mtype (Ncm.FitRunMsgs.SIMPLE) #FULL, SIMPLE
abc.set_rng (rng)

Ncm.func_eval_set_max_threads (2)
Ncm.func_eval_log_pool_stats ()

#abc.set_data_file ("nc_cluster_run_n.fits")

abc.start_run ()
abc.run (1000)
abc.end_run ()

while True:
  abc.start_update ()
  abc.update ()
  abc.end_update ()
  if abc.get_accept_ratio () < 0.002:
    i = 0
    break
