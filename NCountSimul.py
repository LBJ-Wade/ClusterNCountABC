from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm
from numpy import pi
import numpy

from scipy.stats import anderson_ksamp

class NCountSimul:

  def __init__ (self, z_min, z_max, lnM_min, lnM_max, area):
    Ncm.cfg_init ()
    self.cosmo = Nc.HICosmo.new_from_name (Nc.HICosmo, "NcHICosmoDEXcdm")
    dist = Nc.Distance.new (z_max * 1.5)
    wp =  Nc.Window.new_from_name ("NcWindowTophat")
    tf = Nc.TransferFunc.new_from_name ("NcTransferFuncEH")
    vp = Nc.MatterVar.new (Nc.MatterVarStrategy.FFT, wp, tf)
    gf = Nc.GrowthFunc.new ()
    mulf = Nc.MultiplicityFunc.new_from_name ("NcMultiplicityFuncTinkerMean")
    mf = Nc.MassFunction.new (dist, vp, gf, mulf)
    
    cluster_m = Nc.ClusterMass.new_from_name ("NcClusterMassLnnormal{'lnMobs-min':<% 20.15g>, 'lnMobs-max':<% 20.15g>}" % (lnM_min, lnM_max))
    cluster_z = Nc.ClusterRedshift.new_from_name ("NcClusterPhotozGaussGlobal{'pz-min':<%f>, 'pz-max':<%f>, 'z-bias':<0.0>, 'sigma0':<0.03>}" % (z_min, z_max))
    cad = Nc.ClusterAbundance.new (mf, None, cluster_z, cluster_m)

    self.ncdata = Nc.DataClusterNCount.new (cad)

    self.mset = Ncm.MSet ()
    self.mset.set (self.cosmo)
    self.mset.set (cluster_m)

    cluster_m.props.bias = 0.0
    cluster_m.props.sigma = 0.2

    self.rng = Ncm.RNG.pool_get ("example_ca_sampling");

    self.ncdata.init_from_sampling (self.mset, cluster_z, cluster_m, area * (pi / 180.0)**2, self.rng)

    del dist
    del vp
    del gf
    del mulf
    del mf
    del cad
    del cluster_z
    del cluster_m

  def simulation (self, z_max, H0, Omegab, Omegam, OmegaL, Tgamma0, ns, sigma8, w):
    self.cosmo.props.H0      = H0
    self.cosmo.props.Omegab  = Omegab
    self.cosmo.props.Omegac  = Omegam
    self.cosmo.props.Omegax  = OmegaL
    self.cosmo.props.Tgamma0 = Tgamma0
    self.cosmo.props.ns      = ns
    self.cosmo.props.sigma8  = sigma8
    self.cosmo.props.w       = w

    self.ncdata.resample (self.mset, self.rng)

    lnM_true = self.ncdata.get_lnM_true ()
    z_true = self.ncdata.get_z_true ()
    lnM_obs = self.ncdata.get_lnM_obs ()
    z_obs = self.ncdata.get_z_obs ()

    nobjects = self.ncdata.get_len ()

    header = ['z_true', 'Mass_true', 'z_obs', 'Mass_obs']

    sim = []
    for i in range (nobjects):
        #sim.append ([ z_true.get(i), lnM_true.get (i), z_obs.get (i, 0), lnM_obs.get (i, 0)])
        sim.append ([ z_obs.get (i, 0), lnM_obs.get (i, 0)])


    del lnM_true
    del z_true
    del lnM_obs
    del z_obs
    
    return header, sim 

def summary( fid_data, sim_data, mass_bin ):


    #convert mass bins to natural log
    mass = [ numpy.log( item ) for item in mass_bin ]

    #separate fiducial data
    fid_data_bin = numpy.array([ [ item[0] for item in fid_data  if item[1] >= mass[ i ] ] for i in range( len( mass ) ) ])

    #separate simulated data
    sim_data_bin = numpy.array([ [ item[0] for item in sim_data  if item[1] >= mass[ i ] ] for i in range( len( mass ) ) ])

    for ii in range( len( fid_data_bin ) ):
        print 'sim = ' + str( len( sim_data_bin[ ii ] ) ) + '   fid = ' + str( len( fid_data_bin[ ii ] ) )

    #perform anderson darling test
    res = numpy.array( [ [ anderson_ksamp( [ fid_data_bin[ i ], sim_data_bin[ i ] ] )[ j ] for j in [0,2] ]  for i in range( len( fid_data_bin ) ) if ( len( fid_data_bin[ i ] ) > 0 and len( sim_data_bin[ i ] ) > 0 ) ])

    return res 

def deviation( summary_res, choice ):

    #choice -> which diagnostic to use

    valid_res = numpy.array( [ item for item in summary_res if str( item[0] ) != 'nan' and str( item[1] ) != 'nan' ])

    return sum( valid_res[:, choice ] )


