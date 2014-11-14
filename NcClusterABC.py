#!/usr/bin/python3

import threading
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm
from numpy import pi
import numpy as np
import random
import math
import sys

from scipy.stats.mstats import mquantiles

Ncm.cfg_init ()

class NcClusterABC (Ncm.ABC):

  def __init__ (self, mset, prior, dset, quant_list = [ 0.02, 0.09, 0.25, 0.5, 0.75, 0.91, 0.98]):
    super (NcClusterABC, self).__init__ (mset = mset, prior = prior, data_set = dset)
    tkern = Ncm.MSetTransKernGauss.new (0)
    tkern.set_mset (mset)
    tkern.set_cov_from_scale ()
    self.set_trans_kern (tkern)

    self.quant_list = quant_list
    self.epsilon = 1.6e20
    
    del tkern

  def do_data_summary (self):
    ndata = self.props.data_set.get_ndata ()
    ncount = self.props.data_set.get_data (0)
    self.true_data = ncount.using_true_data ()
    assert ndata == 1
    assert type (ncount) is Nc.DataClusterNCount
    
    data = []
    if self.true_data:
      mass_vec = ncount.get_lnM_true ()
      z_vec = ncount.get_z_true ()
      for i in range (mass_vec.len ()):
        data.append ([z_vec.get (i), mass_vec.get (i)])
    else:
      mass_mat = ncount.get_lnM_obs ()
      z_mat = ncount.get_z_obs ()
      for i in range (mass_mat.nrows ()):
        data.append ([z_mat.get (i, 0), mass_mat.get (i, 0)])

    data = np.array (data)
    self.dm_choose = [ min (data[:,1]) ]      
    mq = mquantiles (data[:,1], prob = self.quant_list)
    for it in mq:
      self.dm_choose.append (it)
    self.dm_choose.append (max (data[:,1]))
    
    data_bin = np.array([ [ item[0] for item in data if item[1] >= self.dm_choose[i] ] for i in range (len (self.dm_choose) - 1)])
    self.data_summary = [mquantiles( elem, prob=self.quant_list ) if len( elem ) > 0  else [ 0 for jj in self.quant_list]  for elem in data_bin]

    print (data_bin)
    print (self.data_summary)
    sys.stdout.flush ()
    
    del ncount
    del data
    del mq
    del data_bin
    
    return True
  
  def do_mock_distance (self, dset, theta, thetastar, rng):
    ncount = dset.get_data (0)
    data = []
    if ncount.get_len () > 0:
      if self.true_data:
        mass_vec = ncount.get_lnM_true ()
        z_vec = ncount.get_z_true ()
        for i in range (mass_vec.len ()):
          data.append ([z_vec.get (i), mass_vec.get (i)])
      else:
        mass_mat = ncount.get_lnM_obs ()
        z_mat = ncount.get_z_obs ()
        for i in range (mass_mat.nrows ()):
          data.append ([z_mat.get (i, 0), mass_mat.get (i, 0)])
    data = np.array (data)

    data_bin = np.array([ [ item[0] for item in data if item[1] >= self.dm_choose[i] ] for i in range (len (self.dm_choose) - 1)])
    mock_summary = [ mquantiles( elem, prob=self.quant_list ) if len( elem ) > 0  else [ 0 for jj in self.quant_list]  for elem in data_bin ]
    distance = [ np.sqrt( sum( [ ( self.data_summary[ i ][ j ] -  mock_summary[ i ][ j ] ) ** 2  for j in range( len( self.data_summary[ i ] ) ) ] ) ) for i in range( len( self.data_summary ) ) ]

    del ncount
    del data
    del data_bin
    del mock_summary

    return sum (distance)
  
  def do_update_tkern (self):
    cov = self.mcat.get_covar (None)
    cov.scale (2.0)
    
    self.tkern.set_cov (cov)
    
    dist_75 = self.get_dist_quantile (0.75)
    sys.stdout.write ("# NcClusterABC: ")
    for i in range (8):
      p = (15.0 * i) / 100.0
      if p > 1.0:
        p = 1.0
      sys.stdout.write ("[%2.0f%% %4.2f] " % (p * 100.0, self.get_dist_quantile (p)))
    print (" ")
    print ("# NcClusterABC: epsilon_t = %f, epsilon_tp1 = %f" % (self.epsilon, dist_75))
    self.epsilon = dist_75
    
    del cov
    return
    
  def do_distance_prob (self, distance):
    if distance < self.epsilon:
      return 1.0
    else:
      return 0.0

  def do_get_desc (self):
    return "NcClusterABC"
  
  def do_log_info (self):
    return self.dset.get_info ()

