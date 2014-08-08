import sys 
sys.path.append('/home/emille/Dropbox/WGC/ABC/code')

from math import *
from gi.repository import GObject
import matplotlib.pyplot as plt
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm

from scipy.stats.mstats import mquantiles
from scipy.stats import norm

import numpy
import random 

from NCountSimul import *

import pylab as plt

############################################################
### Fiducial Cosmological parameters: ref. arXiv:1203.5775 (table 5. wCDM CMB+BAO+H0+SNeIa+SPTcl), Tgamma0 and ns are not given.

z_min = 0.3              #minimum redshift
z_max = 1.32             #maximum redshift
H0 = 71.15               #Hubble parameter
Omegam = 0.262           #Dark matter density
Omegab = 0.0439          #Baryon density
Tgamma0 = 2.725          #Radiation temperature today
ns = 0.97                #spectral index 
sigma8 = 0.807           #sigma8
w = -1.01                #Dark energy equation of state     

#mass bin
#dm = [5*10**13, 10**14, 10**14.25, 10**14.5, 10**14.75,  10**15, 10**15.25,  10**15.5, 10**15.75 ]
dm = [10**14.3, 10**14.5, 10**14.7,  10**14.9, 10**15.1, 10**15.3,  10**15.5, 10**15.7 ]

#quantile list
quant_list = [ 0.02, 0.09, 0.25, 0.5, 0.75, 0.91, 0.98]

#sky area
area = 2500

#output file
path1 = 'sequential_ABC_res.dat'


#######
hyper = [ [ 0.3, 0.1], [0.6, 0.3], [ -0.8, 0.5 ] ]


time_steps = 10

N=500		        #particle sample size													
epsilon_ini = 1.0		#Starting tolerance


##############################################################

#new simulation object
ncount = NCountSimul (z_min, z_max, log ( dm[0] ), log ( 10**16 ), area )

#Generate fiducial data
data_fid = numpy.array( ncount.simulation( z_max, H0, Omegab, Omegam, 1-Omegam, Tgamma0, ns, sigma8, w )[1] )

#Calculate summary statistics for fiducial data
summ_fid = summary_quantile( data_fid, dm, quant_list )


#define weights vector   
weights = [numpy.array([ 1.0/N for j in range( N ) ]) for k1 in range( len( hyper ) ) ] 

#generate first set of simulation models
par_surv = choose_surv_par( summ_fid, dm, quant_list, epsilon_ini, N, hyper, 'normal', z_min, z_max, area, ncount  )

#define list to store thresholds
epsilon_list = [ epsilon_ini ]
  

for l in range( time_steps ):

    print 'l = ' + str( l )

    #define new limit for epsilon and keep it
    epsilon_list.append( mquantiles( [ par_surv[ j1 ][-1] for j1 in range( len( par_surv ) ) ], [0.75] )[0]  )

    #calculate the mean of the parameters according to given weights
    mean_prev = [ sum( [ par_surv[ j1 ][ k1 ] * weights[ k1 ][ j1 ] for j1 in range( N ) ] ) for k1 in range( len( hyper ) ) ]  

    #calculate variance for parameters given weights
    var_prev = [ sum( [ ( (par_surv[ j1 ][ l ] - mean_prev[ l ]) ** 2 )*weights[ l ][ j1 ] for j1 in range( N ) ])  for l in range( len( hyper ) ) ]

    #choose one of the instances of simulation set
    indx = random.randint( 0, N - 1 )

    #define mean as the selected instance and standard deviation given by the entire simulation set
    hyper2 = [ [ par_surv[ indx ][ k ], numpy.sqrt( 2*var_prev[ k ] ) ] for k in range( len( hyper ) ) ]
   
    #generate 100 simulations with the required tolerance and given the distribution defined in hyper2
    par_surv2 = choose_surv_par( summ_fid, dm, quant_list, epsilon_list[-1], N, hyper2, 'normal', z_min, z_max, area, ncount  ) 

    #update the weights of the last simulation given the mean and standard deviation from the previous set of simulations
    denominator = [ [ sum( [ weights[ l1 ][ i ]*norm.pdf( par_surv2[ n ][ l1 ], loc=par_surv[ i ][ l1 ], scale=numpy.sqrt( 2*var_prev[ l1 ] )   ) for i in range( N ) ]) for n in range( N ) ] for l1 in range( len( hyper ) ) ]

    numerator = [ [ norm.pdf( par_surv2[ j ][ l1 ], loc=0, scale=hyper[ l1 ][1] )  for j in range( N ) ] for l1 in range( len( hyper ) ) ] 

    weights_new = [ [ numerator[ j1 ][ l1 ]/denominator[ j1 ][ l1 ]  for l1 in range( N ) ] for j1 in range( len( hyper ) ) ] 

    weights = [ [ weights_new[ k ][ l1 ]/sum( weights_new[ k ] ) for l1 in range( N ) ] for k in range( len( hyper ) ) ]

    del par_surv
    
    #update the simulation set
    par_surv = par_surv2

    del par_surv2 
    print '        [om, ol, w, summ_stat] = ' + str( par_surv[0][0] )


#write results to file
op1 = open( path1, 'w')
op1.write( 'om       ol       w      summ_stat\n' )
for elem in par_surv:
    for item in elem:
        op1.write( str( item ) + '    ' )
    op1.write( '\n' )
op1.close()

#plot results
plt.figure()
plt.subplot( 1,2,1)
plt.scatter( numpy.array( par_surv )[:,0], numpy.array( par_surv )[:,1] )
plt.xlabel( 'om' )
plt.ylabel( 'ol' )
plt.xlim( 0,1)
plt.ylim(0,1)

plt.subplot( 1,2,2)
plt.scatter( numpy.array( par_surv )[:,0], numpy.array( par_surv )[:,2] )
plt.xlabel( 'om' )
plt.ylabel( 'w' )
plt.xlim( 0,1)
plt.ylim(-3,0)

plt.show() 

