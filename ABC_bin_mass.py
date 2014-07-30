from math import *
from gi.repository import GObject
import matplotlib.pyplot as plt
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm

import numpy
import random 
import sys 

from NCountSimul import *

############################################################
### Fiducial Cosmological parameters

z_min = 0.0              #minimum redshift
z_max = 2.0               #maximum redshift
H0 = 70.0                #Hubble parameter
Omegam = 0.25            #Dark matter density
OmegaL = 0.7             #Dark energy density
Omegab = 0.05            #Baryon density
Tgamma0 = 2.72           #
ns = 1.0                 #
sigma8 = 0.9             #
w = -1.0                 #Dark energy equation of state     


#Tolerance
epsilon = 20

#Number of tries
N_tries = 100

#mass bin
dm = [10**14, 10**14.25, 10**14.5, 10**14.75,  10**15, 10**15.25,  10**15.5, 10**15.75 ]


#output file
path1 = 'ABC_mass_bin.dat'

#sky area
area = 100

##############################################################

#new simulation object
ncount = NCountSimul (z_min, z_max, log ( dm[0] ), log ( dm[-1]* (10**0.25) ), area )

#Generate fiducial data
data_fid = numpy.array( ncount.simulation( z_max, H0, Omegab, Omegam, 1-Omegam, Tgamma0, ns, sigma8, w )[1] )


#vector to store surviving Om,w values
par_surv = []

op1 = open( path1, 'w')
op1.write( 'Om    w    sum_statistics    sum_pvalue    pvalue_all \n')

for k in range( N_tries ):				
    pvalue = False					
    while pvalue == False:		
               		
        #generate simulated data
        om1 =  0.0
        w1 = 0.0

        while om1 ==0 and w1 == 0.0:   
            om1 = random.uniform(0,1.0)
            w1 = random.uniform(-3.0,0.0)

        print 'om = ' + str( om1 ) + ',  w = ' + str( w1 )
        data_sim = numpy.array( ncount.simulation( z_max, H0, Omegab, om1, 1-om1, Tgamma0, ns, sigma8, w1 )[1] )

        #calculate summary statistics
        result = summary( data_fid, data_sim, dm )

        #calculate deviation
        difference = deviation( result, 0)

        p_value_list = result[:,1]
        print p_value_list

        if numpy.array( p_value_list ).any() < 0.05 or numpy.array( p_value_list ).any() > 1.0:
            pvalue = False
        elif numpy.array( p_value_list ).all() >= 0.05 and numpy.array( p_value_list ).any() <= 1.0:
            pvalue = True

        
	op1.write( str( om1 ) + '    ' + str( w1 ) + '    ' + str( difference ) + '    ' + str( sum( result[:,1] ) ) +  '    ' + str( pvalue ) +'\n' )
        print 'om = ' + str( om1 ) + ',  w = ' + str( w1 )  +  ', p-value[0] = ' + str( result[0][1] )	+ '    diff = ' + str( difference )
    
    if pvalue == True  and difference < epsilon:
        par_surv.append( [ om1, w1 ] )  
        

        print str( k + 1 ) + ', pvalue = ' + str( result[0][1] )

    
        

"""    
op1.close()  				

plt.figure()
plt.scatter( numpy.array( par_surv )[:,0], numpy.array( par_surv )[:,1] )
plt.xlim(0,1)
plt.ylim(-3,0)
plt.xlabel( 'Omega_M' )
plt.ylabel( 'w' )
plt.show()
"""
