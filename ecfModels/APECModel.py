# Using a "TBabs * APEC" Model

# TBabs：Tuebingen-Boulder ISM Absorption Model (Calculating the X-ray absorption by the ISM)
# APEC Emission Spectrum (see http://atomdb.org/)

from xspec import *
import numpy as np
from astropy.table import Table, hstack, vstack

# Cross Section Table set to vern:  Verner, Ferland, Korista, and Yakovlev 1996
Xset.xsect = 'vern'
# Solar Abundance Vector set to wilm:  Wilms, J., Allen, A. & McCray, R. ApJ 542 914 (2000) 
# (abundances are set to zero for those elements not included in the paper).
Xset.abund = 'wilm'

coll = Table(names=('T','nH','Abundance','Redshift','norm','E1','E2','flux','flux_k',))

e1_ = np.array([0.1, 0.5, 0.2])	# Lower Limit of the band
e2_ = np.array([2.4, 2.0, 2.3])	# Upper Limit of the band

nH = np.arange(20,21.01,0.1)				# Log Column Density (in unit of 10^-20 cm^-2)
temp = np.arange(-2,1.1,0.2)			# Temperature (in unit of keV): 0.01 - 10 keV
abundance = [0,0.3]							# The Abundance relative to the solar value
redshift = np.arange(0.0,1.01,0.1)			# Redshift Range: 0.00 - 1.00
norm = 1 									# Normalization


for i in range(len(nH)):
	for j in range(len(temp)):
		for k in range(len(abundance)):
			for h in range(len(redshift)):
				m1 = Model("TBabs*apec")
				m2 = AllModels(1)
				m1.TBabs.nH = 10**(nH[i]-22)
				m1.apec.kT = 10**temp[j]
				m1.apec.Abundanc = abundance[k]
				m1.apec.Redshift = redshift[h]
				m1.apec.norm = norm
				for s in range(len(e1_)):
					flux_range = str(e1_[s])+" "+str(e2_[s])
					AllModels.calcFlux(flux_range)
					s1 = AllModels(1)
					AllModels.calcFlux(flux_range)
					f_s = s1.flux[0]
					llim = np.round(e1_[s]/(1+redshift[h]),2)
					ulim = np.round(e2_[s]/(1+redshift[h]),2)
					flux_range_l = str(llim)+" "+str(ulim)
					AllModels.calcFlux(str(llim)+" "+str(ulim))
					f_s0 = s1.flux[0]
					coll.add_row([np.round(temp[j],1), np.round(nH[i],1), np.round(abundance[k],1), 
								  np.round(redshift[h],1), norm, e1_[s], e2_[s], f_s, f_s0])

coll.write('apec_all.fits', overwrite='True')
