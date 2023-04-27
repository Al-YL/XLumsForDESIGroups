# Using a "TBabs * Power Law" Model

# TBabs：Tuebingen-Boulder ISM Absorption Model (Calculating the X-ray absorption by the ISM)

from xspec import *
import numpy as np
from astropy.table import Table, hstack, vstack

# Cross Section Table set to vern:  Verner, Ferland, Korista, and Yakovlev 1996
Xset.xsect = 'vern'
# Solar Abundance Vector set to wilm:  Wilms, J., Allen, A. & McCray, R. ApJ 542 914 (2000) 
# (abundances are set to zero for those elements not included in the paper).
Xset.abund = 'wilm'

coll = Table(names=('Slope','nH','Redshift','norm','E1','E2','flux','flux_k',))

e1_ = np.array([0.1, 0.5, 0.2])	# Lower Limit of the band
e2_ = np.array([2.4, 2.0, 2.3])	# Upper Limit of the band

phoidx = np.arange(1.4,2.61,0.6)			# Powerlaw Slope: 1.4 - 2.6
nH = np.arange(20,21.01,0.1)				# Log Column Density (in unit of 10^-20 cm^-2)
redshift = np.arange(0.0,1.01,0.1)			# Redshift Range: 0.00 - 1.00
norm = 1 									# Normalization

for i in range(len(nH)):
	for j in range(len(phoidx)):
		for h in range(len(redshift)):
			m1 = Model("TBabs*zpowerlw")
			m1.TBabs.nH = 10**(nH[i]-22)
			m1.zpowerlw.PhoIndex = phoidx[j]
			m1.zpowerlw.norm = norm
			m1.zpowerlw.Redshift = redshift[h]
			fset = FakeitSettings(response='RMF',arf='ARF',exposure=1000000,
				fileName='./FakeSpectrum/pwl_'+
						str(int(nH[i]*10))+'_'+
						str(int(phoidx[j]*10))+'_'+
						str(int(redshift[h]*10))+'.fits')
			AllData.fakeit(1, fset)
