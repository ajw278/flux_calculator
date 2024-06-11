import os


os.environ["PYSYN_CDBS"] = "/home/awinter/Documents/pysyn_data/grp/redcat/trds/"

import stellar_evolution as se
import pysynphot as S
import matplotlib.pyplot as plt
import numpy as np

import warnings


#Lsol in erg/s
Lsol = 3.828e33

#Rsol in cm
Rsol = 6.957e10

def get_spectra(mstar, age, metallicity=0.0):
	
	# Get stellar properties
	Teff, log_g, log_L, R, star_mass = se.fetch_stellar_properties(mstar, age)
	
	# Compute the stellar spectrum using Castelli & Kurucz atmosphere models
	try:
		sp = S.Icat('ck04models', Teff, metallicity, log_g)
	except:
		print('Warning: using blackbody spectrum because stellar parameters outside of atmosphere model range')
		sp = S.BlackBody(Teff)
	#sp = S.Icat('k93models', Teff, metallicity, log_g)
	
	#Renormalize given the stellar luminosity 
	Ltot = np.trapz(sp.flux*np.pi*4.0*R*R*Rsol*Rsol, sp.wave)
	
	Lnorm = Lsol*10.**log_L / Ltot
	
	#if (Lnorm-1.)/Lnorm>0.5:
	#	raise warnings.warn('Luminosity of the spectra very different to the evoluton model')
	
	
	return sp.wave, sp.flux*Lnorm, R*Rsol


def compute_luminosity(wave, flux, Rstar, wavelength_start=0.0, wavelength_end = np.inf):
	"""
	Compute the luminosity of the star between given wavelengths.

	Parameters:
	wave (array): Wavelength array in Angstroms.
	flux (array): Flux array in erg/cm^2/s/Å.
	wavelength_start (float): Starting wavelength in Angstroms.
	wavelength_end (float): Ending wavelength in Angstroms.

	Returns:
	float: Luminosity in erg/s.
	"""
	# Mask to select the wavelength range
	mask = (wave >= wavelength_start) & (wave <= wavelength_end)

	# Integrate the flux over the selected wavelength range
	integrated_flux = np.trapz(flux[mask]*4.*np.pi*Rstar*Rstar, wave[mask])

	# Convert to luminosity (erg/s)
	# Note: The factor of 4πR^2 is already included in the flux normalization, 
	# so we don't need to include it again here.
	luminosity = integrated_flux

	return luminosity

def compute_fuv_luminosity_over_time(mstar, metallicity=0.0, age_start=1e1, age_end=1e7, age_res=40):
	ages = np.logspace(np.log10(age_start), np.log10(age_end), age_res)
	fuv_luminosities = np.zeros(ages.shape)

	for iage, age in enumerate(ages):
		wave, flux, Rstar = get_spectra(mstar, age, metallicity)
		#plt.plot(wave, flux)
		fuv_luminosity = compute_luminosity(wave, flux, Rstar, 910., 2070.0)  # FUV range in Angstroms
		fuv_luminosities[iage] =  fuv_luminosity
	
	#plt.xscale('log')
	return ages, fuv_luminosities

# Define a function to plot the evolution of FUV luminosity
def plot_fuv_luminosity_evolution(stellar_masses, metallicity=0.0):
	plt.figure(figsize=(12, 8))

	for mstar in stellar_masses:
		ages, fuv_luminosities = compute_fuv_luminosity_over_time(mstar, metallicity)
		plt.plot(ages, fuv_luminosities, label=f'{mstar} $M_\\odot$')

	plt.xlabel('Age (years)')
	plt.ylabel('FUV Luminosity (erg/s)')
	plt.yscale('log')
	plt.xscale('log')
	plt.ylim([1e33, 1e40])
	plt.legend()
	plt.grid(True)
	plt.show()
	
def plot_spectrum(mass, age, metallicity=0.0):
	
	wave, flux, R = get_spectra(target_mass, age, metallicity=metallicity)
	# Plot the spectrum
	plt.figure(figsize=(10, 6))
	plt.plot(wave, flux)
	plt.xlabel('Wavelength [Å]')
	plt.ylabel('Flux [erg cm$^2$ s$^{-1}$ Å$^{-1}$]')
	plt.xscale('log')
	plt.yscale('log')
	plt.grid(True)
	plt.show()

if __name__=='__main__':

	stellar_masses = [1, 3, 5, 10, 20, 50, 100]
	plot_fuv_luminosity_evolution(stellar_masses, metallicity=0.0)

	

