# -*- coding: utf-8 -*-
# Try and load some of the data and look at it

from numpy import *
from matplotlib.pyplot import *
from matplotlib.mlab import *
from scipy.stats import *

#import scipy.optimize as optimization
#from scipy.optimize import fmin as simplex
from lmfit import minimize, Parameters, Model, Parameter, report_fit

# Some constants:

k_b = 1.3806488*(10**-3)
m_p = 1.6726178*(10**-7)
G = 6.67384*(10**-1)
M_sol = 1.98855*10
kpc = 3.08567758
H = 70.4  # Hubble parameter

# Load all the data:


n_bins = 50  # Final number of bins per profile
n_profiles = 51  # Total number of profiles

binned_density_dm = empty(shape=(n_bins,n_profiles))
binned_v_r_dm = empty(shape=(n_bins,n_profiles))

binned_v_r_gas = empty(shape=(n_bins,n_profiles))

binned_density_tot = empty(shape=(n_bins,n_profiles))


filelist_dm = ['prof_00001.dark','prof_00003.dark','prof_00004.dark','prof_00005.dark','prof_00006.dark','prof_00008.dark',
              'prof_00009.dark','prof_00010.dark','prof_00011.dark','prof_00012.dark','prof_00013.dark','prof_00014.dark',
              'prof_00015.dark','prof_00016.dark','prof_00017.dark','prof_00018.dark','prof_00020.dark','prof_00021.dark',
              'prof_00022.dark','prof_00023.dark','prof_00024.dark','prof_00025.dark','prof_00026.dark','prof_00027.dark',
              'prof_00028.dark','prof_00029.dark','prof_00030.dark','prof_00032.dark','prof_00033.dark','prof_00034.dark',
              'prof_00035.dark','prof_00036.dark','prof_00037.dark','prof_00038.dark','prof_00039.dark','prof_00040.dark',
              'prof_00043.dark','prof_00049.dark','prof_00054.dark','prof_00057.dark','prof_00061.dark','prof_00065.dark',
              'prof_00068.dark','prof_00073.dark','prof_00076.dark','prof_00079.dark','prof_00081.dark','prof_00085.dark',
              'prof_00088.dark','prof_00097.dark','prof_00099.dark']

filelist_gas = ['prof_00001.gas','prof_00003.gas','prof_00004.gas','prof_00005.gas','prof_00006.gas','prof_00008.gas',
               'prof_00009.gas','prof_00010.gas','prof_00011.gas','prof_00012.gas','prof_00013.gas','prof_00014.gas',
               'prof_00015.gas','prof_00016.gas','prof_00017.gas','prof_00018.gas','prof_00020.gas','prof_00021.gas',
               'prof_00022.gas','prof_00023.gas','prof_00024.gas','prof_00025.gas','prof_00026.gas','prof_00027.gas',
               'prof_00028.gas','prof_00029.gas','prof_00030.gas','prof_00032.gas','prof_00033.gas','prof_00034.gas',
               'prof_00035.gas','prof_00036.gas','prof_00037.gas','prof_00038.gas','prof_00039.gas','prof_00040.gas',
               'prof_00043.gas','prof_00049.gas','prof_00054.gas','prof_00057.gas','prof_00061.gas','prof_00065.gas',
               'prof_00068.gas','prof_00073.gas','prof_00076.gas','prof_00079.gas','prof_00081.gas','prof_00085.gas',
               'prof_00088.gas','prof_00097.gas','prof_00099.gas']



for nfile in range(0,n_profiles):
    
    data_dm = loadtxt(filelist_dm[nfile])
    
    data_gas = loadtxt(filelist_gas[nfile])
    
    radius = data_gas[:,0]
    
    rho_dm = data_dm[:,1]
    v_r_dm = data_dm[:,2]*1000#+H*radius
        
    rho_gas = data_gas[:,1]
    v_r_gas = data_gas[:,3]*1000#+H*radius
    
    bind = len(rho_dm)/n_bins
    rho_tot = rho_dm+rho_gas
    bin_width = (radius[5]-radius[4])
    
    radius_central = empty(n_bins)
    radius_outer = zeros(n_bins+1)
    mass_inside = 0
    v_virial = empty(n_bins)
    
    for newbins in range(0,n_bins):
        
        binned_v_r_dm[newbins,nfile] = sum(v_r_dm[(0+newbins*bind):(bind+newbins*bind)])/bind
        binned_density_dm[newbins,nfile] = sum(rho_dm[(0+newbins*bind):(bind+newbins*bind)])/bind
        
        binned_v_r_gas[newbins,nfile] = sum(v_r_gas[(0+newbins*bind):(bind+newbins*bind)])/bind
        
        binned_density_tot[newbins,nfile] = sum(rho_tot[(0+newbins*bind):(bind+newbins*bind)])/bind
        radius_central[newbins] = sum(radius[(0+newbins*bind):(bind+newbins*bind)])/bind
        radius_outer[newbins+1] = radius[(newbins+1)*bind-1]+0.5*bin_width
        mass_shell = 4*pi/3*(radius_outer[newbins+1]**3-radius_outer[newbins]**3)*binned_density_tot[newbins,nfile]
        mass_inside = mass_inside + mass_shell
        
        if newbins==(int(n_bins/2.)):
            mass_virial = mass_inside
        
        
    v_virial = sqrt(G*mass_virial*M_sol/(radius_central[int(n_bins/2.)]*kpc))
    
    binned_density_dm[:,nfile] = binned_density_dm[:,nfile]/binned_density_dm[int(n_bins/2-0.1),nfile] # Normaliserer tætheden til den værdi den har ved r200
    binned_v_r_dm[:,nfile] = binned_v_r_dm[:,nfile]/v_virial
    binned_v_r_gas[:,nfile] = binned_v_r_gas[:,nfile]/v_virial
    

radius_relative = arange(0.0,2.0,(2.0/n_bins))+(1.0/n_bins)



# Statistisk behandling af indfaldshastighederne:

median_vr_dm = empty(n_bins)
dispersion_vr_dm_upper = empty(n_bins)
dispersion_vr_dm_lower = empty(n_bins)

median_vr_gas = empty(n_bins)
dispersion_vr_gas_upper = empty(n_bins)
dispersion_vr_gas_lower = empty(n_bins)

for nbins3 in range(0,n_bins):
    
    vr_dm_sorted = sorted(binned_v_r_dm[nbins3,:])
    median_vr_dm[nbins3] = vr_dm_sorted[int(len(vr_dm_sorted)*0.5)]
    sigma_file = int(len(vr_dm_sorted)*0.341)
    dispersion_vr_dm_upper[nbins3] = vr_dm_sorted[int(len(vr_dm_sorted)*0.5)+sigma_file] - median_vr_dm[nbins3]
    dispersion_vr_dm_lower[nbins3] = median_vr_dm[nbins3] - vr_dm_sorted[int(len(vr_dm_sorted)*0.5)-sigma_file]
    
    
    vr_gas_sorted = sorted(binned_v_r_gas[nbins3,:])
    median_vr_gas[nbins3] = vr_gas_sorted[int(len(vr_gas_sorted)*0.5)]
    dispersion_vr_gas_upper[nbins3] = vr_gas_sorted[int(len(vr_gas_sorted)*0.5)+sigma_file] - median_vr_gas[nbins3]
    dispersion_vr_gas_lower[nbins3] = median_vr_gas[nbins3] - vr_gas_sorted[int(len(vr_gas_sorted)*0.5)-sigma_file]


# Lav en række funnktioner der ligner dataene lidt:

data_fit = median_vr_dm  # De data der fittes til
error_data_fit = (dispersion_vr_dm_upper+dispersion_vr_dm_lower)/2.

def martinafunc(r,alpha,a,b,c,d):
    
    return -alpha/((r**(-a) + c*(r)**b)**(1/a)-d)


def chi2_martina(params,r,data,dispersion_data):
    alpha = params['alpha'].value
    a = params['a'].value
    b = params['b'].value
    c = params['c'].value
    d = params['d'].value

    model = martinafunc(r,alpha,a,b,c,d)

    return (data-model)/dispersion_data

params = Parameters()

params.add('alpha',value = 0.5)
params.add('a',value = 3.,min=0.)
params.add('b',value = 2.,min=0.)
params.add('c',value = 0.05,min=0.)
params.add('d',value = 0.2,min=0.)

result = minimize(chi2_martina,params,args=(radius_relative,data_fit,error_data_fit))

#report_fit(params)

#print params['a'].value

best_fit_martina_dm = martinafunc(radius_relative,0.1131545,31.9231804,5.1973*10**-7,1.7389*10**-9,0.31414724)

best_fit_martina_dm_grav = martinafunc(radius_relative,0.11255906,12.8245615,0.5*12.8245615,2.2179*10**-6,0.31763796)

best_fit_martina_gas_grav = martinafunc(radius_relative,0.13669927,12.635266,0.5*12.635266,6.3028*10**-6,0.2453202)


#Plot:

velocityplot = figure(65)
axdisp2 = velocityplot.add_subplot(111)
axdisp2.errorbar(radius_relative,median_vr_dm,yerr=[dispersion_vr_dm_lower,dispersion_vr_dm_upper],fmt='o',ecolor='b')
axdisp2.errorbar(radius_relative,median_vr_gas,yerr=[dispersion_vr_gas_lower,dispersion_vr_gas_upper],fmt='ro')
axdisp2.plot(radius_relative,best_fit_martina_dm_grav,'g')
#axdisp2.plot(radius_relative,infall_guessed,'g*')


grid()

velocityplot.show()




