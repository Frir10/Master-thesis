# -*- coding: utf-8 -*-
# Try and load some of the data and look at it

from numpy import *
from matplotlib.pyplot import *
from matplotlib.mlab import *
from scipy.stats import *
from lmfit import minimize, Parameters, Model, Parameter, report_fit


# Some constants:

k_b = 1.3806488*(10**-3)
m_p = 1.6726178*(10**-7)


# Load all the data:


n_bins = 49  # Final number of bins per profile
n_bins1 = n_bins+1
n_profiles = 51  # Total number of profiles

binned_velocitydispersion_r = empty(shape=(n_bins,n_profiles))
binned_velocitydispersion_theta = empty(shape=(n_bins,n_profiles))
binned_velocitydispersion_phi = empty(shape=(n_bins,n_profiles))
#binned_v_r_dm = empty(shape=(n_bins,n_profiles))

binned_T_gas = empty(shape=(n_bins,n_profiles))
#binned_v_r_gas = empty(shape=(n_bins,n_profiles))

beta = empty(shape=(n_bins,n_profiles))
T_dm = empty(shape=(n_bins,n_profiles))
kappa = empty(shape=(n_bins,n_profiles))
temp_difference = empty(shape=(n_bins,n_profiles))

nfile = 0  # Counts the number of the file

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
    
    rho_dm = data_dm[:,1]
    #v_r_dm = data_dm[:,2]
    sigma_r = data_dm[:,3]*1000
    sigma_theta = data_dm[:,4]*1000
    sigma_phi = data_dm[:,5]*1000
    
    rho_gas = data_gas[:,1]
    T_gas = data_gas[:,2]
    #v_r_gas = data_gas[:,3]
    
    binned_density_dm = empty(n_bins)
    binned_density_gas = empty(n_bins)
    binned_density_tot = empty(n_bins)
    
    bind = len(rho_dm)/n_bins
    
    for newbins in range(0,n_bins):
        
        #binned_v_r_dm[newbins,nfile] = sum(v_r_dm[(0+newbins*bind):(bind+newbins*bind)])
        binned_density_dm[newbins] = median(rho_dm[(1+newbins*bind):(1+bind+newbins*bind)])
        binned_velocitydispersion_r[newbins,nfile] = median(sigma_r[(1+newbins*bind):(1+bind+newbins*bind)])
        binned_velocitydispersion_theta[newbins,nfile] = median(sigma_theta[(1+newbins*bind):(1+bind+newbins*bind)])
        binned_velocitydispersion_phi[newbins,nfile] = median(sigma_phi[(1+newbins*bind):(1+bind+newbins*bind)])
        
        binned_density_gas[newbins] = median(rho_gas[(1+newbins*bind):(1+bind+newbins*bind)])
        binned_T_gas[newbins,nfile] = median(T_gas[(1+newbins*bind):(1+bind+newbins*bind)])
        
        beta[newbins,nfile] = 1 - (((binned_velocitydispersion_theta[newbins,nfile])**2+(binned_velocitydispersion_phi[newbins,nfile])**2)/(2*(binned_velocitydispersion_r[newbins,nfile])**2))
        T_dm[newbins,nfile] = m_p/(3*k_b) * (binned_velocitydispersion_r[newbins,nfile]**2+binned_velocitydispersion_theta[newbins,nfile]**2+binned_velocitydispersion_phi[newbins,nfile]**2)
        kappa[newbins,nfile] = T_dm[newbins,nfile]/binned_T_gas[newbins,nfile]
        
    temp_difference[:,nfile] = (T_dm[:,nfile] - binned_T_gas[:,nfile])/binned_T_gas[24,nfile]
    
    #binned_density_dm[:,nfile] = binned_density_dm[:,nfile]/binned_density_dm[int(n_bins/2-0.1),nfile] # Normaliserer tætheden til den værdi den har ved r200
    #binned_v_r_dm[:,nfile] = binned_v_r_dm[:,nfile]/binned_v_r_dm[int(n_bins/2-0.1),nfile]
    


radius_relative = arange((0.0+2.0/len(rho_dm)),2.0-2.0/n_bins1,(2.0/n_bins1))+(1.0/n_bins1)

    
# Statistisk behandling af temperaturerne

median_T_gas = empty(n_bins)
dispersion_T_gas_upper = empty(n_bins)
dispersion_T_gas_lower = empty(n_bins)

median_T_dm = empty(n_bins)
dispersion_T_dm_upper = empty(n_bins)
dispersion_T_dm_lower = empty(n_bins)

for nbins3 in range(0,n_bins):
    
    T_gas_sorted = sorted(binned_T_gas[nbins3,:])
    median_T_gas[nbins3] = T_gas_sorted[int(len(T_gas_sorted)*0.5)]
    sigma_file = int(len(T_gas_sorted)*0.341)
    dispersion_T_gas_upper[nbins3] = T_gas_sorted[int(len(T_gas_sorted)*0.5)+sigma_file] - median_T_gas[nbins3]
    dispersion_T_gas_lower[nbins3] = median_T_gas[nbins3] - T_gas_sorted[int(len(T_gas_sorted)*0.5)-sigma_file]
    
    ##print dispersion_T_gas_upper[nbins3]/median_T_gas[nbins3]
    
    T_dm_sorted = sorted(T_dm[nbins3,:])
    median_T_dm[nbins3] = T_dm_sorted[int(len(T_dm_sorted)*0.5)]
    dispersion_T_dm_upper[nbins3] = T_dm_sorted[int(len(T_dm_sorted)*0.5)+sigma_file] - median_T_dm[nbins3]
    dispersion_T_dm_lower[nbins3] = median_T_dm[nbins3] - T_dm_sorted[int(len(T_dm_sorted)*0.5)-sigma_file]
    

temperatureplot = figure(91)
plot(radius_relative,median_T_gas,'*b',radius_relative,median_T_dm,'*r')
grid()

#temperatureplot.show()

# Plot forskellen i temperatur:

median_diff = empty(n_bins)
dispersion_diff_upper = empty(n_bins)
dispersion_diff_lower = empty(n_bins)

for nbins44 in range(0,n_bins):
    
    diff_sorted = sorted(temp_difference[nbins44,:])
    median_diff[nbins44] = diff_sorted[int(len(diff_sorted)*0.5)]
    
    dispersion_diff_upper[nbins44] = diff_sorted[int(len(diff_sorted)*0.5)+sigma_file] - median_diff[nbins44]
    dispersion_diff_lower[nbins44] = median_diff[nbins44] - diff_sorted[int(len(diff_sorted)*0.5)-sigma_file]
    

differenceplot = figure(11)

errorbar(radius_relative,median_diff,yerr=[dispersion_diff_lower,dispersion_diff_upper],fmt='ob')
grid()

#differenceplot.show()

# Plot kappa:

median_kappa = empty(n_bins)
dispersion_kappa_upper = empty(n_bins)
dispersion_kappa_lower = empty(n_bins)

for nbins4 in range(0,n_bins):
    
    kappa_sorted = sorted(kappa[nbins4,:])
    median_kappa[nbins4] = kappa_sorted[int(len(kappa_sorted)*0.5)]
    #print median_kappa[nbins4]
    dispersion_kappa_upper[nbins4] = kappa_sorted[int(len(T_dm_sorted)*0.5)+sigma_file] - median_kappa[nbins4]
    dispersion_kappa_lower[nbins4] = median_kappa[nbins4] - kappa_sorted[int(len(T_dm_sorted)*0.5)-sigma_file]
    
# Fit en funktion til kappa, som kan bruges som parametrisering:

dispersion_kappa = (dispersion_kappa_upper+dispersion_kappa_lower)/2.

f = 0.84811

def kappafunc(r,a,b,c,d):
    
    return 1/(d*(r/a)**b*(1+(r/a)**f)**(c/f))

def chi2_kappafunc(params_kappa,r,data,dispersion_data):
    
    a = params_kappa['a'].value
    b = params_kappa['b'].value
    c = params_kappa['c'].value
    d = params_kappa['d'].value
    
    model = kappafunc(r,a,b,c,d)
    
    return (data-model)/dispersion_data

params_kappa = Parameters()
params_kappa.add('a',value = 1.)#,max=0.01)
params_kappa.add('b',value = 1.)
params_kappa.add('c',value = 1.)
params_kappa.add('d',value = 1.)

kappa_result = minimize(chi2_kappafunc,params_kappa,args=(radius_relative,median_kappa,dispersion_kappa))
report_fit(params_kappa)
kappa_fit = kappafunc(radius_relative,params_kappa['a'].value,params_kappa['b'].value,params_kappa['c'].value,params_kappa['d'].value)

kappafitfigure = figure(14)
axdisp2 = kappafitfigure.add_subplot(111)
axdisp2.errorbar(radius_relative,median_kappa,yerr=[dispersion_kappa_lower,dispersion_kappa_upper],fmt='o',ecolor='b')
axdisp2.plot(radius_relative,kappa_fit,'g')
grid()
kappafitfigure.show()


# Statistisk behandling af hastighedsdispersionerne

median_sigma_r = empty(n_bins)
dispersion_sigma_r_upper = empty(n_bins)
dispersion_sigma_r_lower = empty(n_bins)

median_sigma_t = empty(n_bins)
dispersion_sigma_t_upper = empty(n_bins)
dispersion_sigma_t_lower = empty(n_bins)

for nbins26 in range(0,n_bins):
    
    sigma_r_sorted = sorted(binned_velocitydispersion_r[nbins26,:])
    median_sigma_r[nbins26] = sigma_r_sorted[int(len(sigma_r_sorted)*0.5)]
    sigma_file = int(len(sigma_r_sorted)*0.341)
    dispersion_sigma_r_upper[nbins26] = sigma_r_sorted[int(len(sigma_r_sorted)*0.5)+sigma_file] - median_sigma_r[nbins26]
    dispersion_sigma_r_lower[nbins26] = median_sigma_r[nbins26] - sigma_r_sorted[int(len(sigma_r_sorted)*0.5)-sigma_file]
    
    ##print dispersion_T_gas_upper[nbins3]/median_T_gas[nbins3]
    
    sigma_t_sorted = sorted(binned_velocitydispersion_theta[nbins26,:])
    median_sigma_t[nbins26] = sigma_t_sorted[int(len(T_dm_sorted)*0.5)]
    dispersion_sigma_t_upper[nbins26] = sigma_t_sorted[int(len(sigma_t_sorted)*0.5)+sigma_file] - median_sigma_t[nbins26]
    dispersion_sigma_t_lower[nbins26] = median_sigma_t[nbins26] - sigma_t_sorted[int(len(sigma_t_sorted)*0.5)-sigma_file]
    

dispersionplot = figure(99)
plot(radius_relative,median_sigma_r,'*r',radius_relative,median_sigma_t,'*b')
grid()

#dispersionplot.show()


# Plot beta:

median_beta = empty(n_bins)
dispersion_beta_upper = empty(n_bins)
dispersion_beta_lower = empty(n_bins)

for nbins5 in range(0,n_bins):
    
    beta_sorted = sorted(beta[nbins5,:])
    median_beta[nbins5] = beta_sorted[int(len(beta_sorted)*0.5)]
    dispersion_beta_upper[nbins5] = beta_sorted[int(len(beta_sorted)*0.5)+sigma_file] - median_beta[nbins5]
    dispersion_beta_lower[nbins5] = median_beta[nbins5] - beta_sorted[int(len(beta_sorted)*0.5)-sigma_file]
    

betaplot = figure(47)
errorbar(radius_relative,median_beta,yerr=[dispersion_beta_lower,dispersion_beta_upper],fmt='o')
grid()

#betaplot.show()


# Fit en funktion til beta for at parametrisere den:

def betafunc(r,a,b,c,d,r0):  # Den funktion vi fitter til data
    
    return b + a*r - c/(1+exp(-d*(r-r0)))
    

def residual_beta(params,r,data,dispersion_data): # Funktion der beregner residualet der skal minimeres
    
    a = params['a'].value
    b = params['b'].value
    c = params['c'].value
    d = params['d'].value
    r0 = params['r0'].value
    
    model = betafunc(r,a,b,c,d,r0)

    return (data-model)/dispersion_data

params = Parameters()

params.add('a',value = 0.1,min=0.)
params.add('b',value = 0.15,min=0.)
params.add('c',value = 1.5,min=0.)
params.add('d',value = 1.,min=0.)
params.add('r0',value = 1.5,min=0.)

result = minimize(residual_beta,params,args=(radius_relative,median_beta,dispersion_beta_upper))

#report_fit(params)

fit_curve_beta = betafunc(radius_relative,params['a'].value,params['b'].value,params['c'].value,params['d'].value,params['r0'].value)


betafitplot = figure(66)
axdisp2 = betafitplot.add_subplot(111)
axdisp2.errorbar(radius_relative,median_beta,yerr=[dispersion_beta_lower,dispersion_beta_upper],fmt='o',ecolor='r')
axdisp2.plot(radius_relative,fit_curve_beta,'g')
grid()

#betafitplot.show()


# Fit en funktion til forskellen mellem T_dm og T_gas:


def difffunc(r,q,p,o,a):
    
    return q - o/((r+1.3)/a)**p


def residual_diff(params_diff,r,data,dispersion_data):
    
    q = params_diff['q'].value
    p = params_diff['p'].value
    o = params_diff['o'].value
    a = params_diff['a'].value
    
    model = difffunc(r,q,p,o,a)

    return (data-model)/dispersion_data

params_diff = Parameters()

params_diff.add('q',value = 2.*10**6,min=0.)
params_diff.add('p',value = 0.)
params_diff.add('o',value = 10**6,min=0.)
params_diff.add('a',value = 1.,min=0.)

result_diff = minimize(residual_diff,params_diff,args=(radius_relative,median_diff,dispersion_diff_upper))

#report_fit(params_diff)

fit_curve_diff = difffunc(radius_relative,params_diff['q'].value,params_diff['p'].value,params_diff['o'].value,params_diff['a'].value)


difffitplot = figure(62)
axdisp2 = difffitplot.add_subplot(111)
axdisp2.errorbar(radius_relative,median_diff,yerr=[dispersion_diff_lower,dispersion_diff_upper],fmt='o',color='b')
axdisp2.plot(radius_relative,fit_curve_diff,'g')

grid()

#difffitplot.show()

print_out = empty(shape=(n_bins,8))

print_out[:,0] = radius_relative

print_out[:,1] = median_kappa
print_out[:,2] = dispersion_kappa_lower
print_out[:,3] = dispersion_kappa_upper

print_out[:,4] = kappa_fit

print_out[:,5] = median_diff
print_out[:,6] = dispersion_diff_lower
print_out[:,7] = dispersion_diff_upper

savetxt('RAMSES_kappa.txt',print_out)

print_out2 = empty(shape=(n_bins,8))

print_out2[:,0] = radius_relative

print_out2[:,1] = median_beta
print_out2[:,2] = dispersion_beta_lower
print_out2[:,3] = dispersion_beta_upper

savetxt('RAMSES_beta.txt',print_out2)
