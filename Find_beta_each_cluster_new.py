# -*- coding: utf-8 -*-
# Try and load some of the data and look at it

from numpy import *
from matplotlib.pyplot import *
from matplotlib.mlab import *
from scipy.stats import *
from lmfit import minimize, Parameters, Model, Parameter, report_fit
from scipy.interpolate import interp1d

from pandas import *

# Some constants:

k_b = 1.3806488
m_p = 1.6726178*(10**-4)
mean_molecular_weight = 1.  # For den temperatur vi har er T/my

G = 6.67384
M_sol = 1.98855
kpc = 3.08567758
H = 70.4  # Hubble parameter
rho_c = 3*H**2/(8*pi*G)*(kpc/M_sol) # Kritisk tæthed af universet


# Definer "Kappa", som du har fra simulationerne:

def difffunc(r,q,p,o,a):
    return q - o/((r+1.3)/a)**p

diff_parameters = [4.0158e+06,16.8419116,6.0739e+05,1.67671123]

# Definer S*r:

def S_term(r,alpha,a,c,d,v_virial,r_virial):
    
    v_p = -alpha/(((r/r_virial)**(-a) + c*(r/r_virial)**(a/2))**(1/a)-d)*v_virial
    
    v_p_derivative = v_virial/r_virial*alpha/a*(a/2*c*(r/r_virial)**(a/2-1)-a*(r/r_virial)**(-a-1))*(((r/r_virial)**(-a)+c*(r/r_virial)**(a/2))**(1/a-1))/(((r/r_virial)**(-a)+c*(r/r_virial)**(a/2))**(1/a) - d)**2
    
    return v_p*v_p_derivative*r + v_p*H*r + v_p_derivative*r**2*H

dm_parameters = [0.13911346,22.4412195,6.4119e-10,0.31590271]
gas_parameters = [0.17613725,16.8621591,1.3588e-07,0.22206032]

chosen_parameters = gas_parameters

# Lav en model for masseprofilen som du kan bruge til at bestemme densiteten med:

def mass_fit_4(r,a,b,c,d):
    
    return 10**14*(c*((r/a)**d*(r/a+1)**-b)**(1/d))

# Lav modeller for gassens densitet og temperatur der gør det lettere at tage logaritmiske afledte:

def dens_fit(r,a,b,c,d):
    
    return 10/(d*(r/a)**b*(1+r/a)**c)

def temp_fit(r,a,b,c,d):
    
    return 10**7/(d*(r/a)**b*(1+r/a)**c)


# Load all the data:


n_bins = 19  # Final number of bins per profile
n_bins1 = n_bins+1
n_profiles = 51  # Total number of profiles

binned_velocitydispersion_r = empty(shape=(n_bins,n_profiles))
mass_actual = empty(shape=(n_bins,n_profiles))
test = empty(shape=(n_bins,n_profiles))

variance_r_actual = empty(shape=(n_bins,n_profiles))
variance_r_calculated = empty(shape=(n_bins,n_profiles))

beta_simulation = empty(shape=(n_bins,n_profiles))
beta_def_matrix = empty(shape=(n_bins,n_profiles))
beta_jeans_matrix = empty(shape=(n_bins,n_profiles))


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
    
    radius = data_dm[:,0]
    rho_dm = data_dm[:,1]
    sigma_r = data_dm[:,3]*1000
    sigma_theta = data_dm[:,4]*1000
    sigma_phi = data_dm[:,5]*1000
    
    rho_gas = data_gas[:,1]
    T_gas = data_gas[:,2]
    
    binned_density_gas = empty(n_bins)
    binned_T_gas = empty(n_bins)
    
    binned_density_dm = empty(n_bins)
    rho_tot = rho_dm+rho_gas
    binned_density_tot = empty(n_bins)
    radius_central = empty(n_bins)
    radius_outer = zeros(n_bins+1)
    mass_inside = 0
    
    
    #binned_velocitydispersion_r = empty(n_bins)
    binned_velocitydispersion_theta = empty(n_bins)
    binned_velocitydispersion_phi = empty(n_bins)
    
    bind = len(rho_dm)/n_bins1
    bin_width = (radius[5]-radius[4])
    
    
    for newbins in range(0,n_bins):
        
        radius_central[newbins] = sum(radius[(1+newbins*bind):(1+bind+newbins*bind)])/bind
        binned_density_dm[newbins] = median(rho_dm[(1+newbins*bind):(1+bind+newbins*bind)])
        
        binned_velocitydispersion_r[newbins,nfile] = median(sigma_r[(1+newbins*bind):(1+bind+newbins*bind)])
        binned_velocitydispersion_theta[newbins] = median(sigma_theta[(1+newbins*bind):(1+bind+newbins*bind)])
        binned_velocitydispersion_phi[newbins] = median(sigma_phi[(1+newbins*bind):(1+bind+newbins*bind)])
        
        binned_density_gas[newbins] = median(rho_gas[(1+newbins*bind):(1+bind+newbins*bind)])
        binned_T_gas[newbins] = median(T_gas[(1+newbins*bind):(1+bind+newbins*bind)])
        
        beta_simulation[newbins,nfile] = 1 - (((binned_velocitydispersion_theta[newbins])**2+(binned_velocitydispersion_phi[newbins])**2)/(2*(binned_velocitydispersion_r[newbins,nfile])**2))
        
        binned_density_tot[newbins] = median(rho_tot[(1+newbins*bind):(1+bind+newbins*bind)]) #*((radius[(1+newbins*bind):(1+bind+newbins*bind)]+bin_width/2.)**3-(radius[(1+newbins*bind):(1+bind+newbins*bind)]-bin_width/2.)**3))/((radius[bind+newbins*bind]+0.5*bin_width)**3-(1+newbins*bind)**3)
        radius_outer[newbins+1] = radius[(newbins+1)*bind]+0.5*bin_width
        mass_shell = 4*pi/3*(radius_outer[newbins+1]**3-radius_outer[newbins]**3)*binned_density_tot[newbins]
        mass_inside = mass_inside + mass_shell
        
        mass_actual[newbins,nfile] = mass_inside
    
    r_virial = radius_central[int(n_bins/2.)]
    v_virial = sqrt(G*mass_actual[int(n_bins/2.),nfile]*M_sol/(r_virial*kpc))
    
    variance_r_actual[:,nfile] = binned_velocitydispersion_r[:,nfile]**2
    
    # Lav fit af densiteten og temperaturen:
    
    def chi2_dens_fit(params_dens,r,data):
        
        a = params_dens['a'].value
        b = params_dens['b'].value
        c = params_dens['c'].value
        d = params_dens['d'].value
        
        model = dens_fit(r,a,b,c,d)
        
        return (data-model)
        
    params_dens = Parameters()
    params_dens.add('a',value = 500.,min=0)
    params_dens.add('b',value = 1.)
    params_dens.add('c',value = 1.)
    params_dens.add('d',value = 1.)
    
    dens_fit_result = minimize(chi2_dens_fit,params_dens,args=(radius_central[2:n_bins],binned_density_gas[2:n_bins]))
    
    density_gas_fit = dens_fit(radius_central,params_dens['a'].value,params_dens['b'].value,params_dens['c'].value,params_dens['d'].value)
    
    def chi2_temp_fit(params_temp,r,data):
        
        a = params_temp['a'].value
        b = params_temp['b'].value
        c = params_temp['c'].value
        d = params_temp['d'].value
        
        model = temp_fit(r,a,b,c,d)
        
        return (data-model)
    
    params_temp = Parameters()
    params_temp.add('a',value = 500.,min=0)
    params_temp.add('b',value = 1.)
    params_temp.add('c',value = 1.)
    params_temp.add('d',value = 1.)
    
    temp_fit_result = minimize(chi2_temp_fit,params_temp,args=(radius_central,binned_T_gas,))
    
    T_gas_fit = temp_fit(radius_central,params_temp['a'].value,params_temp['b'].value,params_temp['c'].value,params_temp['d'].value)
    
    # Udregn de logarimitsk afledte for gassen:
    
    gamma_gas = empty(n_bins) # dlog_rho_dlog_r
    eta_gas = empty(n_bins)   # dlog_T_dlog_r
    
    gamma_gas[1:n_bins-1] = (log(density_gas_fit[1:n_bins-1]/density_gas_fit[0:n_bins-2])/log(radius_central[1:n_bins-1]/radius_central[0:n_bins-2])+log(density_gas_fit[2:n_bins]/density_gas_fit[1:n_bins-1])/log(radius_central[2:n_bins]/radius_central[1:n_bins-1]))/2
    gamma_gas[0] = (3*log(density_gas_fit[1]/density_gas_fit[0])/log(radius_central[1]/radius_central[0])+log(density_gas_fit[2]/density_gas_fit[1])/log(radius_central[2]/radius_central[1]))/4
    gamma_gas[n_bins-1] = (3*log(density_gas_fit[n_bins-1]/density_gas_fit[n_bins-2])/log(radius_central[n_bins-1]/radius_central[n_bins-2])+log(density_gas_fit[n_bins-2]/density_gas_fit[n_bins-3])/log(radius_central[n_bins-2]/radius_central[n_bins-3]))/4
    
    gamma_gas_wolfram = -(params_dens['a'].value*params_dens['b'].value+radius_central*(params_dens['b'].value+params_dens['c'].value))/(params_dens['a'].value+radius_central)
    
    eta_gas[1:n_bins-1] = (log(T_gas_fit[1:n_bins-1]/T_gas_fit[0:n_bins-2])/log(radius_central[1:n_bins-1]/radius_central[0:n_bins-2])+log(T_gas_fit[2:n_bins]/T_gas_fit[1:n_bins-1])/log(radius_central[2:n_bins]/radius_central[1:n_bins-1]))/2
    eta_gas[0] = (3*log(T_gas_fit[1]/T_gas_fit[0])/log(radius_central[1]/radius_central[0])+log(T_gas_fit[2]/T_gas_fit[1])/log(radius_central[2]/radius_central[1]))/4
    eta_gas[n_bins-1] = (3*log(T_gas_fit[n_bins-1]/T_gas_fit[n_bins-2])/log(radius_central[n_bins-1]/radius_central[n_bins-2])+log(T_gas_fit[n_bins-2]/T_gas_fit[n_bins-3])/log(radius_central[n_bins-2]/radius_central[n_bins-3]))/4
    
    eta_gas_wolfram = -(params_temp['a'].value*params_temp['b'].value+radius_central*(params_temp['b'].value+params_temp['c'].value))/(params_temp['a'].value+radius_central)
    
    # Udregn massen ved Euler's ligning:
    
    mass_term = -T_gas_fit*k_b/(m_p)*(gamma_gas_wolfram + eta_gas_wolfram) - S_term(radius_central,chosen_parameters[0],chosen_parameters[1],chosen_parameters[2],chosen_parameters[3],v_virial,r_virial)
    mass_calculated = mass_term*radius_central*kpc/(G*M_sol)
    
    def chi2_mass(params,r,data):
        
        a = params['a'].value
        b = params['b'].value
        c = params['c'].value
        d = params['d'].value
        
        model = mass_fit_4(r,a,b,c,d)
        
        return (data-model)
    
    params = Parameters()
    
    params.add('a',value = 1.,min=0.)
    params.add('b',value = 1,min=0.)
    params.add('c',value = 1.,max=3.)
    params.add('d',value = 1.)#,max=0.)
    
    result = minimize(chi2_mass,params,args=(radius_central,mass_calculated))
    
    #report_fit(params)
    
    mass_fit_result = mass_fit_4(radius_central,params['a'].value,params['b'].value,params['c'].value,params['d'].value)
    
    density_tot_calculated = empty(n_bins)
    density_tot_calculated[1:n_bins-1] = ((mass_fit_result[1:n_bins-1]-mass_fit_result[0:n_bins-2])/(radius_central[1:n_bins-1]-radius_central[0:n_bins-2])+(mass_fit_result[2:n_bins]-mass_fit_result[1:n_bins-1])/(radius_central[2:n_bins]-radius_central[1:n_bins-1]))/(2*4*pi*(radius_central[1:n_bins-1])**2)
    density_tot_calculated[0] = (3*(mass_fit_result[1]-mass_fit_result[0])/(radius_central[1]-radius_central[0])+(mass_fit_result[2]-mass_fit_result[1])/(radius_central[2]-radius_central[1]))/(4*4*pi*(radius_central[0])**2)
    density_tot_calculated[n_bins-1] = (3*(mass_fit_result[n_bins-1]-mass_fit_result[n_bins-2])/(radius_central[n_bins-1]-radius_central[n_bins-2])+(mass_fit_result[n_bins-2]-mass_fit_result[n_bins-3])/(radius_central[n_bins-2]-radius_central[n_bins-3]))/(4*4*pi*(radius_central[n_bins-1])**2)
    density_dm_calculated = density_tot_calculated - density_gas_fit
    
    #print eta_gas/eta_gas_wolfram
    
    #massfigure = figure(2+nfile)
    #plot(radius_central,mass_actual[:,nfile],'b',radius_central,mass_calculated,'g',radius_central,mass_fit_result,'r')
    #grid()
    #title(filelist_dm[nfile])
    #massfigure.show()
    
    # Find Hastighedsdispersionen ved at integrere:
    
    density_dm = density_dm_calculated
    T_dm = T_gas_fit + difffunc(radius_central/r_virial,diff_parameters[0],diff_parameters[1],diff_parameters[2],diff_parameters[3])
    
    psi_integrand = T_gas_fit*k_b/(m_p*mean_molecular_weight)*(gamma_gas_wolfram + eta_gas_wolfram) + 3*k_b/(m_p*mean_molecular_weight)*T_dm
    integrand = 10**-10*psi_integrand*density_dm*radius_central**2
    
    integral_piece = empty(n_bins)
    dispersion_integral = empty(n_bins)
    radial_velocity_variance = empty(n_bins)
    
    integral_piece[0] = 0.5*integrand[0]*radius_central[0]
    dispersion_integral[0] = integral_piece[0]
    radial_velocity_variance[0] = 10**10/(radius_central[0]**3*density_dm[0])*dispersion_integral[0]
    
    for n_integrate in range(1,n_bins):
        
        integral_piece[n_integrate] = 0.5*(integrand[n_integrate]+integrand[n_integrate-1])*(radius_central[n_integrate]-radius_central[n_integrate-1])
        dispersion_integral[n_integrate] = dispersion_integral[n_integrate-1] + integral_piece[n_integrate]
        radial_velocity_variance[n_integrate] = 10**10/(radius_central[n_integrate]**3*density_dm[n_integrate])*dispersion_integral[n_integrate]
    
    beta_inderst = 0.
    variance_r_inderst = (T_dm[0]) * 3*k_b/(mean_molecular_weight*m_p * (3 - 2 * beta_inderst))
    variance_r_yderst = (k_b/m_p*T_dm[n_bins-1])/3.
    radial_velocity_variance = radial_velocity_variance - radial_velocity_variance[n_bins-1] + variance_r_actual[n_bins-1,nfile] #sigma_r_inderst #+ 1.5*10**11#- radial_velocity_variance[0] + median_sigma_r_simulation[0]#+1.7*10**11#
    variance_r_calculated[:,nfile] = radial_velocity_variance
    
    #variancefigure = figure(1+nfile)
    #plot(radius_central,variance_r_actual[:,nfile],'b',radius_central,variance_r_calculated[:,nfile],'g')#,radius_central,variance_fit_result,'r')
    #grid()
    #title(filelist_dm[nfile])
    #variancefigure.show()
    
    
    # Find nu beta på to måder og sammenlign med den fra simulationen:
    
    #Først ved hjælp af definitionen af delta og beta:
    
    beta_definition = 1.5 * (1. - k_b/(mean_molecular_weight*m_p)*T_dm/radial_velocity_variance)
    
    # Så ved hjælp af Jean's ligning:
    
    #density_dm = binned_density_dm
    #radial_velocity_variance = variance_r_actual[:,nfile]
    #mass_fit_result = mass_actual[:,nfile]
    
    gamma_dm = empty(n_bins)
    eta_dm = empty(n_bins)
    
    gamma_dm[1:n_bins-1] = (log(density_dm[1:n_bins-1]/density_dm[0:n_bins-2])/log(radius_central[1:n_bins-1]/radius_central[0:n_bins-2])+log(density_dm[2:n_bins]/density_dm[1:n_bins-1])/log(radius_central[2:n_bins]/radius_central[1:n_bins-1]))/2
    gamma_dm[0] = (3*log(density_dm[1]/density_dm[0])/log(radius_central[1]/radius_central[0])+log(density_dm[2]/density_dm[1])/log(radius_central[2]/radius_central[1]))/4
    gamma_dm[n_bins-1] = (3*log(density_dm[n_bins-1]/density_dm[n_bins-2])/log(radius_central[n_bins-1]/radius_central[n_bins-2])+log(density_dm[n_bins-2]/density_dm[n_bins-3])/log(radius_central[n_bins-2]/radius_central[n_bins-3]))/4
    #gamma_dm = gamma_dm*radius_central/density_dm
    
    eta_dm[1:n_bins-1] = (log(radial_velocity_variance[1:n_bins-1]/radial_velocity_variance[0:n_bins-2])/log(radius_central[1:n_bins-1]/radius_central[0:n_bins-2])+log(radial_velocity_variance[2:n_bins]/radial_velocity_variance[1:n_bins-1])/log(radius_central[2:n_bins]/radius_central[1:n_bins-1]))/2
    eta_dm[0] = (3*log(radial_velocity_variance[1]/radial_velocity_variance[0])/log(radius_central[1]/radius_central[0])+log(radial_velocity_variance[2]/radial_velocity_variance[1])/log(radius_central[2]/radius_central[1]))/4
    eta_dm[n_bins-1] = (3*log(radial_velocity_variance[n_bins-1]/radial_velocity_variance[n_bins-2])/log(radius_central[n_bins-1]/radius_central[n_bins-2])+log(radial_velocity_variance[n_bins-2]/radial_velocity_variance[n_bins-3])/log(radius_central[n_bins-2]/radius_central[n_bins-3]))/4
    #eta_dm = eta_dm*radius_central/radial_velocity_variance
    
    beta_jeans = (-G*mass_fit_result*M_sol/(radius_central*kpc*radial_velocity_variance) - eta_dm - gamma_dm)/2.
    
    #print beta_definition
    
    beta_def_matrix[:,nfile] = beta_definition
    beta_jeans_matrix[:,nfile] = beta_jeans
    
    test[:,nfile] = beta_jeans/beta_simulation[:,nfile]
    
    print nfile +1
    
    eta_T_dm_outer = (3*log(T_dm[n_bins-1]/T_dm[n_bins-2])/log(radius_central[n_bins-1]/radius_central[n_bins-2])+log(T_dm[n_bins-2]/T_dm[n_bins-3])/log(radius_central[n_bins-2]/radius_central[n_bins-3]))/4
    
    #print gamma_dm[n_bins-1]
    
    gamma_outer = gamma_dm[n_bins-1]
    
    #if gamma_outer < -4. :
    #    gamma_outer = -4.
        
    print -(4+gamma_outer)/(1+gamma_outer)*2
    print -((2+eta_T_dm_outer)/(0.+eta_T_dm_outer))
    print -((3+eta_T_dm_outer)/(0.+eta_T_dm_outer))
    print -((2.5+eta_T_dm_outer)/(0.5+eta_T_dm_outer))
    
    
    

radius_relative = arange((0.0+2.0/len(rho_dm)),2.0-2.0/n_bins1,(2.0/n_bins1))+(1.0/n_bins1)

# Statistik:

median_test = empty(n_bins)
dispersion_test_upper = empty(n_bins)
dispersion_test_lower = empty(n_bins)

median_beta_actual = empty(n_bins)
dispersion_beta_actual_upper = empty(n_bins)
dispersion_beta_actual_lower = empty(n_bins)

median_beta_def = empty(n_bins)
dispersion_beta_def_upper = empty(n_bins)
dispersion_beta_def_lower = empty(n_bins)

median_beta_jeans = empty(n_bins)
dispersion_beta_jeans_upper = empty(n_bins)
dispersion_beta_jeans_lower = empty(n_bins)

median_variance_actual = empty(n_bins)
dispersion_variance_actual_upper = empty(n_bins)
dispersion_variance_actual_lower = empty(n_bins)

median_variance_calculated = empty(n_bins)
dispersion_variance_calculated_upper = empty(n_bins)
dispersion_variance_calculated_lower = empty(n_bins)

for n5 in range(0,n_bins):
    
    test_sorted = sorted(test[n5,:])
    median_test[n5] = test_sorted[int(len(test_sorted)*0.5)]
    sigma_file1 = int(len(test_sorted)*0.341)
    dispersion_test_upper[n5] = test_sorted[int(len(test_sorted)*0.5)+sigma_file1] - median_test[n5]
    dispersion_test_lower[n5] = median_test[n5] - test_sorted[int(len(test_sorted)*0.5)-sigma_file1]
    
    beta_actual_sorted = sorted(beta_simulation[n5,:])
    median_beta_actual[n5] = beta_actual_sorted[int(len(beta_actual_sorted)*0.5)]
    sigma_file1 = int(len(beta_actual_sorted)*0.341)
    dispersion_beta_actual_upper[n5] = beta_actual_sorted[int(len(beta_actual_sorted)*0.5)+sigma_file1] - median_beta_actual[n5]
    dispersion_beta_actual_lower[n5] = median_beta_actual[n5] - beta_actual_sorted[int(len(beta_actual_sorted)*0.5)-sigma_file1]
    
    beta_def_sorted = sorted(beta_def_matrix[n5,:])
    median_beta_def[n5] = beta_def_sorted[int(len(beta_def_sorted)*0.5)]
    sigma_file1 = int(len(beta_def_sorted)*0.341)
    dispersion_beta_def_upper[n5] = beta_def_sorted[int(len(beta_def_sorted)*0.5)+sigma_file1] - median_beta_def[n5]
    dispersion_beta_def_lower[n5] = median_beta_def[n5] - beta_def_sorted[int(len(beta_def_sorted)*0.5)-sigma_file1]
    
    beta_jeans_sorted = sorted(beta_jeans_matrix[n5,:])
    median_beta_jeans[n5] = beta_jeans_sorted[int(len(beta_jeans_sorted)*0.5)]
    sigma_file1 = int(len(beta_jeans_sorted)*0.341)
    dispersion_beta_jeans_upper[n5] = beta_jeans_sorted[int(len(beta_jeans_sorted)*0.5)+sigma_file1] - median_beta_jeans[n5]
    dispersion_beta_jeans_lower[n5] = median_beta_jeans[n5] - beta_jeans_sorted[int(len(beta_jeans_sorted)*0.5)-sigma_file1]
    
    variance_actual_sorted = sorted(variance_r_actual[n5,:])
    median_variance_actual[n5] = variance_actual_sorted[int(len(variance_actual_sorted)*0.5)]
    sigma_file1 = int(len(variance_actual_sorted)*0.341)
    dispersion_variance_actual_upper[n5] = variance_actual_sorted[int(len(variance_actual_sorted)*0.5)+sigma_file1] - median_variance_actual[n5]
    dispersion_variance_actual_lower[n5] = median_variance_actual[n5] - variance_actual_sorted[int(len(variance_actual_sorted)*0.5)-sigma_file1]
    
    variance_calculated_sorted = sorted(variance_r_calculated[n5,:])
    median_variance_calculated[n5] = variance_calculated_sorted[int(len(variance_calculated_sorted)*0.5)]
    sigma_file1 = int(len(variance_calculated_sorted)*0.341)
    dispersion_variance_calculated_upper[n5] = variance_calculated_sorted[int(len(variance_calculated_sorted)*0.5)+sigma_file1] - median_variance_calculated[n5]
    dispersion_variance_calculated_lower[n5] = median_variance_calculated[n5] - variance_calculated_sorted[int(len(variance_calculated_sorted)*0.5)-sigma_file1]
 

# Plotning:


varianceplot = figure(94)
axdisp2 = varianceplot.add_subplot(111)
#axdisp2.set_xscale("log", nonposx='clip')
#axdisp2.set_yscale("log", nonposy='clip')
axdisp2.errorbar(radius_relative,median_variance_actual,yerr=[dispersion_variance_actual_lower,dispersion_variance_actual_upper],fmt='ob')
axdisp2.errorbar(radius_relative,median_variance_calculated,yerr=[dispersion_variance_calculated_lower,dispersion_variance_calculated_upper],fmt='or')
grid()
#varianceplot.show()

betaplot = figure(95)
axdisp2 = betaplot.add_subplot(111)
axdisp2.errorbar(radius_relative,median_beta_actual,yerr=[dispersion_beta_actual_lower,dispersion_beta_actual_upper],fmt='ob')
axdisp2.errorbar(radius_relative,median_beta_def,yerr=[dispersion_beta_def_lower,dispersion_beta_def_upper],fmt='or')
axdisp2.errorbar(radius_relative,median_beta_jeans,yerr=[dispersion_beta_jeans_lower,dispersion_beta_jeans_upper],fmt='og')
grid()
betaplot.show()

testplot = figure(96)
errorbar(radius_relative,median_test,yerr=[dispersion_test_lower,dispersion_test_upper],fmt='ob')
grid()
#testplot.show()


# Median aller vægtet gennemsnit? Hmm...

