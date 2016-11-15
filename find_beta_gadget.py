# -*- coding: utf-8 -*-
# Try and load some of the data and look at it

from numpy import *
from matplotlib.pyplot import *
from matplotlib.mlab import *
from scipy.stats import *
from lmfit import minimize, Parameters, Model, Parameter, report_fit
from scipy.interpolate import interp1d


# Some constants:

k_b = 1.3806488
m_p = 1.6726178*(10**-4)
mean_molecular_weight = 1.  # For den temperatur vi har er T/my

G = 6.67384
M_sol = 1.98855
kpc = 3.08567758
H = 72.  # Hubble parameter
rho_c = 3*H**2/(8*pi*G)*(kpc/M_sol) # Kritisk tæthed af universet


# Parametre

n_bins = 36  # Final number of bins per profile
n_bins1 = 50
n_profiles = 29  # Total number of profiles

mass_actual = zeros(shape=(n_bins,n_profiles))

variance_r_actual = zeros(shape=(n_bins,n_profiles))
variance_r_calculated = zeros(shape=(n_bins,n_profiles))

beta_actual = zeros(shape=(n_bins,n_profiles))
beta_def_matrix = zeros(shape=(n_bins,n_profiles))
beta_jeans_matrix = zeros(shape=(n_bins,n_profiles))

beta_def_matrix_less = zeros(shape=(n_bins,n_profiles))
beta_jeans_matrix_less = zeros(shape=(n_bins,n_profiles))

beta_def_matrix_more = zeros(shape=(n_bins,n_profiles))
beta_jeans_matrix_more = zeros(shape=(n_bins,n_profiles))

# Funktioner der skal bruges:

def kappafunc(r,a,b,c,d):
    
    return 1/(d*(r/a)**b*(1+r/a)**c)

kappa_func_parameters = [2.08418728,-0.36177877,2.43918831,0.18572878]

f = 10.

def kappafunc2(r,a,b,c,d):
    
    return 1/(d*(r/a)**b*(1+(r/a)**f)**(c/f))

kappa_func_parameters2 = [0.45329545,-0.16923136,0.53818309,0.48156340]

# Definer S*r:

def S_term(r,alpha,a,c,d,v_virial,r_virial):
    
    v_p = -alpha/(((r/r_virial)**(-a) + c*(r/r_virial)**(a/2))**(1/a)-d)*v_virial
    
    v_p_derivative = v_virial/r_virial*alpha/a*(a/2*c*(r/r_virial)**(a/2-1)-a*(r/r_virial)**(-a-1))*(((r/r_virial)**(-a)+c*(r/r_virial)**(a/2))**(1/a-1))/(((r/r_virial)**(-a)+c*(r/r_virial)**(a/2))**(1/a) - d)**2
    
    return v_p*v_p_derivative*r + v_p*H*r + v_p_derivative*r**2*H

dm_parameters = [0.20988822,19.,1.5417e-05,0.30]
gas_parameters = [0.16249621,10.,0.00069713,0.35478008]

# Lav en model for masseprofilen som du kan bruge til at bestemme densiteten med:

def mass_fit_4(r,a,b,c):
    
    return 10**14*(c*((r/a)*(r/a+1)**b))

# Lav modeller for gassens densitet og temperatur der gør det lettere at tage logaritmiske afledte:

def dens_fit(r,a,b,c,d):
    
    return 10/(d*(r/a)**b*(1+r/a)**c)

def temp_fit(r,a,b,c,d):
    
    return 10**7/(d*(r/a)**b*(1+r/a)**c)

# Data:


filelist_dm = ['prof_01_dm.txt','prof_02_dm.txt','prof_03_dm.txt','prof_04_dm.txt','prof_05_dm.txt','prof_06_dm.txt',
              'prof_07_dm.txt','prof_08_dm.txt','prof_09_dm.txt','prof_10_dm.txt','prof_11_dm.txt','prof_12_dm.txt',
              'prof_13_dm.txt','prof_14_dm.txt','prof_15_dm.txt','prof_16_dm.txt','prof_17_dm.txt','prof_18_dm.txt',
              'prof_19_dm.txt','prof_20_dm.txt','prof_21_dm.txt','prof_22_dm.txt','prof_23_dm.txt','prof_24_dm.txt',
              'prof_25_dm.txt','prof_26_dm.txt','prof_27_dm.txt','prof_28_dm.txt','prof_29_dm.txt']
              

filelist_gas = ['prof_01_gas.txt','prof_02_gas.txt','prof_03_gas.txt','prof_04_gas.txt','prof_05_gas.txt','prof_06_gas.txt',
              'prof_07_gas.txt','prof_08_gas.txt','prof_09_gas.txt','prof_10_gas.txt','prof_11_gas.txt','prof_12_gas.txt',
              'prof_13_gas.txt','prof_14_gas.txt','prof_15_gas.txt','prof_16_gas.txt','prof_17_gas.txt','prof_18_gas.txt',
              'prof_19_gas.txt','prof_20_gas.txt','prof_21_gas.txt','prof_22_gas.txt','prof_23_gas.txt','prof_24_gas.txt',
              'prof_25_gas.txt','prof_26_gas.txt','prof_27_gas.txt','prof_28_gas.txt','prof_29_gas.txt']


inner_limit = 21


for nfile in range(0,n_profiles):
    
    data_dm = loadtxt(filelist_dm[nfile])
    
    data_gas = loadtxt(filelist_gas[nfile])
    
    radius = data_dm[:,0]
    rho_dm = data_dm[:,1]
    sigma_r = data_dm[:,2]*1000#/sqrt(2.)
    sigma_t = data_dm[:,3]*1000*sqrt(2)
    
    rho_gas = data_gas[:,1]
    T_gas = data_gas[:,2]
    
    binned_density_gas = zeros(n_bins)
    binned_density_tot = zeros(n_bins)
    binned_T_gas = zeros(n_bins)
    binned_velocitydispersion_r = empty(n_bins)
    binned_velocitydispersion_t = zeros(n_bins)
    
    bind = len(rho_dm)/n_bins1
    rho_tot = rho_dm+rho_gas
    bin_width = (radius[23+bind]-radius[23])
    
    radius_central = empty(n_bins)
    radius_outer = zeros(n_bins+1)
    mass_inside = zeros(n_bins+1)
    
    #print T_gas[inner_limit+(n_bins-10)*bind]
    
    hej = 0
    
    for newbins in range(0,n_bins):
        
        binned_density_tot[newbins] = median(rho_tot[(inner_limit+newbins*bind):(inner_limit+bind+newbins*bind)])
        
        if binned_density_tot[newbins] == 0:
            binned_density_tot[newbins] = rho_tot[40]
        
        radius_central[newbins] = sum(radius[(inner_limit+newbins*bind):(inner_limit+bind+newbins*bind)])/bind
        radius_outer[newbins+1] = radius_central[newbins]+0.5*bin_width
        mass_shell = 4*pi/3*(radius_outer[newbins+1]**3-radius_outer[newbins]**3)*binned_density_tot[newbins]
        mass_inside[newbins+1] = mass_inside[newbins] + mass_shell
        mass_actual[newbins,nfile] = mass_inside[newbins+1]
        
        if rho_gas[inner_limit+bind+newbins*bind]==0.:
            continue
        
        if rho_gas[inner_limit+newbins*bind]==0.:
            continue
        
        if T_gas[inner_limit+bind+newbins*bind]==0.:
            continue
        
        if T_gas[inner_limit+newbins*bind]==0.:
            continue
        
        if sigma_r[inner_limit+bind+newbins*bind]==0.:
            continue
        
        if sigma_r[inner_limit+newbins*bind]==0.:
            continue
        
        
        binned_density_gas[newbins] = median(rho_gas[(inner_limit+newbins*bind):(inner_limit+bind+newbins*bind)])
        binned_T_gas[newbins] = median(T_gas[(inner_limit+newbins*bind):(inner_limit+bind+newbins*bind)])
        
        
        binned_velocitydispersion_r[newbins] = median(sigma_r[(inner_limit+newbins*bind):(inner_limit+bind+newbins*bind)])
        binned_velocitydispersion_t[newbins] = median(sigma_t[(inner_limit+newbins*bind):(inner_limit+bind+newbins*bind)])
        
        beta_actual[newbins,nfile] = 1 - (((binned_velocitydispersion_t[newbins])**2)/((binned_velocitydispersion_r[newbins])**2))
        
        hej = hej +1
    
    v_virial = sqrt(G*mass_inside[25]*M_sol/(radius_central[24]*kpc))
    r_virial = radius_central[24]
    
    variance_r_actual[:,nfile] = binned_velocitydispersion_r**2
    
    #print nfile+1
    #print hej
    
    # Fjern nu de elementer der er 0, før du fitter:
    
    #if T_gas[inner_limit+(n_bins-1)*bind] == 0.:
     #   print filelist_dm[nfile]
    
    bins_nonzero = [i for i in range(len(binned_T_gas)) if binned_T_gas[i] != 0]
    
    radius_central = radius_central[bins_nonzero]
    binned_density_gas = binned_density_gas[bins_nonzero]
    binned_T_gas = binned_T_gas[bins_nonzero]
    
    n_bins2 = len(binned_T_gas)
    
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
    
    dens_fit_result = minimize(chi2_dens_fit,params_dens,args=(radius_central[2:n_bins2],binned_density_gas[2:n_bins2]))
    
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
    
    gamma_gas = empty(n_bins2) # dlog_rho_dlog_r
    eta_gas = empty(n_bins2)   # dlog_T_dlog_r
    
    gamma_gas[1:n_bins2-1] = (log(density_gas_fit[1:n_bins2-1]/density_gas_fit[0:n_bins2-2])/log(radius_central[1:n_bins2-1]/radius_central[0:n_bins2-2])+log(density_gas_fit[2:n_bins2]/density_gas_fit[1:n_bins2-1])/log(radius_central[2:n_bins2]/radius_central[1:n_bins2-1]))/2
    gamma_gas[0] = (3*log(density_gas_fit[1]/density_gas_fit[0])/log(radius_central[1]/radius_central[0])+log(density_gas_fit[2]/density_gas_fit[1])/log(radius_central[2]/radius_central[1]))/4
    gamma_gas[n_bins2-1] = (3*log(density_gas_fit[n_bins2-1]/density_gas_fit[n_bins2-2])/log(radius_central[n_bins2-1]/radius_central[n_bins2-2])+log(density_gas_fit[n_bins2-2]/density_gas_fit[n_bins2-3])/log(radius_central[n_bins2-2]/radius_central[n_bins2-3]))/4
    
    gamma_gas_wolfram = -(params_dens['a'].value*params_dens['b'].value+radius_central*(params_dens['b'].value+params_dens['c'].value))/(params_dens['a'].value+radius_central)
    
    eta_gas[1:n_bins2-1] = (log(T_gas_fit[1:n_bins2-1]/T_gas_fit[0:n_bins2-2])/log(radius_central[1:n_bins2-1]/radius_central[0:n_bins2-2])+log(T_gas_fit[2:n_bins2]/T_gas_fit[1:n_bins2-1])/log(radius_central[2:n_bins2]/radius_central[1:n_bins2-1]))/2
    eta_gas[0] = (3*log(T_gas_fit[1]/T_gas_fit[0])/log(radius_central[1]/radius_central[0])+log(T_gas_fit[2]/T_gas_fit[1])/log(radius_central[2]/radius_central[1]))/4
    eta_gas[n_bins2-1] = (3*log(T_gas_fit[n_bins2-1]/T_gas_fit[n_bins2-2])/log(radius_central[n_bins2-1]/radius_central[n_bins2-2])+log(T_gas_fit[n_bins2-2]/T_gas_fit[n_bins2-3])/log(radius_central[n_bins2-2]/radius_central[n_bins2-3]))/4
    
    eta_gas_wolfram = -(params_temp['a'].value*params_temp['b'].value+radius_central*(params_temp['b'].value+params_temp['c'].value))/(params_temp['a'].value+radius_central)
    
    
    # Udregn massen ved GIDE:
    
    chosen_parameters = gas_parameters
    
    mass_term = -T_gas_fit*k_b/(m_p)*(gamma_gas_wolfram + eta_gas_wolfram) - S_term(radius_central,chosen_parameters[0],chosen_parameters[1],chosen_parameters[2],chosen_parameters[3],v_virial,r_virial)
    mass_calculated = mass_term*radius_central*kpc/(G*M_sol)
    
    def chi2_mass(params,r,data):
        
        a = params['a'].value
        b = params['b'].value
        c = params['c'].value
        #d = params['d'].value
        
        model = mass_fit_4(r,a,b,c)
        
        return (data-model)
    
    params = Parameters()
    
    params.add('a',value = 1.,min=0.)
    params.add('b',value = 1,min=-1.)
    params.add('c',value = 1.,min=0.)
    #params.add('d',value = 1.)#,max=0.)
    
    result = minimize(chi2_mass,params,args=(radius_central,mass_calculated))
    
    #report_fit(params)
    
    mass_fit_result = mass_fit_4(radius_central,params['a'].value,params['b'].value,params['c'].value)
    mass_fit_result = 0.9*mass_fit_result + 0.1*mass_actual[bins_nonzero,nfile]
    mass_fit_result = mass_actual[bins_nonzero,nfile]
    
    density_tot_calculated = empty(n_bins2)
    density_tot_calculated[1:n_bins2-1] = ((mass_fit_result[1:n_bins2-1]-mass_fit_result[0:n_bins2-2])/(radius_central[1:n_bins2-1]-radius_central[0:n_bins2-2])+(mass_fit_result[2:n_bins2]-mass_fit_result[1:n_bins2-1])/(radius_central[2:n_bins2]-radius_central[1:n_bins2-1]))/(2*4*pi*(radius_central[1:n_bins2-1])**2)
    density_tot_calculated[0] = (3*(mass_fit_result[1]-mass_fit_result[0])/(radius_central[1]-radius_central[0])+(mass_fit_result[2]-mass_fit_result[1])/(radius_central[2]-radius_central[1]))/(4*4*pi*(radius_central[0])**2)
    density_tot_calculated[n_bins2-1] = (3*(mass_fit_result[n_bins2-1]-mass_fit_result[n_bins2-2])/(radius_central[n_bins2-1]-radius_central[n_bins2-2])+(mass_fit_result[n_bins2-2]-mass_fit_result[n_bins2-3])/(radius_central[n_bins2-2]-radius_central[n_bins2-3]))/(4*4*pi*(radius_central[n_bins2-1])**2)
    density_dm_calculated = density_tot_calculated - density_gas_fit
    density_dm = density_dm_calculated
    
    #massfigure = figure(2+nfile)
    #plot(radius_central,mass_actual[bins_nonzero,nfile],'b',radius_central,mass_calculated,'g',radius_central,mass_fit_result,'r')
    #grid()
    #title(filelist_dm[nfile])
    #massfigure.show()
    
    gamma_dm = empty(n_bins2)
    
    gamma_dm[1:n_bins2-1] = (log(density_dm[1:n_bins2-1]/density_dm[0:n_bins2-2])/log(radius_central[1:n_bins2-1]/radius_central[0:n_bins2-2])+log(density_dm[2:n_bins2]/density_dm[1:n_bins2-1])/log(radius_central[2:n_bins2]/radius_central[1:n_bins2-1]))/2
    gamma_dm[0] = (3*log(density_dm[1]/density_dm[0])/log(radius_central[1]/radius_central[0])+log(density_dm[2]/density_dm[1])/log(radius_central[2]/radius_central[1]))/4
    gamma_dm[n_bins2-1] = (3*log(density_dm[n_bins2-1]/density_dm[n_bins2-2])/log(radius_central[n_bins2-1]/radius_central[n_bins2-2])+log(density_dm[n_bins2-2]/density_dm[n_bins2-3])/log(radius_central[n_bins2-2]/radius_central[n_bins2-3]))/4
    
    
    # Find Hastighedsdispersionen ved at integrere:
    
    T_dm = T_gas_fit * kappafunc2(radius_central/r_virial,kappa_func_parameters2[0],kappa_func_parameters2[1],kappa_func_parameters2[2],kappa_func_parameters2[3])
    #T_dm = m_p/(3*k_b) * (binned_velocitydispersion_r[bins_nonzero]**2+2*binned_velocitydispersion_t[bins_nonzero]**2)
    
    beta_inderst = -0.7
    variance_r_inderst = (T_dm[0]) * 3*k_b/(mean_molecular_weight*m_p * (3 - 2 * beta_inderst))
    
    
    psi_integrand = 3*k_b/(m_p*mean_molecular_weight)*T_dm -mass_fit_result*M_sol*G/(radius_central*kpc) - S_term(radius_central,chosen_parameters[0],chosen_parameters[1],chosen_parameters[2],chosen_parameters[3],v_virial,r_virial)
    integrand = psi_integrand*density_dm*radius_central**2
    
    integral_piece = empty(n_bins2)
    dispersion_integral = empty(n_bins2)
    
    integral_piece[0] = 0.1*integrand[0]*radius_central[0]
    dispersion_integral[0] = integral_piece[0] + variance_r_inderst*density_dm[0]*radius_central[0]**3
    
    
    for n_integrate in range(1,n_bins2):
        
        integral_piece[n_integrate] = 0.5*(integrand[n_integrate]+integrand[n_integrate-1])*(radius_central[n_integrate]-radius_central[n_integrate-1])
        dispersion_integral[n_integrate] = dispersion_integral[n_integrate-1] + integral_piece[n_integrate]
    
    radial_velocity_variance = 1./(radius_central**3*density_dm)*dispersion_integral
    
    #variancefigure = figure(1+nfile)
    #plot(radius_central,variance_r_actual[bins_nonzero,nfile],'b',radius_central,variance_r_calculated[bins_nonzero,nfile],'g')#,radius_central,variance_fit_result,'r')
    #grid()
    #title(filelist_dm[nfile])
    #variancefigure.show()
    
    # Find nu beta på to måder og sammenlign med den fra simulationen:
    
    #Først ved hjælp af definitionen af delta og beta:
    
    beta_definition = 1.5 * (1. - k_b/(mean_molecular_weight*m_p)*T_dm/radial_velocity_variance)
    
    # Så ved hjælp af Jean's ligning:
    
    eta_dm = empty(n_bins2)
    
    eta_dm[1:n_bins2-1] = (log(radial_velocity_variance[1:n_bins2-1]/radial_velocity_variance[0:n_bins2-2])/log(radius_central[1:n_bins2-1]/radius_central[0:n_bins2-2])+log(radial_velocity_variance[2:n_bins2]/radial_velocity_variance[1:n_bins2-1])/log(radius_central[2:n_bins2]/radius_central[1:n_bins2-1]))/2
    eta_dm[0] = (3*log(radial_velocity_variance[1]/radial_velocity_variance[0])/log(radius_central[1]/radius_central[0])+log(radial_velocity_variance[2]/radial_velocity_variance[1])/log(radius_central[2]/radius_central[1]))/4
    eta_dm[n_bins2-1] = (3*log(radial_velocity_variance[n_bins2-1]/radial_velocity_variance[n_bins2-2])/log(radius_central[n_bins2-1]/radius_central[n_bins2-2])+log(radial_velocity_variance[n_bins2-2]/radial_velocity_variance[n_bins2-3])/log(radius_central[n_bins2-2]/radius_central[n_bins2-3]))/4
    
    beta_jeans = (-G*mass_fit_result*M_sol/(radius_central*kpc*radial_velocity_variance) - eta_dm - gamma_dm)/2.
    
    
    beta_def_matrix_less[bins_nonzero,nfile] = beta_definition
    beta_jeans_matrix_less[bins_nonzero,nfile] = beta_jeans
    
    # En gang til med anden indre beta:
    
    beta_inderst = 0.
    variance_r_inderst = (T_dm[0]) * 3*k_b/(mean_molecular_weight*m_p * (3 - 2 * beta_inderst))
    
    
    integral_piece[0] = 0.1*integrand[0]*radius_central[0]
    dispersion_integral[0] = integral_piece[0] + variance_r_inderst*density_dm[0]*radius_central[0]**3
    
    
    for n_integrate in range(1,n_bins2):
        
        integral_piece[n_integrate] = 0.5*(integrand[n_integrate]+integrand[n_integrate-1])*(radius_central[n_integrate]-radius_central[n_integrate-1])
        dispersion_integral[n_integrate] = dispersion_integral[n_integrate-1] + integral_piece[n_integrate]
    
    radial_velocity_variance = 1./(radius_central**3*density_dm)*dispersion_integral
    
    #radial_velocity_variance = variance_r_actual[bins_nonzero,nfile]
    
    
    
    # Find nu beta på to måder og sammenlign med den fra simulationen:
    
    #Først ved hjælp af definitionen af delta og beta:
    
    beta_definition = 1.5 * (1. - k_b/(mean_molecular_weight*m_p)*T_dm/radial_velocity_variance)
    
    # Så ved hjælp af Jean's ligning:
    
    eta_dm = empty(n_bins2)
    
    eta_dm[1:n_bins2-1] = (log(radial_velocity_variance[1:n_bins2-1]/radial_velocity_variance[0:n_bins2-2])/log(radius_central[1:n_bins2-1]/radius_central[0:n_bins2-2])+log(radial_velocity_variance[2:n_bins2]/radial_velocity_variance[1:n_bins2-1])/log(radius_central[2:n_bins2]/radius_central[1:n_bins2-1]))/2
    eta_dm[0] = (3*log(radial_velocity_variance[1]/radial_velocity_variance[0])/log(radius_central[1]/radius_central[0])+log(radial_velocity_variance[2]/radial_velocity_variance[1])/log(radius_central[2]/radius_central[1]))/4
    eta_dm[n_bins2-1] = (3*log(radial_velocity_variance[n_bins2-1]/radial_velocity_variance[n_bins2-2])/log(radius_central[n_bins2-1]/radius_central[n_bins2-2])+log(radial_velocity_variance[n_bins2-2]/radial_velocity_variance[n_bins2-3])/log(radius_central[n_bins2-2]/radius_central[n_bins2-3]))/4
    
    beta_jeans = (-G*mass_fit_result*M_sol/(radius_central*kpc*radial_velocity_variance) - eta_dm - gamma_dm)/2.
    
    
    beta_def_matrix[bins_nonzero,nfile] = beta_definition
    beta_jeans_matrix[bins_nonzero,nfile] = beta_jeans
    
    
    # En tredje gang med tredje indre beta:
    
    beta_inderst = 0.5
    variance_r_inderst = (T_dm[0]) * 3*k_b/(mean_molecular_weight*m_p * (3 - 2 * beta_inderst))
    
    
    integral_piece[0] = 0.1*integrand[0]*radius_central[0]
    dispersion_integral[0] = integral_piece[0] + variance_r_inderst*density_dm[0]*radius_central[0]**3
    
    
    for n_integrate in range(1,n_bins2):
        
        integral_piece[n_integrate] = 0.5*(integrand[n_integrate]+integrand[n_integrate-1])*(radius_central[n_integrate]-radius_central[n_integrate-1])
        dispersion_integral[n_integrate] = dispersion_integral[n_integrate-1] + integral_piece[n_integrate]
    
    radial_velocity_variance = 1./(radius_central**3*density_dm)*dispersion_integral
    
    # Find nu beta på to måder og sammenlign med den fra simulationen:
    
    #Først ved hjælp af definitionen af delta og beta:
    
    beta_definition = 1.5 * (1. - k_b/(mean_molecular_weight*m_p)*T_dm/radial_velocity_variance)
    
    # Så ved hjælp af Jean's ligning:
    
    eta_dm = empty(n_bins2)
    
    eta_dm[1:n_bins2-1] = (log(radial_velocity_variance[1:n_bins2-1]/radial_velocity_variance[0:n_bins2-2])/log(radius_central[1:n_bins2-1]/radius_central[0:n_bins2-2])+log(radial_velocity_variance[2:n_bins2]/radial_velocity_variance[1:n_bins2-1])/log(radius_central[2:n_bins2]/radius_central[1:n_bins2-1]))/2
    eta_dm[0] = (3*log(radial_velocity_variance[1]/radial_velocity_variance[0])/log(radius_central[1]/radius_central[0])+log(radial_velocity_variance[2]/radial_velocity_variance[1])/log(radius_central[2]/radius_central[1]))/4
    eta_dm[n_bins2-1] = (3*log(radial_velocity_variance[n_bins2-1]/radial_velocity_variance[n_bins2-2])/log(radius_central[n_bins2-1]/radius_central[n_bins2-2])+log(radial_velocity_variance[n_bins2-2]/radial_velocity_variance[n_bins2-3])/log(radius_central[n_bins2-2]/radius_central[n_bins2-3]))/4
    
    beta_jeans = (-G*mass_fit_result*M_sol/(radius_central*kpc*radial_velocity_variance) - eta_dm - gamma_dm)/2.
    
    beta_def_matrix_more[bins_nonzero,nfile] = beta_definition
    beta_jeans_matrix_more[bins_nonzero,nfile] = beta_jeans
    
    print nfile+1
    

radius_relative1 = arange((0.0+inner_limit*2.0/len(rho_dm)),2.0-2.0/n_bins1,(2.0/n_bins1))+(1.0/n_bins1)

radius_relative = radius_relative1[0:n_bins]



# Statistik:

median_beta_actual = empty(n_bins)
dispersion_beta_actual_upper = empty(n_bins)
dispersion_beta_actual_lower = empty(n_bins)

median_beta_def = empty(n_bins)
dispersion_beta_def_upper = empty(n_bins)
dispersion_beta_def_lower = empty(n_bins)

median_beta_jeans = empty(n_bins)
dispersion_beta_jeans_upper = empty(n_bins)
dispersion_beta_jeans_lower = empty(n_bins)

median_beta_less_def = empty(n_bins)
dispersion_beta_less_def_upper = empty(n_bins)
dispersion_beta_less_def_lower = empty(n_bins)

median_beta_less_jeans = empty(n_bins)
dispersion_beta_less_jeans_upper = empty(n_bins)
dispersion_beta_less_jeans_lower = empty(n_bins)

median_beta_more_def = empty(n_bins)
dispersion_beta_more_def_upper = empty(n_bins)
dispersion_beta_more_def_lower = empty(n_bins)

median_beta_more_jeans = empty(n_bins)
dispersion_beta_more_jeans_upper = empty(n_bins)
dispersion_beta_more_jeans_lower = empty(n_bins)

for nbins3 in range(0,n_bins):
    
    beta_actual_sorted = sorted(beta_actual[nbins3,:])
    beta_actual_sorted = [i for i in beta_actual_sorted if i <> 0.]
    median_beta_actual[nbins3] = beta_actual_sorted[int(len(beta_actual_sorted)*0.5)]
    sigma_file = int(len(beta_actual_sorted)*0.341)
    dispersion_beta_actual_upper[nbins3] = beta_actual_sorted[int(len(beta_actual_sorted)*0.5)+sigma_file] - median_beta_actual[nbins3]
    dispersion_beta_actual_lower[nbins3] = median_beta_actual[nbins3] - beta_actual_sorted[int(len(beta_actual_sorted)*0.5)-sigma_file]
    
    beta_def_sorted = sorted(beta_def_matrix[nbins3,:])
    beta_def_sorted = [i for i in beta_def_sorted if i <> 0.]
    median_beta_def[nbins3] = beta_def_sorted[int(len(beta_def_sorted)*0.5)]
    dispersion_beta_def_upper[nbins3] = beta_def_sorted[int(len(beta_def_sorted)*0.5)+sigma_file] - median_beta_def[nbins3]
    dispersion_beta_def_lower[nbins3] = median_beta_def[nbins3] - beta_def_sorted[int(len(beta_def_sorted)*0.5)-sigma_file]
    
    beta_jeans_sorted = sorted(beta_jeans_matrix[nbins3,:])
    beta_jeans_sorted = [i for i in beta_jeans_sorted if i <> 0.]
    median_beta_jeans[nbins3] = beta_jeans_sorted[int(len(beta_jeans_sorted)*0.5)]
    dispersion_beta_jeans_upper[nbins3] = beta_jeans_sorted[int(len(beta_jeans_sorted)*0.5)+sigma_file] - median_beta_jeans[nbins3]
    dispersion_beta_jeans_lower[nbins3] = median_beta_jeans[nbins3] - beta_jeans_sorted[int(len(beta_jeans_sorted)*0.5)-sigma_file]
    
    beta_less_def_sorted = sorted(beta_def_matrix_less[nbins3,:])
    beta_less_def_sorted = [i for i in beta_less_def_sorted if i <> 0.]
    median_beta_less_def[nbins3] = beta_less_def_sorted[int(len(beta_less_def_sorted)*0.5)]
    sigma_file1 = int(len(beta_less_def_sorted)*0.341)
    dispersion_beta_less_def_upper[nbins3] = beta_less_def_sorted[int(len(beta_less_def_sorted)*0.5)+sigma_file1] - median_beta_less_def[nbins3]
    dispersion_beta_less_def_lower[nbins3] = median_beta_less_def[nbins3] - beta_less_def_sorted[int(len(beta_less_def_sorted)*0.5)-sigma_file1]
    
    beta_less_jeans_sorted = sorted(beta_jeans_matrix_less[nbins3,:])
    beta_less_jeans_sorted = [i for i in beta_less_jeans_sorted if i <> 0.]
    median_beta_less_jeans[nbins3] = beta_less_jeans_sorted[int(len(beta_less_jeans_sorted)*0.5)]
    sigma_file1 = int(len(beta_less_jeans_sorted)*0.341)
    dispersion_beta_less_jeans_upper[nbins3] = beta_less_jeans_sorted[int(len(beta_less_jeans_sorted)*0.5)+sigma_file1] - median_beta_less_jeans[nbins3]
    dispersion_beta_less_jeans_lower[nbins3] = median_beta_less_jeans[nbins3] - beta_less_jeans_sorted[int(len(beta_less_jeans_sorted)*0.5)-sigma_file1]
    
    beta_more_def_sorted = sorted(beta_def_matrix_more[nbins3,:])
    beta_more_def_sorted = [i for i in beta_more_def_sorted if i <> 0.]
    median_beta_more_def[nbins3] = beta_more_def_sorted[int(len(beta_more_def_sorted)*0.5)]
    #print beta_more_def_sorted
    sigma_file1 = int(len(beta_more_def_sorted)*0.341)
    dispersion_beta_more_def_upper[nbins3] = beta_more_def_sorted[int(len(beta_more_def_sorted)*0.5)+sigma_file1] - median_beta_more_def[nbins3]
    dispersion_beta_more_def_lower[nbins3] = median_beta_more_def[nbins3] - beta_more_def_sorted[int(len(beta_more_def_sorted)*0.5)-sigma_file1]
    
    beta_more_jeans_sorted = sorted(beta_jeans_matrix_more[nbins3,:])
    beta_more_jeans_sorted = [i for i in beta_more_jeans_sorted if i <> 0.]
    median_beta_more_jeans[nbins3] = beta_more_jeans_sorted[int(len(beta_more_jeans_sorted)*0.5)]
    sigma_file1 = int(len(beta_more_jeans_sorted)*0.341)
    dispersion_beta_more_jeans_upper[nbins3] = beta_more_jeans_sorted[int(len(beta_more_jeans_sorted)*0.5)+sigma_file1] - median_beta_more_jeans[nbins3]
    dispersion_beta_more_jeans_lower[nbins3] = median_beta_more_jeans[nbins3] - beta_more_jeans_sorted[int(len(beta_more_jeans_sorted)*0.5)-sigma_file1]
    



# Plotning:

betaplot = figure(42)
axdisp2 = betaplot.add_subplot(111)
axdisp2.errorbar(radius_relative,median_beta_actual,yerr=[dispersion_beta_actual_lower,dispersion_beta_actual_upper],fmt='^b',markersize = 7.,label = r'$\beta_{TRUE}$ Gadget')
axdisp2.errorbar(radius_relative-0.005,median_beta_def,yerr=[dispersion_beta_def_lower,dispersion_beta_def_upper],fmt='sr', label = r'$\beta_{DEF}$ Gadget, $\beta (r=0) = 0.0$, $M_{TRUE}$')
axdisp2.errorbar(radius_relative+0.005,median_beta_jeans,yerr=[dispersion_beta_jeans_lower,dispersion_beta_jeans_upper],fmt='pg',markersize = 7., label = r'$\beta_{GJE}$ Gadget, $\beta (r=0) = 0.0$, $M_{TRUE}$')
grid()

legend(loc='upper left',numpoints=1,fontsize=24)
xlabel('$r/r_{200}$', fontsize=22)
ylabel(r' $ \beta $', fontsize=22)

betaplot.show()



defplot = figure(73)
axdisp2 = defplot.add_subplot(111)
axdisp2.errorbar(radius_relative,median_beta_actual,yerr=[dispersion_beta_actual_lower,dispersion_beta_actual_upper],fmt='^b',markersize = 7.,label = r'$\beta_{TRUE}$ Gadget')
axdisp2.errorbar(radius_relative-0.01,median_beta_less_def,yerr=[dispersion_beta_less_def_lower,dispersion_beta_less_def_upper],fmt='sm', label = r'$\beta_{DEF}$ Gadget, $\beta (r=0) = -0.7$, $M_{TRUE}$')
axdisp2.errorbar(radius_relative-0.005,median_beta_def,yerr=[dispersion_beta_def_lower,dispersion_beta_def_upper],fmt='sr', label = r'$\beta_{DEF}$ Gadget, $\beta (r=0) = 0.1$, $M_{TRUE}$')
axdisp2.errorbar(radius_relative+0.005,median_beta_more_def,yerr=[dispersion_beta_more_def_lower,dispersion_beta_more_def_upper],fmt='sg', label = r'$\beta_{DEF}$ Gadget, $\beta (r=0) = 0.5$, $M_{TRUE}$')
grid()
ylim(-0.3,1.4)
legend(loc='upper left',numpoints=1,fontsize=24)
xlabel('$r/r_{200}$', fontsize=22)
ylabel(r' $ \beta $', fontsize=22)

defplot.show()



jeansplot = figure(74)
axdisp2 = jeansplot.add_subplot(111)
axdisp2.errorbar(radius_relative,median_beta_actual,yerr=[dispersion_beta_actual_lower,dispersion_beta_actual_upper],fmt='^b',markersize = 7., label = r'$\beta_{TRUE}$ Gadget')
axdisp2.errorbar(radius_relative-0.01,median_beta_less_jeans,yerr=[dispersion_beta_less_jeans_lower,dispersion_beta_less_jeans_upper],fmt='pm',markersize = 7., label = r'$\beta_{GJE}$ Gadget, $\beta (r=0) = -0.7$, $M_{TRUE}$')
axdisp2.errorbar(radius_relative-0.005,median_beta_jeans,yerr=[dispersion_beta_jeans_lower,dispersion_beta_jeans_upper],fmt='pr',markersize = 7., label = r'$\beta_{GJE}$ Gadget, $\beta (r=0) = 0.1$, $M_{TRUE}$')
axdisp2.errorbar(radius_relative+0.005,median_beta_more_jeans,yerr=[dispersion_beta_more_jeans_lower,dispersion_beta_more_jeans_upper],fmt='pg',markersize = 7., label = r'$\beta_{GJE}$ Gadget, $\beta (r=0) = 0.5$, $M_{TRUE}$')
grid()
ylim(-0.1,1.3)
legend(loc='upper left',numpoints=1,fontsize=24)
xlabel('$r/r_{200}$', fontsize=22)
ylabel(r' $ \beta $', fontsize=22)

jeansplot.show()

