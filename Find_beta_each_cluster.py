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


# Load all the data:


n_bins = 19  # Final number of bins per profile
n_bins1 = n_bins+1
n_profiles = 51  # Total number of profiles

binned_velocitydispersion_r = empty(shape=(n_bins,n_profiles))
mass_simulation = empty(shape=(n_bins,n_profiles))

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
    
    testfigure = figure(nfile+1)
    plot((radius),log10(rho_gas))
    grid()
    testfigure.show()
    
    #binned_velocitydispersion_r = empty(n_bins)
    binned_velocitydispersion_theta = empty(n_bins)
    binned_velocitydispersion_phi = empty(n_bins)
    
    bind = len(rho_dm)/n_bins1
    bin_width = (radius[5]-radius[4])
    
    
    for newbins in range(0,n_bins):
        
        radius_central[newbins] = sum(radius[(1+newbins*bind):(1+bind+newbins*bind)])/bind
        binned_density_dm[newbins] = sum(rho_dm[(1+newbins*bind):(1+bind+newbins*bind)])/bind
        
        binned_velocitydispersion_r[newbins,nfile] = sum(sigma_r[(1+newbins*bind):(1+bind+newbins*bind)])/bind
        binned_velocitydispersion_theta[newbins] = sum(sigma_theta[(1+newbins*bind):(1+bind+newbins*bind)])/bind
        binned_velocitydispersion_phi[newbins] = sum(sigma_phi[(1+newbins*bind):(1+bind+newbins*bind)])/bind
        
        binned_density_gas[newbins] = sum(rho_gas[(1+newbins*bind):(1+bind+newbins*bind)])/bind
        binned_T_gas[newbins] = sum(T_gas[(1+newbins*bind):(1+bind+newbins*bind)])/bind
        
        beta_simulation[newbins,nfile] = 1 - (((binned_velocitydispersion_theta[newbins])**2+(binned_velocitydispersion_phi[newbins])**2)/(2*(binned_velocitydispersion_r[newbins,nfile])**2))
        
        binned_density_tot[newbins] = sum(rho_tot[(1+newbins*bind):(1+bind+newbins*bind)])/bind
        radius_outer[newbins+1] = radius[(newbins+1)*bind]+0.5*bin_width
        mass_shell = 4*pi/3*(radius_outer[newbins+1]**3-radius_outer[newbins]**3)*binned_density_tot[newbins]
        mass_inside = mass_inside + mass_shell
        
        mass_simulation[newbins,nfile] = mass_inside
    
    r_virial = radius[int(len(radius)/2.)]
    
    # Udregn de logarimitsk afledte for gassen:
    
    gamma_gas = empty(n_bins) # dlog_rho_dlog_r
    eta_gas = empty(n_bins)   # dlog_T_dlog_r
    
    gamma_gas[1:n_bins-1] = ((binned_density_gas[1:n_bins-1]-binned_density_gas[0:n_bins-2])/(radius_central[1:n_bins-1]-radius_central[0:n_bins-2])+(binned_density_gas[2:n_bins]-binned_density_gas[1:n_bins-1])/(radius_central[2:n_bins]-radius_central[1:n_bins-1]))/2
    gamma_gas[0] = (3*(binned_density_gas[1]-binned_density_gas[0])/(radius_central[1]-radius_central[0])+(binned_density_gas[2]-binned_density_gas[1])/(radius_central[2]-radius_central[1]))/4
    gamma_gas[n_bins-1] = (3*(binned_density_gas[n_bins-1]-binned_density_gas[n_bins-2])/(radius_central[n_bins-1]-radius_central[n_bins-2])+(binned_density_gas[n_bins-2]-binned_density_gas[n_bins-3])/(radius_central[n_bins-2]-radius_central[n_bins-3]))/4
    gamma_gas = gamma_gas*radius_central/binned_density_gas
    
    eta_gas[1:n_bins-1] = ((binned_T_gas[1:n_bins-1]-binned_T_gas[0:n_bins-2])/(radius_central[1:n_bins-1]-radius_central[0:n_bins-2])+(binned_T_gas[2:n_bins]-binned_T_gas[1:n_bins-1])/(radius_central[2:n_bins]-radius_central[1:n_bins-1]))/2
    eta_gas[0] = (3*(binned_T_gas[1]-binned_T_gas[0])/(radius_central[1]-radius_central[0])+(binned_T_gas[2]-binned_T_gas[1])/(radius_central[2]-radius_central[1]))/4
    eta_gas[n_bins-1] = (3*(binned_T_gas[n_bins-1]-binned_T_gas[n_bins-2])/(radius_central[n_bins-1]-radius_central[n_bins-2])+(binned_T_gas[n_bins-2]-binned_T_gas[n_bins-3])/(radius_central[n_bins-2]-radius_central[n_bins-3]))/4
    eta_gas = eta_gas*radius_central/binned_T_gas
    
    
    # Find Hastighedsdispersionen ved at integrere:
    
    T_dm = binned_T_gas + difffunc(radius_central/r_virial,diff_parameters[0],diff_parameters[1],diff_parameters[2],diff_parameters[3])
    
    psi_integrand = binned_T_gas*k_b/(m_p*mean_molecular_weight)*(gamma_gas + eta_gas) + 3*k_b/(m_p*mean_molecular_weight)*T_dm
    integrand = 10**-10*psi_integrand*binned_density_dm*radius_central**2
    
    #print gamma_gas
    
    integral_piece = empty(n_bins)
    dispersion_integral = empty(n_bins)
    radial_velocity_variance = empty(n_bins)
    
    integral_piece[0] = 0.5*integrand[0]*radius_central[0]
    dispersion_integral[0] = integral_piece[0]
    radial_velocity_variance[0] = 10**10/(radius_central[0]**3*binned_density_dm[0])*dispersion_integral[0]
    
    for n_integrate in range(1,n_bins):
        
        integral_piece[n_integrate] = 0.5*(integrand[n_integrate]+integrand[n_integrate-1])*(radius_central[n_integrate]-radius_central[n_integrate-1])
        dispersion_integral[n_integrate] = dispersion_integral[n_integrate-1] + integral_piece[n_integrate]
        radial_velocity_variance[n_integrate] = 10**10/(radius_central[n_integrate]**3*binned_density_dm[n_integrate])*dispersion_integral[n_integrate]
    
    
    beta_inderst = 0.
    sigma_r_inderst = (T_dm[0]) * 3*k_b/(mean_molecular_weight*m_p * (3 - 2 * beta_inderst))
    radial_velocity_variance = radial_velocity_variance - radial_velocity_variance[0] + sigma_r_inderst #+ 1.5*10**11#- radial_velocity_variance[0] + median_sigma_r_simulation[0]#+1.7*10**11#
    variance_r_calculated[:,nfile] = radial_velocity_variance
    
    variance_r_actual[:,nfile] = binned_velocitydispersion_r[:,nfile]**2
    #print radial_velocity_variance
    
    
    # Find nu beta på to måder og sammenlign med den fra simulationen:
    
    #Først ved hjælp af definitionen af delta og beta:
    
    beta_definition = 1.5 * (1. - k_b/(mean_molecular_weight*m_p)*T_dm/radial_velocity_variance)
    
    # Så ved hjælp af Jean's ligning:
    
    gamma_dm = empty(n_bins)
    eta_dm = empty(n_bins)
    
    gamma_dm[1:n_bins-1] = ((binned_density_dm[1:n_bins-1]-binned_density_dm[0:n_bins-2])/(radius_central[1:n_bins-1]-radius_central[0:n_bins-2])+(binned_density_dm[2:n_bins]-binned_density_dm[1:n_bins-1])/(radius_central[2:n_bins]-radius_central[1:n_bins-1]))/2
    gamma_dm[0] = (3*(binned_density_dm[1]-binned_density_dm[0])/(radius_central[1]-radius_central[0])+(binned_density_dm[2]-binned_density_dm[1])/(radius_central[2]-radius_central[1]))/4
    gamma_dm[n_bins-1] = (3*(binned_density_dm[n_bins-1]-binned_density_dm[n_bins-2])/(radius_central[n_bins-1]-radius_central[n_bins-2])+(binned_density_dm[n_bins-2]-binned_density_dm[n_bins-3])/(radius_central[n_bins-2]-radius_central[n_bins-3]))/4
    gamma_dm = gamma_dm*radius_central/binned_density_dm
    
    eta_dm[1:n_bins-1] = ((radial_velocity_variance[1:n_bins-1]-radial_velocity_variance[0:n_bins-2])/(radius_central[1:n_bins-1]-radius_central[0:n_bins-2])+(radial_velocity_variance[2:n_bins]-radial_velocity_variance[1:n_bins-1])/(radius_central[2:n_bins]-radius_central[1:n_bins-1]))/2
    eta_dm[0] = (3*(radial_velocity_variance[1]-radial_velocity_variance[0])/(radius_central[1]-radius_central[0])+(radial_velocity_variance[2]-radial_velocity_variance[1])/(radius_central[2]-radius_central[1]))/4
    eta_dm[n_bins-1] = (3*(radial_velocity_variance[n_bins-1]-radial_velocity_variance[n_bins-2])/(radius_central[n_bins-1]-radius_central[n_bins-2])+(radial_velocity_variance[n_bins-2]-radial_velocity_variance[n_bins-3])/(radius_central[n_bins-2]-radius_central[n_bins-3]))/4
    eta_dm = eta_dm*radius_central/radial_velocity_variance
    
    beta_jeans = (-G*mass_simulation[:,nfile]*M_sol/(radius_central*kpc*radial_velocity_variance) - eta_dm - gamma_dm)/2.
    
    #print beta_definition
    
    beta_def_matrix[:,nfile] = beta_definition
    beta_jeans_matrix[:,nfile] = beta_jeans
    
    
    

radius_relative = arange((0.0+2.0/len(rho_dm)),2.0-2.0/n_bins1,(2.0/n_bins1))+(1.0/n_bins1)


median_variance_actual = empty(n_bins)
dispersion_variance_actual_upper = empty(n_bins)
dispersion_variance_actual_lower = empty(n_bins)

median_variance_calculated = empty(n_bins)
dispersion_variance_calculated_upper = empty(n_bins)
dispersion_variance_calculated_lower = empty(n_bins)

for n5 in range(0,n_bins):
    
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
 

varianceplot = figure(94)
axdisp2 = varianceplot.add_subplot(111)
#axdisp2.set_xscale("log", nonposx='clip')
#axdisp2.set_yscale("log", nonposy='clip')
axdisp2.errorbar(radius_relative,median_variance_actual,yerr=[dispersion_variance_actual_lower,dispersion_variance_actual_upper],fmt='ob')
axdisp2.errorbar(radius_relative,median_variance_calculated,yerr=[dispersion_variance_calculated_lower,dispersion_variance_calculated_upper],fmt='or')
grid()
#varianceplot.show()


sys.exit(0)

# Statistisk behandling af hastighedsdispersionen:



# Lav et datasæt der svarer til eksperimentel måling af X-ray stråling fra en hob der svarer til alle simulationerne stacked.
# Du har temperatur og tæthed af gassen ved hver radius. 

for nbins3 in range(0,n_bins):
    
    radius_sorted = sorted(radius_central[nbins3,:])
    median_radius[nbins3] = radius_sorted[int(len(radius_sorted)*0.5)]
    sigma_file1 = int(len(radius_sorted)*0.341)
    dispersion_radius_upper[nbins3] = radius_sorted[int(len(radius_sorted)*0.5)+sigma_file1] - median_radius[nbins3]
    dispersion_radius_lower[nbins3] = median_radius[nbins3] - radius_sorted[int(len(radius_sorted)*0.5)-sigma_file1]
    
    density_gas_sorted = sorted(binned_density_gas[nbins3,:])
    median_density_gas[nbins3] = density_gas_sorted[int(len(density_gas_sorted)*0.5)]
    sigma_file2 = int(len(density_gas_sorted)*0.341)
    dispersion_density_gas_upper[nbins3] = density_gas_sorted[int(len(density_gas_sorted)*0.5)+sigma_file2] - median_density_gas[nbins3]
    dispersion_density_gas_lower[nbins3] = median_density_gas[nbins3] - density_gas_sorted[int(len(density_gas_sorted)*0.5)-sigma_file2]
    
    T_gas_sorted = sorted(binned_T_gas[nbins3,:])
    median_T_gas[nbins3] = T_gas_sorted[int(len(T_gas_sorted)*0.5)]
    sigma_file = int(len(T_gas_sorted)*0.341)
    dispersion_T_gas_upper[nbins3] = T_gas_sorted[int(len(T_gas_sorted)*0.5)+sigma_file] - median_T_gas[nbins3]
    dispersion_T_gas_lower[nbins3] = median_T_gas[nbins3] - T_gas_sorted[int(len(T_gas_sorted)*0.5)-sigma_file]
    
# Definer funktionen S (mere specifikt S*r):

def S_term(r,alpha,a,c,d,v_virial,r_virial):
    
    v_p = -alpha/(((r/r_virial)**(-a) + c*(r/r_virial)**(a/2))**(1/a)-d)*v_virial
    
    v_p_derivative = v_virial/r_virial*alpha/a*(a/2*c*(r/r_virial)**(a/2-1)-a*(r/r_virial)**(-a-1))*(((r/r_virial)**(-a)+c*(r/r_virial)**(a/2))**(1/a-1))/(((r/r_virial)**(-a)+c*(r/r_virial)**(a/2))**(1/a) - d)**2
    
    return v_p*v_p_derivative*r + v_p*H*r + v_p_derivative*r**2*H

dm_parameters = [0.11255906,12.8245615,2.2179*10**-6,0.31763796]
gas_parameters = [0.13669927,12.635266,6.3028*10**-6,0.2453202]

chosen_parameters = dm_parameters

t = 1.  #0.928  # Hvor stor en vægt S har (her igennem alpha)

# Definer "Kappa", som du har fra simulationerne:

def difffunc(r,q,p,o,a):
    return q - o/((r+1.3)/a)**p

diff_parameters = [4.0158e+06,16.8419116,6.0739e+05,1.67671123]

# Lav en model for masseprofilen som du kan bruge til at bestemme densiteten med:

def mass_fit_4(r,a,b,c,d):
    
    return 10**14*(c*((r/a)**d*(r/a+1)**-b)**(1/d))

# Lav Monte Carlo sampling:


density_dispersion = (dispersion_density_gas_upper+dispersion_density_gas_lower)/2.
T_dispersion = (dispersion_T_gas_upper+dispersion_T_gas_lower)/2.

sampling_size = 50 #################################################################################

mass_matrix = empty(shape=(n_bins,sampling_size))
beta_def_matrix = empty(shape=(n_bins,sampling_size))
beta_jeans_matrix = empty(shape=(n_bins,sampling_size))
sigma_r_matrix = empty(shape=(n_bins,sampling_size))
density_dm_matrix = empty(shape=(n_bins,sampling_size))


for n1000 in range(0,sampling_size):
    
    # Lav gaussisk fordelte arrays af temperatur og densitet
    
    for n1 in range(0,n_bins):
        density_gas_gauss[n1] = random.normal(median_density_gas_fit[n1],density_dispersion[n1])
        T_gas_gauss[n1] = random.normal(median_T_gas_fit[n1],T_dispersion[n1])
        
    # Udregn logaritmiske hældninger for T og rho:
    
    gamma_gas = empty(n_bins) # dlog_rho_dlog_r
    eta_gas = empty(n_bins)   # dlog_T_dlog_r
    
    gamma_gas[1:n_bins-1] = ((density_gas_gauss[1:n_bins-1]-density_gas_gauss[0:n_bins-2])/(median_radius[1:n_bins-1]-median_radius[0:n_bins-2])+(density_gas_gauss[2:n_bins]-density_gas_gauss[1:n_bins-1])/(median_radius[2:n_bins]-median_radius[1:n_bins-1]))/2
    gamma_gas[0] = (3*(density_gas_gauss[1]-density_gas_gauss[0])/(median_radius[1]-median_radius[0])+(density_gas_gauss[2]-density_gas_gauss[1])/(median_radius[2]-median_radius[1]))/4
    gamma_gas[n_bins-1] = (3*(density_gas_gauss[n_bins-1]-density_gas_gauss[n_bins-2])/(median_radius[n_bins-1]-median_radius[n_bins-2])+(density_gas_gauss[n_bins-2]-density_gas_gauss[n_bins-3])/(median_radius[n_bins-2]-median_radius[n_bins-3]))/4
    gamma_gas = gamma_gas*median_radius/density_gas_gauss
    
    eta_gas[1:n_bins-1] = ((T_gas_gauss[1:n_bins-1]-T_gas_gauss[0:n_bins-2])/(median_radius[1:n_bins-1]-median_radius[0:n_bins-2])+(T_gas_gauss[2:n_bins]-T_gas_gauss[1:n_bins-1])/(median_radius[2:n_bins]-median_radius[1:n_bins-1]))/2
    eta_gas[0] = (3*(T_gas_gauss[1]-T_gas_gauss[0])/(median_radius[1]-median_radius[0])+(T_gas_gauss[2]-T_gas_gauss[1])/(median_radius[2]-median_radius[1]))/4
    eta_gas[n_bins-1] = (3*(T_gas_gauss[n_bins-1]-T_gas_gauss[n_bins-2])/(median_radius[n_bins-1]-median_radius[n_bins-2])+(T_gas_gauss[n_bins-2]-T_gas_gauss[n_bins-3])/(median_radius[n_bins-2]-median_radius[n_bins-3]))/4
    eta_gas = eta_gas*median_radius/T_gas_gauss
    
    # Find masseleddet ud fra hydrostatisk ligevægt:
    
    #mass_term = -T_gas_gauss*k_b/(m_p*mean_molecular_weight)*(gamma_gas + eta_gas)
    #mass = mass_term*median_radius*kpc/(G*M_sol)
    #mass_matrix[:,n1000] = mass
    gamma_gas[1] = gamma_gas[1] + 4.5
    
    # Udregn virialradius ved at interpolere for at finde S, og iterer så:
    
    S_term_list = zeros(n_bins)
    
    iterations = 2
    
    for iteration in range(0,iterations):
        
        mass_term = -T_gas_gauss*k_b/(m_p*mean_molecular_weight)*(gamma_gas + eta_gas) - S_term_list*0.5
        mass = mass_term*median_radius*kpc/(G*M_sol)
        
        mass_interpol_func = interp1d(median_radius,mass)
        mass_interpol = mass_interpol_func(radius_interpol)
        density_average = mass_interpol / (4/3*pi*radius_interpol**3)
        
        find_virial = min(enumerate(density_average), key=lambda x: abs(x[1]-(200*rho_c))) # Finder den indgang hvor densiteten er tættest på 200 gange den kritiske
        
        r_virial = radius_interpol[find_virial[0]]
        mass_virial = mass_interpol[find_virial[0]]
        v_virial = sqrt(G*mass_virial*M_sol/(r_virial*kpc))
        
        S_term_list = S_term(median_radius,chosen_parameters[0],chosen_parameters[1],chosen_parameters[2],chosen_parameters[3],v_virial,r_virial)
        
        #print r_virial
    
    mass_matrix[:,n1000] = mass
    
    
    
    # Udregn den totale densitet og dermed densiteten af mørkt stof:
    
    density_tot = empty(n_bins)
    density_tot[1:n_bins-1] = ((mass[1:n_bins-1]-mass[0:n_bins-2])/(median_radius[1:n_bins-1]-median_radius[0:n_bins-2])+(mass[2:n_bins]-mass[1:n_bins-1])/(median_radius[2:n_bins]-median_radius[1:n_bins-1]))/(2*4*pi*(median_radius[1:n_bins-1])**2)
    density_tot[0] = (3*(mass[1]-mass[0])/(median_radius[1]-median_radius[0])+(mass[2]-mass[1])/(median_radius[2]-median_radius[1]))/(4*4*pi*(median_radius[0])**2)
    density_tot[n_bins-1] = (3*(mass[n_bins-1]-mass[n_bins-2])/(median_radius[n_bins-1]-median_radius[n_bins-2])+(mass[n_bins-2]-mass[n_bins-3])/(median_radius[n_bins-2]-median_radius[n_bins-3]))/(4*4*pi*(median_radius[n_bins-1])**2)
    density_dm = density_tot-density_gas_gauss
    
    density_dm_matrix[:,n1000] = density_dm
    # Udregn nu den radielle hastighedsdispersion for mørkt stof, idet du antager at kende temperatur-relationen:
    

    



    
# Lav statistisk behandling af massen du trækker ud af Monte Carlo:


median_mass_simulation = empty(n_bins)
dispersion_mass_simulation_upper = empty(n_bins)
dispersion_mass_simulation_lower = empty(n_bins)

median_mass = empty(n_bins)
dispersion_mass_upper = empty(n_bins)
dispersion_mass_lower = empty(n_bins)

median_density_dm = empty(n_bins)
dispersion_density_dm_upper = empty(n_bins)
dispersion_density_dm_lower = empty(n_bins)

for nbins3 in range(0,n_bins):
    
    mass_simulation_sorted = sorted(mass_simulation[nbins3,:])
    median_mass_simulation[nbins3] = mass_simulation_sorted[int(len(mass_simulation_sorted)*0.5)]
    sigma_file1 = int(len(mass_simulation_sorted)*0.341)
    dispersion_mass_simulation_upper[nbins3] = mass_simulation_sorted[int(len(mass_simulation_sorted)*0.5)+sigma_file1] - median_mass_simulation[nbins3]
    dispersion_mass_simulation_lower[nbins3] = median_mass_simulation[nbins3] - mass_simulation_sorted[int(len(mass_simulation_sorted)*0.5)-sigma_file1]
    
    mass_sorted = sorted(mass_matrix[nbins3,:])
    median_mass[nbins3] = mass_sorted[int(len(mass_sorted)*0.5)]
    sigma_file1 = int(len(mass_sorted)*0.341)
    dispersion_mass_upper[nbins3] = mass_sorted[int(len(mass_sorted)*0.5)+sigma_file1] - median_mass[nbins3]
    dispersion_mass_lower[nbins3] = median_mass[nbins3] - mass_sorted[int(len(mass_sorted)*0.5)-sigma_file1]
    
    density_dm_sorted = sorted(density_dm_matrix[nbins3,:])
    median_density_dm[nbins3] = density_dm_sorted[int(len(density_dm_sorted)*0.5)]
    sigma_file1 = int(len(density_dm_sorted)*0.341)
    dispersion_density_dm_upper[nbins3] = density_dm_sorted[int(len(density_dm_sorted)*0.5)+sigma_file1] - median_density_dm[nbins3]
    dispersion_density_dm_lower[nbins3] = median_density_dm[nbins3] - density_dm_sorted[int(len(density_dm_sorted)*0.5)-sigma_file1]






densityplot = figure(17)
axdisp2 = densityplot.add_subplot(111)
#axdisp2.set_xscale("log", nonposx='clip')
axdisp2.set_yscale("log", nonposy='clip')
axdisp2.errorbar(median_radius,median_density_dm,yerr=[dispersion_density_dm_lower,dispersion_density_dm_upper],fmt='or')
axdisp2.errorbar(median_radius-20,median_density_dm_simulation,yerr=[dispersion_density_dm_simulation_lower,dispersion_density_dm_simulation_upper],fmt='ob')
grid()
#densityplot.show()







# Statistisk behandling af beta'erne:

median_beta_simulation = empty(n_bins)
dispersion_beta_simulation_upper = empty(n_bins)
dispersion_beta_simulation_lower = empty(n_bins)

median_beta_definition = empty(n_bins)
dispersion_beta_definition_upper = empty(n_bins)
dispersion_beta_definition_lower = empty(n_bins)

median_beta_jeans = empty(n_bins)
dispersion_beta_jeans_upper = empty(n_bins)
dispersion_beta_jeans_lower = empty(n_bins)

for n5 in range(0,n_bins):
    
    beta_simulation_sorted = sorted(beta_simulation[n5,:])
    median_beta_simulation[n5] = beta_simulation_sorted[int(len(beta_simulation_sorted)*0.5)]
    sigma_file1 = int(len(beta_simulation_sorted)*0.341)
    dispersion_beta_simulation_upper[n5] = beta_simulation_sorted[int(len(beta_simulation_sorted)*0.5)+sigma_file1] - median_beta_simulation[n5]
    dispersion_beta_simulation_lower[n5] = median_beta_simulation[n5] - beta_simulation_sorted[int(len(beta_simulation_sorted)*0.5)-sigma_file1]
    
    beta_definition_sorted = sorted(beta_def_matrix[n5,:])
    median_beta_definition[n5] = beta_definition_sorted[int(len(beta_definition_sorted)*0.5)]
    sigma_file1 = int(len(beta_definition_sorted)*0.341)
    dispersion_beta_definition_upper[n5] = beta_definition_sorted[int(len(beta_definition_sorted)*0.5)+sigma_file1] - median_beta_definition[n5]
    dispersion_beta_definition_lower[n5] = median_beta_definition[n5] - beta_definition_sorted[int(len(beta_definition_sorted)*0.5)-sigma_file1]
    
    beta_jeans_sorted = sorted(beta_jeans_matrix[n5,:])
    median_beta_jeans[n5] = beta_jeans_sorted[int(len(beta_jeans_sorted)*0.5)]
    sigma_file1 = int(len(beta_jeans_sorted)*0.341)
    dispersion_beta_jeans_upper[n5] = beta_jeans_sorted[int(len(beta_jeans_sorted)*0.5)+sigma_file1] - median_beta_jeans[n5]
    dispersion_beta_jeans_lower[n5] = median_beta_jeans[n5] - beta_jeans_sorted[int(len(beta_jeans_sorted)*0.5)-sigma_file1]



betaplot = figure(15)
axdisp3 = betaplot.add_subplot(111)
axdisp3.errorbar(median_radius,median_beta_simulation,yerr=[dispersion_beta_simulation_lower,dispersion_beta_simulation_upper],fmt='ob')
axdisp3.errorbar(median_radius+20,median_beta_definition,yerr=[dispersion_beta_definition_lower,dispersion_beta_definition_upper],fmt='or')
#axdisp3.errorbar(median_radius-20,median_beta_jeans,yerr=[dispersion_beta_jeans_lower,dispersion_beta_jeans_upper],fmt='og')
grid()
betaplot.show()




sys.exit(0)




