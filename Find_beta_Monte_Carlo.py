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
H = 70.4  # Hubble parameter
rho_c = 3*H**2/(8*pi*G)*(kpc/M_sol) # Kritisk tæthed af universet


# Load all the data:


n_bins = 19  # Final number of bins per profile
n_bins1 = n_bins+1
n_profiles = 51  # Total number of profiles

radius_central = empty(shape=(n_bins,n_profiles))

binned_density_gas = empty(shape=(n_bins,n_profiles))
binned_T_gas = empty(shape=(n_bins,n_profiles))

beta = empty(shape=(n_bins,n_profiles))



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
    
    binned_density_dm = empty(n_bins)
    binned_density_tot = empty(n_bins)
    
    binned_velocitydispersion_r = empty(n_bins)
    binned_velocitydispersion_theta = empty(n_bins)
    binned_velocitydispersion_phi = empty(n_bins)
    
    bind = len(rho_dm)/n_bins
    
    for newbins in range(0,n_bins):
        
        radius_central[newbins,nfile] = sum(radius[(10+newbins*bind):(10+bind+newbins*bind)])/bind
        binned_density_dm[newbins] = sum(rho_dm[(10+newbins*bind):(10+bind+newbins*bind)])/bind
        
        binned_velocitydispersion_r[newbins] = sum(sigma_r[(10+newbins*bind):(10+bind+newbins*bind)])/bind
        binned_velocitydispersion_theta[newbins] = sum(sigma_theta[(10+newbins*bind):(10+bind+newbins*bind)])/bind
        binned_velocitydispersion_phi[newbins] = sum(sigma_phi[(10+newbins*bind):(10+bind+newbins*bind)])/bind
        
        binned_density_gas[newbins,nfile] = sum(rho_gas[(10+newbins*bind):(10+bind+newbins*bind)])/bind
        binned_T_gas[newbins,nfile] = sum(T_gas[(10+newbins*bind):(10+bind+newbins*bind)])/bind
        
        beta[newbins,nfile] = 1 - (((binned_velocitydispersion_theta[newbins])**2+(binned_velocitydispersion_phi[newbins])**2)/(2*(binned_velocitydispersion_r[newbins])**2))
        
    


radius_relative = arange((0.0+2.0/len(rho_dm)),2.0-2.0/n_bins1,(2.0/n_bins1))+(1.0/n_bins1)

radius_interpol = arange(110,2700,0.1)

        
# Lav et datasæt der svarer til eksperimentel måling af X-ray stråling fra en hob der svarer til alle simulationerne stacked.
# Du har temperatur og tæthed af gassen ved hver radius. 

median_radius = empty(n_bins)
dispersion_radius_upper = empty(n_bins)
dispersion_radius_lower = empty(n_bins)

median_density_gas = empty(n_bins)
dispersion_density_gas_upper = empty(n_bins)
dispersion_density_gas_lower = empty(n_bins)

median_T_gas = empty(n_bins)
dispersion_T_gas_upper = empty(n_bins)
dispersion_T_gas_lower = empty(n_bins)


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



# Lav Monte Carlo sampling:


density_dispersion = (dispersion_density_gas_upper+dispersion_density_gas_lower)/2.
T_dispersion = (dispersion_T_gas_upper+dispersion_T_gas_lower)/2.

density_gas_gauss = zeros(n_bins)
T_gas_gauss = zeros(n_bins)

sampling_size = 5
mass_matrix = empty(shape=(sampling_size,n_bins))
S_matrix = empty(shape=(sampling_size,n_bins))
S_mass_matrix = empty(shape=(sampling_size,n_bins))

for n1000 in range(0,sampling_size):
    
    # Lav gaussisk fordelte arrays af temperatur og densitet
    
    for n1 in range(0,n_bins):
        density_gas_gauss[n1] = random.normal(median_density_gas[n1],density_dispersion[n1])
        T_gas_gauss[n1] = random.normal(median_T_gas[n1],T_dispersion[n1])
        
    # Udregn logaritmiske hældninger for T og rho:
    
    gamma_gas = empty(n_bins) # dlog_rho_dlog_r
    eta_gas = empty(n_bins)   # dlog_T_dlog_r
    
    gamma_gas[1:n_bins-1] = (log(density_gas_gauss[1:n_bins-1]/density_gas_gauss[0:n_bins-2])/(log(median_radius[1:n_bins-1]/median_radius[0:n_bins-2]))+log(density_gas_gauss[2:n_bins]/density_gas_gauss[1:n_bins-1])/(log(median_radius[2:n_bins]/median_radius[1:n_bins-1])))/2
    gamma_gas[0] = (3*log(density_gas_gauss[1]/density_gas_gauss[0])/log(median_radius[1]/median_radius[0])+log(density_gas_gauss[2]/density_gas_gauss[1])/log(median_radius[2]/median_radius[1]))/4
    gamma_gas[n_bins-1] = (3*log(density_gas_gauss[n_bins-1]/density_gas_gauss[n_bins-2])/log(median_radius[n_bins-1]/median_radius[n_bins-2])+log(density_gas_gauss[n_bins-2]/density_gas_gauss[n_bins-3])/log(median_radius[n_bins-2]/median_radius[n_bins-3]))/4
    
    eta_gas[1:n_bins-1] = (log(T_gas_gauss[1:n_bins-1]/T_gas_gauss[0:n_bins-2])/(log(median_radius[1:n_bins-1]/median_radius[0:n_bins-2]))+log(T_gas_gauss[2:n_bins]/T_gas_gauss[1:n_bins-1])/(log(median_radius[2:n_bins]/median_radius[1:n_bins-1])))/2
    eta_gas[0] = (3*log(T_gas_gauss[1]/T_gas_gauss[0])/log(median_radius[1]/median_radius[0])+log(T_gas_gauss[2]/T_gas_gauss[1])/log(median_radius[2]/median_radius[1]))/4
    eta_gas[n_bins-1] = (3*log(T_gas_gauss[n_bins-1]/T_gas_gauss[n_bins-2])/log(median_radius[n_bins-1]/median_radius[n_bins-2])+log(T_gas_gauss[n_bins-2]/T_gas_gauss[n_bins-3])/log(median_radius[n_bins-2]/median_radius[n_bins-3]))/4
    
    
    # Find masseleddet ud fra hydrostatisk ligevægt:
    
    mass_term = -T_gas_gauss*k_b/(m_p*mean_molecular_weight)*(gamma_gas + eta_gas)
    mass = mass_term*median_radius*kpc/(G*M_sol)
    mass_matrix[n1000,:] = mass
    
    
    # Lav interpolation på de udregnede værdier for massen for at finde virialradius:
    
    mass_interpol_func = interp1d(median_radius,mass)
    mass_interpol = mass_interpol_func(radius_interpol)
    density_average = mass_interpol / (4/3*pi*radius_interpol**3)
    
    find_virial = min(enumerate(density_average), key=lambda x: abs(x[1]-(200*rho_c))) # Finder den indgang hvor densiteten er tættest på 200 gange den kritiske
    
    r_virial = radius_interpol[find_virial[0]]
    mass_virial = mass_interpol[find_virial[0]]
    v_virial = sqrt(G*mass_virial*M_sol/(r_virial*kpc))
    
    
    # Udregn nu det samme hvor du gætter på et S:
    
    S_term_list = zeros(n_bins)
    
    iterations = 10
    
    for iteration in range(0,iterations):
        
        mass_term = -T_gas_gauss*k_b/(m_p*mean_molecular_weight)*(gamma_gas + eta_gas) - S_term_list
        mass = mass_term*median_radius*kpc/(G*M_sol)
        
        mass_interpol_func = interp1d(median_radius,mass)
        mass_interpol = mass_interpol_func(radius_interpol)
        density_average = mass_interpol / (4/3*pi*radius_interpol**3)
        
        find_virial = min(enumerate(density_average), key=lambda x: abs(x[1]-(200*rho_c))) # Finder den indgang hvor densiteten er tættest på 200 gange den kritiske
        
        r_virial = radius_interpol[find_virial[0]]
        mass_virial = mass_interpol[find_virial[0]]
        v_virial = sqrt(G*mass_virial*M_sol/(r_virial*kpc))
        
        S_term_list = S_term(median_radius,chosen_parameters[0],chosen_parameters[1],chosen_parameters[2],chosen_parameters[3],v_virial,r_virial)
        
        print r_virial
    
    #S_term_list = S_term(median_radius,t*chosen_parameters[0],chosen_parameters[1],chosen_parameters[2],chosen_parameters[3],v_virial,r_virial)
    #S_matrix[n1000,:] = S_term_list*median_radius*kpc/(G*M_sol)
    
    #mass_term_with_s = -T_gas_gauss*k_b/(m_p*mean_molecular_weight)*(gamma_gas + eta_gas) - S_term_list
    #mass_with_s = mass_term_with_s*median_radius*kpc/(G*M_sol)
    #S_mass_matrix[n1000,:] = mass_with_s
    
    # Udregn den totale densitet og dermed densiteten af mørkt stof:
    
    #density_tot = empty(n_bins)
    #density_tot[1:n_bins-1] = ((mass_with_s[1:n_bins-1]-mass_with_s[0:n_bins-2])/(median_radius[1:n_bins-1]-median_radius[0:n_bins-2])+(mass_with_s[2:n_bins]-mass_with_s[1:n_bins-1])/(median_radius[2:n_bins]-median_radius[1:n_bins-1]))/(2*4*pi*(median_radius[1:n_bins-1])**2)
    #density_tot[0] = (3*(mass_with_s[1]-mass_with_s[0])/(median_radius[1]-median_radius[0])+(mass_with_s[2]-mass_with_s[1])/(median_radius[2]-median_radius[1]))/(4*4*pi*(median_radius[0])**2)
    #density_tot[n_bins-1] = (3*(mass_with_s[n_bins-1]-mass_with_s[n_bins-2])/(median_radius[n_bins-1]-median_radius[n_bins-2])+(mass_with_s[n_bins-2]-mass_with_s[n_bins-3])/(median_radius[n_bins-2]-median_radius[n_bins-3]))/(4*4*pi*(median_radius[n_bins-1])**2)
    #density_dm = density_tot-median_density_gas
    
    #gamma_dm = empty(n_bins)
    #gamma_dm[1:n_bins-1] = (log(density_dm[1:n_bins-1]/density_dm[0:n_bins-2])/(log(median_radius[1:n_bins-1]/median_radius[0:n_bins-2]))+log(density_dm[2:n_bins]/density_dm[1:n_bins-1])/(log(median_radius[2:n_bins]/median_radius[1:n_bins-1])))/2
    #gamma_dm[0] = (3*log(density_dm[1]/density_dm[0])/log(median_radius[1]/median_radius[0])+log(density_dm[2]/density_dm[1])/log(median_radius[2]/median_radius[1]))/4
    #gamma_dm[n_bins-1] = (3*log(density_dm[n_bins-1]/density_dm[n_bins-2])/log(median_radius[n_bins-1]/median_radius[n_bins-2])+log(density_dm[n_bins-2]/density_dm[n_bins-3])/log(median_radius[n_bins-2]/median_radius[n_bins-3]))/4
    
    # Udregn nu den radielle hastighedsdispersion for mørkt stof, idet du antager at kende temperatur-relationen:
    
    #T_dm = T_gas_gauss + difffunc(median_radius/r_virial,diff_parameters[0],diff_parameters[1],diff_parameters[2],diff_parameters[3])
    
    #psi_integrand = T_gas_gauss*k_b/(m_p*mean_molecular_weight)*(gamma_gas + eta_gas) + 3*k_b/(m_p*mean_molecular_weight)*T_dm
    #integrand = 10**-10*psi_integrand*density_dm*median_radius**2
    
    #integral_piece = empty(n_bins)
    #integral_piece[0] = 0.5*integrand[0]*median_radius[0]
    
    #dispersion_integral = empty(n_bins)
    #dispersion_integral[0] = integral_piece[0]
    
    #radial_velocity_variance = empty(n_bins)
    #radial_velocity_variance[0] = 10**10/(median_radius[0]**3*density_dm[0])*dispersion_integral[0]
    
    #for n_integrate in range(1,n_bins):
            
     #   integral_piece[n_integrate] = 0.5*(integrand[n_integrate]+integrand[n_integrate-1])*(median_radius[n_integrate]-median_radius[n_integrate-1])
      #  dispersion_integral[n_integrate] = dispersion_integral[n_integrate-1] + integral_piece[n_integrate]
       # radial_velocity_variance[n_integrate] = 10**10/(median_radius[n_integrate]**3*density_dm[n_integrate])*dispersion_integral[n_integrate]
    
    #print sqrt(radial_velocity_variance)
    

sys.exit(0)


# Udregn logaritmiske hældninger for T og rho:

gamma_gas = empty(n_bins) # dlog_rho_dlog_r
eta_gas = empty(n_bins)   # dlog_T_dlog_r


gamma_gas[1:n_bins-1] = (log(median_density_gas[1:n_bins-1]/median_density_gas[0:n_bins-2])/(log(median_radius[1:n_bins-1]/median_radius[0:n_bins-2]))+log(median_density_gas[2:n_bins]/median_density_gas[1:n_bins-1])/(log(median_radius[2:n_bins]/median_radius[1:n_bins-1])))/2

gamma_gas[0] = (3*log(median_density_gas[1]/median_density_gas[0])/log(median_radius[1]/median_radius[0])+log(median_density_gas[2]/median_density_gas[1])/log(median_radius[2]/median_radius[1]))/4
gamma_gas[n_bins-1] = (3*log(median_density_gas[n_bins-1]/median_density_gas[n_bins-2])/log(median_radius[n_bins-1]/median_radius[n_bins-2])+log(median_density_gas[n_bins-2]/median_density_gas[n_bins-3])/log(median_radius[n_bins-2]/median_radius[n_bins-3]))/4


eta_gas[1:n_bins-1] = (log(median_T_gas[1:n_bins-1]/median_T_gas[0:n_bins-2])/(log(median_radius[1:n_bins-1]/median_radius[0:n_bins-2]))+log(median_T_gas[2:n_bins]/median_T_gas[1:n_bins-1])/(log(median_radius[2:n_bins]/median_radius[1:n_bins-1])))/2

eta_gas[0] = (3*log(median_T_gas[1]/median_T_gas[0])/log(median_radius[1]/median_radius[0])+log(median_T_gas[2]/median_T_gas[1])/log(median_radius[2]/median_radius[1]))/4
eta_gas[n_bins-1] = (3*log(median_T_gas[n_bins-1]/median_T_gas[n_bins-2])/log(median_radius[n_bins-1]/median_radius[n_bins-2])+log(median_T_gas[n_bins-2]/median_T_gas[n_bins-3])/log(median_radius[n_bins-2]/median_radius[n_bins-3]))/4



# Definer funktionen S (mere specifikt S*r):

dm_parameters = [0.11255906,12.8245615,2.2179*10**-6,0.31763796]
gas_parameters = [0.13669927,12.635266,6.3028*10**-6,0.2453202]

chosen_parameters = dm_parameters


def S_term(r,alpha,a,c,d,v_virial,r_virial):
    
    v_p = -alpha/(((r/r_virial)**(-a) + c*(r/r_virial)**(a/2))**(1/a)-d)*v_virial
    
    v_p_derivative = v_virial/r_virial*alpha/a*(a/2*c*(r/r_virial)**(a/2-1)-a*(r/r_virial)**(-a-1))*(((r/r_virial)**(-a)+c*(r/r_virial)**(a/2))**(1/a-1))/(((r/r_virial)**(-a)+c*(r/r_virial)**(a/2))**(1/a) - d)**2
    
    return v_p*v_p_derivative*r + v_p*H*r + v_p_derivative*r**2*H


# Find masseleddet ud fra hydrostatisk ligevægt og iterer:

S_term_list = zeros(n_bins)
radius_more = arange(1.,5001.,1.)*(median_radius[n_bins-1]/5001.)

iterations = 2

for iteration in range(0,iterations):
    
    mass_term = -median_T_gas*k_b/(m_p*mean_molecular_weight)*(gamma_gas + eta_gas) - S_term_list
    mass = mass_term*median_radius*kpc/(G*M_sol)
    
    mass_interpol_func = interp1d(median_radius,mass)
    mass_interpol = mass_interpol_func(radius_more[300:4500])
    density_average = mass_interpol / (4/3*pi*radius_more[300:4500]**3)
    
    find_virial = min(enumerate(density_average), key=lambda x: abs(x[1]-(200*rho_c))) # Finder den indgang hvor densiteten er tættest på 200 gange den kritiske
    
    r_virial = radius_more[find_virial[0]+300]
    mass_virial = mass_interpol[find_virial[0]+300]
    v_virial = sqrt(G*mass_virial*M_sol/(r_virial*kpc))

    
    S_term_list = S_term(median_radius,chosen_parameters[0],chosen_parameters[1],chosen_parameters[2],chosen_parameters[3],v_virial,r_virial)
    
    print r_virial


testplot = figure(81)
plot(median_radius,mass,'*')
grid()
#testplot.show()

# Lav en model for masseprofilen som du kan bruge til at bestemme densiteten med:

def mass_fit(r,a,b,c):
    
    return 10**14*(a*r**2 + b*r + c)

def mass_fit_2(r,mass_tot,a,b,c):
    
    return mass_tot - a/((r+1000.)/b)**c
    

def mass_fit_3(r,a,b,c):
    
    return 10**14*(c*(r/a)*(r/a+1)**-b)

def mass_fit_4(r,a,b,c,d):
    
    return 10**14*(c*((r/a)**d*(r/a+1)**-b)**(1/d))

def chi2_mass(params,r,data):
    a = params['a'].value
    b = params['b'].value
    c = params['c'].value
    #mass_tot = params['mass_tot'].value
    d = params['d'].value
    #e = params['e'].value
    
    model = mass_fit_4(r,a,b,c,d)

    return (data-model)

params = Parameters()

#params.add('mass_tot',value = 10**14.,min=0.)
params.add('a',value = 1.,min=0.)
params.add('b',value = 1,min=0.)
params.add('c',value = 1.,max=3.5)
params.add('d',value = 1.)#,max=0.)
#params.add('e',value = 1.)#,max=0.)


result = minimize(chi2_mass,params,args=(median_radius[0:n_bins-1],mass[0:n_bins-1]))

#report_fit(params)

mass_fit_result = mass_fit_4(median_radius,params['a'].value,params['b'].value,params['c'].value,params['d'].value)

testplot2 = figure(82)
plot(median_radius,mass,'*r',median_radius,mass_fit_result,'-g')
grid()
#testplot2.show()

# Udregn den totale densitet og dermed densiteten af mørkt stof:

density_tot = empty(n_bins)

# Et volumenelement har formen dV = 4*pi*r^2 dr, så differentiationen bliver:
density_tot[1:n_bins-1] = ((mass_fit_result[1:n_bins-1]-mass_fit_result[0:n_bins-2])/(median_radius[1:n_bins-1]-median_radius[0:n_bins-2])+(mass_fit_result[2:n_bins]-mass_fit_result[1:n_bins-1])/(median_radius[2:n_bins]-median_radius[1:n_bins-1]))/(2*4*pi*(median_radius[1:n_bins-1])**2)

density_tot[0] = (3*(mass_fit_result[1]-mass_fit_result[0])/(median_radius[1]-median_radius[0])+(mass_fit_result[2]-mass_fit_result[1])/(median_radius[2]-median_radius[1]))/(4*4*pi*(median_radius[0])**2)
density_tot[n_bins-1] = (3*(mass_fit_result[n_bins-1]-mass_fit_result[n_bins-2])/(median_radius[n_bins-1]-median_radius[n_bins-2])+(mass_fit_result[n_bins-2]-mass_fit_result[n_bins-3])/(median_radius[n_bins-2]-median_radius[n_bins-3]))/(4*4*pi*(median_radius[n_bins-1])**2)

density_dm = density_tot-median_density_gas


# Find den afledte af densiteten for mørkt stof, dvs. gamma_dm:

gamma_dm = empty(n_bins)

gamma_dm[1:n_bins-1] = (log(density_dm[1:n_bins-1]/density_dm[0:n_bins-2])/(log(median_radius[1:n_bins-1]/median_radius[0:n_bins-2]))+log(density_dm[2:n_bins]/density_dm[1:n_bins-1])/(log(median_radius[2:n_bins]/median_radius[1:n_bins-1])))/2

gamma_dm[0] = (3*log(density_dm[1]/density_dm[0])/log(median_radius[1]/median_radius[0])+log(density_dm[2]/density_dm[1])/log(median_radius[2]/median_radius[1]))/4
gamma_dm[n_bins-1] = (3*log(density_dm[n_bins-1]/density_dm[n_bins-2])/log(median_radius[n_bins-1]/median_radius[n_bins-2])+log(density_dm[n_bins-2]/density_dm[n_bins-3])/log(median_radius[n_bins-2]/median_radius[n_bins-3]))/4

#print gamma_dm


# Udregn nu den radielle hastighedsdispersion for mørkt stof, idet du antager at kende temperatur-relationen:

def difffunc(r,q,p,o,a):
    
    return q - o/((r+1.3)/a)**p

diff_parameters = [4.0158e+06,16.8419116,6.0739e+05,1.67671123]

T_dm = median_T_gas + difffunc(median_radius/r_virial,diff_parameters[0],diff_parameters[1],diff_parameters[2],diff_parameters[3])

psi_integrand = median_T_gas*k_b/(m_p*mean_molecular_weight)*(gamma_gas + eta_gas) + 3*k_b/(m_p*mean_molecular_weight)*T_dm

integrand = 10**-10*psi_integrand*density_dm*median_radius**2

integral_piece = empty(n_bins)
integral_piece[0] = 0.5*integrand[0]*median_radius[0]

dispersion_integral = empty(n_bins)
dispersion_integral[0] = integral_piece[0]

radial_velocity_variance = empty(n_bins)
radial_velocity_variance[0] = 10**10/(median_radius[0]**3*density_dm[0])*dispersion_integral[0]

for n_integrate in range(1,n_bins):
    
    integral_piece[n_integrate] = 0.5*(integrand[n_integrate]+integrand[n_integrate-1])*(median_radius[n_integrate]-median_radius[n_integrate-1])
    
    dispersion_integral[n_integrate] = dispersion_integral[n_integrate-1] + integral_piece[n_integrate]
    
    radial_velocity_variance[n_integrate] = 10**10/(median_radius[n_integrate]**3*density_dm[n_integrate])*dispersion_integral[n_integrate]
    
#print integral_piece


# Brug at beta = 0 i centrum, samt definitionerne af beta og kappa, til at bestemme en grænseværdi for sigma_r_^2, så vi kan undgå de negative værdier.

beta_inderst = 0.

sigma_r_inderst = (T_dm[0]) * 3*k_b/(mean_molecular_weight*m_p * (3 - 2 * beta_inderst))
    
#print sigma_r_inderst

radial_velocity_variance = radial_velocity_variance - radial_velocity_variance[0] + sigma_r_inderst 

sigma_r_fundet = sqrt(radial_velocity_variance)


testfigure6 = figure(86)
plot(median_radius,radial_velocity_variance)
grid()
#testfigure6.show()

testfigure = figure(41)
plot(radius_relative,sigma_r_fundet,'*')
grid()
#testfigure.show()


# Find nu beta på to måder og sammenlign med den fra simulationen:


# Først ved hjælp af definitionen af delta og beta:

beta_definition = 1.5 * (1. - k_b/(mean_molecular_weight*m_p)*T_dm/radial_velocity_variance)

#print beta_definition


# Så ved hjælp af Jean's ligning:

eta_dm = empty(n_bins)

eta_dm[1:n_bins-1] = (log(radial_velocity_variance[1:n_bins-1]/radial_velocity_variance[0:n_bins-2])/(log(median_radius[1:n_bins-1]/median_radius[0:n_bins-2]))+log(radial_velocity_variance[2:n_bins]/radial_velocity_variance[1:n_bins-1])/(log(median_radius[2:n_bins]/median_radius[1:n_bins-1])))/2
eta_dm[0] = (3*log(radial_velocity_variance[1]/radial_velocity_variance[0])/log(median_radius[1]/median_radius[0])+log(radial_velocity_variance[2]/radial_velocity_variance[1])/log(median_radius[2]/median_radius[1]))/4
eta_dm[n_bins-1] = (3*log(radial_velocity_variance[n_bins-1]/radial_velocity_variance[n_bins-2])/log(median_radius[n_bins-1]/median_radius[n_bins-2])+log(radial_velocity_variance[n_bins-2]/radial_velocity_variance[n_bins-3])/log(median_radius[n_bins-2]/median_radius[n_bins-3]))/4

beta_jeans = (-G*mass_fit_result*M_sol/(median_radius*kpc*radial_velocity_variance) - eta_dm - gamma_dm)

print eta_dm


# Find beta fra simulationen og plot den sammen med de to andre for at sammenligne:

median_beta = empty(n_bins)
dispersion_beta_upper = empty(n_bins)
dispersion_beta_lower = empty(n_bins)


for nbins9 in range(0,n_bins):
    
    beta_sorted = sorted(beta[nbins9,:])
    median_beta[nbins9] = beta_sorted[int(len(beta_sorted)*0.5)]
    sigma_file1 = int(len(beta_sorted)*0.341)
    dispersion_beta_upper[nbins9] = beta_sorted[int(len(beta_sorted)*0.5)+sigma_file1] - median_beta[nbins9]
    dispersion_beta_lower[nbins9] = median_beta[nbins9] - beta_sorted[int(len(beta_sorted)*0.5)-sigma_file1]


betafigure = figure(33)
plot(median_radius,median_beta,'b',median_radius,beta_jeans,'g',median_radius,beta_definition,'r')
grid()

betafigure.show()

