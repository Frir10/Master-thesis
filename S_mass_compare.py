# -*- coding: utf-8 -*-
# Try and load some of the data and look at it

from numpy import *
from matplotlib.pyplot import *
from matplotlib.mlab import *
from scipy.stats import *


# Some constants:

k_b = 1.3806488*(10**-3)
m_p = 1.6726178*(10**-7)
G = 6.67384
M_sol = 1.98855
kpc = 3.08567758
H = 70.4  # Hubble parameter
rho_c = 3*H**2/(8*pi*G)*(kpc/M_sol) # Kritisk tæthed af universet


n_bins = 50  # Final number of bins per profile
n_profiles = 51  # Total number of profiles

radius_relative = arange(0.0,2.0,(2.0/n_bins))+(1.0/n_bins)

radius_physical = arange(1.,2001.,1.)  # Fysisk radius målt i kiloparsec


# Se Infall_velocity_fit_Martina_new hvis her mangler noget


# Konstruer en Hernquist-profil for det mørke stof:

M_hob = 2. * 10**14  # den totale masse af hoben, i solmasser
a_hob = 500. # skalalængden for profilen, målt i kpc

def hernquist(r,a,M):
    
    return M/((2*pi*a**3)*(r/a)*(r/a+1)**3)

def hernquist_mass(r,a,M):
    
    return M*((r/a)**2/(r/a+1)**2)

density_hernquist = hernquist(radius_physical,a_hob,M_hob)

mass_hernquist = hernquist_mass(radius_physical,a_hob,M_hob)

density_hernquist_average = mass_hernquist/(4/3*pi*radius_physical**3)


find_virial = min(enumerate(density_hernquist_average), key=lambda x: abs(x[1]-(200*rho_c))) # Finder den indgang hvor densiteten er tættest på 200 gange den kritiske

#print find_virial

r_virial = radius_physical[find_virial[0]]
mass_virial = mass_hernquist[find_virial[0]]

v_virial = sqrt(G*mass_virial*M_sol/(r_virial*kpc))

print r_virial


# Definer vores funktion for indfaldshastigheden som fittet af enten det mørke stof eller gassen:

def martinafunc(r,alpha,a,c,d):
    
    return -alpha/(((r/r_virial)**(-a) + c*(r/r_virial)**(a/2))**(1/a)-d)


#best_fit_martina_dm_grav = martinafunc(radius_relative,0.11255906,12.8245615,2.2179*10**-6,0.31763796)

#best_fit_martina_gas_grav = martinafunc(radius_relative,0.13669927,12.635266,6.3028*10**-6,0.2453202)

dm_parameters = [0.14911711,25.7009899,2.4061e-11,0.30949707]
gas_parameters = [0.13669927,12.635266,6.3028*10**-6,0.2453202]

chosen_parameters = dm_parameters

v_p_fit = martinafunc(radius_physical,chosen_parameters[0],chosen_parameters[1],chosen_parameters[2],chosen_parameters[3])

v_p_physical = v_p_fit*v_virial


# Definer en funktion der udregner S ud fra v_p:


# Definer først den afledte af hastigheden, fundet ved diffrentiering i hånden og tjekket med wolfram alpha:

def vp_derivative(r,alpha,a,c,d):
    
    return alpha/a*(a/2*c*(r/r_virial)**(a/2-1)-a*(r/r_virial)**(-a-1))*(((r/r_virial)**(-a)+c*(r/r_virial)**(a/2))**(1/a-1))/(((r/r_virial)**(-a)+c*(r/r_virial)**(a/2))**(1/a) - d)**2

vp_derivative_physical = v_virial/r_virial*vp_derivative(radius_physical,chosen_parameters[0],chosen_parameters[1],chosen_parameters[2],chosen_parameters[3])

# Udregn de tre led der udgør S:

first_term = H*radius_physical*v_p_physical
second_term = vp_derivative_physical*radius_physical**2*H
third_term = v_p_physical*vp_derivative_physical*radius_physical


# Læg leddene sammen og sammenlign med leddet GM/r

S_times_r = first_term + second_term + third_term

mass_term = mass_hernquist*M_sol*G/(radius_physical*kpc)

testplot2 = figure(84)
axdisp1 = testplot2.add_subplot(111)
axdisp1.plot(radius_physical/r_virial,S_times_r/mass_term,'g')#,radius_physical,mass_term,'r')
grid()
#legend(loc='lower left',numpoints=1,fontsize=24)
#xlabel('$r/r_{200}$', fontsize=22)
setp(axdisp1.get_xticklabels(), visible=False)
ylabel('$r^2S/(GM)$', fontsize=22)
testplot2.show()


testplot = figure(83)
axdisp2 = testplot.add_subplot(111)
axdisp2.plot(radius_physical/r_virial,first_term/mass_term,'b',label='$r^2 H_0 \overline{v_p}(r) /(GM(r))$')
axdisp2.plot(radius_physical/r_virial,second_term/mass_term,'r',label='$r^3 H_0 \partial\overline{v_p}/\partial r /(GM(r))$')
axdisp2.plot(radius_physical/r_virial,third_term/mass_term,'m',label='$r^2 \overline{v_p}(r) \partial\overline{v_p}/\partial r/(GM(r))$')

grid()
legend(loc='upper left',numpoints=1,fontsize=24)
xlabel('$r/r_{200}$', fontsize=22)
#ylabel('$r^2S/(GM)$', fontsize=22)
testplot.show()


# Prøv at definere S som en funktion i fysiske enheder

def S_term(r,alpha,a,c,d,v_virial,r_virial):
    
    v_p = -alpha/(((r/r_virial)**(-a) + c*(r/r_virial)**(a/2))**(1/a)-d)*v_virial
    
    v_p_derivative = v_virial/r_virial*alpha/a*(a/2*c*(r/r_virial)**(a/2-1)-a*(r/r_virial)**(-a-1))*(((r/r_virial)**(-a)+c*(r/r_virial)**(a/2))**(1/a-1))/(((r/r_virial)**(-a)+c*(r/r_virial)**(a/2))**(1/a) - d)**2
    
    return v_p*v_p_derivative*r + v_p*H*r + v_p_derivative*r**2*H


#print S_term(radius_physical,chosen_parameters[0],chosen_parameters[1],chosen_parameters[2],chosen_parameters[3],v_virial,r_virial) - S_times_r

testfigure = figure(91)
plot(radius_physical,S_times_r)
grid()
testfigure.show()

