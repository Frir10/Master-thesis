# -*- coding: utf-8 -*-
# Dette program laver forskellige plots ud fra dataene for galaksen
# i Ex.2, del 3

from numpy import *
from matplotlib.pyplot import *


# Load Data:

data = loadtxt('obsweightny.txt')

Radius = data[:,1] #kpc
Density = data[:,2] #10^-25 g/cm^3
T_keV = data[:,3] #keV
Stjernemasse = data[:,4] #Msol

r = Radius*3.09*(10**19)
rho = Density*(10**-22)
T_Joule = T_keV*1.60218*(10**-16)
T_Kelvin = T_Joule/(1.38065*(10**-23))
M_stars = Stjernemasse*1.99*10**40

# Plot temperatur og densitet som funktion af radius

f1 = figure(1)
plot(r,T_Kelvin,'r-')
title('Temperatur i Kelvin')
f2= figure(2)
plot(r,rho,'rx')
ylabel('Densitet som funktion af radius')


# Plot massen af gas som funtkion af radius ved at summe over kugleskaller med 
# fast densitet.

r_average = empty(shape=(15,1))
r_average[14] = r[14]
for Number1 in range(0, 14):
   r_average[Number1] = (r[Number1]+r[Number1+1])/2

m_shell = empty(shape=(15,1))
for n2 in range(1, 14):
    m_shell[n2] = rho[n2] * 4/3 * pi * (r_average[n2]**3 - r_average[n2-1]**3)

m_shell[0] = rho[0] * 4/3 * pi * (r_average[0]**3)
m_shell[14] = rho[14] * 4/3 * pi * (r[14]**3 - r_average[13]**3)


M_shell = zeros(shape=(15, 1))
M_shell[0] = m_shell[0]
for n3 in range(1,15):
    M_shell[n3] = M_shell[n3-1] + m_shell[n3]

f3 = figure(3)
plot(r_average,M_shell)
title('Massen af gas, som funktion af radius')

# Her udregner jeg dlogT/dlogr og dlogrho/dlogr og plotter dem som funktion
# af radius. Jeg tager differenskvotienten mellem to punkter og midler
# derefter over differenskvotienterne på begge sider af et punkt. Jeg plotter 
# værdierne ud for de tilsvarende værdier af r.

T_logdlogr = zeros(shape=(15, 1))
rho_logdlogr = zeros(shape=(15,1))
for n4 in range(1,14):
    T_logdlogr[n4] = (log(T_Kelvin[n4]/T_Kelvin[n4-1])/log(r[n4]/r[n4-1]) + log(T_Kelvin[n4+1]/T_Kelvin[n4])/log(r[n4+1]/r[n4]))/2
    rho_logdlogr[n4] = (log(rho[n4]/rho[n4-1])/log(r[n4]/r[n4-1]) + log(rho[n4+1]/rho[n4])/log(r[n4+1]/r[n4]))/2

f4 = figure(4)
plot(r[1:14],T_logdlogr[1:14])
title('dlog(T)/dlog(r)')

f5 = figure(5)
plot(r[1:14],rho_logdlogr[1:14])
title('dlog(rho)/dlog(r)')

hej = transpose(-rho_logdlogr-T_logdlogr)
M_tot = (transpose(hej*T_Joule*r))/(0.59*6.67*1.67*1.99*(10**-8))
f6 = figure(6)
plot(r[1:14],M_tot[1:14])
# Estimerer massen i den sidste bin, som at have hældning som imellem de sidste to:
M_tot[14] = -(r[14]*T_Joule[14]/(0.59*6.67*1.67*1.99*(10**-8)))*(log(T_Kelvin[14]/T_Kelvin[13])/log(r[14]/r[13])+log(rho[14]/rho[13])/log(r[14]/r[13]))
M_tot[0] = -(r[0]*T_Joule[0]/(0.59*6.67*1.67*1.99*(10**-8)))*(log(T_Kelvin[1]/T_Kelvin[0])/log(r[1]/r[0])+log(rho[1]/rho[0])/log(r[1]/r[0]))

# Lav et plot af stjernemassen som funktion af radius:

f7 = figure(7)
plot(r,M_stars)

# Plot massen af stjerner i solmasser

f8 = figure(8)
plot(r,(M_stars/(1.99*10**30)))

# Plot de tre kurver sammen
f9 = figure(9)
plot(r,M_tot,'r-',r,M_stars/(1.99*10**30),'b-',r,M_shell/(1.99*10**30),'g-')

# Udregning af den logaritmiske hældning af massen i det sidste interval

dMdr_ydre = log(M_tot[14]/M_tot[13])/(log(r[14]/r[13]))

# Plotning af den indre region, under 1 kpc af de forsekllige masser
f10 = figure(10)
plot(r[0:6],M_tot[0:6],'r-',r[0:6],M_stars[0:6]/(1.99*10**30),'b-',r[0:6],M_shell[0:6]/(1.99*10**30),'g-')

# Fit den totale masse med en sum af gas- og stjernemassen og læg en konstant til, for
# at se om der er et sort hul
M_baryonisk = transpose(M_shell[0:6]/(1.99*10**30))+M_stars[0:6]/(1.99*10**30)+10**9
f11 = figure(11)
plot(r[0:6],M_tot[0:6],'r-',r[0:6],M_stars[0:6]/(1.99*10**30),'b-',r[0:6],M_shell[0:6]/(1.99*10**30),'g-',r[0:6],transpose(M_baryonisk),'k-')


# Lav plot af alle masserne logaritmisk
#M_tot_log = log(M_tot)
#M_stars_log = log10(transpose(M_stars))
#M_gas_log = log(M_shell/(1.99*10**30))
#f12 = figure(12)
#plot(r,log(M_tot),'r-',r,log(M_stars/(1.99*10**30)),'b-',r,log(M_shell/(1.99*10**30)),'g-')

# Nyt estimat af den totale masse ud fra drho/dr og dT/dr i stedet:

# Vælg hvilke figurer der skal vises:

#f1.show()
#f2.show()
#f3.show()
#f4.show()
#f5.show()
#f6.show()
#f7.show()
#f8.show()
#f9.show()
#f10.show()
f11.show()
#f12.show()

print M_baryonisk