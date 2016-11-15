# -*- coding: utf-8 -*-
# Lav et progam der kan lave et kunstigt spektrum for CO-molekyle ud fra konstanter
# opnået fra måling af et spektrum

from numpy import *
from matplotlib.pyplot import *

B_e = 1.9319
alpha_e = 0.0175
omega_e = 2.1698136*(10**3)
omega_e_x_e = 13.2883
D_e = 6.12147*(10**(-6))

v_max = 5
J_max = 60

T = 3500
c = 299792458.0
h = 6.626*(10**-34)
k_b = 1.38*(10**-23)


# Definer liste med energier og vægte for værdier af kvantetallene J og v

Energy = zeros(shape=(v_max,J_max))
weight = zeros(shape=(v_max,J_max))

for v in range(0,v_max):
    for J in range(1,J_max):
        
        Energy[v,J] = omega_e*(v+0.5)-omega_e_x_e*((v+0.5)**2)+(B_e-alpha_e*(v+0.5))*J*(J+1)-D_e*J**2*(J+1)**2
        weight[v,J] = (2*J+1)*exp(-Energy[v,J]*100*c*h/(k_b*T)) # Sandsynligheden for at befinde sig i et niveau med den givne
                                                                # energi ganget med udartetheden af det energiniveau


# Lav lister med værdier for overgangene for v -> v+1 (nu for alle v) og J gående 
# både J -> J+1 og J -> J-1

omega_vp1_J_plus1 = zeros(shape=(v_max-1,J_max-2))
omega_vp1_J_minus1 = zeros(shape=(v_max-1,J_max-2))

for v2 in range(0,v_max-1):
    for J2 in range(0,J_max-2):
        
        omega_vp1_J_plus1[v2][J2] = Energy[v2+1][J2+2] - Energy[v2][J2+1]
        omega_vp1_J_minus1[v2][J2] = Energy[v2+1][J2+1] - Energy[v2][J2+2]


# Find nu bølgetallene for overgange svarende til at v går fra v -> v+2 og
# J gør det samme som før

omega_vp2_J_plus1 = zeros(shape=(v_max-2,J_max-2))
omega_vp2_J_minus1 = zeros(shape=(v_max-2,J_max-2))

for v3 in range(0,v_max-2):
    for J3 in range(0,J_max-2):
        
        omega_vp2_J_plus1[v3][J3] = Energy[v3+2][J3+2] - Energy[v3][J3+1]
        omega_vp2_J_minus1[v3][J3] = Energy[v3+2][J3+1] - Energy[v3][J3+2]


# v -> v+3


omega_vp3_J_plus1 = zeros(shape=(v_max-3,J_max-2))
omega_vp3_J_minus1 = zeros(shape=(v_max-3,J_max-2))

for v7 in range(0,v_max-3):
    for J7 in range(0,J_max-3):
        
        omega_vp3_J_plus1[v7][J7] = Energy[v7+3][J7+2] - Energy[v7][J7+1]
        omega_vp3_J_minus1[v7][J7] = Energy[v7+3][J7+1] - Energy[v7][J7+2]


# Nu plottes spektret for alle disse overgange, med sandsynligheden for at være i starttilstanden
# brugt som amplitude for hver overgang.

new_weight = weight/40. # renormaliser aborsptionskoefficienterne så de kan passe ind i et normaliseret spektrum

# Lav en x-akse i bølgetal som kurven kan plottes på

wavenumber = arange(1500,6500,0.01)

spektrum = ones(len(wavenumber))

#print new_weight

# Lav lorentzkurver for hver eneste af overgangene og gang dem med det
# normaliserede kontinuum.


for v4 in range(0,v_max-1):
    for J4 in range(0,(J_max-2)):
        lorentz1 = new_weight[v4][J4+1]*((0.1)**2/(((wavenumber-omega_vp1_J_plus1[v4][J4])**2)+(0.1)**2))
        lorentz2 = new_weight[v4][J4+2]*((0.1)**2/(((wavenumber-omega_vp1_J_minus1[v4][J4])**2)+(0.1)**2))
        spektrum = spektrum * (1-lorentz1) * (1-lorentz2)

# Dette var kun for overgangene med v -> v+1. Nu kommer dem for v -> v+2

for v5 in range(0,v_max-2):
    for J5 in range(0,(J_max-2)):
        lorentz3 = 0.1*new_weight[v5][J5+1]*((0.1)**2/(((wavenumber-omega_vp2_J_plus1[v5][J5])**2)+(0.1)**2))
        lorentz4 = 0.1*new_weight[v5][J5+2]*((0.1)**2/(((wavenumber-omega_vp2_J_minus1[v5][J5])**2)+(0.1)**2))
        spektrum = spektrum * (1-lorentz3) * (1-lorentz4)

# v -> v+3


for v6 in range(0,v_max-3):
    for J6 in range(0,(J_max-2)):
        lorentz5 = 0.01*new_weight[v6][J6+1]*((0.1)**2/(((wavenumber-omega_vp3_J_plus1[v6][J6])**2)+(0.1)**2))
        lorentz6 = 0.01*new_weight[v6][J6+2]*((0.1)**2/(((wavenumber-omega_vp3_J_minus1[v6][J6])**2)+(0.1)**2))
        spektrum = spektrum * (1-lorentz5) * (1-lorentz6)


# Lav nu det fuldstændige spektrum ved at summe over alle lorentzkurverne i et loop

print omega_vp2_J_minus1

spektrumplot = figure(94)
plot(wavenumber,spektrum)
#xlim([1800,2320])
ylim([0.,1.1])

spektrumplot.show()




print "hej"