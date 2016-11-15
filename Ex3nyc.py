# -*- coding: utf-8 -*-


from numpy import *
from matplotlib.pyplot import *

# Definer variable:

rho_nul = 1
a = 1
G = 6.67*(10**(-11))
ds = 0.01
r_ydre = 1000
steps = 100000

r = arange(0,r_ydre,ds)
rho = rho_nul/(r*(1+r)**3)


# Her integrerer jeg densiteten numerisk for mit "datasæt" for at finde massen
# ud til en vis radius

M_num = empty(shape=(steps,1))
M_num[0] = 0
dM = zeros(shape=(steps,1))
for n1 in range(1,steps):
    dM[n1] = 4*pi*rho_nul*r[n1]/((1+r[n1])**3)*ds
    M_num[n1]  = M_num[n1-1] + dM[n1]

# Nu udregner jeg potentialet som funktion af radius. Jeg bruger at potentialet
# er integralet over massen (inden for readius du er ved) fra uendelig ind til
# den radius du er ved. Jeg skal derfor integrere massen numerisk, udefra, og
# bruge at phi(uendelig) = 0.

phi_num = empty(shape=(steps,1))
phi_num[(steps-1)] = 0 # Overvej at starte phi længere udefra, f.eks. phi[100000] = 0

dphi = zeros(shape=(steps,1))

for n3 in range((steps-2),-1,-1):
    dphi[n3] = G*ds*M_num[n3+1]/(r[n3+1]**2)
    phi_num[n3] = phi_num[n3+1] - dphi[n3]


# Her vil jeg plotte de karateristiske hastigheder for systemet, den cirkulære hastighed,
# undslippelseshastigheden og hastighedsdispersionen

r[0] = 10**-10

v_c = empty(shape=(steps,1))
for n4 in range(0,steps):
    v_c[n4] = sqrt(G*M_num[n4]/r[n4])

vcplot = figure(1)
plot(r,v_c)

vclogplot = figure(2)
plot(log10(r),log10(v_c))


v_esc = empty(shape=(steps,1))
for n5 in range(0,steps):
    v_esc[n5] = sqrt(2*abs(phi_num[n5]))

vescplot = figure(3)
plot(r,v_esc)

vesclogplot = figure(4)
plot(log10(r),log10(v_esc))

sigma = zeros(shape=(steps,1))
sigma_squared = zeros(shape=(steps,1))
dI = zeros(shape=(steps,1))
I = zeros(shape=(steps,1))
I[steps-1] = 0
for n6 in range(steps-2,-1,-1):
    dI[n6] = ds*rho[n6]*G*M_num[n6]/(r[n6]**2)
    I[n6] = I[n6+1] + dI[n6]
    sigma_squared[n6] = I[n6]/rho[n6]
    sigma[n6] = sqrt(sigma_squared[n6])

dvplot = figure(5)
plot(r,sigma)

dvlogplot = figure(6)
plot(log10(r),log10(sigma))

# Plot de tre hastigheder sammen

samletplot = figure(7)
plot(log10(r),log10(v_c),'b-',log10(r),log10(v_esc),'r-',log10(r),log10(sigma),'g-')
xlim([-2,3])


#vcplot.show()
#vclogplot.show()
#vescplot.show()
#vesclogplot.show()
#dvplot.show()
#dvlogplot.show()
samletplot.show()



# Plot dlog_rho_dlog_r og dlog_sigma_dlog_r, for at se om det er -3

#rho_log = log10(rho)
#sigma_squared_log = log10(sigma_squared)

#rho_logdlogr = zeros(shape=(steps,1))
#sigma_squared_logdlogr = zeros(shape=(steps,1))

#for n2 in range(1,(steps-1)):
#    rho_logdlogr[n2] = ((rho_log[n2+1]-rho_log[n2])/(log10(r[n2+1]/r[n2]))+(rho_log[n2]-rho_log[n2-1])/(log10(r[n2]/r[n2-1])))/2
#    sigma_squared_logdlogr[n2] = ((sigma_squared_log[n2+1]-sigma_squared_log[n2])/(log10(r[n2+1]/r[n2]))+(sigma_squared_log[n2]-sigma_squared_log[n2-1])/(log10(r[n2]/r[n2-1])))/2

#total_dlogdr = rho_logdlogr + sigma_squared_logdlogr

#logplot = figure(8)
#plot(log10(r[1:(steps-1)]),rho_logdlogr[1:(steps-1)],'r-',log10(r[1:(steps-1)]),sigma_squared_logdlogr[1:(steps-1)],'b-',log10(r[1:(steps-1)]),total_dlogdr[1:(steps-1)],'g-')
#xlim([-1.8,2.3])
#ylim([-10,2])

#logplot.show()

# Udregn den potentielle og den kinetiske energi for systemet:

W = empty(shape=(steps,1))
dW = empty(shape=(steps,1))
W[0] = 0

K = empty(shape=(steps,1))
dK = empty(shape=(steps,1))
K[0] = 0 

for nw in range(1,(steps-1)):
    dW[nw] = -4*pi*G*r[nw]*rho[nw]*M_num[nw]*ds
    W[nw] = W[nw-1] + dW[nw]
    dK[nw] = 6*pi*(r[nw]**2)*rho[nw]*(sigma[nw]**2)*ds
    K[nw] = K[nw-1] + dK[nw]

print K[steps-5]
print W[steps-5]

