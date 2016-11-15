# -*- coding: utf-8 -*-
# Exercise 5


from numpy import *
from matplotlib.pyplot import *
from matplotlib.mlab import *
from scipy.stats import *

# Input:

M_sol = 1.99*(10**30)
G = 6.67*(10**-11)
m_p = 1.67*(10**-27)
mu = 0.59
Mpc = 3.09*10**22

d_center = 99*Mpc
RA_center = radians(180+360.*(59./(24.*60.)+35.7/(24.*60.*60.)))
dec_center = radians(27+57./60.+33./3600.)

data = loadtxt('coma.txt')

RA = radians(data[:,0])
dec = radians(data[:,1])
v_r = data[:,2]

# Man projicerer alle punkter ind på en overflade der er en kugleskal. Man kunne også projicere
# ind på en plan går igennem centrum af clusteren og er vinkelret på afstandsvektoren.

x_skal = d_center*cos(RA)*cos(dec)
y_skal = d_center*sin(RA)*cos(dec)
z_skal = d_center*sin(dec)

x0_skal = d_center*cos(RA_center)*cos(dec_center)
y0_skal = d_center*sin(RA_center)*cos(dec_center)
z0_skal = d_center*sin(dec_center)

r_kugleskal = sqrt((((x_skal-x0_skal)**2)+((y_skal-y0_skal)**2)+((z_skal-z0_skal)**2)))


# Projicer ind på planen i stedet for

x_plan = d_center*tan(RA-RA_center)
y_plan = d_center*tan(dec-dec_center)

r_plan = sqrt(((x_plan**2)+(y_plan**2)))

# Plot den radiale velocity som funktion af radius, op til r = 2.7 Mpc

R = r_kugleskal/Mpc
#R = r_plan/Mpc

velocityplot = figure(33)
plot(R,v_r,'.')
xlim([0,2.7])
#velocityplot.show()

# Lav bins i radius, i intervaller af størrelse 0.3 fra 0 til 2.7 og lav arrays med de
# tilsvarende hastigheder af galakserne.

number_bins = 9 # number of bins

bin_width = 2.7/number_bins
v_r_bins = zeros(shape=(number_bins,1000))# Array containing the line of sight velocities of each bin
n_galaxies = zeros(number_bins) # number of galaxies in each bin

mu_vector = zeros(number_bins)
sigma_vector = zeros(number_bins)
r_bin = zeros(number_bins)

for nbins in range(0,number_bins):
    
    R_bin = R[R<=((nbins+1)*bin_width)]
    v_r_binshej = v_r[R<=((nbins+1)*bin_width)]
    v_r_binshej = v_r_binshej[R_bin>(nbins*bin_width)]
        
    for nelement in range(0,len(v_r_binshej)):
        v_r_bins[nbins][nelement] = v_r_binshej[nelement]
    
    n_galaxies[nbins] = len(v_r_binshej)    
    
    # Få middelværdi og spredning for hver v_r_binshej
    
    (mu_vector[nbins], sigma_vector[nbins]) = norm.fit(v_r_binshej)
    r_bin[nbins] = bin_width*nbins+bin_width*0.5 # r svarende til din bin vi er i

# Plot en af binsne med det gaussiske fit ovenpå    
    
bin_of_interest = 1 # which bin to plot (must be smaller than number_bins

hist_n,hist_bins,hist_patches = hist(v_r_bins[bin_of_interest][0:n_galaxies[bin_of_interest]],30) # hist skaber histogrammet,
# det man kalder dem refererer bare til de objekter kommandoen skaber, jeg har det kun med fordi jeg skal kalde
# dem senere (lige nu faktisk).
gauss_curve = normpdf(hist_bins,mu_vector[bin_of_interest],sigma_vector[bin_of_interest])*(n_galaxies[bin_of_interest]*200)
#hej = figure(4)
velocidistplot = figure(55)
axdisp = velocidistplot.add_subplot(111)
axdisp.plot(hist_bins,gauss_curve,'r-')
axdisp.hist(v_r_bins[bin_of_interest][0:n_galaxies[bin_of_interest]],30)

#velocidistplot.show() 



# Plot hastighedsdispersionen i hver bin som funktion af radius

dispplot = figure(23)
plot(r_bin,sigma_vector,'.')

dispplot.show()

print len(R[R<=2.7])
