

from numpy import *
from matplotlib.pyplot import *
from matplotlib.mlab import *
from scipy.stats import *
from scipy.interpolate import interp1d


keV = 1.60217657/1.3806488*10**7  # Definition af kiloelektronvolt, da det er enheden for temperaturerne


x = linspace(0, 10, 10)
y = cos(-x**2/8.0)
f = interp1d(x, y)
f2 = interp1d(x, y, kind='cubic')



#xnew = linspace(0, 10, 400)

#plot(x,y,'o',xnew,f(xnew),'-', xnew, f2(xnew),'--')
#legend(['data', 'linear', 'cubic'], loc='best')
#show()


n_profiles = 29 # Antallet af filer

n_bins_density_dm = 100  # Antallet af datapunkter i density_dm filer


R_200 = [2416.453857,1646.021729,1780.135986,1609.971313,1152.909790,2323.536865,2315.112305,2397.056641,1042.733398,2275.391113,2152.105469,
        2584.523926,2431.409180,2499.355225,2474.897705,3131.108643,2507.641602,2140.491699,2339.896973,2502.261963,2364.993164,2573.030762,
        2292.840820,2304.907715,2094.564453,2410.086182,2436.733643,2612.871338,2406.322510]

n_print = 0

filelist_density_dm = ['densityDM_g0016649.dat','densityDM_g0052436.dat','densityDM_g0144846.dat','densityDM_g0163178.dat','densityDM_g0168361.dat',
                      'densityDM_g0272097.dat','densityDM_g1212639.dat','densityDM_g1483463.dat','densityDM_g1574117.dat','densityDM_g1657050.dat',
                      'densityDM_g1680241.dat','densityDM_g1987669.dat','densityDM_g2980844.dat','densityDM_g3327821.dat','densityDM_g3346905.dat',
                      'densityDM_g3888703.dat','densityDM_g4425770.dat','densityDM_g4606589.dat','densityDM_g4915399.dat','densityDM_g5265133.dat',
                      'densityDM_g5503149.dat','densityDM_g5699754.dat','densityDM_g6287794.dat','densityDM_g6348555.dat','densityDM_g6802296.dat',
                      'densityDM_g7263961.dat','densityDM_g7358274.dat','densityDM_g7570066.dat','densityDM_g7577931.dat']

filelist_density_gas = ['denG_g0016649.dat','denG_g0052436.dat','denG_g0144846.dat','denG_g0163178.dat','denG_g0168361.dat',
                       'denG_g0272097.dat','denG_g1212639.dat','denG_g1483463.dat','denG_g1574117.dat','denG_g1657050.dat',
                       'denG_g1680241.dat','denG_g1987669.dat','denG_g2980844.dat','denG_g3327821.dat','denG_g3346905.dat',
                       'denG_g3888703.dat','denG_g4425770.dat','denG_g4606589.dat','denG_g4915399.dat','denG_g5265133.dat',
                       'denG_g5503149.dat','denG_g5699754.dat','denG_g6287794.dat','denG_g6348555.dat','denG_g6802296.dat',
                       'denG_g7263961.dat','denG_g7358274.dat','denG_g7570066.dat','denG_g7577931.dat']

filelist_T_gas = ['temp_g0016649.dat','temp_g0052436.dat','temp_g0144846.dat','temp_g0163178.dat','temp_g0168361.dat',
                 'temp_g0272097.dat','temp_g1212639.dat','temp_g1483463.dat','temp_g1574117.dat','temp_g1657050.dat',
                 'temp_g1680241.dat','temp_g1987669.dat','temp_g2980844.dat','temp_g3327821.dat','temp_g3346905.dat',
                 'temp_g3888703.dat','temp_g4425770.dat','temp_g4606589.dat','temp_g4915399.dat','temp_g5265133.dat',
                 'temp_g5503149.dat','temp_g5699754.dat','temp_g6287794.dat','temp_g6348555.dat','temp_g6802296.dat',
                 'temp_g7263961.dat','temp_g7358274.dat','temp_g7570066.dat','temp_g7577931.dat']

filelist_vr_dm = ['vrDM_g0016649.dat','vrDM_g0052436.dat','vrDM_g0144846.dat','vrDM_g0163178.dat','vrDM_g0168361.dat',
                 'vrDM_g0272097.dat','vrDM_g1212639.dat','vrDM_g1483463.dat','vrDM_g1574117.dat','vrDM_g1657050.dat',
                 'vrDM_g1680241.dat','vrDM_g1987669.dat','vrDM_g2980844.dat','vrDM_g3327821.dat','vrDM_g3346905.dat',
                 'vrDM_g3888703.dat','vrDM_g4425770.dat','vrDM_g4606589.dat','vrDM_g4915399.dat','vrDM_g5265133.dat',
                 'vrDM_g5503149.dat','vrDM_g5699754.dat','vrDM_g6287794.dat','vrDM_g6348555.dat','vrDM_g6802296.dat',
                 'vrDM_g7263961.dat','vrDM_g7358274.dat','vrDM_g7570066.dat','vrDM_g7577931.dat']

filelist_vr_gas = ['vrGAS_g0016649.dat','vrGAS_g0052436.dat','vrGAS_g0144846.dat','vrGAS_g0163178.dat','vrGAS_g0168361.dat',
                  'vrGAS_g0272097.dat','vrGAS_g1212639.dat','vrGAS_g1483463.dat','vrGAS_g1574117.dat','vrGAS_g1657050.dat',
                  'vrGAS_g1680241.dat','vrGAS_g1987669.dat','vrGAS_g2980844.dat','vrGAS_g3327821.dat','vrGAS_g3346905.dat',
                  'vrGAS_g3888703.dat','vrGAS_g4425770.dat','vrGAS_g4606589.dat','vrGAS_g4915399.dat','vrGAS_g5265133.dat',
                  'vrGAS_g5503149.dat','vrGAS_g5699754.dat','vrGAS_g6287794.dat','vrGAS_g6348555.dat','vrGAS_g6802296.dat',
                  'vrGAS_g7263961.dat','vrGAS_g7358274.dat','vrGAS_g7570066.dat','vrGAS_g7577931.dat']

filelist_sigma_r = ['sigmarDM_g0016649.dat','sigmarDM_g0052436.dat','sigmarDM_g0144846.dat','sigmarDM_g0163178.dat','sigmarDM_g0168361.dat',
                 'sigmarDM_g0272097.dat','sigmarDM_g1212639.dat','sigmarDM_g1483463.dat','sigmarDM_g1574117.dat','sigmarDM_g1657050.dat',
                 'sigmarDM_g1680241.dat','sigmarDM_g1987669.dat','sigmarDM_g2980844.dat','sigmarDM_g3327821.dat','sigmarDM_g3346905.dat',
                 'sigmarDM_g3888703.dat','sigmarDM_g4425770.dat','sigmarDM_g4606589.dat','sigmarDM_g4915399.dat','sigmarDM_g5265133.dat',
                 'sigmarDM_g5503149.dat','sigmarDM_g5699754.dat','sigmarDM_g6287794.dat','sigmarDM_g6348555.dat','sigmarDM_g6802296.dat',
                 'sigmarDM_g7263961.dat','sigmarDM_g7358274.dat','sigmarDM_g7570066.dat','sigmarDM_g7577931.dat']

filelist_sigma_t = ['sigmatDM_g0016649.dat','sigmatDM_g0052436.dat','sigmatDM_g0144846.dat','sigmatDM_g0163178.dat','sigmatDM_g0168361.dat',
                 'sigmatDM_g0272097.dat','sigmatDM_g1212639.dat','sigmatDM_g1483463.dat','sigmatDM_g1574117.dat','sigmatDM_g1657050.dat',
                 'sigmatDM_g1680241.dat','sigmatDM_g1987669.dat','sigmatDM_g2980844.dat','sigmatDM_g3327821.dat','sigmatDM_g3346905.dat',
                 'sigmatDM_g3888703.dat','sigmatDM_g4425770.dat','sigmatDM_g4606589.dat','sigmatDM_g4915399.dat','sigmatDM_g5265133.dat',
                 'sigmatDM_g5503149.dat','sigmatDM_g5699754.dat','sigmatDM_g6287794.dat','sigmatDM_g6348555.dat','sigmatDM_g6802296.dat',
                 'sigmatDM_g7263961.dat','sigmatDM_g7358274.dat','sigmatDM_g7570066.dat','sigmatDM_g7577931.dat']

save_to_dm = ['prof_01_dm.txt','prof_02_dm.txt','prof_03_dm.txt','prof_04_dm.txt','prof_05_dm.txt','prof_06_dm.txt','prof_07_dm.txt',
             'prof_08_dm.txt','prof_09_dm.txt','prof_10_dm.txt','prof_11_dm.txt','prof_12_dm.txt','prof_13_dm.txt','prof_14_dm.txt',
             'prof_15_dm.txt','prof_16_dm.txt','prof_17_dm.txt','prof_18_dm.txt','prof_19_dm.txt','prof_20_dm.txt','prof_21_dm.txt',
             'prof_22_dm.txt','prof_23_dm.txt','prof_24_dm.txt','prof_25_dm.txt','prof_26_dm.txt','prof_27_dm.txt','prof_28_dm.txt',
             'prof_29_dm.txt']

save_to_gas = ['prof_01_gas.txt','prof_02_gas.txt','prof_03_gas.txt','prof_04_gas.txt','prof_05_gas.txt','prof_06_gas.txt','prof_07_gas.txt',
              'prof_08_gas.txt','prof_09_gas.txt','prof_10_gas.txt','prof_11_gas.txt','prof_12_gas.txt','prof_13_gas.txt','prof_14_gas.txt',
              'prof_15_gas.txt','prof_16_gas.txt','prof_17_gas.txt','prof_18_gas.txt','prof_19_gas.txt','prof_20_gas.txt','prof_21_gas.txt',
              'prof_22_gas.txt','prof_23_gas.txt','prof_24_gas.txt','prof_25_gas.txt','prof_26_gas.txt','prof_27_gas.txt','prof_28_gas.txt',
              'prof_29_gas.txt']


for nfile in range(0,n_profiles):
    
    data_density_dm = loadtxt(filelist_density_dm[nfile])
    radius1 = data_density_dm[:,0]
    density_dm = (data_density_dm[:,1])
    
    density_dm_interpol = interp1d(radius1, density_dm)
    
    
    data_density_gas = loadtxt(filelist_density_gas[nfile])
    radius2 = (data_density_gas[:,1]+data_density_gas[:,2])/2
    density_gas = data_density_gas[:,3]
    
    density_gas_interpol = interp1d(radius2,density_gas)
    
    
    data_T_gas = loadtxt(filelist_T_gas[nfile])
    T_gas = data_T_gas[:,5]*keV
    
    T_gas_interpol = interp1d(radius2,T_gas)
    
    data_sigma_r = loadtxt(filelist_sigma_r[nfile])
    radius3 = data_sigma_r[:,0]
    sigma_r = data_sigma_r[:,1]
    
    sigma_r_interpol = interp1d(radius3,sigma_r)
    
    data_sigma_t = loadtxt(filelist_sigma_t[nfile])
    radius4 = data_sigma_t[:,0]
    sigma_t = data_sigma_t[:,1]
    
    sigma_t_interpol = interp1d(radius4,sigma_t)
    
    
    data_vr_dm = loadtxt(filelist_vr_dm[nfile])
    radius5 = data_vr_dm[:,0]
    vr_dm = data_vr_dm[:,1]
    
    vr_dm_interpol = interp1d(radius5,vr_dm)
    
    
    data_vr_gas = loadtxt(filelist_vr_gas[nfile])
    radius6 = data_vr_gas[:,0]
    vr_gas = data_vr_gas[:,1]
    
    vr_gas_interpol = interp1d(radius6,vr_gas)
    
    output_dm = zeros(shape=(1000,5))
    output_gas = zeros(shape=(1000,4))
    
    radius_virials = linspace(0.,2*R_200[nfile],1000) + R_200[nfile]/2000.
    
    output_dm[:,0] = radius_virials
    
    output_gas[:,0] = radius_virials
    
    hej = 0
    
    for n in range(0,1000):
        
        #output_dm[n,0] = radius_virials[n]
        
        if radius1[0] < radius_virials[n] < radius1[-1]:
            
            output_dm[n,1] = density_dm_interpol(radius_virials[n])
        
        if radius2[0] < radius_virials[n] < radius2[-1]:
            
            output_gas[n,1] = density_gas_interpol(radius_virials[n])
            output_gas[n,2] = T_gas_interpol(radius_virials[n])
        
        if radius3[0] < radius_virials[n] < radius3[-1]:
            
            output_dm[n,2] = sigma_r_interpol(radius_virials[n])
        
        if radius4[0] < radius_virials[n] < radius4[-1]:
            
            output_dm[n,3] = sigma_t_interpol(radius_virials[n])
        
        if radius5[0] < radius_virials[n] < radius5[-1]:
            
            output_dm[n,4] = vr_dm_interpol(radius_virials[n])
        
        if radius6[0] < radius_virials[n] < radius6[-1]:
            
            output_gas[n,3] = vr_gas_interpol(radius_virials[n])
            
    
    savetxt(save_to_dm[nfile],output_dm)
    
    savetxt(save_to_gas[nfile],output_gas)
    
    n_print = n_print+1
    print n_print
    
    #interpol_density_dm = interpol_f(plotradius)
    
    #interpolplot = figure(nfile+1)
    
    #plot(radius,density,'*',plotradius,interpol_density_dm)
    
    #interpolplot.show()



