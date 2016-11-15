# -*- coding: utf-8 -*-
# Dette script plotter de forskellige teoretiske modeller for spektret sammen med det observerede spektrum.
# Den første model er et rent gæt, ud fra de egenskaber vi gætter på stjernen har, og den anden model er
# en forbedring af dette gæt, idet vi har ladet et program fitte vores gæt til det observerede spektrum
# og ladet det variere en del af inputparametrene.

from numpy import *
from matplotlib.pyplot import *


# Få fat i dataene opnået i SME

data_model1 = loadtxt('G0V.txt')
data_model2 = loadtxt('G0V3.txt')

wavelength = data_model1[:,0]
observed = data_model1[:,1]
model1 = data_model1[:,2]
model2 = data_model2[:,2]


# Plot nu de forskellige modeller og det observerede spektrum

modelplot = figure(38)
plot(wavelength,observed,'k',wavelength,model1,'b',wavelength,model2,'r')
modelplot.show()


# Nu plotter jeg spektrummer svarende til modelatmosfærer med forskellig tyngdeacceleration g,
# for at undersøge trykafhængigheden af to linier omkring 5890 Å

data_g3 = loadtxt('t5000g3z0.txt')
data_g4 = loadtxt('t5000g4z0.txt')
data_g5 = loadtxt('t5000g5z0.txt')

wavelength_pressure3 = data_g3[:,0]
model_g3 = data_g3[:,1]
wavelength_pressure4 = data_g4[:,0]
model_g4 = data_g4[:,1]
wavelength_pressure5 = data_g5[:,0]
model_g5 = data_g5[:,1]

print len(model_g5)

# Plot de 3 modeller, og se om der er nogen forskel i intensiteten af de relevante linier

pressureplot = figure(44)
plot(wavelength_pressure3,model_g3,'g',wavelength_pressure4,model_g4,'b',wavelength_pressure5,model_g5,'r')
#pressureplot.show()



