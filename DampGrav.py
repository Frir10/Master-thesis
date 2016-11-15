import scipy
from math import pi

def H(a,x):
    ""
    P  = x**2
    H0 = scipy.exp(-x**2)
    Q  = 1.5/x**2
    return H0 - a / scipy.sqrt(pi) / P * ( H0 * H0 * (4. * P * P + 7. * P + 4. + Q) - Q - 1.0 )
   
def DampedGravity(r, Lambda0, N, b,f,gam,z):
    c  = 2.998e10  #;cm/s
    m_e = 9.1095e-28# ;g
    e  = 4.8032e-10# ;cgs units

    C_a  = scipy.sqrt(pi) * e**2 * f * Lambda0 * 1.e-8 / m_e / c / b
    a = Lambda0 * 1.e-8 * gam / (4.*pi*b)
    
    dl_D = b/c*Lambda0
    x = (Lambda/(z+1.0)-Lambda0)/dl_D+0.01
    
    tau  = C_a*N*H(a,x)
    return scipy.exp(-tau)