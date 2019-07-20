import numpy as np
from scipy.integrate import quad

class constants:
    """Useful constants"""
    speed_of_light = 3e5 # c in km/s
    c_by_H0 = 0.01*speed_of_light #present Hubble radius in Mpc/h

class cosmology(constants):
    """Useful functions for cosmology"""
    def __init__(self,omg_r0,omg_m0,omg_l0):
        self.omg_r0 = omg_r0
        self.omg_m0 = omg_m0
        self.omg_l0 = omg_l0
        self.omg_k0 = 1. - (omg_r0 +omg_m0 + omg_l0)
        print "Intializing cosmology with omega parameters \n omg_r0 = %2.5f\n omg_m0 = %2.5f\n omg_l0 = %2.5f\n omg_k0 = %2.5f"%(self.omg_r0,self.omg_m0,self.omg_l0,self.omg_k0)

    def scale2redshift(self,scale):
        """scale to redshift conversion"""
        redshift = 1.0/(scale) - 1.0
        return redshift

    def redshift2scale(self,redshift):
        """redshift to scale conversion"""
        scale = 1.0/(1.0 + redshift)
        return scale

    def dHub0(self,z):
        """Distance according to Hubble law in units of Mpc/h"""
        dist = self.c_by_H0 * z
        return dist

    def H_by_H0(self,redshift):
        """Friedmann equation H/H0 at a particular redshift"""
        h = self.omg_r0*(1+redshift)**4 + self.omg_m0*(1+redshift)**3 + self.omg_l0 + self.omg_k0*(1+redshift)**2
        h = np.sqrt(h)
        return h

    def redshift2metric(self,z):
        """Gives metric distance at a given redshift"""
        value = quad((lambda z : 1./self.H_by_H0(z)), 0, z)[0]
        value = value * constants.c_by_H0
        k = np.abs(-1.*self.omg_k0/(self.c_by_H0)**2)
        if(self.omg_k0 > 0.0):
            value = np.sinh(np.sqrt(k)*value)/np.sqrt(k)
        elif (self.omg_k0 <  0.0):
            value = np.sin(np.sqrt(k)*value)/np.sqrt(k)
        else:
            value = value
        return value

    def z2d_ang(self,z):
        """Redshift to angular diameter distance (Mpc/h)"""
        return self.redshift2metric(z)/(1.0 + z)

    def z2d_lum(self,z):
        """Redshift to luminosity distance (Mpc/h)"""
        return self.redshift2metric(z)*(1.0 + z)




