# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# ### Components in H<sub>2</sub> Behaviors Models

# <codecell>

from pylab import *

# <codecell>

class Gas:
    def __init__(self, therm, T = None, P = None, rho = None):
        '''class used to describe a gas
        
        Parameters
        ----------
        therm - a thermodynamic class that is used to relate state variables
        T - temperature (K)
        rho - density (kg/m^3)
        '''
        if T != None and P != None:
            self.T = T
            self.rho = therm.rho(T, P)
            self.P = P
        elif T != None and rho != None:
            self.T = T
            self.rho = rho
            self.P = therm.P(T, rho)
        elif P != None and rho != None:
            self.T = therm.T(P, rho)
            self.rho = rho
            self.P = P
        else:
            warnings.warn('system not properly defined')
            self.T, self.P, self.rho = T, P, rho
        self.therm = therm
    def update(self, T = None, P = None, rho = None):
        if T != None and P != None:
            self.T = T
            self.rho = self.therm.rho(T, P)
            self.P = P
        elif T != None and rho != None:
            self.T = T
            self.rho = rho
            self.P = self.therm.P(T, rho)
        elif P != None and rho != None:
            self.T = self.therm.T(P, rho)
            self.rho = rho
            self.P = P
        else:
            warnings.warn('system not properly defined')
            self.T, self.P, self.rho = T, P, rho

# <codecell>

class Orifice:
    def __init__(self, d, Cd = 1.):
        '''
        class used to describe a circular orifice
        
        future versions may be expanded to give effective area for other shapes
        
        Parameters
        ----------
        d - orifice diameter (m)
        Cd - discharge coefficient to account for non-plug flow (always <=1, assumed to be 1 for plug flow)
        
        Contains
        --------
        d - diameter (m)
        Cd- discharge coefficient 
        A - effective area (m^2)
        '''
        self.d, self.Cd, self.A = d, Cd, pi/4*d**2
    
    def mdot(self, rho, v):
        '''
        mass flow rate through orifice given gas density and velocity
        
        Parameters
        ----------
        rho - density (kg/m^3)
        v - velocity (m/s)
        
        Returns
        -------
        mdot - mass flow rate (kg/s)
        '''
        return rho*v*self.A*self.Cd

# <codecell>

class Source:
    '''used to describe a source (tank) that contains gas'''
    def __init__(self, V, gas):
        self.gas = gas
        self.V = V
        self.m = gas.rho*V
    def mdot(self, orifice, Ma = 1):
        '''returns the mass flow rate through an orifice, from the current tank conditions'''
        rho_throat = self.gas.therm.rho_Iflow(self.gas.therm.rho(self.gas.T, self.gas.P), Ma = Ma)
        T_throat = self.gas.therm.T_Iflow(self.gas.T, rho_throat, Ma = Ma)
        return orifice.mdot(rho_throat, self.gas.therm.a(T_throat, rho_throat)*Ma)
    def blowdown(self, t, orifice, Ma = 1):
        '''Returns the mass flow rate over time for a storage tank with an orifice.
        '''
        dt_array = t[1:] - t[:-1]
        mdot = array([self.mdot(orifice, Ma)])
        mtot = self.gas.therm.rho(self.gas.T, self.gas.P)*self.V
        for dt in dt_array:
            mtot -= self.mdot(orifice, Ma)*dt
            T = self.gas.therm.T_IE(self.gas.T, self.gas.therm.rho(self.gas.T, self.gas.P), mtot/self.V)
            P = self.gas.therm.P(self.gas.T, mtot/self.V)
            self.gas.update(T = T, P = P)
            mdot = append(mdot, self.mdot(orifice, Ma))
        return mdot

# <codecell>

class Enclosure:
    '''
    Enclosure used in the overpressure modeling
    '''
    def __init__(self, H, A, ceiling_vent, floor_vent):
        '''
        Describes the enclosure
        
        Parameters
        ----------
        H : encosure height
        A : area of floor and ceiling 
        ceiling_vent : vent class containing vent information for ceiling vent
        floor_vent : vent class containing vent information for floor vent
        '''
        self.H, self.A, self.ceiling_vent, self.floor_vent = H, A, ceiling_vent, floor_vent
        
class Vent:
    '''
    Vent used in overpressure modeling
    '''
    def __init__(self, A, H, Cd = 1, U_wind = 0):
        '''
        Describes the vent
        
        Parameters
        ----------
        A : vent cross-sectional area (m^2)
        H : vent height from floor (m)
        Cd: discharge coefficient of vent
        U_wind : wind velocity (m/s)
        '''
        self.A, self.Cd, self.U_wind = A, Cd, U_wind
        self.Qw = Cd*A*U_wind/sqrt(2) #See Lowesmith et al IJHE 2009

# <codecell>


