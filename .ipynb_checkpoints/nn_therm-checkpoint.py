# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>
from __future__ import print_function, absolute_import, division
import numpy as np
from scipy import optimize, interpolate, integrate
from scipy import constants as const
import pandas as pd
from CoolProp import CoolProp
# <markdowncell>

# ##AbelNoble: 
# 
# - thermodynamic class that Abel-Nobel EOS for various necessary state variables
# - default $b$, $\gamma$, and $M$ set for H<sub>2</sub>, and default ideal gas constant (don't know why anyone would ever want to change this--maybe to check different precisions??)
# - ideal gas law by setting $b = 0$
# - I can envision writing other classes with the same functions (i.e. something that calls refprop to find the same variables)
# 
class EOS:
    # def __init__(self, species = ['hydrogen'], x = [1.]):
        # '''
        # Class that uses CoolProp for equation of state calculations.
        # '''
        # from CoolProp import CoolProp
        # self._cp = CoolProp
        # if len(species) == 1:
            # self.spec = species[0]
        # else:
            # if size(x) != size(species):
                # warnings.warn('mole fractions and species lists not the same length')
            # s = 'REFPROP-MIX:'
            # for g, xspec in zip(species, x):
                # s += g + '[' + str(xspec) + ']' + '&'
            # self.spec = s[:-1]
    def __init__(self, species="hydrogen"):
        
        from CoolProp import CoolProp
        self._cp = CoolProp
        self.spec = species
        #self.MW = self._cp.PropsSI(self.spec, 'molemass') 
    def h(self, T, P):
        '''enthalpy (J/kg) of a gas at temperature T (K) and pressure P (Pa)'''
        return self._cp.PropsSI('H', 'T', T, 'P', P, self.spec)
    def P(self, T, rho):
        '''
        returns the pressure given the temperature and density
        
        Parameters
        ----------
        T - temperature (K)
        rho - density (kg/m^3)
        
        Returns
        -------
        P - pressure (Pa)
        '''
        return self._cp.PropsSI('P', 'T', T, 'D', rho, self.spec)

    def T(self, P, rho):
        '''
        returns the temperature given the temperature and pressure
        
        Parameters
        ----------
        P - pressure (Pa)
        rho - density (kg/m^3)
        
        Returns
        -------
        T - temperature (K)
        '''
        return self._cp.PropsSI('T', 'D', rho, 'P', P, self.spec)
    
    def rho(self, T, P):
        '''
        returns the pressure given the temperature and pressure
        
        Parameters
        ----------
        T - temperature (K)
        P - pressure (Pa)
        
        Returns
        -------
        rho - density (kg/m^3)
        '''
        return self._cp.PropsSI('D', 'T', T, 'P', P, self.spec)
 
    def a(self, T, rho):
        '''
        returns the speed of sound given the temperature and density
        
        Parameters
        ----------
        T - temperature (K)
        P - pressure (Pa)
        
        Returns
        -------
        a - speed of sound (m/s)
        '''
        return self._cp.PropsSI('A', 'T', T, 'D', rho, self.spec)
            
    def errS(self, T0, rho0, T1, rho1):
        '''
        returns the error in entropy (J/kg-K) between 2 states specified by the 
        temperatures and densities.
        '''
        s0 = self._cp.PropsSI('S', 'T', T0, 'D', rho0, self.spec)
        s1 = self._cp.PropsSI('S', 'T', T1, 'D', rho1, self.spec)
        return s1 / s0 - 1.
   
   
    def rho_Iflow(self, rho0, Ma = 1, Ma0 = 0):
        '''
        solves for the density of gas after isentropic flow
        
        Parameters
        ----------
        rho0 - density before exansion (kg/m^3)
        Ma - mach number during epansion 
        
        Returns
        -------
        rho - density after expansion (kg/m^3)
        '''
        from scipy import optimize
        T0 = 273.
        sf = self._cp.State(self.spec, {'T':T0, 'D':rho0})#self._cp
        def errT(T0):
            s0 = self._cp.State(self.spec, {'T':T0, 'D':rho0})
            def errh(h):
                sf.update({'S':s0.get_s(), 'H':h})
                hf = s0.get_h()*1e3 + (Ma0*s0.get_speed_sound())**2/2 - (Ma*sf.get_speed_sound())**2/2.
                return hf/1000. - h
            hf = optimize.root(errh, s0.get_h())
            Tf = self.T_IflowV(T0, 0, s0.p*1e3, Ma*sf.get_speed_sound(), sf.p*1e3)
            return Tf - sf.T
        T0 = optimize.root(errT, T0)
        return sf.get_rho()
       
    def T_Iflow(self, T0, rho, Ma = 1, Ma0 = 0):
        '''
        solves for the temperature of gas after an isentropic flow
        
        Parameters
        ----------
        T0 - temperature before expansion (K)
        rho - density after expansion (kg/m^3)
        Ma - mach number during expansion
        Ma0 - mach number before expansion
        
        Returns
        -------
        T - temperature after expansion(K)
        '''
        from scipy import optimize
        def errS(x):
            [P0, Pf] = x
            P0, Pf = float(P0), float(Pf)
            s0 = self._cp.State(self.spec, {'P':P0/1e3, 'T':T0})
            sf = self._cp.State(self.spec, {'P':Pf/1e3, 'D':rho})
            hf = s0.get_h()*1e3 + (Ma0*s0.get_speed_sound())**2/2. - (Ma*sf.get_speed_sound())**2/2.
            return array([sf.get_s() - s0.get_s(), sf.get_h()*1e3 - hf])
        P0, Pf = optimize.fsolve(errS, array([101325., 101325.]))
        return self._cp.PropsSI('T', 'P', Pf, 'D', rho, self.spec)
    
    def P_Iflow(self, P0, rho, Ma = 1, Ma0 = 0):
        '''
        solves for the temperature of gas after an isentropic flow
        
        Parameters
        ----------
        P0 - Pressure before expansion (Pa)
        rho - density after expansion (kg/m^3)
        Ma - mach number during expansion
        Ma0 - mach number before expansion
        
        Returns
        -------
        P - pressure after expansion(Pa)
        '''
        from scipy import optimize
        def errS(x):
            [T0, Pf] = x
            T0, Pf = float(T0), float(Pf)
            s0 = self._cp.State(self.spec, {'P':P0/1e3, 'T':T0})
            sf = self._cp.State(self.spec, {'P':Pf/1e3, 'D':rho})
            hf = s0.get_h()*1e3 + (Ma0*s0.get_speed_sound())**2/2. - (Ma*sf.get_speed_sound())**2/2.
            return array([sf.get_s() - s0.get_s(), sf.get_h()*1e3 - hf])
        T0, Pf = optimize.fsolve(errS, array([273., 101325.]))
        return Pf
    
    def T_IflowV(self, T0, V0, P0, V, P):
        '''
        solves for the temperature of gas after an isentropic flow
        
        Parameters
        ----------
        T0 - temperature before expansion (K)
        V0 - velocity before expansion (m/s)
        P0 - pressure before expansion (Pa)
        V - velocity after expansion (m/s)
        P - pressure after expansion (Pa)
        
        Returns
        -------
        T - temperature after expansion(K)
        '''
        h = self._cp.PropsSI('H', 'T', T0, 'P', P0, self.spec) + V0**2/2. - V**2/2.
        return self._cp.PropsSI('T', 'P', P, 'H', h, self.spec)
    
    def T_IE(self, T0, rho0, rho):
        from scipy import optimize
        '''
        solves for the temperature of gas after an isentropic expansion, given an old and new 
        density
        
        Parameters
        ----------
        T0 - temperature before expansion (K)
        rho0 - density before expansion (kg/m^3)
        rho - density after expansion (kg/m^3)
        
        Returns
        -------
        T - temperature after expansion(K)
        '''
        S0 = self._cp.PropsSI('S', 'T', T0, 'D', rho0, self.spec)
        def errD(T):
            return rho - self._cp.PropsSI('D', 'T', T, 'S', S0, self.spec)
        return optimize.newton(errD, T0)
    

    
    def T_IflowSonic(self, T0, P0, P):
        from scipy import optimize
        '''
        solves for the temperature of gas after an isentropic flow, assuming that the entrance
        and exit velocities are sonic.
        
        Parameters
        ----------
        T0 - temperature before expansion (K)
        P0 - pressure before expansion (Pa)
        P - pressure after expansion (Pa)
        
        Returns
        -------
        T - temperature after expansion(K)
        '''
        V0 = self.a(T0, self.rho(T0, P0))
        def errH(T):
            V = self.a(T, self.rho(T, P))
            h = self._cp.PropsSI('H', 'T', T0, 'P', P0, self.spec) + V0**2/2. - V**2/2.
            return h - self._cp.PropsSI('H', 'T', T, 'P', P, self.spec)
        T = optimize.newton(errH, T0)
        return T
    
    def T_AdFlow(self, T0, P0, Ma0, P, Ma):
        from scipy import optimize
        '''
        solves for the temperature of gas after an isentropic flow, assuming that the entrance
        and exit velocities are sonic.
        
        Parameters
        ----------
        T0 - temperature before expansion (K)
        P0 - pressure before expansion (Pa)
        P - pressure after expansion (Pa)
        
        Returns
        -------
        T - temperature after expansion(K)
        '''
        V0 = Ma0*self.a(T0, self.rho(T0, P0))
        def errH(T):
            V = Ma*self.a(T, self.rho(T, P))
            h = self._cp.PropsSI('H', 'T', T0, 'P', P0, self.spec) + V0**2/2. - V**2/2.
            return h - self._cp.PropsSI('H', 'T', T, 'P', P, self.spec)
        T = optimize.newton(errH, T0)
        return T
    # def __init__(self, species="hydrogen"):
        
        # from CoolProp import CoolProp
        # self._cp = CoolProp
        # self.spec = species
        # self.MW = self._cp.PropsSI(self.spec, 'molemass') 
        
    # def P(self, T, rho):
        # '''
        # returns the temperature given the pressure and density (and sets the phase)
        
        # Parameters
        # ----------
        # T: float
            # temperature (K)
        # rho: float
            # density (kg/m^3)
        
        # Returns
        # -------
        # P: float
            # pressure (Pa)
        # '''
        # P, Q = self._cp.PropsSI(['P', 'Q'], 'D', rho, 'T', T, self.spec)
        # self.phase = self._cp.PhaseSI('D', rho, 'T', T, self.spec)
        # return P
    
    # def T(self, P, rho):
        # '''
        # returns the temperature given the pressure and density
        
        # Parameters
        # ----------
        # P: float
            # pressure (Pa)
        # rho: float
            # density (kg/m^3)
        
        # Returns
        # -------
        # T: float
            # temperature (K)
        # '''
        # T, Q = self._cp.PropsSI(['T', 'Q'], 'D', rho, 'P', P, self.spec)
        # self.phase = self._cp.PhaseSI('D', rho, 'P', P, self.spec)
        # return T
    
    # def rho(self, T, P):
        # '''
        # returns the denstiy given the temperature and pressure - if at saturation conditions, 
        # requires phase to already be set
        
        # Parameters
        # ----------
        # T: float
            # temperature (K)
        # P: flaot
            # pressure (Pa)
        
        # Returns
        # -------
        # rho:
            # density (kg/m^3)
        # '''
        # try:
            # rho, Q = self._cp.PropsSI(['D', 'Q'], 'T', T, 'P', P, self.spec)
            # self.phase = self._cp.PhaseSI('T', T, 'P', P, self.spec)
            # return rho
        # except:
            # try:
                # warnings.warn('using %s phase information to calculate density'%self.phase)
            # except:
                # warnings.warn('assuming gas phase to calculate density')
                # self.phase = 'gas'
            # return self._cp.PropsSI('D', 'T|'+self.phase, T, 'P', P, self.spec)
            
    # def rho_P(self, T, phase):
        # '''
        # returns the density and pressure given the temperature and phase
        
        # Parameters
        # ----------
        # T: float
            # temperautre (K)
        # phase: string
            # 'gas' or 'liquid'
        
        # Returns
        # -------
        # (rho, P): tuple of floats 
            # rho - density (kg/m^3)
            # P - pressure (Pa)
        # '''
        # rho, P = self._cp.PropsSI(['D', 'P'], 'T', T, 'Q', {'gas':1, 'liquid':0}[phase], self.spec)
        # self.phase = phase
        # return rho, P
        
    # def rho_T(self, P, phase):
        # '''
        # returns the density and temperature given the pressure and phase
        
        # Parameters
        # ----------
        # P: float
            # temperautre (K)
        # phase: string
            # 'gas' or 'liquid'
        
        # Returns
        # -------
        # (rho, P): tuple of floats 
            # rho - density (kg/m^3)
            # P - pressure (Pa)
        # '''
        # rho, T = self._cp.PropsSI(['D', 'T'], 'P', P, 'Q', {'gas':1, 'liquid':0}[phase], self.spec)
        # self.phase = phase
        # return rho, T
        
    # def P_T(self, rho, phase):
        # '''
        # returns the pressure and temperature given the density and phase
        
        # Parameters
        # ----------
        # T: float
            # temperautre (K)
        # phase: string
            # 'gas' or 'liquid'
        
        # Returns
        # -------
        # (rho, P): tuple of floats 
            # rho - density (kg/m^3)
            # P - pressure (Pa)
        # '''
        # P, T = self._cp.PropsSI(['P', 'T'], 'D', rho, 'Q', {'gas':1, 'liquid':0}[phase], self.spec)
        # self.phase = phase
        # return P, T
        
    # def _err_H(self, T1, P1, v1, T2, P2, v2, usePhase1 = False, usePhase2 = False):
        # '''
        # error in total enthalpy (J/kg) for a gas at two different states and velocities
        
        # Parameters
        # ----------
        # T1: float
            # tempearture (K) at state 1
        # P1: float
            # pressure (Pa) at state 1
        # v1: float
            # velocity (m/s) at state 1
        # T2: float
            # tempearture (K) at state 2
        # P2: float
            # pressure (Pa) at state 2
        # v2: float
            # velocity (m/s) at state 2
        
        # Returns
        # -------
        # err_h: float
            # error in enthalpy (J/kg)
        # '''
        # try:
            # h1 = self._cp.PropsSI('H', 'T', T1, 'P', P1, self.spec)
        # except:
            # print('tp')
            # h1 = self._cp.PropsSI('H', 'T|'+self.phase, T1, 'P', P1, self.spec)
        # try:
            # h2 = self._cp.PropsSI('H', 'T', T2, 'P', P2, self.spec)
        # except:
            # print('tp')
            # h2 = self._cp.PropsSI('H', 'T|'+self.phase, T2, 'P', P2, self.spec)
        # return h1 + v1**2/2. - (h2 + v2**2/2.)
        
    # def _err_H_P_rho(self, P1, rho1, v1, P2, rho2, v2):
        # '''
        # error in total enthalpy (J/kg) for a gas at two different states and velocities
        
        # Parameters
        # ----------
        # T1: float
            # tempearture (K) at state 1
        # P1: float
            # pressure (Pa) at state 1
        # v1: float
            # velocity (m/s) at state 1
        # T2: float
            # tempearture (K) at state 2
        # P2: float
            # pressure (Pa) at state 2
        # v2: float
            # velocity (m/s) at state 2
        
        # Returns
        # -------
        # err_h: float
            # error in enthalpy (J/kg)
        # '''
        # try:
            # h1 = self._cp.PropsSI('H', 'P', P1, 'D', rho1, self.spec)
        # except:
            # print('tp')
            # h1 = self._cp.PropsSI('H', 'P|'+self.phase, P1, 'D', rho1, self.spec)
        # try:
            # h2 = self._cp.PropsSI('H', 'P', P2, 'D', rho2, self.spec)
        # except:
            # print('tp')
            # h2 = self._cp.PropsSI('H', 'P|'+self.phase, P2, 'D', rho2, self.spec)
        # return h1 + v1**2/2. - (h2 + v2**2/2.)
    
    # def _err_S(self, T1, P1, T2, P2):
        # '''
        # returns the difference in entropy (J/kg) between 2 states specified by the 
        # temperatures and pressures.
        
        # Parameters
        # ----------
        # T1: float
            # temperature of gas at point 1 (K)
        # P1: float
            # pressure of gas at point 1 (Pa)
        # T2: float
            # temperature of gas at point 2 (K)
        # P2: float
            # Pressure of gas at point 2 (Pa)
            
        # Returns
        # -------
        # err_S: float
            # error in enthalpy between the two different states (J/kg)
        # '''
        # try:
            # s1 = self._cp.PropsSI('S', 'T', T1, 'P', P1, self.spec)
        # except:
            # print('tp')
            # s1 = self._cp.PropsSI('S', 'T|'+self.phase, T1, 'P', P1, self.spec)
        # try:
            # s2 = self._cp.PropsSI('S', 'T', T2, 'P', P2, self.spec)
        # except:
            # print('tp')
            # S2 = self._cp.PropsSI('S', 'T|'+self.phase, T2, 'P', P2, self.spec)
        # return s1 - s2
        
    # def _X(self, Y, other = 'air'):
        # MW = self._cp.PropsSI('M', self.spec)
        # MW_other = self._cp.PropsSI('M', other)
        # return Y/MW/(Y/MW + (1-Y)/MW_other)
    
    # # def _err_adiabatic_out(self, T1, P1, T2, P2, dm_per_m):
        # # '''
        # # returns the difference in internal energy (J/kg) between 2 states specified by the 
        # # temperatures and pressures.
        
        # # Parameters
        # # ----------
        # # T1: float
            # # temperature of gas at point 1 (K)
        # # P1: float
            # # pressure of gas at point 1 (Pa)
        # # T2: float
            # # temperature of gas at point 2 (K)
        # # P2: float
            # # Pressure of gas at point 2 (Pa)
        # # dm_per_m: float
            # # relative mass difference between state 2 and 1 (m2 - m1)/m1
        # # Returns
        # # -------
        # # err_adiabatic: float
            # # error in energy balance between the two different states
        # # '''
        # # U1, H1, rho1 = self._cp.PropsSI(['U', 'H', 'D'], 'T', T1, 'P', P1, self.spec)
        # # U2, rho2 = self._cp.PropsSI(['U', 'D'], 'T', T2, 'P', P2, self.spec)
        # # return (rho2*U2 - rho1*U1) / (dm_per_m*rho1*H1)-1
            
    # def a(self, T, P, S = None):
        # '''
        # returns the speed of sound given the temperature and pressure, or temperature and entropy
        
        # Parameters
        # ----------
        # T: float
            # temperature (K)
        # P: float
            # Pressure (Pa)
        # S: float
            # Entropy (J/K)
        
        # Returns
        # -------
        # a: float
            # speed of sound (m/s)
        # '''
        # try:
            # a = self._cp.PropsSI('A', 'T', T, 'P', P, self.spec)
            # if not np.isfinite(a):
                # P_pterb = 0.1
                # rho = upstream_fluid.therm._cp.PropsSI('D', 'P|'+self.phase,  [P, P + P_pterb], 'T', T, self.spec)
                # a = np.sqrt(P_pterb/np.gradient(rho))[0]
            # return a
        # except:
            # csg, rhog = self._cp.PropsSI(['A', 'D'], 'T|gas', T, 'P', P, self.spec)
            # csl, rhol = self._cp.PropsSI(['A', 'D'], 'T|liquid', T, 'P', P, self.spec)
            # rho, quality = self._cp.PropsSI(['D', 'Q'], 'T', T, 'S', S, self.spec)
             # # calculate gas/liquid void fractions
            # void_frac_gas = rho * quality / rhog
            # void_frac_liquid = 1.0 - void_frac_gas

            # # calculate mixture sound speed
            # term = np.sqrt(rhog * csg ** 2 / (void_frac_liquid * rhog * csg ** 2 + void_frac_gas * rhol * csl ** 2))
            # cs = csg * csl * term / (void_frac_liquid * csg + void_frac_gas * csl * term)
            # return cs
        
    # def h(self, T = None, P = None, rho = None, phase = None):
        # '''
        # enthalpy (J/kg) of a gas at temperature T (K) and pressure P (Pa)
        
        # Parameters
        # ----------
        # T: float
            # tempearture (K)
        # P: float
            # pressure (Pa)
        
        # Returns
        # -------
        # h: float
            # heat capacity (J/kg)
        # '''
        # if phase is None:
            # str = '|not_imposed'
        # else:
            # str = '|' + phase        
        # if T is not None and P is not None:
            # return self._cp.PropsSI('H', 'T' + str, T, 'P', P, self.spec)
        # elif T is not None and rho is not None:
            # return self._cp.PropsSI('H', 'T', T, 'D', rho, self.spec)
        # elif P is not None and rho is not None:
            # return self._cp.PropsSI('H', 'D', rho, 'P', P, self.spec)
        # else:
            # raise warnings.warn('system not properly defined')
            
    # def s(self, T = None, P = None, rho = None, phase = None):
        # '''
        # entropy (J/kg-K) of a fluid at temperature T (K) and pressure P (Pa)
        
        # Parameters
        # ----------
        # T: float
            # tempearture (K)
        # P: float
            # pressure (Pa)
        
        # Returns
        # -------
        # h: float
            # heat capacity (J/kg)
        # '''
        # if phase is None:
            # str = '|not_imposed'
        # else:
            # str = '|' + phase
        # if T is not None and P is not None:
            # return self._cp.PropsSI('S', 'T' + str, T, 'P', P, self.spec)
        # elif T is not None and rho is not None:
            # return self._cp.PropsSI('S', 'T' + str, T, 'D', rho, self.spec)
        # elif P is not None and rho is not None:
            # return self._cp.PropsSI('S', 'D' + str, rho, 'P', P, self.spec)
        # else:
            # raise warnings.warn('system not properly defined')
# # <codecell>



