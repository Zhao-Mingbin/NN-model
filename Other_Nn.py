# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 10:18:06 2014

@author: xueli
"""

#%%
#from pylab import *
from scipy import optimize, interpolate, integrate
from scipy.constants import R
from math import *


from nn_comps import Gas, Orifice, Source
from nn_therm import EOS

# <codecell>

class NotionalNozzle:
    def __init__(self, high_P_gas, orifice, low_P_gas):
        '''Notional nozzle class'''
        self.high_P_gas, self.orifice, self.low_P_gas = high_P_gas, orifice, low_P_gas
        
    def throat(self):
        '''Determine the conditions at the throat'''
        high_P_gas, orifice = self.high_P_gas, self.orifice
        rhot    = high_P_gas.therm.rho_Iflow(high_P_gas.rho) # throat density [kg/m**3]
        Tt      = high_P_gas.therm.T_IE(high_P_gas.T, high_P_gas.rho, rhot) # temperature at throat [K]
        ct      = high_P_gas.therm.a(Tt, rhot) # throat speed of sound [m/s]
        mdot    = orifice.mdot(rhot, ct); # Mass flow rate [kg/s]
        pt      = high_P_gas.therm.P(Tt, rhot) # throat pressure [Pa]
        return (rhot, Tt, ct, mdot, pt)
    
    def calculate(self, model):
        '''
        Calculates the properties after the notional nozzle using one of the following models:
        
        YuceilOtugen
        EwanMoodie
        Birch
        Birch2
        Molkov
        HarstadBellan
        
        Parameters
        ----------
        model - string of model name (default is YuceilOtugen)
        
        Returns
        -------
        tuple of (Gas object, Orifice object, velocity), all at exit of notional nozzle
        '''
        if model == 'Yuceil 2002模型':#YuceilOtugen
            return self._YuceilOtugen()
        elif model == 'Ewan 1986模型':
            return self._EwanMoodie()
        elif model == 'Brich 1984模型':
            return self._Birch()
        elif model == 'Brich 1987模型':
            return self._Birch2()
        elif model == 'Molkov 模型':
            return self._Molkov()
        elif model is 'HarstadBellan':
            return self._HarstadBellan()
        else:
            print ('Error - model by name of ' + model + ' not found.')
            print ('Possible models:')
            return (nan, nan, nan)
        
    
    def _YuceilOtugen(self):
        '''
        #*************************Ma**********************************************
        # This function uses the relationships developed by YUCEIL and OTUGEN 
        # (Physics of Fluids Vol 14 (2002) 4206 - 4215) to compute the effective 
        # area, velocity, and temperature after the Mach disk for an underexpanded 
        # jet. The model is identical to the Birch2 model except that the flow 
        # downstream of the Mach disk is assuemed to be sonic, and the downstream 
        # temperature is assumed to be the throat temperture.
        #
        # An additional feature that has been added to the model is the
        # inclusion of the Abel-Noble equation of state so that non-ideal gas
        # characteristics for HYDROGEN can be captured.
        #
        # Unlike the first Birch and the Ewan and Moodie models, this model 
        # accounts for the conservation of MASS and MOMENTUM; however, the
        # model assumes the flow remains supersonic, and no attempt is made to 
        # account for the change in entropy across the Mach disk. 
        #
        # Created by Isaac Ekoto (March 10, 2011)
        #*************************************************************************
        '''

        # Determine the conditions at the throat
        (rhot, Tt, ct, mdot, pt) = self.throat()
        
        # Simultaneous solution of conservation of mass and momentum (see Birch paper):
        Veff = ct + (pt - self.low_P_gas.P)/(ct*rhot) # Effective Velocity [m/s]
        T2 = self.high_P_gas.therm.T_IflowV(Tt, ct, pt, Veff, self.low_P_gas.P) # Effective Temperature (K)
        rho2 = self.high_P_gas.therm.rho(T2, self.low_P_gas.P)
        gas_eff = Gas(self.high_P_gas.therm, T = T2, rho = rho2)
        
        # Conservation of mass yeilds orifice diameter:
        orifice_eff = Orifice(sqrt(mdot/(rho2*Veff)*4/pi))
        
        return (gas_eff, orifice_eff, Veff)

    def _EwanMoodie(self):
        '''
        #*************************************************************************
        # This function uses the relationships developed by Ewan and Moodie
        # (Combust Sci Tech Vol 45 (1986) 275 - 288) to compute the effective 
        # area, velocity, and temperature after the Mach disk for an underexpanded 
        # jet. 
        # 
        # An additional feature that has been added to the model is the
        # inclusion of the Abel-Noble equation of state so that non-ideal gas
        # characteristics for HYDROGEN can be captured.
        #
        # The model accounts for the conservation of MASS only. The model is 
        # identical to the Birch Model except that the exit temperature is assumed
        # to be the throat temperature instead of the ambient temperature. The 
        # model makes the following assumptions:
        #
        # 1) The effective temperature after the Mach disk is equal to the
        #    throat temperature
        # 2) The Mach number after the Mach disk is assumed to be 1
        # 3) The static pressure at the Mach disk is assumed to be ambient pressure
        # 
        # No attempt is made to account for the change in entropy across the shock
        #
        # Created by Isaac Ekoto (March 10, 2011)
        #*************************************************************************
        '''
        (rhot, Tt, ct, mdot, pt) = self.throat()
        
        # Calculate state at expansion conditions from assumptions 1 & 3:
        T2      = Tt;
        p2      = self.low_P_gas.P;
        rho2    = self.high_P_gas.therm.rho(T2, p2);
        gas_eff = Gas(self.high_P_gas.therm, T = T2, rho = rho2)

        # The Mach number after the Mach disk is assumed to be 1 (this clearly is 
        # a faulty assumption).
        Veff    = self.high_P_gas.therm.a(T2, rho2) # speed of sound at Mach disk [m/s]
        
        # Conservation of mass yeilds orifice diameter:
        orifice_eff = Orifice(sqrt(mdot/(rho2*Veff)*4/pi))
        
        return (gas_eff, orifice_eff, Veff)
    
    def _Birch(self):
        '''
        #*************************************************************************
        # This function uses the relationships developed by Birch et al. 
        # (Combust Sci Tech Vol 36 (1984) 249 - 261) to compute the effective 
        # area, velocity, and temperature after the Mach disk for an underexpanded 
        # jet. 
        # 
        # An additional feature that has been added to the model is the
        # inclusion of the Abel-Noble equation of state so that non-ideal gas
        # characteristics for HYDROGEN can be captured.
        #
        # The model accounts for the conservation of MASS only. The model makes the
        # following assumptions:
        #
        # 1) The effective temperature after the Mach disk is equal to the
        #    stagnation temperature
        # 2) The Mach number after the Mach disk is assumed to be 1
        # 3) The static pressure at the Mach disk is assumed to be ambient pressure
        # 
        # No attempt is made to account for the change in entropy across the shock
        #
        # Created by Isaac Ekoto (March 10, 2011), updated by ESH (May, 2014)
        #*************************************************************************
        '''
        (rhot, Tt, ct, mdot, pt) = self.throat()
        # Calculate state at expansion conditions from assumptions 1 & 3:
        T2      = self.high_P_gas.T; 
        p2      = self.low_P_gas.P;
        rho2    = self.high_P_gas.therm.rho(T2, p2);
        gas_eff = Gas(self.high_P_gas.therm, T = T2, rho = rho2)

        # Calculate speed from assumption 2 (which is clearly faulty):
        Veff    = self.high_P_gas.therm.a(T2, rho2) # speed of sound at Mach disk [m/s]
        
        # Conservation of mass yeilds orifice diameter:
        orifice_eff = Orifice(sqrt(mdot/(rho2*Veff)*4/pi))
        
        return (gas_eff, orifice_eff, Veff)
    
    def _Birch2(self):
        '''
        #*************************************************************************
        # This function uses the second set of relationships developed by Birch et  
        # al. (Combust Sci Tech Vol 52 (1987) 161 - 171) to compute the effective 
        # area, velocity, and temperature after the Mach disk for an underexpanded 
        # jet. 
        # 
        # An additional feature that has been added to the model is the
        # inclusion of the Abel-Noble equation of state so that non-ideal gas
        # characteristics for HYDROGEN can be captured.
        #
        # The model is similar to the Birch1 model, except that it also accounts 
        # for the conservation of Momentum along with the conservation of MASS. 
        # The model makes the following assumptions:
        #
        # 1) The effective temperature after the Mach disk is equal to the
        #    stagnation temperature
        # 2) The static pressure at the Mach disk is assumed to be ambient pressure
        # 
        # No attempt is made to account for the change in entropy across the shock
        #
        # Created by Isaac Ekoto (March 10, 2011), updated by ESH (May, 2014)
        #*************************************************************************
        '''
        # throat conditions:
        (rhot, Tt, ct, mdot, pt) = self.throat()
        
        # Assumptions used to calculate state at expansion conditions:
        T2      = self.high_P_gas.T;# Effective Temperature [K]
        rho2    = self.high_P_gas.therm.rho(T2, self.low_P_gas.P); 
        gas_eff = Gas(self.high_P_gas.therm, T = T2, rho = rho2)
        
        # Simultaneous solution of conservation of mass and momentum (see Birch paper):
        Veff = ct + (pt - self.low_P_gas.P)/(ct*rhot) # Effective Velocity [m/s]
        
        # Conservation of mass yeilds orifice diameter:
        orifice_eff = Orifice(sqrt(mdot/(rho2*Veff)*4/pi))
        
        return (gas_eff, orifice_eff, Veff)
    
    def _Molkov(self):
        '''
        #*************************************************************************
        # This function uses the second set of relationships developed by Molkov et  
        # al. "PHYSICS AND MODELLING OF UNDER-EXPANDED JETS AND HYDROGEN DISPERSION
        # IN ATMOSPHERE" to compute the effective area, velocity, and temperature 
        # after the Mach disk for an underexpanded jet. 
        # 
        # An model includes both Abel-Noble and ideal gas formulations for the EOS.
        #
        # The model is similar to the Birch (1987) model, except that assumes:
        #
        # 1) A sonic velocity at the Mach disk
        # 2) The static pressure at the Mach disk is assumed to be ambient pressure
        # 
        # No attempt is made to account for the change in entropy across the shock
        #
        # Created by Isaac Ekoto (March 10, 2011)
        #*************************************************************************
        '''
        (rhot, Tt, ct, mdot, pt) = self.throat()
        
        Peff  =  self.low_P_gas.P
        Teff  =  self.high_P_gas.therm.T_IflowSonic(Tt, pt, Peff)
        gas_eff = Gas(self.high_P_gas.therm, T = Teff, P = Peff)
        
        Veff    = self.high_P_gas.therm.a(Teff, gas_eff.rho); # Effective Velocity [m/s]
        Aeff    = mdot/(gas_eff.rho*Veff); # Conservation of mass to get effective area [m**2]
        Deff    = sqrt(Aeff*4/pi); # Effective Diameter [m]
        orifice_eff = Orifice(Deff)
        
        return (gas_eff, orifice_eff, Veff)

# # define gas and Orifice
# Dia  = 1 #mm
# Pres = 4.#bar
# Tamb = 293.#K
# Tgas=53#K
# gasname = 'hydrogen'
# therm_gas   = EOS()    # hydrogen or helium, but when there is some
#                                              # problem when set gamma = 1.66, because of
#                                              # fsolve()
#                                              # change estimate to linspace(0.2,30.2,31), Okay!
#
# therm_air   = EOS(species="air")
# gas_high    = Gas(therm_gas, T=Tgas, P= Pres*1e5)
# gas_low     = Gas(therm_air, T=Tamb, P=1.01325e5)
# de          = Dia*1e-3          # exit diameter, m
# orifice     = Orifice(d=de)
# nn = NotionalNozzle(gas_high,orifice,gas_low)
# print(nn.calculate())

