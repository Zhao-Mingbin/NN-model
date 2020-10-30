# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 10:23:11 2014

@author: xueli
"""
#%%
from Other_Nn import *

Dia  = 0 #mm
Pres = 4.#bar
Tamb = 293.#K
Tgas=53#K
gasname = 'hydrogen'

if gasname == 'H2':
    gam = 1.4
    MW_gas = 2.016
    covolume = 0.007691

    
if gasname == 'He':
    gam = 1.667
    MW_gas = 4.003
    covolume = 0.0029

# define gas and Orifice
therm_gas   = EOS(species=gasname)    # hydrogen or helium, but when there is some
                                             # problem when set gamma = 1.66, because of 
                                             # fsolve()
                                             # change estimate to linspace(0.2,30.2,31), Okay!

therm_air   = EOS(species="air")
gas_high    = Gas(therm_gas, T=Tgas, P= Pres*1e5)
gas_low     = Gas(therm_air, T=Tamb, P=1.01325e5)
de          = Dia*1e-3          # exit diameter, m
orifice     = Orifice(d=de)       # exit   


nn = NotionalNozzle(gas_high,orifice,gas_low)

#%%
class effective:
    def __init__(self, gas_eff = None, orifice_eff = None, Veff = None):
        '''class for effective source'''
        self.gas_eff, self.orifice_eff, self.Veff = gas_eff, orifice_eff, Veff
        self.rho = gas_eff.rho
        self.mdot, self.d = orifice_eff.mdot(self.rho,Veff), orifice_eff.d
        self.P, self.T = gas_eff.P, gas_eff.T

#%%    
#model_name=[ "YuceilOtugen","Ewan","Birch84","Birch87","Molkov"]
model_name="YuceilOtugen"
result = effective(*(nn.calculate(model=model_name)))


#%%
string = '''
***************** Yuceil ******************
*     velocity    = %5.4e            *
*     mdot        = %5.4e            *
*     radius      = %5.4e            *
*     Temperature = %5.4e            *
*     Pressure    = %5.4e            *
*******************************************

''' % (result.Veff, result.mdot, result.d/2.,result.T, result.P)

print (string)




          
