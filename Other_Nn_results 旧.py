# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 10:23:11 2014

@author: xueli
"""
#%%
from Other_Nn import *

Dia  = 0 #mm
Pres = 4.#bar

gastype = 'H2'

if gastype == 'H2':
    gam = 1.4
    MW_gas = 2.016
    covolume = 0.007691
    gasname = 'H2_'
    filename = 'nn_'+ gasname + str(Pres)+'bar_'+str(Dia)+'mm.txt'
    
if gastype == 'He':
    gam = 1.667
    MW_gas = 4.003
    covolume = 0.0029
    gasname = 'He_'
    filename = 'nn_' + gasname + str(Pres)+'bar_'+str(Dia)+'mm.txt'
# define gas and Orifice
therm_gas   = EOS(species="hydrogen")    # hydrogen or helium, but when there is some
                                             # problem when set gamma = 1.66, because of 
                                             # fsolve()
                                             # change estimate to linspace(0.2,30.2,31), Okay!
Tamb = 293.
therm_air   = EOS(species="air")
gas_high    = Gas(therm_gas, T=53, P= Pres*1e5)
gas_low     = Gas(therm_air, T=Tamb, P=1.01325e5)
de          = Dia*1e-3          # exit diameter, m
orifice     = Orifice(d=de)       # exit   

#%%     
#therm_hyn = AbelNoble(0.,1.4,2.016)
#therm_air = AbelNoble(0.,1.4,28.97)
#gas_high = Gas(therm_hyn,T=300.,P=30.e5)
#gas_low  = Gas(therm_air,T=300.,P=1.01325e5)
#orifice  = Orifice(d = 0.001)
#filename = 'ref.txt'
#%%
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
Yuceil   = effective(*(nn.calculate(model='YuceilOtugen')))
Ewan     = effective(*(nn.calculate(model='EwanMoodie')))
Birch84  = effective(*(nn.calculate(model='Birch')))
Birch87  = effective(*(nn.calculate(model='Birch2')))
Molkov   = effective(*(nn.calculate(model='Molkov')))
# Harstad  = effective(*(nn.calculate(model='HarstadBellan')))

#%%
string = '''
***************** Yuceil ******************
*     velocity    = %5.4e            *
*     mdot        = %5.4e            *
*     radius      = %5.4e            *
*     Temperature = %5.4e            *
*     Pressure    = %5.4e            *
*******************************************
***************** Ewan   ******************
*     velocity    = %5.4e            *
*     mdot        = %5.4e            *
*     radius      = %5.4e            *
*     Temperature = %5.4e            *
*     Pressure    = %5.4e            *
*******************************************
***************** Birch84 *****************
*     velocity    = %5.4e            *
*     mdot        = %5.4e            *
*     radius      = %5.4e            *
*     Temperature = %5.4e            *
*     Pressure    = %5.4e            *
*******************************************
***************** Birch87 *****************
*     velocity    = %5.4e            *
*     mdot        = %5.4e            *
*     radius      = %5.4e            *
*     Temperature = %5.4e            *
*     Pressure    = %5.4e            *
*******************************************
***************** Molkov  *****************
*     velocity    = %5.4e            *
*     mdot        = %5.4e            *
*     radius      = %5.4e            *
*     Temperature = %5.4e            *
*     Pressure    = %5.4e            *
*******************************************
''' % (Yuceil.Veff, Yuceil.mdot, Yuceil.d/2.,Yuceil.T, Yuceil.P,
       Ewan.Veff,Ewan.mdot,Ewan.d/2.,Ewan.T,Ewan.P,
       Birch84.Veff,Birch84.mdot,Birch84.d/2.,Birch84.T,Birch84.P,
       Birch87.Veff,Birch87.mdot,Birch87.d/2.,Birch87.T,Birch87.P,
       Molkov.Veff,Molkov.mdot,Molkov.d/2.,Molkov.T,Molkov.P)

print (string)




          
