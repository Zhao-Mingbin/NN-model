# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from pylab import *
from scipy import optimize, interpolate, integrate
from scipy.constants import R

# <markdowncell>

# ##AbelNoble: 
# 
# - thermodynamic class that Abel-Nobel EOS for various necessary state variables
# - default $b$, $\gamma$, and $M$ set for H<sub>2</sub>, and default ideal gas constant (don't know why anyone would ever want to change this--maybe to check different precisions??)
# - ideal gas law by setting $b = 0$
# - I can envision writing other classes with the same functions (i.e. something that calls refprop to find the same variables)
# 

class AbelNoble:
    def __init__(self, b = 7.6921e-3, gamma = 1.4, MW = 2.016):
        '''
        Class used to solve the Abel-Nobel equation of state for different thermodynamic parameters.
        See Schefer et al. International Journal of Hydrogen Energy 32 (2007) 2081-2093.
        Note: for ideal gas, use b = 0.
        
        Parameters
        ----------
        b - co-volume constant (7.6921e-3 m^3/kg for H2)
        gamma - specific heat ratio (1.4 for H2)
        MW - molecular weight of gas (2.016 g/mol for H2)
        R - universal gas constant (8.3145 J/K-mol)
        '''
        self.gamma = gamma
        self.b = b
        self.MW = MW
        self.Rg = R*1000/MW
        
    def h(self, T, P):
        '''enthalpy (J/kg) of a gas at temperature T (K) and pressure P (Pa) (not used in this case)'''
        cp = (self.gamma*R*1000./self.MW)/(self.gamma-1)
        return cp*T
    
    def errS(self, T0, rho0, T1, rho1):
        '''
        returns the difference in entropy (J/kg) between 2 states specified by the 
        temperatures and densities.
        '''
        b, g = self.b, self.gamma
        return (T1*((1-b*rho1)/rho1)**(g-1)) / (T0*((1-b*rho0)/rho0)**(g-1)) - 1    
    
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
        g, b = self.gamma, self.b
        cp = (g*R*1000./self.MW)/(g-1)
        return T0 + b*(P0 - P)/cp + 1./2*(V0**2 - V**2)/cp
        
    def rho_Iflow(self, rho0, Ma = 1, Ma0 = 0):
        '''
        solves for the density of gas after isentropic flow
        
        Parameters
        ----------
        rho0 - density before exansion (kg/m^3)
        Ma - mach number after epansion 
        Ma0 - Mach number before expansion
        
        Returns
        -------
        rho - density after expansion (kg/m^3)
        '''
        from scipy import optimize
        g = self.gamma
        b = self.b
        fun = lambda rho: rho0/(1-b*rho0) - rho/(1-b*rho)*((1+((g-1)/(2*(1-b*rho)**2)*Ma**2))/
                                                           (1+((g-1)/(2*(1-b*rho0)**2)*Ma0**2)))**(1/(g-1))
        return optimize.fsolve(fun, rho0)[0]
    
    def T_Iflow(self, T0, rho, Ma = 1):
        '''
        solves for the temperature of gas after an isentropic flow
        
        Parameters
        ----------
        T0 - temperature before expansion (K)
        rho - density after expansion (kg/m^3)
        Ma - mach number after expansion
        
        Returns
        -------
        T - temperature after expansion(K)
        '''
        g = self.gamma; b = self.b
        return T0/(1+((g-1)/(2*(1-b*rho)**2))*Ma**2)
    
    def P_Iflow(self, P0, rho, Ma = 1):
        '''
        solves for the temperature of gas after an isentropic flow
        
        Parameters
        ----------
        T0 - temperature before expansion (K)
        rho - density after expansion (kg/m^3)
        Ma - mach number after expansion
        
        Returns
        -------
        T - temperature after expansion(K)
        '''
        g = self.gamma; b = self.b
        return P0/(1+((g-1)/(2*(1-b*rho)**2))*Ma**2)**(g/(g-1))
    
    def T_IE(self, T0, rho0, rho):
        '''
        solves for the temperature of gas after an isentropic expansion, given an old and new density
        
        Parameters
        ----------
        T0 - temperature before expansion (K)
        rho0 - density before expansion (kg/m^3)
        rho - density after expansion (kg/m^3)
        
        Returns
        -------
        T - temperature after expansion(K)
        '''
        g = self.gamma; b = self.b
        return T0*((1-b*rho0)*rho/((1-b*rho)*rho0))**(g-1)
    
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
        return R*1000/self.MW*rho*T/(1-self.b*rho)
    
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
        
        return P*(1-self.b*rho)/(rho*R*1000/self.MW)
    
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
        return P/(R*1000/self.MW*T+P*self.b)
        
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
        return 1/(1-self.b*rho)*sqrt(self.gamma*R*1000/self.MW*T)
    
#    def T_IflowSonic(self, T0, P0, P):
#        '''
#        solves for the temperature of gas after an isentropic flow, assuming that the entrance
#        and exit velocities are sonic.  This assumes that the final velocity can be calculated using
#        the ideal gas relationship (a = sqrt(gamma*R_gas*T)).
#        
#        Parameters
#        ----------
#        T0 - temperature before expansion (K)
#        P0 - pressure before expansion (Pa)
#        P - pressure after expansion (Pa)
#        
#        Returns
#        -------
#        T - temperature after expansion(K)
#        '''
#        g = self.gamma
#        return 2*T0/(g + 1) + (g-1)/(g+1)*T0/(1-self.rho(T0, P0)*self.b)**2
    
    def T_IflowSonic(self, T0, P0, P):
        from scipy import optimize
        '''
        solves for the temperature of gas after an adiabatic flow, assuming that the entrance
        and exit velocities are sonic
        
        Parameters
        ----------
        T0 - temperature before expansion (K)
        P0 - pressure before expansion (Pa)
        P - pressure after expansion (Pa)
        
        Returns
        -------
        T - temperature after expansion(K)
        '''
        v0 = self.a(T0, self.rho(T0, P0))
        cp = R*1000./self.MW*(self.gamma/(self.gamma - 1))
        def err(T):
            v = self.a(T, self.rho(T, P))
            return cp*T0 + v0**2/2. - cp*T - v**2/2.
        T = optimize.newton(err, T0)
        return T

# <rawcell>

# r = AbelNoble()
# #print r.T_IflowV(160, 100, 1013250, 1000, 101325)
# #print r.T_IE(160, 2, 2.3)
# #print r.rho_Iflow(50, 2)
# #print r.T_Iflow(160., r.rho_Iflow(50, 4), 4)
# print r.T_Iflow(160., .1)
# print r.P(r.T_Iflow(160., r.rho_Iflow(50, 4), 4), r.rho_Iflow(50, 4))

# <rawcell>

# r = RefProp()
# #print r.T_IflowV(160, 100, 1013250, 1000, 101325)
# #print r.T_IE(160, 2, 2.3)
# print r.rho_Iflow(10, 2)
# print r.T_Iflow(160., .1)
# #print r.rho_Iflow3(50, 2)
# 
# #print r.T_Iflow(160., r.rho_Iflow2(50, 4), 4)
# #print r.P(r.T_Iflow(160., r.rho_Iflow2(50, 4), 4), r.rho_Iflow(50, 4))

# <codecell>

class RefProp:
    def __init__(self, species = ['hydrogen'], x = [1.]):
        '''
        Class that uses Refprop for equation of state calculations.
        '''
        import refprop as rp
        self.rp = rp
        if species == 'air':
            self.fld = self.rp.setup(u'def', u'air')
            self.x = self.fld['x']
        else:
            if size(x) != size(species):
                warnings.warn('mole fractions and species lists not the same length')
            self.fld = self.rp.setup(u'def', map(lambda s: unicode(s), species))
            self.x = x

    def h(self, T, P):
        '''enthalpy (J/kg) of a gas at temperature T (K) and pressure P (Pa)'''
        fld = self.flsh(T = T, P = P/1000.)
        return self.perKg(fld['h'])
    
    def errS(self, T0, rho0, T1, rho1):
        '''
        returns the error in entropy (J/kg-K) between 2 states specified by the 
        temperatures and densities.
        '''
        s0 = self.flsh(T = T0, rho = rho0)
        s1 = self.flsh(T = T1, rho = rho1)
        return self.perKg(s1['s']) / self.perKg(s0['s']) - 1
    
    def rho_Iflow(self, rho0, Ma = 1, Ma0 = 0):
        '''
        solves for the density of gas after isentropic flow
        
        Parameters
        ----------
        rho0 - density before exansion (kg/m^3)
        Ma - mach number during epansion 
        Ma0 - mach number before expansion
        
        Returns
        -------
        rho - density after expansion (kg/m^3)
        '''
        from scipy import optimize
        T0 = 273.
        sf = self.flsh(T = T0, rho = rho0)
        def errT(T0):
            T0 = float(T0)
            s0 = self.flsh(T = T0, rho = rho0)
            def errh(h):
                h = float(h)
                sf = self.flsh(S = self.perKg(s0['s']), H = h)
                hf = self.perKg(s0['h'])  + (Ma0*s0['w'])**2/2. - (Ma*sf['w'])**2/2.
                return hf - h
            hf = optimize.root(errh, self.perKg(s0['h']))['x']
            Tf = self.T_IflowV(T0, Ma0*s0['w'], s0['p']*1e3, Ma*self.fld['w'], self.fld['p']*1e3)
            return Tf - self.fld['t']
        T0 = optimize.root(errT, T0)
        return self.fld['D']*self.rp.wmol(self.x)['wmix']

    def T_Iflow(self, T0, rho, Ma = 1):
        '''
        solves for the temperature of gas after an isentropic flow
        
        Parameters
        ----------
        T0 - temperature before expansion (K)
        rho - density after expansion (kg/m^3)
        Ma - mach number during expansion
        
        Returns
        -------
        T - temperature after expansion(K)
        '''
        from scipy import optimize
        def errS(x):
            [P0, Pf] = x
            P0, Pf = float(P0), float(Pf)
            s0 = self.flsh(P = P0, T = T0)
            sf = self.flsh(P = Pf, rho = rho)
            hf = self.perKg(s0['h']) - (Ma*self.fld['w'])**2/2.
            return array([self.fld['s'] - s0['s'], self.perKg(self.fld['h']) - hf])
        P0, Pf = optimize.fsolve(errS, array([101325., 101325.]))
        #P0, Pf = optimize.minimize(errS, array([101325., 101325.]), method = 'SLSQP',#'TNC',#'SLSQP',#'L-BFGS-B',
        #                           bounds = ((1, None), (1, None))
        #                           )['x']
        self.flsh(P = float(Pf), rho = rho)
        return self.fld['t']
    
    def P_Iflow(self, P0, rho, Ma = 1):
        '''
        solves for the pressure of gas after an isentropic flow
        
        Parameters
        ----------
        T0 - temperature before expansion (K)
        rho - density after expansion (kg/m^3)
        Ma - mach number during expansion
        
        Returns
        -------
        P - Pressure after expansion(K)
        '''
        from scipy import optimize
        def errS(x):
            [T0, Pf] = x
            T0, Pf = float(T0), float(Pf)
            s0 = self.flsh(P = P0, T = T0)
            sf = self.flsh(P = Pf, rho = rho)
            hf = self.perKg(s0['h']) - (Ma*self.fld['w'])**2/2.
            return array([self.fld['s'] - s0['s'], self.perKg(self.fld['h']) - hf])
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
        self.flsh(T = T0, P = P0)
        h = self.perKg(self.fld['h']) + V0**2/2. - V**2/2.
        self.fld = self.flsh(P = P, H = h)
        return self.fld['t']
    
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
        self.flsh(T = T0, P = P0)
        V0 = self.fld['w']
        h0 = self.perKg(self.fld['h'])
        def errH(T):
            self.flsh(T = T, P = P)
            V = self.fld['w']
            h = h0 + V0**2/2. - V**2/2.
            return h - self.perKg(self.fld['h'])
        T = optimize.newton(errH, T0)
        return T
    
    def T_IE(self, T0, rho0, rho):
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
        self.flsh(T = T0, rho = rho0)
        self.flsh(rho = rho, S = self.perKg(self.fld['s']))
        return self.fld['t']
    
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
        self.flsh(T =  T, rho = rho)
        return self.fld['p']*1000.
    
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
        self.flsh(P = P, rho = rho)
        return self.fld['t']
    
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
        self.flsh(T = T, P = P)
        return self.fld['D']*self.rp.wmol(self.x)['wmix']
    
    def flsh(self, T = None, P = None, rho = None, H = None, S = None, U = None, Q = None):
        self.rp.resetup(self.fld)
        MW = self.rp.wmol(self.x)['wmix']
        if sum([x != None for x in [T, P, rho, H, S, U, Q]]) != 2:
            raise UserWarning('improperly defined thermodynamic state')
        elif T is not None and P is not None:
            self.fld = self.rp.flsh(u'tp', T, P/1000., self.x)
        elif T is not None and rho is not None:
            self.fld = self.rp.flsh(u'td', T, rho/MW, self.x)
        elif T is not None and H is not None:
            self.fld = self.rp.flsh(u'th', T, H*MW/1000., self.x)
        elif T is not None and S is not None:
            self.fld = self.rp.flsh(u'ts', T, S*MW/1000., self.x)
        elif T is not None and U is not None:
            self.fld = self.rp.flsh(u'te', T, U*MW/1000., self.x)
        elif T is not None and Q is not None:
            self.fld = self.rp.flsh(u'tq', T, Q, self.x)
        elif P is not None and rho is not None:
            self.fld = self.rp.flsh(u'pd', P/1000., rho/MW, self.x)
        elif P is not None and H is not None:
            self.fld = self.rp.flsh(u'ph', P/1000., H*MW/1000., self.x)
        elif P is not None and S is not None:
            self.fld = self.rp.flsh(u'ps', P/1000., S*MW/1000., self.x)
        elif P is not None and U is not None:
            self.fld = self.rp.flsh(u'pe', P/1000., U*MW/1000., self.x)
        elif P is not None and Q is not None:
            self.fld = self.rp.flsh(u'pq', P/1000., Q, self.x)
        elif rho is not None and H is not None:
            self.fld = self.rp.flsh(u'dh', rho/MW, H*MW/1000., self.x)
        elif rho is not None and S is not None:
            self.fld = self.rp.flsh(u'ds', rho/MW, S*MW/1000., self.x)
        elif rho is not None and U is not None:
            self.fld = self.rp.flsh(u'de', rho/MW, U*MW/1000., self.x)
        elif rho is not None and Q is not None:
            self.fld = self.rp.flsh(u'dq', rho/MW, Q, self.x)
        elif H is not None and S is not None:
            self.fld = self.rp.flsh(u'hs', H*MW/1000., S*MW/1000., self.x)
        elif H is not None and U is not None:
            self.fld = self.rp.flsh(u'he', H*MW/1000., U*MW/1000., self.x)
        elif H is not None and Q is not None:
            self.fld = self.rp.flsh(u'hq', H*MW/1000., Q, self.x)
        elif S is not None and U is not None:
            self.fld = self.rp.flsh(u'es', U*MW/1000., S*MW/1000., self.x)
        elif S is not None and Q is not None:
            self.fld = self.rp.flsh(u'sq', S*MW/1000., Q, self.x)
        elif U is not None and Q is not None:
            self.fld = self.rp.flsh(u'eq', U*MW/1000., Q, self.x)
        return self.fld
            
    def perKg(self, x):
        return x/self.rp.wmol(self.x)['wmix']*1000
    
    def perMol(self, x):
        return x*self.rp.wmol(self.x)['wmix']/1000
    
 #   def TH(self, P, Q):
 #       '''returns the temperature and enthalpy given the pressure and quality
 #       
 #       Parameters
 #       P - pressure (Pa)
 #       Q - quality (1 = vapor, 0 = liquid)
 #       
 #       Returns
 #       -------
 #       TH : tuple
 #           temperature (K), enthalpy (J/kg)
 #       '''
 #       self.fld = self.flsh(P = P, Q = Q)
 #       return self.fld['t'], self.perKg(self.fld['h'])
    
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
        self.flsh(T = T, rho = rho)
        return self.fld['w']

# <codecell>

class CoolProp:
    def __init__(self, species = ['hydrogen'], x = [1.]):
        '''
        Class that uses CoolProp for equation of state calculations.
        '''
        from CoolProp import CoolProp
        self.cp = CoolProp
        if len(species) == 1:
            self.spec = species[0]
        else:
            if size(x) != size(species):
                warnings.warn('mole fractions and species lists not the same length')
            s = 'REFPROP-MIX:'
            for g, xspec in zip(species, x):
                s += g + '[' + str(xspec) + ']' + '&'
            self.spec = s[:-1]

    def h(self, T, P):
        '''enthalpy (J/kg) of a gas at temperature T (K) and pressure P (Pa)'''
        return self.cp.PropsSI('H', 'T', T, 'P', P, self.spec)
            
    def errS(self, T0, rho0, T1, rho1):
        '''
        returns the error in entropy (J/kg-K) between 2 states specified by the 
        temperatures and densities.
        '''
        s0 = self.cp.PropsSI('S', 'T', T0, 'D', rho0, self.spec)
        s1 = self.cp.PropsSI('S', 'T', T1, 'D', rho1, self.spec)
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
        sf = self.cp.State(self.spec, {'T':T0, 'D':rho0})
        def errT(T0):
            s0 = self.cp.State(self.spec, {'T':T0, 'D':rho0})
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
            s0 = self.cp.State(self.spec, {'P':P0/1e3, 'T':T0})
            sf = self.cp.State(self.spec, {'P':Pf/1e3, 'D':rho})
            hf = s0.get_h()*1e3 + (Ma0*s0.get_speed_sound())**2/2. - (Ma*sf.get_speed_sound())**2/2.
            return array([sf.get_s() - s0.get_s(), sf.get_h()*1e3 - hf])
        P0, Pf = optimize.fsolve(errS, array([101325., 101325.]))
        return self.cp.PropsSI('T', 'P', Pf, 'D', rho, self.spec)
    
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
            s0 = self.cp.State(self.spec, {'P':P0/1e3, 'T':T0})
            sf = self.cp.State(self.spec, {'P':Pf/1e3, 'D':rho})
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
        h = self.cp.PropsSI('H', 'T', T0, 'P', P0, self.spec) + V0**2/2. - V**2/2.
        return self.cp.PropsSI('T', 'P', P, 'H', h, self.spec)
    
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
        S0 = self.cp.PropsSI('S', 'T', T0, 'D', rho0, self.spec)
        def errD(T):
            return rho - self.cp.PropsSI('D', 'T', T, 'S', S0, self.spec)
        return optimize.newton(errD, T0)
    
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
        return self.cp.PropsSI('P', 'T', T, 'D', rho, self.spec)
    
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
        return self.cp.PropsSI('T', 'D', rho, 'P', P, self.spec)
    
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
        return self.cp.PropsSI('D', 'T', T, 'P', P, self.spec)
 
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
        return self.cp.PropsSI('A', 'T', T, 'D', rho, self.spec)
    
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
            h = self.cp.PropsSI('H', 'T', T0, 'P', P0, self.spec) + V0**2/2. - V**2/2.
            return h - self.cp.PropsSI('H', 'T', T, 'P', P, self.spec)
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
            h = self.cp.PropsSI('H', 'T', T0, 'P', P0, self.spec) + V0**2/2. - V**2/2.
            return h - self.cp.PropsSI('H', 'T', T, 'P', P, self.spec)
        T = optimize.newton(errH, T0)
        return T

# <rawcell>

# cp = CoolProp()
# rp = RefProp()
# an = AbelNoble()
# for mod in [cp, rp, an]:
#     print mod.errS(250., .1, 260, .109999)

# <rawcell>

# import refprop as rp
# rp.flsh()
# #cp.State.State('H2', {'T':293, 'D':1.05})

# <rawcell>

# #cp = CoolProp(['hydrogen', 'oxygen'], [.5, .5])
# #rp = RefProp(['hydrogen', 'oxygen'], [.5, .5])
# 
# an = AbelNoble()
# cp = CoolProp()
# rp = RefProp()
# #for mod in [an, cp, rp]:
# #    print mod.rho_Iflow(20, 2)
# #    % timeit mod.rho_Iflow(20, 2)
#     
# print ''
# for mod in [an, cp, rp]:
#     print mod.T_Iflow(160, .1, 2)
# #    % timeit mod.T_Iflow(160, 100., 2)

# <markdowncell>

# ##H2_Combustion
# 
# A module that performs H<sub>2</sub> combustion calculations.  Initilizes some functions of mass fraction for product temperaure, density, molecular weight, and the derivative of the density with respect to mass fraction.
# 
# 
# H<sub>2</sub> + $\frac{\eta}{2}$ (O<sub>2</sub> + 3.76 N<sub>2</sub>) -> $\max(0, 1-\eta)$ H<sub>2</sub> + $\min(1, \eta)$ H<sub>2</sub>O + $\max\left(0, \frac{\eta-1}{2}\right)$ O<sub>2</sub> + $3.76\frac{\eta}{2}$ N<sub>2</sub>'

# <codecell>

class H2_Combustion:
    def __init__(self, Treac, P = 101325., numpoints = 5000,
                 MW = {'H2':2.016, 'O2':31.999, 'N2':28.013, 'H2O':18.015}, DHc = 118.83e6):
        '''
        Class that performs H2 combustion chemistry calculations.
        
        Initilizes some interpolating functions for T_prod, MW_prod, rho_prod, and drhodf
        
        Parameters
        ----------
        Treac : float
            temperature of reactants (K)
        P : float
            pressure (Pa), default value is 101325.
        numpoints : int
            number of points to solve for temperature to create 
            interpolating functions, default value is 5000
        MW : dict
            molecular weights of species (H2, O2, N2, and H2O)
        DHc : float
            heat of combustion for hydrogen (J/kg), default value of 188.83e6
            
        Contents
        --------
        self.T_prod(f): array_like
            temperature of products (K) at a mass fraction, f
        self.MW_prod:
            mixture averaged molecular weight of products (g/mol) at a mass fraction, f 
            '''
        print ('initilizing chemistry...'), 
        sys.stdout.flush()
        self.MW = MW
        self.DHc = DHc
        self.Treac, self.P = Treac, P
        self.fstoich = 1/(1.5+3.76/2)*MW['H2']/(MW['H2']*1/(1.5+3.76/2) +
                                                MW['O2']*0.5/(1.5+3.76/2) + 
                                                MW['N2']*3.76/2/(1.5+3.76/2) )
        fvals = linspace(0, 1, numpoints);
        T = array(map(lambda fval: self._T_combustion(Treac, fval), fvals))
        MWvals = array(map(lambda f: self._MWmix(self._Yprod(f)), fvals))
        self.MW_prod = lambda f: self._MWmix(self._Yprod(f))
        self.T_prod = lambda f: interpolate.interp1d(fvals, T)(f)
        self.rho_prod = lambda f: P*self.MW_prod(f)/(1000*R*self.T_prod(f))
        self.drhodf = lambda f: interpolate.interp1d(fvals, P/(1000*R*T)*(self._D(fvals, MWvals) - 
                                                                          self._D(fvals, T)/T*MWvals))(f)
        print ('done.')
        sys.stdout.flush()
     
    def reinitilize(self, Treac, P = 101325., numpoints = 1000,
                    MW = {'H2':2.016, 'O2':31.999, 'N2':28.013, 'H2O':18.015}, DHc = 118.83e6):
        '''
        Reinitilizes class to new temperature, pressure, etc.  Can be used rather 
        than creating a new instance taking up additional memory.'''
        self.__init__(Treac, P, numpoints, MW, DHc)
    
    def _D(self, x, y):
        '''numerical derivative of y WRT x'''
        return append(append((y[1] - y[0])/(x[1] -  x[0]),
                             (y[2:] - y[:-2])/(x[2:] - x[:-2])),
                      (y[-2] - y[-1])/(x[-2] - x[-1]))
    
    def _MWtot(self, eta):
        return self.MW['H2'] + eta/2.*(self.MW['O2'] + 3.76*self.MW['N2'])

    def _MWmix(self, Y):
        '''returns the mixture averaged molecular weight from a mole fraction'''
        MWmix = 0
        for spec, Yval in Y.items():
            MWmix += Yval/self.MW[spec]
        MWmix = 1./MWmix
        return MWmix
   
    def _eta(self, f):
        MW = self.MW
        eta = 2*(MW['H2']/(1.e-99*ones_like(f) +f) - MW['H2'])/(MW['O2'] + 3.76*MW['N2'])
        return eta
    
    def _Yprod(self, f):
        '''
        the mass fractions of combustion products as a function of the mixture fraction

        Parameters
        ----------
        f = mixture fraction

        Returns
        -------
        Y - dictionary of mass fractions (kg/kg)
        '''
        eta = self._eta(f)
        return {'H2':(1-eta)*((1-eta) > 0)*self.MW['H2']/self._MWtot(eta),
                'H2O':(eta*(eta < 1) + (eta > 1))*self.MW['H2O']/self._MWtot(eta),
                'O2':((eta-1)/2.*((eta-1) > 0))*self.MW['O2']/self._MWtot(eta),
                'N2':eta/2.*3.76*self.MW['N2']/self._MWtot(eta)}
    
    def _Yreac(self, f):
        '''
        the mass fractions of combustion reactants as a function of the mixture fraction
        (assumes that there is no H2O as a reactant)
        '''
        Y = {'H2':self.MW['H2']/self._MWtot(self._eta(f)), 
             'H2O':0.*f, 
             'O2':self._eta(f)/2.*self.MW['O2']/self._MWtot(self._eta(f)),
             'N2':self._eta(f)/2*3.76*self.MW['N2']/self._MWtot(self._eta(f))}
        return Y

    def _cp_mol(self, T): 
        '''
        heat capactiy of posible species as a function of temperature
        data from NIST chemistry webbook on 5/7/14

        Parameters
        ----------
        T - temperature (K)

        Returns
        -------
        cp - dictionary of heat capacities [J/(mol K)]
        '''
        t = T/1000.
        f = array([1., t, t**2, t**3, 1/t**2])

        if (T > 100.) and (T <= 700):
            O2 = array([31.32234, -20.23531, 57.86644, -36.50624, -0.007374])
        elif (T > 700) and (T <= 2000):
            O2 = array([30.03235, 8.772972, -3.988133, 0.788313, -0.741599])
        elif (T > 2000) and (T <= 6000):
            O2 = array([20.91111, 10.72071, -2.020498, 0.146449, 9.245722])
        else:
            O2 = array(5*[nan])

        if (T > 100.) and (T <= 500):
            N2 = array([28.98641,  1.853978,  -9.647459,  16.63537,  0.000117])
        elif (T > 500) and (T <= 2000.):
            N2 = array([19.50583,  19.88705,  -8.598535,  1.369784,  0.527601])
        elif (T > 2000) and (T <= 6000):
            N2 = array([35.51872,  1.128728,  -0.196103,  0.014662,  -4.55376])
        else:
            N2 = array(5*[nan])

        if (T > 273.) and (T <= 1000):#Webbook says this is only valid to 298, but I lowered it a bit
            H2 = array([33.066178, -11.363417, 11.432816, -2.772874, -0.158558])
        elif (T > 1000) and (T <= 2500):
            H2 = array([18.563083, 12.257357, -2.859786, 0.268238, 1.97799])
        elif (T > 2500) and (T <= 6000):
            H2 = array([43.41356, -4.293079, 1.272428, -0.096876, -20.533862])
        else:
            H2 = array(5*[nan])


        if (T > 500) and (T <= 1700):
            H2O = array([30.092, 6.832514, 6.793435, -2.53448, 0.082139])
        elif (T > 1700) and (T <= 6000):
            H2O = array([41.96426, 8.622053, -1.49978, 0.098119, -11.15764])
        else:
            H2O = array(5*[0.]) # really should be error here

        cp = {'O2': dot(f, O2), 'N2': dot(f, N2), 'H2': dot(f, H2), 'H2O': dot(f, H2O)}
        return cp

    def _cp_mass(self, T, Y):
        '''heat capacity of a mixture J/kg-K

        Parameters
        ----------
        T - Temperature (K)'''
        cp = 0.
        for species, Yval in Y.items():
            if Yval >= 1e-4:
                cp += self._cp_mol(T)[species]/self.MW[species]*1000*Yval
        return cp

    def _T_combustion(self, T_reac, f):
        '''combustion temperature (K)'''
        DHc = self.DHc*self._Yprod(f)['H2O']*self.MW['H2']/self.MW['H2O'] # heat of combustion [J/kg_total]
        H0 = self._cp_mass(T_reac, self._Yreac(f))*T_reac #J/kg
        H0 *= self._MWmix(self._Yreac(f))/self._MWmix(self._Yprod(f)) # ???
        H = H0 + DHc
        T = optimize.brentq(lambda T: T*self._cp_mass(T, self._Yprod(f)) - H, T_reac, 6000)
                            #H/self._cp_mass(T_reac, self._Yreac(f)))[0]
        return T

