# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 09:16:55 2024

@author: User
"""

import numpy as np
import CoolProp as CP
import matplotlib.pyplot as plt
import math

class PTdiagrCO2:
    def plot(self):
        fluid = CP.AbstractState("HEOS", "CO2")
        pc = fluid.keyed_output(CP.iP_critical)
        pmin = fluid.keyed_output(CP.iP_min)
        pmax = pc * 100
        Tc = fluid.keyed_output(CP.iT_critical)
        Tmin = fluid.keyed_output(CP.iT_min)
        Tmax = Tc * 2
        
        fig = plt.figure(figsize = (5,5))
        
        # Saturation curve
        Ts = np.linspace(Tmin, Tc)
        ps = CP.CoolProp.PropsSI('P','T',Ts,'Q',0,'CO2')
        
        # Labels
        plt.plot(Ts,ps,'orange',lw = 3)
        
        # Critical lines
        plt.axvline(Tc, dashes = [2, 2])
        plt.axhline(pc, dashes = [2, 2])
        
        # Labels
        plt.text(Tmax-(Tmax-Tc)/2, pmax-(pmax-pc)/math.log(3), 'supercritical',ha= 'center')
        plt.text(Tc+(Tmax-Tc)/2, pc/5, 'supercritical gas', rotation = 0, ha= 'center')
        plt.text(Tmin+(Tc-Tmin)/2, pmax-(pmax-pc)*.96, 'supercritical liquid', rotation = 90, ha = 'center')
        plt.text(Tmin+(Tc-Tmin)/3, pc/3, 'liquid', rotation = 45,ha= 'center')
        plt.text(Tmin+(Tc-Tmin)/1.5, pc/5, 'gas', rotation = 45,ha= 'center')
        
        plt.ylim(pmin,pmax)
        plt.gca().set_yscale('log')
        plt.gca().set_xlim(Tmin, Tmax)
        plt.ylabel('Pressure [Pa]')
        plt.xlabel('Temperature [K]')
        plt.savefig('PT-CO2.jpg', dpi = 300)
        plt.tight_layout()
        plt.show()