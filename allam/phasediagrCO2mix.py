# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 09:16:55 2024

@author: User
"""

import numpy as np
import pandas as pd
import CoolProp as CP
import matplotlib.pyplot as plt
import math

class PTdiagrmix:
    def __init__(self):
        pd.set_option('display.float_format', '{:.3f}'.format)
        self.p = pd.read_csv('cycle.csv') 
    
    
    def plot(self):
        plt.figure(figsize=(5, 4))
        carbondioxid = CP.AbstractState("HEOS", "CO2")
        pc1 = carbondioxid.keyed_output(CP.iP_critical)
        pmin1 = carbondioxid.keyed_output(CP.iP_min)
        pmax1 = pc1 * 100
        Tc1 = carbondioxid.keyed_output(CP.iT_critical)
        Tmin1 = carbondioxid.keyed_output(CP.iT_min)
        Tmax1 = Tc1 * 2
        
        fig = plt.figure(figsize = (5,5))
        
        # Saturation curve
        Ts1 = np.linspace(Tmin1, Tc1)
        ps1 = CP.CoolProp.PropsSI('P','T',Ts1,'Q',0,'CO2')
        
        # Labels
        plt.plot(Ts1,ps1,'orange', lw=3, label='carbondioxid')
        
        # Critical lines
        plt.axvline(Tc1, dashes = [2, 2], color='orange', lw=1)
        plt.axhline(pc1, dashes = [2, 2], color='orange', lw=1)
        
        water = CP.AbstractState('HEOS', 'Water')
        pc2 = water.keyed_output(CP.iP_critical)
        pmin2 = water.keyed_output(CP.iP_min)
        pmax2 = pc2 * 100
        Tc2 = water.keyed_output(CP.iT_critical)
        Tmin2 = water.keyed_output(CP.iT_min)
        Tmax2 = Tc2 * 2
        
        # Saturation curve
        Ts2 = np.linspace(Tmin2, Tc2, 10)
        ps2 = CP.CoolProp.PropsSI('P','T',Ts2,'Q',0,'water')
        
        # Labels
        plt.plot(Ts2,ps2,'aqua', lw=3, label='water')
        
        # Critical lines
        plt.axvline(Tc2, dashes = [2, 2], lw=1)
        plt.axhline(pc2, dashes = [2, 2], lw=1)
        
        x = np.append(self.p.temp, self.p.temp[0])
        y = np.append(self.p.pres, self.p.pres[0])
        
        # plt.scatter(x, y, color='red') # точка
        plt.plot(x[3:], y[3:], color='blue', label='CO2', lw=1) # точка
        plt.plot(x[:3], y[:3], color='red', label='mix', lw=1) # точка
        
        c = [300, 1000]
        # plt.fill_between(x, y1, y2)
        # plt.fill_between([x, y], 8e6, 18e6, color='purple', alpha= .2, label='cycle')
        # plt.plot(cycle_T, cycle_P, color="red", lw=2, label='cycle', linestyle='--')
        plt.ylim(1e6,1e8)
        plt.gca().set_yscale('log')
        plt.gca().set_xlim(220, 1100)
        plt.ylabel('Pressure [Pa]')
        plt.xlabel('Temperature [K]')
        plt.grid(linestyle='--', linewidth=0.5, color='black') # сетка
        plt.tight_layout()
        plt.legend()
        plt.tight_layout()
        plt.savefig('PTmix.jpg', dpi = 300)
        plt.show()
        
        
        
        
# diagr = PTdiagrmix()
# diagr.plot()