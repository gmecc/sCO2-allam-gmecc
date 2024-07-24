# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 07:45:13 2024

Расчет распределения температуры в рекуператоре
"""

import pandas as pd
import numpy as np
from CoolProp.CoolProp import PropsSI
import matplotlib.pyplot as plt


class HeatEx:
    def __init__(self):
        pd.set_option('display.float_format', '{:.3f}'.format)
        self.p = pd.read_csv('cycle.csv') 
        self.mix = 'CO2[' + str(self.p.CO2[0]) + ']&water[' + str(self.p.H2O[0]) + ']'
        self.r = pd.DataFrame(columns=['x0_cycle','x1_cycle','x0','x1','x0res','x1res'], 
                              index=['hot','cold','dt'])
        
        self.r.loc['hot', 'x0_cycle'] = self.p.temp[1]
        self.r.loc['cold', 'x0_cycle'] = self.p.temp[7]
        self.r.loc['dt', 'x0_cycle'] = self.p.temp[1] - self.p.temp[7]
        
        self.r.loc['hot', 'x1_cycle'] = self.p.temp[2]
        self.r.loc['cold', 'x1_cycle'] = self.p.temp[6]
        self.r.loc['dt', 'x1_cycle'] = self.p.temp[2] - self.p.temp[6]
    
    def solve(self):
        n = 100
        x = np.linspace(0, 1, n)
        th = np.zeros(n)
        dh_hot = (self.p.enth[2] - self.p.enth[1]) / n      
        th[0] = self.p.temp[1]
        for i in range(n-1):
            th[i+1] = th[i] + dh_hot / PropsSI('C','T|supercritical',th[i],'P',self.p.pres[1],self.mix)
        
        x = np.linspace(1, 0, n)
        tc = np.zeros(n)
        dh_cold = (self.p.enth[7] - self.p.enth[6]) / n      
        tc[0] = self.p.temp[6]
        for i in range(n-1):
            tc[i+1] = tc[i] + dh_cold / PropsSI('C','T|supercritical',tc[i],'P',self.p.pres[6],'CO2')
        
        self.r.loc['hot', 'x0'] = th[0]
        self.r.loc['cold', 'x0'] = tc[-1]
        self.r.loc['dt', 'x0'] = th[0] - tc[-1]
        
        self.r.loc['hot', 'x1'] = th[-1]
        self.r.loc['cold', 'x1'] = tc[0]
        self.r.loc['dt', 'x1'] = th[-1] - tc[0]
        
        self.r.x0res =  self.r.x0_cycle - self.r.x0
        self.r.x1res =  self.r.x1_cycle - self.r.x1
        
        # plt.plot(x, th, color='red', label='hot')
        # plt.plot(x, tc, color='blue', label='cold')

        fig, ax = plt.subplots()
        fig.set_figheight(4)
        fig.set_figwidth(5)
        
        ax.plot(x[::-1], th, color='red', label='T hot')
        ax.plot(x, tc, color='blue', label='T cold')
        
        ax.set_xlabel('x')
        ax.set_ylabel('Temperature [K]')
        ax.grid(linestyle='--',linewidth=0.5,color='black')
        ax.set_xlim(0, 1)

        ax1 = ax.twinx()
        
        dt = th[::-1] - tc

        ax1.plot(x, dt, color='green',linestyle='--',label='dT')
        
        ax1.set_ylabel('dT [K]')
        ax1.set_ylim(0,140)
        ax.set_ylim(200,1100)
        ax.legend(loc='lower left')
        ax1.legend(loc='upper right')
        fig.tight_layout()
        fig.savefig('TX_recup.jpg', dpi = 300)
        plt.show()

    def spheat(self):
        plt.figure(figsize=(5, 4))
        n = 100
        x = np.linspace(0, 1, n)
        th = np.linspace(self.p.temp[1], self.p.temp[2], n)
        cph = PropsSI('C','T|supercritical',th,'P',self.p.pres[1],self.mix)
        plt.plot(th, cph, color='red', label='hot')
        
        tc = np.linspace(self.p.temp[7], self.p.temp[6], n)
        cph = PropsSI('C','T|supercritical',th,'P',self.p.pres[6],'CO2')
        plt.plot(tc, cph, color='blue', label='cold')
        
        plt.minorticks_on()
        plt.ylabel('Cp [J/kg]') 
        plt.xlabel('Temperature [K]') 
        plt.grid(linestyle='--', linewidth=0.5, color='black') # сетка
        plt.legend(loc='upper right')
        plt.tight_layout()
        plt.savefig('CpT_recup.jpg', dpi = 300)
        plt.show()
        

    def density(self):
        plt.figure(figsize=(5, 4))
        n = 100
        x = np.linspace(0, 1, n)
        th = np.linspace(self.p.temp[1], self.p.temp[2], n)
        cph = PropsSI('D','T|supercritical',th,'P',self.p.pres[1],self.mix)
        plt.plot(th, cph, color='red', label='hot')
        
        tc = np.linspace(self.p.temp[7], self.p.temp[6], n)
        cph = PropsSI('D','T|supercritical',th,'P',self.p.pres[6],'CO2')
        plt.plot(tc, cph, color='blue', label='cold')
        
        plt.minorticks_on()
        plt.ylabel('Density [J/kg]') 
        plt.xlabel('Temperature [K]') 
        plt.grid(linestyle='--', linewidth=0.5, color='black') # сетка
        plt.legend()
        plt.tight_layout()
        plt.savefig('density_recup.jpg', dpi = 300)
        plt.show()



        