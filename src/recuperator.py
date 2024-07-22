# -*- coding: utf-8 -*-
"""
Created on Sun Jun 16 06:23:30 2024

@author: User
"""
# from __future__ import print_function
import CoolProp
from CoolProp.CoolProp import PropsSI, PhaseSI
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
# fig = plt.figure() 

class RecuperAlam:
    def __init__(self):
        pd.set_option('display.float_format', '{:.3f}'.format)
        self.p = pd.read_csv('cycle.csv') 
        self.mix = 'CO2[' + str(self.p.CO2[0]) + ']&water[' + str(self.p.H2O[0]) + ']'
        
        self.n_mix = [0,1]
        self.n_CO2 = [2,3,4,5,6,7,0]
        
        self.temp_hot_out = self.p.temp[2]
        self.temp_hot_in = self.p.temp[1]
        self.temp_cold_out = self.p.temp[7]
        self.temp_cold_in = self.p.temp[6]
        
    def plot_dt(self):
        x = np.linspace(0, 1)
        th = np.linspace(self.temp_hot_in, self.temp_hot_out)
        tc = np.linspace(self.temp_cold_out, self.temp_cold_in)
        dt = th - tc
        
        fig, ax = plt.subplots()
        
        ax.plot(x, th, color='red', label='T hot')
        ax.plot(x, tc, color='blue', label='T cold')
        
        ax.set_xlabel('x')
        ax.set_ylabel('Temperature [K]')
        ax.grid(linestyle='--', linewidth=0.5)
        ax.set_xlim(0, 1)

        ax1 = ax.twinx()

        ax1.plot(x, dt, color='green',linestyle='--',label='dT')
        
        ax1.set_ylabel('dT [K]')
        ax1.set_ylim(0,200)
        ax.set_ylim(200,1100)
        ax.legend(loc='lower left')
        ax1.legend(loc='upper right')
        fig.tight_layout()

    def plot_dc(self):
        
        # расчет второй производной энтальпии
        x = np.linspace(0, 1)
        th = np.linspace(self.p.temp[1], self.p.temp[2])
        tc = np.linspace(self.p.temp[7], self.p.temp[6])
        ph = np.linspace(self.p.pres[1], self.p.pres[2])
        pc = np.linspace(self.p.pres[7], self.p.pres[6])
        dt = th - tc
        
        d_sp_heat_hot = PropsSI('d(d(Hmass)/d(T)|P)/d(Hmass)|P','P',ph,'T',th,'CO2')
        d_sp_heat_cold = PropsSI('d(d(Hmass)/d(T)|P)/d(Hmass)|P','P',pc,'T',tc,'CO2')
        
        fig, ax = plt.subplots()
        
        ax.plot(x, d_sp_heat_hot, color='red', label='dCp/dT hot')
        ax.plot(x, d_sp_heat_cold, color='b', label='dCp/dT cold')
        
        ax.set_xlabel('x')
        ax.set_ylabel('dCp/dT ')
        ax.grid(linestyle='--', linewidth=0.5)
        ax.set_xlim(0, 1)
        # ax.set_ylim(200,1100)
        
        # расчет фактического температурного напора
        
        dtf = np.zeros(50)
        
        dt0 = self.p.temp[1] - self.p.temp[7]
        
        dtf = dt * d_sp_heat_hot
        
        integral = np.trapz(dtf)
        
        const = dt0 / integral
        
        df = dtf * const
        
        dt = np.zeros(50)
        
        for i in range(50):
            dt[i] = dt0 - np.trapz(df[:i])
        

        ax1 = ax.twinx()

        ax1.plot(x, dt, color='green',linestyle='--',label='dT')
        
        ax1.set_ylabel('dT [K]')
        # ax1.set_ylim(0,200)
        ax.legend(loc='lower left')
        ax1.legend(loc='upper right')
        fig.tight_layout()

    

rec = RecuperAlam()
rec.plot_dc()


    
    # def plotCO2(self, temp, press):
    #     t = np.linspace(temp[0], temp[1])
    #     p = np.linspace(press[0], press[1])

    #     h = PropsSI('H','T',t,'P',p,'CO2')
    #     plt.plot(h,p, color='blue', lw=2) 

    # def CO2(self):
    #     for i in [2,3,4,5,6]:
    #         self.plotCO2([self.p.temp[i], self.p.temp[i+1]],[self.p.pres[i], self.p.pres[i+1]])

    # def plotmix(self, temp, press):
    #     t = np.linspace(temp[0], temp[1])
    #     p = np.linspace(press[0], press[1])

    #     h = PropsSI('H','T|supercritical',t,'P',p,self.mix)
    #     plt.plot(h,p, color='red', lw=2) 

    # def mixCO2(self):
    #     for i in range(2,8):
    #         if i == 7: j = 0
    #         else: j = i + 1
            
    #         self.plotmix([self.p.temp[i], self.p.temp[j]],[self.p.pres[i], self.p.pres[j]])

    # def line(self):
    #     x = self.p.temp
    #     y = self.p.entr
    #     x = np.append(x, x[0])
    #     y = np.append(y, y[0])
        
    #     # plt.plot(y[3:],x[3:], color='red', label='CO2', lw=.5, linestyle='--') 
    #     # plt.plot(x, y, linewidth=1, linestyle='--')
    #     # plt.plot(y[:3],x[:3], color='blue', label='mix', lw=.5, linestyle='--') 
        
    #     plt.minorticks_on()
    #     plt.ylabel('Pressure [Pa]') 
    #     plt.xlabel('Enthalpy [J/kg]') 
    #     plt.grid(linestyle='--', linewidth=0.5, color='black') # сетка
    #     # plt.legend()
    #     plt.tight_layout()
    #     plt.show()

    # def points(self):
    #     x = self.p.temp
    #     y = self.p.entr
    #     plt.scatter(y,x)
    #     for i in [0,1,2,4,5,6,7]:
    #         plt.annotate(i, xy=(y[i]+20,x[i]-40))
        
        
    # def diagr(self):
    #     p = np.linspace(250,320)

    #     h = PropsSI('H','T',t,'Q',1,'CO2')
    #     plt.plot(s,t, color='grey', lw=2) 
    #     s = PropsSI('S','T',t,'Q',0,'CO2')
    #     plt.plot(s,t, color='grey', lw=2) 
        
        
    # def plotdens(self):
    #     self.CO2()
    #     # self.mixCO2()
    #     # self.points()
    #     # self.diagr()
    #     self.line()
        



# print(dens.p)
# plot = dens.CO2()
# plot = dens.mixCO2()
# plot = dens.line()














# for i in range(len(p)):
    
#     cycle_states[i,'P'] = p.pres[i]
#     cycle_states[i,'S'] = p.entr[i]
#     cycle_states[i,'T'] = p.temp[i]
#     cycle_states[i,'H'] = p.enth[i]
#     cycle_states[i,'D'] = p.dens[i]
    
# print(cycle_states)




# plt.rcParams.update({'font.size': 10})
# plot.savefig('ph-cyclopentane.jpg', dpi=300)



# plt.plot(x, f(x), color='red', label='$\\alpha$')
# # plt.minorticks_on()
# plt.xlabel('$x,^{\circ}C$') 
# plt.ylabel('$y_1$') 
# plt.grid(linestyle='--', linewidth=0.5, color='black') # сетка
# # plt.grid(which='minor', linestyle=':')
# plt.title('title', loc='left')
# plt.legend()
# # plt.tight_layout() # оптимизируем поля и расположение объектов
# plt.savefig('pic02.png', dpi = 300)
# plt.show()