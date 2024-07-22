# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 06:32:42 2023
Свойства диоксида углерода
Расчет зависимости теплоемкости от давления
@author: User
"""

import numpy as np
from pyfluids import Fluid, FluidsList, Input
import matplotlib.pyplot as plt

class PropSCO2:
    def __init__(self):
        
        self.fluid = Fluid(FluidsList.CarbonDioxide)
        
    def cp(self, temp, press):
        return self.fluid.with_state(
            Input.pressure(press), Input.temperature(temp)).specific_heat / 1000
    
    def diagramm_pressure(self, pressure, temperature):
        self.pressure = pressure # принимает кортеж
        self.temperature = temperature # принимает список
        
        ps = 200
        self.specific_heat = np.zeros(ps)
        self.pressure = np.linspace(self.pressure[0], self.pressure[1], ps)
        
        self.specific_heat = np.vectorize(self.cp)
        
        plt.figure(figsize=(6,5), dpi=300)
        
        for i in range(len(self.temperature)):
            plt.plot(self.pressure/1e6, self.specific_heat(self.temperature[i], self.pressure), 
                     label= str(self.temperature[i])+'$^{\circ}C$')
        
        plt.minorticks_on()
        plt.xlabel('Pressure, MPa', fontsize=10) 
        plt.ylabel('$C_p$, кДж/(кг.К)]', fontsize=10) 
        plt.grid(linestyle='--', linewidth=0.5, color='black') # сетка
        # plt.title('title', loc='left', fontsize=10)
        plt.legend()
        plt.tight_layout() # оптимизируем поля и расположение объектов
        file_name = 'specific_heat_pressure.png'
        plt.savefig(file_name, dpi = 300)
        plt.show()
        print(file_name)
        
    def diagramm_temperature(self, pressure, temperature):
        self.temperature_set = temperature # принимает кортеж
        self.pressure = pressure # принимает список
        
        ts = 200
        self.specific_heat = np.zeros(ts)
        self.temperature = np.linspace(self.temperature_set[0], self.temperature_set[1], ts)
        
        self.specific_heat = np.vectorize(self.cp)
        
        plt.figure(figsize=(6,5), dpi=300)
        
        for i in range(len(self.pressure)):
            plt.plot(self.temperature, self.specific_heat(self.temperature, self.pressure[i]), 
                     label= str(self.pressure[i]/1e6)+' МПа')
        
        plt.minorticks_on()
        plt.xlabel('Temperature, MPa', fontsize=10) 
        plt.ylabel('$C_p$, кДж/(кг.К)]', fontsize=10) 
        plt.grid(linestyle='--', linewidth=0.5, color='black') # сетка
        # plt.title('title', loc='left', fontsize=10)
        plt.legend()
        plt.tight_layout() # оптимизируем поля и расположение объектов
        file_name = 'specific_heat_temperature.png'
        plt.savefig(file_name, dpi = 300)
        plt.show()
        print(file_name)

    def turbine(self, pressure, pressure_rate, temperature):
        pressure_out = pressure
        pressure_in = pressure_out * pressure_rate
        state_in = self.fluid.with_state(Input.pressure(pressure_in), 
                    Input.temperature(temperature))
        state_out = state_in.isentropic_expansion_to_pressure(pressure_out)
        return abs(state_in.enthalpy - state_out.enthalpy)
    
    def compressor(self, pressure, pressure_rate, temperature):
        pressure_in = pressure
        pressure_out = pressure_in * pressure_rate
        state_in = self.fluid.with_state(Input.pressure(pressure_in), 
                    Input.temperature(temperature))
        state_out = state_in.isentropic_compression_to_pressure(pressure_out)
        return abs(state_in.enthalpy - state_out.enthalpy)

    def work(self, pressure, pressure_rate, temperature):
        # self.pressure = pressure # принимает кортеж
        # self.temperature = temperature # принимает кортеж
        
        ps = 100
        self.pressure = np.linspace(pressure[0], pressure[1], ps)
        self.turbine_diff = np.vectorize(self.turbine)
        self.compressor_diff = np.vectorize(self.compressor)

        plt.figure(figsize=(6,5), dpi=300)
        
        plt.plot(self.pressure/1e6, self.compressor_diff(pressure=self.pressure, 
                                        pressure_rate=pressure_rate, 
                                        temperature=temperature[0])/1000, 'b--', 
                                        label= 'compr')
        
        plt.plot(self.pressure/1e6, self.turbine_diff(pressure=self.pressure, 
                                        pressure_rate=pressure_rate, 
                                        temperature=temperature[1])/1000, 'g--', 
                                        label= 'turbine')
        
        work = (self.turbine_diff(pressure=self.pressure, 
                pressure_rate=pressure_rate, temperature=temperature[1]) -
                self.compressor_diff(pressure=self.pressure, 
                pressure_rate=pressure_rate, temperature=temperature[0])) /1000
        
        plt.plot(self.pressure/1e6, work, 'r', label= 'work')
        
        plt.vlines(7.38, 0, 150, linewidth=1, linestyle='--', color='black', label= '$P_{crit}$') # вертикальная линия
        
        
        plt.minorticks_on()
        plt.xlabel('Давление, MPa', fontsize=10) 
        plt.ylabel('Работа, кДж/кг', fontsize=10) 
        # plt.ylim([0, 200])
        plt.grid(linestyle='--', linewidth=0.5, color='black') # сетка
        title = '$\pi_k$ = ' + str(pressure_rate)
        plt.title(title, loc='left', fontsize=10)
        plt.legend()
        plt.tight_layout() # оптимизируем поля и расположение объектов
        file_name = 'work.png'
        plt.savefig(file_name, dpi = 300)
        plt.show()
        print(file_name)

    def compressibility(self, pressure, temperature):
        ps = 100
        pressure = np.linspace(pressure[0], pressure[1], ps)
        compr_func = lambda press: self.fluid.with_state(Input.pressure(press), 
                    Input.temperature(t)).compressibility
        compr = np.vectorize(compr_func)
        
        plt.figure(figsize=(6,5), dpi=300)
        for temp in temperature:
            t = temp
            plt.plot(pressure/1e6, compr(pressure), label= 'T = '+str(t)+'$^{\circ}C$')
        
        plt.vlines(7.38, 0, 1.1, linewidth=1, linestyle='--', color='black', 
                   label= '$P_{crit}$') # вертикальная линия

        plt.minorticks_on()
        plt.xlabel('Давление, MPa', fontsize=10) 
        plt.ylabel('Сжимаемость', fontsize=10) 
        plt.grid(linestyle='--', linewidth=0.5, color='black') # сетка
        plt.legend()
        plt.tight_layout() # оптимизируем поля и расположение объектов
        file_name = 'compressibility.png'
        plt.savefig(file_name, dpi = 300)
        plt.show()
        print(file_name)



        
        
        
        
    # def recuperator(self, temperature, pressure):

# ts = PropSCO2()
# ts.compressibility(pressure=(5e6, 10e6), temperature=(32, 40, 60, 100, 300, 500))

# ts.work(pressure=(5e6, 10e6), pressure_rate=2.5, temperature=(32, 600))

# temperature = [32, 40, 50, 60, 600]
# ts.diagramm_pressure(pressure=(1e6, 20e6), temperature=temperature)
# ts.diagramm_temperature(pressure=[7.4e6, 15e6], temperature=(32, 300))


