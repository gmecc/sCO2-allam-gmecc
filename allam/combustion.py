# -*- coding: utf-8 -*-
"""
Created on Sun Jan  7 10:38:03 2024
Цикл Аллама
РАСЧЕТ СОСТАВА ГАЗОВОЙ СМЕСИ КАМЕРЫ СГОРАНИЯ ЦИКЛА АЛЛАМА 

Входные данные:
    - температура газа на выходе камеры сгорания
    - температура рециркалируемого СО2
    
Выхадные данные:
    - коэффициент разбавления
    - состав газа на выходе КС

Температура в К
"""

import pandas as pd
import numpy as np
from CoolProp.CoolProp import PropsSI
from scipy.optimize import root
import matplotlib.pyplot as plt
from scipy import constants as cst

pd.set_option('display.float_format', '{:.4f}'.format)
pd.set_option('display.max_columns', None)
pd.options.mode.chained_assignment = None  # default='warn'

class Combust:
    def __init__(self):
        
        # топливо - метан
        self.temp0 = cst.convert_temperature(0, 'C', 'K') # 'C', 'K', 'F', 'R'
        
        self.fuel = pd.Series({'fluid': 'methane', 'calorific': 35800000.})
        
        self.gas = pd.DataFrame(columns=['mol_mass', 'mol', 'mass'], 
                              index=['CO2', 'H2O', 'CO2_recyc'], dtype=float)
        
        self.gas.loc['mol_mass', 'CO2'] = PropsSI('molarmass','CO2')
        self.gas.loc['mol_mass', 'CO2'] = PropsSI('molarmass','CO2')
        self.gas.loc['mol_mass', 'H2O'] = PropsSI('molarmass','H2O')
        
        self.gas_in = pd.DataFrame(columns=['mol_mass', 'mol', 'mass'],
                index=['CH4', 'O2', 'CO2'], dtype=float)
        
        self.gas_in.loc['mol_mass', 'CH4'] = PropsSI('molarmass','CH4')
        self.gas_in.loc['mol_mass', 'O2'] = PropsSI('molarmass','O2')
        self.gas_in.loc['mol_mass', 'CO2'] = PropsSI('molarmass','CO2')

        spHvol_CO2_coef = [1.63479959e+03, 9.75263813e-01, -5.45793612e-04, 1.83324681e-07, -2.67917924e-11]
        spHvol_H2O_coef = [1.49735370e+03, 1.18320983e-01, 1.79154428e-04, -8.92543875e-08, 1.34614261e-11]
        spHvol_O2_coef = [1.30543577e+03, 1.74761949e-01, 6.42329236e-05, -8.82987243e-08, 2.28704857e-11]

        self.spHv_CO2_pol = np.polynomial.Polynomial(spHvol_CO2_coef)
        self.spHv_H2O_pol = np.polynomial.Polynomial(spHvol_H2O_coef)
        self.spHv_O2_pol = np.polynomial.Polynomial(spHvol_O2_coef)
        
        self.g = pd.Series([np.nan], index=['gas_vol'])
        
        
        # РАСЧЕТ СОСТАВА ПРОДУКТОВ СГОРАНИЯ
    def _massFraction(self, k_recyc):         
        self.k_recyc = k_recyc
        
        # Количество кислорода, необходимое для сжигания топлива м3/м3
        self.volume_O2 = 2
        self.volume_CH4 = 1        
        
        # Объем продуктов сгорания при k_recyc=1
        self.volume_H2O = 2
        self.volume_CO2_norm = 1
        self.volume_gas_norm = self.volume_H2O + self.volume_CO2_norm
        
        # Объем продуктов сгорания при k_recyc>1
        self.volume_CO2_recyc = self.volume_gas_norm * self.k_recyc
        self.volume_gas = self.volume_gas_norm * ( self.k_recyc + 1 )
        self.volume_CO2 = self.volume_gas - self.volume_H2O
        
        # Объем продуктов сгорания
        self.g['gas_vol'] = self.volume_gas
        
        # Состав продуктов сгорания - мольные доли
        self.gas.mol.at['CO2'] = self.volume_CO2 / self.g.gas_vol
        self.gas.mol.at['CO2_recyc'] = self.volume_CO2_recyc / self.g.gas_vol
        self.gas.mol.at['H2O'] = self.volume_H2O / self.g.gas_vol
        self.g['mol_mass'] = (self.gas.mol_mass[0:2] * self.gas.mol[0:2]).sum()
        
        # Массовые доли продуктов сгорания
        self.gas.mass = (self.gas.mol * self.gas.mol_mass / self.g.mol_mass).round(4)
        self.g['gas_mass_sum'] = self.gas.mass.sum()
        
        # Мольные доли газа на входе
        self.volume_in = self.volume_O2 + self.volume_CH4 + self.volume_CO2_recyc
        self.gas_in.mol.at['O2'] = self.volume_O2 / self.volume_in
        self.gas_in.mol.at['CH4'] = self.volume_CH4 / self.volume_in
        self.gas_in.mol.at['CO2'] = self.volume_CO2_recyc / self.volume_in 
        
        # Массовые доли газа на входе
        self.gas_in_mol_mass = (self.gas_in.mol_mass * self.gas_in.mol).sum()
        self.gas_in.mass = (self.gas_in.mol * self.gas_in.mol_mass / self.gas_in_mol_mass).round(4)
        
        

        # ФУНКЦИЯ РАСЧЕТА ТЕМПЕРАТУРЫ ГОРЕНИЯ 
    def _tempFunc(self, temp):
        
        sph_CO2 = self.spHv_CO2_pol(temp) 
        sph_H2O = self.spHv_H2O_pol(temp)
        
        temp_O2 = 15
        sph_O2 = self.spHv_O2_pol(temp_O2)
        
        
        q_phisic = (sph_O2 * temp_O2 * self.volume_O2 +
                    self.spHv_CO2_pol(self.temp_recyc) * self.temp_recyc * self.volume_CO2_recyc)
        
        temp_resedual = abs((self.fuel.calorific +  q_phisic)/ 
                (sph_CO2 * self.volume_CO2 + sph_H2O * self.volume_H2O) - temp)
        
        return temp_resedual

        #РАСЧЕТ ТЕМПЕРАТУРЫ ГОРЕНИЯ
    def _burnTemp(self, k_recyc, temp_recyc): # расчет температуры (жаропроизводительность)
        self.temp_recyc = temp_recyc 
        self._massFraction(k_recyc=k_recyc) # расчет фракционного состава продуктов сгорания
        
        sol = root(self._tempFunc, 1000)
        # self.g['temp_gas'] = sol.x[0]
        return sol.x[0]

        # РАСЧЕТ К-АЛЬФА ПО ТЕПМЕРАТУРЕ ГОРЕНИЯ
    def burnAlpha(self, temp_gas, temp_recyc): # расчет k_alpha
        temp_gas = temp_gas - self.temp0
        temp_recyc = temp_recyc - self.temp0
    
        resedual = lambda k_recyc: abs(self._burnTemp(k_recyc=k_recyc, 
                temp_recyc=temp_recyc) - temp_gas)
    
        sol = root(resedual, 1)
        return sol.x[0]
    
        # РАСЧЕТ ЗАВИСИМОСТИ И ПОСТОРЕНИЕ ГРАФИКА
    def tempAlphaPl(self, k_recyc, temp_recyc):
        plt.figure(figsize=(5, 4))
        temp_recyc = temp_recyc - self.temp0        
        n = 50
        ti = np.zeros(n)
        ki = np.linspace(k_recyc[0], k_recyc[1], n)
        
        for i in range(n):
            ti[i] = self._burnTemp(k_recyc=ki[i], temp_recyc=temp_recyc) + self.temp0
        
        plt.plot(ki, ti, color='red')
        plt.minorticks_on()
        plt.xlabel('k_recyc') 
        plt.ylabel('T,K') 
        plt.grid(linestyle='--', linewidth=0.5, color='black') # сетка
        plt.tight_layout() # оптимизируем поля и расположение объектов
        plt.savefig('temp-k_recyc.jpg', dpi = 300)
        plt.show()











