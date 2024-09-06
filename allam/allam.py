'''
ТЕПЛОВОЙ РАСЧЕТ ЦИКЛА АЛЛАМА

Входные данные:
    - давление перед компрессором [Pa]
    - степень сжатия
    - температура (кортеж):
        температура перед комперссором [K]
        температура перед турбиной [K]
        температура рециркулируемого СО2 перед камерой сгорания [K]
    
Выходные данные:
    - .p - (dataframe) параметры в точках
    - .g - (serial) параметры цикла
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
from allam import Combust
from scipy.optimize import root

pd.set_option('display.float_format', '{:.2f}'.format)
# pd.set_option('display.max_columns', None)

class Acycle:
    def __init__(self):
        self.p = pd.DataFrame(columns=['CO2', 'H2O', 'temp', 'dt', 'pres', 'dp_rel',
                'dens', 'entr', 'enth', 'dh', 'sp_heat', 'efc', 'phase'], index=range(8))
        self.p.loc[:, 'dp_rel'] = [.98, np.nan, .95, .99, 1., .95, np.nan, .95] # относительные потери давления
        self.p.loc[:, 'efc'] = (.99,  .9, .85, 1., 1., 1., .86, .85) # КПД
        self.p.loc[:, 'dt'] = (.0, .0, 5., .0, .0, .0, .0, 5.) # температурный напор [град]
        self.g = pd.Series([np.nan], index=['k_recyc'])


        # РАСЧЕТ ТЕМПЕРАТУРЫ ПО ЭНТРОПИИ
    def _tempS(self, pressure, entropy): 
        f = lambda temperature: PropsSI('S','T|supercritical',temperature,'P',pressure,self.fluid_mix)
        def func_resedual(x, param):
            resedual = f(x) - param
            return resedual
        sol = root(func_resedual, 500, args=(entropy))
        return sol.x[0]
    
    
        # РАСЧЕТ ТЕМПЕРАТУРЫ ПО ЭНТАЛЬПИИ
    def _tempH(self, pressure, enthalpy):
        f = lambda temperature: PropsSI('H','T|supercritical',temperature,'P',pressure,self.fluid_mix)
        def func_resedual(x, param):
            resedual = f(x) - param
            return resedual
        sol = root(func_resedual, 1000, args=(enthalpy))
        return sol.x[0]

        # РАСЧЕТ ПАРАМЕТРОВ РТ В ТОЧКАХ ЦИКЛА
    def cycle(self, pressure_min, pressure_rate, temperature, pinch_point=5):
        
        self.pressure_min = pressure_min # давление перед компрессором
        self.pressure_rate = pressure_rate # повышение давления в компрессоре
        self.p.loc[5, 'temp'] = temperature[0] # температура перед компрессором
        self.p.loc[0, 'temp'] = temperature[1] # температура перед турбиной
        self.p.loc[7, 'temp'] = temperature[2] # температура СО2 перед камерой сгорания
        self.pinch_point = pinch_point # пинч-поинт
        

        # давление в контуре
        self.p.loc[5, 'pres'] = self.pressure_min
        self.p.loc[6, 'pres'] = self.p.pres[5] * self.pressure_rate
        self.p.loc[7, 'pres'] = self.p.pres[6] * self.p.dp_rel[7]
        self.p.loc[0, 'pres'] = self.p.pres[7] * self.p.dp_rel[0]
        self.p.loc[4, 'pres'] = self.p.pres[5] / self.p.dp_rel[5]
        self.p.loc[3, 'pres'] = self.p.pres[4] / self.p.dp_rel[4]
        self.p.loc[2, 'pres'] = self.p.pres[3] / self.p.dp_rel[3]
        self.p.loc[1, 'pres'] = self.p.pres[2] / self.p.dp_rel[2]
        
        # фракционный состав рабочего тела
        self.comb = Combust()
        self.g.at['k_recyc'] = self.comb.burnAlpha(temp_gas=self.p.temp[0], 
                    temp_recyc=self.p.temp[7])
        
        molCO2 = round(self.comb.gas.mol['CO2'], 3)
        molH2O = round((1 - molCO2), 3)
        
        # print(molCO2)
        
        self.fluid_mix = 'CO2[' + str(molCO2) + ']&water[' + str(molH2O) + ']'
        
        self.p.loc[0:2, 'CO2'] = molCO2
        self.p.loc[0:2, 'H2O'] = 1 - molCO2
        
        self.p.loc[3:7, 'CO2'] = 1
        self.p.loc[3:7, 'H2O'] = 0
        
        
        # точка 5 - перед компрессором [CO2]
        self.p.loc[5, 'enth'] = PropsSI('H','T',self.p.temp[5],'P',self.p.pres[5],'CO2')
        self.p.loc[5, 'entr'] = PropsSI('S','T',self.p.temp[5],'P',self.p.pres[5],'CO2')
        self.p.loc[5, 'dens'] = PropsSI('D','T',self.p.temp[5],'P',self.p.pres[5],'CO2')
        self.p.loc[5, 'sp_heat'] = PropsSI('C','T',self.p.temp[5],'P',self.p.pres[5],'CO2')
        self.p.loc[5, 'phase'] = PropsSI('Phase','T',self.p.temp[5],'P',self.p.pres[5],'CO2')
        
        # точка 6 - за компрессором / адиабатическое сжатие [CO2]
        enthalpy_isoentr_compr = PropsSI('H','S',self.p.entr[5],'P',self.p.pres[6],'CO2')
        dh_isoentr_compr =  enthalpy_isoentr_compr - self.p.enth[5]
        
        self.p.loc[6, 'enth'] = self.p.enth[5] + dh_isoentr_compr / self.p.efc[6]
        self.p.loc[6, 'temp'] = PropsSI('T','H',self.p.enth[6],'P',self.p.pres[6],'CO2')
        self.p.loc[6, 'entr'] = PropsSI('S','T',self.p.temp[6],'P',self.p.pres[6],'CO2')
        self.p.loc[6, 'dens'] = PropsSI('D','T',self.p.temp[6],'P',self.p.pres[6],'CO2')
        self.p.loc[6, 'sp_heat'] = PropsSI('C','T',self.p.temp[6],'P',self.p.pres[6],'CO2')
        self.p.loc[6, 'dh'] = self.p.enth[6] - self.p.enth[5]  # перепад энтальпий в компрессоре при изоэнтропическом сжатии
        self.p.loc[6, 'dt'] = self.p.temp[6] - self.p.temp[5]
        self.p.loc[6, 'dp_rel'] = self.p.pres[5] / self.p.pres[6]
        self.p.loc[6, 'phase'] = PropsSI('Phase','T',self.p.temp[6],'P',self.p.pres[6],'CO2')
        
        self.p.loc[0, 'enth'] = PropsSI('H','T|supercritical',self.p.temp[0],'P',self.p.pres[0],self.fluid_mix)
        self.p.loc[0, 'entr'] = PropsSI('S','T|supercritical',self.p.temp[0],'P',self.p.pres[0],self.fluid_mix)
        self.p.loc[0, 'dens'] = PropsSI('D','T|supercritical',self.p.temp[0],'P',self.p.pres[0],self.fluid_mix)
        self.p.loc[0, 'sp_heat'] = PropsSI('C','T|supercritical',self.p.temp[0],'P',self.p.pres[0],self.fluid_mix)
        
        self.p.loc[0, 'dt'] = self.p.temp[0] - self.p.temp[2]
        self.p.loc[0, 'dh'] = self.p.enth[0] - self.p.enth[2]
        self.p.loc[0, 'phase'] = PropsSI('Phase','T|supercritical',self.p.temp[0],'P',self.p.pres[0],self.fluid_mix)
        
        
        # точка 1 - за турбиной / адиабатическое расширение [CO2,H2O]
        temp_isoentr_expand = self._tempS(self.p.pres[1], self.p.entr[0])
        enthalpy_isoentr_expand = PropsSI('H','T|supercritical',temp_isoentr_expand,'P',self.p.pres[1],self.fluid_mix)
        dh_isoentr_expand =  enthalpy_isoentr_expand - self.p.enth[0]
        self.p.loc[1, 'enth'] = self.p.enth[0] + dh_isoentr_expand * self.p.efc[1]
        self.p.loc[1, 'temp'] = self._tempH(self.p.pres[1], self.p.enth[1])
        
        self.p.loc[1, 'entr'] = PropsSI('S','T|supercritical',self.p.temp[1],'P',self.p.pres[1],self.fluid_mix)
        self.p.loc[1, 'dens'] = PropsSI('D','T|supercritical',self.p.temp[1],'P',self.p.pres[1],self.fluid_mix)
        self.p.loc[1, 'sp_heat'] = PropsSI('C','T|supercritical',self.p.temp[1],'P',self.p.pres[1],self.fluid_mix)
        
        self.p.loc[1, 'dt'] = self.p.temp[1] - self.p.temp[0]
        self.p.loc[1, 'dh'] = self.p.enth[1] - self.p.enth[0]
        self.p.loc[1, 'dp_rel'] = self.p.pres[0] / self.p.pres[1]
        self.p.loc[1, 'phase'] = PropsSI('Phase','T|supercritical',self.p.temp[1],'P',self.p.pres[1],self.fluid_mix)
        ###################################################
        
        # точка 7 - за рекуператором / холодная часть - нагрев [CO2]
        self.p.loc[7, 'enth'] = PropsSI('H','T',self.p.temp[7],'P',self.p.pres[7],'CO2')
        self.p.loc[7, 'entr'] = PropsSI('S','T',self.p.temp[7],'P',self.p.pres[7],'CO2')
        self.p.loc[7, 'dens'] = PropsSI('D','T',self.p.temp[7],'P',self.p.pres[7],'CO2')
        self.p.loc[7, 'sp_heat'] = PropsSI('C','T',self.p.temp[7],'P',self.p.pres[7],'CO2')
        
        self.p.loc[7, 'dh'] = self.p.enth[7] - self.p.enth[6]
        self.p.loc[7, 'dt'] = self.p.temp[7] - self.p.temp[6]
        self.p.loc[7, 'phase'] = PropsSI('Phase','T',self.p.temp[7],'P',self.p.pres[7],'CO2')
        
        
        # точка 2 - за рекуператором / горячая часть - охлаждение [CO2,H2O]
        
        self.p.loc[2, 'dh'] = - self.comb.gas.mass['CO2_recyc'] * self.p.dh[7]
        self.p.loc[2, 'enth'] = self.p.enth[1] + self.p.dh[2]
        self.p.loc[2, 'temp'] = self._tempH(pressure=self.p.pres[2], enthalpy=self.p.enth[2])
        
        self.p.loc[2, 'entr'] = PropsSI('S','T|supercritical',self.p.temp[2],'P',self.p.pres[2],self.fluid_mix)
        self.p.loc[2, 'dens'] = PropsSI('D','T|supercritical',self.p.temp[2],'P',self.p.pres[2],self.fluid_mix)
        self.p.loc[2, 'sp_heat'] = PropsSI('C','T|supercritical',self.p.temp[2],'P',self.p.pres[2],self.fluid_mix)
        
        self.p.loc[2, 'dp_rel'] = self.p.pres[1] / self.p.pres[2]
        self.p.loc[2, 'dt'] = self.p.temp[2] - self.p.temp[1]
        self.p.loc[2, 'phase'] = PropsSI('Phase','T|supercritical',self.p.temp[2],'P',self.p.pres[2],self.fluid_mix)
        
        ######################################################
        
        # точки 3, 4 - за сепаратором / осушение [CO2]
        
        self.p.loc[3:4, 'temp'] = self.p.temp[2]
        self.p.loc[3:4, 'enth'] = PropsSI('H','T',self.p.temp[3],'P',self.p.pres[3],'CO2')
        self.p.loc[3:4, 'entr'] = PropsSI('S','T',self.p.temp[3],'P',self.p.pres[3],'CO2')
        self.p.loc[3:4, 'dens'] = PropsSI('D','T',self.p.temp[3],'P',self.p.pres[3],'CO2')
        self.p.loc[3:4, 'sp_heat'] = PropsSI('C','T',self.p.temp[3],'P',self.p.pres[3],'CO2')
        
        self.p.loc[3, 'dt'] = self.p.temp[3] - self.p.temp[2]
        self.p.loc[3, 'dh'] = self.p.enth[3] - self.p.enth[2]
        
        self.p.loc[4, 'dt'] = self.p.temp[4] - self.p.temp[3]
        self.p.loc[5, 'dt'] = self.p.temp[5] - self.p.temp[4]
        self.p.loc[4, 'dh'] = self.p.enth[4] - self.p.enth[3]
        self.p.loc[3:4, 'phase'] = PropsSI('Phase','T',self.p.temp[3],'P',self.p.pres[3],'CO2')
        
        
        # точка 0 - за охладителем / охлаждение [CO2,H2O]
        self.p.loc[5, 'dh'] = self.p.enth[5] - self.p.enth[4]
        self.p.loc[0, 'dh'] = self.p.enth[0] - self.p.enth[7] # перепад энтальпий в охладителе
        
        self.p.loc[0, 'dt'] = self.p.temp[0] - self.p.temp[7]
        self.p.loc[0, 'phase'] = PropsSI('Phase','T|supercritical',self.p.temp[0],'P',self.p.pres[0],self.fluid_mix)
        
        # пинч-поинт
        self.g.at['pinch'] = self.p.temp[1] - self.p.temp[7]
        
        # полезная работа цикла
        self.g.at['work_cycle'] = abs(self.p.dh[1] - self.p.dh[6] * self.comb.gas.mass['CO2_recyc'])
        
        # эффективность    
        self.g.at['efc_cycle'] = self.g.work_cycle / self.p.dh[0]

        self.p.to_csv('cycle.csv', index=False)
        self.g.to_csv('cycle_g.csv', index=False)

    def efc_temp_recyc(self, var_temp_recyc, pressure_min, pressure_rate, temperature, pinch_point=5):
        plt.figure(figsize=(5, 4))
        # var_temp - кортеж значений диапазона изменения параметра
        
        temp = np.linspace(var_temp_recyc[0], var_temp_recyc[1], 50)
        efc = np.zeros(50)
        
        for i in range(50):
            self.cycle(pressure_min, pressure_rate, 
                       temperature=(temperature[0], temperature[1], temp[i]), 
                       pinch_point=5)
            efc[i] = self.g.efc_cycle
        
        pol = np.polynomial.Polynomial.fit(temp, efc, deg=3)
        plt.plot(temp, pol(temp))
        plt.minorticks_on()
        plt.ylabel('efc') 
        plt.xlabel('Temperature [K]') 
        plt.grid(linestyle='--', linewidth=0.5, color='black') # сетка
        plt.savefig('efc_temp_recyc.jpg', dpi=300)
        plt.show()
        
        
    def efc_pinch(self, var_temp_recyc, pressure_min, pressure_rate, temperature, pinch_point=5):
        plt.figure(figsize=(5, 4))
        # var_temp - кортеж значений диапазона изменения параметра
        
        temp = np.linspace(var_temp_recyc[0], var_temp_recyc[1], 50)
        pinch = np.zeros(50)
        
        for i in range(50):
            self.cycle(pressure_min, pressure_rate, 
                       temperature=(temperature[0], temperature[1], temp[i]), 
                       pinch_point=5)
            pinch[i] = self.g.pinch
    
        pol = np.polynomial.Polynomial.fit(temp, pinch, deg=3)
        plt.plot(temp, pol(temp))
        plt.minorticks_on()
        plt.ylabel('pinch [K]') 
        plt.xlabel('Temperature recyc [K]') 
        plt.grid(linestyle='--', linewidth=0.5, color='black') # сетка
        plt.show()
        
        
    def efc_temp_heat(self, var_temp_heat, pressure_min, pressure_rate, temperature, pinch_point=5):
        plt.figure(figsize=(5, 4))
        # var_temp - кортеж значений диапазона изменения параметра
        
        temp = np.linspace(var_temp_heat[0], var_temp_heat[1], 50)
        efc = np.zeros(50)
        
        for i in range(50):
            self.cycle(pressure_min, pressure_rate, 
                       temperature=(temperature[0], temp[i], temperature[1]), 
                       pinch_point=5)
            efc[i] = self.g.efc_cycle
        
        pol = np.polynomial.Polynomial.fit(temp, efc, deg=3)
        plt.plot(temp, pol(temp))
        plt.minorticks_on()
        plt.ylabel('efc') 
        plt.xlabel('Temperature [K]') 
        plt.grid(linestyle='--', linewidth=0.5, color='black') # сетка
        plt.savefig('efc_temp_heat.jpg', dpi=300)
        plt.show()
        
        
        
        
        
        
    #####################################################################################
    def power(self, power):
        self.power = power * 1000.
        # self.p['pwr'] = np.nan
        
        self.mfr = self.power / abs(self.p.dh[5] + self.p.dh[2])
        self.mfp = self.mfr * (self.p.temp[4] + 263.15)**.5 / self.p.pres[4] * 1e5
        
        self.p.pwr[0] = self.p.dh[0] * self.mfr # мощность нагревателя
        self.p.pwr[2] = self.p.dh[2] * self.mfr # мощность компрессора
        self.p.pwr[3] = self.p.dh[3] * self.mfr # тепловая мощность экономайзера нагрев
        self.p.pwr[4] = self.p.dh[4] * self.mfr # тепловая мощность нагревателя
        self.p.pwr[5] = self.p.dh[5] * self.mfr # мощность турбины
        self.p.pwr[6] = self.p.dh[6] * self.mfr # тепловая мощность экономайзера охл

    def optim(self, pressure, temperature, pressure_rate):
        press_rate = np.linspace(pressure_rate[0], pressure_rate[1], 50)
        efc_rate = np.zeros(50)
        plt.figure(figsize=(5,5), dpi=300)  
        
        for temp in temperature:
            for i in range(50):
                efc_rate[i] = self.system(pressure_min_rel=pressure,  pressure_rate=press_rate[i], 
                                              temperature=temp, 
                                              pinch_point=10)
            plt.plot(press_rate, efc_rate, label='T '+str(temp))
        
        plt.minorticks_on()
        plt.ylabel('КПД', fontsize=10) 
        plt.xlabel('$\pi_k$', fontsize=10) 
        plt.grid(linestyle='--', linewidth=0.5, color='black') # сетка
        plt.legend()
        plt.tight_layout() # оптимизируем поля и расположение объектов
        file_name = 'optim.png'
        plt.savefig(file_name, dpi = 300)
        plt.show()
        # print(file_name)
        
        
    def specific_speed(self, n):
        self.n = n / 50 # частота вращения ротора турбины
        self.sp_speed = ((2 * np.pi * self.n * (self.mfr / self.p.dens[5])**.5) /
                         (abs(self.p.dh[5])**(3/4)))
        return self.sp_speed        
        
    def speed_optim(self):
        self.n_opt = (0.548 * (abs(self.p.dh[5])**(3/4)) /
                (2 * np.pi * (self.mfr / self.p.dens[5])**.5)  * 50)
        return self.n_opt
        
    
    def ts_diagramm(self):
        # линия насыщения
        temp = np.linspace(-31, self.fluid.critical_temperature)
        entropy_dew = np.vectorize(lambda x: self.fluid.dew_point_at_temperature(x).entropy)
        entropy_bubble = np.vectorize(lambda x: self.fluid.bubble_point_at_temperature(x).entropy)
        
        plt.plot(entropy_dew(temp), temp, 'grey')
        plt.plot(entropy_bubble(temp), temp, 'grey')
        
        entrop = np.vectorize(lambda p, t: self.fluid.with_state(Input.pressure(p), 
                                Input.temperature(t)).entropy)
           
        p = np.linspace(self.p.pres[2], self.p.pres[4])
        t = np.linspace(self.p.temp[2], self.p.temp[4])
        plt.plot(entrop(p,t), t, 'r')
        
        p = np.linspace(self.p.pres[5], self.p.pres[0])
        t = np.linspace(self.p.temp[5], self.p.temp[0])
        plt.plot(entrop(p,t), t, 'b')
        
        p = np.linspace(self.p.pres[4], self.p.pres[5])
        t = np.linspace(self.p.temp[4], self.p.temp[5])
        plt.plot(entrop(p,t), t, 'b')
        
        p = np.linspace(self.p.pres[0], self.p.pres[2])
        t = np.linspace(self.p.temp[0], self.p.temp[2])
        plt.plot(entrop(p,t), t, 'b')
        
        entr = np.delete(self.p.entr, (1,5))
        t = np.delete(self.p.temp, (1,5))
        plt.scatter(entr, t, color='b') # точка
        
        plt.minorticks_on()
        plt.ylabel('$T^\circ C$', fontsize=10) 
        plt.xlabel('$S,Дж/кг$', fontsize=10)    
        plt.grid(linestyle='--', linewidth=0.5, color='black') # сетка
        # plt.legend()
        plt.tight_layout() # оптимизируем поля и расположение объектов
        file_name = 'ts-diagramm.png'
        plt.savefig(file_name, dpi = 300)
        plt.show()
        print(file_name)

'''        
# РАСЧЕТ ПАРАМЕТРОВ ЦИКЛА
temperature[0] # температура перед компрессором [K]
temperature[1] # температура перед камерой сгорания [K]
temperature[2] # температура перед турбиной [K]
'''







