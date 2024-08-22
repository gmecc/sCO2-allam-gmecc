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
                                       'dens', 'entr', 'enth', 'dh', 'sp_heat', 'efc', 'phase', 'mfr'], index=range(8))
        self.p.loc[:, 'dp_rel'] = [.98, np.nan, .95, .99, 1., .95, np.nan, .95]  # относительные потери давления
        self.p.loc[:, 'efc'] = (.99, .9, .85, 1., 1., 1., .86, .85)  # КПД
        self.p.loc[:, 'dt'] = (.0, .0, 5., .0, .0, .0, .0, 5.)  # температурный напор [град]
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
    def cycle(self, pressure_min: object, pressure_rate: object, temperature: object, power: object, pinch_point: object = 5) -> object:
        
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

        # точки 3 - за сепаратором / осушение [CO2]

        self.p.loc[3, 'temp'] = self.p.temp[2]
        self.p.loc[3, 'enth'] = PropsSI('H', 'T', self.p.temp[3], 'P', self.p.pres[3], 'CO2')
        self.p.loc[3, 'entr'] = PropsSI('S', 'T', self.p.temp[3], 'P', self.p.pres[3], 'CO2')
        self.p.loc[3, 'dens'] = PropsSI('D', 'T', self.p.temp[3], 'P', self.p.pres[3], 'CO2')
        self.p.loc[3, 'sp_heat'] = PropsSI('C', 'T', self.p.temp[3], 'P', self.p.pres[3], 'CO2')
        self.p.loc[3, 'dt'] = self.p.temp[3] - self.p.temp[2]
        self.p.loc[3, 'dh'] = self.p.enth[3] - self.p.enth[2]
        self.p.loc[3, 'phase'] = PropsSI('Phase', 'T', self.p.temp[3], 'P', self.p.pres[3], 'CO2')

        # точка 4 - за сепаратором после отвода избытка СО2 [CO2]
        self.p.loc[4, 'temp'] = self.p.temp[3]
        self.p.loc[4, 'enth'] = PropsSI('H', 'T', self.p.temp[4], 'P', self.p.pres[4], 'CO2')
        self.p.loc[4, 'entr'] = PropsSI('S', 'T', self.p.temp[4], 'P', self.p.pres[4], 'CO2')
        self.p.loc[4, 'dens'] = PropsSI('D', 'T', self.p.temp[4], 'P', self.p.pres[4], 'CO2')
        self.p.loc[4, 'sp_heat'] = PropsSI('C', 'T', self.p.temp[4], 'P', self.p.pres[4], 'CO2')

        self.p.loc[4, 'dt'] = self.p.temp[4] - self.p.temp[3]
        self.p.loc[5, 'dt'] = self.p.temp[5] - self.p.temp[4]
        self.p.loc[4, 'dh'] = self.p.enth[4] - self.p.enth[3]
        self.p.loc[4, 'phase'] = PropsSI('Phase', 'T', self.p.temp[4], 'P', self.p.pres[4], 'CO2')

        #  расход

        self.power = power * 1000.

        self.p.loc[0, 'mfr'] = self.power / abs(self.p.dh[5] + self.p.dh[2])
        self.p.loc[1, 'mfr'] = self.p.mfr[0]
        self.p.loc[2, 'mfr'] = self.p.mfr[1]
        self.p.loc[3, 'mfr'] = self.p.mfr[2] * self.comb.gas.mol['CO2']
        self.p.loc[4, 'mfr'] = self.p.mfr[3] * self.comb.gas.mol['CO2_recyc']
        self.p.loc[5, 'mfr'] = self.p.mfr[4]
        self.p.loc[6, 'mfr'] = self.p.mfr[5]
        self.p.loc[7, 'mfr'] = self.p.mfr[6]

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
        
        
    def efc_pinch(self, var_temp_recyc, pressure_min, pressure_rate, temperature,power, pinch_point=5):
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
        entrop = np.vectorize(lambda p, t: PropsSI('S', 'T', t, 'P', p, 'CO2'))
        entrop_mix = np.vectorize(lambda p, t: PropsSI('S', 'T|supercritical', t, 'P', p, self.fluid_mix))

        t = np.linspace(self.p.loc[0, "temp"], self.p.loc[1, "temp"])
        e = np.linspace(self.p.loc[0, "entr"], self.p.loc[1, "entr"])
        plt.plot(e, t, 'r', label='0-1 - адиабатное расширение в турбине')

        p = np.linspace(self.p.loc[1, "pres"], self.p.loc[2, "pres"])
        t = np.linspace(self.p.loc[1, "temp"], self.p.loc[2, "temp"])
        plt.plot(entrop_mix(p, t), t, 'r', label='1-2 - охлаждение в рекуператоре')

        t = np.linspace(self.p.loc[2, "temp"], self.p.loc[3, "temp"])
        # p = np.linspace(self.p.loc[2, "pres"], self.p.loc[3, "pres"])
        # plt.plot(entrop_mix(p,t), t, 'r', label='2-3 - осушение в сепараторе')
        e = np.linspace(self.p.loc[2, "entr"], self.p.loc[3, "entr"])
        plt.plot(e, t, 'r', label='2-3 - осушение в сепараторе')

        t = np.linspace(self.p.loc[3, "temp"], self.p.loc[4, "temp"])
        e = np.linspace(self.p.loc[3, "entr"], self.p.loc[4, "entr"])
        plt.plot(e, t, 'b', label='3-4 - отвод диоксида углерода')

        t = np.linspace(self.p.loc[4, "temp"], self.p.loc[5, "temp"])
        e = np.linspace(self.p.loc[4, "entr"], self.p.loc[5, "entr"])
        plt.plot(e, t, 'b', label='4-5 - нагревание в теплообменнике')

        t = np.linspace(self.p.loc[5, "temp"], self.p.loc[6, "temp"])
        e = np.linspace(self.p.loc[5, "entr"], self.p.loc[6, "entr"])
        plt.plot(e, t, 'b', label='5-6 - сжатие в компрессоре')

        p = np.linspace(self.p.loc[6, "pres"], self.p.loc[7, "pres"])
        t = np.linspace(self.p.loc[6, "temp"], self.p.loc[7, "temp"])
        plt.plot(entrop(p, t), t, 'b', label='6-7 - нагревание в рекуператоре')

        t = np.linspace(self.p.loc[7, "temp"], self.p.loc[0, "temp"])
        e = np.linspace(self.p.loc[7, "entr"], self.p.loc[0, "entr"])
        plt.plot(e, t, 'b', label='7-0 - нагрев в камере сгорания')

        entr = np.delete(self.p.entr, (1, 7))
        t = np.delete(self.p.temp, (1, 7))
        plt.scatter(entr, t, color='b')  # точка

        plt.text(entrop(self.p.loc[0, "pres"], self.p.loc[0, "temp"]), self.p.loc[0, "temp"], "0")
        plt.text(entrop(self.p.loc[1, "pres"], self.p.loc[1, "temp"]), self.p.loc[1, "temp"], "1")
        plt.text(entrop(self.p.loc[2, "pres"], self.p.loc[2, "temp"]) + 60, self.p.loc[2, "temp"] - 30, "2")
        plt.text(entrop(self.p.loc[3, "pres"], self.p.loc[3, "temp"]), self.p.loc[3, "temp"] - 30, "3")
        plt.text(entrop(self.p.loc[4, "pres"], self.p.loc[4, "temp"]), self.p.loc[4, "temp"] + 20, "4")
        plt.text(entrop(self.p.loc[5, "pres"], self.p.loc[5, "temp"]) - 50, self.p.loc[5, "temp"] + 50, "5")
        plt.text(entrop(self.p.loc[6, "pres"], self.p.loc[6, "temp"]) - 50, self.p.loc[6, "temp"] - 70, "6")
        plt.text(entrop(self.p.loc[7, "pres"], self.p.loc[7, "temp"]), self.p.loc[7, "temp"], "7")

        plt.minorticks_on()
        plt.ylabel('$T^\circ C$', fontsize=10)
        plt.xlabel('$S, Дж/кг$', fontsize=10)
        plt.grid(linestyle='--', linewidth=0.5, color='black')  # сетка
        plt.legend()
        plt.tight_layout()  # оптимизируем поля и расположение объектов
        file_name = 'ts-diagramm.png'
        plt.savefig(file_name, dpi=300)
        plt.show()
        print(file_name)

    def ph_diagram(self):

        p = np.linspace(self.p.loc[0, "pres"], self.p.loc[1, "pres"])
        h = np.linspace(self.p.loc[0, "enth"], self.p.loc[1, "enth"])
        plt.plot(h, p, 'r', label='0-1 - адиабатное расширение в турбине')

        p = np.linspace(self.p.loc[1, "pres"], self.p.loc[2, "pres"])
        h = np.linspace(self.p.loc[1, "enth"], self.p.loc[2, "enth"])
        plt.plot(h, p, 'r', label='1-2 - охлаждение в рекуператоре')

        p = np.linspace(self.p.loc[2, "pres"], self.p.loc[3, "pres"])
        h = np.linspace(self.p.loc[2, "enth"], self.p.loc[3, "enth"])
        plt.plot(h, p, 'b', label='2-3 - осушение в сепараторе')

        p = np.linspace(self.p.loc[3, "pres"], self.p.loc[4, "pres"])
        h = np.linspace(self.p.loc[3, "enth"], self.p.loc[4, "enth"])
        plt.plot(h, p, 'b', label='3-4 - отвод диоксида углерода')

        p = np.linspace(self.p.loc[4, "pres"], self.p.loc[5, "pres"])
        h = np.linspace(self.p.loc[4, "enth"], self.p.loc[5, "enth"])
        plt.plot(h, p, 'b', label='4-5 - нагревание в теплообменнике')

        p = np.linspace(self.p.loc[5, "pres"], self.p.loc[6, "pres"])
        h = np.linspace(self.p.loc[5, "enth"], self.p.loc[6, "enth"])
        plt.plot(h, p, 'b', label='5-6 - адиабатное сжатие в компрессоре')

        p = np.linspace(self.p.loc[6, "pres"], self.p.loc[7, "pres"])
        h = np.linspace(self.p.loc[6, "enth"], self.p.loc[7, "enth"])
        plt.plot(h, p, 'b', label='6-7 - нагревание в рекуператоре')

        p = np.linspace(self.p.loc[7, "pres"], self.p.loc[0, "pres"])
        h = np.linspace(self.p.loc[7, "enth"], self.p.loc[0, "enth"])
        plt.plot(h, p, 'b', label='7-0 - нагрев в камере сгорания')

        entalp = np.delete(self.p.enth, (1, 7))
        p = np.delete(self.p.pres, (1, 7))
        plt.scatter(entalp, p, color='b')  # точка

        plt.minorticks_on()
        plt.ylabel('$P, Па$', fontsize=10)
        plt.xlabel('$H, Дж/моль$', fontsize=10)
        plt.grid(linestyle='--', linewidth=0.5, color='black')  # сетка
        plt.legend()
        plt.tight_layout()  # оптимизируем поля и расположение объектов
        file_name = 'ph-diagramm.png'
        plt.savefig(file_name, dpi=600)
        plt.show()

    def hs_diagram(self):


        s = np.linspace(self.p.loc[0, "entr"], self.p.loc[1, "entr"])
        h = np.linspace(self.p.loc[0, "enth"], self.p.loc[1, "enth"])
        plt.plot(s, h, 'r', label='0-1 - адиабатное расширение в турбине')

        s = np.linspace(self.p.loc[1, "entr"], self.p.loc[2, "entr"])
        h = np.linspace(self.p.loc[1, "enth"], self.p.loc[2, "enth"])
        plt.plot(s, h, 'r', label='1-2 - охлаждение в рекуператоре')

        s = np.linspace(self.p.loc[2, "entr"], self.p.loc[3, "entr"])
        h = np.linspace(self.p.loc[2, "enth"], self.p.loc[3, "enth"])
        plt.plot(s, h, 'r', label='2-3 - осушение в сепараторе')

        s = np.linspace(self.p.loc[3, "entr"], self.p.loc[4, "entr"])
        h = np.linspace(self.p.loc[3, "enth"], self.p.loc[4, "enth"])
        plt.plot(s, h, 'b', label='3-4 - отвод диоксида углерода')

        s = np.linspace(self.p.loc[4, "entr"], self.p.loc[5, "entr"])
        h = np.linspace(self.p.loc[4, "enth"], self.p.loc[5, "enth"])
        plt.plot(s, h, 'b', label='4-5 - нагревание в теплообменнике')

        s = np.linspace(self.p.loc[5, "entr"], self.p.loc[6, "entr"])
        h = np.linspace(self.p.loc[5, "enth"], self.p.loc[6, "enth"])
        plt.plot(s, h, 'b', label='5-6 - адиабатное сжатие в компрессоре')

        s = np.linspace(self.p.loc[6, "entr"], self.p.loc[7, "entr"])
        h = np.linspace(self.p.loc[6, "enth"], self.p.loc[7, "enth"])
        plt.plot(s, h, 'b', label='6-7 - нагревание в рекуператоре')

        s = np.linspace(self.p.loc[7, "entr"], self.p.loc[0, "entr"])
        h = np.linspace(self.p.loc[7, "enth"], self.p.loc[0, "enth"])
        plt.plot(s, h, 'b', label='7-0 - нагрев в камере сгорания')

        entalp = np.delete(self.p.enth, (1, 7))
        entr = np.delete(self.p.entr, (1, 7))
        plt.scatter(entr, entalp, color='b')  # точка

        plt.minorticks_on()
        plt.ylabel('$S, Дж/кг$', fontsize=10)
        plt.xlabel('$H, Дж/моль$', fontsize=10)
        plt.grid(linestyle='--', linewidth=0.5, color='black')  # сетка
        # plt.legend()
        plt.tight_layout()  # оптимизируем поля и расположение объектов
        file_name = 'hs-diagramm.png'
        plt.savefig(file_name, dpi=600)
        plt.show()
'''        
# РАСЧЕТ ПАРАМЕТРОВ ЦИКЛА
temperature[0] # температура перед компрессором [K]
temperature[1] # температура перед камерой сгорания [K]
temperature[2] # температура перед турбиной [K]
'''







