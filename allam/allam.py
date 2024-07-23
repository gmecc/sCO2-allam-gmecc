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
from allam.allam import Combust
from scipy.optimize import root

pd.set_option('display.float_format', '{:.2f}'.format)
# pd.set_option('display.max_columns', None)

class Acycle:
    def __init__(self):
        self.p = pd.DataFrame(columns=['CO2', 'H2O', 'temp', 'dt', 'pres', 'dp_rel',
                'dens', 'entr', 'enth', 'dh', 'sp_heat', 'efc', 'phase'], index=range(8))
        self.p.dp_rel[:] = (.98, np.nan, .95, .99, 1., .95, np.nan, .95) # относительные потери давления
        self.p.efc[:] = (.99,  .9, .85, 1., 1., 1., .86, .85) # КПД
        self.p.dt[:] = (.0, .0, 5., .0, .0, .0, .0, 5.) # температурный напор [град]
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
    
    
    # def _temp_cool_out(self):
    #     h = self.p.dh[2] / self.comb.gas.mass['CO2_recyc'] - self.p.enth[6]
    #     f = lambda temperature: PropsSI('H','T',temperature,'P',self.p.pres[5],self.fluid_mix)
        
    #     def func_resedual(x, param):
    #         resedual = f(x) - param
    #         return resedual
    #     sol = root(func_resedual, 1000, args=(h))
    #     return sol.x[0]
        
        
        
        # РАСЧЕТ ПАРАМЕТРОВ РТ В ТОЧКАХ ЦИКЛА
    def cycle(self, pressure_min, pressure_rate, temperature, pinch_point=5):
        
        self.pressure_min = pressure_min # давление перед компрессором
        self.pressure_rate = pressure_rate # повышение давления в компрессоре
        self.p.temp.iloc[5] = temperature[0] # температура перед компрессором
        self.p.temp.iloc[0] = temperature[1] # температура перед турбиной
        self.p.temp.iloc[7] = temperature[2] # температура СО2 перед камерой сгорания
        self.pinch_point = pinch_point # пинч-поинт
        

        # давление в контуре
        self.p.pres.iloc[5] = self.pressure_min 
        self.p.pres.iloc[6] = self.p.pres[5] * self.pressure_rate
        self.p.pres.iloc[7] = self.p.pres[6] * self.p.dp_rel[7]
        self.p.pres.iloc[0] = self.p.pres[7] * self.p.dp_rel[0]
        self.p.pres.iloc[4] = self.p.pres[5] / self.p.dp_rel[5] 
        self.p.pres.iloc[3] = self.p.pres[4] / self.p.dp_rel[4]
        self.p.pres.iloc[2] = self.p.pres[3] / self.p.dp_rel[3]  
        self.p.pres.iloc[1] = self.p.pres[2] / self.p.dp_rel[2]
        
        # фракционный состав рабочего тела
        self.comb = Combust()
        self.g.at['k_recyc'] = self.comb.burnAlpha(temp_gas=self.p.temp[0], 
                    temp_recyc=self.p.temp[7])
        
        molCO2 = round(self.comb.gas.mol['CO2'], 3)
        molH2O = round((1 - molCO2), 3)
        
        # print(molCO2)
        
        self.fluid_mix = 'CO2[' + str(molCO2) + ']&water[' + str(molH2O) + ']'
        
        self.p.CO2.iloc[0:3] = molCO2
        self.p.H2O.iloc[0:3] = 1 - molCO2
        
        self.p.CO2.iloc[3:8] = 1
        self.p.H2O.iloc[3:8] = 0
        
        
        # точка 5 - перед компрессором [CO2]
        self.p.enth.iloc[5] = PropsSI('H','T',self.p.temp[5],'P',self.p.pres[5],'CO2') 
        self.p.entr.iloc[5] = PropsSI('S','T',self.p.temp[5],'P',self.p.pres[5],'CO2')
        self.p.dens.iloc[5] = PropsSI('D','T',self.p.temp[5],'P',self.p.pres[5],'CO2')
        self.p.sp_heat.iloc[5] = PropsSI('C','T',self.p.temp[5],'P',self.p.pres[5],'CO2')
        self.p.phase.iloc[5] = PropsSI('Phase','T',self.p.temp[5],'P',self.p.pres[5],'CO2')
        
        # точка 6 - за компрессором / адиабатическое сжатие [CO2]
        enthalpy_isoentr_compr = PropsSI('H','S',self.p.entr[5],'P',self.p.pres[6],'CO2')
        dh_isoentr_compr =  enthalpy_isoentr_compr - self.p.enth[5]
        
        self.p.enth.iloc[6] = self.p.enth[5] + dh_isoentr_compr / self.p.efc[6]
        self.p.temp.iloc[6] = PropsSI('T','H',self.p.enth[6],'P',self.p.pres[6],'CO2') 
        self.p.entr.iloc[6] = PropsSI('S','T',self.p.temp[6],'P',self.p.pres[6],'CO2')
        self.p.dens.iloc[6] = PropsSI('D','T',self.p.temp[6],'P',self.p.pres[6],'CO2')
        self.p.sp_heat.iloc[6] = PropsSI('C','T',self.p.temp[6],'P',self.p.pres[6],'CO2')
        self.p.dh.iloc[6] = self.p.enth[6] - self.p.enth[5]  # перепад энтальпий в компрессоре при изоэнтропическом сжатии
        self.p.dt.iloc[6] = self.p.temp[6] - self.p.temp[5]
        self.p.dp_rel.iloc[6] = self.p.pres[5] / self.p.pres[6]
        self.p.phase.iloc[6] = PropsSI('Phase','T',self.p.temp[6],'P',self.p.pres[6],'CO2')
        
        ######################################
        # print(self.p)
        # print(self.g)
        
        
        # print("*** program OK ***")
        # sys.exit("*** STOP ***")
        ######################################
        
        # точка 0 - перед турбиной [CO2,H2O]
        # print(self.p.temp[0])
        # print(self.p.pres[0])
        
        
        self.p.enth.iloc[0] = PropsSI('H','T|supercritical',self.p.temp[0],'P',self.p.pres[0],self.fluid_mix) 
        self.p.entr.iloc[0] = PropsSI('S','T|supercritical',self.p.temp[0],'P',self.p.pres[0],self.fluid_mix)
        self.p.dens.iloc[0] = PropsSI('D','T|supercritical',self.p.temp[0],'P',self.p.pres[0],self.fluid_mix)
        self.p.sp_heat.iloc[0] = PropsSI('C','T|supercritical',self.p.temp[0],'P',self.p.pres[0],self.fluid_mix)
        
        self.p.dt.iloc[0] = self.p.temp[0] - self.p.temp[2]
        self.p.dh.iloc[0] = self.p.enth[0] - self.p.enth[2]  
        self.p.phase.iloc[0] = PropsSI('Phase','T|supercritical',self.p.temp[0],'P',self.p.pres[0],self.fluid_mix)
        
        
        # точка 1 - за турбиной / адиабатическое расширение [CO2,H2O]
        temp_isoentr_expand = self._tempS(self.p.pres[1], self.p.entr[0])
        enthalpy_isoentr_expand = PropsSI('H','T|supercritical',temp_isoentr_expand,'P',self.p.pres[1],self.fluid_mix)
        dh_isoentr_expand =  enthalpy_isoentr_expand - self.p.enth[0]
        self.p.enth.iloc[1] = self.p.enth[0] + dh_isoentr_expand * self.p.efc[1]
        self.p.temp.iloc[1] = self._tempH(self.p.pres[1], self.p.enth[1])
        
        self.p.entr.iloc[1] = PropsSI('S','T|supercritical',self.p.temp[1],'P',self.p.pres[1],self.fluid_mix)
        self.p.dens.iloc[1] = PropsSI('D','T|supercritical',self.p.temp[1],'P',self.p.pres[1],self.fluid_mix)
        self.p.sp_heat.iloc[1] = PropsSI('C','T|supercritical',self.p.temp[1],'P',self.p.pres[1],self.fluid_mix)
        
        self.p.dt.iloc[1] = self.p.temp[1] - self.p.temp[0]
        self.p.dh.iloc[1] = self.p.enth[1] - self.p.enth[0]  
        self.p.dp_rel.iloc[1] = self.p.pres[0] / self.p.pres[1]
        self.p.phase.iloc[1] = PropsSI('Phase','T|supercritical',self.p.temp[1],'P',self.p.pres[1],self.fluid_mix)
        
        
        ###################################################
        
        # точка 7 - за рекуператором / холодная часть - нагрев [CO2]
        self.p.enth.iloc[7] = PropsSI('H','T',self.p.temp[7],'P',self.p.pres[7],'CO2') 
        self.p.entr.iloc[7] = PropsSI('S','T',self.p.temp[7],'P',self.p.pres[7],'CO2')
        self.p.dens.iloc[7] = PropsSI('D','T',self.p.temp[7],'P',self.p.pres[7],'CO2')
        self.p.sp_heat.iloc[7] = PropsSI('C','T',self.p.temp[7],'P',self.p.pres[7],'CO2')
        
        self.p.dh.iloc[7] = self.p.enth[7] - self.p.enth[6]  
        self.p.dt.iloc[7] = self.p.temp[7] - self.p.temp[6]
        self.p.phase.iloc[7] = PropsSI('Phase','T',self.p.temp[7],'P',self.p.pres[7],'CO2')
        
        
        # точка 2 - за рекуператором / горячая часть - охлаждение [CO2,H2O]
        
        self.p.dh.iloc[2] = - self.comb.gas.mass['CO2_recyc'] * self.p.dh[7]
        self.p.enth.iloc[2] = self.p.enth[1] + self.p.dh[2]
        self.p.temp.iloc[2] = self._tempH(pressure=self.p.pres[2], enthalpy=self.p.enth[2])
        
        self.p.entr.iloc[2] = PropsSI('S','T|supercritical',self.p.temp[2],'P',self.p.pres[2],self.fluid_mix)
        self.p.dens.iloc[2] = PropsSI('D','T|supercritical',self.p.temp[2],'P',self.p.pres[2],self.fluid_mix)
        self.p.sp_heat.iloc[2] = PropsSI('C','T|supercritical',self.p.temp[2],'P',self.p.pres[2],self.fluid_mix)
        
        self.p.dp_rel.iloc[2] = self.p.pres[1] / self.p.pres[2]
        self.p.dt.iloc[2] = self.p.temp[2] - self.p.temp[1]
        self.p.phase.iloc[2] = PropsSI('Phase','T|supercritical',self.p.temp[2],'P',self.p.pres[2],self.fluid_mix)
        
        
        ######################################################
        
        # точки 3, 4 - за сепаратором / осушение [CO2]
        
        self.p.temp.iloc[3:5] = self.p.temp[2]
        self.p.enth.iloc[3:5] = PropsSI('H','T',self.p.temp[3],'P',self.p.pres[3],'CO2') 
        self.p.entr.iloc[3:5] = PropsSI('S','T',self.p.temp[3],'P',self.p.pres[3],'CO2')
        self.p.dens.iloc[3:5] = PropsSI('D','T',self.p.temp[3],'P',self.p.pres[3],'CO2')
        self.p.sp_heat.iloc[3:5] = PropsSI('C','T',self.p.temp[3],'P',self.p.pres[3],'CO2')
        
        self.p.dt.iloc[3] = self.p.temp[3] - self.p.temp[2]
        self.p.dh.iloc[3] = self.p.enth[3] - self.p.enth[2] 
        
        self.p.dt.iloc[4] = self.p.temp[4] - self.p.temp[3]
        self.p.dt.iloc[5] = self.p.temp[5] - self.p.temp[4]
        self.p.dh.iloc[4] = self.p.enth[4] - self.p.enth[3]  
        self.p.phase.iloc[3:5] = PropsSI('Phase','T',self.p.temp[3],'P',self.p.pres[3],'CO2')
        
        
        # точка 0 - за охладителем / охлаждение [CO2,H2O]
        self.p.dh.iloc[5] = self.p.enth[5] - self.p.enth[4]  
        self.p.dh.iloc[0] = self.p.enth[0] - self.p.enth[7] # перепад энтальпий в охладителе
        
        self.p.dt.iloc[0] = self.p.temp[0] - self.p.temp[7]
        self.p.phase.iloc[0] = PropsSI('Phase','T|supercritical',self.p.temp[0],'P',self.p.pres[0],self.fluid_mix)
        
        # пинч-поинт
        self.g.at['pinch'] = self.p.temp[1] - self.p.temp[7]
        
        # полезная работа цикла
        self.g.at['work_cycle'] = abs(self.p.dh[1] - self.p.dh[6] * self.comb.gas.mass['CO2_recyc'])
        
        # эффективность    
        self.g.at['efc_cycle'] = self.g.work_cycle / self.p.dh[0]
        
        
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


# sco = Acycle()
# sco.cycle(pressure_min=8e6, pressure_rate=2.2, temperature=(310,1000, 760))

# print(sco.p)
# print(sco.g)
# print(f'{sco.g.efc_cycle = :.3f}')






# sco.ts_diagramm()

# # sco.optim(pressure=64., temperature=((32,500), (32, 500), (32, 600)), pressure_rate=(1.2,3.5))

# sco.power()
# sco.n_optim_power()

######################################
# print(self.p)
# print(self.g)


# print("*** program OK ***")
# sys.exit("*** STOP ***")
######################################








