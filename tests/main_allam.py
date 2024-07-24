# Calculation of cycle parameters
from allam import Acycle
sco = Acycle()
sco.cycle(pressure_min=8e6, pressure_rate=2.2, temperature=(310,1073,900))
print(sco.p)
print(sco.g)

sco.efc_temp_recyc(var_temp_recyc=(500,900), pressure_min=8e6, pressure_rate=2.2, temperature=(310,1073))
# sco.efc_temp_heat(var_temp_heat=(700,1000), pressure_min=8e6, pressure_rate=2.2, temperature=(310,600))