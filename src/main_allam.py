# Calculation of cycle parameters
from allam import Acycle
sco = Acycle()
sco.cycle(pressure_min=8e6, pressure_rate=2.2, temperature=(310,1073,900))
print(sco.p)
print(sco.g)