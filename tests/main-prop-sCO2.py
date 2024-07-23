from allam import PropSCO2
temperature = (32, 40, 50, 60, 600)
ts = PropSCO2()

# Dependence of specific isobaric heat capacity $sCO_2$ on pressure
# at different temperatures
ts.diagramm_pressure(pressure=(1e6, 20e6), temperature=temperature)

# Dependence of specific isobaric heat capacity on pressure at different temperatures
ts.diagramm_temperature(pressure=(7.4e6, 10e6, 15e6, 20e6), temperature=(20, 200))

# Compressibility
ts.compressibility(pressure=(5e6, 10e6), temperature=(32, 40, 60, 100, 300, 500))