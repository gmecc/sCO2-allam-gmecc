from allam import PropSCO2
temperature = (32, 40, 50, 60, 600)
ts = PropSCO2()
ts.diagramm_pressure(pressure=(1e6, 20e6), temperature=temperature)