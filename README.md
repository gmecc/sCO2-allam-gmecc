# Allam cycle calculation
The Allam cycle, patented in 2013 by Rodney Allam, 
uses oxy-combustion and a supercritical CO2 stream 
as the working fluid to produce high-purity liquid pipeline CO2.

The program is based on the liquid properties library `coolProp`.

## Navigation



## Project structure



## Allam cycle diagram

![Allam cycle diagram](image/allam-scheme.jpg)

### Units systems
- temperature - Kelvins _[K]_;
- absolute pressure _[Pa]_

## List of variables
- `fluid` - name of gas (class instance `FluidList`)
- `mol_mass` - molar mass _[kg/mol]_
- `mol` - mole fraction
- `mass` - mass fraction

## List of methods

## Examples


## Parameters
The fractional composition of the fuel is specified 
in the form of a dictionary. Key - faction name (str). Properties tuple:
- `fluid`
- `mol`
- volumetric heat of combustion at 0 deg. C _[J/m3]_

## Returns:



## Термодинамические свойства sCO2. Теоретическая полезная работа цикла


## Тепловая модель камеры сгорания


## Расчет цикла


### Phase diagram CO2
```python
from allam import PTdiagrCO2
diagr = PTdiagrCO2()
diagr.plot()
```

## Исследование рекуператора


## About the authors
Sergey Besedin, dr. of sc., prof.

Andry Fydorov, engineer