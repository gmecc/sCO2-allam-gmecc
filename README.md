# Allam cycle calculation

## Allama cycle diagram

![Allama cycle diagram](image/alam-sceme.jpg)

# Calculation of fuel combustion processes in the combustion chamber of a gas turbine engine
Program for calculating the fuel combustion process in the combustion chambers of gas turbine plants

The program is based on the liquid properties library `coolProp`.

### Units
- temperature - degree Celsius _[°C]_;
- absolute pressure _[Pa]_

## Variables
- `fluid` - name of gas (class instance `FluidList`)
- `mol_mass` - molar mass _[kg/mol]_
- `mol` - mole fraction
- `mass` - mass fraction

## Initial data
The fractional composition of the fuel is specified 
in the form of a dictionary. Key - faction name (str). Properties tuple:
- `fluid`
- `mol`
- volumetric heat of combustion at 0 deg. C _[J/m3]_






## About the author
Sergey Besedin, dr. of sc., prof.

Andry Fydorov, engineer