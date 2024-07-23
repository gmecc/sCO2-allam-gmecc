from allam import Combust
comb = Combust()

# Calorimetric combustion temperature and gas composition calculation
temp = comb._burnTemp(k_recyc=20, temp_recyc=900)
print(f'Calorimetric combustion temperature {temp:.1f} °C')

# Calculation of recirculation coefficient
k_recyc = comb.burnAlpha(temp_gas=1000, temp_recyc=580)
print(f'Recirculation coefficient {k_recyc:.2f}')

print(comb.gas)
print(comb.gas_in)

comb.tempAlphaPl(k_recyc=(5, 25), temp_recyc=580)