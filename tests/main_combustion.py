from allam import Combust
comb = Combust()
temp = comb._burnTemp(k_recyc=20, temp_recyc=900)
print(f'Calorimetric combustion temperature {temp:.1f} °C')