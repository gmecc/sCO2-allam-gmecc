from allam import Combust
comb = Combust()
temp = comb._burnTemp(k_recyc=20, temp_recyc=900)
print(f'Калориметрическая температура горения {temp:.1f} С')