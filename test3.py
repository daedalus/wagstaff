base = 2
power = 5

#if power == 5:
tmp = 0
for i in range(0, power):
  tmp += ((-1) ** i) *  base ** i
  
print(base + 1, base ** 4 - base ** 3 + base ** 2 - base + 1)
print(tmp)

#return factor_special_forms(p) + factor_special_forms(q)

