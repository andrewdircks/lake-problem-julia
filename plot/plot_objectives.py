import numpy as np
import matplotlib.pyplot as plt
solution_f = "solutions/solution_shallow.txt"
solution_f2 = "solutions/solution_cayuga.txt"

nvars = 100
nobjs = 2

# for _f, label in [(solution_f, 'Shallow'), (solution_f2, 'Cayuga')]:
for _f, label in [(solution_f2, 'Cayuga')]:
    with open(_f, 'r') as f:
        obj1 = []
        obj2 = []
        for line in f:
            sol = line.split()
            obj1.append(-1*float(sol[nvars]))
            obj2.append(float(sol[nvars+1]))    
        plt.scatter(obj1, obj2, label=label, c='#ff7f0e')

plt.xlabel('Economic Benefit (maximize)')
plt.ylabel('Environmental Impact (minimize)')
plt.title('Lake Problem - Pareto Front')
plt.legend()
plt.show()

# 3  -> environmental
# 42 -> economic
# 30 -> compromise