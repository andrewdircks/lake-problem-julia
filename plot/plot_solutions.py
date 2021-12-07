import random
import matplotlib.pyplot as plt
from lake import Lake

# solution_f = "solutions/solution_shallow.txt"
solution_f = "solutions/solution_cayuga.txt"

threshold = 0.54454

def plot_soln(p_release, label=None, color=None):
    years = list(range(101))
    lake_state = Lake(random.randint(0,1000000), p_release)
    plt.plot(years, lake_state, label=label, c=color)


nvars = 100
nobjs = 2

labels = {
    3: 'Environmental',
    30: 'Compromise',
    42: 'Ecomonic'
}
colors = {
    3: 'green',
    30: 'blue',
    42: 'red'
}


with open(solution_f, 'r') as f:
    for i, line in enumerate(f):
        sol = line.split()
        economic = -1*float(sol[100])
        # if i in {3, 42, 30}:
        if True:
            plot_soln([float(x) for x in sol[:100]])


# plt.plot(list(range(101)), [threshold]*101, label='Lake Eutrophication Threshold', c='#ff7f0e')
plt.xlabel('Year', size=14)
plt.ylabel('Lake Phosphorus Concentration', size=13)
# plt.legend()
plt.title('Lake Cayuga - Solutions', size=13)
plt.show()
