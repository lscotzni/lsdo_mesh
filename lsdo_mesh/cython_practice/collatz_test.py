from collatz_conjecture import collatz_steps
import random
import numpy as np
import matplotlib.pyplot as plt


rand_num = 10
rand_int = []

plt.figure(1, figsize=(10,10))

for i in range(rand_num):
    print('Working on random integer {}'.format(i+1))
    rand_int.append(random.randrange(2,20))
    path = collatz_steps(rand_int[i])
    plt.plot(path[0])
    plt.xlabel('Iterations')
    plt.ylabel('Function Value')

plt.show()
