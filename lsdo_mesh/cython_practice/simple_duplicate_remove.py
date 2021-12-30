from simple_duplicate_remove import simple_duplicate_remove
import numpy as np

max_val = 10
num_val = 178760
full_list = []
for i in range(num_val):
    full_list.append(np.random.randint(1, max_val))

full_list = np.array(full_list, dtype=int)

test = simple_duplicate_remove(full_list)
# print(test[0])
# print(test[1])
# print(test[2])