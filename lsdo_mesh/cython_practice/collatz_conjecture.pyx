# if even: divide by 2
# if odd: multiply by 3 and add 1
import numpy as np
cimport numpy

def collatz_steps(val):

    if val == 1:
        raise ValueError('Input value must be greater than 1.')
    elif not isinstance(val, int):
        raise TypeError('Input value must be a whole number.')

    odd_steps, even_steps = 0, 0
    active_val = val
    path = []
    counter = 0
    while True:

        if np.linalg.norm(active_val - 1.0) < 1e-8: # FINAL STEP
            print(np.linalg.norm(active_val - 1.0))
            total_steps = even_steps + odd_steps
            break

        if active_val % 2 == 0.: # EVEN NUMBER
            active_val /= 2
            even_steps +=  1 
            
        elif active_val % 2 == 1.: # ODD NUMBER
            active_val *= 3
            active_val += 1
            odd_steps +=  1

        
        
        counter += 1
        print(active_val)
        print('Working on step {}'.format(counter))
        path.append(active_val)
        
    return path, even_steps, odd_steps, active_val