import numpy as np
from HH_brian2 import *
from scipy.io import savemat


# adjacancy matrix
A = np.array([
            [0, 0, 1, 0],
            [0, 0, 0, 0],
            [0]*4,
            [0]*4
            ])

#input stimulus defined as [neuron index, start time, duration (ms), amplitude]
stimulus = np.array([
                [0, 2, .5, .04],
                [1, 2, .5, .04],
                [0, 10, .5, .15],
                [1, 25, .5, .15],
                [0, 35, .5, .15],
                [1, 45, .5, .15]
            ])

#simulation duration in ms
duration = 80


[t, v, I_ext] = simulate_HH_network(A, stimulus, duration, dt=.01, post_noise_sigma=0.05)

out = v/volt

dic = {"out": out}

savemat("testmat.mat",dic)