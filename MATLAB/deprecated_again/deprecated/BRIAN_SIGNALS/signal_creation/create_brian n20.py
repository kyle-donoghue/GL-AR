import numpy as np
from HH_brian2 import *
from scipy.io import savemat
import matlab.engine

# adjacancy matrix
eng = matlab.engine.start_matlab()
paths = eng.genpath("../graph_creation")
eng.addpath(paths)
A = eng.construct_graph(20,'er',0.2)
A = eng.to_directed(A)
A = np.array(A)
print(A)

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


[t, v, I_ext] = simulate_HH_network(A, stimulus, duration, dt=.01, post_noise_sigma=0.00)

out = v/volt

dic = {"out": out, "A": A}

savemat("n20.mat",dic)