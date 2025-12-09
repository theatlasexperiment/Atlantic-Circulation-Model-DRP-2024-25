import numpy as np
from time_step import dState

State_0 = np.array([10., 10., 20., 35., 34., 34., 0.]) #[T_u, T_l, T_e, S_u, S_l, S_e, h_ice]

current_state = State_0
dt = 1 #day

for day in range(1, 1000):
    current_state[:-1] += dState(current_state, day)[:-1] * dt #simplest possible Euler integration
    current_state[-1] = dState(current_state, day)[-1]
    print(current_state)
