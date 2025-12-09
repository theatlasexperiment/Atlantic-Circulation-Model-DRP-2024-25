from constants import vertical_gamma, calc_q, equator_equilibrium_sal, icegrowth_per_day, sealevelrise
import pandas as pd
import numpy as np


def read_eq_atm(timestep, filedir='/Users/achereque/Documents/DRP/Project/'):
    forcing = pd.read_csv(filedir + 'simulation.ice_age.csv')

    # if you don't have pandas
    #forcing = np.genfromtxt(filedir + 'simulation.getting_colder.csv', delimiter=',', skip_header=1)
    #return forcing[(timestep-1)%365, 2]

    return forcing.loc[(timestep-1)%365].T_atm_eq

def read_u_atm(timestep, filedir='/Users/achereque/Documents/DRP/Project/'):
    forcing = pd.read_csv(filedir + 'simulation.ice_age.csv')
    
    # if you don't have pandas
    #forcing = np.genfromtxt(filedir + 'simulation.getting_colder.csv', delimiter=',', skip_header=1)
    #return forcing[(timestep-1)%365, 1]

    return forcing.loc[(timestep-1)%365].T_atm_pole

def dState(current_state,timestep):

    T_u,T_l,T_e,S_u,S_l,S_e,h_ice = current_state

    # temperature equilibration with atmosphere
    a = 1 / 365 #[1 year timescale]

    # salinity equilibration
    c = 1 / (10*365) #[10 year timescale]

    # sea level rise rate
    sealevelrise_rate = 0.0001

    # overturning circulation
    qu, qb = calc_q(T_u, T_l, T_e, S_u, S_l, S_e)

    # vertical mixing strength
    gamma = vertical_gamma(T_u,T_l,T_e,S_u,S_l, S_e)

    # boundary conditions
    Teq_atm = read_eq_atm(timestep)
    Tu_atm = read_u_atm(timestep)

    # ice 
    d_Sw, h_ice = icegrowth_per_day(T_u, Tu_atm, S_u, h_ice, timestep, sealevelrise_rate)

    # sea level rise dS
    d_S_slr = sealevelrise(Tu_atm, S_u, timestep, sealevelrise_rate)

    dTeq = a*(Teq_atm - T_e) - qu * (T_u - T_e) + qb * (T_l - T_e)
    dTu = a * max(0, 1-h_ice/5) * (Tu_atm - T_u) + qu * (T_u - T_e) + gamma * (T_l - T_u)
    dTb = -qb * (T_l - T_e) - gamma * (T_l - T_u)

    dSeq = c*(equator_equilibrium_sal - S_e) - qu * (S_u - S_e) + qb * (S_l - S_e)
    dSu = qu * (S_u - S_e) + gamma * (S_l - S_u) + d_Sw + d_S_slr
    dSl = -qb * (S_l - S_e) - gamma * (S_l - S_u)

    return np.array([dTu, dTb, dTeq, dSu, dSl, dSeq, h_ice])