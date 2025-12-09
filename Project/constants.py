import numpy as np

# equilibrium conditions for basins #
equator_equilibrium_sal = 35.2 #[psu]

def del_rho(T1,T2,S1,S2):
    '''
    Calculate density difference using T and S difference

    delT(float): Tupper-Tlower
    delS(float): Supper-Slower
    '''

    # physical constants #
    alpha = 0.2 #[kg m-3 degC-1]
    beta = 1 #[kg m-3 psu-1]

    return -alpha*(T1-T2) + beta*(S1-S2)

def density_water(T, S, p=1.01325):
    # Density of pure water as a function of the temperature (kg/m3)
    # water salinity = 0 ppt
    rhow = 999.842594 + (6.793952E-2 * T) - (9.095290E-3 * T**2) + \
            (1.001685E-4 * T**3) - (1.120083E-6 * T**4) + (6.536332E-9 * T**5)
    
    # Water density (kg/m3) at one standard atmosphere, p = 0. 
    rhost0 = rhow + (S * (0.824493 - (4.0899E-3 * T) + (7.6438E-5 * T**2) - \
            (8.2467E-7 * T**3) + (5.3875E-9 * T**4))) + ( S**(3/2) * \
            (-5.72466E-3 + (1.0227E-4 * T) - (1.6546E-6 * T**2))) + \
            ( 4.8314E-4 * S**2)
    
    # Water pure secant bulk modulus
    kw = 19652.21 + (148.4206 * T) - (2.327105 * T**2) + \
        (1.360477E-2 * T**3) - (5.155288E-5 * T**4)
    
    # Water secant bulk modulus at one standard atmosphere (p = 0)
    kst0 = kw + (S * (54.6746 - (0.603459 * T) + (1.09987E-2 * T**2) - \
        (6.1670E-5 * T**3))) + (S**(3/2) * (7.944E-2 + (1.6483E-2 * T) - \
        (5.3009E-4 * T**2)))
    
    # Water secant bulk modulus at pressure values, p
    p=p/10   # use decibar (~1m depth) instead of bar (1 atm)

    kstp = kst0 + \
        (p*(3.239908 + (1.43713E-3 * T) + (1.16092E-4 * T**2) -(5.77905E-7* T**3))) + \
        ((p*S) *(2.2838E-3 - (1.0981E-5 * T) - (1.6078E-6 * T**2))) + \
        (1.91075E-4 * p * S**(3/2)) + \
        (p**2 * (8.50935E-5 - (6.12293E-6 * T) + (5.2787E-8 * T**2))) + \
        ((p**2*S) * (-9.9348E-7 + (2.0816E-8 * T) + (9.1697E-10 * T**2)) )
    
    # Water density at any pressure (kg/m3)
    rho = rhost0 / (1-(p/kstp))
    return rho

# strength of vertical mixing #
def vertical_gamma(T_u,T_l,T_e,S_u,S_l,S_e):
    '''
    Returns a mixing time scale (units [time-1]).

    del_rho(float): density_upper - density_lower
    '''

    drho = del_rho(T_u, T_l, S_u, S_l)

    rho_b = density_water(T_l, S_l, p=101)
    
    m = ((400-10) / (0.1*rho_b))
    tau = (m * drho) + (m * 0.05 * rho_b)
    
    gamma = 1/(2*tau)

    return gamma
    
def calc_q(T_u,T_l,T_e,S_u,S_l,S_e):
    #takes in the upper and lower and returns the q value calcuated differently for upper and lower 
    
    k = 1.5e-6 #rate constant from overturning circulation reference 
    rho_o = 1e3 #average density 
    alpha = 1.5e-4 #1/deg C
    beta = 8e-4 #1/psu

    q_u = k*(alpha*(T_u - T_e) - beta*(S_u - S_e)) 
    q_l = k*(alpha*(T_l - T_e) - beta*(S_l - S_e)) 
    
    print(q_u, q_l)

    return [q_u, q_l]

# ice formation #
def icegrowth_per_day(T_u, Tatm, S_u, h_init, N_timestep, sealevelrise_rate):
    '''
    Determines the growth rate of the sea ice
    given air and water temperature.
    
    T_u(float): temperature of upper polar box
    Tatm(float): temperature of air above upper polar box
    hinit(float): initial ice thickness
    '''

    #constants:

    k_ice = 2.2 #thermal conductivity of ice in W/m·K
    latentheat = 334000 #latent heat of fusion of ice in J/kg
    rho_ice = 916 #density of ice in kg/m^3 at 0 degrees C
    ice_S = 5 #salinity of ice in psu
    Tmelt = 0 # could be -1.8C
    D = 100 # ice box depth, arbitrary 
    #note: thermal conductivity actually depends on temperature; latent heat on salinity

    D_current = D + N_timestep * sealevelrise_rate

    h_total = h_init
    const = 2*k_ice/(latentheat*rho_ice)

    if Tatm < Tmelt:
        G = np.sqrt(const * abs(Tmelt - Tatm - T_u) * 60 * 60 * 24)  #G is growth in ice thickness per day
    elif Tatm >= Tmelt:
        G = -np.sqrt(const * abs(Tmelt - Tatm) * 60 * 60 * 24) # G is melt in ice thickness per day

    h_total_new = h_total + G # ice thickness after growth or melt

    h_total_new = min(20., h_total_new) # ice stops growing after 10 m
    
    h_total_new = max(0, h_total_new) # ice thickness not negative
    
    G_effective = h_total_new - h_total # diagnose how much ice grew/melted, given the chosen limits

    waterdepth = D_current - h_total_new # new depth of unfrozen water

    S_new = (waterdepth * S_u + G_effective * (S_u - ice_S)) / waterdepth
    
    d_Sw = S_new - S_u
    
    return d_Sw, h_total_new

# ice formation #
def sealevelrise(Tatm, S_u,  N_timestep, sealevelrise_rate):

    #constants:
    fresh_S = 0.5

    D = 100
    D_current = D + N_timestep * sealevelrise_rate

    S_new = ((D_current * S_u) + (sealevelrise_rate * fresh_S)) / (D_current + sealevelrise_rate)

    # if Tatm <= 0:
    #     d_Sw = 0
    # else:
    d_Sw = S_new - S_u

    return d_Sw
