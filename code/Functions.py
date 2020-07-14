"""
    Contains all functions and classes needed to run the notebooks and
    simulations for my Masters thesis on ontogenetic growth models done as part
    of the "Computational Methods in Ecology and Evolution" MSc at Imperial
    College London 2019/2020
"""


###### Imports ######
from numpy import arange, array
##optimisation
from numpy import unravel_index, argmax, isnan, nan_to_num, zeros
## maths functions
from scipy.integrate import odeint 
from numpy import exp, sin, pi, log10, log
## plotting
import matplotlib.pyplot as plt

###### Functions ######

## Hou et al 2011 ##

def hou_mt(m0, M, Em, B0):
    """
    to calculate m(t) to feed into equation 1.5 from hou et al 2011 / the `hou-dmdt` function.

    Args:
        m0 (int): initial mass of organism
        M (int): adult mass
        Em (int): metabolic energy required to synthesize one unit of biomass
        B0 (int): Taxon specific constant

    Returns:
        int: mt - the mass of an organism at time t
    """

    return ((1 - (1 - ((m0/M)**0.25)) * exp(-B0 / (4 * Em * (M**0.25))))**4) * M

def hou_AS(f, Ec, Em, B0, mt, M):
    """
    Calculates the assimilation (A) of an organism as described in Hou et al 2011.

    Args:
        f (int): activity scope
        Ec (int): energy stored in newly synthesised biomass
        Em (int): metabolic energy required to synthesize one unit of biomass
        B0 (int): Taxon specific constant
        mt (int): mass at time t
        M (int): adult mass

    Returns:
        int: Assimilation of organism of given size at time t
    """
    return (((f + (Ec / Em)) * B0 * (mt**0.75)) - ((Ec / Em) * B0 * (M**-0.25) * mt))

def hou_dmdt(m0, t, beta, f, Ec, Em, B0, B0FR,  M):
    """[summary]

    Args:
        m0 (int): Initial mass of organism
        beta (int): Proportion of feeding restriciton
        f (int): Activity scope
        Ec (int): Energy stored in newly synthesised biomass
        Em (int): Metabolic energy required to synthesize one unit of biomass
        B0 (int): Taxon specific constant ad libitum
        B0FR (int): Taxon specific constant under feeding restriction
        M (int): Adult mass

    Returns:
        int: Change in mass dm/dt
    """
    mt = hou_mt(m0, M, Em, B0)
    mFR = m0
    beta_AS = beta * hou_AS(f, Ec, Em, B0, mt, M)
    B_totFR = f * (B0*(mFR**0.75))

    dmdt = (beta_AS - B_totFR) / Ec

    return dmdt

def hou_integrate(m0, time, params):
    """
    A function to simulate the growth of an organism over time based on 
    Hou et al 2011.

    Args:
        m0 (int): Initial mass of organism
        time (int): The amouont of time to integrate over 
        params (dict): A dictionary of parameters:
            beta (int): Proportion of feeding restriciton
            f (int): Activity scope
            Ec (int): Energy stored in newly synthesised biomass
            Em (int): Metabolic energy required to synthesize one unit of 
                        biomass
            B0 (int): Taxon specific constant ad libitum
            B0FR (int): Taxon specific constant under feeding restriction
            M (int): Adult mass

    Returns:
        array: The growth of the organism over time
    """    
    
    t = arange(0, time, 1)
    beta = params["beta"]
    f = params["f"]
    Ec = params["Ec"]
    Em = params["Em"]
    B0 = params["B0"]
    B0FR = params["B0FR"]
    M = params["M"]

    return odeint(hou_dmdt, m0, t, args = (beta, f, Ec, Em, B0, B0FR,  M))


def plot_hou(m0, time, params):
    """
    A function to simulate the growth of an organism over time and plot the result based on Hou et al 2011.

    Args:
        m0 (int): Initial mass of organism
        time (int): The amouont of time to integrate over 
        params (dict): A dictionary of parameters:
            beta (int): Proportion of feeding restriciton
            f (int): Activity scope
            Ec (int): Energy stored in newly synthesised biomass
            Em (int): Metabolic energy required to synthesize one unit of 
                        biomass
            B0 (int): Taxon specific constant ad libitum
            B0FR (int): Taxon specific constant under feeding restriction
            M (int): Adult mass

    Returns:
        array: The growth of the organism over time
    """ 

    m = hou_integrate(m0, time, params)[:,0] #change dimensions from col to row
    t = arange(0, time, 1)

    plt.figure()
    plt.plot(t, m, label="Mass")
    plt.xlabel("Time")
    plt.ylabel("Mass")
    plt.legend()
    plt.show()

    return m


def hou_simple(m, t, B0, Em, M):
    """
    Function for calculating growth based on Eq 1.2 from Hou et al 2011.

    Args:
        m ([type]): [description]
        t ([type]): [description]
        B0 ([type]): [description]
        Em ([type]): [description]
        M ([type]): [description]

    Returns:
        [type]: [description]
    """
    Bm = B0 * M**-0.25

    dmdt = ((B0 * m**0.75) - (Bm * m)) / Em

    return dmdt

def hou_integrate_simple(m0, time, params):

    B0 = params["B0"]
    Em = params["Em"]
    M = params["M"]

    t = arange(0, time, 1)

    mass =  odeint(hou_simple, m0, t, args=(B0, Em, M))

    return mass

def plot_hou_simple(m0, time, params):


    m = hou_integrate_simple(m0, time, params)[:,0] #change dimensions from col to row
    t = arange(0, time, 1)

    plt.figure()
    plt.plot(t, m, label="Mass")
    plt.xlabel("Time")
    plt.ylabel("Mass")
    plt.legend()
    plt.show()
    return m



## Supply Model ##
def am(m, dimensionality = "3D"):
    """
    Calculates mass specific search rate in a functional response as derived in 
    Pawar et al 2012.

    Args:
        m (float): Mass of individual
        dimensionality (str, optional): Dimensionality of the functional 
        response. Defaults to "3D".

    Returns:
        float: Mass specific search rate (a)
    """
    if dimensionality == "3D":
        a0 = 10**-1.77  
        gamma = 1.05 # scarce resources
        return (m**gamma) * a0

    if dimensionality == "2D":
        a0 = 10**-3.08 
        gamma = 0.68 # scarce resources
        return (m**gamma) * a0

def hm(m, dimensionality = "3D"):
    """
    Calculates mass specific handling time in a functional response as derived
     in Pawar et al 2012 SI. 
    Scaling exponent is 0.75 due to caveats stated in the paper to do with the 
    gathering of data for saturated resources.

    Args:
        m (float): Mass of individual
        dimensionality (str, optional): Dimensionality of the functional response. Defaults to "3D".

    Returns:
        float: Mass specific handling time (h)
    """

    if dimensionality == "3D":
        tk0 = 10**3.04  
        beta = 0.75 
        return (m**-beta) * tk0

    if dimensionality == "2D":
        tk0 = 10**3.95 
        beta = 0.75 
        return (m**-beta) * tk0

def Xrt(t, amp, centre, period = 365):
    
    """
    To simulate the fluctuation of resource density in a functional response 
    through time according to a sine wave.

    Args:
        t (int): time passed (is converted to radians in function)
        amp (float): The amplitude of the sin wave as a percentage of centre. 
                    1 = 100% flucuation, i.e. from 0 to 2*centre
        centre (float): The value around which resource density fluctuates.
        period (int): Period of the wave in time. Defaults to 365


    Returns:
        float: Resource density
    """   

    x = t * (2 * pi / period) 

    return (amp*centre * sin(x)) + centre

def Fun_Resp(m, Xr, dimensionality = "3D"):
    """
    Calculates the functional response of an organism dependent on mass.

    Args:
        m (float): Mass of individual
        ## a0 ([type]): mass dependent search rate 
        R (float): Resource density
        ## h ([type]): [description]
        dimensionality (str): Used to determine how  serach rate and handling rate are calculated. See functions for details.

    Returns:
        [float]: consumption rate of the organism
    """    
    a = am(m, dimensionality)  # find mass dependent search rate
    h = hm(m, dimensionality) # find mass dependent handling time
    
    f = (a *Xr) / (1 + a*h*Xr)
    return f

def Bm (m, delta, proportion = 0.05, dimensionality = "2D"):
    """
    Calculated mass specific metabolic cost based on some percentage 
    of effective intake rate i.e. search rate

    Args:
        m (float): Mass of individual
        proportion (float): The proportion of intake rate that will be used
        dimensionality (str): Used to determine how  serach rate and handling rate are calculated. See functions for details.

    Returns:
        [float]: Mass specific metabolic cost
    """    

    if dimensionality == "3D":
        a0 = 10**-1.77  
        delta = 0.75 
        return proportion *  m**delta + a0

    if dimensionality == "2D":
        a0 = 10**-3.08 
        delta = 0.75
        return proportion *  m**delta + a0

# def Bm(m):
#     """
#     metabolic rate based on Gillooly et al 2001 "Effects of size and temperature on metabolic rate".
#     Uses Boltzmann's factor to approximate for the cost of all reactions in the body.

#     Args:
#         m (float): Mass of individual

#     Returns:
#         [float]: Mass specific metabolic cost

#     """
#     Ei = 0.6
#     k = 1.3807 * 10**-23 # boltzmann constant
#     T = 293.15 # 20C in Kelvin
#     M = m**-0.25 # so as it is m**0.75 when multiplied by mass
#     return exp(-Ei/(k*T))*M

def L(t, k = 0.01):
    """
    Suvival Function for reproduction modelled as an exponetially increasing number through time.
    Based on thinking that most organisms will not live long enough for factors such as 
    reproductive senescence to be a factor.

    Args:
        t (int): Time
        k (float, optional): Reproductive senescence. Defaults to 0.01.

    Returns:
        [type]: [description]
    """
    return exp(-k*t)

def reproduction(t, c, m, rho, alpha, k = 0.01):
    """
    Calculates the reproductive output of an organisms at time `t` in terms of biomass

    Args:
        t (int): Time
        c (float): Reproductive scaling constant
        m (float): Mass of individual
        rho (float): Reproductive scaling exponent
        k (float, optional): Reproductive senescence. Defaults to 0.01.]
        
    Returns:
        float: Reproductive output in terms of biomass

    """    
    Q = L(t-alpha) # mortality

    return Q * c * (m**rho)

def metabolic_cost(m):
    """
    Calculates the metabolic cost of an organism in term of mass/time from Barneche et al 2014.

    Args:
        m (float): Mass of individual (units: mass)
        metabolic_rate (float): The standard metabolic rate or resting metabolic rate of the organism (units: energy * mass / time)
        conversion_factor ([type]): Value for how much energy is in a unit of mass (units: energy / mass)

    Returns:
        float: The "mass cost" of the organism at the given mass
    """

    alpha = 0.76
    intercept = exp(-5.71)
    return intercept * m**alpha

def dmdt(mR0, t, 
         alpha, epsilon, norm_const, meta_prop, meta_exp, 
         c, rho, 
         Xr, amp, period, dimensionality = "3D"):
    """
    Calculates the instantaneous change in mass at time `t`. 

    Args:
        m (float): Mass of individual
        t (int): time
        alpha (int): maturation time
        epsilon (float): Efficiency term
        norm_const (float) : Normalisation constant, 
                            i.e. Functional response value for a 1kg organism
        meta_prop (float) : Proportion of optimum intake rate that 
                            is assigned to metabolism
        meta_exp (float) : Scaling exponent for metabolism
        c (float): Metabolic cost constant
        rho (float): Metabolic cost exponent
        Xr (float): The expected median value for resource density
        amp (float): amplitude of resource fluctuation around `Xr`
        period (int): the period (duration) of the resource cycle
        dimensionality (str): See `Func_Resp`


    Returns:
        float: change in mass (dm/dt) at time t
    """
    
    # check if individual is at/past maturation
    m, R = mR0
    k = 0.01 #  reproductive senesence
    if t < alpha:
        R = 0 # reproductive cost
    repro = 0 # reproductive output
    if t >= alpha:
        R = c * norm_const * (m**rho) # kg/d
        repro = repro_out = reproduction(t, c*norm_const, m, rho, alpha, k = 0.1)

    # Gain
    Xr_t = Xrt(t, amp, Xr, period)
    gain = epsilon * Fun_Resp(m, Xr, dimensionality) # kg/s
    gain = gain * (60 * 60 * 24)# /s -> /min -> /hour -> /day
    # Loss
    B_m  = norm_const * meta_prop * m**meta_exp  
    loss = B_m + R
    
    dmdt =  gain - loss
    
    # check for shrinking
    if dmdt + m < 0:
        dmdt = -m

        
    return array([dmdt, repro])
    
def dmdt_integrate(m0, R0, time, params):
    """
    integrates dmdt to return a growth curve.

    Args:
        m0 (float): Initial mass of individual
        R0 (float): Initial reproductive output
        time (int): Time
        params (dict): see `dmdt` function for needed names of parameters

    Returns:
        array: growth curve of organism, and reproductive output of organism
    """    

    # Organise Parameters for integration
    t = arange(0, time, 1)   
    mR0 = array([m0, R0])
    arg = (params["alpha"], params["epsilon"], 
            params["norm_const"], params["meta_prop"], params["meta_exp"], 
            params["c"], params["rho"], 
            params["Xr"], params["amp"], params["period"], 
            params["dimensionality"])

    # Simulate growth
    mR = odeint(func=dmdt, y0=mR0, t=t, args=arg)

    return mR

def plot_supply(m0, R0, time, params):
    """
    Plots the growth curve of a bottom up supply based model in the form:
    dmdt = [gain - loss]m ;  where gain is modeled as a functional response 
    scaled by an efficiency constant and loss is metabolic and reproductive 
    costs.

    Args:
        m0 (float): Initial mass of individual
        R0 (float): Initial reproductive output
        time (int): Time
        params (dict): see `dmdt` function for needed names of parameters

    Returns:
        array: growth curve of organism, and reproductive output of organism
    """

    # Simulate growth
    mR = dmdt_integrate(m0, R0, time, params) 

    # unpack results
    m = mR[:,0]
    repro = mR[:,1]

    t = arange(0, time, 1)

    plt.figure()
    plt.plot(t, m, label="Mass") #change dimensions from col to row
    plt.plot(t, repro, label="Reproductive Output") 
    plt.xlabel("Time")
    plt.ylabel("Mass")
    plt.legend()
    plt.show()
    
    return mR
       

## Optimisation


def find_max(arr):
    """
    A function to find the maximum of an array and return its indices

    Arguments:
        arr {np.array} -- the matrix to find the max value of 

    Returns:
        {tuple} --  the indices of the maximum value

    """

    max_ind = unravel_index(argmax(array, axis=None), arr.shape)
    return max_ind

def find_optimum(c_vec, rho_vec, m0, R0, time, params):
# do the c and rho vectors as a meshgrid 
# this should speed it up 
# from there the same as the notebook with all the debugging left to do

    
    # array to store final reproduction values
    # `c` will be row and `rho` columns, ith val is ith val in c_vec or jth val is jth in rho_vec
    repro_array = zeros((len(c_vec), len(rho_vec)))

    for i, c in enumerate(c_vec):
        params["c"] = c
        for j, rho in enumerate(rho_vec):
            params["rho"] = rho
            result = dmdt_integrate(m0, R0, time, params)
            mass = result[:, 0]
            repro = result[:, 1]
            total_repro = repro[-1]
            if mass[-1] < mass[params["alpha"]]:
                total_repro = -1 # if the fish shrinks
            repro_array[i, j] = total_repro
    repro_array = nan_to_num(repro_array) # replace `nan` with 0
    max_ind = find_max(repro_array)
    max_repro = repro_array[max_ind]
    c_opt = c_vec[max_ind[0]]
    rho_opt = rho_vec[max_ind[1]]

    return array([c_opt, rho_opt])
###### Classes ######



###### Testing ######
# a section to test functions and other functionality
# trying to see what is causing m to drop to 0 in most cases
# from numpy import linspace

###### Notes / To Do ######
