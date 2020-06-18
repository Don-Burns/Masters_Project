"""
    Contains all functions and classes needed to run the notebooks and
    simulations for my Masters thesis on ontogenetic growth models done as part
    of the "Computational Methods in Ecology and Evolution" MSc at Imperial
    College London 2019/2020
"""


###### Imports ######
from scipy import arange as arange

## maths functions
from scipy.integrate import odeint as odeint
from scipy import exp as exp
from scipy import sin as sin
from scipy import pi as pi

## plotting
import matplotlib.pyplot as plt

###### Functions ######

## Hou et al 2011 ##

def hou_mt(m0, M, Em, B0):
    """
    to calculate m(t) to feed into equation 1.5 from hou et al 2011 / the ___ function.

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
            Em (int): Metabolic energy required to synthesize one unit of biomass
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
            Em (int): Metabolic energy required to synthesize one unit of biomass
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
def Fun_Resp(m, a0, R, h):
    """[summary]

    Args:
        m ([type]): mass of indicidual
        a0 ([type]): mass dependent search rate 
        R ([type]): [description]
        h ([type]): [description]

    Returns:
        [type]: [description]
    """    
    f = ((a0 * m**0.75)*R) / (1 + (a0 * m**0.75)*h*R)
    return f

def dmdt(m, t, epsilon, L_B, L_R, a0, R, h, amp, period):
    """[summary]

    Args:
        m (float): [description]
        t (int): [description]
        epsilon (float): [description]
        L_B (float): [description]
        L_R (float): [description]
        a (float): [description]
        R (float): The expected median value for resource density
        amp (float): [description]
        period (int): [description]

    Returns:
        float: change in mass (dm/dt) at time t
    """
    # put params as a dict?

    R_t = Rt(t, amp, period, centre=R)

    gain = epsilon * Fun_Resp(m, a0, R_t, h)
    loss = (L_B  ) + (L_R)
    dmdt = ((gain) - loss) * m # `gain` is times m**0.75
    
    return dmdt

def dmdt_integrate(m0, time, params):
    """
    integrates dmdt to return a growth curve.

    Args:
        m0 (float): Initial mass of individual
        time (int): Time
        params (dict): see `dmdt` function for needed names of parameters

    Returns:
        [type]: [description]
    """    
    t = arange(0, time, 1)   
    arg = (params["epsilon"], params["L_B"], params["L_R"], params["a0"], params["R"], params["h"], params["amp"], params["period"])

    mass = odeint(dmdt, m0, t, args=arg)

    return mass

def plot_supply(m0, time, params):
    """
    Plots the growth curve of a bottom up supply based model in the form:
    dmdt = [gain - loss]m ;  where gain is modeled as a functional response 
    scaled by an efficiency constant and loss is metabolic and reproductive 
    costs.

    Args:
        m0 (float): Initial mass of individual
        time (int): Time
        params (dict): see `dmdt` function for needed names of parameters
    """

    m = dmdt_integrate(m0, time, params)[:,0] #change dimensions from col to row
    t = arange(0, time, 1)

    plt.figure()
    plt.plot(t, m, label="Mass")
    plt.xlabel("Time")
    plt.ylabel("Mass")
    plt.legend(loc="upper left")
    plt.show()
    
    return m

def Rt(t, amp, period, centre):
    """
    To simulate the fluctuation of resource density in a functional response 
    through time according to a sine wave.

    Args:
        t (int): time passed (is converted to radians in function)
        amp (float): The amplitude of the sin wave
        period (int): Period of the wave in time
        centre (float): The value around which resource density fluctuates

    Returns:
        float: Resource density
    """   

    x = t * (2 * pi / period) 

    return amp * sin(x) + centre



###### Classes ######



###### Testing ######
# a section to test functions and other functionality


###### Notes / To Do ######
# may want to define B0 and B0FR based on b0 and calculate in function?
