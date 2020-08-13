"""
Author: Donal Burns
Data: 09/08/2020
This file is to generate the points for the sensitivity analysis.
Since output is messy progress report are save to report_plots.log
"""


## Imports
import Functions as F
import csv
from numpy import arange, around, array, linspace, isnan

from copy import deepcopy as deepcopy
# plotting
import matplotlib.pyplot as plt
import seaborn as sb
    # for in picture in picture graphs
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from time import time as timer
import multiprocessing as mp
import logging




## plot sizes
a4_sheet = array([8.3, 11.7]) # in inches
a4_width = a4_sheet[0]
a4_heigth = a4_sheet[1]
a4_box = array([8.3, 8.3])
titlefont = 10
labelfont = 5




###functions
def sensitivity(key, values, c_vec, rho_vec, m0, R0, time, params, shrinkage):
    """
        key (str) : the name of the param being examined in the parameter dictionary
        values (list) : the values for the key 
        params (dict) : the parameter dictionary
    """
    start = timer()
    c_list = []
    rho_list = []
    dim = params["dimensionality"]
    for i, val in enumerate(values):
        if key == "shrinkage": # in case i am looking at shrinking
            shrinkage = val
        else:
            params[key] = val

        c, rho = F.find_optimum(c_vec, rho_vec, m0, R0, time, params, shrinkage)

        # handle results with reproduction which have all values as optima so have 0 as optimum
        if c ==0  or rho == 0:
            c = None
            rho = None 


        c_list.append(c)
        rho_list.append(rho)

        if i % 3 == 0:
            progress_string = str(round(i/len(values)*100)) + "% finished " + key + dim
            logging.info(progress_string)
            logging.info(f"running {dim} {key} for {timer() -  start}")
    logging.info(f"{dim} {key} took: {timer() - start}")
    return array([c_list, rho_list])


def sens_and_save_res (dim, key, values, c_vec, rho_vec, m0, R0, time, params, shrinkage, prefix):
    """
    Finds the optimal rho and c values which lead to maximum reproductive output. Saves the results and plots both a line and scatter plot sensitivity analysis style plot.

    Args:
        dim (str): "2D" or "3D"
        key (str): The value in param dict being changed
        values (list): values for the key to be tested
        c_vec ([list): vector of possible c values
        rho_vec (list): vector of possible rho values
        m0 (float): starting mass
        R0 (float): starting reproductive output
        time (int): length of time to simulate over in days
        params (dict): parameter dictionary for dmdt see Functions.py `dmdt`
        shrinkage (float): amount of shrinking allowed at maturation
        prefix (str): label naming files
    Returns:
        [type]: [description]
    """

    temp_params = deepcopy(params)
    temp_params["dimensionality"] = dim
    c_list, rho_list = sensitivity(key, values, c_vec, rho_vec, m0, R0, 
                time, temp_params, shrinkage)
    
    with open(f"../results/sensitivity{prefix}{key}{dim}.csv", "w") as file:
        writer = csv.writer(file)
        writer.writerow(["c", "rho", key, "Xr", "metabolic exponent"])
        for i in range(len(c_list)):
            writer.writerow([c_list[i], rho_list[i], values[i], temp_params["Xr"], temp_params["meta_exp"]])
    
    return 0


### MAIN ####
if __name__ == "__main__":

    logging.basicConfig(filename="report_plots.log", level=logging.INFO)

        #### Parameters #####
    m0 = 1**-3 #1mg#10**-4 #1g
    R0 = 0
    time = 365*10 #10**6
    shrinkage = 0.05
    resolution = 0.01
    c_vec = around(arange(0, 0.4, resolution), decimals=4)
    rho_vec = around(arange(0, 2, resolution), decimals=4)

    params = {"alpha" : 365, "epsilon" : 0.7,
          "norm_const" : None, "meta_prop" : None, "meta_exp" : 1,
          "c" : None, "rho" : None, 
          "Xr" : 0.01, "amp" : 0, "period" : 365, "dimensionality" : None
}

    
    prefix = "LowResources_2D_metaexp1" # for naming files
    dimensions = ["2D"]#["2D", "3D"]
    keys = ["alpha", "shrinkage"]#["Xr", "meta_exp", "alpha", "shrinkage"]





    resources = linspace(0.01, 0.5, 25)#linspace(0.01, 0.5, 50)#linspace(0.5, 30, 50)#
    exponents = arange(0.75, 1.001, 0.01)#arange(0.75, 1.001, 0.01)
    alphas = around(linspace(1, 20, 20, dtype=int), decimals=0)#around(linspace(1, 2000, 100, dtype=int), decimals=0)
    shrink = linspace(0, 0.3, 20) # shrink percentages
    processes = []
    for key in keys:
        if key == "Xr":
            values = resources
        elif key == "meta_exp":
            values = exponents
        elif key == "alpha":
            values = alphas
        elif key == "shrinkage":
            values = shrink
        else:
            logging.info("key incorrect.  Key = ", key)
            ValueError
        for dim in dimensions:
            logging.info(f"{key} {dim} being assigned process")
            p = mp.Process(target=sens_and_save_res, args=(dim, key, values, c_vec, rho_vec, m0, R0, time, params, shrinkage, prefix))
            p.start()
            processes.append(p)

    for process in processes:
        logging.info("Joining Processes")
        process.join()
        process.close()


