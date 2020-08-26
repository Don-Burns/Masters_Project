<img src = "./writing/images/IC_crest.svg" width = 100 height = 100 />


# Masters Thesis
My masters thesis for the 2019/2020 academic year as part of the Computational Methods in Ecology and Evolution at Imperial College London.

The project looks at ontogenetic models (OGMs) and using them to explore hyper-allomentric scaling in fish reproduction and how this may possibly relate to previously found hyper-allomentric scaling in consumption rate of consumers in 3D space.

## Getting started

Once the repository has been cloned locally, all code to produce the results found in the project can be found in the [`code`](https://github.com/Don-Burns/Masters_Project/tree/master/code) directory. All results will be output to the [`results`](https://github.com/Don-Burns/Masters_Project/tree/master/results) directory.  `data` contains the results which were obtained at the time of writting the thesis.  As such for a fresh `clone` the file will need to be copied to [`results`](https://github.com/Don-Burns/Masters_Project/tree/master/results)  if the user does not wish to generate them locally.

# [`code`](https://github.com/Don-Burns/Masters_Project/tree/master/code)
The project began looking at OGMs following the framework of [West et al. 2001](https://www.researchgate.net/publication/11676588_A_general_model_for_ontogenic_growth) and later changed to a resource supply based approach.  As such the notebooks within [`code`](https://github.com/Don-Burns/Masters_Project/tree/master/code) can be divided into these two periods.
[`Functions.py`](https://github.com/Don-Burns/Masters_Project/blob/master/code/Functions.py) contains functions which are used across multiple of the notebooks.

## OGMs
### Notebooks
[`Hou_Exploration.ipynb`](https://github.com/Don-Burns/Masters_Project/blob/master/code/Hou_Exploration.ipynb) :
An exploration of trying to include restricted resource supply following the work of Hou et al. [2008](doi.org/10.1126/science.1162302) and [2011](https://royalsocietypublishing.org/doi/10.1098/rspb.2011.0047).  

[`Model_Functions.ipynb`](https://github.com/Don-Burns/Masters_Project/blob/master/code/Model_Functions.ipynb) : A script using West et al.'s and exploring the effects of each parameter on the shape of the growth curve.

[`Model_Implementation_Exploration.ipynb`](https://github.com/Don-Burns/Masters_Project/blob/master/code/Model_Implementation_Exploration.ipynb) :  The initial exploration of using West et al.. Breaking down the various parameters and their derivations.  With some initial function development for integration over the equation.

[`Optimisation.ipynb`](https://github.com/Don-Burns/Masters_Project/blob/master/code/Optimisation.ipynb) : An initial attempt at finding values which yield optimal reproduction within the model.  Has an attempt at using optimisation algorithms before implementing the methods used in 2019 by Luke Vassor.


[`Other_Angles.ipynb`](https://github.com/Don-Burns/Masters_Project/blob/master/code/Other_Angles.ipynb) : Breaks down the methodology of Hou et al. and what the parameters in both the 2011 feeding restricted and 2008 ad libitum model mean.

## Resource Supply 
### Notebooks
[`Functional_Response_Model.ipynb`] : An initial attempt at building the supply model which incorrectly used the search and handling scaling from [Pawar et al. 2012](https://www.researchgate.net/publication/227857329_Dimensionality_of_consumer_seach_space_drives_tropic_interaction_strengths).

[`Optimise2.ipynb`](https://github.com/Don-Burns/Masters_Project/blob/master/code/Optimise2.ipynb) : Optimisation of the supply model and viewing the effects of shifting parameter values on the produced heatmaps and growth curves.

[`Report_Plots.ipynb`](https://github.com/Don-Burns/Masters_Project/blob/master/code/Report_Plots.ipynb) : Contains the code to produce the heatmaps used in the thesis and some prototyping for `report_plots.py`.

[`Search_and_Handling_Check.ipynb`](https://github.com/Don-Burns/Masters_Project/blob/master/code/Search_and_Handling_Check.ipynb) : A simple check of whether the seach rate and handling time functions match what is seen in Pawar et al. 2012.


[`Sensitivity_Plots.ipynb`](https://github.com/Don-Burns/Masters_Project/blob/master/code/Sensitivity_Plots.ipynb) : Contains the sensitivity analysis plots for the thesis and other miscellaneous plots such as the scaling plot.


[`Supply_Model.ipynb`](https://github.com/Don-Burns/Masters_Project/blob/master/code/Supply_Model.ipynb) : More detailed development steps of the supply model with different routes that were attempted.  

### Scripts
[`report_plots.py`](https://github.com/Don-Burns/Masters_Project/blob/master/code/report_plots.py) : Script used to generate the results for the sensitivity analysis for the supply model.  The file prefix and parameter value ranges will need to be specified in `__main__` depending on what is being analysed.  Will generate a `.log` file so a script progress can be monitored and will use all availible cores to produce results.  

# [`writing`](https://github.com/Don-Burns/Masters_Project/blob/master/writing)
Contains the project proposal and thesis .tex and pdf files along with the images needed for them which are not generated from the code.

# [`data`](https://github.com/Don-Burns/Masters_Project/blob/master/data)
Contains the outputs of `report_plots.py` for various paramter ranges which were used in `Sensitivity_Plots.ipynb` to generate the sensitivity plots.  The contents will need to be moved to `results` to be used with `Sensitivity_Plots.ipynb`.