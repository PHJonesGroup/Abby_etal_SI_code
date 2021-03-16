Simulations of clonal dynamics for Abby et al.

These use a Wright-Fisher style process constrained to a two-dimensional grid to simulate clonal growth in the basal layer of esophageal epithelium.   

The scripts require Python 3 and use the code from https://github.com/michaelhall28/clone-competition-simulation to run simulations.

Other required packages (numpy, scipy, pandas and matplotlib are also required for the installation of the clone-competition-simulation code):
-	numpy  
-	scipy  
-	pandas  
-	matplotlib  
-	pyabc
-	openpyxl
-	jupyter


Run the fitting of the model to data using   
``` python abc_fitting.py data_type ```   
where `data_type` is one of 'wt', 'het_ctl', 'hom_ctl', 'het' or 'hom'.  

The fitting results will be output into a database, which can be explored using PyABC functions (see pyabc.readthedocs.io). The figures from the paper can be recreated using the notebook Fitting_results.ipynb. Note there may be some small changes due to the random nature of the simulations.  

Simulations demonstrating the impact of the haploinsufficiency of Notch1 can be run using   
```python haploinsufficiency_sims.py```   
