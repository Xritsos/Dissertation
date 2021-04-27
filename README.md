# Joule Heating Error Analysis
[![forthebadge made-with-python](http://ForTheBadge.com/images/badges/made-with-python.svg)](https://www.python.org/)

[![PyPI pyversions](https://upload.wikimedia.org/wikipedia/commons/3/34/Blue_Python_3.6_Shield_Badge.svg)](https://pypi.python.org/pypi/ansicolortags/)

The code in this repository is used in my dissertation: "Error analysis of Joule heating calculations in the Upper Atmosphere from in-situ satellite measurements". It is implemented under the Daedalus project, a mission proposed to ESA for the Earth Observation programme’s 10th Earth Explorer. The version one corresponds to the error propagation analysis, while version two deals with random errors or else noise. The measurements are simulated using data from the TIE-GCM(2.0) and IGRF(12) models. 
The data are not provided here beacause of large file sizes but anyone can request a run for the models in the links below:

* https://ccmc.gsfc.nasa.gov/requests/IT/TIE-GCM/tiegcm_user_registration.php

* http://www.geomag.bgs.ac.uk/data_service/models_compass/igrf_calc.html

In this case the data for the IGRF model where not request but imported with the pyglow library (including more atmospheric models), which can be found here: https://github.com/timduly4/pyglow.git
The operating system used is Ubuntu 18.04.5 LTS. More prerequisites can be found at `requirements.txt`
