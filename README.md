# Joule Heating Error Analysis
[![forthebadge made-with-python](http://ForTheBadge.com/images/badges/made-with-python.svg)](https://www.python.org/)

[![PyPI pyversions](https://upload.wikimedia.org/wikipedia/commons/3/34/Blue_Python_3.6_Shield_Badge.svg)](https://pypi.python.org/pypi/ansicolortags/)

## About the project

The code in this repository is used in my dissertation: "Error analysis of Joule heating calculations in the Upper Atmosphere from in-situ satellite measurements". It is implemented under the Daedalus project, a mission proposed to ESA for the Earth Observation programme’s 10th Earth Explorer. This work is used to demonstrate the effects of error in satellite measurements. Although, the data used do not come from in-situ measurements, but atmospheric models are used to simulate the measured data instead. The models used are the TIE-GCM (V2.0) and IGRF-12. More information about the Daedalus mission and the errors used can be found at the mission's web page: https://daedalus.earth.
In this repository two versions of the code can be found in the folders:
* `Joule Heating Error Propagation` (version 1)
* `Joule Heating Statistical Error` (version 2)

In the first version, the error propagation method is used to determine the systematic error, based on the errors proposed by the Daedalus team. The second version deals with random error or else refered as noise or statistical error. Both versions calculate some key ionospheric quantities and their errors, such as:
* Heating rates
* Collision frequencies
* Conductivities
* Currents and
* Cross sections.
 
The output plots depict the quantities and errors as height, latitude-longitude and altitude-latitude depended.

## Prerequisites

The data used in the mission are not provided here beacause of large file sizes but anyone can request a run for the models in the links below:

* https://ccmc.gsfc.nasa.gov/requests/IT/TIE-GCM/tiegcm_user_registration.php

* http://www.geomag.bgs.ac.uk/data_service/models_compass/igrf_calc.html.

In this case the data for the IGRF model where not requested but imported with the pyglow library (including more atmospheric models), which can be found here: https://github.com/timduly4/pyglow.git. In this link the installation guide for the pyglow is provided.

## Built with

The operating system used is Ubuntu 18.04.5 LTS and the python version is 3.6.9. More about module versions used can be found at `requirements.txt`.

## License

Distributed under the GPL-3.0 License. See more at `LICENSE`.

## Contact

I would be glad to answer your questions and listen to your suggestions. Also any correction would be highly appreciated.

Psychalas Christos - kobecris24@yahoo.gr
Project Link: https://github.com/Xritsos/Dissertation.git
