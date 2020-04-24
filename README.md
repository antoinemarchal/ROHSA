# ROHSA

ROHSA (Regularized Optimization for Hyper-Spectral Analysis) was developped by the Hyperstars collaboration at Paris-Saclay University (IAS/CEA) to study the statistical properties of interstellar gas through atomic and molecular lines.

This code is a Gaussian decomposition algorithm designed to decompose any kind of hyper-spectral observations into a sum of coherent Gaussian. It is written in Fortran 90 and can be run on a single CPU. A user-friendly ROHSApy python interface can be used to run the code easily.

## Publication
[Marchal et al., A&A 626, A101 (2019)](https://ui.adsabs.harvard.edu/abs/2019A%26A...626A.101M/abstract)
Please cite this publication if you are using ROHSA. 

## Replicating the results
All codes used in this work are available [here](https://github.com/antoinemarchal/ROHSA/tree/master/publication)

## Documentation
Check out the [online documentation](https://antoinemarchal.github.io/ROHSA/lists/files.html) generated by [FORD](https://github.com/cmacmackin/ford).

## Getting started
A google colab notebook is available [here!](https://github.com/antoinemarchal/ROHSA/blob/master/ROHSApy.ipynb) It shows how to install, compile, run and read ROHSA's output.

Please note that if you are not using ROHSA for separating phases of 21cm observation, the hyper-paramters _"lambda_var_sig"_** should be set to 0. 

## CONTACT 
If you have any queries or questions related to ROHSA, please do not hesitate to contact us [here](amarchal@cita.utoronto.ca). Your feedback is welcome and important to us. 


