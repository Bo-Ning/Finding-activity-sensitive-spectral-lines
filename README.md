# R code for finding activity sensitive spectral lines

Stellar activity including spots and plages may lead to false discoveries or poor mass estimates for small exoplanets when using radial velocity (RV) techniques. In this file, we provide the R code for a recently proposed method to identifying activity sensitive spectral linesa -- using a Bayesian variable selection approach to fast and automatically  search  for  activity-sensitive lines. 

The code file contains all the functions for the Bayesian variable selection algorithm. The "main.R" is the main function. The "example.R" provides an example on how to use the "main.R" function using the data from alpha Centrauri B. 

The data folder provides the wavelengths, the activity indices, and Julian dates of the alpha Centrauri B dataset. We are unable to upload the normalized spectrum data as the file size is too large. However, a copy is avaiable to obtain by sending us a request at bo.ning@yale.edu or aww@udel.edu.
