# R code for finding activity sensitive spectral lines

Stellar activity including spots and plages can create signals that are similarto planetary RV signals, making them particularly problematic in planet searches. As some pixels (or absorption lines) are more sensitive to stellar activity than others, identifying and isolating these pixels from an analysis could help usmeasure RV shifts that are truly due to exoplanets. In this file, we provide the R code for our proposed method to identifying activity sensitive spectral linesa -- using a Bayesian variable selection approach to fast and automatically  search  for  activity-sensitive lines. 

The code folder contains all the functions for the Bayesian variable selection algorithm. The "main.R" is the main function. The "example.R" provides an example on how to use the "main.R" function using the data from alpha Centrauri B. 

The results folder includes the pixels selected using S-index and H alpha, NaD, BIS, FWHM indices. 

The data folder provides the wavelengths, the activity indices, and Julian dates of the alpha Centrauri B dataset. We are unable to upload the normalized spectrum data as the file size is too large. However, a copy is avaliable by sending us a request at bo.ning@yale.edu or aww@udel.edu.
