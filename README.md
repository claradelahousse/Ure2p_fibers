# Ure2p_fibers
internship April-June 2023 LPTMS - modelling fiber formation during self assembly of Ure2p prion.


Organisation :

    • Initial_configuration file : contains the coordinates of the residues of the initial protein, and the interactions. To change configurations modify only this file.
    • Neighbors file : contains all the functions used to generate the coordinates of neighbor proteins
    • Plotting file : contains all the functions necessary to plot the outcome configuration. It is necessary to copy and paste these functions into a jupyter notebook if you want to interact with the figures.
    • Energy file : contains the Energy function that computes the energy of the initial protein in a given fiber configuration, and all the functions it calls. 
    •  Gradient file : contains a Numerical_Gradient function and a Gradient function. Usefull only when the minimization is carried on with a conjugate gradient method. The Numercical_Gradient function only serves to test the implementation of the Gradient function.
    • Launcher file : used to run the code on the cluster. 
    • Main file : loops that test all the fiber configurations chosen, minimization of the continuous parameters either in a Nelder Mead or a Conjugate Gradient method, returns the optimized parameters and minimum energy for each fiber configuration tested.

If you want to remove the effect of one interaction you can simply comment it out in the energy and gradient functions, without modifying the intial_configuration file.

To run the code : python3 -u launcher.py
