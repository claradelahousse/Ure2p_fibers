# Ure2p_fibers
internship April-June 2023 LPTMS - modelling fiber formation during self assembly of Ure2p prion.


Organisation of one folder:

    • Initial_configuration file : contains the coordinates of the residues of the initial protein, and the interactions. To change configurations modify only this file.
    • Neighbors file : contains all the functions used to generate the coordinates of neighbor proteins
    • Plotting file : contains all the functions necessary to plot the outcome configuration. It is necessary to copy and paste these functions into a jupyter notebook if you want to interact with the figures.
    • Energy file : contains the Energy function that computes the energy of the initial protein in a given fiber configuration, and all the functions it calls. 
    •  Gradient file : contains a Numerical_Gradient function and a Gradient function. Usefull only when the minimization is carried on with a conjugate gradient method. The Numercical_Gradient function only serves to test the implementation of the Gradient function.
    • Launcher file : used to run the code on the cluster. 
    • Main file : loops that test all the fiber configurations chosen, minimization of the continuous parameters either in a Nelder Mead or a Conjugate Gradient method, returns the optimized parameters and minimum energy for each fiber configuration tested.

the different folders :

    • dimer_modelling_NM : folder used to get the dimer conformation out of 2 native proteins. Uses a Nerlde Mead method for optimisation of the continuous parameters.
    • fiber_modelling_dim_NM : folder used to get fibers out of the dimers found with the file above. Uses a Nerlde Mead method for optimisation of the continuous parameters.
    • fiber_modelling_mono : folder used to get fibers out of the native monomers. Uses a Conjugate Gradient method for optimisation of the continuous parameters.
    • fiber_modelling_mono_NM : folder used to get fibers out of the native monomers. Uses a Nelder Mead method for optimisation of the continuous parameters.

If you want to remove the effect of one interaction you can simply comment it out in the energy and gradient functions, without modifying the intial_configuration file.

To run the code in the cluster : python3 -u launcher.py

To run the code locally : python3 Main.py
