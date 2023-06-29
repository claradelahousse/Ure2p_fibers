import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import scipy.optimize as sp
from tqdm import tqdm
#from initial_configuration_mono_to_di import *
from initial_configuration_di_to_fi import *
from Energy import *
from Plotting_functions import *
#from Gradient import *
from Neighbors import *

for n in range(-1,2) :
    for m in range(0,2) :
        nb_copies = 1
    
        nb_interprot_interactions =  (nb_inter_attra + nb_inter_rep)*nb_copies
        particles_positions = find_protein_0(parameters,particles_positions_i,n,m,)
    
        print(particles_positions)
    
        possible_indexes = possible_interprot_indexes(nb_interprot_interactions,1)
    
        print(possible_indexes)
    

        minimisation_parameters = mini_para(particles_positions,nb_copies,parameters)
        interprot_indexes = [1 for k in range(nb_interprot_interactions)]
        arg = (interprot_indexes,nb_copies,n,m)
        print('initial_energy')
        print(total_energy(minimisation_parameters,*arg))
    
        energy_list = []
        configuration_list = []
        l = len(possible_indexes)
        for i in tqdm(range(l)) :
        #for i in tqdm(range(nb_trials)) :
            interprot_indexes = possible_indexes[i]
            #print('interprot', interprot_indexes)
            #interprot_indexes = NN_indexes(i,(nb_units-1),nb_interprot_interactions)
            #print(interprot_indexes)
            arg = (interprot_indexes,nb_copies,n,m)#,fprime=Gradient
            res_full = sp.minimize(total_energy, minimisation_parameters, args=arg,method='Nelder-Mead', options={'tol': 1e-6, 'maxiter' : 10e6})
            configuration_list.append(res_full.x)
            energy_list.append(res_full.fun)
            #print(res_full[1])
            #print(res_full[0])
        
            if res_full.status != 0 :
                print('ERROR')
    
        index_min = np.argmin(energy_list)
        print('final_energy')
        print(energy_list[index_min])
        final_configuration = configuration_list[index_min]
        final_particles_positions = positions_from_mini_para(final_configuration,nb_copies,n,m)
        print(final_particles_positions)
        final_parameters = final_configuration[-3:]
        print(final_parameters)
        final_interprot_indexes = possible_indexes[index_min]
        #final_interprot_indexes = NN_indexes(index_min,(nb_units),nb_interprot_interactions)
        print(final_interprot_indexes)
        print(n,m)

