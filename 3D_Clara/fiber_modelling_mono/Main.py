import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import scipy.optimize as sp
from tqdm import tqdm
#from initial_configuration_mono_to_di import *
from initial_configuration_di_to_fi import *
from Energy import *
from Plotting_functions import *
from Gradient import *
from Neighbors import *

n=0
m=1

for nb_copies in range(1,2) :
    
    nb_interprot_interactions =  (nb_inter_attra + nb_inter_rep)*nb_copies
    #print(nb_interprot_interactions)
    particles_positions = find_protein_0(parameters,particles_positions_i,n,m,)
    
    print(particles_positions)
    
    possible_indexes = possible_interprot_indexes(nb_interprot_interactions,1)
    
    print(possible_indexes)
    
#     nb_trials = ((nb_units-1))**(nb_interprot_interactions) 
#     print(nb_trials)

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
        print('interprot', interprot_indexes)
        #interprot_indexes = NN_indexes(i,(nb_units-1),nb_interprot_interactions)
        #print(interprot_indexes)
        arg = (interprot_indexes,nb_copies,n,m)#,fprime=Gradient
        res_full = sp.fmin_cg(total_energy, minimisation_parameters,fprime=Gradient, args=arg,gtol=1e-06, full_output=1, disp=0)
        configuration_list.append(res_full[0])
        energy_list.append(res_full[1])
        #print(res_full[1])
        #print(res_full[0])
        
        if res_full[4] != 0 :
            g = Gradient(configuration_list[-1],*arg)
            g_num = Gradient_Numerical(configuration_list[-1],*arg)
            if np.max(np.abs(g-g_num))>1e-3 :
                print(g)
                print(g_num)
                print('error')
                print(res_full[4])
                print(i)
    
    index_min = np.argmin(energy_list)
    print('final_energy')
    print(energy_list[index_min])
    final_configuration = configuration_list[index_min]
    print(final_configuration)
    final_particles_positions = positions_from_mini_para(final_configuration,nb_copies)
    final_parameters = final_configuration[-3:]
    print(final_parameters)
    final_interprot_indexes = possible_indexes[index_min]
    #final_interprot_indexes = NN_indexes(index_min,(nb_units),nb_interprot_interactions)
    print(final_interprot_indexes)

def protein_rpz(e,a,title) :
    # Create the figure and 3D axes
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #ax.set_box_aspect([1, 1, 1]) 
    # ax.set_xlim([-3,3])
    # ax.set_ylim([-3,3])
    # ax.set_zlim([-3,3])
    
    #c_all = ['pink','magenta','purple','deeppink']

    interaction_unit_core,interaction_unit_core_length,interaction_unit_core_strength = super_interaction(nb_copies,interaction_unit_core_i,interaction_unit_core_length_i,interaction_unit_core_strength_i)
    interaction_inter_attra,interaction_inter_attra_length,interaction_inter_attra_strength = super_interaction(nb_copies,interaction_inter_attra_i,interaction_inter_attra_length_i,interaction_inter_attra_strength_i)
    interaction_inter_rep,interaction_inter_rep_length,interaction_inter_rep_strength = super_interaction(nb_copies,interaction_inter_rep_i,interaction_inter_rep_length_i,interaction_inter_rep_strength_i)

    
    for u in range(nb_units) :
        for p in range(nb_particles*nb_copies) :
            c = c_all #[p]
            carte = cyl_to_carte(find_neighbor_coordinates(final_parameters,final_particles_positions,p,u,n,m))
            X,Y,Z = surface_map(carte,paticle_radius)
            ax.plot_surface(X,Y,Z, color=c,shade = True)
        
        l = len(interaction_unit_core_length)
        for i in range(l) :
            point1 = cyl_to_carte(find_neighbor_coordinates(final_parameters,final_particles_positions,int(interaction_unit_core[i,0]),u,n,m))
            point2 = cyl_to_carte(find_neighbor_coordinates(final_parameters,final_particles_positions,int(interaction_unit_core[i,1]),u,n,m))
            r = np.linalg.norm(np.array(point1)-np.array(point2))
            linex = np.linspace(point1[0],point2[0],10)
            liney = np.linspace(point1[1],point2[1],10)
            linez = np.linspace(point1[2],point2[2],10)
            ax.plot(linex,liney,linez,'grey')
            #print(r)

    l = len(interaction_inter_attra_length)
    print(l)
    for i in range(l) :
        point1 = cyl_to_carte(find_neighbor_coordinates(final_parameters,final_particles_positions,int(interaction_inter_attra[i,0]),0,n,m))
        point2 = cyl_to_carte(find_neighbor_coordinates(final_parameters,final_particles_positions,int(interaction_inter_attra[i,1]),final_interprot_indexes[i],n,m))#-nb_units//2))
        r = np.linalg.norm(np.array(point1)-np.array(point2))
        linex = np.linspace(point1[0],point2[0],10)
        liney = np.linspace(point1[1],point2[1],10)
        linez = np.linspace(point1[2],point2[2],10)
        ax.plot(linex,liney,linez,'tomato')
        #print(r)

    l = len(interaction_inter_rep_length)
    print(l)
    for i in range(l) :
        point1 = cyl_to_carte(find_neighbor_coordinates(final_parameters,final_particles_positions,int(interaction_inter_rep[i,0]),0,n,m))
        point2 = cyl_to_carte(find_neighbor_coordinates(final_parameters,final_particles_positions,int(interaction_inter_rep[i,1]),final_interprot_indexes[i+nb_inter_attra],n,m))#-nb_units//2))
        r = np.linalg.norm(np.array(point1)-np.array(point2))
        linex = np.linspace(point1[0],point2[0],10)
        liney = np.linspace(point1[1],point2[1],10)
        linez = np.linspace(point1[2],point2[2],10)
        ax.plot(linex,liney,linez,'cornflowerblue')
        #print(r)

    # Set the plot labels
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # Set the plot title
    ax.set_title(title)
    
    # Set the initial view
    elev = e
    azim = a
    rotate_view(elev, azim,ax,fig)
    
    plt.show()

dossier = 'fibers'
protein_rpz(30,30,'dimer')
plt.savefig(dossier+'/fig'+str(0)+'.png')