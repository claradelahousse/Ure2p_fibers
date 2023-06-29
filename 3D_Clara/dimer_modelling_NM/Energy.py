# Importing packages and files
import numpy as np
from initial_configuration_mono_to_di import *
#from initial_configuration_di_to_fi import *
from Neighbors import *

# defining a function that returns the distance betzeen two points in cylindrical coordinates :

def cyl_to_carte(cyl) :
    r = cyl[0]
    theta = cyl[1]
    z = cyl[2]
    x = r*np.cos(theta)
    y = r*np.sin(theta)
    return([x,y,z])


def cyl_distance(c1,c2) :
    dr = c1[0] - c2[0]
    dtheta = c1[1] - c2[1]
    dz = c1[2] - c2[2]
    res = c1[0]**2 + c2[0]**2 - 2*c1[0]*c2[0]*np.cos(dtheta) + dz**2
    return(np.sqrt(res))

#Defining interaction functions :

#spring interaction with arguments : coordinates particle 1, coordinates particle 2, interaction distance.
def interaction_spring(c1,c2,d0,k_spring) :
    d = np.linalg.norm(c1-c2)
    e = (1/2)*k_spring*(d-d0)**2
    return(e)

#harmonic attractive interaction with arguments : coordinates particle 1, coordinates particle 2, interaction cutoff.
def interaction_attra(c1,c2,d0,k_attra) :
    d = np.linalg.norm(c1-c2)
    if d > d0 :
        e = (1/2)*k_attra*(d-d0)**2
    else :
        e = 0
    return(e)

#harmonic repulsive interaction with arguments : coordinates particle 1, coordinates particle 2, interaction cutoff.
def interaction_rep(c1,c2,d0,k_rep) :
    d = np.linalg.norm(c1-c2)
    #print(d)
    if d < d0 :
        e = (1/2)*k_rep*(d-d0)**2
    else :
        e = 0
    return(e)


def force_spring(c1,c2,d0,k_spring) :
    d = np.linalg.norm(c1-c2)
    f = k_spring*(d-d0)/d
    F12 = f*(c1-c2)
    return(F12)

def force_attra(c1,c2,d0,k_attra) :
    d = np.linalg.norm(c1-c2)
    f = k_attra*(d-d0)/d
    if d > d0 :
        F12 = f*(c1-c2)
    else :
        F12 = np.zeros(3)
    return(F12)

def force_rep(c1,c2,d0,k_rep) :
    d = np.linalg.norm(c1-c2)
    f = k_rep*(d-d0)/d
    if d < d0 :
        F12 = f*(c1-c2)
    else :
        F12 = np.zeros(3)
    return(F12)


#Defining function that returns all the minimisation parameters in a 1D array. used for gradient as well

def mini_para(particles_positions,nb_copies,parameters) :
    positions_1D = np.reshape(particles_positions[1:,:],((nb_particles-1)*nb_copies*3))
    #print(positions_1D)
    minimisation_parameters = np.concatenate((positions_1D ,parameters))
    return(minimisation_parameters)

#defininf function that returns particles_positions from minimisation parameters

def positions_from_mini_para(mimisation_parameters,nb_copies,n,m) :
    a1 = mimisation_parameters[-3]
    b1 = mimisation_parameters[-2]
    b2 = mimisation_parameters[-1]
    N = np.sqrt((m*b2)**2 + (n*a1 + m*b1)**2) #N/(2*np.pi)
    R = N/(2*np.pi)
    particle_0 = np.array([R,0,0])
    particles_positions = np.concatenate((particle_0,mimisation_parameters))
    particles_positions = np.reshape(particles_positions[:3*nb_particles*nb_copies],(nb_particles*nb_copies,3))
    return(particles_positions)

def total_energy(minimisation_parameters,*args) :
    (interprot_indexes,nb_copies,n,m) = args
    #Defining the interactions in the super cells
    interaction_unit_core,interaction_unit_core_length,interaction_unit_core_strength = super_interaction(nb_copies,interaction_unit_core_i,interaction_unit_core_length_i,interaction_unit_core_strength_i)
    #interaction_unit_attra,interaction_unit_attra_length,interaction_unit_attra_strength = super_interaction(nb_copies,interaction_unit_attra_i,interaction_unit_attra_length_i,interaction_unit_attra_strength_i)
    #interaction_unit_rep,interaction_unit_rep_length,interaction_unit_rep_strength = super_interaction(nb_copies,interaction_unit_rep_i,interaction_unit_rep_length_i,interaction_unit_rep_strength_i)
    interaction_inter_attra,interaction_inter_attra_length,interaction_inter_attra_strength = super_interaction(nb_copies,interaction_inter_attra_i,interaction_inter_attra_length_i,interaction_inter_attra_strength_i)
    interaction_inter_rep,interaction_inter_rep_length,interaction_inter_rep_strength = super_interaction(nb_copies,interaction_inter_rep_i,interaction_inter_rep_length_i,interaction_inter_rep_strength_i)
    interaction_glob_rep,interaction_glob_rep_length,interaction_glob_rep_strength = super_interaction(nb_copies,interaction_glob_rep_i,interaction_glob_rep_length_i,interaction_glob_rep_strength_i)
    #print(interaction_unit_core)
    #Getting particles positions and parameters from minimisation parameters
    particles_positions = positions_from_mini_para(minimisation_parameters,nb_copies,n,m)
    #print(particles_positions)
    parameters = minimisation_parameters[-3:]
    #COMPUTE ALL protein positions. 
    interacting = np.unique(interprot_indexes)
    list_proteins = [particles_positions]
    for k in interacting :
        list_proteins.append(find_protein_x(parameters,particles_positions,k,n,m))
    en = 0
    #intraprotein
    #unit interactions 
    l = np.size(interaction_unit_core_length)
    for i in range(l) :
        p1 = int(interaction_unit_core[i,0]) #first particle in interaction ...
        p2 = int(interaction_unit_core[i,1]) #... with second particle
        #print(p1)
        c1 = particles_positions[p1,:]
        c2 = particles_positions[p2,:]
        #print(i)
        d0 = interaction_unit_core_length[i]
        k_spring = interaction_unit_core_strength[i]
        en += interaction_spring(c1,c2,d0,k_spring)
    # #attractive interactions
    # l = np.size(interaction_unit_attra_length)
    # for i in range(l) :
    #     p1 = int(interaction_unit_attra[i,0]) #first particle in interaction ...
    #     p2 = int(interaction_unit_attra[i,1]) #... with second particle
    #     #print(p1)
    #     c1 = particles_positions[p1,:]
    #     c2 = particles_positions[p2,:]
    #     d0 = interaction_unit_attra_length[i]
    #     k_attra = interaction_unit_attra_strength[i]
    #     en += interaction_attra(c1,c2,d0,k_attra)
    # #repulsive interactions
    # l = np.size(interaction_unit_rep_length)
    # for i in range(l) :
    #     p1 = int(interaction_unit_rep[i,0]) #first particle in interaction ...
    #     p2 = int(interaction_unit_rep[i,1]) #... with second particle
    #     #print(p1)
    #     c1 = particles_positions[p1,:]
    #     c2 = particles_positions[p2,:]
    #     d0 = interaction_unit_rep_length[i]
    #     k_rep = interaction_unit_rep_strength[i]
    #     en += interaction_rep(c1,c2,d0,k_rep)
    # # interprotein
    count = 0 # counting the number of interprotein interactions taken care of in order to une the interprot_indexes list.
    #attractive interactions
    l = np.size(interaction_inter_attra_length)
    #print(count)
    for i in range(l) :
        p1 = int(interaction_inter_attra[i,0]) #first particle in interaction ...
        p2 = int(interaction_inter_attra[i,1]) #... with second particle
        #print(count)
        unit_index = interprot_indexes[count] #- nb_units//2
        c2 = list_proteins[unit_index][p2,:]
        c1 = particles_positions[p1,:]
        d0 = interaction_inter_attra_length[i]
        k_attra = interaction_inter_attra_strength[i]
        count += 1
        #print(interaction_attra(c1,cNN,d0))
        en += (1/2)*interaction_attra(c1,c2,d0,k_attra)
    # # # #repulsive interactions
    count = 5
    l = np.size(interaction_inter_rep_length)
    for i in range(l) :
        p1 = int(interaction_inter_rep[i,0]) #first particle in interaction ...
        p2 = int(interaction_inter_rep[i,1]) #... with second particle
        unit_index = interprot_indexes[count] #- nb_units//2 ###ATTENTION REMETTRE COUNT
        c2 = list_proteins[unit_index][p2,:]
        c1 = particles_positions[p1,:]
        d0 = interaction_inter_rep_length[i]
        k_rep = interaction_inter_rep_strength[i]
        count += 1
        en += (1/2)*interaction_rep(c1,c2,d0,k_rep)
    # GLOBAL rep
    interacting = np.unique(interprot_indexes)
    for k in interacting :
        l = np.size(interaction_glob_rep[:,0])
        for i in range(l) :
            p1 = int(interaction_glob_rep[i,0]) #first particle in interaction ...
            p2 = int(interaction_glob_rep[i,1]) #... with second particle
            unit_index = k
            c2 = list_proteins[unit_index][p2,:]
            c1 = particles_positions[p1,:]
            d0 = interaction_glob_rep_length[i]
            k_rep = interaction_glob_rep_strength[i]
            en += (1/2)*interaction_rep(c1,c2,d0,k_rep)
        #print(interaction_rep(c1,c2,d0,k_rep))
    # #print(e)
    return(en)