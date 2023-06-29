# Importing packages and files
import numpy as np
import itertools
from initial_configuration_mono_to_di import *
#from initial_configuration_di_to_fi import *

#Neighbors numberss can either be positive or negative and describe how the helical transformation is applied.

#Defining function that gives out 3D rotation matrix in xy plane.

ab = np.zeros((43,2))
ab[0,:] = [0,0]#[0,0]
ab[1,:] = [1,0]#[1,0]
ab[2,:] = [-1,0]#[0,1] #[0,-1]
ab[3,:] = [0,1]#[1,1]#[1,0]
ab[4,:] = [0,-1]#[2,0]#[-1,0]
ab[5,:] = [1,1]#[0,2]#[1,1]
ab[6,:] = [-1,-1]#[2,1]#[1,-1]
ab[7,:] = [1,-1]#[1,2]#[-1,-1]
ab[8,:] = [-1,1]#[2,2]#[-1,1]
ab[9,:] = [2,0]#[3,0]#[2,0]
ab[10,:] = [0,2]#[0,3]#[-2,0]
ab[11,:] = [-2,0]#[1,3]#[0,2]
ab[12,:] = [0,-2]#[3,1]#[2,1]
ab[13,:] = [2,-1]#[2,-1]
ab[14,:] = [-2,1]
ab[15,:] = [-2,-2]
ab[16,:] = [1,2]
ab[17,:] = [-1,2]
ab[18,:] = [1,-2]
ab[19,:] = [-1,-2]
ab[20,:] = [2,2]
ab[21,:] = [2,-2]
ab[22,:] = [-2,3]
ab[23,:] = [-2,2]
ab[24,:] = [3,0]
ab[25,:] = [-3,0]
ab[26,:] = [0,3]
ab[27,:] = [0,-3]
ab[28,:] = [3,-1]
ab[29,:] = [-3,1]
ab[30,:] = [-3,-1]
ab[31,:] = [3,1]
ab[32,:] = [-1,3]
ab[33,:] = [1,-3]
ab[34,:] = [-1,-3]
ab[35,:] = [1,3]
ab[36,:] = [3,-2]
ab[37,:] = [-2,3]
ab[38,:] = [2,3]
ab[39,:] = [-3,2]
ab[40,:] = [-3,-2]
ab[41,:] = [-2,-3]
ab[42,:] = [2,-3]

def number_to_coordinates(nb_units,unit_index) :
    combinations = list(itertools.product([k for k in range(0,nb_units%4+1,1)], repeat=2))
    #i,j = combinations[unit_index]
    #print(i,j)
    i = ab[unit_index,0]
    j = ab[unit_index,1]
    return(i,j)

def cylinder_parameters(parameters,n,m,i,j) :
    a1 = parameters[0]
    b1 = parameters[1]
    b2 = parameters[2]
    N = np.sqrt((m*b2)**2 + (n*a1 + m*b1)**2)
    phi = ((i*a1+j*b1)*(n*a1+m*b1)+m*j*b2**2)*2*np.pi/N**2  #%(2*np.pi)
    e = (-(i*a1+j*b1)*m*b2+j*b2*(n*a1+m*b1))/N
    return(phi,e)

def transformation(phi,e) :  
    R_matrix = np.reshape(np.array([np.cos(phi),np.sin(phi),0,0,-np.sin(phi),np.cos(phi),0,0,0,0,1,e,0,0,0,1]),(4,4)) # rotation matrix around axix z + translation of direction z
    return(R_matrix)


#Defines a function that computes the coordinates of a particle in a neighbor unit from the unit number and the helix parameters in cylindric coordinates 


def find_protein_0(parameters,particles_positions_i,n,m,) :
    a1 = parameters[0]
    b1 = parameters[1]
    b2 = parameters[2]
    N = np.sqrt((m*b2)**2 + (n*a1 + m*b1)**2) #N/(2*np.pi)
    R = N/(2*np.pi)
    protein_0_positions = np.zeros((nb_particles,3))
    x0 = particles_positions_i[0,0]
    y0 = particles_positions_i[0,1]
    z0 = particles_positions_i[0,2]
    for i in range(nb_particles) :
        protein_0_positions[i,:] = particles_positions_i[i,:] + np.array([R-x0,-y0,-z0])
    return(protein_0_positions)


def find_neighbor_coordinates(parameters,protein_0_positions, residue, unit_index,n,m) :
    i,j = number_to_coordinates(nb_units, unit_index)
    a1 = parameters[0]
    b1 = parameters[1]
    b2 = parameters[2]
    N = np.sqrt((m*b2)**2 + (n*a1 + m*b1)**2) 
    phi,e = cylinder_parameters(parameters,n,m,i,j)
    #print('e',e)
    v = np.concatenate((protein_0_positions[residue,:],np.array([1])))
    neighbor_position = np.dot(transformation(phi,e),v)
    #print('v',v)
    #print(neighbor_position)
    return(neighbor_position[:-1])

def find_protein_x(parameters,protein_0_positions,unit_index,n,m) :
    protein_x_positions = np.zeros((nb_particles,3))
    for residue in range(nb_particles) :
        protein_x_positions[residue:] = find_neighbor_coordinates(parameters,protein_0_positions, residue, unit_index,n,m)
    if unit_index == 0 :
        return(protein_0_positions)
    else :  
        return(protein_x_positions)

# def find_protein_coordinates(parameters,particles_positions, unit_index,n,m) :
#     protein_position = np.zeros((nb_particles,3))
#     i,j = number_to_coordinates(nb_units, unit_index)
#     a1 = parameters[0]
#     b1 = parameters[1]
#     b2 = parameters[2]
#     N = np.sqrt((m*b2)**2 + (n*a1 + m*b1)**2)
#     N0 = np.sqrt((m*b20)**2 + (n*a10+ m*b10)**2)
#     for i in range(nb_particles) : #N/(2*np.pi)
#         protein_position[i,:] = particles_positions[i,:] + np.array([N/(2*np.pi) - N0/(2*np.pi),((i*a1+j*b1)*(n*a1+m*b1)+m*j*b2**2)*2*np.pi/N**2,(b2*a1*(j*n-i*m))/N])
#     return(protein_position)

#Defining a function that creates the supercell by stacking the coordinates of each particle for each copy of the protein in the supercell

# ATTENTION marche pour 2 copies mais pas encore teste pour plus, choses a changer a partir de boucle for j in range(2,nb_copies)

def create_super_cell(particles_positions_i,nb_copies) : #a : distance between 2 copies
    if nb_copies == 1 :
        return(particles_positions_i)
    sym_particles_positions = np.zeros((nb_particles,3))
    for i in range(nb_particles) :
            sym_particles_positions[i,:] = particles_positions_i[i,:]*np.array([1,-1,-1]) #+ np.array([0,2*np.pi/nb_copies,np.max(particles_positions_i[:,2])])
    particles_positions = np.concatenate((particles_positions_i,sym_particles_positions))
    return(particles_positions)
    
#Define a function that returns the indexes of the proteins to be considered neighbors.

def possible_interprot_indexes(nb_interprot_interactions, farthest_neighbor) :
    combinations1 = list(itertools.product([k  for k in range(1,farthest_neighbor + 1)], repeat=nb_interprot_interactions))
    #combinations2 = list(itertools.product([k for k in range(-farthest_neighbor,0)], repeat=nb_interprot_interactions))
    #combinations = combinations1 + combinations2
    combinations = combinations1
    return(combinations)