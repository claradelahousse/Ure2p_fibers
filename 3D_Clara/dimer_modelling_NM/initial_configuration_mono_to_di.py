# Importing packages and files
import numpy as np

#switching between cylindrical and cartesian coordinates

def carte_to_cyl(carte) :
    x = carte[0]
    y = carte[1]
    z = carte[2]
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y,x)
    return(r,theta,z)

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

# Define number of particles and their positions for initial unit

nb_particles = 18 #number of particles in the unit ie n of subparts of the protein, here number of residues we have inforation about *2 because base motive of fiber = dimer
paticle_radius = 1.9 #used in plot and in gradient
nb_units = 2 #number of copies of the protein, including the original one
nb_max_copies = 4 #maximal number of copies in a super cell to be tested in the first for loop of the Main code


#parameters : p the pitch of the helix and a the roation angle betweeen one protein and the rest, e the rotation angle between two helices

#initialisation :

a10 = 50
b10 = 30
b20 = 50


parameters = [a10,b10,b20]

# Definie intial pqrticle positions

#initialisation :

#N = np.sqrt((m*b20)**2 + (n*a10+ m*b10)**2)

particles_positions_i = np.zeros((nb_particles,3)) #cylindric coordinates r,theta,z

#18 residues

# #mono1

#Nterminal domain
particles_positions_i[0,:] = [10,32,130] #+ np.array([N/(2*np.pi),0,0]) #redisidue 1
particles_positions_i[1,:] = [10,32,120] # + np.array([N/(2*np.pi),0,0]) #redisidue 6
particles_positions_i[2,:] = [10,32,110] #+ np.array([N/(2*np.pi),0,0])  #redisidue 10
particles_positions_i[3,:] = [10,32,100] #+ np.array([N/(2*np.pi),0,0])   #redisidue 13
particles_positions_i[4,:] = [10,32,90] #+ np.array([N/(2*np.pi),0,0])  #redisidue 30
particles_positions_i[5,:] = [10,32,80] #+ np.array([N/(2*np.pi),0,0])  #redisidue 33
particles_positions_i[6,:] = [10,32,70] #+ np.array([N/(2*np.pi),0,0])  #redisidue 52
particles_positions_i[7,:] = [10,32,60] #+ np.array([N/(2*np.pi),0,0])  #redisidue 68
particles_positions_i[8,:] = [10,32,50] #+ np.array([N/(2*np.pi),0,0])  #redisidue 78
particles_positions_i[9,:] = [10,32,40] #+ np.array([N/(2*np.pi),0,0])  #redisidue 92
particles_positions_i[10,:] = [10,32,30] #+ np.array([N/(2*np.pi),0,0])  #redisidue 95


#Cterminal domain
particles_positions_i[11,:] = [9.57,31.83,26.19] # + np.array([N/(2*np.pi),0,0]) #redisidue 104
particles_positions_i[12,:] = [13.24,-16.89,16.61] #+ np.array([N/(2*np.pi),0,0])  #redisidue 121
particles_positions_i[13,:] = [-8.66,-13.10,18.41]  #+ np.array([N/(2*np.pi),0,0]) #redisidue 137
particles_positions_i[14,:] = [1.20,-12.46,-1.43] #+ np.array([N/(2*np.pi),0,0])  #redisidue 158
particles_positions_i[15,:] = [4.33,-6.35,14.06] #+ np.array([N/(2*np.pi),0,0])  #redisidue 181
particles_positions_i[16,:] = [11.80,1.99,15.11] # + np.array([N/(2*np.pi),0,0]) #redisidue 221
particles_positions_i[17,:] = [35.36,-9.33,16.58] #+ np.array([N/(2*np.pi),0,0])  #redisidue 240

#k_spring = 10 #strength of spring interaction
#k_attra = 10 #strength of attractive interaction
#k_rep = 10 #strength of repulsive interaction

li = [5,4,3,17,3,19,16,10,14,3,9]

#Define interactions in the case of one single protein per super cell

# Define interactions inside a protein

#triangle configuration :
interaction_unit_core_i = np.zeros((26,2),dtype = int)
interaction_unit_core_i[:,0] = [0,1,2,3,4,5,6,7,8,9,10,  11,11,11,11,11,11,  12,12,12,12,12,   13,13,13,13] #particles that are in interactions with ...
interaction_unit_core_i[:,1] = [1,2,3,4,5,6,7,8,9,10,11, 12,13,14,15,16,17,  13,14,15,16,17,   14,15,16,17] #these particles
interaction_unit_core_length_i = li #corresponds to zero energy length of vertices when only minimizing a single protein
l_N = 95 #length of the N-terminal in terms of aa
k_N = 1/np.sqrt(l_N)
k_C = 10 #rqndom, change later
interaction_unit_core_strength_i = [k_N,k_N,k_N,k_N,k_N,k_N,k_N,k_N,k_N,k_N,k_N,k_N,] #corresponds to the spring strength in each vertices

#computing the distances in the undefrmer particle and setting the; as d0 for the core interactions in the C-terminal domaim

for i in range(11,14) :
    for j in range(i+1,18) :
        d = np.linalg.norm(particles_positions_i[i,:]-particles_positions_i[j,:])
        interaction_unit_core_length_i.append(d)
        interaction_unit_core_strength_i.append(k_C)


#Define interactions between proteins

#attractive interactions
nb_inter_attra = 15
interaction_inter_attra_i = np.zeros((nb_inter_attra,2))
interaction_inter_attra_i[:,0] = [0,0, 1, 2,2, 3,3, 4, 5, 6, 8,8,8, 10, 11] 
interaction_inter_attra_i[:,1] = [11,17, 13, 12,13, 6,13, 6, 13, 9, 10,12,17, 12, 14] 
interaction_inter_attra_length_i = [12,12, 4, 4,4, 4,4, 4, 4, 4, 12,12,12, 12, 12]
k_attra = 10
interaction_inter_attra_strength_i = [k_attra,k_attra,k_attra,k_attra,k_attra,k_attra,k_attra,k_attra,k_attra,k_attra,k_attra,k_attra,k_attra,k_attra,k_attra]

#repulsive interactions
nb_inter_rep = 1
interaction_inter_rep_i = np.zeros((nb_inter_rep,2))
interaction_inter_rep_i[:,0] = [2]
interaction_inter_rep_i[:,1] = [15]
interaction_inter_rep_length_i = [4]
k_rep = 10
interaction_inter_rep_strength_i = [k_rep]  #NOT USED

#global repulsive interactions
nb_glob_rep = 42
interaction_glob_rep_i = np.zeros((nb_glob_rep,2))
interaction_glob_rep_i[:,0] = [11,11,11,11,11,11,   12,12,12,12,12,12,   13,13,13,13,13,13,   14,14,14,14,14,14,   15,15,15,15,15,15,   16,16,16,16,16,16,   17,17,17,17,17,17]
interaction_glob_rep_i[:,1] = [12,13,14,15,16,17,   11,13,14,15,16,17,   11,12,14,15,16,17,   11,12,13,15,16,17,   11,12,13,14,16,17,   11,12,13,14,15,17,   11,12,13,14,15,16]
g = 4
interaction_glob_rep_length_i = [g,g,g,g,g,g,   g,g,g,g,g,g,   g,g,g,g,g,g,   g,g,g,g,g,g,   g,g,g,g,g,g,   g,g,g,g,g,g,   g,g,g,g,g,g]
k_glob = 100
interaction_glob_rep_strength_i = [k_glob,k_glob,k_glob,k_glob,k_glob,k_glob,   k_glob,k_glob,k_glob,k_glob,k_glob,k_glob,   k_glob,k_glob,k_glob,k_glob,k_glob,k_glob,   k_glob,k_glob,k_glob,k_glob,k_glob,k_glob,   k_glob,k_glob,k_glob,k_glob,k_glob,k_glob,   k_glob,k_glob,k_glob,k_glob,k_glob,k_glob,   k_glob,k_glob,k_glob,k_glob,k_glob,k_glob,]



#Defining the interactions for a super cell using the interactions for a single protein

def super_interaction(nb_copies, interaction_i, interaction_length_i, interaction_strength_i) :    
    l = len(interaction_length_i)
    interaction = np.zeros((l*nb_copies,2))
    interaction_length = []
    interaction_strength = []
    for j in range(nb_copies) :
        interaction[j*l:(j+1)*l,:] = interaction_i + nb_particles*j
        interaction_length = interaction_length + interaction_length_i
        interaction_strength = interaction_strength + interaction_strength_i
    return(interaction,interaction_length,interaction_strength)