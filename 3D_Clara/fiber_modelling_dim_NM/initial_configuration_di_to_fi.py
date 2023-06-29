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

nb_particles = 2*18 #number of particles in the unit ie n of subparts of the protein, here number of residues we have inforation about *2 because base motive of fiber = dimer
paticle_radius = 1.9 #used in plot and in gradient
nb_units = 13 #number of copies of the protein, including the original one
nb_max_copies = 4 #maximal number of copies in a super cell to be tested in the first for loop of the Main code


#parameters : p the pitch of the helix and a the roation angle betweeen one protein and the rest, e the rotation angle between two helices

#initialisation :

a10 = 100
b10 = 50
b20 = 100


parameters = [a10,b10,b20]

# Definie intial pqrticle positions

#initialisation :

#N = np.sqrt((m*b20)**2 + (n*a10+ m*b10)**2)

particles_positions_i = np.zeros((nb_particles,3)) #cylindric coordinates r,theta,z

#18 residues

# #mono1

#Nterminal domain
particles_positions_i[0,:] = [3.23426869e+01,  0.00000000e+00,  0.00000000e+00] #+ np.array([N/(2*np.pi),0,0]) #redisidue 1
particles_positions_i[1,:] = [4.39277942e+00 , 4.37539301e-04 ,-2.42526313e+01] # + np.array([N/(2*np.pi),0,0]) #redisidue 6
particles_positions_i[2,:] = [9.57533818e+00  ,4.85677544e-03, -2.47669076e+01] #+ np.array([N/(2*np.pi),0,0])  #redisidue 10
particles_positions_i[3,:] = [6.23037814e+00,  1.64823145e-02 ,-1.93994429e+01] #+ np.array([N/(2*np.pi),0,0])   #redisidue 13
particles_positions_i[4,:] = [1.51541105e+01 , 2.90124186e-03, -1.38233874e+01] #+ np.array([N/(2*np.pi),0,0])  #redisidue 30
particles_positions_i[5,:] = [5.37838706e+00  ,1.35593688e-02 ,-2.15331379e+01] #+ np.array([N/(2*np.pi),0,0])  #redisidue 33
particles_positions_i[6,:] = [1.08964875e+01,  1.22310210e-02, -4.81757710e+01] #+ np.array([N/(2*np.pi),0,0])  #redisidue 52
particles_positions_i[7,:] = [-9.48926710e+00, -1.18168037e-03, -2.59740350e+01] #+ np.array([N/(2*np.pi),0,0])  #redisidue 68
particles_positions_i[8,:] = [2.82158068e+01  ,1.42470351e-02, -1.86147668e+01] #+ np.array([N/(2*np.pi),0,0])  #redisidue 78
particles_positions_i[9,:] = [1.33602638e+01, -1.28611363e-02 ,-7.68770294e+01] #+ np.array([N/(2*np.pi),0,0])  #redisidue 92
particles_positions_i[10,:] = [2.18508798e+01,  2.18409926e-02 ,-4.04983795e+01] #+ np.array([N/(2*np.pi),0,0])  #redisidue 95


#Cterminal domain
particles_positions_i[11,:] = [2.36874414e+01 , 3.26311943e-01, -1.34982683e+01] # + np.array([N/(2*np.pi),0,0]) #redisidue 104
particles_positions_i[12,:] = [1.99894734e+01  ,5.46784176e-01 ,-6.08408478e+01] #+ np.array([N/(2*np.pi),0,0])  #redisidue 121
particles_positions_i[13,:] = [-3.44339542e-01, -1.69247442e+00 ,-5.60243631e+01]  #+ np.array([N/(2*np.pi),0,0]) #redisidue 137
particles_positions_i[14,:] = [1.56613350e+1, -1.85098470e+01, -5.87020871e+01] #+ np.array([N/(2*np.pi),0,0])  #redisidue 158
particles_positions_i[15,:] = [1.26494814e+01 ,-7.48788495e+00, -5.17110278e+01] #+ np.array([N/(2*np.pi),0,0])  #redisidue 181
particles_positions_i[16,:] = [1.99269003e+01  ,8.92149449e+00 ,-4.39817430e+01] # + np.array([N/(2*np.pi),0,0]) #redisidue 221
particles_positions_i[17,:] = [4.04536792e+01, -1.11435213e+01, -5.33063288e+01] #+ np.array([N/(2*np.pi),0,0])  #redisidue 240

#Nterminal domain
particles_positions_i[18,:] = [3.16230962e+01 , 6.78448072e+00,  3.34469594e+01] #+ np.array([N/(2*np.pi),0,0]) #redisidue 1
particles_positions_i[19,:] = [4.29495293e+00  ,9.21895071e-01 , 9.19432816e+00] # + np.array([N/(2*np.pi),0,0]) #redisidue 6
particles_positions_i[20,:] = [9.36127822e+00,  2.01335417e+00  ,8.68005185e+00] #+ np.array([N/(2*np.pi),0,0])  #redisidue 10
particles_positions_i[21,:] = [6.08830133e+00 , 1.32305341e+00,  1.40475166e+01] #+ np.array([N/(2*np.pi),0,0])   #redisidue 13
particles_positions_i[22,:] = [1.48163389e+01  ,3.18169350e+00 , 1.96235720e+01] #+ np.array([N/(2*np.pi),0,0])  #redisidue 30
particles_positions_i[23,:] = [5.25587929e+00,  1.14147451e+00  ,1.19138215e+01] #+ np.array([N/(2*np.pi),0,0])  #redisidue 33
particles_positions_i[24,:] = [1.06514865e+01 , 2.29770001e+00 ,-1.47288116e+01] #+ np.array([N/(2*np.pi),0,0])  #redisidue 52
particles_positions_i[25,:] = [-9.27789305e+00 ,-1.99170583e+00 , 7.47292442e+00] #+ np.array([N/(2*np.pi),0,0])  #redisidue 68
particles_positions_i[26,:] = [2.75850463e+01,  5.93272085e+00  ,1.48321926e+01] #+ np.array([N/(2*np.pi),0,0])  #redisidue 78
particles_positions_i[27,:] = [1.30657099e+01 , 2.78998908e+00, -4.34300699e+01] #+ np.array([N/(2*np.pi),0,0])  #redisidue 92
particles_positions_i[28,:] = [2.13601392e+01  ,4.60498390e+00 ,-7.05142009e+00] #+ np.array([N/(2*np.pi),0,0])  #redisidue 95


#Cterminal domain
particles_positions_i[29,:] = [2.30919709e+01,  5.28793369e+00,  1.99486911e+01] # + np.array([N/(2*np.pi),0,0]) #redisidue 104
particles_positions_i[30,:] = [1.94300305e+01 , 4.72778301e+00 ,-2.73938883e+01] #+ np.array([N/(2*np.pi),0,0])  #redisidue 121
particles_positions_i[31,:] = [1.83496689e-02, -1.72705028e+00, -2.25774036e+01]  #+ np.array([N/(2*np.pi),0,0]) #redisidue 137
particles_positions_i[32,:] = [1.91956718e+01 ,-1.48127658e+01 ,-2.52551277e+01] #+ np.array([N/(2*np.pi),0,0])  #redisidue 158
particles_positions_i[33,:] = [1.39387671e+01, -4.66782318e+00, -1.82640684e+01] #+ np.array([N/(2*np.pi),0,0])  #redisidue 181
particles_positions_i[34,:] = [1.76120982e+01 , 1.29030390e+01 ,-1.05347836e+01] # + np.array([N/(2*np.pi),0,0]) #redisidue 221
particles_positions_i[35,:] = [4.18911886e+01, -2.40967732e+00, -1.98593693e+01] #+ np.array([N/(2*np.pi),0,0])  #redisidue 240

#k_spring = 10 #strength of spring interaction
k_attra = 10 #strength of attractive interaction
k_rep = 10 #strength of repulsive interaction

li = [5,4,3,17,3,19,16,10,14,3,9]

#Define interactions in the case of one single protein per super cell

# Define interactions inside a protein

#triangle configuration :
interaction_unit_core_i = np.zeros((52,2),dtype = int)
interaction_unit_core_i[:,0] = [0,1,2,3,4,5,6,7,8,9,10,  11,11,11,11,11,11,  12,12,12,12,12,   13,13,13,13,   18,19,20,21,22,23,24,25,26,27,28,   29,29,29,29,29,29,   30,30,30,30,30,   31,31,31,31] #particles that are in interactions with ...
interaction_unit_core_i[:,1] = [1,2,3,4,5,6,7,8,9,10,11, 12,13,14,15,16,17,  13,14,15,16,17,   14,15,16,17,   19,20,21,22,23,24,25,26,27,28,29,   30,31,32,33,34,35,   31,32,33,34,35,   32,33,34,35] #these particles
interaction_unit_core_length_i = li #corresponds to zero energy length of vertices when only minimizing a single protein
l_N = 95 #length of the N-terminal in terms of aa
k_N = 0.1 #1/np.sqrt(l_N)
k_C = 10 #rqndom, change later
interaction_unit_core_strength_i = [k_N,k_N,k_N,k_N,k_N,k_N,k_N,k_N,k_N,k_N,k_N,k_N,] #corresponds to the spring strength in each vertices

#computing the distances in the undefrmer particle and setting the; as d0 for the core interactions in the C-terminal domaim

for i in range(11,14) :
    for j in range(i+1,18) :
        d = np.linalg.norm(particles_positions_i[i,:]-particles_positions_i[j,:])
        interaction_unit_core_length_i.append(d)
        interaction_unit_core_strength_i.append(k_C)

for i in range(11) :
    interaction_unit_core_length_i.append(li[i])
    interaction_unit_core_strength_i.append(k_N)

for i in range(29,32) :
    for j in range(i+1,36) :
        d = np.linalg.norm(particles_positions_i[i,:]-particles_positions_i[j,:])
        interaction_unit_core_length_i.append(d)
        interaction_unit_core_strength_i.append(k_C)


#attractive interactions
interaction_unit_attra_i = np.zeros((15,2))
interaction_unit_attra_i[:,0] = [0,0, 1, 2,2, 3,3, 4, 5, 6, 8,8,8, 10, 11] #,7]
interaction_unit_attra_i[:,1] = [29,35, 31, 30,31, 24,31,24, 31, 27, 28,30,35, 30, 32]#[11,27, 15, 13,15, 6,15, 6, 15, 9, 10,13,27, 13, 18] #,8]
interaction_unit_attra_length_i = [12,12, 4, 4,4, 4,4, 4, 4, 4, 12,12,12, 12, 12] #,3])
k_attra = 10
interaction_unit_attra_strength_i = [k_attra,k_attra,k_attra,k_attra,k_attra,k_attra,k_attra,k_attra,k_attra,k_attra,k_attra,k_attra,k_attra,k_attra,k_attra] #NOT USED, Just to give every interaction the same structure for the super_interaction function>

#repulsive interactions
interaction_unit_rep_i = np.zeros((1,2))
interaction_unit_rep_i[:,0] = [2] #,9]
interaction_unit_rep_i[:,1] = [33]#[21] #,6]
interaction_unit_rep_length_i = [4] #,0.9])
k_rep = 10
interaction_unit_rep_strength_i = [k_rep] #NOT USED

#Define interactions between proteins

#attractive interactions
nb_inter_attra = 10
interaction_inter_attra_i = np.zeros((nb_inter_attra,2),dtype=int)
interaction_inter_attra_i[:,0] = [7,7, 8, 10, 16, 25,25,26,28,34] #45,45,46,48,63,] #particles from the initial unit interactig whth ...
interaction_inter_attra_i[:,1] = [7,8, 8, 10, 16,25,26,26,28,34] #7,8, 8, 10, 25,] # ... particles from other units (not necessarily on the same other unit)
interaction_inter_attra_length_i = [4,4,4,4,4,4,4,4,4,4,]
k_attra = 10
interaction_inter_attra_strength_i = [k_attra,k_attra,k_attra,k_attra,k_attra,k_attra,k_attra,k_attra,k_attra,k_attra] #NOT USED

#repulsive interactions
nb_inter_rep = 2
interaction_inter_rep_i = np.zeros((nb_inter_rep,2),dtype=int)
interaction_inter_rep_i[:,0] = [1,19]#39} #,11] #particles from the initial unit interactig whth ...
interaction_inter_rep_i[:,1] = [1,19]#] #,11] # ... particles from other units (not necessarily on the same other unit)
interaction_inter_rep_length_i = [4,4]#,4] #in A
k_rep = 10
interaction_inter_rep_strength_i = [k_rep,k_rep] #NOT USED

#global repulsive interactions
nb_glob_rep = 182
interaction_glob_rep_i = np.zeros((nb_glob_rep,2))
interaction_glob_rep_i[:,0] = [11,11,11,11,11,11,11,11,11,11,11,11,11,   12,12,12,12,12,12,12,12,12,12,12,12,12,  13,13,13,13,13,13,13,13,13,13,13,13,13,   14,14,14,14,14,14,14,14,14,14,14,14,14,   15,15,15,15,15,15,15,15,15,15,15,15,15,   16,16,16,16,16,16,16,16,16,16,16,16,16,   17,17,17,17,17,17,17,17,17,17,17,17,17,   29,29,29,29,29,29,29,29,29,29,29,29,29,   30,30,30,30,30,30,30,30,30,30,30,30,30,   31,31,31,31,31,31,31,31,31,31,31,31,31,   32,32,32,32,32,32,32,32,32,32,32,32,32,   33,33,33,33,33,33,33,33,33,33,33,33,33,   34,34,34,34,34,34,34,34,34,34,34,34,34,   35,35,35,35,35,35,35,35,35,35,35,35,35]
interaction_glob_rep_i[:,1] = [12,13,14,15,16,17,29,30,31,32,33,34,35,   11,13,14,15,16,17,29,30,31,32,33,34,35,  11,12,14,15,16,17,29,30,31,32,33,34,35,   11,12,13,15,16,17,29,30,31,32,33,34,35,   11,12,13,14,16,17,29,30,31,32,33,34,35,   11,12,13,14,15,17,29,30,31,32,33,34,35,   11,12,13,14,15,16,29,30,31,32,33,34,35,   11,12,13,14,15,16,17,30,31,32,33,34,35,   11,12,13,14,15,16,17,29,31,32,33,34,35,   11,12,13,14,15,16,17,29,30,32,33,34,35,   11,12,13,14,15,16,17,29,30,31,33,34,35,   11,12,13,14,15,16,17,20,30,31,32,34,35,   11,12,13,14,15,16,17,29,30,31,32,33,35,   11,12,13,14,15,16,17,29,30,31,32,33,34]
g = 20
interaction_glob_rep_length_i = [g,g,g,g,g,g,g,g,g,g,g,g,g,   g,g,g,g,g,g,g,g,g,g,g,g,g,   g,g,g,g,g,g,g,g,g,g,g,g,g,   g,g,g,g,g,g,g,g,g,g,g,g,g,   g,g,g,g,g,g,g,g,g,g,g,g,g,   g,g,g,g,g,g,g,g,g,g,g,g,g,   g,g,g,g,g,g,g,g,g,g,g,g,g,   g,g,g,g,g,g,g,g,g,g,g,g,g,   g,g,g,g,g,g,g,g,g,g,g,g,g,   g,g,g,g,g,g,g,g,g,g,g,g,g,   g,g,g,g,g,g,g,g,g,g,g,g,g,   g,g,g,g,g,g,g,g,g,g,g,g,g,   g,g,g,g,g,g,g,g,g,g,g,g,g,   g,g,g,g,g,g,g,g,g,g,g,g,g]
k_glob = 10
interaction_glob_rep_strength_i = [10,10,10,10,10,10,10,10,10,10,10,10,10,   10,10,10,10,10,10,10,10,10,10,10,10,10,   10,10,10,10,10,10,10,10,10,10,10,10,10,   10,10,10,10,10,10,10,10,10,10,10,10,10,   10,10,10,10,10,10,10,10,10,10,10,10,10,   10,10,10,10,10,10,10,10,10,10,10,10,10,   10,10,10,10,10,10,10,10,10,10,10,10,10,   10,10,10,10,10,10,10,10,10,10,10,10,10,   10,10,10,10,10,10,10,10,10,10,10,10,10,   10,10,10,10,10,10,10,10,10,10,10,10,10,   10,10,10,10,10,10,10,10,10,10,10,10,10,   10,10,10,10,10,10,10,10,10,10,10,10,10,   10,10,10,10,10,10,10,10,10,10,10,10,10,   10,10,10,10,10,10,10,10,10,10,10,10,10]


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
