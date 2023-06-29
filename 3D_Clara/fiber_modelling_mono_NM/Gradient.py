# Importing packages and files
import numpy as np
from Energy import *
# from initial_configuration_mono_to_di import *
from initial_configuration_di_to_fi import *
from Neighbors import *



def derivate_parameters(parameters,unit_index,n,m) :
    i,j = number_to_coordinates(nb_units, unit_index)
    a1 = parameters[0]
    b1 = parameters[1]
    b2 = parameters[2]
    N = np.sqrt((m*b2)**2 + (n*a1 + m*b1)**2)

    # dpara = np.zeros((2,3))
    # dpara[0,0] = 2*np.pi*((n*(b1*j + i*a1) + i*(a1*n + b1*m))/((a1*n + b1*m)**2 + b2**2*m**2) - (2*n*(a1*n + b1*m)*(b2**2*j*m + (b1*j + i*a1)*(a1*n + b1*m)))/((a1*n + b1*m)**2 + b2**2*m**2)**2) #dphi/da1
    # dpara[0,1] = 2*np.pi*((j*(a1*n + b1*m) + m*(b1*j + i*a1))/((a1*n + b1*m)**2 + b2**2*m**2) - (2*m*(a1*n + b1*m)*(b2**2*j*m + (b1*j + i*a1)*(a1*n + b1*m)))/((a1*n + b1*m)**2 + b2**2*m**2)**2) #dphi/db1
    # dpara[0,2] = 2*np.pi*(2*a1*b2*m*(j*n - i*m)*(a1*n + b1*m))/((a1*n + b1*m)**2 + b2**2*m**2)**2 #dphi/db2
    # dpara[1,0] = (b2*m*(j*n - i*m)*(a1*b1*n + m*(b1**2 + b2**2)))/((a1*n + b1*m)**2 + b2**2*m**2)**(3/2) #de/da1
    # dpara[1,1] = (-(n*b2*j*a1-m*b2*i*a1)*(m*(n*a1+m*b1))/(N))/(N**2)#de/db1
    # dpara[1,2] = (a1*(j*n - i*m)*(a1*n + b1*m)**2)/((a1*n + b1*m)**2 + b2**2*m**2)**(3/2) #de/db2

    # dpara = np.zeros((2,3))
    # dpara[0,0] = 2*np.pi*(i*(m + i*j*n)*(2*a1*m*n*(b1**2 + b2**2) + a1**2*b1*n**2 + b1*m**2*(b1**2 + b2**2)))/(2*a1*b1*m*n + a1**2*n**2 + m**2*(b1**2 + b2**2))**2 #dphi/da1
    # dpara[0,1] = 2*np.pi*(a1*(j*n - i*m)*(a1**2*n**2 + 2*a1*b1*m*n + b1**2*m**2 - b2**2*m**2))/(a1**2*n**2 + 2*a1*b1*m*n + b1**2*m**2 + b2**2*m**2)**2 #dphi/db1
    # dpara[0,2] = 2*np.pi*(2*a1*b2*m*(j*n - i*m)*(a1*n + b1*m))/((a1*n + b1*m)**2 + b2**2*m**2)**2 #dphi/db2
    # dpara[1,0] = (b2*m*(j*n - i*m)*(a1*b1*n + m*(b1**2 + b2**2)))/((a1*n + b1*m)**2 + b2**2*m**2)**(3/2) #de/da1
    # dpara[1,1] = (i*a1*b2*m*(m + i*j*n)*(a1*n + b1*m))/((a1*n + b1*m)**2 + b2**2*m**2)**(3/2)#de/db1
    # dpara[1,2] = (a1*(j*n - i*m)*(a1*n + b1*m)**2)/((a1*n + b1*m)**2 + b2**2*m**2)**(3/2) #de/db2

    dpara = np.zeros((3,3))
    dpara[0,0] = (n*(a1*n + b1*m))/(2*np.pi*np.sqrt((a1*n + b1*m)**2 + b2**2*m**2)) #dr/da1 
    dpara[0,1] = (m*(a1*n + b1*m))/(2*np.pi*np.sqrt((a1*n + b1*m)**2 + b2**2*m**2))  #dr/db1
    dpara[0,2] = (b2*m**2)/(2*np.pi*np.sqrt((a1*n + b1*m)**2 + b2**2*m**2)) #dr/db2
    dpara[1,0] = 2*np.pi*((n*(b1*j + i*a1) + i*(a1*n + b1*m))/((a1*n + b1*m)**2 + b2**2*m**2) - (2*n*(a1*n + b1*m)*(b2**2*j*m + (b1*j + i*a1)*(a1*n + b1*m)))/((a1*n + b1*m)**2 + b2**2*m**2)**2) #dtheta/da1
    dpara[1,1] = 2*np.pi*((j*(a1*n + b1*m) + m*(b1*j + i*a1))/((a1*n + b1*m)**2 + b2**2*m**2) - (2*m*(a1*n + b1*m)*(b2**2*j*m + (b1*j + i*a1)*(a1*n + b1*m)))/((a1*n + b1*m)**2 + b2**2*m**2)**2) #dtheta/db1
    dpara[1,2] = 2*np.pi*(2*a1*b2*m*(j*n - i*m)*(a1*n + b1*m))/((a1*n + b1*m)**2 + b2**2*m**2)**2 #dtheta/db2
    dpara[2,0] = (b2*m*(j*n - i*m)*(a1*b1*n + m*(b1**2 + b2**2)))/((a1*n + b1*m)**2 + b2**2*m**2)**(3/2) #dz/da1
    dpara[2,1] = (-(n*b2*j*a1-m*b2*i*a1)*(m*(n*a1+m*b1))/(N))/(N**2)#dz/db1
    dpara[2,2] = (a1*(j*n - i*m)*(a1*n + b1*m)**2)/((a1*n + b1*m)**2 + b2**2*m**2)**(3/2) #dz/db2

    return(dpara)

def parameters_derivative(parameters,unit_index,n,m,c1,c2,param) :
    i,j = number_to_coordinates(nb_units, unit_index)
    x1 = c1[0]
    y1 = c1[1]
    z1 = c1[2]
    x2 = c2[0]
    y2 = c2[1]
    z2 = c2[2]
    phi,e = cylinder_parameters(parameters,n,m,i,j)
    dpara = derivate_parameters(parameters,unit_index,n,m)
    dphi = dpara[0,param]
    de = dpara[1,param]
    #d = np.linalg.norm(c1-cNN)
    #f = k*(d-d0)/d*{c1-cNN}
    der = np.array([-np.sin(phi)*x2*dphi - np.cos(phi)*y2*dphi,-np.sin(phi)*y2*dphi + np.cos(phi)*x2*dphi,z2*de])
    return(der)


def Gradient_Numerical(minimisation_parameters,*args) :
    (interprot_indexes) = args
    l = np.size(minimisation_parameters)
    grad_num = np.zeros(l)
    E0 = total_energy(minimisation_parameters,*args)
    for i in range(l) :
        p1 = np.zeros(l)
        p2 = np.zeros(l)
        for j in range(l) :
            p1[j] = minimisation_parameters[j]
            p2[j] = minimisation_parameters[j]
        p1[i] = p1[i] + 0.0000001/2
        p2[i] = p2[i] - 0.0000001/2
        E1 = total_energy(p1,*args)
        E2 = total_energy(p2,*args)
        grad_num[i] = (E1 - E2)/0.0000001
    #print(grad_num)
    return(grad_num)
    

def Gradient(minimisation_parameters,*args) :  
    (interprot_indexes,nb_copies,n,m) = args
    #Defining the interactions in the super cells
    interaction_unit_core,interaction_unit_core_length,interaction_unit_core_strength = super_interaction(nb_copies,interaction_unit_core_i,interaction_unit_core_length_i,interaction_unit_core_strength_i)
    interaction_unit_attra,interaction_unit_attra_length,interaction_unit_attra_strength = super_interaction(nb_copies,interaction_unit_attra_i,interaction_unit_attra_length_i,interaction_unit_attra_strength_i)
    interaction_unit_rep,interaction_unit_rep_length,interaction_unit_rep_strength = super_interaction(nb_copies,interaction_unit_rep_i,interaction_unit_rep_length_i,interaction_unit_rep_strength_i)
    interaction_inter_attra,interaction_inter_attra_length,interaction_inter_attra_strength = super_interaction(nb_copies,interaction_inter_attra_i,interaction_inter_attra_length_i,interaction_inter_attra_strength_i)
    interaction_inter_rep,interaction_inter_rep_length,interaction_inter_rep_strength = super_interaction(nb_copies,interaction_inter_rep_i,interaction_inter_rep_length_i,interaction_inter_rep_strength_i)
    interaction_glob_rep,interaction_glob_rep_length,interaction_glob_rep_strength = super_interaction(nb_copies,interaction_glob_rep_i,interaction_glob_rep_length_i,interaction_glob_rep_strength_i)
    #print(interaction_unit_attra)
    #Getting particles positions and parameters from minimisation parameters
    particles_positions = positions_from_mini_para(minimisation_parameters,nb_copies,n,m)
    parameters = minimisation_parameters[-3:]
    interacting = np.unique(interprot_indexes)
    list_proteins = [particles_positions]
    for k in interacting :
        list_proteins.append(find_protein_x(parameters,particles_positions,k,n,m))
    #print(parameters)
    total_force = np.zeros((nb_particles*nb_copies,3))
    grad_para = np.zeros((3))
    #intraprotein interactions
    # unit core interactions
    l = np.size(interaction_unit_core_length)
    #print(l)
    for i in range(l) :
        #print(i)
        p1 = int(interaction_unit_core[i,0])
        #print(p1)
        p2 = int(interaction_unit_core[i,1])
        #print(p2)
        c1 = particles_positions[p1,:]
        c2 = particles_positions[p2,:]
        d0 = interaction_unit_core_length[i]
        k_spring = interaction_unit_core_strength[i]
        total_force[p1,:] = total_force[p1,:] + force_spring(c1,c2,d0,k_spring)
        total_force[p2,:] = total_force[p2,:] - force_spring(c1,c2,d0,k_spring) #*np.array([-1,1,1])
        #print('f',total_force)
    #unit attractive interactions
    l = np.size(interaction_unit_attra_length)
    for i in range(l) :
        p1 = int(interaction_unit_attra[i,0])
        p2 = int(interaction_unit_attra[i,1])
        c1 = particles_positions[p1,:]
        c2 = particles_positions[p2,:]
        d0 = interaction_unit_attra_length[i]
        k_attra = interaction_unit_attra_strength[i]
        total_force[p1,:] = total_force[p1,:] + force_attra(c1,c2,d0,k_attra)
        total_force[p2,:] = total_force[p2,:] - force_attra(c1,c2,d0,k_attra) #*np.array([-1,1,1])
        #print('f',force_attra(c1,c2,d0,1))
    #unit repulsive ineractions
    l = np.size(interaction_unit_rep_length)
    for i in range(l) :
        p1 = int(interaction_unit_rep[i,0])
        p2 = int(interaction_unit_rep[i,1])
        c1 = particles_positions[p1,:]
        c2 = particles_positions[p2,:]
        d0 = interaction_unit_rep_length[i]
        k_rep = interaction_unit_rep_strength[i]
        total_force[p1,:] = total_force[p1,:] + force_rep(c1,c2,d0,k_rep)
        total_force[p2,:] = total_force[p2,:] - force_rep(c1,c2,d0,k_rep) #*np.array([-1,1,1])
        #print('f',total_force)
    #interprotein interactions
    count = 0 #to know which interactions have already been treated.
    #attractive interactions
    #print('inter_attra')
    l = np.size(interaction_inter_attra_length)
    #print(l)
    for i in range(l) :
        p1 = int(interaction_inter_attra[i,0])
        p2 = int(interaction_inter_attra[i,1])
        #print(p1,p2)
        c1 = particles_positions[p1,:]
        c2 = particles_positions[p2,:]
        unit_index = interprot_indexes[count] #- nb_units//2
        cNN = list_proteins[unit_index][p2,:] #coordinates of te second particle in the interacting unit
        dpara = derivate_parameters(parameters,unit_index,n,m)
        d0 = interaction_inter_attra_length[i]
        k_attra = interaction_inter_attra_strength[i]
        phi,e = cylinder_parameters(parameters,n,m,i,j)
        force = force_attra(c1,cNN,d0,k_attra)
        #print('force',force_attra(c1,cNN,d0,k_attra))
        if p1==p2 :
            total_force[p1,:] = total_force[p1,:] + (1/2)*force*(1-np.cos(phi)) 
            total_force[p2,:] = total_force[p2,:] - (1/2)*force*(np.cos(phi)-1)
        else :
            total_force[p1,:] = total_force[p1,:] + (1/2)*force 
            total_force[p2,:] = total_force[p2,:] - (1/2)*force
        count += 1
        dera1 = parameters_derivative(parameters,unit_index,n,m,c1,c2,0)
        derb1 = parameters_derivative(parameters,unit_index,n,m,c1,c2,1)
        derb2 = parameters_derivative(parameters,unit_index,n,m,c1,c2,2)
        if p1 == 0 and p2 == 0 :
            grad_para[0] = grad_para[0] -1/2*(force[0]*(dpara[0,0] + np.sin(phi)*c2[0]*dpara[1,0] - np.cos(phi)*dpara[0,0]) + force[1]*(0 - np.cos(phi)*dpara[1,0]*c2[0] - np.sin(phi)*dpara[0,0]) + force_attra(c1,cNN,d0,k_attra)[2]*(-dpara[2,0]))
            grad_para[1] = grad_para[1] -1/2*(force[0]*(dpara[0,1] + np.sin(phi)*c2[0]*dpara[1,1] - np.cos(phi)*dpara[0,1]) + force[1]*(0 - np.cos(phi)*dpara[1,1]*c2[0] - np.sin(phi)*dpara[0,1]) + force_attra(c1,cNN,d0,k_attra)[2]*(-dpara[2,1]))
            grad_para[2] = grad_para[2] -1/2*(force[0]*(dpara[0,2] + np.sin(phi)*c2[0]*dpara[1,2] - np.cos(phi)*dpara[0,2]) + force[1]*(0 - np.cos(phi)*dpara[1,2]*c2[0] - np.sin(phi)*dpara[0,2]) + force_attra(c1,cNN,d0,k_attra)[2]*(-dpara[2,2]))
        elif p1 == 0 and p2 != 0 :
            grad_para[0] = grad_para[0] -1/2*(force[0]*(dpara[0,0] + np.sin(phi)*c2[0]*dpara[1,0] + np.cos(phi)*c2[1]*dpara[1,0]) + force[1]*(0 - np.cos(phi)*dpara[1,0]*c2[0] + np.sin(phi)*dpara[1,0]*c2[1]) + force[2]*(-dpara[2,0]))
            grad_para[1] = grad_para[1] -1/2*(force[0]*(dpara[0,1] + np.sin(phi)*c2[0]*dpara[1,1] + np.cos(phi)*c2[1]*dpara[1,1]) + force[1]*(0 - np.cos(phi)*dpara[1,1]*c2[0] + np.sin(phi)*dpara[1,1]*c2[1]) + force[2]*(-dpara[2,0]))
            grad_para[2] = grad_para[2] -1/2*(force[0]*(dpara[0,2] + np.sin(phi)*c2[0]*dpara[1,2] + np.cos(phi)*c2[1]*dpara[1,2]) + force[1]*(0 - np.cos(phi)*dpara[1,2]*c2[0] + np.sin(phi)*dpara[1,2]*c2[1]) + force[2]*(-dpara[2,0]))
        elif p1 != 0 and p2 == 0 :
            grad_para[0] = grad_para[0] -1/2*(force[0]*(0 + np.sin(phi)*c2[0]*dpara[1,0] - np.cos(phi*dpara[0,0])) + force[1]*(0 - np.cos(phi)*dpara[1,0]*c2[0] - np.sin(phi)*dpara[0,0]) + force[2]*(-dpara[2,0]))
            grad_para[1] = grad_para[1] -1/2*(force[0]*(0 + np.sin(phi)*c2[0]*dpara[1,1] - np.cos(phi*dpara[0,1])) + force[1]*(0 - np.cos(phi)*dpara[1,1]*c2[0] - np.sin(phi)*dpara[0,1]) + force[2]*(-dpara[2,1]))
            grad_para[2] = grad_para[2] -1/2*(force[0]*(0 + np.sin(phi)*c2[0]*dpara[1,2] - np.cos(phi*dpara[0,2])) + force[1]*(0 - np.cos(phi)*dpara[1,2]*c2[0] - np.sin(phi)*dpara[0,2]) + force[2]*(-dpara[2,2]))
        elif p1 != 0 and p2 !=0 :
            grad_para[0] = grad_para[0] -1/2*(force[0]*(0 + np.sin(phi)*c2[0]*dpara[1,0] + np.cos(phi)*c2[1]*dpara[1,0]) + force[1]*(0 - np.cos(phi)*dpara[1,0]*c2[0] + np.sin(phi)*dpara[1,0]*c2[1]) + force[2]*(-dpara[2,0]))
            grad_para[1] = grad_para[1] -1/2*(force[0]*(0 + np.sin(phi)*c2[0]*dpara[1,1] + np.cos(phi)*c2[1]*dpara[1,1]) + force[1]*(0 - np.cos(phi)*dpara[1,1]*c2[0] + np.sin(phi)*dpara[1,1]*c2[1]) + force[2]*(-dpara[2,1]))
            grad_para[2] = grad_para[2] -1/2*(force[0]*(0 + np.sin(phi)*c2[0]*dpara[1,2] + np.cos(phi)*c2[1]*dpara[1,2]) + force[1]*(0 - np.cos(phi)*dpara[1,2]*c2[0] + np.sin(phi)*dpara[1,2]*c2[1]) + force[2]*(-dpara[2,2]))
            #print('force' , force_attra(c1,cNN,d0,k_attra)[0]*(0 + np.sin(phi)*c2[0]*dpara[1,0] + np.cos(phi)*c2[1]*dpara[1,0]) + force_attra(c1,cNN,d0,k_attra)[1]*(0 - np.cos(phi)*dpara[1,0]*c2[0] + np.sin(phi)*dpara[1,0]*c2[1]) + force_attra(c1,cNN,d0,k_attra)[2]*(-dpara[2,0]))
        #print(grad_para[0])

        #grad_para = grad_para - 1/2*np.array([np.dot(force_attra(c1,cNN,d0,k_attra),dera1),np.dot(force_attra(c1,cNN,d0,k_attra),derb1),np.dot(force_attra(c1,cNN,d0,k_attra),derb2)])
#     # #repulsive interactions
#     #count = 5
#     l = np.size(interaction_inter_rep_length)
#     #print(l)
#     for i in range(l) :
#         p1 = int(interaction_inter_rep[i,0])
#         p2 = int(interaction_inter_rep[i,1])
#         c1 = particles_positions[p1,:]
#         c2 = particles_positions[p2,:]
#         unit_index = interprot_indexes[count] #- nb_units//2 
#         #print(unit_index)
#         #print(number_to_coordinates(nb_units,unit_index))
#         cNN = list_proteins[unit_index][p2,:] #coordinates of the second particle in the interacting unit
#         dpara = derivate_parameters(parameters,unit_index,n,m)
#         #print(dpara)
#         #print(dpara)
#         d0 = interaction_inter_rep_length[i]
#         k_rep = interaction_inter_rep_strength[i]
#         total_force[p1,:] = total_force[p1,:] + (1/2)*force_rep(c1,cNN,d0,k_rep)
#         total_force[p2,:] = total_force[p2,:] - (1/2)*force_rep(c1,cNN,d0,k_rep) #*np.array([-1,1,1])
#         count += 1
#         dera1 = parameters_derivative(parameters,unit_index,n,m,c1,c2,0)
#         derb1 = parameters_derivative(parameters,unit_index,n,m,c1,c2,1)
#         derb2 = parameters_derivative(parameters,unit_index,n,m,c1,c2,2)
#         grad_para = grad_para - 1/2*np.array([np.dot(force_rep(c1,cNN,d0,k_rep),dera1),np.dot(force_rep(c1,cNN,d0,k_rep),derb1),np.dot(force_rep(c1,cNN,d0,k_rep),derb2)])
# # GLOBAL repulsion
#     interacting = np.unique(interprot_indexes)
#     #print(interacting)
#     for k in interacting :
#         #print(k)
#         l = np.size(interaction_glob_rep_length)
#         #print(l)
#         for i in range(l) :
#             p1 = int(interaction_glob_rep[i,0])
#             p2 = int(interaction_glob_rep[i,1])
#             c1 = particles_positions[p1,:]
#             c2 = particles_positions[p2,:]
#             unit_index = k
#             #print(unit_index)
#             #print(number_to_coordinates(nb_units,unit_index))
#             cNN = list_proteins[unit_index][p2,:] #coordinates of the second particle in the interacting unit
#             #print(dpara)
#             #print(dpara)
#             d0 = interaction_glob_rep_length[i]
#             k_glob= interaction_glob_rep_strength[i]
#             total_force[p1,:] = total_force[p1,:] + (1/2)*force_rep(c1,cNN,d0,k_glob)
#             total_force[p2,:] = total_force[p2,:] - (1/2)*force_rep(c1,cNN,d0,k_glob) #*np.array([-1,1,1])
#             dera1 = parameters_derivative(parameters,unit_index,n,m,c1,c2,0)
#             derb1 = parameters_derivative(parameters,unit_index,n,m,c1,c2,1)
#             derb2 = parameters_derivative(parameters,unit_index,n,m,c1,c2,2)
#             grad_para = grad_para - 1/2*np.array([np.dot(force_rep(c1,cNN,d0,k_glob),dera1),np.dot(force_rep(c1,cNN,d0,k_glob),derb1),np.dot(force_rep(c1,cNN,d0,k_glob),derb2)])
#     #creating a 1D gradient :
    grad1D = np.reshape(total_force[1:,:],((nb_particles-1)*nb_copies*3)) #remove the derivatives with respect to coordinates of particle 0
    total_grad1D = np.concatenate((grad1D ,grad_para))
    #print('grad',total_grad1D)
    return(total_grad1D)