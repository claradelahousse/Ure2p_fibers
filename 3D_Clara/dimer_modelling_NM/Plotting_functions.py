######### First cell

#%matplotlib notebook #Don't forget to remove it

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sp
from tqdm import tqdm
from initial_configuration_mono_to_di import *
#from initial_configuration_di_to_fi import *
from Energy import *
from Gradient import *
from Neighbors import *

def surface_map(carte,r) :
    u,v = np.mgrid[0:2*np.pi:200j,0:np.pi:100j]
    X = carte[0] + r*np.sin(u)*np.cos(v)
    Y = carte[1] + r*np.sin(u)*np.sin(v)
    Z = carte[2] + r*np.cos(u)
    return(X,Y,Z)

def color_couples() :
    l = len(final_interprot_indexes)
    n = nb_inter_attra*(index_min+1) #index_min + 1 corresponds to the number of copies in a super cell in the minimal energy configuration
    r = nb_inter_rep*(index_min+1)
    attra_couples_list = [(0,k) for k in interaction_inter_attra[:,1]]
    rep_couples_list = [(0,k) for k in interaction_inter_rep[:,1]]
    for i in range(n) :
        attra_couples_list.append((final_interprot_indexes[i]- nb_units//2,interaction_inter_attra[i,1]))
        rep_couples_list.append((final_interprot_indexes[i]- nb_units//2,interaction_inter_rep[i,1]))
    return(attra_couples_list, rep_couples_list)


#defining ineraction color caracteristics.
c_all = 'mediumseagreen'
c_attra = 'hotpink'
c_rep = 'mediumpurple'

# Define the rotation function
def rotate_view(elev, azim,ax,fig):
    ax.view_init(elev=elev, azim=azim)
    fig.canvas.draw()

# # Set the initial view
# elev = 30
# azim = -60
# rotate_view(elev, azim,ax,fig)

# Define the key press function
def on_key_press(event):
    global elev, azim
    if event.key == 'up':
        elev += 10
    elif event.key == 'down':
        elev -= 10
    elif event.key == 'left':
        azim -= 10
    elif event.key == 'right':
        azim += 10
    rotate_view(elev, azim)

#fig.canvas.mpl_connect('key_press_event', on_key_press)

