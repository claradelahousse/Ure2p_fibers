import subprocess
import numpy as np





"""Choose the queue"""
# myqueue = 'q-1jour'
#myqueue = 'q-1mois'
myqueue = 'q-1sem'
# myqueue = 'q-2sem'



"""Info on the queues"""
# q-1jour mem<3GB
# q-1sem mem<3.3GB
# q-2sem mem<3GB
# q-1mois mem<15GB
choose_type = 'pbs' #'sbatch', 'pbs'

"""Choose name of the job"""
name = 'fiber_mono'


"""Choose the code you want to run"""
directory= '/home/cdelahousse/code/v2_fiber_modelling_dim_NM_g'
codeToExecute = 'Main.py'
# codeToExecute = 'FindLEL.py'
# codeToExecute = 'ComputeFlow.py'
# codeToExecute = 'PatchyHexagons.py'


"""Write the instructions in a pbs file that will be read by titan"""
extension='.pbs'
f = open("launcher.pbs","w")
f.write("#!/bin/bash\n")
f.write("#PBS -V\n")
f.write("#PBS -q " + myqueue + " \n")
 #f.write("#PBS -l cput=100:00:00\n")
 #f.write("#PBS -l mem=2GB\n")
f.write("#PBS -N " + name + " \n")
f.write("#PBS -m ae \n")
f.write("#PBS -M clara.delahousse@ens-paris-saclay.fr\n")
f.write("cd "+directory+" \n")

# f.write("make clean")
# f.write("make")
f.write("python3 -u "+codeToExecute+"\n") 

f.close()



"""Launch the code"""
terminal_input="sbatch launcher"+extension
proc = subprocess.Popen([terminal_input], stdout=subprocess.PIPE,shell=True)
(out, err) = proc.communicate()









