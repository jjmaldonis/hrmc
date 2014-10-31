import sys, os, shutil
sys.path.append('/home/jjmaldonis/model_analysis/scripts')
from model import Model
from subproc import run_subproc

def collect(m,modelfile):
    # Get some paramfile with the correct format
    content = open('/home/jjmaldonis/development/mostupdated_hrmc_for_dev/hrmc/param_file.in').readlines()
    content[1] = modelfile+'\n'
    # After changing the modelfile to the one that was given, save a temporary paramfile
    open('paramfile.temp.in','w').writelines(content)
    path = os.getcwd()
    # Get the eam file if it isnt here already
    if( not os.path.isfile('ZrCuAl2011.eam.alloy') ):
        copy = True
        shutil.copyfile('/home/jjmaldonis/development/mostupdated_hrmc_for_dev/hrmc/ZrCuAl2011.eam.alloy',path+'/ZrCuAl2011.eam.alloy')
    else:
        copy = False
    # Run the energy_per_atom code and keep the output
    output = run_subproc('/home/jjmaldonis/bin/energy_per_atom {0}/paramfile.temp.in'.format(path))
    # Remove temporary paramfile and eam file if it wasnt here
    os.remove('paramfile.temp.in')
    if(copy):
        os.remove(path+'/ZrCuAl2011.eam.alloy')
    output = output.strip().split('\n')
    # Collect and save the energies per atom into atom.energy
    while('Avg energy per atom' not in output[0]):
        output.pop(0)
    output.pop(0)
    for line in output:
        line = line.strip().split()
        at = int(line[0]) - 1
        energy = float(line[1])
        m.atoms[at].energy = energy
    # Return the updated model, but it should be changed upon return anyways
    return m

def main():
    modelfile = sys.argv[1]
    m = Model(modelfile)
    collect(m,modelfile)


if __name__ == "__main__":
    main()
