import sys, os, shutil
sys.path.append('/home/jjmaldonis/model_analysis/scripts')
from model import Model
from subproc import run_subproc

def collect(m,modelfile):
    content = open('/home/jjmaldonis/development/mostupdated_hrmc_for_dev/hrmc/param_file.in').readlines()
    content[1] = modelfile+'\n'
    open('paramfile.temp.in','w').writelines(content)
    path = os.getcwd()
    if( not os.path.isfile('ZrCuAl2011.eam.alloy') ):
        copy = True
        shutil.copyfile('/home/jjmaldonis/development/mostupdated_hrmc_for_dev/hrmc/ZrCuAl2011.eam.alloy',path+'/ZrCuAl2011.eam.alloy')
    else:
        copy = False
    output = run_subproc('/home/jjmaldonis/bin/energy_per_atom {0}/paramfile.temp.in'.format(path))
    os.remove('paramfile.temp.in')
    if(copy):
        os.remove(path+'/ZrCuAl2011.eam.alloy')
    output = output.strip().split('\n')
    while('Avg energy per atom' not in output[0]):
        output.pop(0)
    output.pop(0)
    for line in output:
        line = line.strip().split()
        at = int(line[0]) - 1
        energy = float(line[1])
        m.atoms[at].energy = energy
    return m

def main():
    modelfile = sys.argv[1]
    m = Model(modelfile)
    collect(m,modelfile)


if __name__ == "__main__":
    main()
