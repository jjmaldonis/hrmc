import sys, os
import shlex, subprocess, shutil

def main():
    ts = 24000000 # total number of steps in the simulation
    js = 400000 # number of steps per job
    ss = 0 # starting step, this number gets incremented by js on each submit
    temp = 50 # starting temperature
    mm = 1.5 # starting maxmove

    for i in range(0,ts/js):
        if(i > 0):
            src = os.getcwd() + '/param_file.in'
            dst = os.getcwd() + '/param_file.{0}.in'.format(jobid)
            shutil.copyfile(src,dst)
            #shutil.copyfile('paramfile.in','paramfile.{0}.in'.format(jobid))
            args = "sed -i '5s/.*/{0}   !startingstep/' param_file.{1}.in".format(ss,jobid)
        else:
            args = "sed -i '5s/.*/{0}   !startingstep/' param_file.in".format(ss)
        args = shlex.split(args)
        print(args)
        p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        poutput = p.stdout.read()
        perr = p.stderr.read()
        print(poutput,perr)
        preturncode = p.wait()

        if(preturncode != 0):
            raise Exception("sed failed!")
        if(i > 0):
            args = "sbatch --dependency=afterok:{0} submits/slurm_noomp.sh param_file.{0}.in".format(jobid)
            #args = "qsub --dependency=afterok:{0} submits/slurm_noomp.sh param_file.{0}.in".format(jobid)
        else:
            args = "sbatch submits/slurm_noomp.sh param_file.in"
            #args = "qsub submits/slurm_noomp.sh param_file.in"
        args = shlex.split(args)
        print(args)
        p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        pcomm = p.communicate()

        print('')
        print(pcomm[0].strip())
        print(pcomm[1].strip())
        #jobid = int(pcomm[0].strip().split('\n').pop().split()[2])
        #print(jobid)
        print('')

        preturncode = p.wait()
        if(preturncode != 0):
            raise Exception("sbatch failed!")
        jobid = int(pcomm[0].strip().split('\n').pop().split()[-1])
        ss = ss + js


if __name__ == '__main__':
    main()
