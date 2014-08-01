import sys, os
import shlex, subprocess, shutil

def run_subproc(args):
    """ Run subprocess given args as a command line string"""
    print(args)
    args = shlex.split(args)
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    pcomm = p.communicate()
    poutput = pcomm[0] #p.stdout.read()
    perr = pcomm[1] #p.stderr.read()
    print(poutput,perr)
    preturncode = p.wait()
    if(preturncode != 0):
        raise Exception("{0} failed!".format(args[0]))
    return pcomm

def submit_job(paramfile,prev_jobid,ss,modelfile=None):
    """ Pass None to prev_jobid if using param_file.in as paramfile """
    # Call sed to change starting step number to ss
    args = "sed -i '5s/.*/{0}   !startingstep/' {1}".format(ss,paramfile)
    pcomm = run_subproc(args)

    # Call sed to change modelfile if modelfile was given
    if(modelfile != None):
        args = "sed -i '2s/.*/model_final_{0}.txt/' {1}".format(jobid,paramfile)
        pcomm = run_subproc(args)

    # Do the submitting
    if(prev_jobid != None):
        args = "sbatch --dependency=afterok:{0} submits/slurm_noomp.sh {1}".format(jobid,paramfile)
    else:
        args = "sbatch submits/slurm_noomp.sh {0}".format(paramfile)
    pcomm = run_subproc(args)
    jobid = int(pcomm[0].strip().split('\n').pop().split()[-1])
    return jobid

def main():
    """ If you want to resume a job, all you have to do is change ss and starting_param_filename """
    ts = 24000000 # total number of steps in the simulation
    js = 400000 # number of steps per job
    ss = 0 # starting step, this number gets changed in the loop ~ incremented by js on each submit
    temp = 50 # starting temperature
    mm = 1.5 # starting maxmove
    starting_param_filename = 'param_file.in'

    for s in range(ss, ts, js):
        if(s > ss):
            # Copy the starting parameter file to another filename which includes the step being started
            src = os.getcwd() + '/' + starting_param_filename
            dst = os.getcwd() + '/param_file.{0}.in'.format(s)
            shutil.copyfile(src,dst)
            # Submit the job using that parameter file, but change it first
            jobid = submit_job(dst,jobid,ss,'model_final_{0}.txt'.format(jobid))
        else:
            # On the first run, just submit the starting parameter file
            dst = starting_param_filename
            jobid = submit_job(dst,None,ss)


if __name__ == '__main__':
    main()
