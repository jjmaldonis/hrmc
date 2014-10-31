import sys, os, math
import shlex, subprocess, shutil

def run_subproc(args):
    """ Run subprocess given args as a command line string"""
    print(args)
    args = shlex.split(args)
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    pcomm = p.communicate()
    poutput = pcomm[0] #p.stdout.read()
    perr = pcomm[1] #p.stderr.read()
    print(poutput)
    print(perr)
    preturncode = p.wait()
    if(preturncode != 0):
        raise Exception("{0} failed!".format(args[0]))
    return pcomm

def submit_job(paramfile, prev_jobid, s, starting_temp, starting_maxmove, modelfile=None):
    """ Pass None to prev_jobid if using param_file.in as paramfile """
    # Call sed to change starting step number to s
    args = "sed -i '5s/.*/{0}   !startingstep/' {1}".format(s,paramfile)
    pcomm = run_subproc(args)
    # Call sed to change temp
    args = "sed -i '6s/.*/{0}   !starting temp/' {1}".format(starting_temp*math.sqrt(0.7)**(s/200000),paramfile)
    pcomm = run_subproc(args)
    # Call sed to change maxmove
    args = "sed -i '7s/.*/{0}   !starting maxmove/' {1}".format(starting_maxmove*math.sqrt(0.94)**(s/200000),paramfile)
    pcomm = run_subproc(args)

    # Call sed to change modelfile if modelfile was given
    if(modelfile != None):
        args = "sed -i '2s/.*/model_final_{0}.txt/' {1}".format(prev_jobid,paramfile)
        pcomm = run_subproc(args)

    # Do the submitting
    if(prev_jobid != None):
        args = "sbatch --dependency=afterok:{0} submits/stampede_submit.sh {1}".format(prev_jobid,paramfile)
    else:
        args = "sbatch submits/stampede_submit.sh {0}".format(paramfile)
    pcomm = run_subproc(args)
    jobid = int(pcomm[0].strip().split('\n').pop().split()[-1])
    return jobid

def main():
    """ If you want to resume a job, all you have to do is change ss and starting_param_filename """
    ts = 4000000 # total number of steps in the simulation
    js = 2000000 # number of steps per job
    ss = 0 # starting step, this number gets changed in the loop ~ incremented by js on each submit
    temp = 226 # starting temperature
    mm = 1.5 # starting maxmove
    starting_param_filename = 'param_file.in'

    for s in range(ss, ts, js):
        if(s > ss):
            # Copy the starting parameter file to another filename which includes the step being started
            src = os.getcwd() + '/' + starting_param_filename
            dst = os.getcwd() + '/param_file.{0}.in'.format(s)
            shutil.copyfile(src,dst)
            # Submit the job using that parameter file, but change it first
            jobid = submit_job(dst,jobid,s,temp,mm,'model_final_{0}.txt'.format(jobid))
        else:
            # On the first run, just submit the starting parameter file
            dst = starting_param_filename
            jobid = submit_job(dst,None,s,temp,mm)
        print("Jobid: {0}".format(jobid))


if __name__ == '__main__':
    main()
