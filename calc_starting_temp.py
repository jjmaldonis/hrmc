import sys, math
import numpy as np

def main():
    # Import a file containing values for del-chi for the first ~1000 moves
    acceptance = 0.90
    delfile = sys.argv[1]
    content = open(delfile).readlines()
    content = [float(line.strip()) for line in content]

    values = [abs(x) for x in content]
    values = np.array(values)
    h,xx = np.histogram(values,bins=100)
    xx = [(xx[i-1]+xx[i])/2.0 for i in range(1,len(xx))]
    remove = []
    for i,x in enumerate(h):
        if x == 0:
            remove.append(i)
    xx = np.delete(xx,remove)
    h = [math.log(x) for x in h if x!=0]
    tau = - np.polyfit(xx,h,1)[0]
    temp = 1.0/tau*math.log(1-acceptance/2.0) / (8.6171e-05*math.log(acceptance))
    print("Try Temperature = {0} if you want a starting acceptance rate of {1}%".format(temp,acceptance*100))

    ## Use this to check if you have enough points
    ## Plot jason.txt in Igor and see if it has converged.
    #f = open('jason.txt','w')
    #for i in range(100,len(content),100):
    #    values = [abs(x) for x in content[0:i]]
    #    values = np.array(values)
    #    h,xx = np.histogram(values,bins=100)
    #    xx = [(xx[i-1]+xx[i])/2.0 for i in range(1,len(xx))]
    #    remove = []
    #    for i,x in enumerate(h):
    #        if x == 0:
    #            remove.append(i)
    #    xx = np.delete(xx,remove)
    #    h = [math.log(x) for x in h if x!=0]
    #    tau = - np.polyfit(xx,h,1)[0]
    #    temp = 1.0/tau*math.log(1-acceptance/2.0) / (8.6171e-05*math.log(acceptance))
    #    f.write('{0}\n'.format(temp))
    #f.close()



if __name__ == '__main__':
    main()
