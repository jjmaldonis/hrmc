import sys, os

def reformat(f):
    of = 'reformatted/'+f
    with open(f, 'r') as input, open(of, 'w') as output:
        line = input.readline() # Comment 1
        output.write(line)
        line = input.readline() # Comment 2
        output.write(line)
        line = input.readline() # Comment 3
        output.write(line)

        line = input.readline() # Number of elements and elements
        output.write(line)
        sline = line.strip().split()
        while '' in sline:
            sline.remove('')
        nelements = int(sline[0])
        
        line = input.readline() # nrho drho nr dr cutoff
        output.write(line)
        sline = line.strip().split()
        while '' in sline:
            sline.remove('')
        nrho = int(sline[0])
        drho = float(sline[1])
        nr = int(sline[2])
        dr = float(sline[3])
        cutoff = float(sline[4])

        for elem in range(nelements):
            line = input.readline()
            output.write(line)
            i = 0
            while i < nrho:
                line = input.readline().strip().split()
                while '' in line:
                    line.remove('')
                #for num in [float(x) for x in line]:
                    #output.write(str(num)+'\n')
                for num in line:
                    output.write(num+'\n')
                    i += 1
            i = 0
            while i < nr:
                line = input.readline().strip().split()
                while '' in line:
                    line.remove('')
                #for num in [float(x) for x in line]:
                    #output.write(str(num)+'\n')
                for num in line:
                    output.write(num+'\n')
                    i += 1

        for elem1 in range(nelements):
            for elem2 in range(nelements):
                if elem1 >= elem2:
                    i = 0
                    while i < nr:
                        line = input.readline().strip().split()
                        while '' in line:
                            line.remove('')
                        #for num in [float(x) for x in line]:
                            #output.write(str(num)+'\n')
                        for num in line:
                            output.write(num+'\n')
                            i += 1


def main():
    f = sys.argv[1]
    if not os.path.exists(os.path.split(os.path.abspath(f))[0] + '/reformatted'):
        os.makedirs(os.path.split(os.path.abspath(f))[0] + '/reformatted')
    reformat(f)

if __name__ == '__main__':
    main()
