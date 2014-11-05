import sys, math

def main():
    # Import a file containing values for del-chi for the first ~1000 moves
    delfile = sys.argv[1]
    content = open(delfile).readlines()
    content = [float(line.strip()) for line in content]
    avg = sum(content)/len(content)
    acceptance = 0.875
    temp = -avg/(8.6171e-05*math.log(acceptance))
    print("Try Temperature = {0} if you want a starting acceptance rate of {1}%".format(temp,acceptance*100))

    temp10 = -avg/(8.6171e-05*math.log(0.1))
    s = math.log(temp10/temp)/math.log(math.sqrt(0.7))*200000.0
    print("Then the acceptance rate will reach 10% when the number of steps reaches {0}".format(s))
    temp10 = -avg/(8.6171e-05*math.log(0.01))
    s = math.log(temp10/temp)/math.log(math.sqrt(0.7))*200000.0
    print("                     and will reach 1%  when the number of steps reaches {0}".format(s))



if __name__ == '__main__':
    main()
