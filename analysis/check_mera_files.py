import numpy as np
from sys import argv
from os.path import exists

def get_mera_fname(W, e, dis):
    data_dir = "../data/L_{}_l_{}/".format(L, l)
    return "{}W_{}/epsilon_{:02}_{:02}.txt".format(data_dir, W, e, dis)

def check_last_line(fname):
    if not exists(fname):
        return False, "file doesn't exist"
    with open(fname, 'r') as fp:
        data = fp.readlines()
        if len(data) == 0:
            return False, "file empty"
        last_line = np.array(data[-1].split())
        if len(last_line) != 6:
           return  False, "last line error: {}".format(last_line)
    return True, ""

if __name__ == "__main__":

    L = int(argv[1])
    l = int(argv[2])

<<<<<<< Updated upstream
    #W_list = [0.0001] + list(range(1,11))
=======
>>>>>>> Stashed changes
    W_list = [0.33, 0.66, 1.33, 1.66, 2.33, 2.66, 3.33, 3.66, 4.33, 4.66, 5.33]
    num_W = len(W_list)
    num_dis = 100
    num_energies = 11


    succ = 0
    fail = 0
    f = open("rerun_L_{}_l_{}.txt".format(L,l), 'w')
    for W_idx, W in enumerate(W_list):
        for dis in range(num_dis):
            for e in range(num_energies):
                mera_fn = get_mera_fname(W, e, dis)
                good_file, msg = check_last_line(mera_fn)
                if good_file:
                    succ +=1
                else:
                    print("W = {}, dis = {}, e = {}: {}".format(W, dis, e, msg))
                    f.write("{} {} {}\n".format(W, dis, e))
                    fail += 1
    print("successes: {}, failures: {}".format(succ, fail))
    f.close()

