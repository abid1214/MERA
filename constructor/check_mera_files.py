import numpy as np
from sys import argv
from os.path import exists

def get_mera_fname(W, e, dis):
    data_dir = "../data/L_{}_l_{}/".format(L, l)
    return "{}W_{}/epsilon_{:02}_{:02}.txt".format(data_dir, W, e, dis)

def get_dmrg_fname(W, dis):
    data_dir = "../data/dmrg/L_{}/".format(L)
    return "{}W_{}/epsilon_{:02}.txt".format(data_dir, W, dis)


def check_last_line(fname):
    if not exists(fname):
        return False, "file doesn't exist"
    with open(fname, 'r') as fp:
        data = fp.readlines()
        if len(data) == 0:
            return False, "file empty"
        last_line = np.array(data[-1].split())
        if len(last_line) != 6 or "Info:" not in last_line:
            return  False, "last line error: {}".format(last_line[:10])
    return True, ""

def check_dmrg_last_line(fname):
    try:
        with open(fname, 'r') as fp:
            data = fp.readlines()[-1].split()
            minE, maxE = float(data[0]), float(data[1])
    except Exception as e:
        return False, e
    return True, ""

if __name__ == "__main__":

    L = int(argv[1])

    W_list = [0.0001] + list(range(1,11))
    W_list += [0.33, 0.66, 1.33, 1.66, 2.33, 2.66, 3.33, 3.66, 4.33, 4.66, 5.33]
    num_W = len(W_list)
    num_dis = 100
    num_energies = 11

    if len(argv) >2:
        l = int(argv[2])
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

    succ = 0
    fail = 0
    f = open("dmrg_L_{}.txt".format(L), 'w')
    for W_idx, W in enumerate(W_list):
        for dis in range(num_dis):
            mera_fn = get_dmrg_fname(W, dis)
            good_file, msg = check_dmrg_last_line(mera_fn)
            if good_file:
                succ +=1
            else:
                print("W = {}, dis = {}: {}".format(W, dis, msg))
                f.write("{} {}\n".format(W, dis))
                fail += 1
    print("successes: {}, failures: {}".format(succ, fail))
    f.close()

