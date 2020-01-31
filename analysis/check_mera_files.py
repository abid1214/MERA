import numpy as np
from sys import argv

def get_mera_fname(W, e, dis):
    data_dir = "../data/L_{}_l_{}/".format(L, l)
    return "{}W_{}/epsilon_{:02}_{:02}.txt".format(data_dir, W, e, dis)

def get_mera_minmax_energies(W, dis):
    E_min, _,_,_,_ = load_diagnostics(get_mera_fname(W, 0, dis))
    E_max, _,_,_,_ = load_diagnostics(get_mera_fname(W, num_energies-1, dis))
    return E_min, E_max

def load_diagnostics(fname):
    E_dict = {}
    with open(fname, 'r') as fp:
        data = fp.readlines()
        assert len(data) > 10, "lines in file = {}".format(len(data))
        energies = np.array(data[-10].split()[2:]).astype(float)
        e_diffs = energies[1:] - energies[:-1]
        spacings = [min(e_diffs[i:i+2])/max(e_diffs[i:i+2]) for i in range(len(e_diffs)-2)]
        ll = data[-1].split()
        E = float(ll[2])
        var = float(ll[3])
        EE = float(ll[4])
        Sz = float(ll[5])
        r = np.mean(spacings)
    return E, var, EE, r, Sz

if __name__ == "__main__":

    L = int(argv[1])
    l = int(argv[2])

    W_list = [0.0001] + list(range(1,11))
    num_W = len(W_list)
    num_dis = 100
    num_energies = 11


    succ = 0
    fail = 0
    f = open("rerun_L_{}_l_{}.txt".format(L,l), 'w')
    for W_idx, W in enumerate(W_list):
        for dis in range(num_dis):
            for e in range(num_energies):
                try:
                    mera_fn = get_mera_fname(W, e, dis)
                    E, var, EE, r, Sz = load_diagnostics(mera_fn)
                    E_min_mera, E_max_mera = get_mera_minmax_energies(W, dis)
                    ep_mera = (E - E_min_mera)/(E_max_mera - E_min_mera)
                    succ += 1
                except Exception as ex:
                    print("W = {}, dis = {}, e = {}: {}".format(W, dis, e, ex))
                    f.write("{} {} {}\n".format(W, dis, e))
                    fail += 1
    print("successes: {}, failures: {}".format(succ, fail))
    f.close()

