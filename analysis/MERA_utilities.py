import numpy as np
from scipy.signal import convolve2d

E_IDX, VAR_IDX, LOGVAR_IDX, EE_IDX, R_IDX, SZ_IDX, EP_MERA_IDX, EP_DMRG_IDX, \
        W_IDX, DIS_IDX, EP_IDX, L_IDX, LL_IDX = list(range(13))

def get_mera_fname(L, l, W, e, dis):
    ''' grabs the file name for MERA config (L, l, W, e, dis) '''
    data_dir = "../data/L_{}_l_{}/".format(L, l)
    return "{}W_{}/epsilon_{:02}_{:02}.txt".format(data_dir, W, e, dis)


def get_dmrg_fname(L, W, dis):
    ''' grabs the file name for DMRG config (L, W, dis) '''
    data_dir = "../data/dmrg/L_{}/".format(L)
    return "{}W_{}/epsilon_{:02}.txt".format(data_dir, W, dis)


def load_diagnostics(L, l, W, e, dis):
    ''' loads MERA config (L, l, W, e, dis) data '''
    fname = get_mera_fname(L, l, W, e, dis)
    with open(fname, 'r') as fp:
        data = fp.readlines()
        energies = np.array(data[-4].split()[2:]).astype(float)
        e_diffs = energies[1:] - energies[:-1]
        r = np.mean([min(e_diffs[i:i+2])/max(e_diffs[i:i+2]) for i in \
                range(len(e_diffs)-2)])
        ll = data[-1].split()
        E, var, EE, Sz = [float(ll[i]) for i in range(2,6)]
    return E, var, EE, r, Sz


def get_dmrg_minmax_energies(L, W, dis):
    ''' returns the min and max energies of a DMRG config (L, W, dis) '''
    fname = get_dmrg_fname(L, W, dis)
    with open(fname, 'r') as fp:
        data = fp.readlines()[-1].split()
        return float(data[0]), float(data[1])


def get_mera_minmax_energies(L, l, W, dis, num_energies):
    ''' returns the min and max energies of a MERA config (L, l, W, dis) '''
    E_min, _,_,_,_ = load_diagnostics(L, l, W, 0, dis)
    E_max, _,_,_,_ = load_diagnostics(L, l, W, num_energies-1, dis)
    return E_min, E_max


def load_all_data(L_list, l_list, W_list, num_dis, num_energies):
    ''' gets MERA and DMRG data for all W, disorders and energies for a
        given (L, l)
    '''
    num_L, num_l, num_W, num_params = len(L_list), len(l_list), len(W_list), 13
    all_data = np.zeros((num_L, num_l, num_W, num_dis, num_energies, num_params))
    for L_idx, L in enumerate(L_list):
        for l_idx, l in enumerate(l_list):
            print("L = {}, l = {}".format(L, l))
            for W_idx, W in enumerate(W_list):
                for dis in range(num_dis):
                    for e in range(num_energies):
                        E, var, EE, r, Sz = load_diagnostics(L, l, W, e, dis)
                        E_min_dmrg, E_max_dmrg = get_dmrg_minmax_energies(L, W, dis)
                        E_min_mera, E_max_mera = get_mera_minmax_energies(L, l, W, \
                                dis, num_energies)
                        ep_mera = (E - E_min_mera)/(E_max_mera - E_min_mera)
                        ep_dmrg = (E - E_min_dmrg)/(E_max_dmrg - E_min_dmrg)

                        all_data[L_idx][l_idx][W_idx][dis][e][:] = \
                                np.array([E, var, np.log10(var), EE, r, Sz, ep_mera, \
                                ep_dmrg, W, dis, e/(num_energies - 1), L, l])
    return all_data


def dis_avg_data(all_data):
    ''' returns a disorder average of the data'''
    return np.mean(all_data, axis=3)


def get_W_data(data, W_idx):
    ''' gets flattened array of data for a given W '''
    num_L, num_l, num_W, num_dis, num_energies, num_params = all_data.shape
    return data[:,:,W_idx,:,:,:].reshape(num_L*num_l*num_dis*num_energies, num_params)


def get_dis_data(all_data, dis):
    ''' gets flattened array of data for a given disorder '''
    num_L, num_l, num_W, num_dis, num_energies, num_params = all_data.shape
    return all_data[:,:,:,dis,:,:].reshape(num_L*num_l*num_W*num_energies, num_params)


def get_L_data(all_data, L_idx):
    ''' gets flattened array of data for a given L '''
    num_L, num_l, num_W, num_dis, num_energies, num_params = all_data.shape
    return all_data[L_idx,:,:,:,:,:].reshape(num_l*num_W*num_dis*num_energies, num_params)


def get_l_data(all_data, l_idx):
    ''' gets flattened array of data for a given L '''
    num_L, num_l, num_W, num_dis, num_energies, num_params = all_data.shape
    return all_data[:,l_idx,:,:,:,:].reshape(num_L*num_W*num_dis*num_energies, num_params)


def flatten_data(all_data):
    num_params = all_data.shape[-1]
    return all_data.reshape(-1, num_params)


def get_ep_data(flattened_data, ep, tol = 0.05):
    ''' gets flattened array of data for a given epsilon '''
    d = []
    for i in range(len(flattened_data)):
        data = flattened_data[i]
        if data[EP_DMRG_IDX] <= ep + tol and data[EP_DMRG_IDX] >= ep - tol:
            d.append(data)
    return np.array(d)


def get_lists(flattened_data):
    ''' separates the flattened data into lists of parameters '''
    num_params = flattened_data.shape[-1]
    return [flattened_data[:,i] for i in range(num_params)]

def get_sorted_data(flattened_data, x_idx):
    ''' given a flattened array , sorted the data by the x parameter'''
    list_arr = get_lists(flattened_data)
    sorted_arr = [list_arr[x_idx]] + list_arr
    sorted_arr = list(zip(*sorted(zip(*sorted_arr))))
    return sorted_arr[1:]

def smooth_data(arr, N=10):
    ''' smoothes 1d data '''
    return np.convolve(arr, np.ones((N,))/N, mode='valid')

def smooth_data2d(arr, N=3):
    ''' smoothes 2d data '''
    return convolve2d(arr, np.ones((N,N))/(N*N), mode='valid')


