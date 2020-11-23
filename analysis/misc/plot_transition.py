from sys import argv
import matplotlib.pyplot as plt

W_list = ['1.0000', '2.3300', '2.6600', '3.0000',  '3.3300', '3.6600', '4.0000', '4.3300', '4.6600', \
          '5.0000', '6.0000', '7.0000', '8.0000', '9.0000', '10.0000']
thresh = 0.0001
x_list = []
zero_list = []
for W in W_list:
    fname = 'dis_avg_overlap_W_{}_L_32_l_3_e_0.50.txt'.format(W)
    f = open(fname, 'r')
    data = []
    for line in f:
        data.append(float(line))
    L = len(data)
    z = 0
    for i in range(int(L/2), 0, -1):
        l = .5*(data[i] + data[L-1-i])
        if l < thresh:
            z = i
            break
    zero_list.append(L/2 - i)
    x_list.append(float(W))



plt.plot(x_list, zero_list, 'o-')
plt.xlabel("W")
plt.ylabel("r")
#plt.legend(loc='best')
plt.show()

