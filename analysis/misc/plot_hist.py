from sys import argv
import matplotlib.pyplot as plt

for fname in argv[1:]:
    f = open(fname, 'r')
    data = []
    for line in f:
        data.append(float(line))
    plt.hist(data, bins=50, histtype=u'step', label=fname, alpha=0.5)

plt.legend(loc='best')
plt.show()

