
import os
import re
import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append(r"C:\IDE\Anaconda")

def parse_data():
    numeric_const_pattern = r"""[-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?"""
    _momentum_list = list()

    for line in open('data.txt', encoding='utf-8', errors = 'ignore'):
        _find = re.findall(numeric_const_pattern, line)
        if len(_find) > 0:
            if len(_find) == 3:
                _momentum_list.append(float(_find[2]))
            elif len(_find) == 2:
                _momentum_list.append(float(_find[1]))

    return _momentum_list


def plot_data(xs, momentum):

    try:

        plt.plot(xs, momentum_list)
        #plt.axis([0, np.max(xs), 0, np.max(momentum_list)])
        plt.title("Cylindrical shell momentum")
        plt.xlabel("X")
        plt.ylabel("M[X]")
        plt.show()

    except:

        print('Maybe you got quazi singular matrix')


    return 0


momentum_list = parse_data()

xs = [i for i in range(0, len(momentum_list))]


print('\n<X> = {}'.format(np.mean(xs)))
print('<M[X]> = {}'.format(np.mean(momentum_list)))

plot_data(xs, momentum_list)




