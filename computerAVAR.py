# -*- coding: utf-8 -*-
import math
import numpy as np


def computer_allan_variance(ydot, dt):
    pts = 500    # % Arbitrary size of log spaced vector
    # print(ydot.shape)
    N = ydot.shape[0]
    M = ydot.shape[1]
    # determine largest bin size
    n = [2 ** i for i in range(math.floor(np.log2(N / 2)) + 1)]
    maxN = n[-1]
    endLogINC = np.log10(maxN)
    # create log spaced vector average factor
    m = np.unique(np.ceil(np.logspace(0, endLogINC, pts)).astype(int))
    t = m * dt
    # % integration of samples over time to obtain output angle
    theta = np.cumsum(ydot, axis=0) * dt
    # % loop over the various cluster sizes
    AV = np.zeros((len(t), M))
    theta2 = theta * 2
    for i in range(len(m)):
        AV[i, :] = np.sum((theta[2 * m[i]:, :] - theta2[m[i]:N - m[i], :] + theta[0:N - 2 * m[i], :]) ** 2, axis=0)
    # Allan Variance
    AV = AV / (np.tile(2 * (N - 2 * m) * t ** 2, (M, 1))).T
    ad = np.sqrt(AV)
    return ad,  t


def read_csv(filename):
    data = np.genfromtxt(open(filename, "rb"), delimiter=',')
    return data


if __name__ == "__main__":
    # import tkinter
    # print(tkinter.__version__)
    x = np.arange(6).reshape(2, 3) + 2
    print(x)
    x1 = x *2
    print(x1)
    print(x/x1)
    print(x/np.linalg.inv(x1))

    # curr_t = None
    # if curr_t is None:
    #     print('None')
    # print(np.isnan(curr_t))
    # filename = "data/wz.csv"
    # data = read_csv(filename)
    # dt = 1/8
    # Y = data[:, 1]
    # Y = Y.reshape(len(Y), 1)
    # ad, t = computer_allan_variance(Y, dt)
    # print(type(ad), type(t), ad.shape, t.shape)
    # # 水平叠加
    # ad_t = np.column_stack((ad, t))
    # np.savetxt("data/computeravar_correct.csv", ad_t, delimiter=',')

