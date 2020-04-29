# -*- coding: utf-8 -*- 
"""
Project: Allan_Project
Creator: Administrator
Create time: 2020-02-19 16:25
IDE: PyCharm
Introduction:
"""
import matplotlib.pyplot as plt
import numpy as np
import math
import time


def allan(omega, fs, maxNumM):
    """

    :param omega: 输入样本集合
    :param fs: 采样频率
    :param maxNumM: 最大量化平均因子数
    :return: Tau
             allan variance
    """
    t0 = 1/fs
    N = omega.shape[0]
    M = omega.shape[1]
    n = [2 ** i for i in range(math.floor(np.log2(N / 2)) + 1)]
    maxN = n[-1]
    log_min = np.log10(1)
    log_max = np.log10(maxN)
    # 去重，取整
    m = np.unique(np.ceil(np.logspace(log_min, log_max, maxNumM)).astype(int))
    T = m * t0
    theta = np.cumsum(omega, axis=0) * t0
    sigma2 = np.zeros((len(T), M))

    theta2 = theta * 2
    for i in range(len(m)):
        sigma2[i, :] = np.sum((theta[2 * m[i]:, :] - theta2[m[i]:N - m[i], :] + theta[0:N - 2 * m[i], :]) ** 2, axis=0)

    sigma2 = sigma2 / (np.tile(2 * (N - 2 * m) * T ** 2, (M, 1))).T
    sigma = np.sqrt(sigma2)

    return T, sigma


if __name__ == "__main__":
    t1 = time.time()
    matfn = "D:\\GyroAllan-master\\Software\\data.txt"
    data = np.loadtxt(matfn)
    data = data[:, 2:5]*3600
    print(data.shape)
    [T, sigma] = allan(data, 100, 100)
    plt.loglog(T, sigma, ".")
    plt.xlabel('time:sec')
    plt.ylabel('Sigma: deg/h')
    plt.legend(['X axis', 'Y axis', "Z axis"])
    plt.savefig("pictures/Allan_variance.png")
    plt.show()
    print(time.time()-t1)