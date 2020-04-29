# % Allan Variance analyzer: This script takes in computed Allan Deviation,
# % averaging time (T), along with desired slopes and averaging times of
# % interest to extract noise coefficients of interest.
# % [sigmasOut,f,h] = AnalyzeAVAR(AD,dt,T,Slopes,Ts,plots,method)
# %
# % Outputs:
# % sigmasOut: Vector of noise coefficients of interest
# % f: figure handles to Allan Deviation and Slope plots
# % h: axis handles to plotted elements
# %
# % Inputs:
# % AD: An Allan deviation vector
# % T: The corresponding averaging time vector for Allan deviation
# % Slopes: A vector containing Allan deviation slopes of interest
# % Ts: The corresponding vector containing averaging times of interest for
# % extrapolation. Use NaN for slopes with no averaging time of interest such
# % as bias instability.
# % plots: 0 for no plots, 1 for plots
# % method: 1 for left-to-right slope search, 2 for min distance search
# %
import numpy as np
import math
import matplotlib.pyplot as plt

# def logdiff(X):
#     # m = size(X,dim);%返回矩阵的行数或列数，dim=1返回行数，dim=2返回列数
#     N = X.shape[0]
#     # print('X.shape: ', X.shape)
#     logdiff = np.diff(np.log10(X), axis=0)
#
#     # logdiff = np.zeros((N-1, 1))
#     # for ii in range(N-1):
#     #     logdiff[ii, 0] = np.log10(X[ii+1, 0]) - np.log10(X[ii, 0])
#     return logdiff


def AnalyzeAVAR(AD, T, Slopes, Ts, plots, method):
    """

    :param AD: nx1 ndarray
    :param T:  nx1 ndarray
    :param Slopes: list lenght=5
    :param Ts: list length=5
    :param plots: 1
    :param method:1
    :return:
    """
    h1_list = []
    h2_list = []
    if len(Slopes) != len(Ts):
        raise Exception("Each slope must have a corresponding Ts. Enter NaN if no Ts")

    # % Next, compute Log - Log slope of input Allan deviation signal.
    # Log - Log slope is: (Log10(Y(t + 1)) - Log10(Y(t))) / ((Log10(X(t + 1)) - Log10(X(t)))
    slope = np.diff(np.log10(AD), axis=0) / np.diff(np.log10(T), axis=0)   # nx1 ndarray
    # Create basic figure shells if plots were requested
    if plots:
        # h = np.zeros((2, len(Slopes)+1))
        # 数据图
        fig1 = plt.figure(num='fig1')
        ax1, = plt.loglog(T, AD, LineWidth=2)
        h1_list.append(ax1)
        plt.xlabel('Averaging Time')
        plt.ylabel('Allan Deviation')
        plt.title("Allan Deviation")
        # figure 2
        fig2 = plt.figure(num='fig2')
        ax2, = plt.semilogx(T[1:, 0], slope)
        h2_list.append(ax2)
        plt.xlabel('Averaging Time')
        plt.ylabel('Log-Log Slope')
        plt.title("Slope of Allan Deviation")
        # plt.savefig("Pictures/Allan_variance.png")
        c = ['r', 'b', 'k', 'g', 'y']
        f = ['fig1', 'fig2']
    else:
        h = 0
        f = 0
    # % Creaate output vector shell
    sigmasOut = np.zeros((len(Slopes), 1))
    Tbias = []
    # Now, step throgh each slope of interest, find the slope in the Allan
    # deviation signal,extrapolate to corresponding Ts of interest and record the coefficient found at that location.
    for ii in range(len(Slopes)):
        curr_slope = Slopes[ii]
        curr_t = Ts[ii]
        # method 1: use knowledge of Allan deviation shape, starting from the left(small T)
        # look occurence of desired slope.
        # method 2: use min distance method to locate point with closet matching slope to the one desired
        if method == 1:
            if math.copysign(1, slope[0, 0]) == -1.0:
                # k = find(X, n, direction)（其中 direction为'last'）查找与X中的非零元素对应的最后n个索引。direction
                # 的默认值为'first'，即查找与非零元素对应的前n个索引。
                # idx = find(slope > curr_slope, 1, 'first');
                if np.argwhere(slope > curr_slope).shape[0] != 0:  # # 找到了
                    idx = np.argwhere(slope > curr_slope)[0, 0]
                else:
                    idx = ""
            else:
                if np.argwhere(slope > curr_slope).shape[0] != 0:  # # 找到了
                    idx = np.argwhere(slope < curr_slope)[0, 0]
                else:
                    idx = ""
        elif method == 2:
            dist = (slope - curr_slope) ** 2
            idx = np.min(dist, axis=0)
        if idx != "":
            if curr_t is None:
                sigmasOut[ii, :] = AD[idx, :]
                Tbias = T[idx, 0]
                # linefun = @(t) AD(idx) + 0 * t; 匿名函数
                linefun = lambda t: AD[idx, 0] + 0*t
            else:
                linefun = lambda t: 10 ** (curr_slope * (np.log10(t)-np.log10(T[idx, 0])) + np.log10(AD[idx, 0]))
                sigmasOut[ii, :] = linefun(curr_t)

            if plots:
                plt.figure(num='fig1')
                ax3, = plt.plot(T[idx, :], AD[idx, :], marker='s', color=c[ii])
                h1_list.append(ax3)
                # print(linefun(T))
                plt.plot(T, linefun(T), linestyle='--', lineWidth=1.2, color=c[ii])
                plt.figure(num='fig2')
                ax4, = plt.plot(T[idx+1, :], slope[idx, :], marker='s', color=c[ii])
                h2_list.append(ax4)
    return sigmasOut, Tbias,  f, [h1_list, h2_list]







