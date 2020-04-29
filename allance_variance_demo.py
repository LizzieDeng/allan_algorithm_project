# -*- coding: utf-8 -*-
import math
import numpy as np
from AnalyzeAVAR import AnalyzeAVAR
from computerAVAR import computer_allan_variance
import matplotlib.pyplot as plt
import time
# % Allan Variance: Data simulation and analysis demo script.
# %if data is short :plot_them = 0


def main():
    plot_them = 1
    filename = "data/wz.csv"
    data = np.genfromtxt(open(filename, "rb"), delimiter=',')
    dt = 1 / 8
    Y = data[0:, 1].reshape(len(data[0:, 1]), 1)
    AllanSigma, T = computer_allan_variance(Y, dt)
    T = T.reshape(len(T), 1)
    # % AllanSigma is expressed in deg/s (gyro) or m/s^2 (accel), we'll convert
    # % the denominator to hours prior to anlayzing so we don't have to scale the
    # % outputs later.
    AllanSigma = AllanSigma * 3600     # % Now in deg/hr (gyros) or m/s^2/hr(accels)
    T = T/3600  # % Now in hrs
    # % Finally, feed AllanSigma and T to analyzer script to extract noise
    # % parameters and compare to those we first introducted into Y.
    # % The script takes in Log-Log slopes of interest along with respective T
    # % values to extrapolate to prior to reading desired Allan Deviation. Per
    # % the literature, the slopes of interest are: -1, -1/2, 0 and 1/2;
    slopes = [-1, -0.5, 0, 0.5, 1]
    # slopes = [-1, -0.497, 0, 0.5, 1]
    # % Also per literature, the corresponding T values to extrapolate to are:
    # % T=sqrt(3), T=1, N/A, and T=3. Note we've already converted T to hours at
    # % this point, so there is no need to multiply these values by 3600.
    Ts = [math.sqrt(3), 1, None, 3, math.sqrt(2)]
    # % The script will analyse AllanSigma, look for the slopes, extrapolate to
    # % Ts and report the appropriate AllanSigma found. The last parameter is a
    # % plot option.
    [sigmasOut, Tbias, figs, hs] = AnalyzeAVAR(AllanSigma, T, slopes, Ts, 1, 1)
    # % We finally take the output of the analyzer and apply any other scaling
    # % per the literature. The only scaling needed at this point is for bias
    # % instability, where we divide by sqrt((2*log(2)/pi)).
    sigmasOut[2, 0] = sigmasOut[2, 0] / math.sqrt((2*np.log(2)/math.pi))
    print(' Quantization:%0.2e [deg]\n', sigmasOut[0, :])
    print('Random Walk:%0.2e [deg/sqrt{hr}]\n', sigmasOut[1, :])
    print('Bias Instability:%0.2e [deg/hr]\n', sigmasOut[2, :])
    print('Rate Random Walk:%0.2e [deg/hr/sqrt{hr}]\n', sigmasOut[3, :])
    print('Rate Ramp:%0.2e [deg/hr/hr]\n', sigmasOut[4, :])
    # % sigmasOut(1) % Quantization  deg (gyros) OR m/s (accels)
    # % sigmasOut(2) %  Random Walk  deg/sqrt(hr) (gyros) OR m/s/sqrt(hr) (accels)
    # % sigmasOut(3)%  Bias Instability  deg/hr (gyros) OR m/s/hr (accels)
    # % sigmasOut(4)% Rate Random Walk deg/hr/sqrt(hr) (gyros) OR m/s/hr/sqrt(hr)(accels)
    # % sigmasOut(5) % Rate Ramp deg/hr/hr
    if plot_them:
        sigmaQ = 'Quantization:%0.2e [deg]'% sigmasOut[0, :]
        sigmaRW = 'Random Walk:%0.2e [deg/sqrt(hr)]' % sigmasOut[1, :]
        sigmaBias = 'Bias Instability:%0.2e [deg/hr]' % sigmasOut[2, :]
        sigmaRRW = 'Rate Random Walk:%0.2e [deg/hr/sqrt(hr)]' % sigmasOut[3, :]
        sigmaRR = 'Rate Ramp:%0.2e [deg/hr/hr]' % sigmasOut[4, :]
        plt.figure(num=figs[0])
        plt.title("Simulated Data: Allan Deviation")
        plt.legend(hs[0][1:], [sigmaQ, sigmaRW, sigmaBias, sigmaRRW, sigmaRR],  loc='lower right')
        plt.xlabel("Averaging Time, $\tau$ (hr)")
        plt.ylabel("Allan Deviation, $\sigma(\tau)$ (deg/hr)")
        plt.figure(num=figs[1])
        plt.title("Simulated Data: Allan Deviation Slope")
        plt.legend(hs[1][1:], ['-1: Quantization', '-1/2: Random Walk','0: Bias Instability', '+1/2: Rate Random Walk', '+1: Rate Ramp'], loc='upper left')
        plt.xlabel('Averaging Time, $\tau$ (hr)')
        plt.ylabel("Allan Deviation Slope")
        plt.show()

    # % % Bias stability and Bias
    N = 10 # % Ns
    num = int(N / dt) # % % m每组内数据个数
    m = int(math.floor(len(Y) / num))  # % % 共可以分成m组数据
    gx = np.zeros((m, 1))
    Xbais_all = np.mean(Y)
    for i in range(m-1):
        gx[i, :] = np.mean(Y[i * num: (i+1) * num, 0])
    gx = gx[1:len(gx)-1, :]
    # 矩阵标准差
    Xbais_stability = np.std(gx) * 3600
    print('Bias Stability 10s:%0.2e [deg/hr]\n', Xbais_stability)


if __name__ == "__main__":
    t1 = time.time()
    main()
    print("main:", time.time()-t1)