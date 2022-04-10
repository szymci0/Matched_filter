import random
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal
from scipy import signal
from scipy.signal import find_peaks
from peakdetect import peakdetect


def periodic_square(length):
    t = np.linspace(0, length, 1, endpoint=False)
    sig = np.sin(2 * np.pi * t)
    pwm = signal.square(2 * np.pi * 30 * t, duty=(sig + 1) / 2)

    return pwm


def rand_square_gen(length, num_of_pulses):
    # inicjalizacja szeregu o danej dlugosci
    square_sig = np.zeros(length)
    t = np.linspace(0, length, 1)
    random_list = []

    # wytwarzanie losowych przedziałów
    for i in range(num_of_pulses * 2):
        n = random.randint(length / num_of_pulses * 0.5, length * 0.9)
        random_list.append(n)

    random_list.sort()
    random_comps = list(zip(random_list, random_list[1:]))[::2]

    # przydzielanie wartosci 1 do losowych przedziałów
    for comp in random_comps:
        square_sig[comp[0] : comp[1]] = 1

    return (signal.square(t, square_sig) + 1) / 2


def matched_filtering(sig, sig_noise):
    flipped = sig[::-1]
    sig_lengthened = np.append(flipped, np.zeros(LENGTH))

    return signal.lfilter(sig_noise, 1, sig_lengthened) / LENGTH
    # return signal.decimate(signal.lfilter(sig_noise, 1, sig_lengthened), 2)


def squares(length, przedzialy):
    square_sig = np.zeros(length)
    t = np.linspace(0, length, 1)

    przedzialyt = list(zip(przedzialy, przedzialy[1:]))[::2]

    for comp in przedzialyt:
        square_sig[comp[0]: comp[1]] = 1

    return (signal.square(t, square_sig) + 1) / 2


LENGTH = 1000

# generacja sygnalu
sygnal = squares(LENGTH, [150, 250, 350, 450, 550, 650])

# zaszumianie
noise = np.random.normal(0, 0.2, LENGTH)
sygnal_szum = sygnal + noise

# aplikacja filtru dopasowanego
square_sig = np.zeros(LENGTH)
t = np.linspace(0, LENGTH, 1)
square_sig[0:100] = 1

sig_lengthened = (signal.square(t, square_sig) + 1) / 2
matched_filter = signal.lfilter(sygnal_szum, 1, sig_lengthened) / 100


# znajdowanie szczytów

end, _ = find_peaks(matched_filter, prominence=0.5)
begin = []
for peak in end:
    begin.append(peak-99)

# pika = peakdetect(matched_filter, lookahead=1, delta=0.5)
# higherPeaks = np.array(pika[0])
# lowerPeaks = np.array(pika[1])


# rysowanie wykresów
fig, (f_orig, f_peaks, f_corr, f_smooth, f_peaki) = plt.subplots(5, 1, sharex=True)

f_orig.plot(sygnal)
f_smooth.plot(sig_lengthened)
f_corr.plot(matched_filter)
f_peaks.plot(sygnal_szum)
f_peaki.plot(begin,sygnal_szum[begin], "xr")
f_peaki.plot(end,sygnal_szum[end], "xr")
f_peaki.plot(sygnal_szum)
# f_pika.plot(higherPeaks[:,0], higherPeaks[:,1], 'xr')
# f_pika.plot(lowerPeaks[:,0], lowerPeaks[:,1], 'xk')
# f_pika.plot(sygnal)
plt.show()
