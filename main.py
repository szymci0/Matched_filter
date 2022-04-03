import random
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal


def square_gen(length, num_of_pulses):
    # inicjalizacja szeregu o danej dlugosci
    square_sig = np.zeros(length)
    t = np.linspace(0, length, 1)
    random_list = []

    # wytwarzanie losowych przedziałów
    for i in range(num_of_pulses * 2):
        n = random.randint(length / num_of_pulses * 0.1, length * 0.9)
        random_list.append(n)

    random_list.sort()
    random_comps = list(zip(random_list, random_list[1:]))[::2]

    # przydzielanie wartosci 1 do losowych przedziałów
    for comp in random_comps:
        square_sig[comp[0] : comp[1]] = 1

    return (signal.square(t, square_sig) + 1) / 2


def matched_filtering(sig, sig_noise):
    sig_lengthened = np.append(np.zeros(LENGTH), sig)
    flipped = sig_lengthened[::-1]

    return signal.decimate(signal.lfilter(sig_noise, 1, flipped), 2)


LENGTH = 1000
NUM_OF_PULS = 5

# generacja sygnalu
sygnal = square_gen(LENGTH, NUM_OF_PULS)

# zaszumianie
noise = np.random.normal(0, 0.8, LENGTH)
sygnal_szum = sygnal + noise

# aplikacja filtru dopasowanego
matched_filter = matched_filtering(sygnal, sygnal_szum)


# rysowanie wykresów
fig, (f_orig, f_smooth, f_corr, f_peaks) = plt.subplots(4, 1, sharex=False)

f_orig.plot(sygnal)
f_smooth.plot(sygnal)
f_corr.plot(matched_filter)
f_peaks.plot(sygnal_szum)
plt.show()
