import math
import numpy as np
import matplotlib.pyplot as plt


class Hodgkin:
    def __init__(self):
        self.dt = 0.01
        self.t_max = 100

        self.v = -65
        self.v_t = []

        self.n = 0.3177
        self.n_t = []

        self.m = 0.0529
        self.m_t = []

        self.h = 0.5961
        self.h_t = []

        self.external_current = 0.200
        self.g_l = 0.003
        self.g_k = 0.36
        self.g_na = 1.2
        self.e_l = -54.387
        self.e_k = -77
        self.e_na = 50
        self.c_m = 0.010

    def step(self):

        self.v_t.append(self.v)
        self.n_t.append(self.n)
        self.m_t.append(self.m)
        self.h_t.append(self.h)
        self.v += self.dv()
        self.n += self.dn()
        self.m += self.dm()
        self.h += self.dh()

    def dv(self):
        return (-self.i_m() + self.external_current) / self.c_m * self.dt

    def i_m(self):
        return self.g_l * (self.v - self.e_l) + self.g_k * math.pow(self.n, 4) * (self.v - self.e_k)\
            + self.g_na * math.pow(self.m, 3) * self.h * (self.v - self.e_na)

    def dn(self):
        return (self.n_inf() - self.n) / self.tau_n() * self.dt

    def tau_n(self):
        return 1 / (self.alpha_n() + self.beta_n())

    def n_inf(self):
        return self.alpha_n() / (self.alpha_n() + self.beta_n())

    def alpha_n(self):
        return 0.01 * (self.v + 55) / (1 - math.exp(-0.1 * (self.v + 55)))

    def beta_n(self):
        return 0.125 * math.exp(-0.0125 * (self.v + 65))

    def tau_m(self):
        return 1 / (self.alpha_m() + self.beta_m())

    def tau_h(self):
        return 1 / (self.alpha_h() + self.beta_h())

    def dm(self):
        return (self.m_inf() - self.m) / self.tau_m() * self.dt

    def dh(self):
        return (self.h_inf() - self.h) / self.tau_h() * self.dt

    def m_inf(self):
        return self.alpha_m() / (self.alpha_m() + self.beta_m())

    def h_inf(self):
        return self.alpha_h() / (self.alpha_h() + self.beta_h())

    def alpha_m(self):
        return 0.1 * (self.v + 40) / (1 - math.exp(-0.1 * (self.v + 40)))

    def beta_m(self):
        return 4 * math.exp(-0.0556 * (self.v + 65))

    def alpha_h(self):
        return 0.07 * math.exp(-0.05 * (self.v + 65))

    def beta_h(self):
        return 1 / (1 + math.exp(-0.1 * (self.v + 35)))


def a():
    h = Hodgkin()
    for step in range(int(h.t_max / h.dt)):
        h.step()

    for data in [(h.v_t, "Voltage", "221", "mV"), (h.n_t, "n", "222"), (h.m_t, "m", "223"), (h.h_t, "h", "224")]:
        plt.subplot(data[2])
        plt.plot(np.arange(0, len(data[0]) * h.dt, h.dt), data[0])
        plt.xlabel("t/ms")
        plt.ylabel("" if len(data) <= 3 else data[3])

        plt.title(data[1])

    plt.show()


def b():
    firing_rates = []
    external_currents = list(range(0, 500, 10))
    for external in external_currents:
        h = Hodgkin()
        h.external_current = external / 1000
        for step in range(int(h.t_max / h.dt)):
            h.step()

        num_spikes = 0
        is_greater_zero = False
        for v in h.v_t:
            if v > 0 and not is_greater_zero:
                num_spikes += 1
                is_greater_zero = True
            elif v < 0:
                is_greater_zero = False

        firing_rates.append(num_spikes / h.t_max)

    plt.plot(external_currents, firing_rates)
    plt.xlabel("I_e/A / nA/mm^2")
    plt.ylabel("Firing rate/Hz")
    plt.show()


def c():
    h = Hodgkin()
    h.external_current = -50 / 1000
    for step in range(int(h.t_max / h.dt)):
        h.step()
        if step > 500:
            h.external_current = 0

    for data in [(h.v_t, "Voltage", "221", "Voltage/mV"), (h.n_t, "n", "222"), (h.m_t, "m", "223"), (h.h_t, "h", "224")]:
        plt.subplot(data[2])
        plt.plot(np.arange(0, len(data[0]) * h.dt, h.dt), data[0])
        plt.xlabel("t/ms")
        plt.ylabel("" if len(data) <= 3 else data[3])

        plt.title(data[1])

    plt.show()


print("External current = 200nA/mm^2...")
a()
print("Plot of firing rate vs external current...")
b()
print("Negative pulse...")
c()
