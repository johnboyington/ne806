import numpy as np
import matplotlib.pyplot as plt

Sig_c, Sig_s, Sig_f = 0.011437, 0.05, 0.013
Sig_t = Sig_c + Sig_s + Sig_f


# diffusion
def y(x, B):
    return np.pi/(2*B)


x = np.linspace(-1, 1, 64)
z = 0.7104 / Sig_s

Thickness = 10 * (1 / Sig_t)
T = np.linspace(-Thickness, Thickness, 64)
B = np.pi/(2*(Thickness+z))
diff10 = np.zeros(len(T))
for i, x_i in enumerate(T):
    diff10[i] = np.sin(B * (x_i + y(x_i, B)))

Thickness = 0.5 * (1 / Sig_t)
T = np.linspace(-Thickness, Thickness, 64)
B = np.pi/(2*(Thickness+z))
diff05 = np.zeros(len(T))
for i, x_i in enumerate(T):
    diff05[i] = np.sin(B * (x_i + y(x_i, B)))

diff05 = diff05 / np.sum(diff05)
diff10 = diff10 / np.sum(diff10)


# power iteration
phi05 = np.load('p2_05.npy')
phi10 = np.load('p2_10.npy')

phi05 = phi05 / np.sum(phi05)
phi10 = phi10 / np.sum(phi10)


# 5 mfp
fig = plt.figure(0)
ax = fig.add_subplot(111)
ax.set_xlabel('x (cm)')
ax.set_ylabel('Flux Profile')
ax.plot(x, phi05, 'k', label='Power Iteration')
ax.plot(x, diff05, 'k--', label='Diffusion Approximation')

ax.legend()
fig.savefig('p4_05.png', dpi=250)


# 10 mfp
fig = plt.figure(1)
ax = fig.add_subplot(111)
ax.set_xlabel('x (cm)')
ax.set_ylabel('Flux Profile')
ax.plot(x, phi10, 'k', label='Power Iteration')
ax.plot(x, diff10, 'k--', label='Diffusion Approximation')

ax.legend()
fig.savefig('p4_10.png', dpi=250)
