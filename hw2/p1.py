import numpy as np

# (a)

I = 3e6
flux_density = 4 * I
print('(a) Flux Density:  {:4.2e}'.format(flux_density))


# (b)

x = I * (1 + 0.5)
y = I * (-1 + 0.7071)
z = I * (1 + 0.5)
print('(b) Current Density:  {:4.2e} i_x + {:4.2e} i_y + {:4.2e} i_z'.format(x, y, z))

# (c)
current_density_magnitude = np.sqrt(x**2 + y**2 + z**2)
print('(c) Current Density Magnitude:  {:4.2e}'.format(current_density_magnitude))

# (d)
print('(d) Neutrons Crossing y-direction:  {:4.2e}'.format(abs(y)))
