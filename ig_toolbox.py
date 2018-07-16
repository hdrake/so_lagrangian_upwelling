import numpy as np
######################################

# Inverse Gaussian (scaled)
def inversegaussian(x, Gamma, Delta, Scaling):
    inversegaussian = []
    for i in range(x.size):
        inversegaussian += [Scaling*(((Gamma**3)/(4 * np.pi * (Delta**2) * (x[i]**3)))**(0.5))
                            * np.exp(((-Gamma)*((x[i]-Gamma)**2))/(4*(Delta**2)*x[i]))]
    return np.array(inversegaussian)

######################################
# Residual Calculation for least squares
def res(p, y, x):
    Gam, Del, Sca = p
    y_fit = inversegaussian(x, Gam, Del, Sca)
    err = y - y_fit
    return err

######################################
