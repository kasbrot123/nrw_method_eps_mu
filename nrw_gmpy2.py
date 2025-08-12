
import numpy as np
from scipy import constants as const
from matplotlib import pyplot as plt
import gmpy2 as g



# float/complex precision
ctx = g.get_context()
ctx.precision = 30

################################################################################
# helper functions

def db(S):
    return 20 * np.log10(np.abs(S))
    return [20*g.log10(g.sqrt(g.norm(S[i]))) for i in range(len(S))]

def real(S):
    return [S[i].real for i in range(len(S))]

def imag(S):
    return [S[i].imag for i in range(len(S))]

def norm(S):
    return [float(g.sqrt(g.norm(S[i]))) for i in range(len(S))]

def unwrap(Phase):
    Phase_diff = [Phase[i+1]-Phase[i] for i in range(len(Phase)-1)]
    Phase_unwrap = [Phase[0]]
    add = 0
    for i in range(len(Phase_diff)):
        if g.sqrt((Phase_diff[i])**2) > g.const_pi()/2:
            if Phase_diff[i] > 0:
                add -= 2*g.const_pi()
            else:
                add += 2*g.const_pi()

        Phase_unwrap.append(Phase[i+1] + add)
    return Phase_unwrap

########################################################


def reverse_nrw2(f, Length, eps, mu, lambda_cut=g.inf()):

    # check if every input is ok
    if False:
        raise Exception('Input Error')

    f = [g.mpc(i) for i in f]
    lambda_cut = g.mpc(lambda_cut)
    mu = [g.mpc(i) for i in mu]
    eps = [g.mpc(i) for i in eps]
    lambda_free = [const.c / i for i in f]
    Length = g.mpc(Length)

    T_minus = [g.exp(-g.sqrt(1/lambda_cut**2-eps[i]*mu[i]/lambda_free[i]**2)*2*g.const_pi()*Length) for i in range(len(lambda_free))]
    T_plus = [g.exp(+g.sqrt(1/lambda_cut**2-eps[i]*mu[i]/lambda_free[i]**2)*2*g.const_pi()*Length) for i in range(len(lambda_free))]

    T = []
    for i in range(len(T_minus)):
        # first plus value, then minus
        # minus first did not work, also T = 1 is plus
        if g.sqrt(g.norm(T_plus[i])) <= 1:
            T.append(T_plus[i])
        else:
            T.append(T_minus[i])

    Lambda = [1/g.sqrt(eps[i]*mu[i]/lambda_free[i]**2-1/lambda_cut**2) for i in range(len(lambda_free))]

    Gamma_plus  = [(mu[i] * Lambda[i] * g.sqrt(1/lambda_free[i]**2-1/lambda_cut**2) - 1) / (mu[i] * Lambda[i] * g.sqrt(1/lambda_free[i]**2-1/lambda_cut**2) + 1) for i in range(len(Lambda))]
    Gamma_minus  = [(mu[i] * -Lambda[i] * g.sqrt(1/lambda_free[i]**2-1/lambda_cut**2) - 1) / (mu[i] * Lambda[i] * g.sqrt(1/lambda_free[i]**2-1/lambda_cut**2) + 1) for i in range(len(Lambda))]
    Gamma = []
    for i in range(len(Gamma_minus)):
        if g.sqrt(g.norm(Gamma_plus[i])) <= 1:
            Gamma.append(Gamma_plus[i])
        else:
            Gamma.append(Gamma_minus[i])

    S11 = [Gamma[i] * (1 - T[i]**2) / (1-Gamma[i]**2*T[i]**2) for i in range(len(Gamma))]
    S21 = [T[i] * (1 - Gamma[i]**2) / (1-Gamma[i]**2*T[i]**2) for i in range(len(Gamma))]

    S11 = np.array(S11, dtype=np.complex128)
    S21 = np.array(S21, dtype=np.complex128)

    return S11, S21





################################################################################
# START OF NRW


def nrw2(f, Length, S11, S21, lambda_cut=g.inf()):

    f = [g.mpc(i) for i in f]
    lambda_cut = g.mpc(lambda_cut)
    Length = g.mpc(Length)
    S11 = [g.mpc(i) for i in S11]
    S21 = [g.mpc(i) for i in S21]


    X_nrw = [(1 - (S21[i]**2 - S11[i]**2))/(2 * S11[i]) for i in range(len(S11))]

    lambda_free_nrw = [const.c / f[i] for i in range(len(f))]
    Gamma_plus_nrw = [X_nrw[i] + (X_nrw[i]**2 - 1)**0.5 for i in range(len(X_nrw))]
    Gamma_minus_nrw = [X_nrw[i] - (X_nrw[i]**2 - 1)**0.5 for i in range(len(X_nrw))]
    Gamma_nrw = []
    for i in range(len(Gamma_minus_nrw)):
        if g.sqrt(g.norm(Gamma_minus_nrw[i])) < 1:
            Gamma_nrw.append(Gamma_minus_nrw[i])
        else:
            Gamma_nrw.append(Gamma_plus_nrw[i])


    T_nrw = [(S11[i] + S21[i] - Gamma_nrw[i])/(1 - (S11[i] + S21[i])*Gamma_nrw[i]) for i in range(len(S11))]
    # T.real = -T.real # reverse engineering, else did not work
    # because the flow in comsol was in the wrong direction


    Phase_nrw = [g.phase(T_nrw[i]) for i in range(len(T_nrw))]
    Phase_unwrap_nrw = unwrap(Phase_nrw)


    Dphi_df_nrw = [(Phase_unwrap_nrw[i+1]-Phase_unwrap_nrw[i])/(f[i+1]-f[i]) for i in range(len(f)-1)]
    Dphi_df_nrw.append(Dphi_df_nrw[-1])
    Tau_meas_nrw = [-1/(2*g.const_pi()) * Dphi_df_nrw[i] for i in range(len(Dphi_df_nrw))]

    T_log_nrw = [g.log(1/T_nrw[i]) for i in range(len(T_nrw))]
    N_freq = len(T_log_nrw)
    N_matrix = []
    N_ambig = 10

    for i in range(-N_ambig,N_ambig+1,1):
        Helper = []
        for j in range(N_freq):
            Helper.append(i)
        N_matrix.append(Helper)

    T_log_imag_nrw = [T_log_nrw[i].imag for i in range(len(T_log_nrw))]
    T_log_imag_unwrap_nrw = unwrap(T_log_imag_nrw)


    AMBIGUOUS = []
    for i in range(2*N_ambig+1):
        Helper = []
        for j in range(N_freq):
            Helper.append(g.mpc(T_log_nrw[j].real + g.mpc(1j*T_log_imag_unwrap_nrw[j]) + 1j*2*g.const_pi()*N_matrix[i][j]))
        AMBIGUOUS.append(Helper)


    LAMBDA_SQUARE_INV = []
    for i in range(2*N_ambig+1):
        Helper = []
        for j in range(N_freq):
            Helper.append(-(1/2/g.const_pi()/Length*AMBIGUOUS[i][j])**2)
        LAMBDA_SQUARE_INV.append(Helper)

    LAMBDA = []
    for i in range(2*N_ambig+1):
        Helper = []
        for j in range(N_freq):
            # this is not doing anything
            value = g.sqrt(LAMBDA_SQUARE_INV[i][j])
            if value.real > 0:
                Helper.append(1/value)
            else:
                # print(j)
                # print('lambda not positive')
                Helper.append(-1/value)

        LAMBDA.append(Helper)

    GAMMA = []
    for i in range(2*N_ambig+1):
        Helper = []
        for j in range(N_freq):
            Helper.append(Gamma_nrw[j])
        GAMMA.append(Helper)

    LAMBDA_FREE = []
    for i in range(2*N_ambig+1):
        Helper = []
        for j in range(N_freq):
            Helper.append(lambda_free_nrw[j])
        LAMBDA_FREE.append(Helper)


    MU = []
    for i in range(2*N_ambig+1):
        Helper = []
        for j in range(N_freq):
            Helper.append((1 + GAMMA[i][j])/LAMBDA[i][j]/(1 - GAMMA[i][j])/(g.sqrt(1/LAMBDA_FREE[i][j]**2-1/lambda_cut**2)))
        MU.append(Helper)


    EPS = []
    for i in range(2*N_ambig+1):
        Helper = []
        for j in range(N_freq):
            Helper.append(LAMBDA_FREE[i][j]**2 / MU[i][j] * ( 1/lambda_cut**2 + LAMBDA_SQUARE_INV[i][j] ))
        EPS.append(Helper)

    MU_EPS = []
    for i in range(2*N_ambig+1):
        Helper = []
        for j in range(N_freq):
            Helper.append(MU[i][j]*EPS[i][j])
        MU_EPS.append(Helper)

    F = []
    for i in range(2*N_ambig+1):
        Helper = []
        for j in range(N_freq):
            Helper.append(f[j])
        F.append(Helper)

    GRAD = []
    for i in range(2*N_ambig+1):
        Helper = []
        for j in range(N_freq-1):
            Helper.append((MU_EPS[i][j+1]-MU_EPS[i][j])/(f[j+1]-f[j]))
        Helper.append(Helper[-1])
        GRAD.append(Helper)

    TAU_CALC = []
    for i in range(2*N_ambig+1):
        Helper = []
        for j in range(N_freq):
            # Helper.append(1/const.c**2*(F[i][j]*MU_EPS[i][j] + F[i][j]**2/2*GRAD[i][j])*Length/(g.sqrt(MU_EPS[i][j]*F[i][j]**2/const.c**2-1/lambda_cut**2)))
            Helper.append((1/const.c/LAMBDA_FREE[i][j]*MU_EPS[i][j] + 1/LAMBDA_FREE[i][j]**2/2*GRAD[i][j])*Length/(g.sqrt(MU_EPS[i][j]/LAMBDA_FREE[i][j]**2-1/lambda_cut**2)))
        TAU_CALC.append(Helper)

    MINIMUM = []
    for i in range(2*N_ambig+1):
        Helper = []
        for j in range(N_freq):
            value = TAU_CALC[i][j].real + Tau_meas_nrw[j].real
            # abs value
            if value > 0:
                Helper.append(value)
            else:
                Helper.append(-value)
        MINIMUM.append(Helper)

    Argmin = []
    for j in range(N_freq):
        Helper = []
        for i in range(2*N_ambig+1):
            Helper.append(MINIMUM[i][j])
        Argmin.append(Helper.index(min(Helper)))
        # Argmin.append(7)


    mu_nrw = []
    eps_nrw = []
    n_best_list = []
    for j in range(N_freq):
        a = Argmin[j]
        n_best_list.append(N_matrix[a][j])
        mu_nrw.append(MU[a][j])
        eps_nrw.append(EPS[a][j])


    mu_nrw = np.array(mu_nrw, dtype=np.complex128)
    eps_nrw = np.array(eps_nrw, dtype=np.complex128)

    return eps_nrw, mu_nrw

################################################################################
# Create noise

def noise(S, A=0.1, PHI=0.01):
    """
    typical errors for VNA measurement:
        ampltiude error < +-0.5 dB (+-5%)
        phase error < 1°

    A = 0.1 # dB
    PHI = 0.01 # Grad

    """
    a = 10**(2*A*(np.random.rand(len(S)) - 0.5) / 20)
    phi = 2*np.pi/180*(2*PHI*(np.random.rand(len(S))-0.5))

    return a*S*np.exp(1j*phi)




################################################################################


if __name__ == '__main__':

    print('Test example:')

    f = np.linspace(60e9, 90e9, 100)
    eps = np.ones(f.shape)*(8.5+0.2j)
    mu = np.ones(f.shape)
    Length = 0.01 # m
    lambda_cut = const.c / 48.373e9 # wr12 waveguide
    lambda_cut = 1e100


    # reverse calculation
    S11, S21 = reverse_nrw2(f, Length, eps, mu, lambda_cut)

    # noise
    S11 = noise(S11, 0.1, 0.01)
    S21 = noise(S21, 0.1, 0.01)

    # nrw calculation
    eps_nrw, mu_nrw = nrw2(f, Length, S11, S21, lambda_cut)



    plt.figure()
    plt.plot(f, eps_nrw.real, label='eps_nrw real')
    plt.plot(f, eps_nrw.imag, label='eps_nrw imag')
    plt.plot(f, eps.real, '--', label='eps_start real')
    plt.plot(f, eps.imag, '--', label='eps_start imag')
    plt.xlabel('frequency [Hz]')
    plt.ylabel('Permittivity [1]')
    plt.title('Permittivity')
    plt.legend()
    plt.tight_layout()

    plt.figure()
    plt.plot(f, mu_nrw.real, label='mu_nrw real')
    plt.plot(f, mu_nrw.imag, label='mu_nrw imag')
    plt.plot(f, mu.real, '--', label='mu_start real')
    plt.plot(f, mu.imag, '--', label='mu_start imag')
    plt.xlabel('frequency [Hz]')
    plt.ylabel('Permeability [1]')
    plt.title('Permeability')
    plt.legend()
    plt.tight_layout()

    fig, ax = plt.subplots(2)
    ax[0].plot(f, db(S11), label='|S11|')
    ax[0].plot(f, db(S21), label='|S21|')
    ax[1].plot(f, np.unwrap(np.angle(S11)), label='angle(S11)')
    ax[1].plot(f, np.unwrap(np.angle(S21)), label='angle(S21)')
    ax[0].set_xlabel('frequency [Hz]')
    ax[0].set_ylabel('Amplitude [dB]')
    ax[1].set_xlabel('frequency [Hz]')
    ax[1].set_ylabel('Phase [rad]')
    fig.tight_layout()


    plt.show()



