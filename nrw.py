
import numpy as np
from matplotlib import pyplot as plt
from scipy import constants as const



DTYPE = np.complex128



def nrw(S11, S21, f, d, lambda_cut, n_best=None):

    f = f.astype(DTYPE)
    d = DTYPE(d)
    S11 = S11.astype(DTYPE)
    S21 = S21.astype(DTYPE)


    X = (1 - (S21**2 - S11**2))/(2 * S11)

    lambda_free = const.c / f
    Gamma = X + (X**2 - 1)**0.5
    L = np.abs(Gamma) > 1
    Gamma[L] = X[L] - (X[L]**2 - 1)**0.5


    T = (S11 + S21 - Gamma)/(1 - (S11 + S21)*Gamma)
    # T.real = -T.real # reverse engineering, else did not work
    # because the flow in comsol was in the wrong direction

    # do not use S21 for phase calculation
    phase = np.unwrap(np.angle(T))
    dphi_df = np.gradient(phase, f) # should have same length
    tau_meas = -1/(2*np.pi) * dphi_df


    T_term = np.log(1/T)
    N_freq = len(T_term)
    if n_best is None:
        N_ambig = 10 # number of wavelength ambiguities
        N_matrix = np.ones((2*N_ambig+1, N_freq)) * np.arange(-N_ambig,N_ambig+1,1).reshape(2*N_ambig+1,1)
    else:
        N_ambig = 0
        N_matrix = np.ones((1, N_freq))*n_best

    Ambiguous = T_term.real.reshape(1,N_freq) + 1j*np.unwrap(T_term.imag).reshape(1,N_freq) + 1j*2*np.pi*N_matrix

    LAMBDA_SQUARE_INV = -(1/2/np.pi/d*Ambiguous)**2
    LAMBDA = 1 / ( LAMBDA_SQUARE_INV )**0.5

    GAMMA = np.ones((2*N_ambig+1, N_freq)) * Gamma.reshape(1, N_freq)
    LAMBDA_FREE = np.ones((2*N_ambig+1, N_freq)) * lambda_free.reshape(1, N_freq)

    MU = (1 + GAMMA)/LAMBDA/(1-GAMMA)/(np.sqrt(1/LAMBDA_FREE**2-1/lambda_cut**2))
    EPS = LAMBDA_FREE**2 / MU * ( 1/lambda_cut**2 + LAMBDA_SQUARE_INV )
    MU_EPS = MU*EPS

    F = np.ones((2*N_ambig+1, N_freq))*f.reshape(1,N_freq)
    # ignore the fact that you differentiate complex values
    GRAD = np.gradient(MU_EPS, f, axis=1)

    # either f or lambda
    # tau_calc = 1/const.c**2*(F*MU_EPS + F**2/2*GRAD)*d/(np.sqrt(MU_EPS*F**2/const.c**2-1/lambda_cut**2))
    tau_calc = (1/(const.c*LAMBDA_FREE)*MU_EPS + 1/LAMBDA_FREE**2/2*GRAD)*d/(np.sqrt(MU_EPS/LAMBDA_FREE**2-1/lambda_cut**2))

    Minimum = np.abs(tau_calc.real + tau_meas.real.reshape(1,N_freq))
    Argmin = Minimum.argmin(axis=0)

    # print('best n values')
    # print(N_matrix[Argmin,np.arange(N_freq)])

    mu = MU[Argmin,np.arange(N_freq)]
    eps = EPS[Argmin,np.arange(N_freq)]
    refractive_index = np.sqrt(mu*eps)

    return eps, mu, refractive_index





################################################################################
# start reverse calculation


def nrw_reverse(f, Length, eps, mu, lambda_cut):

    mu = np.ones(len(f),dtype=DTYPE)*mu
    eps = np.ones(len(f),dtype=DTYPE)*eps
    lambda_cut = DTYPE(lambda_cut)

    lambda_free = const.c/f
    T_minus = np.exp(-np.sqrt(1/lambda_cut**2-eps*mu/lambda_free**2)*2*np.pi*Length)
    T_plus = np.exp(+np.sqrt(1/lambda_cut**2-eps*mu/lambda_free**2)*2*np.pi*Length)

    T = T_plus
    Logical = np.abs(T_plus) <= 1
    T[Logical] = T_minus[Logical]

    # minus für die lösung in der Anleitung
    # plus für alle anderen lösungen mit mu null
    Lambda = -1/np.sqrt(eps*mu/lambda_free**2-1/lambda_cut**2)
    Lambda_inv_square = 1/Lambda**2

    Gamma_plus  = (mu * Lambda * np.sqrt(1/lambda_free**2-1/lambda_cut**2) - 1) / (mu * Lambda * np.sqrt(1/lambda_free**2-1/lambda_cut**2) + 1)
    Gamma_minus  = (mu * (-Lambda) * np.sqrt(1/lambda_free**2-1/lambda_cut**2) - 1) / (mu * (-Lambda) * np.sqrt(1/lambda_free**2-1/lambda_cut**2) + 1)
    Gamma = Gamma_plus

    Logical = np.abs(Gamma) <= 1
    Gamma[Logical] = Gamma_minus[Logical]

    S11 = Gamma * (1 - T**2) / (1-Gamma**2*T**2)
    S21 = T * (1 - Gamma**2) / (1-Gamma**2*T**2)

    return S11, S21


################################################################################

def db(S):
    return 20*np.log10(np.abs(S))


def save_s2p(filename, S11, S21):
    data = np.array([f,S11.real,S11.imag,S21.real,S21.imag,S21.real,S21.imag,S11.real,S11.imag])
    np.savetxt(filename, data.T, header='# Hz S RI R 50')



################################################################################


if __name__ == '__main__':
    f = np.linspace(60e9,90e9,100, dtype=DTYPE)

    mu = 1
    # mu = 1.1
    # mu = (3 + 3j)

    eps = 2.5*(1+0.100j)
    eps = 8.4*(1+0.0500j)
    # eps = 8.4*(1+0.0500j)*(60e9)**0.5/f**0.5

    lambda_cut = const.c / (48.373e9)
    # lambda_cut = 1e100
    Length = 0.010

    S11, S21 = nrw_reverse(f, Length, eps, mu, lambda_cut)

    ## save S11/S21 in .s2p file
    # save_s2p('reverse.s2p', S11, S21)
    # data = np.loadtxt('./reverse.s2p', comments=['!', '#'], dtype=DTYPE)
    # f = data[:,0]
    # S11 = data[:,1] + 1j*data[:,2]
    # S22 = data[:,7] + 1j*data[:,8]

    eps_nrw, mu_nrw, N_nrw = nrw(S11, S21, f, Length, lambda_cut)

    f_ghz = f / 1e9
    eps = eps * np.ones(f.shape)
    mu = mu * np.ones(f.shape)

    f1 = plt.figure()
    plt.plot(f_ghz, db(S11), label='S11')
    plt.plot(f_ghz, db(S21), label='S12')
    plt.xlabel('frequency [GHz]')
    plt.ylabel('S parameter [dB]')
    plt.legend()
    plt.title('s-parameter')

    f2 = plt.figure()
    plt.plot(f_ghz, eps_nrw.real, label='eps real (nrw)')
    plt.plot(f_ghz, eps_nrw.imag, label='eps imag (nrw)')
    plt.plot(f_ghz, eps.real, '--', label='eps real')
    plt.plot(f_ghz, eps.imag, '--', label='eps imag')
    plt.xlabel('frequency [GHz]')
    plt.ylabel('Permittivity [1]')
    plt.legend()
    plt.title('permittivity')

    f3 = plt.figure()
    plt.plot(f_ghz, mu_nrw.real, label='mu real (nrw)')
    plt.plot(f_ghz, mu_nrw.imag, label='mu imag (nrw)')
    plt.plot(f_ghz, mu.real, '--', label='mu real')
    plt.plot(f_ghz, mu.imag, '--', label='mu imag')
    plt.legend(['real', 'imag'])
    plt.ylabel('Permeability [1]')
    plt.title('permeability')

    f4 = plt.figure()
    plt.plot(f_ghz, N_nrw.real)
    plt.plot(f_ghz, N_nrw.imag)
    plt.legend(['real', 'imag'])
    plt.title('refractive index')

    plt.show()


