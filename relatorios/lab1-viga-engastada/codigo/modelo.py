"""Modelo analítico numérico: equações características, frequências naturais, modos
normalizados e FRF (acelerância) para a base engastada e para a base elástica.
Notação da Seção 1: lambda = beta L, mu = m/(rho A L), kappa_R = k_R L/EI,
kappa_L = k_L L^3/EI. Base engastada: kappa_R = kappa_L = inf."""
import numpy as np
from scipy.optimize import brentq

from parametros import L, rhoA, lambda_para_freq

INF = np.inf


# <<eq_caracteristica>>
def eq_caracteristica(lam, kR=INF, kL=INF, mu=0.0):
    """Equação característica dividida por cosh(lambda) (e, nas molas, por
    1 + alpha_L alpha_R > 0), o que não muda as raízes e evita overflow."""
    s, c, ch = np.sin(lam), np.cos(lam), np.cosh(lam)
    th = np.tanh(lam)
    if np.isinf(kR) and np.isinf(kL):               # engaste, Eq. (car_engaste)
        return 1 / ch + c + mu * lam * (c * th - s)
    aR, aL = kR / lam, kL / lam**3                  # alpha_R, alpha_L
    aa = aR * aL
    D = (1 + aa) / ch + (aa - 1) * c + (aL - aR) * c * th - (aL + aR) * s \
        + mu * lam * ((1 - aa) * (s - c * th) - 2 * aL * s * th - 2 * aR * c)
    return D / (1 + aa)                             # Eq. (car_molas) / [(1+aa) cosh]


def raizes(kR=INF, kL=INF, mu=0.0, n=8, passo=1e-3):
    """As n primeiras raízes lambda_r: varredura de sinais + bisseção (Brent)."""
    grade = np.arange(passo, np.pi * (n + 1), passo)
    v = eq_caracteristica(grade, kR, kL, mu)
    trocas = np.nonzero(np.sign(v[:-1]) != np.sign(v[1:]))[0][:n]
    return np.array([brentq(eq_caracteristica, grade[i], grade[i + 1],
                            args=(kR, kL, mu), xtol=1e-14) for i in trocas])
# <<fim_eq_caracteristica>>


def frequencias(kR=INF, kL=INF, mu=0.0, n=8):
    """Frequências naturais f_r [Hz]."""
    return lambda_para_freq(raizes(kR, kL, mu, n))


# <<modos>>
def forma_modal(xi, lam, kR=INF, kL=INF):
    """F_r(xi), xi = x/L, não normalizada: Eqs. (F_engaste) e (F_molas).
    No caso com molas, phi_1 e phi_2 são divididas por 1 + alpha_L alpha_R."""
    z = lam * np.asarray(xi)
    if np.isinf(kR) and np.isinf(kL):
        gama = -(np.cos(lam) + np.cosh(lam)) / (np.sin(lam) + np.sinh(lam))
        return np.cos(z) - np.cosh(z) + gama * (np.sin(z) - np.sinh(z))
    aR, aL = kR / lam, kL / lam**3
    p, q = 1.0, (1 - aR * aL) / (1 + aR * aL)
    gL, gR = 2 * aL / (1 + aR * aL), 2 * aR / (1 + aR * aL)
    B11 = -p * np.cos(lam) + q * np.cosh(lam) - gL * np.sinh(lam)
    B12 = -p * np.sin(lam) + q * np.sinh(lam) + gR * np.cosh(lam)
    phi1 = p * np.cos(z) + q * np.cosh(z) - gL * np.sinh(z)
    phi2 = p * np.sin(z) + q * np.sinh(z) + gR * np.cosh(z)
    return phi1 - B11 / B12 * phi2


def modos_normalizados(x, lams, kR=INF, kL=INF, mu=0.0):
    """W_r(x) normalizados pela massa, Eq. (Wr_norm); uma linha por modo."""
    xi = np.linspace(0, 1, 4001)
    W = []
    for lam in lams:
        F = forma_modal(xi, lam, kR, kL)
        massa = rhoA * L * (np.trapezoid(F**2, xi) + mu * F[-1] ** 2)
        W.append(forma_modal(np.asarray(x) / L, lam, kR, kL) / np.sqrt(massa))
    return np.array(W)
# <<fim_modos>>


# <<frf>>
def acelerancia(freq, x_e, x_i, zetas, kR=INF, kL=INF, mu=0.0, n=8):
    """Acelerância A = -Omega^2 sum_r W_r(x_e) W_r(x_i) /
    (omega_r^2 - Omega^2 + 2 i zeta_r omega_r Omega), Eqs. (H) e (acelerancia)."""
    lams = raizes(kR, kL, mu, n)
    w_r = 2 * np.pi * lambda_para_freq(lams)
    W = modos_normalizados([x_e, x_i], lams, kR, kL, mu)
    Om = 2 * np.pi * np.asarray(freq)[:, None]
    H = W[:, 0] * W[:, 1] / (w_r**2 - Om**2 + 2j * np.asarray(zetas) * w_r * Om)
    return -Om[:, 0] ** 2 * H.sum(axis=1)
# <<fim_frf>>
