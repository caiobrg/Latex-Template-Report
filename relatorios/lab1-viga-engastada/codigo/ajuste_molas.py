"""Ajuste das rigidezes adimensionais da base (kappa_R, kappa_L) às duas primeiras
frequências naturais medidas e validação com os modos restantes."""
import numpy as np
from scipy.optimize import least_squares

from modelo import frequencias
from parametros import EI, L


# <<ajuste>>
def ajusta_molas(f_medidas, mu=0.0, chute=(1.0, 3.0)):
    """Resolve f_1(kR, kL) = f_1,med e f_2(kR, kL) = f_2,med.

    As incógnitas são log10(kappa_R) e log10(kappa_L): as rigidezes variam por
    ordens de grandeza e ficam positivas. Os resíduos são erros relativos."""
    alvo = np.asarray(f_medidas[:2], float)

    def residuos(p):
        return frequencias(10 ** p[0], 10 ** p[1], mu, n=2) / alvo - 1

    sol = least_squares(residuos, chute, xtol=1e-12, ftol=1e-12, gtol=1e-12)
    kR, kL = 10 ** sol.x
    return kR, kL, sol


def rigidezes_fisicas(kR, kL):
    """k_R = kappa_R EI / L [N m/rad] e k_L = kappa_L EI / L^3 [N/m]."""
    return kR * EI / L, kL * EI / L**3


def valida(f_medidas, kR, kL, mu=0.0):
    """Frequências previstas e erro relativo para todos os modos medidos."""
    f_med = np.asarray(f_medidas, float)
    f_mod = frequencias(kR, kL, mu, n=len(f_med))
    return f_mod, (f_mod - f_med) / f_med
# <<fim_ajuste>>


def sensibilidade(f_medidas, delta=0.01, mu=0.0):
    """Refaz o ajuste com f_1 ou f_2 perturbada em +-delta e recalcula os erros
    em relação às frequências medidas originais."""
    casos = []
    for i in (0, 1):
        for sinal in (-1, 1):
            alvo = np.array(f_medidas[:2], float)
            alvo[i] *= 1 + sinal * delta
            kR, kL, _ = ajusta_molas(alvo, mu)
            casos.append((i + 1, sinal * delta, kR, kL, valida(f_medidas, kR, kL, mu)[1]))
    return casos


def mapa_frequencias(log_kR, log_kL, mu=0.0):
    """f_1 e f_2 numa grade de (log10 kappa_R, log10 kappa_L), para a figura."""
    F1 = np.empty((len(log_kL), len(log_kR)))
    F2 = np.empty_like(F1)
    for i, a in enumerate(log_kL):
        for j, b in enumerate(log_kR):
            F1[i, j], F2[i, j] = frequencias(10**b, 10**a, mu, n=2)
    return F1, F2
