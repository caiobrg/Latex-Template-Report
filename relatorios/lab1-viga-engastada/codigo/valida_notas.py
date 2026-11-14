"""Confronta as expressões das notas manuscritas com as deduções corrigidas
(ver errata_notas.md). Não entra no relatório como listing."""
import numpy as np
import sympy as sp
from scipy.optimize import brentq
from solucao_geral import C1, C2, C3, C4
from contorno_ponta import lam, mu
from caso_engastado import Amat, car_eng
from caso_molas import car_mol, aL, aR

s, c, sh, ch = sp.sin(lam), sp.cos(lam), sp.sinh(lam), sp.cosh(lam)
trig2 = {s**2: 1 - c**2, sh**2: ch**2 - 1}

print("== Caso engastado ==")
# Matriz das notas (p. 3): A12 = cos - cosh + mu lam (sin - sinh); 2a linha vinda de (4b)
A_notas = sp.Matrix([[s - sh + mu * lam * (c - ch), c - ch + mu * lam * (s - sh)],
                     [c + ch, s + sh]])
A_certa = sp.Matrix([Amat.row(1), -Amat.row(0)])  # mesma ordem de linhas das notas
print("A_notas - A_certa =", sp.simplify(A_notas - A_certa))
print("det(A_notas) - det(A_certa) =",
      sp.simplify(sp.expand(A_notas.det() - A_certa.det()).subs(trig2)))
eq_notas = mu * lam * (c * sh - s * c) + c * ch + 1           # como escrita (p. 3)
print("eq. final das notas - correta =", sp.simplify(eq_notas - car_eng))

print("\n== Caso com molas ==")
# (14a)-(14b) das notas: alpha = k/EI e sinal perdido em -EI W'''(0) = k_L W(0)
sol = sp.solve([C3 - (aR * (C2 + C4) + C1), C4 - (aL * (C1 + C3) + C2)], [C3, C4])
print("(15) das notas reproduzida:", {k: sp.factor(v) for k, v in sol.items()})
M4 = sp.symbols("M4")  # m beta^4/(rho A), como nas notas (p. 8)
eq_notas_m = (aL * aR - 1 + (aR - aL) * s * ch - (aL + aR) * c * sh + (aL * aR + 1) * c * ch
              + M4 * ((-1 - aL * aR) * s * ch + (1 + aL * aR) * c * sh
                      - 2 * aL * s * sh + 2 * aR * c * ch))
# Desfazendo os dois erros (k_L -> -k_L, M4 -> mu lam) a eq. das notas vira -car_mol,
# exceto pelo termo em cos*sinh (erro de sinal na expansão do determinante):
dif = sp.expand(eq_notas_m.subs(aL, -aL).subs(M4, mu * lam) + car_mol)
print("notas(k_L -> -k_L, M4 -> mu lam) + correta =", sp.factor(dif))

print("\n== lambda_r com kappa_R = 1e6, mu = 0 (engastado: 1.8751, 4.6941) ==")
kR, kL = sp.symbols("kappa_R kappa_L", positive=True)


def raizes(expr, kLv):
    e = (expr * lam**4).subs({aR: kR / lam, aL: kL / lam**3})
    f = sp.lambdify(lam, e.subs({kR: 1e6, kL: kLv, mu: 0, M4: 0}))
    g = np.linspace(0.05, 6, 6000)
    v = f(g)
    i = np.nonzero(np.sign(v[:-1]) != np.sign(v[1:]))[0]
    return [round(brentq(f, g[j], g[j + 1]), 4) for j in i[:3]]


for kLv in (1, 10, 100, 1e3, 1e5):
    print(f"kappa_L = {kLv:>8g}: correta {raizes(car_mol, kLv)}"
          f"  notas {raizes(eq_notas_m, kLv)}")
