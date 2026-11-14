"""Caso 1: base engastada + ponta livre com massa concentrada."""
import sympy as sp
from scipy.optimize import brentq
from solucao_geral import x, beta, C1, C2, C3, C4, W, dW
from contorno_ponta import lam, mu, bc_M, bc_Q

s, c, sh, ch = sp.sin(lam), sp.cos(lam), sp.sinh(lam), sp.cosh(lam)

# --- Engaste: W(0) = 0 e W'(0) = 0  ->  C3 = -C1, C4 = -C2
base = sp.solve([dW(0, 0), dW(1, 0)], [C3, C4])

# --- Sistema A [C1, C2]^T = 0 a partir das condições na ponta
eqs = [sp.expand(e.subs(base)) for e in (bc_M, bc_Q)]
Amat = sp.Matrix([[e.coeff(Ck) for Ck in (C1, C2)] for e in eqs])

# --- Equação característica: det(A) = 0
det = sp.expand(Amat.det()).subs({s**2: 1 - c**2, sh**2: ch**2 - 1})
car_eng = sp.expand(det / 2)
esperado = 1 + c * ch + mu * lam * (c * sh - s * ch)
assert sp.simplify(car_eng - esperado) == 0

# --- Forma modal não normalizada (C2/C1 tirado da 1a linha, W''(L) = 0)
razao = sp.simplify(-Amat[0, 0] / Amat[0, 1])  # -(cos + cosh)/(sin + sinh)
xi = sp.symbols("xi", positive=True)           # xi = x/L  ->  beta x = lambda xi
F_eng = (W.subs(base).subs(C2, razao * C1) / C1).subs(beta * x, lam * xi)
F_eng = sp.simplify(F_eng)

if __name__ == "__main__":
    print("C3, C4 =", base)
    sp.pprint(Amat)
    print("eq. característica:", esperado, "= 0")
    print("C2/C1 =", razao)
    print("F(xi) =", F_eng)
    # Verificação: sem massa na ponta (mu = 0) -> 1.8751, 4.6941, 7.8548
    f = sp.lambdify(lam, car_eng.subs(mu, 0))
    print("lambda_r (mu = 0):",
          [round(brentq(f, a, a + 1.5), 4) for a in (1.0, 4.0, 7.0)])
