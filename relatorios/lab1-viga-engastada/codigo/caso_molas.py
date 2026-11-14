"""Caso 2: base com molas torcional (k_R) e transversal (k_L) + ponta com massa."""
import sympy as sp
from solucao_geral import x, E, I, beta, C1, C2, C3, C4, W, dW
from contorno_ponta import lam, mu, bc_M, bc_Q

s, c, sh, ch = sp.sin(lam), sp.cos(lam), sp.sinh(lam), sp.cosh(lam)
kR, kL = sp.symbols("k_R k_L", positive=True)
aR, aL = sp.symbols("alpha_R alpha_L", positive=True)

# --- Base elástica:  EI W''(0) = k_R W'(0)   e   -EI W'''(0) = k_L W(0)
# alpha_R = k_R/(EI beta) = kappa_R/lam,  alpha_L = k_L/(EI beta^3) = kappa_L/lam^3
bc_R = sp.expand((E * I * dW(2, 0) - kR * dW(1, 0)) / (E * I * beta**2))
bc_L = sp.expand((-E * I * dW(3, 0) - kL * dW(0, 0)) / (E * I * beta**3))
bc_R, bc_L = (sp.expand(e.subs({kR: aR * E * I * beta, kL: aL * E * I * beta**3}))
              for e in (bc_R, bc_L))
base = {k: sp.factor(v) for k, v in sp.solve([bc_R, bc_L], [C3, C4]).items()}
D = 1 + aL * aR  # denominador comum (sempre > 0)

# --- Sistema B [C1, C2]^T = 0 (multiplicado por D)
eqs = [sp.expand(sp.cancel(e.subs(base) * D)) for e in (bc_M, bc_Q)]
Bmat = sp.Matrix([[sp.collect(e.coeff(Ck), mu) for Ck in (C1, C2)] for e in eqs])

# --- Equação característica: det(B) / (2 D) = 0
det = sp.expand(sp.expand(Bmat.det()).subs({s**2: 1 - c**2, sh**2: ch**2 - 1}))
car_mol = sp.expand(sp.cancel(det / (2 * D)))
esperado = (1 + aL * aR + (aL * aR - 1) * c * ch + (aL - aR) * c * sh
            - (aL + aR) * s * ch
            + mu * lam * ((1 - aL * aR) * (s * ch - c * sh)
                          - 2 * aL * s * sh - 2 * aR * c * ch))
assert sp.expand(car_mol - esperado) == 0

# --- Limite de molas infinitamente rígidas -> volta ao caso engastado
eps = sp.symbols("epsilon", positive=True)
limite = sp.limit(car_mol.subs({aL: 1 / eps, aR: 1 / eps}) * eps**2, eps, 0)
assert sp.simplify(limite - (1 + c * ch + mu * lam * (c * sh - s * ch))) == 0

# --- Forma modal não normalizada (C2/C1 da 1a linha, W''(L) = 0)
razao = sp.simplify(-Bmat[0, 0] / Bmat[0, 1])
xi = sp.symbols("xi", positive=True)
F_mol = sp.expand(W.subs(base) * D / C1).subs(C2, razao * C1).subs(beta * x, lam * xi)

if __name__ == "__main__":
    print("C3 =", base[C3], "\nC4 =", base[C4])
    for i, j in [(0, 0), (0, 1), (1, 0), (1, 1)]:
        print(f"B[{i + 1},{j + 1}] =", Bmat[i, j])
    print("Delta =", sp.collect(car_mol, [mu, c * ch, s * ch, c * sh]))
    print("C2/C1 =", razao)
