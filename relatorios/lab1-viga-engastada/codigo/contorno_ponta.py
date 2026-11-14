"""Condições de contorno na extremidade livre com massa concentrada (x = L),
comuns aos dois modelos de base."""
import sympy as sp
from solucao_geral import x, L, E, I, rho, A, beta, omega, W, dW

m = sp.symbols("m", positive=True)       # massa concentrada na ponta
lam, mu = sp.symbols("lambda mu", positive=True)

# Momento nulo:            EI W''(L) = 0
# Cortante = inércia:      EI W'''(L) = -m omega^2 W(L)
ponta_M = E * I * dW(2, L)
ponta_Q = E * I * dW(3, L) + m * omega**2 * dW(0, L)

# omega^2 = EI beta^4 / (rho A)  e  m = mu rho A L  ->  divide por EI beta^2 e EI beta^3
adim = {omega: sp.sqrt(E * I / (rho * A)) * beta**2, m: mu * rho * A * L}
bc_M = sp.expand(ponta_M / (E * I * beta**2))
bc_Q = sp.expand(ponta_Q.subs(adim) / (E * I * beta**3))

# Por fim, beta L = lambda
bc_M, bc_Q = (sp.expand(e.subs(L, lam / beta)) for e in (bc_M, bc_Q))

if __name__ == "__main__":
    print("W''(L)/beta^2           =", bc_M, "= 0")
    print("W'''(L)/beta^3 + mu lam W(L) =", bc_Q, "= 0")
