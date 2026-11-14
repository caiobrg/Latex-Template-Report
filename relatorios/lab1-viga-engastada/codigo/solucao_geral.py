"""Separação de variáveis e solução geral da viga de Euler-Bernoulli."""
import sympy as sp

x, t, L = sp.symbols("x t L", positive=True)
E, I, rho, A = sp.symbols("E I rho A", positive=True)
beta, omega = sp.symbols("beta omega", positive=True)
C1, C2, C3, C4 = sp.symbols("C1:5")

# --- Separação de variáveis: w(x, t) = W(x) g(t) em EI w'''' + rho A w_tt = 0
Wf, g = sp.Function("W"), sp.Function("g")
edp = E * I * sp.diff(Wf(x) * g(t), x, 4) + rho * A * sp.diff(Wf(x) * g(t), t, 2)
separada = sp.expand(edp / (rho * A * Wf(x) * g(t)))  # (EI/rhoA) W''''/W + g''/g

# Cada parcela é igual a uma constante (omega^2) -> duas EDOs independentes
sol_g = sp.dsolve(g(t).diff(t, 2) + omega**2 * g(t)).rhs

# beta^4 = rho A omega^2 / (EI)  ->  W'''' - beta^4 W = 0, com W = e^{s x}
s = sp.symbols("s")
raizes = sp.solve(s**4 - beta**4, s)  # +-beta, +-i beta

# --- Solução geral na forma trigonométrica/hiperbólica, Eq. (2)
W = C1 * sp.cos(beta * x) + C2 * sp.sin(beta * x) \
    + C3 * sp.cosh(beta * x) + C4 * sp.sinh(beta * x)
assert sp.simplify(W.diff(x, 4) - beta**4 * W) == 0


def dW(n, x0):
    """n-ésima derivada de W(x) avaliada em x = x0."""
    return W.diff(x, n).subs(x, x0)


if __name__ == "__main__":
    print("EDP / (rho A W g) =", separada)
    print("g(t) =", sol_g)
    print("raízes de s^4 - beta^4:", raizes)
    for n in range(4):
        print(f"W^({n})(x) =", W.diff(x, n))
