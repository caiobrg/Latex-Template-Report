"""Equação modal amortecida, resposta harmônica e FRF analítica da viga."""
import sympy as sp

t, x, L = sp.symbols("t x L", positive=True)
x_i, x_e = sp.symbols("x_i x_e", positive=True)   # entrada (força) e saída (resposta)
Om, w_r, z_r, F0 = sp.symbols("Omega omega_r zeta_r F_0", positive=True)
r, N = sp.symbols("r N", integer=True, positive=True)
Wr = sp.Function("W_r")
eta = sp.Function("eta")

# --- Força generalizada de uma carga pontual f = F0 delta(x - x_i) e^{i Omega t}
N_r = sp.integrate(sp.DiracDelta(x - x_i) * Wr(x), (x, 0, L))
N_r = N_r.subs(sp.Heaviside(L - x_i), 1)  # 0 < x_i < L
print("N_r / (F0 e^{i Omega t}) =", N_r)  # W_r(x_i)

# --- Equação modal com amortecimento modal zeta_r
edo = eta(t).diff(t, 2) + 2 * z_r * w_r * eta(t).diff(t) + w_r**2 * eta(t) \
    - N_r * F0 * sp.exp(sp.I * Om * t)

# Regime permanente: eta = Y e^{i Omega t}
Y = sp.symbols("Y")
Ysol = sp.solve(sp.simplify(edo.subs(eta(t), Y * sp.exp(sp.I * Om * t)).doit()
                            / sp.exp(sp.I * Om * t)), Y)[0]
G_r = 1 / (1 - (Om / w_r)**2 + 2 * sp.I * z_r * Om / w_r)
assert sp.simplify(Ysol - N_r * F0 / w_r**2 * G_r) == 0

# --- Receptância H = w(x_e)/f(x_i): w(x_e) = sum W_r(x_e) eta_r, truncada em N modos
termo = Wr(x_e) * Ysol / F0
assert sp.simplify(termo - Wr(x_e) * Wr(x_i)
                   / (w_r**2 - Om**2 + 2 * sp.I * z_r * w_r * Om)) == 0
H = sp.Sum(termo, (r, 1, N))

# --- Mobilidade e acelerância (o ensaio mede aceleração): d/dt -> i Omega
Y_mob = sp.I * Om * H
A_acc = -Om**2 * H

# --- Módulo e fase de um termo modal
mod_G = sp.simplify(sp.Abs(G_r))
fase_G = -sp.atan2(2 * z_r * Om / w_r, 1 - (Om / w_r)**2)
assert sp.simplify(sp.tan(fase_G) - sp.im(G_r) / sp.re(G_r)) == 0

if __name__ == "__main__":
    print("eta_r,p / e^{i Omega t} =", sp.simplify(Ysol))
    print("H(x_e, x_i, Omega) ="); sp.pprint(H)
    print("A(x_e, x_i, Omega) ="); sp.pprint(A_acc)
    print("|G_r| =", mod_G)
    print("fase G_r =", sp.simplify(fase_G))
