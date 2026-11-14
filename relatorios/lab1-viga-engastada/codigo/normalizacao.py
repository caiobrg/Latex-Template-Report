"""Normalização pela massa (viga + massa na ponta) e verificação de ortogonalidade."""
import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.optimize import brentq
from contorno_ponta import lam, mu
from caso_engastado import car_eng, F_eng, xi
from caso_molas import car_mol, F_mol, aL, aR

kapR, kapL, rhoAL = sp.symbols("kappa_R kappa_L rhoAL", positive=True)

# W_r(x) = C_r F_r(xi):  rho A L C_r^2 [ int_0^1 F_r^2 dxi + mu F_r(1)^2 ] = 1
F = sp.Function("F")(xi)
C_r = 1 / sp.sqrt(rhoAL * (sp.Integral(F**2, (xi, 0, 1)) + mu * F.subs(xi, 1)**2))


def modos(car, Fr, par, n=3):
    """Raízes lambda_r de car = 0 e modos normalizados phi_r(xi) = sqrt(rhoAL) W_r."""
    f = sp.lambdify(lam, (car * lam**4).subs(par))
    grade = np.linspace(0.3, 15, 3000)
    v = f(grade)
    idx = np.nonzero(np.sign(v[:-1]) != np.sign(v[1:]))[0][:n]
    saida = []
    for i in idx:
        lr = brentq(f, grade[i], grade[i + 1])
        Fx = Fr.subs(par).subs(lam, lr)
        d = [sp.lambdify(xi, Fx.diff(xi, k)) for k in range(3)]
        N = quad(lambda z: d[0](z)**2, 0, 1)[0] + par[mu] * d[0](1)**2
        saida.append((lr, [lambda z, dk=dk, N=N: dk(z) / np.sqrt(N) for dk in d]))
    return saida


def matrizes(ms, par):
    """Massa: int phi_r phi_s + mu phi_r(1) phi_s(1);
    Rigidez: int phi_r'' phi_s'' + kR phi_r'(0) phi_s'(0) + kL phi_r(0) phi_s(0)."""
    n = len(ms)
    M, K = np.zeros((n, n)), np.zeros((n, n))
    for r, (_, p) in enumerate(ms):
        for s, (_, q) in enumerate(ms):
            M[r, s] = (quad(lambda z: p[0](z) * q[0](z), 0, 1)[0]
                       + par[mu] * p[0](1) * q[0](1))
            K[r, s] = (quad(lambda z: p[2](z) * q[2](z), 0, 1)[0]
                       + par.get(kapR, 0) * p[1](0) * q[1](0)
                       + par.get(kapL, 0) * p[0](0) * q[0](0))
    return M, K


if __name__ == "__main__":
    sp.pprint(C_r)
    casos = {
        "engastado": (car_eng, F_eng, {mu: 0.3}),
        "molas": (car_mol.subs({aR: kapR / lam, aL: kapL / lam**3}),
                  F_mol.subs({aR: kapR / lam, aL: kapL / lam**3}),
                  {mu: 0.3, kapR: 20.0, kapL: 500.0}),
    }
    for nome, (car, Fr, par) in casos.items():
        ms = modos(car, Fr, par)
        M, K = matrizes(ms, par)
        lr = np.array([m[0] for m in ms])
        print(nome, "lambda_r =", lr.round(4))
        assert np.allclose(M, np.eye(len(ms)), atol=1e-6)            # M = I
        assert np.allclose(K, np.diag(lr**4), rtol=1e-5, atol=1e-5)  # K = lambda^4
    print("ortonormalidade verificada")
