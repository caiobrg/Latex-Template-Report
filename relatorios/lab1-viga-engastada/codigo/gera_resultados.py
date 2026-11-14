"""Processa os ensaios, ajusta o modelo e gera as figuras (../figures), as tabelas
(../tabelas) e os valores citados no texto (../tabelas/valores.tex).
Uso: python gera_resultados.py   (dentro da pasta codigo/)"""
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scienceplots  # noqa: F401  (registra os estilos)
from matplotlib.ticker import FuncFormatter, MultipleLocator, NullLocator
from scipy.optimize import brentq

import parametros as P
from ajuste_molas import (ajusta_molas, mapa_frequencias, rigidezes_fisicas,
                          sensibilidade, valida)
from dados import frf_impacto, ler_martelo_dewesoft, ler_martelo_serie, ler_shaker
from identificacao import meia_potencia
from modelo import acelerancia, forma_modal, frequencias, modos_normalizados, raizes

RAIZ = Path(__file__).resolve().parent.parent
FIG, TAB = RAIZ / "figures", RAIZ / "tabelas"
FIG.mkdir(exist_ok=True)
TAB.mkdir(exist_ok=True)

# ================================================================ estilo
plt.style.use(["science", "grid"])
plt.rcParams.update({
    "text.latex.preamble": r"\usepackage[T1]{fontenc}\usepackage{libertine}"
                           r"\usepackage[libertine]{newtxmath}\usepackage{amsmath}",
    "font.size": 9, "axes.labelsize": 9, "legend.fontsize": 8,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "grid.alpha": 0.35, "legend.frameon": True, "legend.framealpha": 0.9,
    "legend.edgecolor": "0.8", "savefig.dpi": 300,
})
LARGURA = 5.9  # largura do texto no relatório [pol]
COR = {"medido": "k", "h2": "#9e9e9e", "dewesoft": "#FF9500", "shaker": "#00B945",
       "engaste": "#FF2C00", "molas": "#0C5DA5",
       "engaste_massa": "#FF9500", "molas_massa": "#845B97"}
DESLOC_FASE = 0.45  # fase em [-pi + 0,45; pi + 0,45): o patamar em pi não salta


def br(v, casas=None):
    """Número com vírgula decimal (texto LaTeX)."""
    s = f"{v:.{casas}f}" if casas is not None else f"{v:g}"
    return s.replace(".", "{,}")


def fmt_virgula(ax, eixo="y"):
    f = FuncFormatter(lambda v, _: f"${br(round(v, 6))}$")
    (ax.yaxis if eixo == "y" else ax.xaxis).set_major_formatter(f)


def fmt_freq_log(ax):
    ax.set_xscale("log")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"${v:g}$"))


def fase(H):
    return np.mod(np.angle(H) + np.pi - DESLOC_FASE, 2 * np.pi) - np.pi + DESLOC_FASE


def rotulo_pi(valor, _pos):
    n = round(2 * valor / np.pi)
    if n == 0:
        return "$0$"
    num, den = (n // 2, "") if n % 2 == 0 else (n, "/2")
    coef = {1: "", -1: "-"}.get(num, str(num))
    return rf"${coef}\pi{den}$"


def eixo_fase(ax):
    ax.set_ylim(-np.pi + DESLOC_FASE - 0.15, np.pi + DESLOC_FASE + 0.15)
    ax.yaxis.set_major_locator(MultipleLocator(np.pi / 2))
    ax.yaxis.set_major_formatter(FuncFormatter(rotulo_pi))
    ax.set_ylabel("Fase [rad]")


def decima(freq, mag, n=3000):
    """Índices que preservam máximos e mínimos de |H| em n faixas logarítmicas
    (só para desenhar; a análise usa todos os pontos)."""
    if len(freq) < n:
        return np.arange(len(freq))
    bordas = np.geomspace(max(freq[1], 1e-3), freq[-1], n)
    faixa = np.digitize(freq, bordas)
    idx = []
    for b in np.unique(faixa):
        sel = np.nonzero(faixa == b)[0]
        idx += [sel[np.argmax(mag[sel])], sel[np.argmin(mag[sel])]]
    return np.unique(idx)


def numera_picos(ax, picos, assoc):
    for j, p in enumerate(picos):
        modo = [r for r, jj in assoc.items() if jj == j]
        rot = f"{modo[0] + 1}" if modo else "?"
        ax.annotate(rot, (p["f_r"], p["mag"]), xytext=(0, 3),
                    textcoords="offset points", ha="center", fontsize=7)


def salva(fig, nome):
    fig.savefig(FIG / f"{nome}.pdf")
    plt.close(fig)


# ================================================================ dados
t, F, a = ler_martelo_serie()
fm, H1m, H2m, cohm, t_impactos = frf_impacto(t, F, a)
fd, H1d, H2d, cohd = ler_martelo_dewesoft()
fs, Ts, cohs = ler_shaker()

picos_m = meia_potencia(fm, H1m, cohm)
picos_d = meia_potencia(fd, H1d, cohd)
picos_s = meia_potencia(fs, Ts, cohs, faixa=(10, 2000))

# ======================================= ajuste das molas (picos 1 e 2)
f_med = np.array([p["f_r"] for p in picos_m])
kR, kL, sol = ajusta_molas(f_med[:2])
kR_s, kL_s, _ = ajusta_molas([p["f_r"] for p in picos_s[:2]])
k_R, k_L = rigidezes_fisicas(kR, kL)
N_MOD = 8
f_eng = frequencias(n=N_MOD)
f_mol = frequencias(kR, kL, n=N_MOD)
f_eng_m = frequencias(mu=P.mu_ponta, n=N_MOD)
f_mol_m = frequencias(kR, kL, mu=P.mu_ponta, n=N_MOD)


def associa(picos, f_modelo, tol=0.2):
    """Cada modo do modelo recebe o pico medido mais próximo (até tol)."""
    fp = np.array([p["f_r"] for p in picos])
    assoc = {}
    for r, fr in enumerate(f_modelo):
        j = int(np.argmin(np.abs(fp - fr)))
        if abs(fp[j] - fr) / fr < tol and j not in assoc.values():
            assoc[r] = j
    return assoc


assoc_m = associa(picos_m, f_mol)
assoc_s = associa(picos_s, f_mol)
assoc_d = associa(picos_d, f_mol)
n_val = len(assoc_m)  # modos de flexão identificados no martelo
f_flex = np.array([picos_m[assoc_m[r]]["f_r"] for r in range(n_val)])
_, erro_mol = valida(f_flex, kR, kL)
erro_eng = f_eng[:n_val] / f_flex - 1
casos_sens = sensibilidade(f_flex)


def zetas_modelo(picos, assoc, reserva):
    """zeta_r medido para os modos identificados; 'reserva' para os demais."""
    return np.array([picos[assoc[r]]["zeta"] if r in assoc else reserva[r]
                     for r in range(N_MOD)])


# Verificações citadas na discussão:
# (a) só mola torcional (kappa_L -> inf) ajustada a f_1
kR_t = brentq(lambda k: frequencias(k, 1e12, n=1)[0] - f_flex[0], 0.5, 1e3)
f_tor = frequencias(kR_t, 1e12, n=n_val)
# (b) massa na ponta que levaria f_1 da base elástica ao valor do shaker
mu_sh = brentq(lambda m: frequencias(kR, kL, m, n=1)[0] - picos_s[assoc_s[0]]["f_r"], 0, 0.5)
f_mu_sh = frequencias(kR, kL, mu_sh, n=4)

# (c) engaste com sqrt(EI/rho A)/L^2 corrigido por um único fator c (mínimos
# quadrados em f_r/f_med - 1): mede quanto da diferença é uma escala uniforme
razao = f_eng[:n_val] / f_flex
c_escala = razao.sum() / (razao**2).sum()
erro_escala = c_escala * razao - 1
# (d) sinal da constante modal W_r(x_e) W_r(x_i): fase +pi/2 ou -pi/2 no pico
W_e = modos_normalizados([P.x_acel, P.x_impacto], raizes(n=n_val))
W_m = modos_normalizados([P.x_acel, P.x_impacto], raizes(kR, kL, n=n_val), kR, kL)
const_eng, const_mol = W_e[:, 0] * W_e[:, 1], W_m[:, 0] * W_m[:, 1]
# (e) parcela da energia de deformação de cada modo guardada nas molas da base:
# [kR phi'(0)^2 + kL phi(0)^2] / lambda^4, com int_0^1 phi^2 dxi = 1 (mu = 0)
xi = np.linspace(0, 1, 20001)
energia_molas = []
for lam in raizes(kR, kL, n=n_val):
    phi = forma_modal(xi, lam, kR, kL)
    phi /= np.sqrt(np.trapezoid(phi**2, xi))
    dphi0 = np.gradient(phi, xi, edge_order=2)[0]
    energia_molas.append([kR * dphi0**2 / lam**4, kL * phi[0] ** 2 / lam**4])
energia_molas = np.array(energia_molas)  # colunas: mola torcional, mola transversal

zeta_res = np.full(N_MOD, np.mean([picos_m[assoc_m[r]]["zeta"] for r in range(2, n_val)]))
zeta_m = zetas_modelo(picos_m, assoc_m, zeta_res)
zeta_s = zetas_modelo(picos_s, assoc_s, zeta_m)

# ================================================================ FIGURAS
# --- 1. Série temporal do martelo
fig, ax = plt.subplots(2, 2, figsize=(LARGURA, 3.6), layout="constrained",
                       gridspec_kw={"width_ratios": [2.2, 1]})
ax[0, 0].plot(t, F, color=COR["medido"], lw=0.5)
ax[1, 0].plot(t, a, color=COR["molas"], lw=0.4)
for k, ti in enumerate(t_impactos):
    sel = (t >= ti - 0.002) & (t <= ti + 0.008)
    ax[0, 1].plot((t[sel] - ti) * 1e3, F[sel], lw=0.8, label=f"Impacto {k + 1}")
    sel = (t >= ti) & (t <= ti + 1.5)
    ax[1, 1].semilogy(t[sel] - ti, np.abs(a[sel]) + 1e-3, lw=0.3)
ax[0, 0].set_ylabel("Força [N]")
ax[1, 0].set_ylabel(r"Aceleração [m\,s$^{-2}$]")
ax[1, 0].set_xlabel("Tempo [s]")
ax[0, 1].set_xlabel(r"$t - t_k$ [ms]")
ax[1, 1].set_xlabel(r"$t - t_k$ [s]")
ax[0, 1].legend(fontsize=7)
ax[1, 1].set_ylabel(r"$|a|$ [m\,s$^{-2}$]")
for axx in (ax[0, 0], ax[1, 0], ax[1, 1]):
    fmt_virgula(axx, "x")
salva(fig, "martelo_serie")


# --- 2 e 3. FRFs medidas
def fig_frf(freq, H, coh, cor, rotulo, extra=(), xlim=(5, 2500)):
    fig, (a1, a2, a3) = plt.subplots(3, 1, figsize=(LARGURA, 5.2), sharex=True,
                                     height_ratios=[3, 2, 1.2], layout="constrained")
    for f2, H2, c2, kw in extra:          # curvas de referência, por baixo
        j = decima(f2, np.abs(H2))
        a1.plot(f2[j], np.abs(H2[j]), zorder=1, **kw)
        if c2 is not None:
            a3.plot(f2[j], c2[j], color=kw["color"], lw=kw.get("lw", 0.7), zorder=1)
    i = decima(freq, np.abs(H))
    a1.plot(freq[i], np.abs(H[i]), color=cor, lw=0.7, label=rotulo, zorder=2)
    a2.plot(freq[i], fase(H[i]), color=cor, lw=0.6)
    a3.plot(freq[i], coh[i], color=cor, lw=0.6, zorder=2)
    a1.set_yscale("log")
    a1.set_ylabel(r"$|\mathcal{A}|$ [m\,s$^{-2}$\,N$^{-1}$]")
    eixo_fase(a2)
    a3.set_ylim(0, 1.05)
    a3.set_ylabel(r"$\gamma^2$")
    fmt_virgula(a3)
    a3.set_xlabel("Frequência [Hz]")
    fmt_freq_log(a3)
    a3.set_xlim(*xlim)
    fig.align_ylabels()
    return fig, (a1, a2, a3)


fig, (a1, a2, a3) = fig_frf(
    fm, H1m, cohm, COR["medido"], r"$H_1$ (série temporal, $T = 3$\,s)",
    extra=[(fm, H2m, None, dict(color=COR["h2"], lw=0.6, ls="--",
                                 label=r"$H_2$ (série temporal)")),
           (fd, H1d, cohd, dict(color=COR["dewesoft"], lw=1.0, alpha=0.8,
                                label=r"$H_1$ Dewesoft ($\Delta f = 2{,}44$\,Hz)"))])
a1.set_ylim(1e-2, 2e3)
h, lab = a1.get_legend_handles_labels()
a1.legend([h[2], h[0], h[1]], [lab[2], lab[0], lab[1]], loc="upper left")
numera_picos(a1, picos_m, assoc_m)
salva(fig, "frf_martelo")

fig, (a1, a2, a3) = fig_frf(fs, Ts, cohs, COR["shaker"],
                            r"$\ddot w(L)/f(L)$ (shaker)", xlim=(5, 2000))
a1.plot(fs, np.abs(Ts), "o", color=COR["shaker"], ms=1.2)
numera_picos(a1, picos_s, assoc_s)
a1.set_ylim(5e-2, 3e2)
a1.legend(loc="upper left")
salva(fig, "frf_shaker")

# --- 4. Efeito da massa na ponta: f_r em função de mu
mus = np.linspace(0, 1.0, 81)
fig, ax = plt.subplots(figsize=(LARGURA, 3.4), layout="constrained")
curv_e = np.array([frequencias(mu=m, n=4) for m in mus])
curv_m = np.array([frequencias(kR, kL, mu=m, n=4) for m in mus])
for r in range(4):
    ax.plot(mus, curv_e[:, r], color=COR["engaste"], ls="--",
            label="Base engastada" if r == 0 else None)
    ax.plot(mus, curv_m[:, r], color=COR["molas"],
            label="Base elástica (ajustada)" if r == 0 else None)
    ax.annotate(f"modo {r + 1}", (1.0, curv_e[-1, r]), xytext=(3, 0),
                textcoords="offset points", va="center", fontsize=7)
ax.plot(np.zeros(4), f_flex[:4], "ko", ms=4, clip_on=False, zorder=5,
        label="Medido (martelo, sem massa)")
ax.plot([P.mu_ponta] * 4, f_mol_m[:4], "o", color=COR["molas"], mfc="w", ms=4,
        zorder=5, label=r"Previsto com $m = 0{,}445$\,kg")
ax.axvline(P.mu_ponta, color="0.4", lw=0.8, ls=":")
ax.text(P.mu_ponta + 0.01, 1350, rf"$m = {br(P.m_ponta)}$\,kg ($\mu = {br(P.mu_ponta, 3)}$)",
        fontsize=7, va="center")
ax.set_yscale("log")
ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"${v:g}$"))
ax.set_xlim(0, 1.0)
ax.set_ylim(10, 1600)
ax.set_xlabel(r"Razão de massas $\mu = m/(\rho A L)$")
ax.set_ylabel("Frequência natural [Hz]")
fmt_virgula(ax, "x")
ax.legend(loc="upper left", bbox_to_anchor=(0.03, 0.53))
salva(fig, "massa_ponta")

# --- 5. Ajuste das molas: curvas f1 = f1,med e f2 = f2,med no plano (kR, kL)
lgR, lgL = np.linspace(-0.5, 3.5, 81), np.linspace(1, 6, 91)
F1, F2 = mapa_frequencias(lgR, lgL)
fig, ax = plt.subplots(figsize=(LARGURA, 3.6), layout="constrained")
cf = ax.contourf(lgR, lgL, F1 / f_eng[0], levels=np.linspace(0, 1, 11),
                 cmap="Blues", alpha=0.45)
cb = fig.colorbar(cf, ax=ax, label=r"$f_1 / f_{1,\mathrm{eng}}$")
cb.ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"${br(round(v, 2))}$"))
ax.contour(lgR, lgL, F1, levels=[f_med[0]], colors=[COR["molas"]], linewidths=1.4)
ax.contour(lgR, lgL, F2, levels=[f_med[1]], colors=[COR["engaste"]],
           linewidths=1.4, linestyles="--")
ax.contour(lgR, lgL, F2, levels=[0.99 * f_med[1], 1.01 * f_med[1]],
           colors=[COR["engaste"]], linewidths=0.6, linestyles=":")
ax.plot(np.log10(kR), np.log10(kL), "k*", ms=9, label="Ajuste (martelo)")
ax.plot(np.log10(kR_s), np.log10(kL_s), "kx", ms=6, label="Ajuste (shaker)")
ax.plot([], [], color=COR["molas"], lw=1.4, label=rf"$f_1 = {br(f_med[0], 2)}$\,Hz")
ax.plot([], [], color=COR["engaste"], lw=1.4, ls="--",
        label=rf"$f_2 = {br(f_med[1], 2)}$\,Hz")
ax.plot([], [], color=COR["engaste"], lw=0.6, ls=":", label=r"$f_2 \pm 1\,\%$")
ax.set_xlabel(r"$\log_{10}\kappa_R$")
ax.set_ylabel(r"$\log_{10}\kappa_L$")
fmt_virgula(ax, "x")
ax.legend(loc="lower right")
ax.grid(False)
salva(fig, "ajuste_molas")

# --- 6. Validação: erro relativo por modo
fig, ax = plt.subplots(figsize=(LARGURA, 2.8), layout="constrained")
modos = np.arange(1, n_val + 1)
faixa = np.array([c[4] for c in casos_sens])         # reajustes com +-1 %
ax.fill_between(modos, 100 * faixa.min(axis=0), 100 * faixa.max(axis=0),
                color=COR["molas"], alpha=0.15, lw=0,
                label=r"Base elástica, $f_1$ ou $f_2$ $\pm 1\,\%$")
ax.axvspan(0.5, 2.5, color="0.92", lw=0, zorder=0)
ax.text(1.5, 24.5, "ajuste", ha="center", fontsize=7)
ax.text(4.5, 24.5, "validação", ha="center", fontsize=7)
ax.plot(modos, 100 * erro_eng, "s--", color=COR["engaste"], ms=4, label="Base engastada")
ax.plot(modos, 100 * erro_mol, "o-", color=COR["molas"], ms=4,
        label="Base elástica (ajustada)")
ax.axhline(0, color="k", lw=0.6)
ax.set_xticks(modos)
ax.xaxis.set_minor_locator(NullLocator())
ax.set_xlim(0.5, n_val + 0.5)
ax.set_ylim(-16, 28)
ax.set_xlabel("Modo $r$")
ax.set_ylabel(r"$(f_{r,\mathrm{mod}} - f_{r,\mathrm{med}})/f_{r,\mathrm{med}}$ [\%]")
ax.legend(loc="center right", bbox_to_anchor=(1.0, 0.62))
salva(fig, "validacao")


# --- 7 e 8. FRF analítica x medida
def fig_comparacao(freq_med, H_med, x_e, x_i, zetas, nome, rot_med, cor_med, xlim):
    fmod = np.geomspace(xlim[0], xlim[1], 6000)
    extra = [fr * (1 + np.linspace(-6, 6, 121) * z) for fv in (f_eng, f_mol)
             for fr, z in zip(fv, zetas)]
    fmod = np.unique(np.r_[fmod, np.concatenate(extra)])
    fmod = fmod[(fmod >= xlim[0]) & (fmod <= xlim[1])]
    A_eng = acelerancia(fmod, x_e, x_i, zetas, n=N_MOD)
    A_mol = acelerancia(fmod, x_e, x_i, zetas, kR, kL, n=N_MOD)
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(LARGURA, 4.4), sharex=True,
                                 height_ratios=[3, 2], layout="constrained")
    i = decima(freq_med, np.abs(H_med))
    for axx, fun in ((a1, np.abs), (a2, fase)):
        axx.plot(freq_med[i], fun(H_med[i]), color=cor_med, lw=0.8, label=rot_med)
        axx.plot(fmod, fun(A_eng), color=COR["engaste"], lw=0.9, ls="--",
                 label="Modelo: base engastada")
        axx.plot(fmod, fun(A_mol), color=COR["molas"], lw=0.9,
                 label="Modelo: base elástica ajustada")
    a1.set_yscale("log")
    a1.set_ylabel(r"$|\mathcal{A}|$ [m\,s$^{-2}$\,N$^{-1}$]")
    a1.legend(loc="upper left")
    eixo_fase(a2)
    a2.set_xlabel("Frequência [Hz]")
    fmt_freq_log(a2)
    a2.set_xlim(*xlim)
    fig.align_ylabels()
    salva(fig, nome)


fig_comparacao(fm, H1m, P.x_acel, P.x_impacto, zeta_m, "frf_martelo_modelos",
               r"Medido: martelo ($H_1$)", COR["medido"], (5, 2500))
fig_comparacao(fs, Ts, P.x_shaker, P.x_shaker, zeta_s, "frf_shaker_modelos",
               "Medido: shaker", COR["shaker"], (5, 2000))

# --- 9. Formas modais normalizadas
x = np.linspace(0, P.L, 400)
W_e = modos_normalizados(x, raizes(n=4))
W_m = modos_normalizados(x, raizes(kR, kL, n=4), kR, kL)
fig, axs = plt.subplots(2, 2, figsize=(LARGURA, 3.9), sharex=True, layout="constrained")
for r, axx in enumerate(axs.flat):
    axx.plot(x, np.sign(W_e[r, -1]) * W_e[r], color=COR["engaste"], ls="--",
             label="Base engastada")
    axx.plot(x, np.sign(W_m[r, -1]) * W_m[r], color=COR["molas"],
             label="Base elástica")      # mesma convenção de sinal: W_r(L) > 0
    axx.axhline(0, color="k", lw=0.5)
    for xp in (P.x_impacto, P.x_acel):
        axx.axvline(xp, color="0.25", lw=0.8, ls="-.")
    axx.set_title(rf"Modo {r + 1}: $f_{{\mathrm{{eng}}}} = {br(f_eng[r], 1)}$\,Hz, "
                  rf"$f_{{\mathrm{{el}}}} = {br(f_mol[r], 1)}$\,Hz", fontsize=8)
    fmt_virgula(axx, "x")
    fmt_virgula(axx, "y")
axs[0, 0].legend(loc="upper left", fontsize=7)
for axx in axs[1]:
    axx.set_xlabel("$x$ [m]")
for axx in axs[:, 0]:
    axx.set_ylabel(r"$W_r(x)$ [kg$^{-1/2}$]")
y0 = axs[0, 0].get_ylim()[0]
axs[0, 0].text(P.x_impacto - 0.005, y0 + 0.1, "impacto", fontsize=6, ha="right", va="bottom")
axs[0, 0].text(P.x_acel - 0.005, y0 + 0.1, "acelerômetro", fontsize=6, ha="right",
               va="bottom")
salva(fig, "modos")

# ================================================================ TABELAS
n = r"\num"


def linha(*c):
    return " & ".join(c) + r" \\"


# Martelo: Dewesoft x série temporal
L_ = []
for r in range(n_val):
    pm = picos_m[assoc_m[r]]
    pd_ = picos_d[assoc_d[r]] if r in assoc_d else None
    cel_d = ([f"{n}{{{pd_['f_r']:.2f}}}", f"{n}{{{100 * pd_['zeta']:.2f}}}",
              f"{pd_['pontos']}"] if pd_ else ["--"] * 3)
    L_.append(linha(f"{r + 1}", *cel_d, f"{n}{{{pm['f_r']:.2f}}}",
                    f"{n}{{{100 * pm['zeta']:.2f}}}", f"{pm['pontos']}",
                    f"{n}{{{pm['coerencia']:.3f}}}"))
(TAB / "martelo_resolucao.tex").write_text("\n".join(L_) + "\n")

# Martelo x shaker
L_ = []
for r in range(n_val):
    pm = picos_m[assoc_m[r]]
    cel = [f"{r + 1}", f"{n}{{{pm['f_r']:.2f}}}", f"{n}{{{100 * pm['zeta']:.2f}}}"]
    if r in assoc_s:
        ps = picos_s[assoc_s[r]]
        cel += [f"{n}{{{ps['f_r']:.2f}}}", f"{n}{{{100 * ps['zeta']:.2f}}}",
                f"{ps['pontos']}", f"{n}{{{ps['coerencia']:.3f}}}",
                f"{n}{{{100 * (ps['f_r'] / pm['f_r'] - 1):+.1f}}}"]
    else:
        cel += ["--"] * 5
    L_.append(linha(*cel))
for j, ps in enumerate(picos_s):
    if j not in assoc_s.values():
        L_.append(linha("?", "--", "--", f"{n}{{{ps['f_r']:.2f}}}",
                        f"{n}{{{100 * ps['zeta']:.2f}}}", f"{ps['pontos']}",
                        f"{n}{{{ps['coerencia']:.3f}}}", "--"))
(TAB / "martelo_shaker.tex").write_text("\n".join(L_) + "\n")

# Efeito da massa na ponta
L_ = []
for r in range(6):
    fmed = f"{n}{{{f_flex[r]:.1f}}}" if r < n_val else "--"
    L_.append(linha(f"{r + 1}", fmed, f"{n}{{{f_eng[r]:.1f}}}", f"{n}{{{f_eng_m[r]:.1f}}}",
                    f"{n}{{{100 * (f_eng_m[r] / f_eng[r] - 1):.1f}}}",
                    f"{n}{{{f_mol[r]:.1f}}}", f"{n}{{{f_mol_m[r]:.1f}}}",
                    f"{n}{{{100 * (f_mol_m[r] / f_mol[r] - 1):.1f}}}"))
(TAB / "massa_ponta.tex").write_text("\n".join(L_) + "\n")

# Validação
L_ = []
for r in range(n_val):
    L_.append(linha(f"{r + 1}", f"{n}{{{f_flex[r]:.1f}}}", f"{n}{{{f_eng[r]:.1f}}}",
                    f"{n}{{{100 * erro_eng[r]:+.1f}}}", f"{n}{{{f_mol[r]:.1f}}}",
                    f"{n}{{{100 * erro_mol[r]:+.1f}}}",
                    f"{n}{{{f_flex[r] / f_eng[r]:.3f}}}",
                    "ajuste" if r < 2 else "validação"))
(TAB / "validacao.tex").write_text("\n".join(L_) + "\n")

# Sensibilidade do ajuste
L_ = []
for i, d, a_, b_, e in casos_sens:
    L_.append(linha(rf"$f_{i}$ ${'+' if d > 0 else '-'}\,1\,\%$", f"{n}{{{a_:.2f}}}",
                    f"{n}{{{b_:.0f}}}",
                    *[f"{n}{{{100 * v:+.1f}}}" for v in e[2:]]))
(TAB / "sensibilidade.tex").write_text("\n".join(L_) + "\n")

# Sinal das constantes modais x fase medida no pico
L_ = []
for r in range(n_val):
    ph = picos_m[assoc_m[r]]["fase"]
    L_.append(linha(f"{r + 1}", f"{n}{{{ph:+.2f}}}", "$+$" if ph > 0 else "$-$",
                    f"{n}{{{const_eng[r]:+.3f}}}", f"{n}{{{const_mol[r]:+.3f}}}"))
(TAB / "sinais.tex").write_text("\n".join(L_) + "\n")

# Valores citados no texto
pm = {r: picos_m[assoc_m[r]] for r in range(n_val)}
val = {
    "kapR": f"{kR:.2f}", "kapL": f"{kL:.0f}", "kR": f"{k_R:.0f}", "kL": f"{k_L / 1e6:.2f}",
    "kapRs": f"{kR_s:.2f}", "kapLs": f"{kL_s:.0f}",
    "mviga": f"{P.m_viga:.3f}", "muponta": f"{P.mu_ponta:.3f}",
    "fum": f"{pm[0]['f_r']:.2f}", "fdois": f"{pm[1]['f_r']:.2f}",
    "zum": f"{100 * pm[0]['zeta']:.2f}",
}
macros = [rf"\newcommand{{\v{k}}}{{\num{{{v}}}}}" for k, v in val.items()]
(TAB / "valores.tex").write_text("\n".join(macros) + "\n")

# Resumo para conferência
resumo = {
    "impactos_s": t_impactos.round(4).tolist(),
    "martelo": [{k: round(float(v), 5) for k, v in p.items()} for p in picos_m],
    "dewesoft": [{k: round(float(v), 5) for k, v in p.items()} for p in picos_d],
    "shaker": [{k: round(float(v), 5) for k, v in p.items()} for p in picos_s],
    "kR": kR, "kL": kL, "k_R": k_R, "k_L": k_L, "kR_s": kR_s, "kL_s": kL_s,
    "f_eng": f_eng.round(2).tolist(), "f_mol": f_mol.round(2).tolist(),
    "f_eng_massa": f_eng_m.round(2).tolist(), "f_mol_massa": f_mol_m.round(2).tolist(),
    "f_mol_ajuste_shaker": frequencias(kR_s, kL_s, n=n_val).round(2).tolist(),
    "erro_eng_%": (100 * erro_eng).round(2).tolist(),
    "erro_mol_%": (100 * erro_mol).round(2).tolist(),
    "zeta_m": zeta_m.round(5).tolist(), "zeta_s": zeta_s.round(5).tolist(),
    "so_torcional": {"kR": kR_t, "f": f_tor.round(2).tolist(),
                     "erro_%": (100 * (f_tor / f_flex - 1)).round(2).tolist()},
    "massa_eq_shaker": {"mu": mu_sh, "m_g": 1e3 * mu_sh * P.m_viga,
                        "f": f_mu_sh.round(2).tolist(),
                        "var_%": (100 * (f_mu_sh / f_mol[:4] - 1)).round(2).tolist()},
    "jacobiano_ajuste": sol.jac.round(4).tolist(),
    "escala_engaste": {"c": c_escala, "EI_rhoA_rel": c_escala**2,
                       "L_ef": P.L / np.sqrt(c_escala),
                       "erro_%": (100 * erro_escala).round(2).tolist()},
    "energia_molas_%": {"torcional": (100 * energia_molas[:, 0]).round(1).tolist(),
                        "transversal": (100 * energia_molas[:, 1]).round(1).tolist()},
}
(TAB / "resumo.json").write_text(json.dumps(resumo, indent=1))
print(json.dumps({k: v for k, v in resumo.items()
                  if k not in ("martelo", "dewesoft", "shaker")}, indent=1))
