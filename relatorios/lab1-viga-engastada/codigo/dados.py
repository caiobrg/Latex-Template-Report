"""Leitura dos arquivos dos ensaios e cálculo da FRF do martelo a partir da série
temporal. Arquivos originais:
  martelo_frf.csv              <- FRFs_martelo.csv (exportação Dewesoft, H1/H2)
  martelo_serie_temporal.csv   <- ensaio_martelo.csv (força e aceleração, 5 kHz)
  shaker_transmissibilidade.txt <- transmissibility1,0(f) Oct 20, 2026 15-47-20.txt
  shaker_coerencia.txt          <- coherence1,0(f) Oct 20, 2026 15-47-20.txt
"""
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks

DADOS = Path(__file__).resolve().parent.parent / "dados"

# Nos dois ensaios o sentido positivo da força é oposto ao do acelerômetro: abaixo
# da 1a ressonância a fase medida é ~0, quando deveria ser pi (a = -Omega^2 w, com
# w em fase com F). As FRFs medidas são multiplicadas por -1.
SINAL = -1


def ler_martelo_dewesoft():
    """FRF exportada pelo Dewesoft (blocos de 2048 pontos, df = 2,44 Hz)."""
    linhas = [ln.split() for ln in open(DADOS / "martelo_frf.csv")]
    d = {ln[0]: np.array(ln[1:], float) for ln in linhas}
    f = d["Frequency(Hz)"]
    H1 = d["TM/FRF_H1_1Z+/1Z-_Ampl(m/s2/N)"] * np.exp(1j * d["TM/FRF_H1_1Z+/1Z-_Phase(rad)"])
    H2 = d["TM/FRF_H2_1Z+/1Z-_Ampl(m/s2/N)"] * np.exp(1j * d["TM/FRF_H2_1Z+/1Z-_Phase(rad)"])
    return f, SINAL * H1, SINAL * H2, d["TM/Coherence_1Z+/1Z-"]


def ler_martelo_serie():
    """Série temporal: tempo [s], força do martelo [N] e aceleração [m/s^2]."""
    texto = open(DADOS / "martelo_serie_temporal.csv", encoding="latin-1").read()
    linhas = texto.splitlines()
    i0 = next(i for i, ln in enumerate(linhas) if ln.startswith("Time(s)"))
    D = np.array([ln.split()[:3] for ln in linhas[i0 + 1:] if ln.strip()], float)
    return D[:, 0], D[:, 1], D[:, 2]


def ler_dactron(nome):
    """Exportação ASCII do Dactron: 32 linhas de cabeçalho + 3 de rótulos."""
    return np.loadtxt(DADOS / nome, skiprows=35)


def ler_shaker():
    """Transmissibilidade a/F na ponta (varredura senoidal, passo logarítmico)."""
    T = ler_dactron("shaker_transmissibilidade.txt")
    coh = ler_dactron("shaker_coerencia.txt")
    return T[:, 0], SINAL * (T[:, 1] + 1j * T[:, 2]), coh[:, 1]


# <<frf_serie>>
def frf_impacto(t, F, a, T_janela=3.0, t_pre=0.01, t_forca=0.02, nfft=2**17):
    """H1, H2 e coerência a partir dos impactos registrados na série temporal.

    Cada impacto gera um registro de T_janela segundos (a resposta decai até o
    nível de ruído antes do fim, então não há vazamento). A força é zerada
    fora dos primeiros t_forca segundos (janela de força), e os registros são
    completados com zeros até nfft pontos, o que apenas interpola o espectro.
    """
    fs = 1 / np.median(np.diff(t))
    impactos, _ = find_peaks(F, height=0.1 * F.max(), distance=int(fs))
    n, n_pre, n_forca = int(T_janela * fs), int(t_pre * fs), int(t_forca * fs)

    Sff = Saa = Sfa = 0
    for k in impactos:
        f_reg = F[k - n_pre:k - n_pre + n].copy()
        a_reg = a[k - n_pre:k - n_pre + n]
        f_reg[n_pre + n_forca:] = 0.0                     # janela de força
        Fw, Aw = np.fft.rfft(f_reg, nfft), np.fft.rfft(a_reg, nfft)
        Sff = Sff + np.abs(Fw) ** 2                       # autoespectros e
        Saa = Saa + np.abs(Aw) ** 2                       # espectro cruzado,
        Sfa = Sfa + np.conj(Fw) * Aw                      # médios nos impactos
    freq = np.fft.rfftfreq(nfft, 1 / fs)
    H1 = Sfa / Sff                                        # ruído na saída
    H2 = Saa / np.conj(Sfa)                               # ruído na entrada
    coerencia = np.abs(Sfa) ** 2 / (Sff * Saa)
    return freq, SINAL * H1, SINAL * H2, coerencia, t[impactos]
# <<fim_frf_serie>>
