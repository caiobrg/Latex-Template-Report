"""Identificação de frequências naturais e fatores de amortecimento a partir das FRFs
medidas: picos por proeminência e amortecimento pelo método da meia potência."""
import numpy as np
from scipy.signal import find_peaks


# <<meia_potencia>>
def meia_potencia(freq, H, coerencia, faixa=(10, 2400), proeminencia=0.5):
    """Para cada pico de |H| (proeminência em décadas), retorna f_r, as frequências
    de meia potência f_1 e f_2, em que |H| = |H|_max / sqrt(2), e
    zeta_r = (f_2 - f_1) / (2 f_r)."""
    freq, mag = np.asarray(freq), np.abs(H)
    picos, _ = find_peaks(np.log10(mag), prominence=proeminencia)
    picos = [k for k in picos if faixa[0] <= freq[k] <= faixa[1]]

    resultado = []
    for k in picos:
        alvo = mag[k] / np.sqrt(2)
        f_meia = []
        for passo in (-1, 1):
            j = k                    # anda a partir do pico até cruzar o alvo
            while 0 <= j + passo < len(mag) and mag[j + passo] >= alvo:
                j += passo
            par = [j + passo, j]     # interpola em log f entre os pontos vizinhos
            f_meia.append(np.exp(np.interp(alvo, mag[par], np.log(freq[par]))))
        f1, f2 = f_meia
        resultado.append({
            "f_r": freq[k], "f_1": f1, "f_2": f2,
            "zeta": (f2 - f1) / (2 * freq[k]),
            "pontos": int(np.sum((freq > f1) & (freq < f2))),  # pontos na banda
            "coerencia": float(np.asarray(coerencia)[k]),
            "mag": mag[k], "fase": float(np.angle(H[k])),
        })
    return resultado
# <<fim_meia_potencia>>
