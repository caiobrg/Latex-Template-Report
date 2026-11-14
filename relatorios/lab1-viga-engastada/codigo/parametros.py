"""Parâmetros da viga e da instrumentação (Tabela 1 e Seção 2.1 do roteiro)."""
import numpy as np

# --- Geometria e material
b = 0.026        # largura [m]
h = 0.007        # espessura [m]
L = 0.43         # comprimento livre [m]
E = 210e9        # módulo de elasticidade [Pa]
rho = 8069.4     # massa específica [kg/m^3]

A = b * h               # área da seção [m^2]
I = b * h**3 / 12       # momento de inércia de área [m^4]
EI = E * I              # rigidez à flexão [N m^2]
rhoA = rho * A          # massa por unidade de comprimento [kg/m]
m_viga = rhoA * L       # massa da viga [kg]

# --- Massa concentrada da Figura 1 do roteiro (m_end). O roteiro a lista como
# m_b = 0,445 kg, mas ela não é a massa da viga (rho A L = 0,632 kg). Os ensaios
# foram feitos SEM essa massa; ela só entra na comparação da Seção 3.3.
m_ponta = 0.445                 # [kg]
mu_ponta = m_ponta / m_viga     # razão de massas mu = m/(rho A L)

# --- Pontos de medição
x_impacto = 0.200       # martelo: ponto de impacto (entrada) [m]
x_acel = 0.356          # martelo: acelerômetro (saída) [m]
x_shaker = L            # shaker: cabeça de impedância na ponta (FRF pontual)


# <<conversao>>
def lambda_para_freq(lam):
    """f_r = omega_r / 2 pi, com omega_r = (lambda_r/L)^2 sqrt(EI / rho A)  [Hz]."""
    return (np.asarray(lam) / L) ** 2 * np.sqrt(EI / rhoA) / (2 * np.pi)


def freq_para_lambda(f):
    """Inversa de lambda_para_freq."""
    return L * np.sqrt(2 * np.pi * np.asarray(f) / np.sqrt(EI / rhoA))
# <<fim_conversao>>
