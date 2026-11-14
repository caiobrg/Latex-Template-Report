# Errata das notas manuscritas (`Laboratório 1.pdf`)

Correções feitas ao transformar as notas na introdução teórica
(`content/01_introduction.tex`). Cada item diz **onde estava**, **como estava**, **como
ficou** e **como foi verificado**. Todas as verificações podem ser reproduzidas com os
scripts em `codigo/`, em especial `codigo/valida_notas.py`, que monta as versões das
notas e as compara com as deduções corrigidas.

Notação usada no relatório: $\lambda=\beta L$, $\mu = m/(\rho A L)$,
$\kappa_R = k_R L/EI$, $\kappa_L = k_L L^3/EI$.

---

## A. Caso engastado

### A1. Eq. (5a), coeficiente de $C_2$ (p. 3)
- **Estava:** $C_2\left[\cos\lambda - \cosh\lambda + \tfrac{m\lambda}{\rho A L}(\sin\lambda-\sinh\lambda)\right]$
- **Correto:** $C_2\left[-(\cos\lambda + \cosh\lambda) + \mu\lambda(\sin\lambda-\sinh\lambda)\right]$
- **Origem:** a parcela $-\beta^3C_2(\cos\beta L+\cosh\beta L)$ de (4c) foi transcrita como
  $\cos\lambda-\cosh\lambda$.
- **Verificação:** `caso_engastado.py` (matriz `Amat`); `valida_notas.py` mostra
  `A_notas - A_certa = [[0, 2cos λ], [0, 0]]`.
- O mesmo erro está no elemento $A_{12}$ da matriz $A$ da p. 3.

### A2. Eq. (5b) (p. 3)
- **Estava:** $0 = C_1(\cos\lambda-\cosh\lambda) + C_2(\sin\lambda-\sinh\lambda)$
- **Correto:** $0 = C_1(\cos\lambda+\cosh\lambda) + C_2(\sin\lambda+\sinh\lambda)$ (é a própria
  (4b) com $\beta L=\lambda$).
- **Observação:** a matriz $A$ e a razão $C_2/C_1$ da p. 4 já usavam a forma correta,
  então o erro não se propagou.

### A3. Determinante de $A$ (p. 3)
- Com o $A_{12}$ errado (A1), $\det A$ difere do correto por
  $-2\cos\lambda(\cos\lambda+\cosh\lambda)$ (`valida_notas.py`). A equação final
  escrita nas notas não vem dessa matriz; ela está certa exceto pelo item A4.

### A4. Equação característica (p. 3)
- **Estava:** $\frac{m}{\rho A L}\lambda(\cos\lambda\sinh\lambda - \sin\lambda\cos\lambda) + \cos\lambda\cosh\lambda + 1 = 0$

- **Correto:** $1 + \cos\lambda\cosh\lambda + \mu\lambda(\cos\lambda\sinh\lambda -
  \sin\lambda\cosh\lambda) = 0$ (erro de digitação: $\sin\lambda\cos\lambda \to
  \sin\lambda\cosh\lambda$).
- **Verificação:** `caso_engastado.py` (assert) e, com $\mu=0$, as raízes
  $1{,}8751;\ 4{,}6941;\ 7{,}8548$ da viga engastada-livre clássica.

### A5. Referência de equação (p. 2, última linha)
- **Estava:** "Substituindo (4b) e (4c) em (3c)". **Correto:** "... em **(3d)**".

## B. Normalização

### B1. Constante de normalização, Eq. (8a) e $C_r$ (p. 4)
- **Estava:** $\left(\rho A\int_0^L F_r^2\,dx + m F_r^2(x)\right)^{1/2}$
- **Correto:** $\left(\rho A\int_0^L F_r^2\,dx + m F_r^2(L)\right)^{1/2}$ — a massa
  concentrada entra com o valor do modo **na ponta**, não como função de $x$.
- **Verificação:** `normalizacao.py` confirma numericamente
  $\int_0^1\phi_r\phi_s\,d\xi+\mu\phi_r(1)\phi_s(1)=\delta_{rs}$ para os dois casos.

## C. Resposta forçada e FRF

### C1. Função $G_r$, Eq. (9) (p. 5)
- **Estava:** $G_r(i\Omega) = \dfrac{1}{1+(\Omega/\omega_r)^2 + 2i\zeta_r\Omega/\omega_r}$
- **Correto:** $G_r(i\Omega) = \dfrac{1}{1-(\Omega/\omega_r)^2 + 2i\zeta_r\Omega/\omega_r}$
- **Consequência:** com o sinal $+$ não há ressonância (o denominador nunca se anula
  para $\zeta_r\to0$).
- **Verificação:** `frf.py` resolve a EDO modal com $\eta_r = Y e^{i\Omega t}$ e
  compara com $G_r$ (assert).

### C2. Denominadores em $\eta_{r,p}$, (11) e (12) (p. 5 e 6)
- **Estava:** $\omega_r^2 + \Omega^2 + 2i\zeta_r\Omega\omega_r$
- **Correto:** $\omega_r^2 - \Omega^2 + 2i\zeta_r\Omega\omega_r$ (mesma origem de C1).

### C3. Eq. (12): ponto de saída (p. 6)
- **Estava:** $H(\Omega)=\sum \dfrac{W_r(x)W_r(x_i)}{\cdots}$
- **Correto:** $H(\Omega)=\sum \dfrac{W_r(x_e)W_r(x_i)}{\cdots}$

### C4. Limite inferior do somatório (p. 5)
- **Estava:** $w(x,t)=\sum_{r=0}^{\infty}W_r\eta_r$. **Correto:** $\sum_{r=1}^{\infty}$
  (como nas linhas seguintes da mesma página).

### C5. Notação
- As notas alternam $i$ e $j$ como unidade imaginária; no relatório usa-se só $i$.

## D. Caso com molas na base

### D1. Texto (p. 6)
- **Estava:** "a base da **torre** é modelada ... com uma mola torcional e outra
  **axial**". **Correto:** base da **viga**, com mola torcional e mola
  **transversal** (linear, na direção de $w$). $k_L$ não é rigidez axial; a própria p. 1
  diz "rigidez torcional e transversal".

### D2. Fator $\beta$ omitido em (14a) e na definição de $\alpha_R$ (p. 6 e 7)
- **Estava:** $-C_1 + C_3 = \frac{k_R}{EI}(C_2 + C_4)$, com $\alpha_R = k_R/EI$
- **Correto:** $EI\beta^2(-C_1+C_3) = k_R\beta(C_2+C_4)$, logo
  $-C_1 + C_3 = \alpha_R(C_2 + C_4)$ com $\alpha_R = \dfrac{k_R}{EI\beta} =
  \dfrac{\kappa_R}{\lambda}$.
- **Origem:** $W''$ traz $\beta^2$ e $W'$ traz $\beta$; a razão deixa $1/\beta$.
- **Consequência:** $k_R/EI$ tem unidade de $\text{m}^{-1}$, então somá-la a termos adimensionais não faz sentido. Com a correção, $\alpha_R$ é adimensional e **depende de $\lambda$**.

### D3. Sinal e fator $\beta^3$ em (14b) (p. 7)
- **Estava:** $EIW'''(0)\,g = k_L W(0)\,g \Rightarrow -C_2 + C_4 = \frac{k_L}{EI}(C_1+C_3)$
- **Correto:** de (13b), $-EIW'''(0) = k_L W(0)$:
  $-EI\beta^3(-C_2+C_4) = k_L(C_1+C_3) \Rightarrow C_4 = C_2 - \alpha_L(C_1+C_3)$,
  com $\alpha_L = \dfrac{k_L}{EI\beta^3} = \dfrac{\kappa_L}{\lambda^3}$.
- **Origem:** (13b) estava certa ($Q(0) = -EIw'''(0) = k_Lw(0)$, confirmada por
  princípio variacional), mas o sinal de menos se perdeu ao aplicá-la; além disso faltou
  o $\beta^3$ (como em D2).
- **Consequência:** trocar o sinal equivale a usar uma mola de **rigidez negativa**
  $-k_L$.

### D4. Expressões para $C_3$ e $C_4$, Eq. (15a) e seguinte (p. 7)
- **Estava (herdado de D3):**
  $C_3 = \dfrac{(1+\alpha_L\alpha_R)C_1 + 2\alpha_RC_2}{1-\alpha_L\alpha_R}$,
  $C_4 = \dfrac{2\alpha_LC_1 + (1+\alpha_L\alpha_R)C_2}{1-\alpha_L\alpha_R}$
- **Correto:**
  $C_3 = \dfrac{(1-\alpha_L\alpha_R)C_1 + 2\alpha_RC_2}{1+\alpha_L\alpha_R}$,
  $C_4 = \dfrac{-2\alpha_LC_1 + (1-\alpha_L\alpha_R)C_2}{1+\alpha_L\alpha_R}$
- **Observação:** o denominador correto, $1+\alpha_L\alpha_R$, nunca se anula; o
  "caso singular" $\alpha_L\alpha_R = 1$ das notas (e do `../simbólico.py`) é um
  artefato do erro de sinal. A álgebra de (14) → (15) nas notas está certa para as
  equações (14) que foram escritas (`valida_notas.py` reproduz (15)).

### D5. Rótulos das condições na ponta (p. 8)
- **Estava:** "Para a condição (3a)" / "(3b)". **Correto:** **(3c)** e **(3d)** (as
  condições de ponta; (3a) e (3b) são as do engaste).

### D6. Lado direito de (3d) aplicado (p. 8)
- **Estava:** $\cdots = -\dfrac{m}{\rho A}\beta^4\{\cdots\}$, com o lado esquerdo já
  dividido por $\beta^3$
- **Correto:** $\cdots = -\dfrac{m\beta}{\rho A}\{\cdots\} = -\mu\lambda\{\cdots\}$ (o
  mesmo que foi feito corretamente em (5a) no caso engastado).
- **Consequência:** a equação característica das notas tem $m\beta^4/\rho A$ (dimensão
  de m$^{-3}$) no lugar de $\mu\lambda$ (adimensional). O `../simbólico.py` também usa
  `mu = m*beta**4/(rho*A)`.

### D7. Typo em (16b), coeficiente de $C_2$ (p. 8)
- **Estava:** $-\alpha_L\alpha_R\sinh(\lambda) + \alpha_L\alpha_R\sinh(\lambda)$ (os
  termos se anulariam). **Correto (dado D3):** $-\alpha_L\alpha_R\sin\lambda +
  \alpha_L\alpha_R\sinh\lambda$.
- Fora isso, (16a) e (16b) estão de acordo com as equações (15) das notas (verificado
  com sympy).

### D8. Equação característica (p. 8)
- **Estava:**
  $\alpha_L\alpha_R - 1 + \frac{m\beta^4}{\rho A}\big[(-1-\alpha_L\alpha_R)\sin\lambda\cosh\lambda + (1+\alpha_L\alpha_R)\cos\lambda\sinh\lambda - 2\alpha_L\sin\lambda\sinh\lambda + 2\alpha_R\cos\lambda\cosh\lambda\big] + (\alpha_R-\alpha_L)\sin\lambda\cosh\lambda - (\alpha_L+\alpha_R)\cos\lambda\sinh\lambda + (\alpha_L\alpha_R+1)\cos\lambda\cosh\lambda = 0$
- **Correto:**
  $1+\alpha_L\alpha_R + (\alpha_L\alpha_R-1)\cos\lambda\cosh\lambda + (\alpha_L-\alpha_R)\cos\lambda\sinh\lambda - (\alpha_L+\alpha_R)\sin\lambda\cosh\lambda + \mu\lambda\big[(1-\alpha_L\alpha_R)(\sin\lambda\cosh\lambda - \cos\lambda\sinh\lambda) - 2\alpha_L\sin\lambda\sinh\lambda - 2\alpha_R\cos\lambda\cosh\lambda\big] = 0$
- **Origem:** D3 e D6, mais um erro de sinal na expansão do determinante: mesmo
  partindo das (15) das notas, o termo em $\cos\lambda\sinh\lambda$ sai
  $+(\alpha_L+\alpha_R)$, e não $-(\alpha_L+\alpha_R)$.
- **Verificação:**
  1. `caso_molas.py`: dedução direta + assert; o limite $\alpha_L,\alpha_R\to\infty$
     recupera exatamente a equação do engaste (A4); com $\alpha=0$ e $\mu=0$ recupera
     a viga livre-livre, $\cos\lambda\cosh\lambda = 1$.
  2. `normalizacao.py`: os modos obtidos satisfazem
     $\int_0^1\phi_r''\phi_s''\,d\xi + \kappa_R\phi_r'(0)\phi_s'(0) + \kappa_L\phi_r(0)\phi_s(0)
     = \lambda_r^4\delta_{rs}$ (energia de deformação com molas positivas). Com o sinal
     das notas essa identidade falha (ex.: $\kappa_R=20$, $\kappa_L=500$, $\mu=0{,}3$:
     diagonal $4{,}87;\ 360{,}7;\ 3872$ contra $\lambda_r^4 = 4{,}79;\ 296{,}8;\ 3044$).
  3. `valida_notas.py`: com $\kappa_R = 10^6$, $\mu=0$, a primeira raiz correta cresce
     monotonicamente com $\kappa_L$ até $1{,}8751$ (engaste); a equação das notas dá
     valores sem sentido físico para molas flexíveis (p.ex. $\kappa_L=1$: primeira raiz
     $3{,}92$, acima até da viga engastada).

---

## Itens verificados e mantidos

- Eq. de Euler–Bernoulli, separação de variáveis, (1a), (1b), (2) e (6).
- Condições (3a)–(3d), incluindo o sinal de (3d): $W'''(L) = -\frac{m}{\rho A}\beta^4W(L)$.
- (4a), (4b), (4c), a razão $C_2/C_1$ e a forma modal (7)/(8b) do caso engastado.
- Normalização pela massa (com a correção B1) e a projeção $N_r = F_0W_r(x_i)$.
- (13a) e (13b) do caso com molas (o erro está só na aplicação, D3).
