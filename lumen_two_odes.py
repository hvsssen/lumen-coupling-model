import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# ---- Parametres du probleme ----
C12 = 1.0        # conductivité de la jonction entre les deux lumens
gamma1 = 1.0     # tension de surface du lumen 1
gamma2 = 1.0     # tension de surface du lumen 2
gammaC = 1.0     # tension de la jonction cellule-cellule (fixe)
R1_0 = 4.0       # rayon initial du lumen 1
R2_0 = 4.0     # rayon initial du lumen 2
t_max = 50.0     # temps maximal de simulation

# Le lumen est une lentille formee de deux segments circulaires :
# A(R) = R^2 * (2*theta - sin(2*theta))
# theta ne depend que des tensions (fixes) donc theta est CONSTANT dans le temps
theta1 = np.arccos(gammaC / (2 * gamma1))
theta2 = np.arccos(gammaC / (2 * gamma2))
g1 = 2 * theta1 - np.sin(2 * theta1)
g2 = 2 * theta2 - np.sin(2 * theta2)

#A = (R^2)*(2*téta - 2sin(2*téta)) = R^2 * g
def aire(R, g):
    return R**2 * g

# En derivant A_i = R_i^2 * g_i (g_i constant) on a 
# dA_i/dt = 2*g_i* R_i* dR_i/dt
# donc en remplacant dans dA1/dt = -J12 et dA2/dt = +J12 on obtient directement
# les equations differentielles sur R1 et R2 :
#dR1/dt = C12·(γ1/R1 - γ2/R2) / (2R1·(2θ1-sin2θ1))
#dR2/dt = C12·(γ2/R2 - γ1/R1) / (2R2·(2θ2-sin2θ2))

def equations(t, y):
    R1, R2 = y
    P1 = gamma1 / R1
    P2 = gamma2 / R2
    J12 = C12 * (P1 - P2) # = C12((gamma1/R1) - (gamma2/R2))
    dR1_dt = -J12 / (2 * g1 * R1)
    dR2_dt = J12 / (2 * g2 * R2)
    return [dR1_dt, dR2_dt]

# Un lumen ne peut pas avoir un rayon negatif : des que l'un des deux est
# (presque) vide on arrete le calcul, sinon P = gamma/R explose
def lumen_vide(t, y):
    R1, R2 = y
    return min(R1, R2) - 1e-4
lumen_vide.terminal = True

# ---- Resolution numerique ----
solution = solve_ivp(
    equations,
    (0, t_max),
    [R1_0, R2_0],
    events=lumen_vide,
    max_step=0.1,
    dense_output=True,
)

t = np.linspace(0, solution.t[-1], 400)
R1, R2 = solution.sol(t)
P1, P2 = gamma1 / R1, gamma2 / R2
A1, A2 = aire(R1, g1), aire(R2, g2)

if solution.t_events[0].size > 0:
    print(f"Un des deux lumens se vide a t = {solution.t_events[0][0]:.2f}")
print(f"Verification de la conservation A1+A2 : {A1[0]+A2[0]:.3f} -> {A1[-1]+A2[-1]:.3f}")

# ---- Trace des resultats ----
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(6, 9), sharex=True)

ax1.plot(t, R1, label="R1")
ax1.plot(t, R2, label="R2")
ax1.set_ylabel("Rayon")
ax1.legend()

ax2.plot(t, A1, label="A1")
ax2.plot(t, A2, label="A2")
ax2.set_ylabel("Aire")
ax2.legend()

ax3.plot(t[:-10], P1[:-10], label="P1")
ax3.plot(t[:-10], P2[:-10], label="P2")
ax3.set_ylabel("Pression")
ax3.set_xlabel("t")
ax3.legend()

plt.tight_layout()
plt.savefig("lumen_gamma1=gamma2_R1_0=2.0_R2_0=4.0.png", dpi=150)
plt.show()
