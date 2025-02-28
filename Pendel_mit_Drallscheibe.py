
# In[0]

from ELH_mod import EL_eq, EL_Vars, Prog_EH_py, Prog_EL_py
from ELH_mod import Prog_EH_opt_py, Prog_EL_opt_py
import sympy as sp
from IPython.display import display
sp.init_printing()

# In[1]

def Dm(x):
    """.

    Drehmatrix im R3

    Parameters
    ----------
    x : Winkel
    Returns
    -------
    TYPE : Drehmatrix
    """
    return sp.Matrix([[sp.cos(x), -sp.sin(x)],
                      [sp.sin(x), sp.cos(x)]])

# In[2]

Sys = EL_Vars('t', ['\\varphi', '\\theta'])

# m1 = Stab, m2 = Scheibe
m1, m2, g, l, C = sp.symbols('m_{Stab} m_{Scheibe} g l C')

# Ortsvektor - Schwerpunkt Stab
r_SStab = Dm(Sys.x[0]) * sp.Matrix([l/2, 0])
v_SStab = r_SStab.applyfunc(Sys.ttd)

# Ortsvektor - Schwerpunkt scheibe
r_SScheibe = Dm(Sys.x[0]) * sp.Matrix([l, 0])
v_SScheibe = r_SScheibe.applyfunc(Sys.ttd)

# Kinetische Energie
E_kinStab = m1/2 * (v_SStab.T * v_SStab)[0] + 1/2 * 1/12 * m1 * l**2 * Sys.v[0]**2
E_kinScheibe = m2/2 * (v_SScheibe.T * v_SScheibe)[0] + 1/2 * C * (Sys.v[0]+Sys.v[1])**2

E_kin = E_kinStab + E_kinScheibe
E_kin = sp.simplify(E_kin)

# Potentielle Energie
g_vec = sp.Matrix([0, -g]).T
E_potStab = - (m1 * g_vec * r_SStab)[0]
E_potScheibe = - (m2 * g_vec * r_SScheibe)[0]

E_pot = E_potStab + E_potScheibe
E_pot = sp.simplify(E_pot)

# Lagrangian und Bewegungsgleichung
L = E_kin - E_pot
res = EL_eq(L, Sys)
display(res)


# Programm erstellen und einsetzen der Größen
M, RhsE = sp.linear_eq_to_matrix(res, Sys.a)
Prog_EL_py(M, RhsE, Sys, 'x')

# Scheibe mit d=15cm und 1cm dicke:
# m = r^2 * pi * h*rho
m2_ = (0.075**2 * sp.pi * 0.01 * 2700).evalf()
# I = 1/2 * m * r^2
C_ = 1/2 * m2_ * 0.075**2
# display(m2_, C_)


sdic = {m1: .25, g: 9.81, l: .3, m2: m2_, C: C_}

with open('EL_func.py', 'w') as f:
    Prog_EL_py(M.subs(sdic), RhsE.subs(sdic), Sys, 'x', f=f)

# In[3]

from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from numpy import arange


# In[4]

# Plotten der Bewegung
from EL_func import EL_func
# from EH_func import EH_func

Te = 10
DT = 0.1
X0 = [-0, 0, 0, 0]
EL_func(0, X0)

Sol_L = solve_ivp(EL_func, [0, Te], X0, rtol=1e-7,
                  t_eval=arange(0, Te + DT, DT))
display(len(Sol_L))
plt.plot(Sol_L.t, Sol_L.y[0], 'r', label='$\\varphi(t)$')
plt.plot(Sol_L.t, Sol_L.y[1], 'g', label='$\\theta(t)$')
plt.plot(Sol_L.t, Sol_L.y[2], 'b', label='$\dot{\\varphi}(t)$')
plt.plot(Sol_L.t, Sol_L.y[3], 'c', label='$\dot{\\theta}(t)$')
plt.xlabel('t [s]')
plt.legend()
plt.title('Demo Simulation EL')
plt.grid(True)
plt.show()

# In[5]

# Ruhelagenberechnung

# Erstelle Substitutionsvariable für Ruhelage (v = a = 0)
subs_RL = {}
for key in Sys.a:
    subs_RL[key] = 0
for key in Sys.v:
    subs_RL[key] = 0

# Parameter in Bewegungsgleichung substituieren und lösen
res_RL = res.subs(subs_RL)

sol_RL = sp.solve(res_RL, Sys.x)
print("Ruhelagen:")
for sol in sol_RL:
    display(sp.Eq(sp.Matrix(Sys.x), sp.Matrix(sol), evaluate=False))
# display(res_RL, sol_RL)

display(f"Es ist ersichtlich, dass die Ruhelagen von theta unabhängig sind (theta ist invariant)")

# In[6]

# Linearisierung
# display(El_func(0, x0))