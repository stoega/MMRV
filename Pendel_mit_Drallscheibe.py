
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

M_ = sp.symbols('M')

# Ortsvektor - Schwerpunkt Stab
r_SStab = Dm(Sys.x[0]) * sp.Matrix([l/2, 0])
v_SStab = r_SStab.applyfunc(Sys.ttd)

# Ortsvektor - Schwerpunkt scheibe
r_SScheibe = Dm(Sys.x[0]) * sp.Matrix([l, 0])
v_SScheibe = r_SScheibe.applyfunc(Sys.ttd)

# Kinetische Energie
E_kinStab = m1/2 * (v_SStab.T * v_SStab)[0] + 1/2 * 1/12 * m1 * l**2 * Sys.v[0]**2
E_kinScheibe = m2/2 * (v_SScheibe.T * v_SScheibe)[0] + 1/2 * C * (Sys.v[0] + Sys.v[1])**2

E_kin = E_kinStab + E_kinScheibe
E_kin = sp.simplify(E_kin)

# Potentielle Energie
g_vec = sp.Matrix([0, -g]).T
E_potStab = - (m1 * g_vec * r_SStab)[0]
E_potScheibe = - (m2 * g_vec * r_SScheibe)[0]

E_pot = E_potStab + E_potScheibe
E_pot = sp.simplify(E_pot)

# Lagrangian und Bewegungsgleichung
L = E_kin - E_pot + M_*Sys.x[1] # Hier wird Motormoment eingebracht
res = EL_eq(L, Sys).subs(1.0, 1)
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
display(m2_)

# In[]

# M = 0 für Simulation
sdic = {m1: .1, g: 9.81, l: .3, m2: m2_, C: C_, M_: 0}

with open('EL_func.py', 'w') as f:
     Prog_EL_py(M.subs(sdic), RhsE.subs(sdic), Sys, 'x', f=f)

# In[3]

from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from numpy import arange
from importlib import reload


# In[4]

# Plotten der Bewegung
import EL_func as EL
reload(EL)  # Needed to catch changes

Te = 10
DT = 0.1
# X0 = [-3.1415/2, 0, 0, 0]
X0 = [0, 0, 0, 0]

Sol_L = solve_ivp(EL.EL_func, [0, Te], X0, rtol=1e-7,
                  t_eval=arange(0, Te + DT, DT))
# display(Sol_L)
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

# Linearisierung um [pi/2, 0]
subs_RL[Sys.x[0]] = sol_RL[0][0]
subs_RL[Sys.x[1]] = sol_RL[0][1]
display(subs_RL)

# Aufstellen x_dot = f(x, u)
f = (M.inv()*RhsE).subs(1.0, 1)
display(f)

# Linearisierungsparameter
delta_phi, delta_theta = sp.symbols('\\Delta{\\varphi} \\Delta{\\theta}') 
delta_x = sp.Matrix([delta_phi, delta_theta])

# f_lin = f(x0) + f'(x0) * delta_x
f_lin = f.subs(subs_RL) + f.jacobian(Sys.x).subs(subs_RL) * delta_x
display(f_lin)

# Zustandsraumdarstellung
A = f_lin.jacobian(delta_x)
A = sp.simplify(A)
b = f_lin.jacobian(sp.Matrix([M_]))
b = sp.simplify(b)

display(A, b)

# In[7]

#Eigenvalues
evs = (M.inv()*RhsE.jacobian(Sys.x)).eigenvals()
for ev in evs:
    print(sp.latex(ev.simplify()))

# haltemoment
M_halt = sp.solve((M.inv()*RhsE)[0], M_)[0]
print(sp.latex(M_halt))


