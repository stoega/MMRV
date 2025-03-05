
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
m1, m2, g, l, C = sp.symbols('m_1 m_2 g l C')

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
display(m2_ ,C_)

# In[]

# M = 0 für Simulation
params = {m1: .1, g: 9.81, l: .3, m2: m2_, C: C_, M_: 0}

with open('EL_func.py', 'w') as f:
     Prog_EL_py(M.subs(params), RhsE.subs(params), Sys, 'x', f=f)

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
X0 = [-3.1415/2, 0, 0, 0]
# X0 = [0, 0, 0, 0]

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
upper = True
subs_RL[Sys.x[0]] = sol_RL[0][0] if upper else sol_RL[1][0]
subs_RL[Sys.x[1]] = sol_RL[0][1] if upper else sol_RL[1][1]
display(subs_RL)

# Aufstellen x_dot = f(x, u)
x_vec = sp.Matrix([sp.Matrix(Sys.x), sp.Matrix(Sys.v)])
x_ddot = (M.inv()*RhsE).subs(1.0, 1)

f = sp.Matrix([sp.Matrix(Sys.v), x_ddot])
display(sp.Eq(sp.S('f(x, u)'), f, evaluate=False))

# Linearisierungsparameter
delta_phi, delta_theta, delta_phi_dot, delta_theta_dot = sp.symbols('\\Delta{\\varphi} \\Delta{\\theta} \\Delta{\\dot{\\varphi}} \\Delta{\\dot{\\theta}}') 
delta_x = sp.Matrix([delta_phi, delta_theta, delta_phi_dot, delta_theta_dot])
# delta_x = sp.Matrix([delta_phi, delta_theta])

# f_lin = f(x0) + f'(x0) * delta_x
f_lin = f.subs(subs_RL) + f.jacobian(x_vec).subs(subs_RL) * delta_x
display(sp.Eq(sp.symbols('f(\\Delta x, \\Delta u)'), f_lin, evaluate=False))

# Zustandsraumdarstellung
A = f_lin.jacobian(delta_x)
A = sp.simplify(A)
b = f_lin.jacobian(sp.Matrix([M_]))
b = sp.simplify(b)

display(A, b)

# In[7]

#Eigenvalues
evs = A.eigenvals()
display(evs)
for ev in evs:
    # print(sp.latex(ev.simplify()))
    display('EW: ', ev.simplify())
    display('Vielfachheit: ', evs[ev])

# haltemoment
M_halt = sp.solve((M.inv()*RhsE)[0], M_)[0]
# print(sp.latex(M_halt))
display(M_halt)



# %%

# Abtastsystem

from support import Equ, Sol
import numpy as np
import control as ctrl
from scipy.optimize import linprog

n = A.shape[0]

sys = ctrl.ss(sp.matrix2numpy(A.subs(params)), sp.matrix2numpy(b.subs(params)), np.eye(n), np.zeros((n, 1)))

x0 = np.array(n*[[0]])
xN = np.array(n*[[0]], np.float64)
xN[0, 0] = 10*np.pi/180.0    # TODO

# Abtastsystem

T = 0.1     # Abtastzeit
sysd = ctrl.sample_system(sys, T, method='zoh')

# Schrittanzahl

N = 300      # Anzahl der Schritte

# l_1 mit u Begrenzung

M = Equ(np.array(sysd.A), np.array(sysd.B), N)

# Sympole für l_\inf und l1 Lösungen

U = [sp.Symbol('u_' + str(i)) for i in range(N)]
Up = [sp.Symbol('u⁺_' + str(i)) for i in range(N)]
Um = [sp.Symbol('u⁻_' + str(i)) for i in range(N)]
umax = sp.Symbol('u_max')

# Gleichung

Gln_xn = list(M @ U - xN[:, 0])

# Gleichungen für theta_dot
V_theta = [sp.Symbol('vtheta_' + str(i)) for i in range(1, N)]

# Betrachten nur 4.ten Eintrag von x_k=0..N -> M[3,...] * U liefert skalar 
Gln_v_theta = [(M[3:, -i - 2: -1] @ U[: i + 1])[0] - v for i, v in enumerate(V_theta)]

V_phi = [sp.Symbol('vphi_' + str(i)) for i in range(1, N)]

# Betrachten nur 4.ten Eintrag von x_k=0..N -> M[3,...] * U liefert skalar 
Gln_v_phi = [(M[2:3, -i - 2: -1] @ U[: i + 1])[0] - v for i, v in enumerate(V_phi)]

# Liste der Optimierungsvariablen

var = U + Up + Um + [umax] + V_theta

# Gleichungsbeschränkungen

HH = [Up[i] - Um[i] - u for i, u in enumerate(U)]
Aeq, beq = sp.linear_eq_to_matrix(HH + Gln_v_theta + Gln_xn, var)

# Gütefunktional

c, dummy = sp.linear_eq_to_matrix([umax], var)

# Ungleichungen für u_i

Ugln = [Up[i] + Um[i] - umax for i, u in enumerate(U)]
Auq, buq = sp.linear_eq_to_matrix(Ugln, var)

# Schranken der Optimierungsvariablen

Limit_u = 10 # Nm
Limit_v_theta = 100000 / 60 * 2 * np.pi   # rad/s
bounds = [(-Limit_u, Limit_u) for i in U] + [(0, None) for i in Up + Um] +\
            [(0, None)] + [(-Limit_v_theta, Limit_v_theta) for i in V_theta]

# Lineares Programm und Berechnung der optimalen Lösung

res = linprog(c, A_eq=Aeq, b_eq=beq, A_ub=Auq, b_ub=buq, bounds=bounds)

# Optimale Stellfolge und X-Folge

Uopt = res.x[0:N]
Xopt = Sol(np.array(sysd.A), np.array(sysd.B), Uopt, x0)

# In[10] Plotten

col = ['g^', 'r^', 'b^', 'm^', 'c^', 'y^', 'k^']

plt.figure()
plt.plot(range(N), Uopt, col[0], label='$u$')
plt.xlabel(f'n  ($T_a = ${str(T)})')
plt.ylabel('$u$')
plt.legend()
plt.title('$l_\\infty$' + '-optimal mit v Begrenzung')
plt.grid(True)
plt.xlim(0, N)
plt.show()

plt.figure()
for i in range(n//2):
    plt.plot(range(N + 1), Xopt[i], col[i + 1], label=f'$x_{i}$')
plt.xlabel('n  (' + '$T_a = $' + str(T) + ')')
plt.ylabel('X')
plt.legend()
plt.title('$l_\\infty$' + '-optimal mit v Begrenzung')
plt.grid(True)
plt.xlim(0, N)
plt.show()

plt.figure()
for i in range(n//2,n):
    plt.plot(range(N + 1), Xopt[i], col[i + 1], label=f'$x_{i}$')
plt.xlabel('n  (' + '$T_a = $' + str(T) + ')')
plt.ylabel('V')
plt.legend()
plt.title('$l_\\infty$' + '-optimal mit v Begrenzung')
plt.grid(True)
plt.xlim(0, N)
plt.show()
# %%
