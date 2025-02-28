""".

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

@author: Kurt
"""

# In[0]

from sys import stdout
from sympy import cse, diff, Matrix, Symbol, var

# In[1]


class EL_Vars ():
    """._."""

    def __init__(self, t, x):
        """.

        Parameters
        ----------
        t : Zeit
        x : [Variablennamen] 'x' nicht indiziert, ('x', n) indizier
        Returns
        -------
        None.
        """
        self.t = var(t)
        self.x = []
        self.v = []
        self.a = []
        for y in x:
            if isinstance(y, tuple):
                self.x += [Symbol(f'{y[0]}_{i}')
                           for i in range(1, y[1] + 1)]
                self.v += [Symbol(f'\\dot {y[0]}_{i}')
                           for i in range(1, y[1] + 1)]
                self.a += [Symbol(f'\\ddot {y[0]}_{i}')
                           for i in range(1, y[1] + 1)]
            else:
                self.x += [Symbol(f'{y}')]
                self.v += [Symbol(f'\\dot {y}')]
                self.a += [Symbol(f'\\ddot {y}')]

    def ttd(self, x):
        """.

        Parameters
        ----------
        x : Funktion der EL Variablen
        Returns
        -------
        res : totale Zeitableitung
        """
        res = diff(x, self.t)
        res += sum([i[1] * diff(x, i[0])
                    for i in list(zip(self.x, self.v))])
        res += sum([i[1] * diff(x, i[0])
                    for i in list(zip(self.v, self.a))])
        return res


def EL_eq(L, Vars):
    """.

    Parameters
    ----------
    L : EL Funktion
    Vars : EL Variablen
    Returns
    -------
    list : EL-Gleichungen
    """
    LL = Matrix([L])
    return LL.jacobian(Vars.v).applyfunc(Vars.ttd) - LL.jacobian(Vars.x)

# In[2]


def Prog_EL_py(M, R, Vars, x, f=stdout):
    """.

    Parameters
    ----------
    M : Massenmatric
    R : Rechte Seite
    Vars : EL-Variablen
    x : Name der Lage
    f : TYPE, Ausgabedatei
    Returns
    -------
    None.
    """
    n = len(Vars.x)
    sdic1 = {y: Symbol(f'x[{i}]') for i, y in enumerate(Vars.x)}
    sdic2 = {y: Symbol(f'x[{i + n}]') for i, y in enumerate(Vars.v)}
    sdic = sdic1 | sdic2

    print('from numpy import *\nfrom scipy import linalg\n', file=f)
    print(f'def EL_func (t, {x}):', file=f)
    print(f'    M = zeros(({n}, {n}))', file=f)
    print(f'    R = zeros({n})', file=f)
    for i in range(n):
        for j in range(i, n):
            print(f'    M[{i}, {j}] = {M[i,j].subs(sdic1)}', file=f)
    for i in range(0, n):
        print(f'    R[{i}] = {R[i].subs(sdic)}', file=f)
    print(f"""    return hstack((x[{n}:], \
linalg.solve(M, R, assume_a='pos')))""", file=f)


def Prog_EH_py(M, R, Vars, x, v, f=stdout):
    """.

    Parameters
    ----------
    M : Massenmatric
    R : Rechte Seite
    Vars : EL-Variablen
    x : Name des Lage
    v : Name der Lageableitung
    f : TYPE, Ausgabedatei
    Returns
    -------
    None.
    """
    n = len(Vars.x)
    sdic1 = {y: Symbol(f'{x}[{i}]') for i, y in enumerate(Vars.x)}
    sdic2 = {y: Symbol(f'{v}[{i}]') for i, y in enumerate(Vars.v)}
    sdic = sdic1 | sdic2
    print('from numpy import *\nfrom scipy import linalg\n', file=f)
    print(f'def EH_func (t, {x}):', file=f)
    print(f'    M = zeros(({n}, {n}))', file=f)
    print(f'    R = zeros({n})', file=f)
    for i in range(n):
        for j in range(i, n):
            print(f'    M[{i}, {j}] = {M[i,j].subs(sdic1)}', file=f)
    print(f"    {v} = linalg.solve(M, x[{n}:], assume_a='pos')", file=f)
    for i in range(0, n):
        print(f'    R[{i}] = {R[i].subs(sdic)}', file=f)
    print("    return hstack((V, R))",  file=f)

# In[3]


def Prog_EL_opt_py(M, R, Vars, x, f=stdout):
    """.

    Parameters
    ----------
    M : Massenmatric
    R : Rechte Seite
    Vars : EL-Variablen
    x : Name der Lage
    f : TYPE, Ausgabedatei
    Returns
    -------
    None.
    """
    n = len(Vars.x)
    sdic1 = {y: Symbol(f'x[{i}]') for i, y in enumerate(Vars.x)}
    sdic2 = {y: Symbol(f'x[{i + n}]') for i, y in enumerate(Vars.v)}
    sdic = sdic1 | sdic2

    cmds = cse([M, R], optimizations='basic')
    M1 = cmds[1][0]
    R1 = cmds[1][1]
    print('from numpy import *\nfrom scipy import linalg\n', file=f)
    print(f'def EL_func (t, {x}):', file=f)
    print(f'    M = zeros(({n}, {n}))', file=f)
    print(f'    R = zeros({n})', file=f)
    for i in cmds[0]:
        print(f'    {i[0]} = {i[1].subs(sdic)}', file=f)
    for i in range(n):
        for j in range(i, n):
            print(f'    M[{i}, {j}] = {M1[i,j].subs(sdic)}', file=f)
    for i in range(0, n):
        print(f'    R[{i}] = {R1[i].subs(sdic)}', file=f)
    print(f"""    return hstack((x[{n}:], \
linalg.solve(M, R, assume_a='pos')))""", file=f)


def Prog_EH_opt_py(M, R, Vars, x, v, f=stdout):
    """.

    Parameters
    ----------
    M : Massenmatric
    R : Rechte Seite
    Vars : EL-Variablen
    x : Name des Lage
    v : Name der Lageableitung
    f : TYPE, Ausgabedatei
    Returns
    -------
    None.
    """
    n = len(Vars.x)
    sdic1 = {y: Symbol(f'{x}[{i}]') for i, y in enumerate(Vars.x)}
    sdic2 = {y: Symbol(f'{v}[{i}]') for i, y in enumerate(Vars.v)}
    sdic = sdic1 | sdic2
    cmds = cse([M, R], optimizations='basic')
    Setv = set(Vars.v)
    cmds1 = []
    cmds2 = []
    for i in cmds[0]:
        if i[1].free_symbols.isdisjoint(Setv):
            cmds1.append(i)
        else:
            cmds2.append(i)
            Setv.add(i[0])
    M1 = cmds[1][0]
    R1 = cmds[1][1]
    print('from numpy import *\nfrom scipy import linalg\n', file=f)
    print(f'def EH_func (t, {x}):', file=f)
    print(f'    M = zeros(({n}, {n}))', file=f)
    print(f'    R = zeros({n})', file=f)
    for i in cmds1:
        print(f'    {i[0]} = {i[1].subs(sdic)}', file=f)
    for i in range(n):
        for j in range(i, n):
            print(f'    M[{i}, {j}] = {M1[i,j].subs(sdic)}', file=f)
    print(f"    {v} = linalg.solve(M, x[{n}:], assume_a='pos')", file=f)
    for i in cmds2:
        print(f'    {i[0]} = {i[1].subs(sdic)}', file=f)
    for i in range(0, n):
        print(f'    R[{i}] = {R1[i].subs(sdic)}', file=f)
    print("    return hstack((V, R))",  file=f)


