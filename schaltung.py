""".

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

@author: Kurt
"""

# In[0]

import sympy as sympy
from IPython.display import display
sympy.init_printing()

# In[1] Netzwerksverwaltung


class Netzwerk():
    """_Basisklasse_."""

    def __init__(self):
        self.Namen = dict()

    def Eleverbindungen(self):
        """_Verbinden der Komponenten_."""
        res = True
        Kanten = set()
        for x in self.Namen.values():
            if isinstance(x, Kante):
                Kanten.add(x)
        s = set()
        for x in Kanten:
            s.add(x.Kp)
            s.add(x.Km)
        for x in s:
            y = self.Namen.get(x)
            if not y:
                self.Namen[x] = Knoten(self, x)
            elif not isinstance(y, Knoten):
                res = False
                print(f'Knotename: {y.Name} ist bereits Vergeben')
        for x in Kanten:
            x.Kp = self.Namen[x.Kp]
            x.Km = self.Namen[x.Km]
        for x in self.Namen.values():
            if isinstance(x, (Iiquelle, Uuquelle)):
                y = self.Namen.get(x.In)
                if y and isinstance(y, Kante):
                    x.In = y
                else:
                    print(f'Verbindungsfehler: {x.Name}')
        return res

    def Gleichungen(self):
        """_Aufstellen des Gleichungssystem_."""
        self.aVars = set()
        for x in self.Namen.values():
            self.aVars = self.aVars.union(x.aVar())
        for x in self.Namen.values():
            if isinstance(x, Kante):
                x.Kp.Im.add(x.I)
                x.Km.Ip.add(x.I)
        self.GlSys = set()
        for x in self.Namen.values():
            self.GlSys = self.GlSys.union(x.Gleichung())
        self.Res = sympy.solve(self.GlSys, self.aVars)
        return self.Res

    def DynGL(self):
        """_Gleichungen der dynamischen Element_."""
        S = set()
        for x in self.Namen.values():
            if isinstance(x, (Induktivität, Kapazität)):
                S = S.union(x.DGl())
        return {x.subs(self.Res) for x in S}


class Objekt(Netzwerk):
    """_Subklasse der Netzwerkskomponenten_."""

    def __init__(self, Netzwerk, Name):
        """.

        Parameters
        ----------
        Netzwerk : Instanze Netzwerk
        Name : Name

        Returns
        -------
        None.
        """
        self.Name = Name
        x = Netzwerk.Namen.get(Name)
        if x:
            print(f'Objektname: {Name} ist bereits Vergeben')
        Netzwerk.Namen[Name] = self

# In[2] Knoten


class Knoten(Objekt):
    """_Klasse Knoten_."""

    def __init__(self, Netzwerk, Name):
        """.

        Parameters
        ----------
        Netzwerk : Instanz Netzwert
        Name :Name

        Returns
        -------
        None.
        """
        self.Ip = set()
        self.Im = set()
        super(Knoten, self).__init__(Netzwerk, Name)

    def aVar(self):
        """_Variablenname_."""
        self.aVar = sympy.Symbol(f'u_{{{str(self.Name)}}}')
        self.U = self.aVar
        return {self.aVar}

    def Gleichung(self):
        """_Gleichung_."""
        return {sum(self.Ip) - sum(self.Im)}


class KnotenP(Knoten):
    """_Klasse Knoten mit Potential_."""

    def __init__(self, Netzwerk, Name, Spannung='U'):
        """.

        Parameters
        ----------
        Netzwerk : Instanz Netzwerk
        Name : Name
        Spannung : Wert der Spannung

        Returns
        -------
        None.
        """
        if isinstance(Spannung, str):
            self.Spannung = sympy.Symbol(f'Spannung + {{{str(Name)}}}')
        else:
            self.Spannung = Spannung
        super(KnotenP, self).__init__(Netzwerk, Name)

    def aVar(self):
        """_Variablenname_."""
        self.aVar = sympy.Symbol(f'i_{{{str(self.Name)}}}')
        self.U = self.Spannung
        return {self.aVar}

    def Gleichung(self):
        """_Gleichung_."""
        return {self.aVar + sum(self.Ip) - sum(self.Im)}

# In[3] Kante und 1-Zore


class Kante(Objekt):
    """_Klasse Kante_."""

    def __init__(self, Netzwerk, Name, Kp, Km):
        """.

        Parameters
        ----------
        Netzwerk : Instanz Netzwerk
        Name : Name
        Kp : Knotenname +
        Km : Knotenname -

        Returns
        -------
        None.
        """
        self.Kp = Kp
        self.Km = Km
        super(Kante, self).__init__(Netzwerk, Name)


class Widerstand(Kante):
    """_Klasse Kante_."""

    def __init__(self, Netzwerk, Name, Kp, Km, R='R'):
        """.

        Parameters
        ----------
        Netzwerk : Instanz Netzwerk
        Name : Name
        Kp : Knotenname +
        Km : Knotenname -
        R : Wert des Widerstands

        Returns
        -------
        None.
        """
        if isinstance(R, str):
            self.R = sympy.Symbol(f'{R}_{{{str(Name)}}}')
        else:
            self.R = R
        super(Widerstand, self).__init__(Netzwerk, Name, Kp, Km)

    def aVar(self):
        """_Variablenname_."""
        self.aVar = sympy.Symbol(f'i_{{{str(self.Name)}}}')
        self.I = self.aVar
        return {self.aVar}

    def Gleichung(self):
        """_Gleichung_."""
        return {self.Kp.U - self.R*self.I - self.Km.U}


class Uquelle(Kante):
    """_Klasse Spannungsquelle_."""

    def __init__(self, Netzwerk, Name, Kp, Km, Uq='U'):
        """.

        Parameters
        ----------
        Netzwerk : Instanz Netzwerk
        Name : Name
        Kp : Knotenname +
        Km : Knotenname -
        Uq : Wert der Spannungsquelle

        Returns
        -------
        None.
        """
        if isinstance(Uq, str):
            self.Uq = sympy.Symbol(f'Uq_{{{str(Name)}}}')
        else:
            self.Uq = Uq
        super(Uquelle, self).__init__(Netzwerk, Name, Kp, Km)

    def aVar(self):
        """_Variablenname_."""
        self.aVar = sympy.Symbol(f'i_{{{str(self.Name)}}}')
        self.I = self.aVar
        return {self.aVar}

    def Gleichung(self):
        """_Gleichung_."""
        return {self.Kp.U + self.Uq - self.Km.U}


class Iquelle(Kante):
    """_Klasse Spannungsquelle_."""

    def __init__(self, Netzwerk, Name, Kp, Km, Iq='I'):
        """.

        Parameters
        ----------
        Netzwerk : Instanz Netzwerk
        Name : Name
        Kp : Knotenname +
        Km : Knotenname -
        Iq : Wert der Stromquelle

        Returns
        -------
        None.
        """
        if isinstance(Iq, str):
            self.Iq = sympy.Symbol(f'Iq_{{{str(Name)}}}')
        else:
            self.Iq = Iq
        super(Iquelle, self).__init__(Netzwerk, Name, Kp, Km)

    def aVar(self):
        """_Variablenname_."""
        self.I = self.Iq
        return {}

    def Gleichung(self):
        """_Gleichung_."""
        return {}

# in[4] Lineare dynamische 1-Tore


class Kapazität(Kante):
    """_Klasse Kapazität_."""

    def __init__(self, Netzwerk, Name, Kp, Km, C='C'):
        """.

        Parameters
        ----------
        Netzwerk : Instanz Netzwerk
        Name : Name
        Kp : Knotenname +
        Km : Knotenname -
        C : Wert der Kapazität

        Returns
        -------
        None.
        """
        if isinstance(C, str):
            self.C = sympy.Symbol(C + '_{' + str(Name) + '}')
        else:
            self.C = C
        super(Kapazität, self).__init__(Netzwerk, Name, Kp, Km)

    def aVar(self):
        """_Variablennamen_."""
        self.aVar = sympy.Symbol(f'i_{{{str(self.Name)}}}')
        self.sVar = sympy.Symbol(f'u_{{{str(self.Name)}}}')
        self.sdVar = sympy.Symbol(f'\\dot{{u}}_{{{str(self.Name)}}}')
        self.I = self.aVar
        return {self.aVar}

    def Gleichung(self):
        """_Gleichung_."""
        return {self.Kp.U - self.sVar - self.Km.U}

    def DGl(self):
        """_Diff-Gleichung_."""
        return {self.C*self.sdVar - self.aVar}


class Induktivität(Kante):
    """_Klasse Induktivität_."""

    def __init__(self, Netzwerk, Name, Kp, Km, L='L'):
        """.

        Parameters
        ----------
        Netzwerk : Instanz Netzwerk
        Name : Name
        Kp : Knotenname +
        Km : Knotenname -
        L : Wert der Induktivität

        Returns
        -------
        None.
        """
        if isinstance(L, str):
            self.L = sympy.Symbol(f'L_{{{str(Name)}}}')
        else:
            self.L = L
        super(Induktivität, self).__init__(Netzwerk, Name, Kp, Km)

    def aVar(self):
        """_Variablennamen_."""
        self.sVar = sympy.Symbol(f'i_{{{str(self.Name)}}}')
        self.sdVar = sympy.Symbol(f'\\dot{{i}}_{{{str(self.Name)}}}')
        self.I = self.sVar
        return {}

    def Gleichung(self):
        """_Gleichung_."""
        return {}

    def DGl(self):
        """_Diff-Gleichung_."""
        return {self.L*self.sdVar - self.Kp.U + self.Km.U}

# In[5] 2-Tore


class Uuquelle(Kante):
    """_Spannungsgesteuerte Spannungsquelle_."""

    def __init__(self, Netzwerk, Name, Kp, Km, In, ku='ku'):
        """.

        Parameters
        ----------
        Netzwerk : Instanz Netzwerk
        Name : Name
        Kp : Knotenname +
        Km : Knotenname -
        In : Name der steuernden Kante
        ku : Wert des Verstärkungsfaktors

        Returns
        -------
        None.
        """
        if isinstance(ku, str):
            self.ku = sympy.Symbol(f'ku_{{{str(Name)}}}')
        else:
            self.ku = ku
        self.In = In
        super(Uuquelle, self).__init__(Netzwerk, Name, Kp, Km)

    def aVar(self):
        """_Variablennamen_."""
        self.aVar = sympy.Symbol(f'i_{{{str(self.Name)}}}')
        self.I = self.aVar
        return {self.aVar}

    def Gleichung(self):
        """_Gleichung_."""
        return {self.Kp.U - self.Km.U -
                self.ku * (self.In.Kp.U - self.In.Km.U)}


class Iiquelle(Kante):
    """_Stromgesteuerte Stromquelle_."""

    def __init__(self, Netzwerk, Name, Kp, Km, In, ki='ki'):
        """.

        Parameters
        ----------
        Netzwerk : Instanz Netzwerk
        Name : Name
        Kp : Knotenname +
        Km : Knotenname -
        In : Name der steuernden Kante
        ki : Wert des Verstärkungsfaktors

        Returns
        -------
        None.
        """
        if isinstance(ki, str):
            self.ki = sympy.Symbol(f'ki_{{{str(Name)}}}')
        else:
            self.ki = ki
        self.In = In
        super(Iiquelle, self).__init__(Netzwerk, Name, Kp, Km)

    def aVar(self):
        """_Variablennamen_."""
        self.aVar = sympy.Symbol(f'I_{{{str(self.Name)}}}')
        self.I = self.aVar
        return {self.aVar}

    def Gleichung(self):
        """_Gleichung_."""
        return {self.I - self.ki * self.In.I}


class Iwandler(Objekt):
    """_Idealer Wandler_."""

    def __init__(self, Netzwerk, Name, Klp, Klm, Krp, Krm, ü='ü'):
        """.

        Parameters
        ----------
        Netzwerk : Instanz Netzwerk
        Name : Name
        Klp : Knotenname + links
        Klm : Knotenname - links
        Krp : Knotenname + rechts
        Krm : Knotenname - rechts
        ü : Wert der ÜbersetzungYPE, optional

        Returns
        -------
        None.
        """
        if isinstance(ü, str):
            self.ü = sympy.Symbol(f'ü _{{{str(Name)}}}')
        else:
            self.ü = ü
        NameL = '\\_L' + Name
        NameR = '\\_R' + Name
        self.Ltor = Iiquelle(Netzwerk, NameL, Klp, Klm, NameR, 1/self.ü)
        self.Rtor = Uuquelle(Netzwerk, NameR, Krp, Krm, NameL, self.ü)
        super(Iwandler, self).__init__(Netzwerk, Name)

    def aVar(self):
        """."""
        return {}

    def Gleichung(self):
        """."""
        return {}


class Iverstärker(Objekt):
    """_Idealer Operationsverstärker_."""

    def __init__(self, Netzwerk, Name, Klp, Klm, Krp, Krm):
        """.

        Parameters
        ----------
        Netzwerk : Instanz Netzwerk
        Name : Name
        Klp : Knotenname + links
        Klm : Knotenname - links
        Krp : Knotenname + rechts
        Krm : Knotenname - rechts
        ü : Wert der ÜbersetzungYPE, optional

        Returns
        -------
        None.
        """
        self.Ltor = Iquelle(Netzwerk, '\\_L' + Name, Klp, Klm, 0)
        self.Rtor = Uquelle(Netzwerk, '\\_R' + Name, Krp, Krm)
        super(Iverstärker, self).__init__(Netzwerk, Name)

    def aVar(self):
        """_Variablenname_."""
        return {self.Rtor.Uq}

    def Gleichung(self):
        """_Gleichung_."""
        return {self.Ltor.Kp.U - self.Ltor.Km.U}


# In[6]
from support import Sol
import control as ctrl

N = Netzwerk()
KnotenP(N, 'M', 0)
Knoten(N, 'E')
Uquelle(N, 'In', 'M', 'E')

Knoten(N, 'K')
Knoten(N, 'A')

Widerstand(N, 'R1', 'E', 'K')
Widerstand(N, 'R2', 'K', 'A')

Iverstärker(N, 'OPV', 'K', 'M', 'M', 'A')

# Kapazität(N, 'L1', 'X', 'K')
# Kapazität(N, 'L2', 'Y', 'A')

# Kapazität(N, 'C1', 'E', 'K')
# Kapazität(N, 'C2', 'K', 'A')
Induktivität(N, 'L1', 'E', 'K')
Induktivität(N, 'L2', 'K', 'A')

res = N.Eleverbindungen()
display(res)
# if res:
# display(N.Gleichungen())
# display(N.DynGL())

import sympy as sp
import control as ctrl
import support
if res:
    display(N.Gleichungen())
    gl = N.DynGL()
    display(gl)
    LL = ['L1', 'L2']
    l1 = [N.Namen[x].sVar for x in LL]
    l1d = [N.Namen[x].sdVar for x in LL]
    gli = sympy.solve(gl, l1d)
    display(gli)
    A, b = sympy.linear_eq_to_matrix([x.subs(gli) for x in l1d], l1)
    display(A)
    display(b)
    
    R2 = sp.Symbol("R_{R2}")
    R1 = sp.Symbol("R_{R1}")
    
    c = sp.Matrix([-R2, R2])
    d = sp.Matrix([-R2/ R1])
    
    subs = {R2:1000, R1: 100, sp.Symbol("L_{L1}"): 1e-3, 
						sp.Symbol("L_{L2}"):10e-3,
						sp.Symbol("Uq_{In}"):1}
    # sys = ctrl.StateSpace(A.subs(subs), b.subs(subs), c.subs(subs), d.subs(subs))
    sys = ctrl.ss(sp.matrix2numpy(A.subs(subs)),
									sp.matrix2numpy(b.subs(subs)),
									sp.matrix2numpy(c.subs(subs)),
									sp.matrix2numpy(d.subs(subs)))
		
    ctrl.bode_plot(sys)

