#!/usr/bin/env python
#
# @purpose
#   a collection of physical constants
#
# @Author Peng Lian
# @Email penglian518@gmail.com
# @Date Sep 30 2016
#

# Physical constants, derived quantities, conversion factors, etc.

from math import *

pi      = acos(-1.0)
e       = exp(1.0)
clight  = 2.99792458e8          # Speed of light (m s^-1)
hartree = 2.0*2.1798736e-18     # J (atomic unit of energy = 2*Rydberg constant)
# M. Nooijen values
h       = 6.6260755e-34         # Planck's constant (J s) MN
kb      = 1.38064852e-23          # Boltzmann constant (J K^-1) MN
a0      = 5.29177249e-11        # Bohr radius (m) MN
me      = 9.1093897e-31         # kg MN
amu     = 1.6605402e-27         # kg MN
mp      = 1.6726231e-27         # kg MN
eV      = 1.60217733e-19        # J MN
h2kcal  = 627.5095              # Hartree to kcal/mol
cal     = 4.184                 # 1 cal = 4.184 J
kcal    = 4184.0                #
#NA      = 6.0221367e23          # Avogadro's number MN
NA      = 6.022140857e23        #

hbar    = h/(2.0*pi)
kcm     = kb/(h*100.*clight)    # speed of light in wave numbers (cm^-1)
R       = NA*kb                 # gas constant (J K^-1 mol^-1)
T       = 298.15                # temperature in K
RT      = kb*NA*T/kcal          # gas constant in kcal/mol
P       = 101325.0              # atmospheric pressure in Pa ( 101325 Pa = 1 atm = 760 mmHg )

faraday = 23.061                # kcal mol^-1 eV^-1
SCE = 4.52                      # Reference electrode (Saturated Calomel Electrode)

kb_kJmolK = kb * NA / 1000.     # Boltzmann constant in kJ/(mol*K)
kb_kcalmolK = kb_kJmolK/cal     # Boltzmann constant in kcal/(mol*K)

#gas2star = 0.0
#gas2star = RT*log(24.46/0.002)                    # Gas to liquid standard state correction: 1 mol/atm --> 1 mol/L = +1.89 kcal/mol
gas2star = RT*log(24.46/1.0)                    # Gas to liquid standard state correction: 1 mol/atm --> 1 mol/L = +1.89 kcal/mol
liquid2star = 1.0*RT*log(55.34/1.0)         # Standard state correction for ion solvation = 2.38 kcal/mol for n_water = 1
#liquid2star = 0.0         # Standard state correction for ion solvation = 2.38 kcal/mol for n_water = 1

G0_proton = -6.28/h2kcal                    # Gas-phase free energy of H+ (Sackur-Tetrode entropy plus translational enthalpy at 298 K at 1 atm)
Gs_proton = (-6.28 + gas2star)/h2kcal       # Gas-phase free energy of H+ (Sackur-Tetrode entropy plus translational enthalpy at 298 K at 1 M)
G0_electron = -0.868/h2kcal                 # Gas-phase free energy of e- at 1 atm (from Fermi-Dirac statistical mechanics, see Bartmess, JPC, 1994, 98, 6420)
Gs_electron = (-0.868 + gas2star)/h2kcal    # Gas-phase free energy of e- at 1 M (from Fermi-Dirac statistical mechanics, see Bartmess, JPC, 1994, 98, 6420)
dG0solv_proton = -264.0/h2kcal              # 1 mol/atm solvation free energy of H+ at 298 K (see Tissandier, 1998)
dGssolv_proton = (-264.0 - gas2star)/h2kcal # 1 mol/L solvation free energy of H+ at 298 K (see Tissandier, 1998)
#dG0solv_electron = -35.5/h2kcal              # 1 mol/atm solvation free energy of e- at 298 K (see Zhang, JPCB, 2003, 107, 4403)
#dGssolv_electron = (-35.5 - gas2star)/h2kcal # 1 mol/K solvation free energy of e- at 298 K (see Zhang, JPCB, 2003, 107, 4403)
