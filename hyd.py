#!/usr/bin/env python3
"""
hyd.py  -  steady-state, pseudo-homogeneous, 1-D non-isothermal model of the laboratory
steam-methane-reforming (SMR) reactor of

    Costamagna et al., "Modeling of Laboratory Steam Methane Reforming and CO2 Methanation
    Reactors", Energies 2020, 13, 2624.

REACTOR CHAIN (simulated in this order, each section feeds the next)
    PBR   12 mm  catalytic packed bed inside the furnace (kinetics + energy balance)
    P-R1  50 mm  empty glass tube inside the furnace (gas-phase water-gas-shift only)
    P-R2  40 mm  empty AISI-316L pipe, wall held at 523 K (wall-catalysed WGS, gas cools down)
  Paper Fig. 4 shows only the PBR and the first 3 mm of P-R2; P-R1 is simulated but not drawn.

MODEL
    Kinetics : Xu & Froment (1989): SMR, WGS, GRR (Langmuir-Hinshelwood), catalyst-specific
               adjustment factors zeta1..zeta6 (paper Table 3).
    Balances : dF_i/dV = sum_j nu_ij r_j                              (species)
               sum_i F_i cp_i dT/dV = sum_j (-dH_j) r_j + U a (Tw-T)  (energy)

HOW TO RUN  (keep hyd.py, validation.py, hyd_contours.py, compare_h_length.py and the data/ folder together)
    python hyd.py                                  # 3 catalysts, 853 K, 120 mL/min -> plots + csv
    python hyd.py --catalyst HMMC --T 853 --flow 120
    python hyd.py --seed-study                     # H2-seed sensitivity (numerical regularisation only)
    python hyd.py --validate                       # also runs validation.py (needs data/*.csv)
    python hyd.py --h-length particle              # sensitivity: particle diameter in the PBR Nusselt number
  All figures/CSVs go to --out (default: current folder).  Full instructions: README.md
"""
# ---- imports -------------------------------------------------------------------------------
import argparse                      # command-line options in main()
from dataclasses import dataclass    # Reactor geometry container

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

# ---- physical constants and species bookkeeping ---------------------------------------------
R = 8.314462618            # J/mol/K
P_BAR = 1.01325            # total pressure [bar]  (ASSUMPTION: 1 atm, no pressure drop)
P_PA = P_BAR * 1e5

SPECIES = ["CH4", "H2O", "CO", "CO2", "H2", "He"]
iCH4, iH2O, iCO, iCO2, iH2, iHe = range(6)
MW = np.array([16.043, 18.015, 28.010, 44.010, 2.016, 4.0026])   # g/mol

# Stoichiometric matrix NU: rows = species, columns = reactions [SMR, WGS, GRR]
#   SMR: CH4 + H2O  <-> CO  + 3 H2     WGS: CO + H2O <-> CO2 + H2     GRR: CH4 + 2 H2O <-> CO2 + 4 H2
NU = np.zeros((6, 3))
NU[[iCH4, iH2O, iCO, iH2], 0] = [-1, -1, 1, 3]
NU[[iCO, iH2O, iCO2, iH2], 1] = [-1, -1, 1, 1]
NU[[iCH4, iH2O, iCO2, iH2], 2] = [-1, -2, 1, 4]

# ----------------------------------------------------------------------------
# Thermodynamics: NIST Shomate polynomials  (t = T/1000)
#   cp [J/mol/K] = A + B t + C t^2 + D t^3 + E / t^2
#   Hf(T) [kJ/mol] = A t + B t^2/2 + C t^3/3 + D t^4/4 - E/t + F
# each entry: (Tmax, A, B, C, D, E, F)
# ----------------------------------------------------------------------------
SHOMATE = {
    "CH4": [(1300, -0.703029, 108.4773, -42.52157, 5.862788, 0.678565, -76.84376),
            (6000, 85.81217, 11.26467, -2.114146, 0.138190, -26.42221, -153.5327)],
    "H2O": [(1700, 30.09200, 6.832514, 6.793435, -2.534480, 0.082139, -250.8810),
            (6000, 41.96426, 8.622053, -1.499780, 0.098119, -11.15764, -272.1797)],
    "CO":  [(1300, 25.56759, 6.096130, 4.054656, -2.671301, 0.131021, -118.0089),
            (6000, 35.15070, 1.300095, -0.205921, 0.013550, -3.282780, -127.8375)],
    "CO2": [(1200, 24.99735, 55.18696, -33.69137, 7.948387, -0.136638, -403.6075),
            (6000, 58.16639, 2.720074, -0.492289, 0.038844, -6.447293, -425.9186)],
    "H2":  [(1000, 33.066178, -11.363417, 11.432816, -2.772874, -0.158558, -9.980797),
            (2500, 18.563083, 12.257357, -2.859786, 0.268238, 1.977990, -1.147438)],
    "He":  [(6000, 20.786, 0.0, 0.0, 0.0, 0.0, -20.786 * 0.29815)],
}
_SH = [SHOMATE[s] for s in SPECIES]


def _coef(i, T):
    """pick the Shomate coefficient row valid at temperature T for species i"""
    for row in _SH[i]:
        if T <= row[0]:
            return row[1:]
    return _SH[i][-1][1:]


def cp_i(T):
    """molar heat capacities [J/mol/K], array over species"""
    t = T / 1000.0                       # Shomate polynomials use t = T/1000
    out = np.empty(6)
    for i in range(6):                   # loop over the six species
        A, B, C, D, E, F = _coef(i, T)
        out[i] = A + B * t + C * t**2 + D * t**3 + E / t**2
    return out


def hf_i(T):
    """enthalpy of formation at T [J/mol], array over species"""
    t = T / 1000.0
    out = np.empty(6)
    for i in range(6):
        A, B, C, D, E, F = _coef(i, T)
        out[i] = 1000.0 * (A * t + B * t**2 / 2 + C * t**3 / 3 + D * t**4 / 4 - E / t + F)
    return out


def dH_rxn(T):
    """reaction enthalpies [J/mol] for [SMR, WGS, GRR] at temperature T"""
    return NU.T @ hf_i(T)                # dH_rxn = sum_i nu_i * Hf_i  (one value per reaction)


# Equilibrium constants.  Paper Table 4 (A*exp(-EK/RT)); in the PDF the two EK
# values appear swapped between columns, the physically consistent assignment
# is used here (checked against Xu-Froment K1 = exp(-26830/T + 30.114) bar^2).
ATM2BAR2 = 1.01325**2        # converts a constant in atm^2 to bar^2


KP_MODEL = "D"   # equilibrium-constant model: A literal Table 4 | B Pa^2 reading | D Table 4 with E-row swapped (default) | C Hou & Hughes


def Kp_SMR(T):
    """SMR equilibrium constant [bar^2]; the branch is chosen by KP_MODEL (default D)"""
    if KP_MODEL == "A":      # literal Table 4: A=7.846e12 atm^2, E=-37.72 kJ/mol
        return 7.846e12 * np.exp(37.72e3 / (R * T)) * ATM2BAR2
    if KP_MODEL == "B":      # literal Table 4 with the Pa^2 reading of the unit (x1e-10 -> bar^2)
        return 7.846e12 * np.exp(37.72e3 / (R * T)) * 1e-10
    if KP_MODEL == "C":      # Hou & Hughes (2001), bar^2
        return np.exp(-26830.0 / T + 30.114)
    return 7.846e12 * np.exp(-220.2e3 / (R * T)) * ATM2BAR2     # model D, bar^2


def Kp_WGS(T):
    """WGS equilibrium constant [-]; same KP_MODEL switch"""
    if KP_MODEL in ("A", "B"):   # literal Table 4: A=1.412e-2, E=220.2 kJ/mol
        return 1.412e-2 * np.exp(-220.2e3 / (R * T))
    if KP_MODEL == "C":
        return np.exp(4400.0 / T - 4.036)
    return 1.412e-2 * np.exp(37.72e3 / (R * T))                  # model D


def Kp_GRR(T):
    """GRR equilibrium constant = Kp_SMR * Kp_WGS (paper Eq. 33)"""
    return Kp_SMR(T) * Kp_WGS(T)                                 # bar^2


# ----------------------------------------------------------------------------
# Transport properties (Chapman-Enskog + modified Eucken, Wilke mixing)
# ----------------------------------------------------------------------------
# Lennard-Jones parameters used by the Chapman-Enskog viscosity formula
LJ = {  # sigma [A], eps/k [K]
    "CH4": (3.758, 148.6), "H2O": (2.641, 809.1), "CO": (3.690, 91.7),
    "CO2": (3.941, 195.2), "H2": (2.827, 59.7), "He": (2.551, 10.22),
}


def _omega(Ts):
    """collision integral Omega_mu(T*) (Neufeld correlation)"""
    return 1.16145 * Ts**-0.14874 + 0.52487 * np.exp(-0.7732 * Ts) + 2.16178 * np.exp(-2.43787 * Ts)


def mu_pure(T):
    """pure-species dynamic viscosities [Pa s] (Chapman-Enskog), array over species"""
    out = np.empty(6)
    for i, s in enumerate(SPECIES):
        sig, ek = LJ[s]
        out[i] = 2.6693e-6 * np.sqrt(MW[i] * T) / (sig**2 * _omega(T / ek))
    return out


def k_pure(T):
    """pure-species thermal conductivities [W/m/K] (modified Eucken: k = mu/M (1.32 cv + 1.77 R))"""
    mu = mu_pure(T)
    cv = cp_i(T) - R
    return mu / (MW * 1e-3) * (1.32 * cv + 1.77 * R)


def mix_props(y, T):
    """returns mu [Pa s], kf [W/m/K] (mole-weighted, as in paper), cp_mass [J/kg/K], rho [kg/m3]"""
    y = np.clip(y, 0, None)
    y = y / y.sum()
    mu_i = mu_pure(T)
    # Wilke mixing rule for the mixture viscosity
    phi = np.empty((6, 6))
    for i in range(6):
        for j in range(6):
            phi[i, j] = (1 + np.sqrt(mu_i[i] / mu_i[j]) * (MW[j] / MW[i])**0.25)**2 / \
                        np.sqrt(8 * (1 + MW[i] / MW[j]))
    mu = np.sum(y * mu_i / (phi @ y))                 # mixture viscosity
    kf = np.sum(y * k_pure(T))                         # mole-weighted conductivity (as stated in the paper)
    Mmix = np.sum(y * MW)                              # mean molar mass [g/mol]
    cp_mass = np.sum(y * cp_i(T)) / (Mmix * 1e-3)      # mass heat capacity [J/kg/K]
    rho = P_PA * Mmix * 1e-3 / (R * T)                 # ideal-gas density [kg/m3]
    return mu, kf, cp_mass, rho


# solid (silica glass) conductivity: quadratic fit through the ks values of paper Table 5
_ks_fit = np.polyfit([473.0, 773.0, 1173.0], [1.57, 1.99, 3.39], 2)


def k_solid(T):
    """silica-glass conductivity [W/m/K] at temperature T"""
    return float(np.polyval(_ks_fit, T))


def k_radial_ZS(kf, ks, eps=0.4):
    """Zehner-Schlunder effective radial conductivity of a packed bed (spheres)."""
    kap = ks / kf                                     # solid/gas conductivity ratio kappa
    B = 1.25 * ((1 - eps) / eps)**(10.0 / 9.0)
    s = np.sqrt(1 - eps)                              # sqrt(1-eps)
    d = 1 - B / kap
    # bracket term of the Zehner-Schlunder equation (paper Eq. 14); the -(B-1)/2 term was a sign bug before
    br = ((1 - 1 / kap) * B / d**2) * np.log(kap / B) - (B - 1) / d - (B - 1) / 2
    return kf * (1 - s + 2 * s / d * br)


def k_axial_votruba(kf, ks, eps=0.4):
    """auxiliary reporting calculation only (compare with paper Table 5); the 1-D energy balance has no axial
    conduction/dispersion term, as in the paper"""
    kap = ks / kf
    return kf * kap**(0.28 - 0.757 * np.log(eps) - 0.057 * np.log(kap))   # natural log (reproduces Table 5)


# ----------------------------------------------------------------------------
# Kinetics (Xu & Froment) with adjustment factors (paper Tables 2 and 3)
# ----------------------------------------------------------------------------
# Catalyst adjustment factors zeta1..zeta6 (paper Table 3). z1-z3 scale the rate constants of SMR/WGS/GRR,
# z4-z5 scale the activation energies of SMR/WGS, z6 scales the water adsorption constant.
ZETA = {  # zeta1..zeta6
    "HMMC":  dict(z1=15, z2=15,  z3=15, z4=1.00, z5=1.0, z6=1),
    "57-4Q": dict(z1=25, z2=100, z3=10, z4=1.17, z5=1.5, z6=90),
    "25-4Q": dict(z1=8,  z2=130, z3=10, z4=1.22, z5=1.7, z6=95),
    "base":  dict(z1=1,  z2=1,   z3=1,  z4=1.00, z5=1.0, z6=1),    # untouched Xu-Froment
}


def _k_arr(kTr, E, T, Tr):
    """re-parameterised Arrhenius law: k(T) = k(Tr) exp(-E/R (1/T - 1/Tr))"""
    return kTr * np.exp(-E / R * (1 / T - 1 / Tr))


def rates_catalyst(p, T, z):
    """p: partial pressures [bar] array over SPECIES. Returns [rSMR, rWGS, rGRR] in kmol/kgcat/h."""
    pCH4, pH2O, pCO, pCO2, pH2 = p[iCH4], p[iH2O], p[iCO], p[iCO2], max(p[iH2], 1e-12)   # floor avoids 0 division
    # rate constants of the three reactions, scaled by the catalyst factors
    kS = z["z1"] * _k_arr(1.842e-4, 240.1e3 * z["z4"], T, 648.0)
    kW = z["z2"] * _k_arr(7.558,    67.13e3 * z["z5"], T, 648.0)
    kG = z["z3"] * _k_arr(2.193e-5, 243.9e3,           T, 648.0)
    # adsorption equilibrium constants (van't Hoff form) of CH4, CO, H2 and H2O
    KCH4 = _k_arr(0.1791, 0.0, T, 823.0) * np.exp(38.28e3 / R * (1 / T - 1 / 823.0))
    KCO = 40.91 * np.exp(70.65e3 / R * (1 / T - 1 / 648.0))
    KH2 = 0.02960 * np.exp(82.90e3 / R * (1 / T - 1 / 648.0))
    KH2O = z["z6"] * 0.4152 * np.exp(-88.68e3 / R * (1 / T - 1 / 823.0))
    DEN = 1 + KCO * pCO + KH2 * pH2 + KCH4 * pCH4 + KH2O * pH2O / pH2   # Langmuir-Hinshelwood denominator
    # Xu-Froment rate expressions: forward driving force minus reverse term (p^n / Kp)
    rS = kS / pH2**2.5 * (pCH4 * pH2O - pH2**3 * pCO / Kp_SMR(T)) / DEN**2
    rW = kW / pH2 * (pCO * pH2O - pH2 * pCO2 / Kp_WGS(T)) / DEN**2
    rG = kG / pH2**3.5 * (pCH4 * pH2O**2 - pH2**4 * pCO2 / Kp_GRR(T)) / DEN**2
    return np.array([rS, rW, rG])


def rate_wgs_wall(p, T, kWGS_Tr=4.7):
    """WGS on metal pipe wall: Xu-Froment form, kWGS in kmol/bar/m2(inner surface)/h  -> kmol/m2/h"""
    pCH4, pH2O, pCO, pCO2, pH2 = p[iCH4], p[iH2O], p[iCO], p[iCO2], max(p[iH2], 1e-12)
    kW = _k_arr(kWGS_Tr, 67.13e3, T, 648.0)          # wall rate constant, k_WGS = 4.7 at 648 K (paper, area basis)
    KCH4 = 0.1791 * np.exp(38.28e3 / R * (1 / T - 1 / 823.0))
    KCO = 40.91 * np.exp(70.65e3 / R * (1 / T - 1 / 648.0))
    KH2 = 0.02960 * np.exp(82.90e3 / R * (1 / T - 1 / 648.0))
    KH2O = 0.4152 * np.exp(-88.68e3 / R * (1 / T - 1 / 823.0))
    DEN = 1 + KCO * pCO + KH2 * pH2 + KCH4 * pCH4 + KH2O * pH2O / pH2
    return kW / pH2 * (pCO * pH2O - pH2 * pCO2 / Kp_WGS(T)) / DEN**2


def rate_wgs_gas(y, T):
    """homogeneous gas-phase WGS (Bustamante), returns mol/m3/s"""
    # concentrations of the ideal gas, and the Bustamante rate law  r = k sqrt(C_CO) C_H2O (1 - approach to equilibrium)
    C = np.clip(y, 1e-14, None) * P_PA / (R * T)          # mol/m3
    Ccm = C * 1e-6                                          # mol/cm3
    k = 7.4e11 * np.exp(-288.3e3 / (R * T))                 # (cm3/mol)^0.5 / s
    pp = np.clip(y, 1e-14, None) * P_BAR                    # partial pressures [bar]
    approach = 1 - pp[iH2] * pp[iCO2] / (pp[iCO] * pp[iH2O] * Kp_WGS(T))
    return k * np.sqrt(Ccm[iCO]) * Ccm[iH2O] * approach * 1e6


# ----------------------------------------------------------------------------
# Geometry & operating conditions
# ----------------------------------------------------------------------------
# (Re characteristic length, Nu-to-h length, Pr conductivity).
# The legacy h_length field remains effective when transport_config is None.
TRANSPORT_CONFIGS = {
    "v3_reproduction": ("tube", "tube", "gas"),
    "v3_h_particle_only": ("tube", "particle", "gas"),
    "particle_lengths_gas_pr": ("particle", "particle", "gas"),
    "printed_equations": ("particle", "particle", "bed"),
}


@dataclass
class Reactor:
    """geometry and options of the laboratory reactor (SI units); defaults = paper values"""
    d: float = 6e-3          # inner diameter [m]
    s: float = 1e-3          # wall thickness [m]
    L_pbr: float = 12e-3
    L_pr1: float = 50e-3
    L_pr2: float = 40e-3
    m_cat: float = 88.2e-6   # kg
    kw_glass: float = 1.6
    kw_steel: float = 16.3   # AISI 316L
    T_pr2_wall: float = 523.0
    h_length: str = "tube"   # characteristic length in Nu->h: 'tube' or 'particle'
    dp: float = 0.2e-3
    U_scale: float = 1.0     # sensitivity knob: multiply PBR U
    transport_config: str = None  # explicit named configuration overrides legacy h_length


def effective_transport(r):
    """Resolve named or legacy options. Lengths [m]; conductivity choice is explicit."""
    name = r.transport_config
    if name is None:
        if r.h_length not in ("tube", "particle"):
            raise ValueError("h_length must be tube or particle")
        name = "v3_reproduction" if r.h_length == "tube" else "v3_h_particle_only"
    if name not in TRANSPORT_CONFIGS:
        raise ValueError(f"Unknown transport configuration: {name}")
    re_length, h_length, pr_conductivity = TRANSPORT_CONFIGS[name]
    return dict(configuration=name, Re_length_m=r.d if re_length == "tube" else r.dp,
                Nu_length_m=r.d if h_length == "tube" else r.dp,
                Pr_conductivity=pr_conductivity)


def transport_from_properties(r, *, rho, us, mu, cp, kf, ks):
    """PBR closure for supplied SI properties; supports tabulated-property diagnostics.

    rho [kg/m3], us [m/s], mu [Pa s], cp [J/kg/K], kf/ks [W/m/K].
    Returns U [W/m2/K] and local properties, definitions and dimensionless groups.
    """
    cfg = effective_transport(r)
    kr = k_radial_ZS(kf, ks)
    k_pr = kf if cfg["Pr_conductivity"] == "gas" else kr
    Re = rho * us * cfg["Re_length_m"] / mu
    Pr = cp * mu / k_pr
    Nu = 2.67 + 0.53 * Re**0.77 * Pr**0.53
    h = Nu * kr / cfg["Nu_length_m"]
    U = r.U_scale / (1 / h + r.s / r.kw_glass)
    return U, dict(Re=Re, Pr=Pr, Nu=Nu, kr=kr, kf=kf, us=us, mu=mu,
                   rho=rho, cp=cp, h=h, ks=ks, k_Pr_W_mK=k_pr,
                   kw_W_mK=r.kw_glass, U_W_m2K=U, **cfg)


SMR_REACTOR = Reactor()      # SMR mode: furnace + 523 K P-R2 wall
# MCO2_REACTOR: geometry of the paper's CO2-methanation reactor, retained for completeness; NOT used in this SMR work
MCO2_REACTOR = Reactor(L_pbr=26e-3, L_pr1=150e-3, L_pr2=300e-3, T_pr2_wall=298.0)


def area(r):
    """tube cross-section [m2]"""
    return np.pi * r.d**2 / 4


def U_pbr(F, T, r):
    """overall PBR-to-furnace heat-transfer coefficient U [W/m2/K] from Nu (Singhal), k_r (Zehner-Schlunder) and wall"""
    y = F / F.sum()
    mu, kf, cpm, rho = mix_props(y, T)
    Ftot = F.sum()
    us = Ftot * R * T / (P_PA * area(r))
    return transport_from_properties(r, rho=rho, us=us, mu=mu, cp=cpm,
                                     kf=kf, ks=k_solid(T))


def U_pipe(F, T, r, kw):
    """U for the empty tubes P-R1/P-R2: laminar flow, Nu = 3.66, in series with the wall resistance s/kw"""
    y = F / F.sum()
    mu, kf, cpm, rho = mix_props(y, T)
    h = 3.66 * kf / r.d                    # laminar, constant wall T
    return 1 / (1 / h + r.s / kw)


# ----------------------------------------------------------------------------
# ODE right-hand sides   state = [F_CH4, F_H2O, F_CO, F_CO2, F_H2, F_He, T]
# ----------------------------------------------------------------------------
def rhs_pbr(V, state, r, Tf, z):
    """ODE right-hand side of the PBR: species molar flows and temperature versus bed volume V"""
    F = np.clip(state[:6], 0, None) + 1e-30          # molar flows [mol/s] (clipped: no negative flows)
    T = state[6]
    y = F / F.sum()                                    # mole fractions
    p = y * P_BAR                                      # partial pressures [bar]
    rho_bed = r.m_cat / (area(r) * r.L_pbr)                     # kg_cat / m3_bed
    rj = rates_catalyst(p, T, z) * rho_bed * 1000.0 / 3600.0   # mol/m3/s
    dF = NU @ rj                                       # species balance: dF_i/dV = sum_j nu_ij r_j
    U, _ = U_pbr(F, T, r)
    a = 4.0 / r.d                                      # heat-exchange area per volume of tube
    q = -(dH_rxn(T) @ rj) + U * a * (Tf - T)           # reaction heat + exchange with the furnace [W/m3]
    dT = q / np.sum(F * cp_i(T))                       # energy balance
    return np.append(dF, dT)


def rhs_pr1(V, state, r, Tf):
    """ODE right-hand side of P-R1 (empty glass tube in the furnace): gas-phase WGS only, wall at furnace T"""
    F = np.clip(state[:6], 0, None) + 1e-30
    T = state[6]
    y = F / F.sum()
    rw = rate_wgs_gas(y, T)                            # homogeneous WGS rate [mol/m3/s]
    dF = NU[:, 1] * rw                                 # only the WGS column of NU contributes
    U = U_pipe(F, T, r, r.kw_glass)                    # glass wall
    a = 4.0 / r.d
    q = -(dH_rxn(T)[1] * rw) + U * a * (Tf - T)        # WGS heat + exchange with furnace
    return np.append(dF, q / np.sum(F * cp_i(T)))


def rhs_pr2(V, state, r, Tw):
    """ODE right-hand side of P-R2 (steel pipe, wall temperature Tw = 523 K): wall-catalysed WGS + cooling"""
    F = np.clip(state[:6], 0, None) + 1e-30
    T = state[6]
    y = F / F.sum()
    p = y * P_BAR
    a = 4.0 / r.d
    rw = rate_wgs_wall(p, T) * 1000.0 / 3600.0 * a            # mol/m3/s (surface rate x area/volume)
    dF = NU[:, 1] * rw                                 # WGS stoichiometry
    U = U_pipe(F, T, r, r.kw_steel)                    # steel wall
    q = -(dH_rxn(T)[1] * rw) + U * a * (Tw - T)        # gas is cooled toward the 523 K wall
    return np.append(dF, q / np.sum(F * cp_i(T)))


# ----------------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------------
def feed_moles(flow_ml_min, y_feed):          # block: convert a volumetric flow at RTP into molar flows
    """molar flows [mol/s] from flow at RTP (273.15 K, 101.325 kPa)"""
    ntot = 101325.0 * (flow_ml_min * 1e-6 / 60.0) / (R * 273.15)
    return ntot * np.asarray(y_feed)


# numerical settings (the paper used MATLAB ode45; LSODA is used here for robustness)
SOLVER = "LSODA"     # LSODA (default) | RK45 | BDF | Radau
TOL_SCALE = 1.0      # multiplies rtol/atol (0.5 for the tolerance-halving check)


def integrate(rhs, state0, L, args, n=400):
    """integrate one reactor section over its length L; returns z [m] and the 7 states on n points"""
    A = args[0]
    Vend = area(A) * L                                 # integrate over volume, convert back to length below
    sol = solve_ivp(rhs, (0, Vend), state0, args=args, method=SOLVER,
                    rtol=1e-8 * TOL_SCALE, atol=np.r_[np.full(6, 1e-16), 1e-6] * TOL_SCALE,
                    t_eval=np.linspace(0, Vend, n))
    if not sol.success:                                # never continue silently after a solver failure
        raise RuntimeError(sol.message)
    z = sol.t / area(A)
    return z, sol.y


def run(catalyst="HMMC", Tf=853.0, flow=120.0, y_in=(0.05, 0.20, 0, 0, 0, 0.75),
        reactor=SMR_REACTOR, h2_seed=1e-6, zeta=None, n=400):
    """Integrate PBR -> P-R1 -> P-R2. Returns dict with profiles per section."""
    # feed composition (5% CH4, 20% H2O, 75% He); a 1e-6 H2 seed removes the p_H2 = 0 singularity of the
    # Xu-Froment rates. It is a NUMERICAL regularisation, not a physical feed component.
    y0 = np.array(y_in, float)
    y0[iH2] = max(y0[iH2], h2_seed)           # avoid the p_H2 = 0 singularity of Xu-Froment
    y0 /= y0.sum()
    F0 = feed_moles(flow, y0)                          # inlet molar flows [mol/s]
    state0 = np.append(F0, Tf)                         # gas enters at furnace temperature
    zz = ZETA[catalyst] if zeta is None else zeta      # catalyst adjustment factors
    out = {"F0": F0}

    # 1) PBR: kinetics + energy balance, bed wall at furnace temperature Tf
    zb, yb = integrate(rhs_pbr, state0, reactor.L_pbr, (reactor, Tf, zz), n)
    out["PBR"] = (zb, yb)
    # 2) P-R1: starts from the PBR outlet; gas-phase WGS only (simulated, but NOT drawn in Fig. 4)
    z1, y1 = integrate(rhs_pr1, yb[:, -1], reactor.L_pr1, (reactor, Tf), n)
    out["PR1"] = (z1, y1)
    # 3) P-R2: full 40 mm at 523 K wall; 2000 points resolve the steep cooling in the first mm
    z2, y2 = integrate(rhs_pr2, y1[:, -1], reactor.L_pr2, (reactor, reactor.T_pr2_wall), 2000)
    out["PR2"] = (z2, y2)
    return out


def mole_fractions(Y):
    """molar flows (rows 0-5 of Y) -> mole fractions, column by column"""
    F = Y[:6]
    return F / F.sum(axis=0)


def summarize(res):
    """exit compositions/temperatures of each section, CH4 conversion and atom-balance errors"""
    rows = {}
    for sec in ("PBR", "PR1", "PR2"):
        z, Y = res[sec]
        y = mole_fractions(Y)[:, -1]
        rows[sec] = dict(y=y, T=Y[6, -1])
    F0 = res["F0"]
    Fe = res["PR2"][1][:6, -1]
    rows["X_CH4"] = 1 - Fe[iCH4] / F0[iCH4]
    rows["atom_err"] = atom_balance(F0, Fe)
    return rows


def atom_balance(F0, Fe):
    """relative C/H/O atom-balance error between inlet and outlet flows (should be ~1e-15)"""
    C = np.array([1, 0, 1, 1, 0, 0]); H = np.array([4, 2, 0, 0, 2, 0]); O = np.array([0, 1, 1, 2, 0, 0])
    return {n: float((v @ Fe - v @ F0) / (v @ F0)) for n, v in (("C", C), ("H", H), ("O", O))}


# ----------------------------------------------------------------------------
# Thermodynamic equilibrium (SMR + WGS) for validation
# ----------------------------------------------------------------------------
def equilibrium(y_in, T):
    """equilibrium composition of the SMR + WGS system at temperature T (1 atm), by least squares on the two
    extents of reaction; used to check that HMMC reaches equilibrium within ~2 mm"""
    n0 = np.asarray(y_in, float)

    def comp(x):
        n = n0 + NU[:, 0] * x[0] + NU[:, 1] * x[1]
        return n

    def res(x):
        n = comp(x)
        if np.any(n[:5] <= 0):
            return np.array([1e3, 1e3])
        y = n / n.sum()
        p = y * P_BAR
        return np.array([np.log(p[iCO] * p[iH2]**3 / (p[iCH4] * p[iH2O])) - np.log(Kp_SMR(T)),
                         np.log(p[iCO2] * p[iH2] / (p[iCO] * p[iH2O])) - np.log(Kp_WGS(T))])

    best = None
    for x0 in ([0.04, 0.02], [0.02, 0.01], [0.049, 0.03]):
        s = least_squares(res, x0, bounds=([1e-9, -0.05], [n0[iCH4] - 1e-9, 0.2]), xtol=1e-14, ftol=1e-14)
        if best is None or s.cost < best.cost:
            best = s
    n = comp(best.x)
    return n / n.sum()


# ----------------------------------------------------------------------------
# Post-processing: profiles, figures, validation, summary
# ----------------------------------------------------------------------------
# colours as in the paper: HMMC red, 57-4Q black, 25-4Q green
COLORS = {"HMMC": "tab:red", "57-4Q": "k", "25-4Q": "tab:green"}
PLOT_SPECIES = [("CH4", iCH4), ("H2", iH2), ("CO", iCO), ("CO2", iCO2), ("H2O", iH2O), ("He", iHe)]
PR2_VALIDATION_MM = 3.0     # paper Fig. 4 integrates/reports only the first 3 mm of P-R2


def profiles(res, pr2_mm=PR2_VALIDATION_MM, with_pr1=False):
    # stitches PBR (0-12 mm) and the first pr2_mm of P-R2 onto one axis; P-R1 is skipped unless with_pr1=True
    """Concatenate PBR [-> P-R1] -> first pr2_mm of P-R2. Returns z[mm], y[species, n], T, boundaries[mm]."""
    zb, Yb = res["PBR"]
    zs, Ys = [zb], [Yb]
    off = zb[-1]
    bounds = {"PBR_end": off * 1e3}
    if with_pr1:
        z1, Y1 = res["PR1"]
        zs.append(off + z1)
        Ys.append(Y1)
        off += z1[-1]
        bounds["PR1_end"] = off * 1e3
    z2, Y2 = res["PR2"]
    m = z2 <= pr2_mm * 1e-3
    zs.append(off + z2[m])
    Ys.append(Y2[:, m])
    Y = np.hstack(Ys)
    return np.concatenate(zs) * 1e3, mole_fractions(Y), Y[6], bounds


def _axlabel(with_pr1, pr2_mm):
    """x-axis label text"""
    return ("PBR | P-R1 | " if with_pr1 else "PBR | ") + f"first {pr2_mm:g} mm of P-R2"


def _mpl():
    """import matplotlib lazily with a non-interactive backend (works on servers)"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def plot_results(results, fname, title="SMR laboratory reactor", pr2_mm=PR2_VALIDATION_MM, with_pr1=False):
    """Species (CH4,H2,CO,CO2,H2O) + temperature, same layout as paper Fig. 4 (simulation lines only;
    validation.py adds the paper curves and experimental markers)."""
    plt = _mpl()
    fig, axs = plt.subplots(3, 2, figsize=(10, 10))
    axs = axs.ravel()
    for cat, res in results.items():                   # one line per catalyst in every panel
        z, y, T, b = profiles(res, pr2_mm, with_pr1)
        for ax, (nm, idx) in zip(axs[:5], PLOT_SPECIES[:5]):
            ax.plot(z, y[idx], color=COLORS.get(cat), label=cat)
            ax.set_title(nm)
        axs[5].plot(z, T, color=COLORS.get(cat), label=cat)
    axs[5].set_title("Temperature [K]")
    for ax in axs:
        for xb in b.values():
            ax.axvline(xb, ls="--", c="grey", lw=0.8)   # dashed line = PBR / P-R2 boundary
        ymax = ax.get_ylim()[1]
        ax.text(b["PBR_end"] / 2, ymax, "PBR", ha="center", va="top", fontsize=9)
        ax.text(b["PBR_end"] + (axs[0].get_xlim()[1] - b["PBR_end"]) / 2, ymax, "P-R2", ha="center", va="top", fontsize=9)
        ax.set_xlabel(f"axial position [mm]  ({_axlabel(with_pr1, pr2_mm)})")
        ax.grid(alpha=0.3)
    axs[0].legend()
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(fname, dpi=150)
    plt.close(fig)


def plot_axial_heatmaps(res, fname, catalyst, pr2_mm=PR2_VALIDATION_MM, with_pr1=False):
    """1-D axial field strips. NOT radial contours: the reference model is 1-D (axial only)."""
    plt = _mpl()
    z0, y, T, b = profiles(res, pr2_mm, with_pr1)
    # imshow needs a uniform grid, so resample every profile onto 3000 equidistant points
    z = np.linspace(z0[0], z0[-1], 3000)          # imshow needs a uniform grid; sections are sampled non-uniformly
    y = np.array([np.interp(z, z0, row) for row in y])
    T = np.interp(z, z0, T)
    # one colour strip per species plus temperature
    fields = [(nm, y[i], "mole fraction") for nm, i in PLOT_SPECIES] + [("Temperature", T, "K")]
    fig, axs = plt.subplots(len(fields), 1, figsize=(10, 1.5 * len(fields)), sharex=True)
    for ax, (nm, f, unit) in zip(axs, fields):
        im = ax.imshow(f[None, :], aspect="auto", extent=[z[0], z[-1], 0, 1], origin="lower")
        for xb in b.values():
            ax.axvline(xb, ls="--", c="white", lw=0.8)
        ax.set_yticks([])
        ax.set_ylabel(nm, rotation=0, ha="right", va="center")
        fig.colorbar(im, ax=ax, label=unit, pad=0.01)
    axs[-1].set_xlabel(f"axial position [mm]  ({_axlabel(with_pr1, pr2_mm)})")
    fig.suptitle(f"{catalyst}: 1-D axial fields (not radial contours; model is 1-D)")
    fig.tight_layout()
    fig.savefig(fname, dpi=150)
    plt.close(fig)


def point_at_pr2(res, mm):
    # state of the gas mm millimetres into P-R2 (used for the 3 mm validation point)
    """(mole fractions, T, molar flows) at mm into P-R2"""
    z2, Y2 = res["PR2"]
    k = min(int(np.searchsorted(z2, mm * 1e-3)), len(z2) - 1)
    return mole_fractions(Y2)[:, k], Y2[6, k], Y2[:6, k]


def write_summary(results, fname, pr2_mm=PR2_VALIDATION_MM):
    """CSV with compositions, T and X_CH4 at PBR exit, P-R1 exit, P-R2 at 3 mm and P-R2 outlet"""
    import csv
    with open(fname, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["catalyst", "section"] + SPECIES + ["T_K", "X_CH4"])
        for cat, res in results.items():
            F0 = res["F0"]
            rows = []
            for sec, lab in (("PBR", "PBR exit"), ("PR1", "P-R1 exit")):
                Y = res[sec][1]
                rows.append((lab, mole_fractions(Y)[:, -1], Y[6, -1], Y[:6, -1]))
            y, T, F = point_at_pr2(res, pr2_mm)
            rows.append((f"P-R2 at {pr2_mm:g} mm (validation point)", y, T, F))
            Y = res["PR2"][1]
            rows.append(("P-R2 full 40 mm (extrapolation beyond paper Fig. 4)",
                         mole_fractions(Y)[:, -1], Y[6, -1], Y[:6, -1]))
            for lab, y, T, F in rows:
                w.writerow([cat, lab] + [f"{v:.6f}" for v in y] + [f"{T:.2f}", f"{1 - F[iCH4] / F0[iCH4]:.4f}"])


def seed_sensitivity(catalyst="HMMC", **kw):
    """Re-run the chain for several H2 seeds; max relative change of PBR-exit and 3 mm compositions vs smallest seed."""
    seeds = [1e-6, 1e-8, 1e-10, 1e-12]
    pts = []
    for e in seeds:
        r = run(catalyst, h2_seed=e, n=200, **kw)
        pts.append(np.r_[mole_fractions(r["PBR"][1])[:5, -1], point_at_pr2(r, PR2_VALIDATION_MM)[0][:5]])
    print(f"\nH2-seed sensitivity ({catalyst}). The seed is a numerical regularisation only: "
          "the paper's feed has 0% H2 and the seed is not a physical component.")
    for e, p in zip(seeds, pts):
        print(f"  eps={e:.0e}  max rel. change vs eps=1e-12: {np.max(abs(p - pts[-1]) / abs(pts[-1])):.2e}")


# ----------------------------------------------------------------------------
# Main: command-line driver
# ----------------------------------------------------------------------------
def main():
    import os
    # ---- command-line options -------------------------------------------------------------
    ap = argparse.ArgumentParser(description="Costamagna 2020 SMR laboratory reactor model")
    ap.add_argument("--catalyst", default=None, choices=["HMMC", "57-4Q", "25-4Q"], help="default: all three")
    ap.add_argument("--T", type=float, default=853.0, help="furnace temperature [K]")
    ap.add_argument("--flow", type=float, default=120.0, help="total flow at RTP [mL/min]")
    ap.add_argument("--pr2-mm", type=float, default=PR2_VALIDATION_MM, help="P-R2 length shown/validated [mm]")
    ap.add_argument("--with-pr1", action="store_true", help="include the 50 mm P-R1 in plots (paper Fig. 4 omits it)")
    ap.add_argument("--h-length", choices=["tube", "particle"], default=None,
                    help="length scale in PBR Nu->h. Paper text says particle; default 'tube' reproduces Fig. 4 better (see compare_h_length.py)")
    ap.add_argument("--transport-config", choices=list(TRANSPORT_CONFIGS), default=None,
                    help="named PBR transport interpretation; default is v3_reproduction")
    ap.add_argument("--seed-study", action="store_true", help="run the H2-seed sensitivity study")
    ap.add_argument("--validate", action="store_true",
                    help="also run validation.py (Fig. 4 overlay + error table; needs the data/ folder)")
    ap.add_argument("--out", default=".", help="output directory")
    a = ap.parse_args()
    if a.transport_config is not None and a.h_length is not None:
        if TRANSPORT_CONFIGS[a.transport_config][1] != a.h_length:
            ap.error("--h-length conflicts with --transport-config Nu-to-h length")
    # The old validator always reruns defaults; prevent new alternatives from
    # being silently presented as validated. General custom-case repair is Phase 3.
    if a.validate and a.transport_config not in (None, "v3_reproduction"):
        ap.error("alternative transport validation: use tools/compare_transport.py")
    os.makedirs(a.out, exist_ok=True)

    def P(f):                                          # helper: path inside the output directory
        return os.path.join(a.out, f)

    # ---- simulate -------------------------------------------------------------------------
    reactor = Reactor(h_length=a.h_length or "tube", transport_config=a.transport_config)
    cats = [a.catalyst] if a.catalyst else ["HMMC", "57-4Q", "25-4Q"]
    results = {c: run(c, Tf=a.T, flow=a.flow, reactor=reactor) for c in cats}

    # ---- console report: compositions at PBR exit, P-R1 exit and P-R2 validation point ----
    print(f"Feed: CH4 5%, H2O 20%, He 75% | furnace {a.T:.0f} K | {a.flow:.0f} mL/min (RTP) | "
          f"PBR transport: {effective_transport(reactor)['configuration']} | H2 seed 1e-6 (numerical only)")
    for c, res in results.items():
        s = summarize(res)
        print(f"\n{c}")
        y3, T3, _ = point_at_pr2(res, a.pr2_mm)
        for sec, y, T in (("PBR exit", s["PBR"]["y"], s["PBR"]["T"]), ("P-R1 exit", s["PR1"]["y"], s["PR1"]["T"]),
                          (f"P-R2 @{a.pr2_mm:g} mm", y3, T3)):
            print(f"  {sec:14s} CH4={y[iCH4]:.4f} H2O={y[iH2O]:.4f} CO={y[iCO]:.4f} CO2={y[iCO2]:.4f} "
                  f"H2={y[iH2]:.4f} He={y[iHe]:.4f} T={T:.1f} K")
        print(f"  X_CH4={s['X_CH4']:.3f}  max|atom err|={max(abs(v) for v in s['atom_err'].values()):.1e}  "
              f"Tmin(PBR)={res['PBR'][1][6].min():.1f} K")

    # ---- figures and summary table --------------------------------------------------------
    plot_results(results, P("01_smr_results.png"), f"SMR, furnace {a.T:.0f} K, {a.flow:.0f} mL/min",
                 a.pr2_mm, a.with_pr1)
    for c, res in results.items():
        plot_axial_heatmaps(res, P(f"02_axial_fields_{c}.png"), c, a.pr2_mm, a.with_pr1)
    write_summary(results, P("smr_summary.csv"), a.pr2_mm)

    # ---- optional studies -----------------------------------------------------------------
    if a.seed_study:
        seed_sensitivity(cats[0], Tf=a.T, flow=a.flow, reactor=reactor)
    if a.validate:                                     # validation is a separate, mandatory step for the assignment
        import validation
        validation.main(out=a.out)
    else:
        print("\nNote: quantitative validation vs the paper is done by validation.py "
              "(or run this script with --validate).")
    print(f"\nFigures/CSV written to {os.path.abspath(a.out)}")


if __name__ == "__main__":
    main()
