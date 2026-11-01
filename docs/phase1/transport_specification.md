# V3 transport-model scientific specification

Scope: existing `v3-baseline` calculations only. No alternatives are implemented
in Phase 1. Source references below use the unchanged model's line numbering.

## Sources and governing assumptions

Costamagna et al., Energies 2020, 13, 2624, DOI 10.3390/en13102624;
local reference: `../energies-13-02624.pdf` relative to repository root.
Relevant pages: 6 (balances and resistance model), 7 (Eqs. 11–16),
8 (dispersion assumptions), 9–10 (post-reactors), 11 (Table 5), 13–14 (Fig. 4).
The paper is outside this repository; its SHA-256 is recorded in the run report.

Steady, ideal-gas, pseudo-homogeneous 1-D plug flow; constant pressure 1.01325 bar;
solid/gas temperature equality; effectiveness factors treated as unity; no axial
mass/thermal dispersion, pressure-drop equation or radial field solution.
Furnace radiation is treated as ideal through the wall-temperature assumption.
Holders are omitted, matching the paper. These assumptions are preserved, not
independently demonstrated across the 700–1100 K operating sweep in this phase.

## Exact implementation

`hyd.py:143–216,284–332` defines the transport path. At every local state:

- `y_i=F_i/sum(F)`; mixture molar mass M is in kg/mol after conversion from g/mol.
- `rho=P M/(R T)` [kg/m3]. `us=sum(F) R T/(P A)` [m/s], A=pi d²/4 [m2].
- Pure viscosity: Chapman–Enskog with Lennard-Jones parameters and Neufeld
  collision integral (`mu_pure`); units Pa s. Mixture viscosity uses Wilke mixing.
- Pure conductivity: modified Eucken, `ki=mui/Mi*(1.32*(cpi-R)+1.77*R)` [W/m/K].
  `kf=sum(yi*ki)` is mole-weighted **gas** conductivity [W/m/K].
- `cp=sum(yi*cpi)/M` [J/kg/K]; species Shomate cp values are in J/mol/K.
- `ks(T)` [W/m/K] is a quadratic interpolant through (473,1.57), (773,1.99),
  (1173,3.39). Wall conductivity `kw_glass=1.6` is separately constant.
- `kappa=ks/kf`; void fraction eps=0.4; `B=1.25*((1-eps)/eps)^(10/9)`;
  `D=1-B/kappa`; `S=sqrt(1-eps)`.
  `kr=kf*[1-S+2*S/D*(B*(1-1/kappa)/D²*ln(kappa/B)
                     -(B-1)/D-(B-1)/2)]` [W/m/K].
  This matches the structure/sign of paper Eq. 14. Exponent 10/9 is the code's
  approximation to the paper's printed 1.11; it is unchanged.
- **Re**=`rho*us*d/mu`, with tube diameter d=0.006 m, dimensionless.
- **Pr**=`cp*mu/kf`, with gas conductivity, dimensionless.
- **Nu**=`2.67+0.53*Re^0.77*Pr^0.53`, dimensionless (paper Eq. 16).
- Default **h**=`Nu*kr/d` [W/m2/K]. Existing `--h-length particle` changes only
  this denominator to dp=0.0002 m. It does not change Re or Pr.
- **U**=`U_scale/(1/h+s/kw_glass)` [W/m2/K], s=0.001 m, default U_scale=1.
  The thin-wall resistance and `a=4/d` convention match paper Eqs. 9–10;
  no extra cylindrical-wall correction is introduced.

P-R1/P-R2 use Nu=3.66 and `h=3.66*kf/d`; U uses the same series resistance,
with glass kw=1.6 or steel kw=16.3 W/m/K respectively. P-R1 wall is at furnace
T; P-R2 wall is at 523 K. Surface WGS is converted to volumetric rate with 4/d.
`k_axial_votruba` is auxiliary only, uses natural logarithms and is not an axial
conduction term in the energy equation.

Species rates and the energy balance are integrated versus volume V=A*z.
Solver is LSODA, rtol=1e-8, atol=[1e-16 mol/s]*6 plus 1e-6 K; PBR/P-R1
have 400 saved samples, P-R2 2000. Output sampling is not adaptive-step control.
The H2 seed is 1e-6 before feed renormalization. Catalyst factors and equilibrium
model D remain identical to V3.

## Paper ambiguities and reconciliation limits

| Item | Printed equations/text | Table 5 | V3 interpretation |
|---|---|---|---|
| Re length | Eq. 11 L, text says particle | Row explicitly rho us d/mu | Tube |
| Pr conductivity | Eq. 12 kr (bed) | Row explicitly mu cp/kf (gas) | Gas |
| Nu-to-h | Eq. 13 kr/L, text says particle L | No explicit conversion-length row | kr/tube |
| ks | Page 7 text approximately 1.6 | 1.57,1.99,3.39 at listed temperatures | Temperature fit |
| kw | Page 7 glass approximately 1.6 | Not independently tabulated | Constant 1.6 |

Consequently, V3 is a stated reproduction interpretation, not a literal
implementation of every printed definition. The old h-only particle result is
not evidence that a complete particle-length convention has been evaluated.
Lowest Figure 4 RMSE is not a criterion for establishing physical correctness.

Table 5 lists SMR at 773 and 1173 K with GHSV 2.4e4/h; Figure 4 is at 853 K,
GHSV 2.1e4/h (120 mL/min RTP). They must not share an assumed inlet case.
Even Table 5's rounded rows are not wholly internally reproducible: at 773 K,
its rho=.0474, us=.200, d=.006 and mu=3.56e-5 give Re about 1.598, consistent
with 1.600; cp=4520 and kf=.242 give Pr about .665. Substituting those rounded
Re/Pr into Eq. 16 yields Nu about 3.28, rather than the tabulated 2.91.
This mismatch cannot be explained as a model-property discrepancy alone.
It does not identify the authors' actual implementation or establish an erratum.

Phase 2 should distinguish: (a) algebra using tabulated properties, (b) properties
calculated from stated feed/conditions, and (c) full Figure 4 simulations. Keep
uncertain Table 5 operating inputs explicit rather than fitting them to U.

The parameter-row swap in Table 4 is an existing equilibrium interpretation;
Phase 1 neither changes nor independently validates it. The thermodynamic
consistency check shares equilibrium constants with the kinetic model.

## Verification versus proposed improvements

Executed Phase 1 evidence belongs in `baseline_reproduction_report.md` and the
versioned run directory. No full solver convergence, independent solver,
energy-integral check, new transport formulation, or physical validation at
80 mL/min is claimed. These are subsequent approved phases.

Proposed Phase 2 configurations remain the four previously agreed names:
`v3_reproduction`, `v3_h_particle_only`, `particle_lengths_gas_pr`, and
`printed_equations`. Preserve the ks fit initially; consider constant ks as a
separate sensitivity only after explaining the source ambiguity. Do not retune
kinetics or change the default according to fitting error.
