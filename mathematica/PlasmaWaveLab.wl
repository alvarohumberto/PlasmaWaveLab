(* ::Package:: *)

(* ::Title:: *)
(*PlasmaWaveLab (Wolfram Language edition)*)


(* ::Text:: *)
(*Wolfram Language companion to the PlasmaWaveLab Python package (v0.3.0).*)
(*It implements the same models and numerical methods:*)
(*  - Sagdeev pseudopotential of the cold-ion / Boltzmann-electron plasma,*)
(*    critical Mach number, soliton amplitude and profile;*)
(*  - the small-amplitude (Washimi-Taniuti) KdV limit;*)
(*  - a periodic Fourier pseudospectral KdV solver with 2/3 de-aliasing and*)
(*    integrating-factor fourth-order Runge-Kutta (IF-RK4) time stepping;*)
(*  - spectral mass / L2 / Hamiltonian diagnostics.*)
(*All quantities are dimensionless (see docs/THEORY.md).*)
(**)
(*Authors: Alvaro H. Salas, Weaam Alhejaili, Samir A. El-Tantawy. License: MIT.*)
(*Usage:  Get["PlasmaWaveLab.wl"]   (or Needs after placing it on $Path)*)


BeginPackage["PlasmaWaveLab`"];


SagdeevPotential::usage =
"SagdeevPotential[phi, M] gives S(phi;M) = 1 - Exp[phi] + M^2 (1 - Sqrt[1 - 2 phi/M^2]) \
for the cold-ion/Boltzmann-electron plasma. Numeric arguments are checked (M > 1, phi < M^2/2); \
symbolic arguments return the expression.";

SagdeevPrime::usage =
"SagdeevPrime[phi, M] gives dS/dphi = (1 - 2 phi/M^2)^(-1/2) - Exp[phi].";

CriticalMach::usage =
"CriticalMach[] gives the upper Mach limit M_c of the soliton branch, the root of \
1 - Exp[M^2/2] + M^2 = 0 (M_c ~ 1.5852010652).";

PeakPotential::usage =
"PeakPotential[M] gives the soliton amplitude phi_m, the nonzero root of S(phi;M) = 0 below M^2/2. \
Returns $Failed (with a message) when no soliton exists (M <= 1 or M >= M_c). \
Options: \"Samples\" (2500) bracketing points, \"AccuracyGoal\" (13).";

SolitaryProfile::usage =
"SolitaryProfile[M] integrates phi'' = -S'(phi;M) from the crest (phi_m, 0) and returns an Association \
with keys \"Xi\", \"Phi\", \"Amplitude\", \"XiEnd\" and \"Function\". The profile is symmetric and is \
set to zero beyond the point where phi falls to Tail*phi_m. \
Options: \"XiMax\" (45), \"Points\" (1201), \"Tail\" (10^-5), \"PrecisionGoal\" (12), \"AccuracyGoal\" (14).";

SmallAmplitudeSoliton::usage =
"SmallAmplitudeSoliton[xi, M] gives the leading-order KdV approximation 3 delta Sech[Sqrt[delta/2] xi]^2, \
delta = M - 1, of the Sagdeev soliton.";

IonAcousticKdVCoefficients::usage =
"IonAcousticKdVCoefficients[] gives <|\"Alpha\" -> 1, \"Beta\" -> 1/2|>, the coefficients of the \
ion-acoustic KdV equation u_t + Alpha u u_x + Beta u_xxx = 0 in the present normalisation.";

KdVSoliton::usage =
"KdVSoliton[x, t] gives the exact soliton A Sech[Sqrt[alpha A/(12 beta)] (x - x0 - alpha A t/3)]^2 of \
u_t + alpha u u_x + beta u_xxx = 0. Options: \"Amplitude\" (1), \"Alpha\" (1), \"Beta\" (1), \"Center\" (0).";

KdVPropagate::usage =
"KdVPropagate[u0] evolves the periodic grid function u0 (even length N >= 16) with Fourier \
differentiation, 2/3 de-aliasing of u^2 and IF-RK4. Returns {times, states}. \
Options: \"Length\" (120), \"FinalTime\" (2), \"TimeStep\" (0.002), \"Alpha\" (1), \"Beta\" (1), \"SaveEvery\" (100).";

KdVInvariants::usage =
"KdVInvariants[u, dx] gives <|\"Mass\", \"L2\", \"Hamiltonian\"|> = dx*Total of u, u^2 and \
beta u_x^2/2 - alpha u^3/6, with u_x computed spectrally. Options: \"Alpha\" (1), \"Beta\" (1).";

KdVTwoSolitonPhaseShifts::usage =
"KdVTwoSolitonPhaseShifts[A1, A2] gives the exact phase shifts {d1, d2} of the faster (A1) and slower (A2) \
soliton after an overtaking collision. Options: \"Alpha\" (1), \"Beta\" (1).";

PlasmaWaveLabVersion::usage = "PlasmaWaveLabVersion gives the version string.";


Begin["`Private`"];


PlasmaWaveLabVersion = "0.3.0";

(* numpy-compatible FFT: forward without normalisation, exp(-i ...); inverse with 1/N *)
fft[u_] := Fourier[u, FourierParameters -> {1, -1}];
ifft[u_] := InverseFourier[u, FourierParameters -> {1, -1}];

(* integer wavenumber indices in numpy.fft.fftfreq order, N even *)
modeIndices[n_Integer] := Join[Range[0, n/2 - 1], Range[-n/2, -1]];


(* ::Section:: *)
(*Sagdeev model*)


SagdeevPotential::mach = "Supersonic Mach number M > 1 required (M = `1`).";
SagdeevPotential::dom = "Ion-density singularity: phi = `1` must be < M^2/2 = `2`.";

SagdeevPotential[phi_?NumericQ, mach_?NumericQ] := Which[
  mach <= 1, Message[SagdeevPotential::mach, mach]; $Failed,
  phi >= mach^2/2, Message[SagdeevPotential::dom, phi, mach^2/2]; $Failed,
  True, 1 - Exp[phi] + mach^2 (1 - Sqrt[1 - 2 phi/mach^2])
];
SagdeevPotential[phi_, mach_] := 1 - Exp[phi] + mach^2 (1 - Sqrt[1 - 2 phi/mach^2]);

SagdeevPrime[phi_, mach_] := (1 - 2 phi/mach^2)^(-1/2) - Exp[phi];


CriticalMach[] := CriticalMach[] = N[
  m /. FindRoot[1 - Exp[m^2/2] + m^2 == 0, {m, 3/2, 11/10, 2},
    WorkingPrecision -> 30, AccuracyGoal -> 25, PrecisionGoal -> 25]
];


PeakPotential::mach = "M must exceed 1 (M = `1`).";
PeakPotential::nowell = "No resolved positive-polarity Sagdeev well for M = `1`.";
PeakPotential::noroot = "No positive turning point before the ion-density limit (M = `1` >= critical Mach `2`).";

Options[PeakPotential] = {"Samples" -> 2500, "AccuracyGoal" -> 13};

PeakPotential[mach_?NumericQ, OptionsPattern[]] := Module[
  {top, grid, vals, neg, last, pos, upper, p},
  If[mach <= 1, Message[PeakPotential::mach, mach]; Return[$Failed]];
  top = N[mach^2/2 (1 - 10^-11)];
  grid = N @ Subdivide[Max[10^-8, top 10^-5], top, OptionValue["Samples"] - 1];
  vals = 1 - Exp[grid] + mach^2 (1 - Sqrt[1 - 2 grid/mach^2]);
  (* ignore round-off near phi = 0: last clearly negative sample followed by a positive one *)
  neg = Flatten @ Position[vals, _?(# < -10^-12 &), {1}, Heads -> False];
  If[neg === {}, Message[PeakPotential::nowell, mach]; Return[$Failed]];
  last = Last[neg];
  pos = Flatten @ Position[vals[[last + 1 ;;]], _?(# > 0 &), {1}, Heads -> False];
  If[pos === {}, Message[PeakPotential::noroot, mach, CriticalMach[]]; Return[$Failed]];
  upper = last + First[pos];
  p /. FindRoot[1 - Exp[p] + mach^2 (1 - Sqrt[1 - 2 p/mach^2]),
    {p, grid[[last]], grid[[upper]]}, Method -> "Brent",
    AccuracyGoal -> OptionValue["AccuracyGoal"], PrecisionGoal -> OptionValue["AccuracyGoal"]]
];


Options[SolitaryProfile] = {"XiMax" -> 45, "Points" -> 1201, "Tail" -> 10^-5,
  "PrecisionGoal" -> 12, "AccuracyGoal" -> 14};

SolitaryProfile[mach_?NumericQ, OptionsPattern[]] := Module[
  {amp, ximax, tail, sol, phiFn, end, xs, prof, phi, x},
  amp = PeakPotential[mach];
  If[FailureQ[amp], Return[$Failed]];
  ximax = N @ OptionValue["XiMax"]; tail = N @ OptionValue["Tail"];
  (* the far field phi = 0 is a saddle: stop near it instead of integrating through it *)
  sol = NDSolve[{
      phi''[x] == -((1 - 2 phi[x]/mach^2)^(-1/2) - Exp[phi[x]]),
      phi[0] == amp, phi'[0] == 0,
      WhenEvent[phi[x] < tail amp, "StopIntegration"],
      (* safeguard: stop if round-off pushes the orbit off the separatrix and it turns back *)
      WhenEvent[phi'[x] > 0, "StopIntegration"]},
    phi, {x, 0, ximax},
    Method -> {"ExplicitRungeKutta", "DifferenceOrder" -> 5},
    PrecisionGoal -> OptionValue["PrecisionGoal"],
    AccuracyGoal -> OptionValue["AccuracyGoal"], MaxStepSize -> 0.05];
  phiFn = phi /. First[sol];
  end = phiFn["Domain"][[1, 2]];
  xs = N @ Subdivide[-ximax, ximax, OptionValue["Points"] - 1];
  prof = If[Abs[#] <= end, phiFn[Abs[#]], 0.] & /@ xs;
  <|"Xi" -> xs, "Phi" -> prof, "Amplitude" -> amp, "XiEnd" -> end, "Function" -> phiFn|>
];


IonAcousticKdVCoefficients[] := <|"Alpha" -> 1, "Beta" -> 1/2|>;

SmallAmplitudeSoliton[xi_, mach_?NumericQ] /; mach > 1 :=
  3 (mach - 1) Sech[Sqrt[(mach - 1)/2] xi]^2;


(* ::Section:: *)
(*KdV model*)


Options[KdVSoliton] = {"Amplitude" -> 1, "Alpha" -> 1, "Beta" -> 1, "Center" -> 0};

KdVSoliton[x_, t_, OptionsPattern[]] := With[
  {A = OptionValue["Amplitude"], a = OptionValue["Alpha"], b = OptionValue["Beta"],
   c = OptionValue["Center"]},
  A Sech[Sqrt[a A/(12 b)] (x - c - a A t/3)]^2
];


KdVPropagate::grid = "Require an even number of grid points N >= 16 and positive Length, TimeStep and SaveEvery.";

Options[KdVPropagate] = {"Length" -> 120, "FinalTime" -> 2, "TimeStep" -> 0.002,
  "Alpha" -> 1, "Beta" -> 1, "SaveEvery" -> 100};

KdVPropagate[u0_List, OptionsPattern[]] := Module[
  {u, n, len, tf, dt, alpha, beta, save, k, lin, mask, nonlin, steps, h, eFull, eHalf,
   U, a, b, c, d, out},
  u = N[u0]; n = Length[u];
  {len, tf, dt, alpha, beta, save} = N @ OptionValue[
    {"Length", "FinalTime", "TimeStep", "Alpha", "Beta", "SaveEvery"}];
  If[n < 16 || OddQ[n] || len <= 0 || dt <= 0 || tf < 0 || save < 1,
    Message[KdVPropagate::grid]; Return[$Failed]];
  save = Round[save];
  k = N[2 Pi/len modeIndices[n]];
  lin = I beta k^3;
  mask = N @ UnitStep[Quotient[n, 3] - Abs[modeIndices[n]]];   (* keep |m| <= N/3 *)
  nonlin[W_] := -0.5 I alpha k (fft[Re[ifft[W]]^2] mask);
  steps = Ceiling[tf/dt];
  If[steps == 0, Return[{{0.}, {u}}]];
  h = tf/steps;
  eFull = Exp[lin h]; eHalf = Exp[lin h/2];
  U = fft[u];
  out = Reap[
    Sow[{0., u}];
    Do[
      a = nonlin[U];
      b = nonlin[eHalf (U + h a/2)]/eHalf;
      c = nonlin[eHalf (U + h b/2)]/eHalf;
      d = nonlin[eFull (U + h c)]/eFull;
      U = eFull (U + h (a + 2 b + 2 c + d)/6);
      If[Mod[s, save] == 0 || s == steps, Sow[{s h, Re[ifft[U]]}]],
      {s, 1, steps}]
  ][[2, 1]];
  {out[[All, 1]], out[[All, 2]]}
];


Options[KdVInvariants] = {"Alpha" -> 1, "Beta" -> 1};

KdVInvariants[u_List, dx_?NumericQ, OptionsPattern[]] := Module[
  {n = Length[u], k, ux, alpha = N @ OptionValue["Alpha"], beta = N @ OptionValue["Beta"]},
  k = N[2 Pi/(n dx) modeIndices[n]];
  ux = Re[ifft[I k fft[N[u]]]];
  <|"Mass" -> dx Total[u], "L2" -> dx Total[u^2],
    "Hamiltonian" -> dx Total[beta ux^2/2 - alpha u^3/6]|>
];


Options[KdVTwoSolitonPhaseShifts] = {"Alpha" -> 1, "Beta" -> 1};

KdVTwoSolitonPhaseShifts[a1_?NumericQ, a2_?NumericQ, OptionsPattern[]] := Module[
  {k1, k2, lg, al = OptionValue["Alpha"], be = OptionValue["Beta"]},
  k1 = Sqrt[al a1/(12 be)]; k2 = Sqrt[al a2/(12 be)];
  lg = Log[(k1 + k2)/(k1 - k2)];
  N[{lg/k1, -lg/k2}]
];


End[];

EndPackage[];
