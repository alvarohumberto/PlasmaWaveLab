#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Title:: *)
(*PlasmaWaveLab - Wolfram Language verification*)


(* ::Text:: *)
(*Reproduces the verification results of the SoftwareX manuscript with the Wolfram Language*)
(*implementation and compares them with the reference values produced by the Python package*)
(*(figures/verification/summary.json). Run either as*)
(*    wolframscript -file Verification.wls*)
(*or open it in Mathematica and evaluate all cells. Figures are written to mathematica/output/.*)


$dir = Which[
  StringQ[$InputFileName] && $InputFileName =!= "", DirectoryName[$InputFileName],
  $FrontEnd =!= Null, Quiet @ Check[NotebookDirectory[], Directory[]],
  True, Directory[]];
Get[FileNameJoin[{$dir, "PlasmaWaveLab.wl"}]];
$out = FileNameJoin[{$dir, "output"}];
If[!DirectoryQ[$out], CreateDirectory[$out]];

$results = {};
check[name_String, value_, ok_] := (
  AppendTo[$results, {name, value, If[TrueQ[ok], "PASS", "FAIL"]}];
  Print[StringPadRight[name, 46],
    With[{s = ToString[value, InputForm]}, StringPadRight[s, Max[34, StringLength[s] + 1]]],
    If[TrueQ[ok], "PASS", "FAIL"]]);

Print["PlasmaWaveLab (Wolfram Language) ", PlasmaWaveLabVersion, " on ", $Version];
t0 = AbsoluteTime[];

(* Python reference values (v0.3.0, figures/verification/summary.json) *)
ref = <|
  "Mc" -> 1.5852010652445132,
  "phi12" -> 0.5243854849359129,
  "smallDelta" -> {0.005, 0.01, 0.02, 0.04, 0.08, 0.16},
  "smallProfErr" -> {0.006421356981016364, 0.012781337863045777, 0.02532042413282467,
                     0.0496914914034968, 0.09574515819613129, 0.17804954986314603},
  "spatialN" -> {48, 64, 96, 128, 160, 192, 256, 320, 384},
  "spatialErr" -> {0.07487601010581213, 0.018659606312392274, 0.0007451000011412523,
                   5.308104438039862*^-5, 3.461502785628845*^-6, 1.5943156128628289*^-7,
                   7.559088910369383*^-10, 1.7502665983215593*^-12, 2.195647863291481*^-14},
  "dts" -> {0.08, 0.04, 0.02, 0.01, 0.005},
  "temporalErr" -> {3.750341497732457*^-5, 1.0659862925699315*^-6, 7.118111722270726*^-8,
                    4.530812613978696*^-9, 2.849900315737841*^-10},
  "phaseNumerical" -> {2.6982747395833293, -5.381705729166669}
|>;

relClose[a_, b_, tol_] := Abs[a - b] <= tol Max[Abs[b], 10^-300];


(* ::Section:: *)
(*1. Sagdeev model: critical Mach number, amplitude, existence*)


mc = CriticalMach[];
check["Critical Mach number M_c", mc, Abs[mc - ref["Mc"]] < 10^-9];

phi12 = PeakPotential[1.2];
check["Amplitude phi_m(M = 1.2)", phi12, Abs[phi12 - ref["phi12"]] < 10^-9];

check["No soliton above M_c (M = 1.59)", "$Failed",
  FailureQ[Quiet @ PeakPotential[1.59]]];

machs = Subdivide[1.005, mc - 10^-4, 119];
amps = PeakPotential /@ machs;
check["phi_m increasing on 1 < M < M_c", "monotone", And @@ Positive[Differences[amps]]];

sagdeevPlot = GraphicsRow[{
   Quiet @ Plot[Evaluate[SagdeevPotential[p, #] & /@ {1.1, 1.2, 1.4, 1.58}],
     {p, 0, 1.3}, PlotRange -> {-0.4, 0.06}, AxesLabel -> {"\[Phi]", "S(\[Phi];M)"},
     PlotLegends -> {"M = 1.1", "M = 1.2", "M = 1.4", "M = 1.58"},
     PlotStyle -> {Automatic, Dashed, DotDashed, Dotted},
     PlotLabel -> "(a) Sagdeev pseudopotential",
     RegionFunction -> Function[{p, s}, s < 0.06]],
   ListLinePlot[
     Table[With[{r = SolitaryProfile[m, "XiMax" -> 25, "Points" -> 1001]},
       Transpose[{r["Xi"], r["Phi"]}]], {m, {1.1, 1.2, 1.4, 1.58}}],
     PlotRange -> {{-15, 15}, All}, AxesLabel -> {"\[Xi]", "\[Phi]"},
     PlotStyle -> {Automatic, Dashed, DotDashed, Dotted},
     PlotLabel -> "(b) solitary-wave profiles"]}, ImageSize -> 900];
Export[FileNameJoin[{$out, "fig_sagdeev.pdf"}], sagdeevPlot];


(* ::Section:: *)
(*2. Small-amplitude (KdV) limit*)


smallErr = Table[
   With[{m = 1 + d, xm = Max[60, 12/Sqrt[d/2]]},
     With[{r = SolitaryProfile[m, "XiMax" -> xm, "Points" -> 4001]},
       Max[Abs[r["Phi"] - SmallAmplitudeSoliton[r["Xi"], m]]]/(3 d)]],
   {d, ref["smallDelta"]}];
slope = Last @ LeastSquares[
    Transpose[{ConstantArray[1., Length[smallErr]], Log[ref["smallDelta"]]}], Log[smallErr]];
check["Small-amplitude profile errors vs Python", smallErr,
  And @@ MapThread[relClose[#1, #2, 10^-3] &, {smallErr, ref["smallProfErr"]}]];
check["Small-amplitude error order (~0.96)", slope, Abs[slope - 0.96] < 0.02];

existencePlot = GraphicsRow[{
   ListLinePlot[{Transpose[{machs, amps}], Transpose[{machs, 3 (machs - 1)}],
      Transpose[{machs, machs^2/2}]}, AxesLabel -> {"M", "amplitude"},
     PlotLegends -> {"\[Phi]_m", "KdV limit 3(M-1)", "ion limit M^2/2"},
     PlotStyle -> {Automatic, Dashed, Dotted}, PlotLabel -> "(a) existence domain"],
   ListLogLogPlot[Transpose[{ref["smallDelta"], smallErr}], Joined -> True,
     PlotMarkers -> Automatic, AxesLabel -> {"\[Delta] = M-1", "relative error"},
     PlotLabel -> "(b) small-amplitude limit"]}, ImageSize -> 900];
Export[FileNameJoin[{$out, "fig_existence.pdf"}], existencePlot];


(* ::Section:: *)
(*3. KdV solver: spatial and temporal convergence*)


L = 120.;
spatialErr = Table[
   With[{x = Most @ Subdivide[-L/2, L/2, n]},
     Max @ Abs[Last[KdVPropagate[KdVSoliton[x, 0.], "Length" -> L, "FinalTime" -> 2.,
         "TimeStep" -> 0.001][[2]]] - KdVSoliton[x, 2.]]],
   {n, ref["spatialN"]}];
check["Spatial error at N = 256 (~7.6e-10)", spatialErr[[7]],
  relClose[spatialErr[[7]], ref["spatialErr"][[7]], 0.05]];
check["Spatial error at N = 384 (round-off)", spatialErr[[9]], spatialErr[[9]] < 10^-12];

xT = Most @ Subdivide[-L/2, L/2, 512];
temporalErr = Table[
   Max @ Abs[Last[KdVPropagate[KdVSoliton[xT, 0., "Amplitude" -> 2.], "Length" -> L,
        "FinalTime" -> 4., "TimeStep" -> dt][[2]]] - KdVSoliton[xT, 4., "Amplitude" -> 2.]],
   {dt, ref["dts"]}];
orders = Log2[Most[temporalErr]/Rest[temporalErr]];
check["Temporal errors vs Python", temporalErr,
  And @@ MapThread[relClose[#1, #2, 0.02] &, {temporalErr, ref["temporalErr"]}]];
check["Observed temporal order (finest pair)", Last[orders], Abs[Last[orders] - 4] < 0.1];

convPlot = GraphicsRow[{
   ListLogPlot[Transpose[{ref["spatialN"], spatialErr}], Joined -> True,
     PlotMarkers -> Automatic, AxesLabel -> {"N", "max error"},
     PlotLabel -> "(a) spectral convergence"],
   ListLogLogPlot[{Transpose[{ref["dts"], temporalErr}],
      Transpose[{ref["dts"], Last[temporalErr] (ref["dts"]/Last[ref["dts"]])^4}]},
     Joined -> True, PlotStyle -> {Automatic, Dashed}, PlotMarkers -> {Automatic, None},
     PlotLegends -> {"IF-RK4", "slope 4"}, AxesLabel -> {"\[CapitalDelta]t", "max error"},
     PlotLabel -> "(b) fourth order in time"]}, ImageSize -> 900];
Export[FileNameJoin[{$out, "fig_convergence.pdf"}], convPlot];


(* ::Section:: *)
(*4. Two-soliton collision, phase shifts and invariants*)


Lc = 150.; nc = 512; T = 100.;
xc = Most @ Subdivide[-Lc/2, Lc/2, nc];
u0 = KdVSoliton[xc, 0., "Amplitude" -> 2., "Center" -> -50.] +
     KdVSoliton[xc, 0., "Amplitude" -> 0.5, "Center" -> -25.];
tc = AbsoluteTime[];
{times, states} = KdVPropagate[u0, "Length" -> Lc, "FinalTime" -> T, "TimeStep" -> 0.005,
   "SaveEvery" -> 100];
collisionSeconds = AbsoluteTime[] - tc;

inv = KdVInvariants[#, Lc/nc] & /@ states;
drift[key_] := Max[Abs[(#[key] - First[inv][key])/First[inv][key]] & /@ inv];
check["Mass drift, 2x10^4 steps (round-off)", drift["Mass"], drift["Mass"] < 10^-14];
check["L2 drift, 2x10^4 steps (~5.6e-11)", drift["L2"], drift["L2"] < 2*^-10];
check["Hamiltonian drift, 2x10^4 steps (~1.0e-10)", drift["Hamiltonian"],
  drift["Hamiltonian"] < 3*^-10];

(* sub-grid crest position from the band-limited (trigonometric) interpolant *)
crest[u_, guess_] := Module[{n = Length[u], uh, k, f, y},
  uh = Fourier[u, FourierParameters -> {1, -1}];
  k = 2 Pi/Lc Join[Range[0, n/2 - 1], Range[-n/2, -1]];
  f[y_?NumericQ] := Re[Total[uh Exp[I k (y - First[xc])]]]/n;
  y /. Last[FindMaximum[f[y], {y, guess}]]];
uT = Last[states];
iBig = First @ Ordering[uT, -1];
xBig = crest[uT, xc[[iBig]]];
iSmall = First @ Ordering[uT (UnitStep[(xBig - 10) - xc]), -1];
xSmall = crest[uT, xc[[iSmall]]];
shiftsNum = {xBig - (-50 + 2 T/3), xSmall - (-25 + 0.5 T/3)};
shiftsExact = KdVTwoSolitonPhaseShifts[2., 0.5];
check["Phase shifts (numerical)", shiftsNum,
  Max[Abs[shiftsNum - ref["phaseNumerical"]]] < 10^-3];
check["Phase shifts vs exact {2.691, -5.382}", shiftsExact,
  Max[Abs[shiftsNum - shiftsExact]] < 0.01];

collisionPlot = GraphicsRow[{
   ListDensityPlot[states[[;; ;; 2, ;; ;; 2]], DataRange -> {{-Lc/2, Lc/2 - 2 Lc/nc}, {0, T}},
     ColorFunction -> "BlueGreenYellow", FrameLabel -> {"x", "t"},
     PlotLabel -> "(a) two-soliton collision", PlotLegends -> Automatic],
   ListLogPlot[{
      Transpose[{times, Max[#, 10^-17] & /@ FoldList[Max, Abs[(#["Mass"] - inv[[1]]["Mass"])/inv[[1]]["Mass"]] & /@ inv]}],
      Transpose[{times, Max[#, 10^-17] & /@ FoldList[Max, Abs[(#["L2"] - inv[[1]]["L2"])/inv[[1]]["L2"]] & /@ inv]}],
      Transpose[{times, Max[#, 10^-17] & /@ FoldList[Max, Abs[(#["Hamiltonian"] - inv[[1]]["Hamiltonian"])/inv[[1]]["Hamiltonian"]] & /@ inv]}]},
     Joined -> True, PlotRange -> {10^-17, 10^-6}, PlotStyle -> {Automatic, Dashed, Dotted},
     PlotLegends -> {"mass", "L2 norm", "Hamiltonian"}, AxesLabel -> {"t", "max relative drift"},
     PlotLabel -> "(b) invariant drift"]}, ImageSize -> 900];
Export[FileNameJoin[{$out, "fig_collision.pdf"}], collisionPlot];


(* ::Section:: *)
(*Summary*)


nFail = Count[$results[[All, 3]], "FAIL"];
Print["\n", Length[$results] - nFail, " passed, ", nFail, " failed; collision run ",
  NumberForm[collisionSeconds, {6, 1}], " s; total ", NumberForm[AbsoluteTime[] - t0, {6, 1}], " s"];
Export[FileNameJoin[{$out, "verification_results.csv"}],
  Prepend[MapAt[ToString[#, InputForm] &, $results, {All, 2}], {"check", "value", "status"}]];
Export[FileNameJoin[{$out, "summary.json"}], <|
   "version" -> PlasmaWaveLabVersion, "wolfram" -> $Version,
   "critical_mach" -> mc, "phi_m_1.2" -> phi12, "small_amplitude_error" -> smallErr,
   "spatial_error" -> spatialErr, "temporal_error" -> temporalErr, "temporal_orders" -> orders,
   "phase_shift_numerical" -> shiftsNum, "phase_shift_exact" -> shiftsExact,
   "max_rel_drift" -> <|"mass" -> drift["Mass"], "l2" -> drift["L2"],
      "hamiltonian" -> drift["Hamiltonian"]|>, "passed" -> Length[$results] - nFail,
   "failed" -> nFail|>, "JSON"];
If[!$Notebooks, Exit[If[nFail == 0, 0, 1]]];
