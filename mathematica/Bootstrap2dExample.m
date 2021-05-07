(* ::Package:: *)

(*Write here the path of the folder containing sdpb and pvm2sdp*)
executablesPath = "../build";


<<"SDPB.m";

half = SetPrecision[1/2, prec];

rho[z_] := z/(1+Sqrt[1-z])^2;

(* An SL_2 conformal block in the rho coordinate *)
chiralBlock[x_, rho_] := rho^(x/2) Hypergeometric2F1[1/2,x/2,(x+1)/2,rho^2];

(* Turn a series a_0 + a_1 x + ... into a list of rules { name[0] ->
   0! a_0, name[1] -> 1! a_1, ... } *)
seriesToRules[name_, series_] := Table[
    (* series[[3]] is the list of coefficients *)
    name[n] -> n! series[[3, n+1]],
    {n, 0, Length[series[[3]]] - 1}];

(* A list of replacement rules zDeriv[n] -> sum of rhoDeriv[k]'s
giving the derivatives of a function with respect to z around 1/2 in
terms of derivatives with respect to rho around rho[1/2] *)
zDerivTable[order_] :=
    (zDerivTable[order] =
     Module[{dRhoSeries = Series[rho[half+dz] - rho[half], {dz, 0, order}]},
            seriesToRules[zDeriv, Sum[rhoDeriv[n] dRhoSeries^n / n!, {n, 0, order}]]]);

(* A table of the form { prefactor, { zDeriv[n] -> polynomial_n(x),
... } }, where prefactor * polynomial_n(x) approximates the n-th
derivative of a chiral conformal bock around z=1/2 *)
chiralBlockTable[derivativeOrder_, keptPoleOrder_] :=
    (chiralBlockTable[derivativeOrder, keptPoleOrder] =
     Module[
         {
             numerator,
             prefactor,
             poles,
             derivTable
         },
         poles = Table[n, {n, 1, keptPoleOrder-1, 2}];
         numerator = (Series[chiralBlock[x,rho], {rho, 0, keptPoleOrder}] Product[x+n, {n, poles}]) // Normal;
         prefactor = DampedRational[1,Table[-n,{n,poles}],rhoCrossing^(1/2),x];
         (* We compute derivatives with respect to rho and use
         zDerivTable to convert to zDerivatives *)
         derivTable = zDerivTable[derivativeOrder] /. Table[
             rhoDeriv[n] -> ((rho^(-x/2) D[numerator, {rho,n}] // Expand) /. rho->rhoCrossing),
             {n, 0, derivativeOrder}];
         {prefactor,Expand[derivTable]}]);

(* Replacement rules giving derivatives of (1-z)^deltaPhi f[z] around
z=1/2 in terms of derivatives of f[z] *)
withDeltaPhiDerivTable[deltaPhi_,order_] :=
    seriesToRules[
        withDeltaPhiDeriv,
        Series[(half-dz)^deltaPhi f[half+dz], {dz, 0, order}] /. {
            Derivative[j_][f][_] :> zDeriv[j],
            f[_] :> zDeriv[0]}];

(* Derivatives such that m >= n, m+n <= derivativeOrder, and m+n is odd *)
oddDerivs[derivativeOrder_] :=
    Flatten[Table[zzbDeriv[m,n]/(m! n!), 
                  {m,0,derivativeOrder},
                  {n, 1 - Mod[m,2], Min[m, derivativeOrder-m], 2}]
            ,1];

(* Test whether the point (deltaPhi, deltaPhiSqLowPrecision) is
allowed in a Z_2-symmetric CFT.  derivativeOrder gives the number of
derivatives to use, keptPoleOrder controls the accuracy of the
conformal block approximation, and Lmax sets the number of included
spins *)
singletAllowed2d[deltaPhiLowPrecision_, deltaPhiSqLowPrecision_, derivativeOrder_, keptPoleOrder_, Lmax_] :=
    Module[
        {
            deltaPhi  = SetPrecision[deltaPhiLowPrecision,prec],
            deltaPhiSq = SetPrecision[deltaPhiSqLowPrecision,prec],
            chiralBlocksPrefactor,
            chiralBlockPols,
            chiralBlocksWithDeltaPhi,
            unitOp,
            pols,
            norm,
            obj
        },
        {chiralBlocksPrefactor, chiralBlockPols} = chiralBlockTable[derivativeOrder, keptPoleOrder];
        chiralBlocksWithDeltaPhi = withDeltaPhiDerivTable[deltaPhi, derivativeOrder] /. chiralBlockPols;
        pols=Table[
            PositiveMatrixWithPrefactor[
                (chiralBlocksPrefactor/.x->x+2L) (chiralBlocksPrefactor),
                (* These are 1x1 matrices of polynomial vectors *)
                {{oddDerivs[derivativeOrder]/.{
                    zzbDeriv[m_,n_] :> (
                        (withDeltaPhiDeriv[m]/.chiralBlocksWithDeltaPhi /. x -> x+2L)
                        (withDeltaPhiDeriv[n]/.chiralBlocksWithDeltaPhi) + 
                        (withDeltaPhiDeriv[n]/.chiralBlocksWithDeltaPhi /. x -> x+2L)
                        (withDeltaPhiDeriv[m]/.chiralBlocksWithDeltaPhi))
                                              }//Expand}}],
            {L, 0, Lmax, 2}];
        (* Replace x -> x + deltaPhiSq for scalar operators *)
        pols = MapAt[# /. x -> x + deltaPhiSq &, pols, 1];
        unitOp = oddDerivs[derivativeOrder] /. {
            zzbDeriv[m_,n_] :> 2 withDeltaPhiDeriv[m] withDeltaPhiDeriv[n]
                                               } /. withDeltaPhiDerivTable[deltaPhi, derivativeOrder] /. {
                                                   zDeriv[0] :> 1, zDeriv[_] :> 0};
        norm = unitOp;
        obj = 0*norm;
        SolveBootstrapSDP[SDP[obj,norm,pols]]];

(* This is not a recommended long-term solution for evaluating SDPs
for the following reasons: 1) different sdpFiles should have unique
names so that they can be solved in parallel and so that their
checkpoints and output files don't overwrite each other; 2) The
Run[...] command forces Mathematica to be running until sdpb
finishes. It is better to use WriteBootstrapSDP instead of
SolveBootstrapSDP and run sdpb by hand or with an external script. *)
SolveBootstrapSDP[sdp_] := Module[
    {
        sdpFile = "mySDP.xml",
        outFolder = "mySDP_out",
        sdpFolder = "mySDP",
        fullResult
    },
    WriteBootstrapSDP[sdpFile, sdp];
    (* Most of the defaults are way over the top for this size
    problem, but we'll use them because it's easy. If you want speed,
    try fiddling with some of the parameters. *)
    WriteString["stdout","Converting xml to pvm format\n"];
    Run[StringRiffle[{
        FileNameJoin[{executablesPath,"pvm2sdp"}],
        ToString[prec], sdpFile, sdpFolder
    }," "]];
    WriteString["stdout","Running sdpb\n"];
    Run[StringRiffle[{
        "mpirun",
        FileNameJoin[{executablesPath,"sdpb"}],
        "-s", sdpFolder, "-o", outFolder,
        "--procsPerNode", "6",
        "--precision", ToString[prec],
        "--findPrimalFeasible",
        "--findDualFeasible",
        "--noFinalCheckpoint"
    }," "]];
    (* Careful! Simply evaluating the output file will bring a bunch
    of variables into scope! *)
    fullResult = Import[FileNameJoin[{outFolder,"out.txt"}], "String"];
    terminateReason = StringCases[fullResult,"terminateReason = \""~~Shortest[a__]~~"\";\n" :> a][[1]];
    Switch[
        terminateReason,
        "found primal feasible solution", True,
        "found dual feasible solution", False,
        _, Throw[terminateReason]]];

(* For a binary-valued function f, find the value of x in the interval
(true, false) where f changes from True to False, to within an error
of thresh. Returns the closest false value. *)
binarySearch[f_, true_, false_, thresh_] := Module[
    {test = (true + false)/2},
    If[Abs[true - false] < thresh,
       false,
       WriteString["stdout", "> trying: ", test, "\n"];
       If[f[test],
          binarySearch[f, test, false, thresh],
          binarySearch[f, true, test,  thresh]]]];

(* Compute an approximate upper bound on deltaPhi^2, as a function of deltaPhi. *)
bootstrapBound2d[derivativeOrder_, keptPoleOrder_, Lmax_] :=
    Table[
        WriteString["stdout", "> deltaPhi = ", deltaPhi, "\n"];
        {deltaPhi,
         binarySearch[
             singletAllowed2d[deltaPhi, #, derivativeOrder, keptPoleOrder, Lmax] &, 
             0.1,
             2,
             0.01]},
        {deltaPhi, 0.005,0.255,0.01}];

(* Uncomment this line to run the above computation (~8
minutes). Increase derivativeOrder to make a stronger bound, and be
sure to increase keptPoleOrder and Lmax until the results
stabilize. *)

(*Print[myBootstrapBound = bootstrapBound2d[7, 10, 15]];*)

(* Uncomment this line to plot the result.  It should display a kink
that moves closer to the 2d Ising point as derivativeOrder is
increased. *)

(* ListPlot[myBootstrapBound, GridLines -> {{1/8}, {1}}] *)
