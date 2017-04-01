(* ::Package:: *)

(* Examples *)
SetDirectory[NotebookDirectory[]];
<<"SDPBschmudgen.m";

(* Simple Schmudgen example *)
(* Want to maxmimize b for x(x+b)(x-2) \[GreaterEqual]0 for all 1>=x>=0 and x\[GreaterEqual]2 *)
(* Should find primal-dual optimal value of -1 *)

(* Maximize {a,b}.{0, 1} = b over {a,b} such that {a,b}.{1,0}=a=1 and 

E^(-x)(a(x^3- 2x^2) + b( x^2 - 2x)) >= 0 for all 1>=x>=0 and x\[GreaterEqual]2

Equivalently,

x^3- 2x^2 + b(x^2 - 2x) >= 0 for all 1>=x>=0 and x\[GreaterEqual]2

The prefactor DampedRational[1,{},1/E,x] doesn't affect the answer,
but does affect the choice of sample scalings and bilinear basis.
*)

testSchmudgenSDP[datFile_] := Module[
    {
        pols = {PositiveMatrixWithPrefactor[DampedRational[1,{}, 1/E,x], {{{x^3- 2x^2, x^2 - 2x}}},{1,2}]},
		(*pols = {PositiveMatrixWithPrefactor[DampedRational[1,{}, 1/E,x], {{{x^3- x^2, -x^2+ x}}},{0.4,1.8}]},*)
        norm = {1, 0},
        obj  = {0, 1}
    },

    WriteBootstrapSDP[datFile, SDP[obj, norm, pols]];];



testSchmudgenSDP["testSchmudgen.xml"]
