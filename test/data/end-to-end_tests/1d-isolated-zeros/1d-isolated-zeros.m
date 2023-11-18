(* ::Package:: *)

SetDirectory@NotebookDirectory[];
<<"../../../../mathematica/SDPB.m";

(* The following is the modified example from the manual *)
(* Maximize {a,b}.{0,-1} = -b over {a,b} such that {a,b}.{1,0}=a=1 and 

E^(-x)(a(1+x^4) + b(x^4/12 + x^2)) >= 0 for all x>=2, and for x=2/3, x=4/3

Equivalently,

1+x^4 + b(x^4/12 + x^2) >= 0 for all x>=2, and for x=2/3, x=4/3

The prefactor DampedRational[1,{},1/E,x] doesn't affect the answer,
but does affect the choice of sample scalings and bilinear basis.

The resulting polynomial should have two zeros, one of which should be found by spectrum
*)
Module[
 {
  pols = {
    PositiveMatrixWithPrefactor[
     DampedRational[1, {}, 1/E, x],
     {{{1 + x^4, x^4/12 + x^2}}} /. x -> x + 2
     ],
    PositiveMatrixWithPrefactor[
     DampedRational[1, {}, 1/E, x],
     {{{1 + x^4, x^4/12 + x^2}}} /. x -> 2/3
     ],
    PositiveMatrixWithPrefactor[
     DampedRational[1, {}, 1/E, x],
     {{{1 + x^4, x^4/12 + x^2}}} /. x -> 4/3
     ]
    },
  norm = {1, 0},
  obj = {0, -1}
  },
 WritePmpJson["input/pmp.json", SDP[obj, norm, pols]]
 ]
