(* ::Package:: *)

SetDirectory@NotebookDirectory[];
<<"../../../../mathematica/SDPB.m";

If[DirectoryQ@"input",DeleteDirectory["input",DeleteContents->True]];
CreateDirectory@"input";
SetDirectory@"input";

(* The following is the modified example from the manual *)
(* Maximize {a,b}.{0,-1} = -b over {a,b} such that {a,b}.{1,0}=a=1 and 

E^(-x)(a(1+x^4) + b(x^4/12 + x^2)) >= 0 for all x>=2, and for x=2/3, x=4/3

Equivalently,

1+x^4 + b(x^4/12 + x^2) >= 0 for all x>=2, and for x=2/3, x=4/3

The prefactor DampedRational[1,{},1/E,x] doesn't affect the answer,
but does affect the choice of sample scalings and bilinear basis.

The resulting polynomial should have two zeros, one of which should be found by spectrum.

For testing purposes, we add several diagonal SDP matrices (1x1, 2x2, 3x3)
that represent the same constraint, f(x) >= 0 at x = 4/3.
*)
Module[{
prec=200
,norm = {1, 0}
,obj = {0, -1}
,poly0 = 1 + x^4
,poly1 = x^4 / 12 + x^2
,polyVector
,makeBlock43
,namedBlocks
,pmpNsv
}
,
polyVector = {poly0, poly1};

(*
Build (dim x dim) matrix = diag(f, 2*f,...rank*f,1,1...)
where f = 1+x^4 + y*(x^4/12+x^2) at x=4/3
and the rest are equal to 1.
spectrum should find eigenvectors corresponding to the first K=rank columns.
TODO: change basis to make it non-diagonal?
*)
makeBlock43[dim_,rank_]:={
  "constant_x=4%3_dim="<>ToString[dim]<>"_rank="<>ToString[rank]
  , DiagonalMatrix[Range[dim]]
  }/.{0->{0,0}}/.Table[
    i->If[i<=rank
    (*make different eigenvalues {1,2,3...}*const *)
    ,i*polyVector/.x->4/3
    ,{1,0}
    ]
  ,{i,dim}
];

namedBlocks={
{"continuum_x_gt_2",{{polyVector}}/.x->x+2}
,{"constant_x=2%3",{{polyVector}}/. x->2/3}
(*Different matrices describing the same isolated zero at x=4/3*)
,makeBlock43[1,1]
,makeBlock43[2,1]
,makeBlock43[2,2]
,makeBlock43[3,1]
,makeBlock43[3,2]
};
namedBlocks={
"json/"<>#[[1]]<>".json"
,#[[2]]
}&/@namedBlocks;

(*Write null-separated list of files*)
pmpNsv=OpenWrite@"pmp.nsv";
Table[
WriteString[pmpNsv,path];
WriteString[pmpNsv,FromCharacterCode[0]];
,{path,namedBlocks[[All,1]]}
];
Close@pmpNsv;

(*Write each block to a separate JSON file*)
Table[
WritePmpJson[
block[[1]]
,SDP[obj,norm,
{PositiveMatrixWithPrefactor[<|
"prefactor"->DampedRational[1,{}, 1/E,x],
"polynomials"->block[[2]]
|>]
}]
,prec
]
,{block,namedBlocks}
]
]



