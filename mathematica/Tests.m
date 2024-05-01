(* Examples *)

<<"SDPB.m";

(* The following is the example in the manual *)
(* Maximize {a,b}.{0,-1} = -b over {a,b} such that {a,b}.{1,0}=a=1 and 

E^(-x)(a(1+x^4) + b(x^4/12 + x^2)) >= 0 for all x>=0

Equivalently,

1+x^4 + b(x^4/12 + x^2) >= 0 for all x >= 0

The prefactor DampedRational[1,{},1/E,x] doesn't affect the answer,
but does affect the choice of sample scalings and bilinear basis.

*)
testSDP[jsonFile_, prec_:200] := Module[
    {
        pols = {PositiveMatrixWithPrefactor[<|
        "prefactor"->DampedRational[1,{}, 1/E,x],
        "polynomials"->{{{1 + x^4, x^4/12 + x^2}}}
        |>]},
        norm = {1, 0},
        obj  = {0, -1}
    },
    
    (*
    If you want to specify sample points, sample scalings and/or bilinear bases explicitly,
    you may provide a function computing this data.
    See SDPB.m, getAnalyticSampleData[PositiveMatrixWithPrefactor[pmp_?AssociationQ],prec_]
    *)
    WritePmpJson[jsonFile, SDP[obj, norm, pols], prec(*, getAnalyticSampleData*)]
    ];

(* A similar computation to the above, except with nontrivial matrix semidefiniteness constraints *)
testSDPMatrix[jsonFile_, prec_:200] := Module[
    {
        pols = {
            PositiveMatrixWithPrefactor[<|
                "prefactor"->DampedRational[1, {}, 1/E, x],
                "polynomials"->
                {{{1 + x^4, 1 + x^4/12 + x^2}, {x^2,     x/5}},
                 {{x^2,     x/5},              {2 + x^4, x^4/3 + 2*x^2}}}
                 |>],
            PositiveMatrixWithPrefactor[<|
                "prefactor"->DampedRational[1, {}, 1/E, x],
                "polynomials"->
                {{{1 + 3x^4/4, 1 + x^4/12 + x^2}, {x^2,     1/2 + x/5}},
                 {{x^2,     1/2 + x/5},        {2 + 3x^4/5, x^4/3 + 2*x^2}}}
                 |>]
        },
        norm = {1, 0},
        obj  = {0, -1}
    },

    WritePmpJson[jsonFile, SDP[obj, norm, pols], prec(*, getAnalyticSampleData*)]
];
