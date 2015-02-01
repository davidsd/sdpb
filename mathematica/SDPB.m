(* ::Package:: *)

(*Setup*)
prec = 200;

(* A matrix with constant anti-diagonals given by the list bs *)
antiBandMatrix[bs_] := Module[
    {n = Ceiling[Length[bs]/2]},
    Reverse[Normal[
        SparseArray[
            Join[
                Table[Band[{i, 1}] -> bs[[n - i + 1]], {i, n}],
                Table[Band[{1, i}] -> bs[[n + i - 1]], {i, 2, n}]],
            {n, n}]]]];

(* DampedRational[c, {p1, p2, ...}, b, x] stands for c b^x / ((x-p1)(x-p2)...) *)
(* It satisfies the following identities *)

DampedRational[const_, poles_, base_, x + a_] := 
    DampedRational[base^a const, # - a & /@ poles, base, x];

DampedRational[const_, poles_, base_, a_ /; FreeQ[a, x]] := 
    const base^a/Product[a - p, {p, poles}];

DampedRational/:x DampedRational[const_, poles_ /; MemberQ[poles, 0], base_, x] :=
    DampedRational[const, DeleteCases[poles, 0], base, x];


(* bilinearForm[f, m] = Integral[x^m f[x], {x, 0, Infinity}] *)
(* The special case when f[x] has no poles *)
bilinearForm[DampedRational[const_, {}, base_, x], m_] :=
    const Gamma[1+m] (-Log[base])^(-1-m);

memoizeGamma[a_,b_]:=memoizeGamma[a,b]=Gamma[a,b];

(* The general DampedRational case *)
bilinearForm[DampedRational[const_, poles_, base_, x], m_] := 
    const Sum[
        ((-poles[[i]])^m) ( base^poles[[i]]) Gamma[1 + m] memoizeGamma[-m, poles[[i]] Log[base]]/
        Product[poles[[i]] - p, {p, Delete[poles, i]}],
        {i, Length[poles]}];

(* orthogonalPolynomials[f, n] is a set of polynomials with degree 0
through n which are orthogonal with respect to the measure f[x] dx *)
orthogonalPolynomials[const_ /; FreeQ[const, x], 0] := {1/Sqrt[const]};

orthogonalPolynomials[const_ /; FreeQ[const, x], degree_] := 
    error["can't get orthogonal polynomials of nonzero degree for constant measure"];

orthogonalPolynomials[DampedRational[const_, poles_, base_, x], degree_] := 
    Table[x^m, {m, 0, degree}] . Inverse[
        CholeskyDecomposition[
            antiBandMatrix[
                Table[bilinearForm[DampedRational[const, Select[poles, # < 0&], base, x], m],
                      {m, 0, 2 degree}]]]];

(* Preparing SDP for Export *)
rho = SetPrecision[3-2 Sqrt[2], prec];

rescaledLaguerreSamplePoints[n_] := Table[
    SetPrecision[\[Pi]^2 (-1+4k)^2/(-64n Log[rho]), prec],
    {k,0,n-1}];

maxIndexBy[l_,f_] := SortBy[
    Transpose[{l,Range[Length[l]]}],
    -f[First[#]]&][[1,2]];

(* finds v' such that a . v = First[v'] + a' . Rest[v'] when normalization . a == 1, where a' is a vector of length one less than a *)
reshuffleWithNormalization[normalization_, v_] := Module[
    {j = maxIndexBy[normalization, Abs], const},
    const = v[[j]]/normalization[[j]];
    Prepend[Delete[v - normalization*const, j], const]];

(* XML Exporting *)
nf[x_Integer] := x;
nf[x_] := NumberForm[SetPrecision[x,prec],prec,ExponentFunction->(Null&)];

safeCoefficientList[p_, x_] := Module[
    {coeffs = CoefficientList[p, x]},
    If[Length[coeffs] > 0, coeffs, {0}]];

WriteBootstrapSDP[file_, SDP[objective_, normalization_, positiveMatricesWithPrefactors_]] := Module[
    {
        stream = OpenWrite[file],
        node, real, int, vector, polynomial,
        polynomialVector, polynomialVectorMatrix,
        affineObjective, polynomialVectorMatrices
    },

    (* write a single XML node to file.  children is a routine that writes child nodes when run. *)
    node[name_, children_] := (
        WriteString[stream, "<", name, ">"];
        children[];
        WriteString[stream, "</", name, ">\n"];
    );

    real[r_][] := WriteString[stream, nf[r]];
    int[i_][] := WriteString[stream, i];
    vector[v_][] := Do[node["elt", real[c]], {c, v}];
    polynomial[p_][] := Do[node["coeff", real[c]], {c, safeCoefficientList[p,x]}];
    polynomialVector[v_][] := Do[node["polynomial", polynomial[p]], {p, v}];

    polynomialVectorMatrix[PositiveMatrixWithPrefactor[prefactor_, m_]][] := Module[
        {degree = Max[Exponent[m, x]], samplePoints, sampleScalings, bilinearBasis},

        samplePoints   = rescaledLaguerreSamplePoints[degree + 1];
        sampleScalings = Table[prefactor /. x -> a, {a, samplePoints}];
        bilinearBasis  = orthogonalPolynomials[prefactor, Floor[degree/2]];

        node["rows", int[Length[m]]];
        node["cols", int[Length[First[m]]]];
        node["elements", Function[
            {},
            Do[node[
                "polynomialVector",
                polynomialVector[reshuffleWithNormalization[normalization,pv]]],
               {row, m}, {pv, row}]]];
        node["samplePoints", vector[samplePoints]];
        node["sampleScalings", vector[sampleScalings]];
        node["bilinearBasis", polynomialVector[bilinearBasis]];
    ];

    node["sdp", Function[
        {},
        node["objective", vector[reshuffleWithNormalization[normalization, objective]]];
        node["polynomialVectorMatrices", Function[
            {},
            Do[node["polynomialVectorMatrix", polynomialVectorMatrix[pvm]], {pvm, positiveMatricesWithPrefactors}];
        ]];
    ]];                                          

    Close[stream];
];

(*Test Computations*)

(* Maximize {a,b}.{0,-1} = -b over {a,b} such that {a,b}.{1,0}=a=1 and 

E^(-x)(a(1+x^4) + b(x^4/12 + x^2)) >= 0 for all x>=0

Equivalently,

1+x^4 + b(x^4/12 + x^2) >= 0 for all x >= 0

The prefactor DampedRational[1,{},1/E,x] doesn't affect the answer,
but does affect the choice of sample scalings and bilinear basis.

*)
testSDP[datFile_] := Module[
    {
        pols = {PositiveMatrixWithPrefactor[DampedRational[1,{}, 1/E,x], {{{1 + x^4, x^4/12 + x^2}}}]},
        norm = {1, 0},
        obj  = {0, -1}
    },

    WriteBootstrapSDP[datFile, SDP[obj, norm, pols]];];

(* A similar computation to the above, except with nontrivial matrix semidefiniteness constraints *)
testSDPMatrix[datFile_] := Module[
    {
        pols = {
            PositiveMatrixWithPrefactor[
                DampedRational[1, {}, 1/E, x],
                {{{1 + x^4, 1 + x^4/12 + x^2}, {x^2,     x/5}},
                 {{x^2,     x/5},              {2 + x^4, x^4/3 + 2*x^2}}}],
            PositiveMatrixWithPrefactor[
                DampedRational[1, {}, 1/E, x],
                {{{1 + 3x^4/4, 1 + x^4/12 + x^2}, {x^2,     1/2 + x/5}},
                 {{x^2,     1/2 + x/5},        {2 + 3x^4/5, x^4/3 + 2*x^2}}}]
        },
        norm = {1, 0},
        obj  = {0, -1}
    },

    WriteBootstrapSDP[datFile, SDP[obj, norm, pols]];];

(*2d Bootstrap Example*)

tableDir="/data/dsd/sdpa-multicorrelators/tables-zzb-singlespin";
spacetimeDim=3;

prec=200;

error[msg_] := (Print[msg]; Exit[];);

floatToString[x_]:=If[Abs[x]<10^(-10),"0", ToString[CForm[SetAccuracy[x,10]]]];

(* Fail with an error message if a file doesn't exist *)
safeGet[file_] := If[FileExistsQ[file], Get[file], error[file <> " does not exist. Exiting."];];

(* Fail with an error message if a file doesn't exist *)
safeImport[file_] := If[FileExistsQ[file], Import[file], error[file <> " does not exist. Exiting."];];

(* Memoize rcDerivPolTable so we don't fetch files multiple times *)
(* This is the old basic table with thresholds instead of keptPoleOrder *)
rcDerivPolTable[nmax_] := (
    rcDerivPolTable[nmax] = safeGet[
        $HomeDirectory<>"/Dropbox/SDP3D/rho-expansion/rcDerivTable-nmax"<>ToString[nmax]<>"-thresh10E-2-order60-shifting-allL.m"]);

(* These are the newer tables with keptPoleOrder *)
derivPolTable[d_, delta12_, delta34_, L_, nmax_, keptPoleOrder_, order_] :=
    safeImport[
        StringJoin[
            tableDir,
            "/nmax",            ToString[nmax],
            "/zzbDerivTable-d", ToString[N[d]],
            "-delta12-",        floatToString[delta12],
            "-delta34-",        floatToString[delta34],
            "-L",               ToString[L],
            "-nmax",            ToString[nmax],
            "-keptPoleOrder",   ToString[keptPoleOrder],
            "-order",           ToString[order], 
            ".mx"]];

(* Computation Setup *)

seriesDataToRules[ruleName_,series_] := Module[
    {toRule},
    toRule[expr_,{m_,n_}] := ruleName[m-1,n-1]->(m-1)!(n-1)!expr;
    toRule[expr_,_]       := expr;
    Flatten[MapIndexed[toRule,Normal[series],2]]
];

(* Here we use that z and zb are symmetric *)
withDeltaPhiToWithoutDeltaPhi[DeltaPhi_][mMax_,nMax_] := Module[
    {series, half},
    half = SetPrecision[1/2,prec];
    series = Series[
        ((1-z) (1-zb))^DeltaPhi f[z,zb]/.{z->half+dz,zb->half+dzb},
    {dz,0,mMax},
    {dzb,0,nMax}
    ] // Normal;
seriesDataToRules[withDeltaPhiDeriv,
    CoefficientList[series,{dz,dzb}]/.{
        Derivative[j_,k_][f][_,_]:>zzbDeriv[Max[j,k],Min[j,k]],
        f[_,_]:>zzbDeriv[0,0]
                                      }
    ]
];

leadingCoefficient[pols_,x_]:=Coefficient[pols,x^Max[Exponent[pols,x]]];

oddDerivs[nMax_]  := Flatten[Table[withDeltaPhiDeriv[n+i-1,n-i]/((n+i-1)!(n-i)!),{n,1,nMax},{i,1,n}]];
evenDerivs[nMax_] := Flatten[Table[withDeltaPhiDeriv[n+i-2,n-i]/((n+i-2)!(n-i)!),{n,1,nMax},{i,1,n}]];

singletSpectrumDisallowed[datFile_][deltaSigLowPrec_,deltaEpsLowPrec_,nmax_,Lmax_] := Module[
    {
        coeffs,unitOp,norm,obj,
        DeltaScalarShift=SetPrecision[deltaEpsLowPrec-1,prec],
        DeltaPhi=SetPrecision[deltaSigLowPrec,prec],
        addDeltaPhi,
        pols
    },

    addDeltaPhi=withDeltaPhiToWithoutDeltaPhi[DeltaPhi][2nmax,2nmax];

    (* FIXME! *)
    pols=Table[
        PositiveMatrixWithPrefactor[DampedRational[1,{},rho,x], {{
            Expand[
                Expand[
                    Expand[
                        oddDerivs[nmax]/.addDeltaPhi
                          ]/.zzbToRcTable[nmax]
                      ]/.rcDerivPolTable[nmax][[L+1]]
                  ]
         }}],
        {L,0,Lmax,2}
    ];

    unitOp=oddDerivs[nmax]/.addDeltaPhi/.{
        zzbDeriv[0,0]->1,
        zzbDeriv[_,_]->0
    };

    pols = MapAt[#/.x->x+DeltaScalarShift&,pols,1];
    norm = unitOp;
    obj = 0*unitOp;
    WriteBootstrapSDP[datFile,SDP[obj,norm,pols]];
];
