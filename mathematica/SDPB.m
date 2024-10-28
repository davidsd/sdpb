(* ::Package:: *)

(* ::Section:: *)
(*Setup*)


(*TODO get rid of global prec, pass as an argument*)
(*Now it is used only for obsolete XML output*)
prec = 200;

(* DampedRational[c, {p1, p2, ...}, b, x] stands for c b^x / ((x-p1)(x-p2)...) *)
(* It satisfies the following identities *)

DampedRational[const_, poles_, base_, x + a_] := 
    DampedRational[base^a const, # - a & /@ poles, base, x];

DampedRational[const_, poles_, base_, a_ /; FreeQ[a, x]] := 
    const base^a/Product[a - p, {p, poles}];

DampedRational/:x DampedRational[const_, poles_ /; MemberQ[poles, 0], base_, x] :=
    DampedRational[const, DeleteCases[poles, 0], base, x];

DampedRational/:DampedRational[c1_,p1_,b1_,x] DampedRational[c2_,p2_,b2_,x] :=
    DampedRational[c1 c2, Join[p1, p2], b1 b2, x];

DampedRational/:a_ DampedRational[c_,p_,b_,x] /; FreeQ[a,x] := DampedRational[a c, p, b, x];

evalDampedRationalRegulated[DampedRational[c_,poles_,b_,x],x0_,minPoleDistance_]:=c b^x0/Product[Max[x0-p,minPoleDistance],{p,poles}];
evalDampedRationalRegulated[const_?NumericQ,x0_,minPoleDistance_]:=const;

nf[x_Integer, prec___] := x;
nf[x_, prec_:prec] := NumberForm[SetPrecision[x,prec],prec,ExponentFunction->(Null&)];
nf[x_] := nf[x ,prec];

safeCoefficientList[p_, x_] := Module[
    {coeffs = CoefficientList[p, x]},
    If[Length[coeffs] > 0, coeffs, {0}]];

(*WritePmpXml is equivalent to old WriteBootstrapSDP, writing XML output.
We want to enforce transition from XML to more effective and compact JSON,
thus we set WriteBootstrapSDP := WritePmpJson.
Note that users should change file extension from .xml to .json*)

(*WriteBootstrapSDP := WritePmpXml*)
(*WriteBootstrapSDP := Throw["WriteBootstrapSDP is deprecated, use WritePmpJson instead."];*)
WriteBootstrapSDP := WritePmpJson


(* ::Section:: *)
(*JSON export*)


toJsonNumber[x_, prec_] := ToString@nf[x,prec];

toJsonNumberArray[xs_List, prec_] := toJsonNumber[#,prec]& /@ xs;
toJsonNumberArray[xs_?MissingQ, args___] := xs;

toJsonObject[value_?MissingQ, args___]:=value;

bilinearBasisToJson[value_?MissingQ,args___]:=value;
bilinearBasisToJson[value_List,prec_]:=toJsonNumberArray[safeCoefficientList[#, x],prec]&/@value;

toJsonDampedRational[DampedRational[constant_, poles_List, base_, x],prec_] := <|
   "base" -> toJsonNumber[base,prec],
   "constant" -> toJsonNumber[constant,prec],
   "poles" -> toJsonNumberArray[poles,prec]
   |>;
   
toJsonDampedRational[constant_?NumericQ,prec_] := <|
   "base" -> toJsonNumber[1,prec],
   "constant" -> toJsonNumber[constant,prec],
   "poles" -> toJsonNumberArray[{},prec]
   |>;
toJsonDampedRational[value_?MissingQ, args___]:=value;


(*With default getSampleDataFn, "samplePoints" and next fields will be missing*)
toJsonObject[PositiveMatrixWithPrefactor[pmp_?AssociationQ], prec_, getSampleDataFn_:Function[<||>]] :=
Module[
{sampleData=getSampleDataFn[PositiveMatrixWithPrefactor[pmp],prec]}
,
<|
  "prefactor" -> toJsonDampedRational[pmp[["prefactor"]],prec],
  "reducedPrefactor" -> toJsonDampedRational[pmp[["reducedPrefactor"]],prec],
  "polynomials" -> Map[toJsonNumberArray[safeCoefficientList[#, x],prec] &, pmp[["polynomials"]], {3}],
  "samplePoints" -> toJsonNumberArray[sampleData[["samplePoints"]],prec],
  "sampleScalings" -> toJsonNumberArray[sampleData[["sampleScalings"]],prec],
  "reducedSampleScalings" -> toJsonNumberArray[sampleData[["reducedSampleScalings"]],prec],
  "bilinearBasis" -> bilinearBasisToJson[sampleData[["bilinearBasis"]],prec],
  "bilinearBasis_0" -> bilinearBasisToJson[sampleData[["bilinearBasis_0"]],prec],
  "bilinearBasis_1" -> bilinearBasisToJson[sampleData[["bilinearBasis_1"]],prec]
|>//DeleteMissing
  ];

toJsonObject[PositiveMatrixWithPrefactor[prefactor_, m_],args___] := 
toJsonObject[PositiveMatrixWithPrefactor[<|
"prefactor"->prefactor,
"polynomials"->m
|>],
args
];

toJsonObject[
  SDP[objective_List, normalization_List,
   positiveMatricesWithPrefactors_List], prec_, getSampleDataFn_:Function[<||>]
   ] := <|
  "objective" -> toJsonNumberArray[objective, prec],
  "normalization" -> toJsonNumberArray[normalization, prec],
  "PositiveMatrixWithPrefactorArray" ->
   Table[toJsonObject[pmp, prec, getSampleDataFn],{pmp,positiveMatricesWithPrefactors}]
  |>;

exportJson[file_,expr_]:=If[
  StringEndsQ[file,".json"],
  Export[file,expr,"JSON"],
  Throw["Expected .json extension: "<>ToString[file]]
];

(*getSampleDataFn computes sample points, sample scalings and bilinear bases.
getSampleDataFn[PositiveMatrixWithPrefactor[pmp],prec] should return an association
with some of the following (optional) fields:
<|
  "samplePoints" \[Rule]...
  "sampleScalings" \[Rule] ...
  "reducedSampleScalings" \[Rule] ...
  "bilinearBasis" \[Rule] ...
  "bilinearBasis_0" \[Rule] ...
  "bilinearBasis_1" \[Rule]...
|>
See getAnalyticSampleData[] below for an example.
By default, sampling data is not computed (pmp2sdp will compute it automatically).
*)
WritePmpJson[
  file_,
  SDP[objective_, normalization_, positiveMatricesWithPrefactors_],
  prec_, getSampleDataFn_:Function[<||>]
  ]:=exportJson[
    file,
    toJsonObject[SDP[objective, normalization, positiveMatricesWithPrefactors], prec, getSampleDataFn]
    ];


(* ::Section:: *)
(*Compute sample points, sample scalings and bilinear bases (same algorithm as in pmp2sdp)*)


(* ::Subsection:: *)
(*getAnalyticSampleData (main function)*)


(*
PMP input:
PolynomialMatrixWithPrefactor[<|
"prefactor"->DampedRational[...],
"reducedPrefactor"->DampedRational[...],
"polynomials"->{...}
|>]

Output:
<|
"samplePoints"\[Rule]...,
"sampleScalings"\[Rule]...,
"reducedSampleScalings"\[Rule]...,
"bilinearBasis_0"\[Rule]...,
"bilinearBasis_1"\[Rule]...
|>
*)

getAnalyticSampleData[PositiveMatrixWithPrefactor[pmp_?AssociationQ],prec_]:=Module[
{
numeratorDegOld,
prefactorOld,
prefactorNew
},
If[MissingQ[pmp[["polynomials"]]], Return[$Failed]];
numeratorDegOld=Max[Exponent[pmp[["polynomials"]], x]];

prefactorOld=pmp[["prefactor"]];
(*Set default prefactor to 1 for constant constraints and to Exp[-x] otherwise*)
If[MissingQ@prefactorOld,
  prefactorOld=If[numeratorDegOld==0,
    1,
    DampedRational[1,{},1/E,x]
  ]
];

prefactorNew=pmp[["reducedPrefactor"]];
If[MissingQ@prefactorNew, prefactorNew=prefactorOld];

getAnalyticSampleData[numeratorDegOld,prefactorOld,prefactorNew,prec]
];


(* ::Subsection:: *)
(*Implementation: sample points*)


(*
Compute Table[root of f=1/2+n, {n,nmin,NN-1}]. First we find the root with n=nmin, then we use that as an initial condition to find the next root, and so on. If we set nmin=0, then the total number of roots returned is NN.
*)
findBohrSommerfeldRoots[f_,nmin_,NN_,x_,x0_,prec_]:=findBohrSommerfeldHelper[f,NN,x,{},prec][nmin,x0];
(* Here, we use Mathematica's FindRoot. In the C++ implementation, perhaps we should use Newton's method (since we will know analytic formula for both the function f and its derivative)? *)
findBohrSommerfeldHelper[f_,NN_,x_,xs_,prec_][n_,x0_]:=If[
n>=NN,
xs,
Module[
{xPrime=x/.FindRoot[f==1/2+n,{x,x0},WorkingPrecision->prec]},
findBohrSommerfeldHelper[f,NN,x,Append[xs,xPrime],prec][n+1,xPrime]
]
];


(* Compute (nearly) optimal sample points for the given DampedRational. We use smallPoleThreshold to decide whether a pole is 0 or very close to 0. In that case, we include 0 as a sample point, and compute the remaining sample points starting from n=1. *)
getSamplePoints[const_?NumericQ,numSamplePoints_,smallPoleThreshold_,prec_]:=
If[numSamplePoints==1,{0},error["numSamplePoints is not equal to 1:",numSamplePoints]];

getSamplePoints[DampedRational[_,poles_,base_,x_],numSamplePoints_,smallPoleThreshold_,prec_]:=Module[
{
bVar,
bEquation,
bGuess,
b,
integratedDensity,
z,
z0,
nmin,
numSmallRoots,
smallRoots,
smallRootEnd,
bohrSommerfeldRoots,
(* Mathematica seems to require some quantities to be higher precision in order for FindRoot to give answers with precision prec *)
highPrec=2*prec,
result
},

If[numSamplePoints==1,Return[{0}]];

bEquation=Sum[1-Sqrt[-p/(bVar-p)],{p,poles}]-1/2 bVar Log[base]-numSamplePoints;
bGuess=-((2 numSamplePoints)/Log[base]);
(* Petr explains that to get rid of highPrec, we need to treat the p=0 case analytically in this sum *)
b=bVar/.FindRoot[bEquation,{bVar,bGuess},WorkingPrecision->highPrec];
(*TODO choose a better way to print warnings?*)
If[b<smallPoleThreshold,
  Print["b is too small, setting b=smallPoleThreshold"];
  b=smallPoleThreshold;
];

integratedDensity=Sum[ 1 /\[Pi] ( ArcCos[1-(2z(b-p))/(b (z-p))]- Sqrt[-p/(b-p)] ArcCos[1-(2 z)/b]),{p,poles}]-Log[base]/\[Pi] (Sqrt[(b-z) z]+ b/2 ArcCos[1-(2 z)/b]);

numSmallRoots=Count[poles,_?(Abs[#]<smallPoleThreshold&)];
numSmallRoots=Min[numSmallRoots,numSamplePoints];


z0=SetPrecision[smallPoleThreshold+(b-smallPoleThreshold)/(numSamplePoints-numSmallRoots+1.0),highPrec];

bohrSommerfeldRoots=findBohrSommerfeldRoots[integratedDensity,numSmallRoots,numSamplePoints,z,z0,prec];

smallRootEnd=If[numSmallRoots==numSamplePoints,b,bohrSommerfeldRoots[[1]]];


smallRoots=Table[
smallRootEnd*(i-1)/numSmallRoots
,{i,numSmallRoots}
];

result=Join[
smallRoots,
findBohrSommerfeldRoots[integratedDensity,numSmallRoots,numSamplePoints,z,z0,prec]
];

Table[
  Assert[result[[i+1]]>result[[i]]];
  ,{i,numSamplePoints-1}
];

result
];


(* ::Subsection:: *)
(*Implementation: sample scalings and bilinear bases*)


(* A matrix with constant anti-diagonals given by the list bs *)
antiBandMatrix[bs_] := Module[
    {n = Ceiling[Length[bs]/2]},
    Reverse[Normal[
        SparseArray[
            Join[
                Table[Band[{i, 1}] -> bs[[n - i + 1]], {i, n}],
                Table[Band[{1, i}] -> bs[[n + i - 1]], {i, 2, n}]],
            {n, n}]]]];


poleDegree[HoldPattern[DampedRational[_,poles_,_,_]]]:=Length[poles];
poleDegree[constant_?NumericQ]:=0;


(*getAnalyticSampleData[] implementation*)
getAnalyticSampleData[numeratorDegOld_,prefactorOld_,prefactorNew_,prec_]:=Module[
{
numeratorDeg,
numSamplePoints,
(* TODO: Is this a good value? *)
smallPoleThreshold=10^-10,
(* TODO: Is this a good value? *)
minPoleDistance=10^-16,
samplePoints,
sampleScalings,
reducedSampleScalings,
prefactorNewAtZero,
integrateMeasure1,
integrateMeasure2,
\[Delta]1,
\[Delta]2,
invL1,
invL2,
bilinearBasis1,
bilinearBasis2
},
numeratorDeg=numeratorDegOld-poleDegree[prefactorOld]+poleDegree[prefactorNew];
numSamplePoints=numeratorDeg+1;
samplePoints=getSamplePoints[prefactorNew,numSamplePoints,smallPoleThreshold,prec];

sampleScalings=Table[
evalDampedRationalRegulated[prefactorOld,xx,minPoleDistance],
{xx,samplePoints}
];

reducedSampleScalings=Table[
evalDampedRationalRegulated[prefactorNew,xx,minPoleDistance],
{xx,samplePoints}
];
integrateMeasure1[f_]:=Sum[
reducedSampleScalings[[i]](f/.x->samplePoints[[i]]),
{i,Length[samplePoints]}
];

integrateMeasure2[f_]:=Sum[
reducedSampleScalings[[i]](x* f/.x->samplePoints[[i]]),
{i,Length[samplePoints]}
];

\[Delta]1=Floor[numeratorDeg/2];
\[Delta]2=Floor[(numeratorDeg-1)/2];

invL1=Inverse[CholeskyDecomposition[antiBandMatrix[Table[
integrateMeasure1[x^n],
{n,0,2\[Delta]1}
]]]];
bilinearBasis1=Transpose[invL1] . Table[x^n,{n,0,\[Delta]1}];

If[\[Delta]2<0,
bilinearBasis2={};
,
invL2=
Inverse[CholeskyDecomposition[antiBandMatrix[Table[
integrateMeasure2[x^n],
{n,0,2\[Delta]2}
]]]];
bilinearBasis2=Transpose[invL2] . Table[x^n,{n,0,\[Delta]2}];
];

<|
"samplePoints"->samplePoints,
"sampleScalings"->sampleScalings,
"reducedSampleScalings"->reducedSampleScalings,
"bilinearBasis_0"->bilinearBasis1,
"bilinearBasis_1"->bilinearBasis2
|>
];


(* ::Section:: *)
(*XML export (obsolete)*)


(* ::Subsection:: *)
(*Helper functions*)


(* bilinearForm[f, m] = Integral[x^m f[x], {x, 0, Infinity}] *)
(* The special case when f[x] has no poles *)
bilinearForm[DampedRational[const_, {}, base_, x], m_] :=
    const Gamma[1+m] (-Log[base])^(-1-m);

(*memoizeGamma[a_,b_]:=memoizeGamma[a,b]=Gamma[a,b];*)

(* The case where f[x] has only single poles *)
(*bilinearForm[DampedRational[const_, poles_, base_, x], m_] := 
    const Sum[
        ((-poles[[i]])^m) ( base^poles[[i]]) Gamma[1 + m] memoizeGamma[-m, poles[[i]] Log[base]]/
        Product[poles[[i]] - p, {p, Delete[poles, i]}],
        {i, Length[poles]}];*)

(* The case where f[x] can have single or double poles *)
(*bilinearForm[DampedRational[c_, poles_, b_, x_], m_] := Module[
    {
        gatheredPoles = Gather[poles],
        quotientCoeffs = CoefficientList[PolynomialQuotient[x^m, Product[x-p, {p, poles}], x], x],
        integral, p, rest
    },
    integral[a_,1] := b^a Gamma[0, a Log[b]];
    integral[a_,2] := -1/a + b^a Gamma[0, a Log[b]] Log[b];
    c (Sum[
        p = gatheredPoles[[n,1]];
        rest = x^m / Product[x-q, {q, Join@@Delete[gatheredPoles, n]}];
        Switch[Length[gatheredPoles[[n]]],
               1, integral[p,1] rest /. x->p,
               2, integral[p,2] rest + integral[p,1] D[rest, x] /. x->p],
        {n, Length[gatheredPoles]}] + 
       Sum[
           quotientCoeffs[[n+1]] Gamma[1+n] (-Log[b])^(-1-n),
           {n, 0, Length[quotientCoeffs]-1}])];*)

(* A bilinearForm that allows for arbitrary collections of poles *)
bilinearForm[DampedRational[c_, poles_, b_, x_], m_] := Module[
    {
        gatheredPoles = GatherBy[poles, Round[#, 10^(-prec/2)]&],
        quotientCoeffs = CoefficientList[PolynomialQuotient[x^m,Product[x-p,{p,poles}],x],x],
        integral,
        rest,
        logRest,
        p,
        otherPoles,
        l
    },
    integral[0] := b^x Gamma[0,x Log[b]];
    integral[k_] := integral[k] = Simplify[D[integral[k-1],x]/k];
    dExp[k_] := dExp[k] = Expand[E^(-f[x])D[E^(f[x]),{x,k}]];
    c (
        Sum[
            p = gatheredPoles[[n,1]];
            l = Length[gatheredPoles[[n]]];
            Clear[otherPoles, logRest, rest];
            otherPoles = Join@@Delete[gatheredPoles, n];
            logRest[k_] := logRest[k] = (k-1)! (-1)^(k-1) (m / p^k - Sum[1/(p - q)^k, {q, otherPoles}]);
            rest[k_] := (dExp[k] /. { Derivative[n_][f][x] :> logRest[n] }) / k!;
            p^m / Product[p-q, {q, otherPoles}] * Sum[(integral[l-k-1]/.x->p) * rest[k], {k,0,l-1}],
            {n, Length[gatheredPoles]}] +
        Sum[quotientCoeffs[[n+1]] Gamma[1+n] (-Log[b])^(-1-n),
            {n,0,Length[quotientCoeffs]-1}]
      )];

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
rhoCrossing = SetPrecision[3-2 Sqrt[2], prec];

rescaledLaguerreSamplePoints[n_] := Table[
    SetPrecision[\[Pi]^2 (-1+4k)^2/(-64n Log[rhoCrossing]), prec],
    {k,0,n-1}];

maxIndexBy[l_,f_] := SortBy[
    Transpose[{l,Range[Length[l]]}],
    -f[First[#]]&][[1,2]];

(* finds v' such that a . v = First[v'] + a' . Rest[v'] when normalization . a == 1, where a' is a vector of length one less than a *)
reshuffleWithNormalization[normalization_, v_] := Module[
    {j = maxIndexBy[normalization, Abs], const},
    const = v[[j]]/normalization[[j]];
    Prepend[Delete[v - normalization*const, j], const]];


(* ::Subsection::Closed:: *)
(*Export to XML file*)


WritePmpXml[
    file_,
    SDP[objective_, normalization_, positiveMatricesWithPrefactors_],
    samplePointsFn_ : rescaledLaguerreSamplePoints] := Module[
    {
        stream = OpenWrite[file],
        node, real, int, vector, polynomial,
        polynomialVector, polynomialVectorMatrix,
        affineObjective, polynomialVectorMatrices
    },
    If[!StringEndsQ[file,".xml"],
      Throw["Expected .xml extension: "<>ToString[file]];
      ];

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

        samplePoints   = samplePointsFn[degree + 1];
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
