(* Examples *)

<<"SDPB.m";

(* 2d Bootstrap Example *)

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
