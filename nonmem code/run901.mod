$PROBLEM CARBOPLATIN ADULTE 188
$INPUT ID TIME EVID AMT CMT RATE MDV DV BSA S2_SAMPLING
$DATA data_to_fit901.csv IGNORE=@
$SUBROUTINES ADVAN3 TRANS4

$PK 
iETA1 = 0
iETA2 = 0
iETA3 = 0
TVCL   = THETA(1)
TVV1   = THETA(2)
TVV2   = THETA(3)
TVQ    = THETA(4)

ETCL   = ETA(1) + iETA1
ETV1   = ETA(2) + iETA2
ETV2   = ETA(3) + iETA3

CL     = TVCL  * EXP(ETCL)
V1     = TVV1  * EXP(ETV1) * BSA
V2     = TVV2  * EXP(ETV2)
Q      = TVQ

S1 = V1

$ERROR
IPRED=F
W= F
Y = F + W * EPS(1)
IRES=DV-IPRED
IWRES=IRES/(W+0.001)

$THETA
(6.38) FIX ; 1  CL
(9.50) FIX ; 2  V1
(9.56) FIX ; 3  V2
(2.70) FIX ; 4  Q

$OMEGA 
0.1050 FIX ; 1 CL
0.0378 FIX ; 2 V1
0.0260 FIX ; 3 V2

$SIGMA 
0.0381 FIX ; ERR.PROP 

$EST METHOD=1 INTER NOABORT MAXEVAL=0 FORMAT= s1PE16.8E3
$COV UNCONDITIONAL


$TABLE ID TIME EVID AMT CMT RATE MDV DV BSA S2_SAMPLING PRED IPRED NOPRINT NOAPPEND ONEHEADER FORMAT=s1PE16.8E3 FILE=run901.tab