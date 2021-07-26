$PROBLEM EXAMPLINIB1
$INPUT ID TIME EVID AMT CMT II ADDL MDV NDV DV S2_SAMPLING
$DATA data_to_fit407.csv IGNORE=@
$SUBROUTINES ADVAN13 TOL=8
$MODEL
COMP=(DEPOT) ; 1
COMP=(CENTRAL) ; 2

$PK 
iETA1 = 0
iETA2 = 0
iETA3 = 0
ETCL  = ETA(1) + iETA1 ; CL
ETVC  = ETA(2) + iETA2 ; VC
ETKA  = ETA(3) + iETA3 ; KA

TVCL   =  4.0
TVVC   =  70.0
TVKA   =  1.0

CL    = TVCL  * EXP(ETCL )
VC    = TVVC  * EXP(ETVC )
KA    = TVKA  * EXP(ETKA )
K20 = CL / VC

S2 = VC ; dv in mg/l ; amt in mg

$DES
DADT(1) = -KA * A(1)
DADT(2) =  KA * A(1) - K20*A(2)
$OMEGA
0.2 ; CL
0.2 ; VC
0.2 ; KA 
$SIGMA
0.05 FIX ; log additive

$ERROR
FLAG=0
IF(AMT.NE.0)FLAG=1  ;dosing records only
IPRED=LOG(F+FLAG)   ;transform the prediction to the log of the prediction
; IPRED=log(f) for concentration records and
; IPRED=log(f+1) for dose records
; be careful : IPRED value ok for event record only (not for adm lines). Please use NIPRED instead
W=1                 ;additive error model
Y= IPRED + W*EPS(1)
NIPRED = F

$EST METHOD=1 INTER NOABORT MAXEVAL=0 FORMAT= s1PE16.8E3
$COV UNCONDITIONAL


$TABLE ID TIME EVID AMT CMT II ADDL MDV NDV DV S2_SAMPLING NIPRED PRED IPRED NOPRINT NOAPPEND ONEHEADER FORMAT=s1PE16.8E3 FILE=run407.tab

