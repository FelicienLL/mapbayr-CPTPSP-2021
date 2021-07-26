$PROBLEM EXAMPLINIB1
$INPUT ID TIME EVID AMT CMT II ADDL RATE MDV DV S2_SAMPLING
$DATA data_to_fit202.csv IGNORE=@
$SUBROUTINES ADVAN13 TOL=8
$MODEL
COMP=(DEPOT) ; 1
COMP=(CENTRAL) ; 2

$PK 
iETA1 = 0
iETA2 = 0
iETA3 = 0
iETA4 = 0
ETVMAX  = ETA(1) + iETA1 ; VMAX
ETVC  = ETA(2) + iETA2 ; VC
ETKA  = ETA(3) + iETA3 ; KA
ETKM  = ETA(4) + iETA4 ; KM

TVVMAX = 10000
TVVC   = 70.0
TVKA   = 1.0
TVKM   = 2500

VMAX  = TVVMAX* EXP(ETVMAX)
VC    = TVVC  * EXP(ETVC)
KA    = TVKA  * EXP(ETKA)
KM    = TVKM * EXP(ETKM)

S2 = VC ; dv in mg/l ; amt in mg

$DES
DADT(1) = -KA * A(1)
DADT(2) =  KA * A(1) - VMAX*(A(2)/S2)/(KM+(A(2)/S2))

$OMEGA
0.2 ; VMAX
0.2 ; VC
0.2 ; KA 
0.2 ; KM 
$SIGMA
0.05 FIX    ; prop
0 FIX ; additive

$ERROR
DEL=0
IF (F.EQ.0) DEL=1
W=F+DEL
Y=F+W*EPS(1)
IPRED=F
IRES=DV-IPRED
IWRES=IRES/W

$EST METHOD=1 INTER NOABORT MAXEVAL=0 FORMAT= s1PE16.8E3
$COV UNCONDITIONAL


$TABLE ID TIME EVID AMT CMT II ADDL RATE MDV DV S2_SAMPLING PRED IPRED NOPRINT NOAPPEND ONEHEADER FORMAT=s1PE16.8E3 FILE=run202.tab

