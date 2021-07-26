$PROBLEM CABOZANTINIB NGUYEN
$INPUT ID TIME EVID AMT CMT II ADDL RATE MDV NDV DV AOLA AGE SEX WT RCC CRPC MTC GB OTHER HCC S2_SAMPLING
; NDV   : "normal" concentration mg/L
; DV    : log conc
; AMT   : mg
; AGE   : 60 yo
; SEX   : 0=male, 1=female
; WT    : 80 kg
; RCC   : renal cell cancer (0=no, 1=yes)
; CRPC  : castration resistant prostate cancer (0=no, 1=yes)
; MTC   : metast med thyroid cancer(0=no, 1=yes)
; GB    : glioblastoma (0=no, 1=yes)
; OTHER : other malignancies(0=no, 1=yes)
; HCC   : hepatocellular carcin (0=no, 1=yes)
; AOLA  : amont last administr (60 mg)
; LAG_AOLA  : lagged amont last administr (60 mg)

$DATA data_to_fit931.csv IGNORE=@

$SUBROUTINES ADVAN13 TOL=8
$MODEL
COMP=(DEPOT1) ; 1
COMP=(CENTRAL) ; 2
COMP=(PERIPH) ; 3

$PK 
iETA1 = 0
iETA2 = 0
iETA3 = 0
iETA4 = 0
ETKA   = ETA(1) + iETA1
ETCL   = ETA(2) + iETA2
ETV2   = ETA(3) + iETA3
ETFR   = ETA(4) + iETA4

TVKA    = THETA(1)
TVALAG1 = THETA(2)
TVCL    = THETA(3)
TVV2    = THETA(4)
TVQ     = THETA(5)
TVV3    = THETA(6)
LF1     = LOG(THETA(7) / ( 1 - THETA(7) )); Logit transform for F1
TVD2    = THETA(8)

KA     = TVKA    * EXP(ETKA) * ((AOLA/60)**THETA(9))
ALAG1  = TVALAG1

TVCL_COV =  TVCL  * ((AGE/60)**THETA(10)) * ((THETA(11))**SEX) * ((WT/80)**THETA(12)) * ((THETA(13))**RCC) * ((THETA(14))**CRPC) * ((THETA(15))**MTC) * ((THETA(16))**GB) * ((THETA(17))**OTHER) * ((THETA(18))**HCC)
CL     = TVCL_COV * EXP(ETCL)

TVV2_COV = TVV2 * ((AGE/60)**THETA(19)) * ((THETA(20))**SEX) * ((WT/80)**THETA(21)) * ((THETA(22))**RCC) * ((THETA(23))**CRPC) * ((THETA(24))**MTC) * ((THETA(25))**GB) * ((THETA(26))**OTHER) * ((THETA(27))**HCC)
V2     = TVV2_COV * EXP(ETV2)
Q      = TVQ
V3     = TVV3
F1     = EXP(LF1+ETFR)/(1+EXP(LF1+ETFR)) ; Fraction of dose absorbed via 1st order (from cmt1)
F2     = 1-F1 ; Fraction of dose absorbed via 0 order (from cmt2)
D2     = TVD2

S2 = V2 ; dv in mg/l ; amt in mg

K20 = CL / V2
K23 = Q  / V2
K32 = Q  / V3

ALAG2 = 0.000001

$DES

DADT(1) = -KA*A(1)
DADT(2) =  KA*A(1) + K32*A(3) - K23*A(2) - K20*A(2)
DADT(3) = -K32*A(3) + K23*A(2)

$THETA
(1.24    FIX) ;1  KA
(0.821   FIX) ;2  ALAG1
(2.48    FIX) ;3  CL
(212.0   FIX) ;4  V2
(30.0    FIX) ;5  Q
(177.0   FIX) ;6  V3
(0.83    FIX) ;7  F1
(2.48    FIX) ;8  D2
(0.734   FIX) ;9  KA_DOSE
(-0.157  FIX) ;10 CL_AGE
(0.76    FIX) ;11 CL_SEX
(-0.0393 FIX) ;12 CL_WT
(0.87    FIX) ;13 CL_RCC
(0.989   FIX) ;14 CL_CRPC
(1.9     FIX) ;15 CL_MTC
(1.2     FIX) ;16 CL_GB
(1.19    FIX) ;17 CL_OTHER
(0.878   FIX) ;18 CL_HCC
(0.0644  FIX) ;19 VC_AGE
(1.1     FIX) ;20 VC_SEX
(1.19    FIX) ;21 VC_WT
(0.656   FIX) ;22 VC_RCC
(0.743   FIX) ;23 VC_CRPC
(0.936   FIX) ;24 VC_MTC
(0.479   FIX) ;25 VC_GB
(0.762   FIX) ;26 VC_OTHER
(0.847   FIX) ;27 VC_HCC

$OMEGA
2.02 FIX        ;1 KA
$OMEGA BLOCK(2)
0.213 FIX       ;2 CL
0.211 0.443     ;3 V2
$OMEGA
2.55 FIX        ;4 F1

$SIGMA 0.127 FIX

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


$TABLE ID TIME EVID AMT CMT II ADDL RATE MDV NDV DV AOLA AGE SEX WT RCC CRPC MTC GB OTHER HCC S2_SAMPLING NIPRED PRED IPRED NOPRINT NOAPPEND ONEHEADER FORMAT=s1PE16.8E3 FILE=run931.tab

