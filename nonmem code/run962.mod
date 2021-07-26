$PROBLEM Voriconazole, Liu et al, 2014
$INPUT ID TIME EVID AMT CMT II ADDL MDV NDV DV WT S2_SAMPLING
$DATA data_to_fit962.csv IGNORE=@
$SUBROUTINES ADVAN13 TOL=8
$MODEL
COMP=(DEPOT) ; 1
COMP=(CENTRAL) ; 2
COMP=(PERIPH) ; 3

$PK 
iETA1 = 0
iETA2 = 0
iETA3 = 0
iETA4 = 0
iETA5 = 0
iETA6 = 0
iETA7 = 0
ETVMAXKM = ETA(1) + iETA1 ; Vmax and KM
ETCL     = ETA(2) + iETA2 ; CL
ETV2     = ETA(3) + iETA3 ; V2
ETV3     = ETA(4) + iETA4 ; V3
ETQ      = ETA(5) + iETA5 ; Q
ETF1     = ETA(6) + iETA6 ; F1
ETKA     = ETA(7) + iETA7 ; KA

TVKM      = 1150   ; KM constant (microg/L) // original 1.15 microg/mL
TVVMAX    = 114000 ; VMAX constant (microg/h/70kg)
VMAXINH   = 1.50   ; Maximum fract of Vmax inhibition ()
VMAXSCALE = 0.584  ; Scaling on random effect ()
T50       = 2.41   ; Time of half inhibition (h)
TVCL      = 6.16   ; Clearance (L/h/70kg)
TVV2      = 79.0   ; Central volume (L/70kg)
TVV3      = 103.0  ; Peripheral volume (L/70kg)
TVQ       = 15.5   ; Intercompartmental clearance (L/h/70kg)
TVF1      = 0.585  ; Typical bioavailability ()
BCF1      = 0.367  ; Box-Cox parameter F1 ()
TVKA      = 100  ; Absorption rate (h-1)
ALAG1     = 0.949  ; Lag time (h)


KM = TVKM * EXP(iETA1 + ETA(1)) ;
iVMAX = TVVMAX * EXP(VMAXSCALE * (iETA1 + ETA(1))) * ((WT/70) ** 0.75)
iVMAXINH  = EXP(VMAXINH) / (1 + EXP(VMAXINH)) 
CL        = TVCL * EXP(iETA2 + ETA(2)) * ((WT/70) ** 0.75) ;
V2        = TVV2 * EXP(iETA3 + ETA(3)) * WT/70 ;
V3        = TVV3 * EXP(iETA4 + ETA(4)) * WT/70 ;
Q         = TVQ  * EXP(iETA5 + ETA(5)) * ((WT/70) ** 0.75) ;
ETATR     = (EXP((ETA(6) + iETA6) * BCF1) - 1) / BCF1 ;
LTVF1     = LOG(TVF1 / (1 - TVF1)) ;
F1        = EXP(LTVF1 + ETATR) / (1 + EXP(LTVF1 + ETATR)) ;
KA        = TVKA * EXP(ETKA)

K20 = CL / V2 ;
K23 = Q / V2 ;
K32 = Q / V3 ;

S2 = V2 ; dv in microg/L ; amt in microg

$DES
VMAX = iVMAX * (1 - iVMAXINH * ((T - 1)/ ((T-1) + (T50 - 1)))) ;
CONC  = A(2) / V2 ;

DADT(1) = - KA * A(1) ;
DADT(2) =   KA * A(1) - (K20 + K23) * A(2) + K32 * A(3) - (VMAX * CONC)/(CONC + KM) ;
DADT(3) = K23 * A(2) - K32 * A(3) ;

$OMEGA
1.36  FIX ; VMAX_KM
0.239 FIX ; CL
0.136 FIX ; V2
0.769 FIX ; V3
0.424 FIX ; Q
0.686 FIX ; F1
0.898 FIX ; KA

$SIGMA
0 FIX ; err prop
0.017424 FIX ; err log additive

$ERROR
FLAG=0
IF(AMT.NE.0)FLAG=1  ;dosing records only
IPRED=LOG(F+FLAG)   ;transform the prediction to the log of the prediction
; IPRED=log(f) for concentration records and
; IPRED=log(f+1) for dose records
; be careful : IPRED value ok for event record only (not for adm lines). Please use NIPRED instead
W=1                 ;additive error model
Y= IPRED + W*EPS(2)
NIPRED = F

$EST METHOD=1 INTER NOABORT MAXEVAL=0 FORMAT= s1PE16.8E3
$COV UNCONDITIONAL


$TABLE ID TIME EVID AMT CMT II ADDL MDV NDV DV WT S2_SAMPLING NIPRED PRED IPRED NOPRINT NOAPPEND ONEHEADER FORMAT=s1PE16.8E3 FILE=run962.tab

