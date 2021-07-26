$PROB LAG model 

$PARAM @annotated
TVCL   : 4.00 : Clearance (L/h)
TVVC   : 70.0 : Central volume of distribution (L)
TVKA1  : 1.00 : Absorption rate 1 (h-1)
TVKA3  : .25  : Absorption rate 3 (h-1)
FR     : 0.2  : Fraction absorbed from Depot 1 ()

ETA1 : 0 : CL
ETA2 : 0 : VC
ETA3 : 0 : KA1
ETA4 : 0 : KA3

$OMEGA
0.2 // CL
0.2 // VC
0.2 // KA1
0.2 // KA3

$SIGMA 
0.05 // err prop
0   //  err additive 


$CMT @annotated
DEPOT1 : Depot 1 () [ADM]
CENTRAL : Central () [OBS]
DEPOT3 : Depot 3 () [ADM]

$TABLE
double DV = (CENTRAL / VC) * (1 + EPS(1)) + EPS(2) ;

$MAIN
double CL  = TVCL   * exp(ETA(1) + ETA1) ; 
double VC  = TVVC   * exp(ETA(2) + ETA2) ;
double KA1 = TVKA1  * exp(ETA(3) + ETA3) ;
double KA3 = TVKA3  * exp(ETA(4) + ETA4) ;
double K20 = CL / VC ;

F_DEPOT1 = FR ;
F_DEPOT3 = 1 - FR ; 

$ODE
dxdt_DEPOT1 = - KA1 * DEPOT1 ;
dxdt_CENTRAL = - K20 * CENTRAL + KA1 * DEPOT1 + KA3 * DEPOT3 ;
dxdt_DEPOT3 = - KA3 * DEPOT3 ;

$CAPTURE DV