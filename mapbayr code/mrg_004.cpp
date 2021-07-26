$PROB Reference model 

$PARAM @annotated
TVCL   : 4.00 : Clearance (L/h)
TVVC   : 70.0 : Central volume of distribution (L)
TVD2   : 4.00 : Absorption rate (h-1)

ETA1 : 0 : CL
ETA2 : 0 : VC
ETA3 : 0 : D2

$OMEGA
0.2 // CL
0.2 // VC
0.2 // D2 

$SIGMA 
0.05 // err prop
0   //  err additive 


$CMT @annotated
DEPOT   : Depot () []
CENTRAL : Central () [ADM, OBS]

$TABLE
double DV = (CENTRAL / VC) * (1 + EPS(1)) + EPS(2) ;

$MAIN
double CL  = TVCL  * exp(ETA(1) + ETA1) ; 
double VC  = TVVC  * exp(ETA(2) + ETA2) ;
double D2  = TVD2  * exp(ETA(3) + ETA3) ;

double K20 = CL / VC ;

D_CENTRAL = D2 ; 

$ODE
dxdt_DEPOT   = 0 ;
dxdt_CENTRAL = - K20 * CENTRAL ;

$CAPTURE DV