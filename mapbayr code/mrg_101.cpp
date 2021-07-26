$PROB Reference model 

$PARAM @annotated
TVCL   : 4.00 : Clearance (L/h)
TVVC   : 70.0 : Central volume of distribution (L)
TVKA   : 1.00 : Absorption rate (h-1)
TVVP   : 50.0 : Peripheral volume ()
Q      : 4.00 : Intercompartmental clearance (L/h)

ETA1 : 0 : CL
ETA2 : 0 : VC
ETA3 : 0 : KA
ETA4 : 0 : VP

$OMEGA
0.2 // CL
0.2 // VC
0.2 // KA 
0.2 // VP

$SIGMA 
0.05 // err prop
0   //  err additive 

$CMT @annotated
DEPOT   : Depot () [ADM]
CENTRAL : Central () [OBS]
PERIPH  : Peripheral ()

$TABLE
double DV = (CENTRAL / VC) * (1 + EPS(1)) + EPS(2) ;

$MAIN
double CL  = TVCL  * exp(ETA(1) + ETA1 ) ; 
double VC  = TVVC  * exp(ETA(2) + ETA2 ) ;
double KA  = TVKA  * exp(ETA(3) + ETA3 ) ;
double VP  = TVVP  * exp(ETA(4) + ETA4 ) ;

double K23 = Q / VC ; 
double K32 = Q / VP ;
double K20 = CL / VC ;

$ODE
dxdt_DEPOT  = - KA * DEPOT ;
dxdt_CENTRAL = - (K20 + K23) * CENTRAL + K32 * PERIPH + KA * DEPOT ;
dxdt_PERIPH = K23 * CENTRAL - K32 * PERIPH ;

$CAPTURE DV