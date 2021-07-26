$PROB Reference model 

$PARAM @annotated
TVVMAX : 10000 : Maximum elimination (mg/h)
TVVC   : 70.0 : Central volume of distribution (L)
TVKA   : 1.00 : Absorption rate (h-1)
TVKM   : 2500 : MM constant (mg/L)

ETA1 : 0 : VMAX
ETA2 : 0 : VC
ETA3 : 0 : KA
ETA4 : 0 : KM

$OMEGA
0.2 // VMAX 
0.2 // VC
0.2 // KA 
0.2 // KM 

$SIGMA 
0.05 // err prop
0   //  err additive 


$CMT @annotated
DEPOT   : Depot () []
CENTRAL : Central () [ADM, OBS]

$TABLE
double DV = (CENTRAL / VC) * (1 + EPS(1)) + EPS(2) ;

$MAIN
double VMAX= TVVMAX* exp(ETA(1) + ETA1 ) ; 
double VC  = TVVC  * exp(ETA(2) + ETA2 ) ;
double KA  = TVKA  * exp(ETA(3) + ETA3 ) ;
double KM  = TVKM  * exp(ETA(4) + ETA4 ) ;

$ODE
dxdt_DEPOT   = - KA * DEPOT ;
dxdt_CENTRAL =   KA * DEPOT - VMAX * (CENTRAL/VC) / (KM + (CENTRAL/VC)) ;

$CAPTURE DV