$PROB Reference model 

$PARAM @annotated
TVCL : 4.0 : Clearance (L/h)
TVVC : 70.0 : Central volume of distribution (L)
TVKA : 1.00 : Absorption rate (h-1)
TVVMAX : 10000 : VMAX (mg/L)

ETA1 : 0 : CL
ETA2 : 0 : VC
ETA3 : 0 : KA
ETA4 : 0 : VMAX

$OMEGA
0.2 // CL 
0.2 // VC
0.2 // KA 
0.2 // VMAX

$SIGMA 
0.05 // err prop
0   //  err additive 


$CMT @annotated
DEPOT   : Depot () []
CENTRAL : Central () [ADM,OBS]

$TABLE
double DV = (CENTRAL / VC) * (1 + EPS(1)) + EPS(2) ;

$MAIN
double CL  = TVCL * exp(ETA(1) + ETA1 ) ; 
double VC  = TVVC  * exp(ETA(2) + ETA2 ) ;
double KA  = TVKA  * exp(ETA(3) + ETA3 ) ;
double VMAX  = TVVMAX  * exp(ETA(4) + ETA4 ) ;
double KE  = CL / VC ; 
double KM = 2500 ; 

$ODE
dxdt_DEPOT   = - KA * DEPOT ;
dxdt_CENTRAL =   KA * DEPOT - VMAX * (CENTRAL/VC) / (KM + (CENTRAL/VC)) - KE * CENTRAL ;

$CAPTURE DV