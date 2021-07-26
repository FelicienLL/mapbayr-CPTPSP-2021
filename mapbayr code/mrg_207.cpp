$PROB Reference model 

$PARAM @annotated
TVVMAX : 10000 : Maximum elimination (mg/h)
TVVC   : 70.0 : Central volume of distribution (L)
TVKA   : 1.00 : Absorption rate (h-1)
TVKM   : 2500 : MM constant (mg/L)
TVCL   : 4 : Clearance (L/h)

ETA1 : 0 : VMAX
ETA2 : 0 : VC
ETA3 : 0 : KA
ETA4 : 0 : KM
ETA5 : 0 : CL

$OMEGA
0.2 // VMAX 
0.2 // VC
0.2 // KA 
0.2 // KM 
0.2 // CL

$SIGMA 
0.05 // err prop
0   //  err additive 


$CMT @annotated
DEPOT   : Depot () [ADM]
CENTRAL : Central () [OBS]

$TABLE
double DV = (CENTRAL / VC) * (1 + EPS(1)) + EPS(2) ;

$MAIN
double VMAX= TVVMAX* exp(ETA(1) + ETA1 ) ; 
double VC  = TVVC  * exp(ETA(2) + ETA2 ) ;
double KA  = TVKA  * exp(ETA(3) + ETA3 ) ;
double KM  = TVKM  * exp(ETA(4) + ETA4 ) ;
double CL  = TVCL * exp(ETA(5) + ETA5) ;
double KE = CL / VC ; 

$ODE
dxdt_DEPOT   = - KA * DEPOT ;
dxdt_CENTRAL =   KA * DEPOT - VMAX * (CENTRAL/VC) / (KM + (CENTRAL/VC)) - KE * CENTRAL ;

$CAPTURE DV