[[PROB]]
https://www.ncbi.nlm.nih.gov/pubmed/25870266


[[SET]] end = 365


[[ PARAM ]]

epsilon = 0.9    
pv2 = 0.65          
kdlta = 0.005    
n = 2              
c = 23             
delta = 1          
N = 2000           
           
alpha_L = 1e-6     
rho = 0.0045       
d_L = 0.004        
a = 0.001          
lambda = 1e4       
d = 0.01           
beta = 1.5e-8      
lambda_e = 1       
b_e = 1 
k_b = 0.1
d_e = 2
k_d = 5  
mu = 2  
m = 0.42 


[[ INIT ]]
art = 0, T = 874320,  L = 3, I = 1102, V = 95826, E = 0.3334

[[ MAIN ]]
double p = N*delta;

[[ ODE ]]

dxdt_art = 0;
double eps = art*epsilon;
dxdt_T = lambda - d*T - (1.0-eps)*beta*V*T;
dxdt_L = (alpha_L)*(1.0-eps)*beta*V*T + (rho - a - d_L)*L;
dxdt_I = (1.0-alpha_L)*(1.0-eps)*beta*V*T - delta*I + a*L - m*E*I;
dxdt_V = p*I  - c*V;
dxdt_E = lambda_e + (b_e)*(I/(k_b+I))*E - d_e*(I/(k_d+I))*E - mu*E;

[[ TABLE ]]
capture logV = log10(V);
capture logT = log10(T);
capture logI = log10(I);
capture logE = log10(E);
capture logL = log10(L);
capture year = TIME/365;

