$title Indonesia-Energy Standard CGE Model in Ch. 6

$onText
Make an Indonesia-Energy Standard CGE Model
based on Chapter 6 of Hosoe, N, Gasawa, K, and Hashimoto, H
Handbook of Computible General Equilibrium Modeling
University of Tokyo Press, Tokyo, Japan, 2004
SAM in billion Rupiah unit

by Hadi Prasojo - for MSc Thesis - Economic Analysis, Corvinus University of Budapest, 2020
$offText

Set
   u    'SAM entry' / AFF, OIL, EMS, PIN, UGW, CON, VTI, OSV, CAP, LAB, IDT, HOH, GOV, INV, EXT /
   i(u) 'goods'     / AFF, OIL, EMS, PIN, UGW, CON, VTI, OSV                                    /
   h(u) 'factor'    /                                         CAP, LAB                          /;

Alias (u,v), (i,j), (h,k);

Table SAM(u,v) 'social accounting matrix'
           AFF      OIL      EMS      PIN      UGW      CON       VTI      OSV      CAP      LAB      IDT      HOH      GOV      INV      EXT
   AFF   64002        0       77   697881        0    37162     81724    14110                              393323        0   156249    38317
   OIL       0    16337    13900   193044       69       10        87      168                                   0        0     3735   101182
   EMS       3    23352    46609   244971    38301    85320       136     1377                                   0        0     8884   301574
   PIN   99971     7951    54082  1460592    47165   899197    599266   235701                             1962103    15154   447763  1081838
   UGW     543     1368      484    49165   145355     3798     21380    10788                               54048      228       57      592
   CON   22721     6578    15921     9269      521     8604     45373    44336                                   0        0  1600541     4696
   VTI    4258     5545    14543    85588     3552    48450    220081   155738                              811469      221     9842  1528924
   OSV   17695    18015    26507    77654     6184    48310    177198   167353                              637879   602575    29864    49042
   CAP  724052   151748   418485  1066983    70777   351335   1118344   554375
   LAB  246581    18280    95450   472878    20745   237196    494084   584862
   IDT    8325     1536     6459   195191   -50977    33778     25813    17833
   HOH                                                                          4456099  2170076
   GOV                                                                                             237958   385626
   INV                                                                                                     2381727     5406           -130198
   EXT  294694    77822    58010  2357567     6114     5400    104725    71635                                                                ;


* Loading the initial values ===================================================
Parameter
   Y0(j)   'composite factor'
   F0(h,j) 'the h-th factor input by the j-th firm'
   X0(i,j) 'intermediate input'
   Z0(j)   'output of the j-th good'
   Xp0(i)  'household consumption of the i-th good '
   Xg0(i)  'government consumption'
   Xv0(i)  'investment demand'
   E0(i)   'exports'
   M0(i)   'imports'
   Q0(i)   "Armington's composite good"
   D0(i)   'domestic good'
   Sp0     'private saving'
   Sg0     'government saving'
   Td0     'direct tax'
   Tz0(j)  'production tax'
   FF(h)   'factor endowment of the h-th factor'
   Sf      'foreign saving in US dollars'
   pWe(i)  'export price in US dollars'
   pWm(i)  'import price in US dollars'
   tauz(i) 'production tax rate'
;
Td0     = SAM("GOV","HOH");
Tz0(j)  = SAM("IDT",j);

F0(h,j) = SAM(h,j);
Y0(j)   = sum(h, F0(h,j));
X0(i,j) = SAM(i,j);
Z0(j)   = Y0(j) +sum(i, X0(i,j));
M0(i)   = SAM("EXT",i);

tauz(j) = Tz0(j)/Z0(j);

Xp0(i)  = SAM(i,"HOH");
FF(h)   = SAM("HOH",h);

Xg0(i)  = SAM(i,"GOV");
Xv0(i)  = SAM(i,"INV");
E0(i)   = SAM(i,"EXT");
Q0(i)   = Xp0(i)+Xg0(i)+Xv0(i)+sum(j, X0(i,j));
D0(i)   = (1+tauz(i))*Z0(i)-E0(i);
Sp0     = SAM("INV","HOH");
Sg0     = SAM("INV","GOV");
Sf      = SAM("INV","EXT");

pWe(i)  = 1;
pWm(i)  = 1;

display Y0, F0, X0, Z0, Xp0, Xg0, Xv0, E0, M0, Q0, D0, Sp0, Sg0, Td0, Tz0,
        FF, Sf, tauz;

* Calibration ==================================================================
Parameter
   sigma(i) 'elasticity of substitution'
   psi(i)   'elasticity of transformation'
   eta(i)   'substitution elasticity parameter'
   phi(i)   'transformation elasticity parameter';

sigma(i) =  2;
psi(i)   =  2;
eta(i)   = (sigma(i) - 1)/sigma(i);
phi(i)   = (psi(i) + 1)/psi(i);


Parameter
   alpha(i)  'share parameter in utility func.'
   beta(h,j) 'share parameter in production func.'
   b(j)      'scale parameter in production func.'
   ax(i,j)   'intermediate input requirement coeff.'
   ay(j)     'composite fact. input req. coeff.'
   mu(i)     'government consumption share'
   lambda(i) 'investment demand share'
   deltam(i) 'share par. in Armington func.'
   deltad(i) 'share par. in Armington func.'
   gamma(i)  'scale par. in Armington func.'
   xid(i)    'share par. in transformation func.'
   xie(i)    'share par. in transformation func.'
   theta(i)  'scale par. in transformation func.'
   ssp       'average propensity for private saving'
   ssg       'average propensity for gov. saving'
   taud      'direct tax rate';

alpha(i)  =  Xp0(i)/sum(j, Xp0(j));
beta(h,j) =  F0(h,j)/sum(k, F0(k,j));
b(j)      =  Y0(j)/prod(h, F0(h,j)**beta(h,j));

ax(i,j)   =  X0(i,j)/Z0(j);
ay(j)     =  Y0(j)/Z0(j);
mu(i)     =  Xg0(i)/sum(j, Xg0(j));
lambda(i) =  Xv0(i)/(Sp0+Sg0+Sf);

deltam(i) =  M0(i)**(1-eta(i))/(M0(i)**(1-eta(i)) + D0(i)**(1-eta(i)));
deltad(i) =  D0(i)**(1-eta(i))/(M0(i)**(1-eta(i)) + D0(i)**(1-eta(i)));
gamma(i)  =  Q0(i)/(deltam(i)*M0(i)**eta(i)+deltad(i)*D0(i)**eta(i))**(1/eta(i));

xie(i)    =  E0(i)**(1-phi(i))/(E0(i)**(1-phi(i))+D0(i)**(1-phi(i)));
xid(i)    =  D0(i)**(1-phi(i))/(E0(i)**(1-phi(i))+D0(i)**(1-phi(i)));
theta(i)  =  Z0(i)/(xie(i)*E0(i)**phi(i)+xid(i)*D0(i)**phi(i))**(1/phi(i));

ssp       =  Sp0/sum(h, FF(h));
ssg       =  Sg0/(Td0+sum(j, Tz0(j)));
taud      =  Td0/sum(h, FF(h));

display alpha, beta,  b,   ax,  ay, mu, lambda, deltam, deltad, gamma, xie
        xid,   theta, ssp, ssg, taud;

Variable
   Y(j)    'composite factor'
   F(h,j)  'the h-th factor input by the j-th firm'
   X(i,j)  'intermediate input'
   Z(j)    'output of the j-th good'
   Xp(i)   'household consumption of the i-th good'
   Xg(i)   'government consumption'
   Xv(i)   'investment demand'
   E(i)    'exports'
   M(i)    'imports'
   Q(i)    "Armington's composite good"
   D(i)    'domestic good'
   pf(h)   'the h-th factor price'
   py(j)   'composite factor price'
   pz(j)   'supply price of the i-th good'
   pq(i)   "Armington's composite good price"
   pe(i)   'export price in local currency'
   pm(i)   'import price in local currency'
   pd(i)   'the i-th domestic good price'
   epsilon 'exchange rate'
   Sp      'private saving'
   Sg      'government saving'
   Td      'direct tax'
   Tz(j)   'production tax'
   UU      'utility [fictitious]';

Equation
   eqpy(j)   'composite factor agg. func.'
   eqF(h,j)  'factor demand function'
   eqX(i,j)  'intermediate demand function'
   eqY(j)    'composite factor demand function'
   eqpzs(j)  'unit cost function'
   eqTd      'direct tax revenue function'
   eqTz(j)   'production tax revenue function'
   eqXg(i)   'government demand function'
   eqXv(i)   'investment demand function'
   eqSp      'private saving function'
   eqSg      'government saving function'
   eqXp(i)   'household demand function'
   eqpe(i)   'world export price equation'
   eqpm(i)   'world import price equation'
   eqepsilon 'balance of payments'
   eqpqs(i)  'Armington function'
   eqM(i)    'import demand function'
   eqD(i)    'domestic good demand function'
   eqpzd(i)  'transformation function'
   eqDs(i)   'domestic good supply function'
   eqE(i)    'export supply function'
   eqpqd(i)  'market clearing cond. for comp. good'
   eqpf(h)   'factor market clearing condition'
   obj       'utility function [fictitious]';

* domestic production
eqpy(j)..   Y(j)   =e= b(j)*prod(h, F(h,j)**beta(h,j));

eqF(h,j)..  F(h,j) =e= beta(h,j)*py(j)*Y(j)/pf(h);

eqX(i,j)..  X(i,j) =e= ax(i,j)*Z(j);

eqY(j)..    Y(j)   =e= ay(j)*Z(j);

eqpzs(j)..  pz(j)  =e= ay(j)*py(j) + sum(i, ax(i,j)*pq(i));

* government behavior
eqTd..      Td     =e= taud*sum(h, pf(h)*FF(h));

eqTz(j)..   Tz(j)  =e= tauz(j)*pz(j)*Z(j);

eqXg(i)..   Xg(i)  =e= mu(i)*(Td + sum(j, Tz(j)) - Sg)/pq(i);

* investment behavior
eqXv(i)..   Xv(i)  =e= lambda(i)*(Sp + Sg + epsilon*Sf)/pq(i);

* savings
eqSp..      Sp     =e= ssp*sum(h, pf(h)*FF(h));

eqSg..      Sg     =e= ssg*(Td + sum(j, Tz(j)));

* household consumption
eqXp(i)..   Xp(i)  =e= alpha(i)*(sum(h, pf(h)*FF(h)) - Sp - Td)/pq(i);

* international trade
eqpe(i)..   pe(i)  =e= epsilon*pWe(i);

eqpm(i)..   pm(i)  =e= epsilon*pWm(i);

eqepsilon.. sum(i, pWe(i)*E(i)) + Sf =e= sum(i, pWm(i)*M(i));

* Armington function
eqpqs(i)..  Q(i)   =e=  gamma(i)*(deltam(i)*M(i)**eta(i) + deltad(i)*D(i)**eta(i))**(1/eta(i));

eqM(i)..    M(i)   =e= (gamma(i)**eta(i)*deltam(i)*pq(i)/pm(i))**(1/(1-eta(i)))*Q(i);

eqD(i)..    D(i)   =e= (gamma(i)**eta(i)*deltad(i)*pq(i)/pd(i))**(1/(1-eta(i)))*Q(i);

* transformation function
eqpzd(i)..  Z(i)   =e=  theta(i)*(xie(i)*E(i)**phi(i)+xid(i)*D(i)**phi(i))**(1/phi(i));

eqE(i)..    E(i)   =e= (theta(i)**phi(i)*xie(i)*(1+tauz(i))*pz(i)/pe(i))**(1/(1-phi(i)))*Z(i);

eqDs(i)..   D(i)   =e= (theta(i)**phi(i)*xid(i)*(1+tauz(i))*pz(i)/pd(i))**(1/(1-phi(i)))*Z(i);

* market clearing condition
eqpqd(i)..  Q(i)   =e= Xp(i) + Xg(i) + Xv(i) + sum(j, X(i,j));

eqpf(h)..   sum(j, F(h,j)) =e= FF(h);

* fictitious objective function
obj..       UU     =e= prod(i, Xp(i)**alpha(i));

* Initializing variables
Y.l(j)    = Y0(j);
F.l(h,j)  = F0(h,j);
X.l(i,j)  = X0(i,j);
Z.l(j)    = Z0(j);
Xp.l(i)   = Xp0(i);
Xg.l(i)   = Xg0(i);
Xv.l(i)   = Xv0(i);
E.l(i)    = E0(i);
M.l(i)    = M0(i);
Q.l(i)    = Q0(i);
D.l(i)    = D0(i);
pf.l(h)   = 1;
py.l(j)   = 1;
pz.l(j)   = 1;
pq.l(i)   = 1;
pe.l(i)   = 1;
pm.l(i)   = 1;
pd.l(i)   = 1;
epsilon.l = 1;
Sp.l      = Sp0;
Sg.l      = Sg0;
Td.l      = Td0;
Tz.l(j)   = Tz0(j);

* Setting lower bounds to avoid division by zero
Y.lo(j)    = 0.00001;
F.lo(h,j)  = 0.00001;
X.lo(i,j)  = 0.0000;
Z.lo(j)    = 0.00001;
Xp.lo(i)   = 0.0000;
Xg.lo(i)   = 0.0000;
Xv.lo(i)   = 0.00001;
E.lo(i)    = 0.00001;
M.lo(i)    = 0.00001;
Q.lo(i)    = 0.00001;
D.lo(i)    = 0.00001;
pf.lo(h)   = 0.00001;
py.lo(j)   = 0.00001;
pz.lo(j)   = 0.00001;
pq.lo(i)   = 0.00001;
pe.lo(i)   = 0.00001;
pm.lo(i)   = 0.00001;
pd.lo(i)   = 0.00001;
epsilon.lo = 0.00001;
Sp.lo      = 0.00001;
Sg.lo      = 0.00001;
Td.lo      = 0.00001;
* Tz.lo(j)   = 0.0000;

* numeraire
pf.fx("LAB") = 1;

Model INAencge / all /;


** Simulation run: without any shocks (base run equilibrium) ===================
solve INAencge maximizing UU using nlp;

** Hypothetical simulation run: export oil price shocks (counterfactual eq)=====
pWe('OIL') = 0.7;
option bRatio = 1;

solve INAencge maximizing UU using nlp;


* Display of changes ===========================================================
Parameter
   dY(j), dF(h,j), dX(i,j), dZ(j),    dXp(i),     dXg(i), dXv(i)
   dE(i), dM(i),   dQ(i),   dD(i),    dpf(h),     dpy(j), dpz(i), dpq(i)
   dpe(i),dpm(i),  dpd(i),  depsilon, dTd,dTz(i), dSp,    dSg;

dY(j)    = (Y.l(j)  /Y0(j)  -1)*100;
dF(h,j)  = (F.l(h,j)/F0(h,j)-1)*100;
dX(i,j)$ (X0(i,j) <> 0)  = (X.l(i,j)/X0(i,j)-1)*100;
dZ(j)    = (Z.l(j)  /Z0(j)  -1)*100;
dXp(i) $ (Xp0(i) <> 0)  = (Xp.l(i) /Xp0(i) -1)*100;
dXg(i) $ (Xg0(i) <> 0)  = (Xg.l(i) /Xg0(i) -1)*100;
dXv(i)   = (Xv.l(i) /Xv0(i) -1)*100;
dE(i)    = (E.l(i)  /E0(i)  -1)*100;
dM(i)    = (M.l(i)  /M0(i)  -1)*100;
dQ(i)    = (Q.l(i)  /Q0(i)  -1)*100;
dD(i)    = (D.l(i)  /D0(i)  -1)*100;
dpf(h)   = (pf.l(h) /1 -1)*100;
dpy(j)   = (py.l(j) /1 -1)*100;
dpz(j)   = (pz.l(j) /1 -1)*100;
dpq(i)   = (pq.l(i) /1 -1)*100;
dpe(i)   = (pe.l(i) /1 -1)*100;
dpm(i)   = (pm.l(i) /1 -1)*100;
dpd(i)   = (pd.l(i) /1 -1)*100;
depsilon = (epsilon.l/1 -1)*100;
dTd      = (Td.l    /Td0    -1)*100;
dTz(j)   = (Tz.l(j) /Tz0(j) -1)*100;
dSp      = (Sp.l    /Sp0    -1)*100;
dSg      = (Sg.l    /Sg0    -1)*100;

display dY,  dF,  dX,  dZ,  dXp,      dXg, dXv, dE,  dM,  dQ, dD, dpf, dpy, dpz
        dpq, dpe, dpm, dpd, depsilon, dTd, dTz, dSp, dSg;


* Welfare measure: Hicksian equivalent variations ==============================
Parameter
   UU0 'utility level in the base run eq.'
   ep0 'expenditure func. in the base run eq.'
   ep1 'expenditure func. in the C-f eq.'
   EV  'Hicksian equivalent variations';

UU0 = prod(i, Xp0(i)**alpha(i));
ep0 = UU0 /prod(i, (alpha(i)/1)**alpha(i));
ep1 = UU.l/prod(i, (alpha(i)/1)**alpha(i));
EV  = ep1 - ep0;

display UU0, ep0, ep1, EV;
