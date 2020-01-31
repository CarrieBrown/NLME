
data simulated;

n_effects = 5;
x_min = 1;
x_max = 200;
x_int = 5;

a = 50;
b = 0.5;
c = 0.01;

ai = 5;
bi = 0.01;
ci = 0.001;

var_ai = .5;
var_bi = 0.001;
var_ci = 0.00001;
var_eu = 1;

do id = 1 to n_effects;
 ai_x = ai + sqrt(var_ai)*rand("Normal");
 bi_x = bi + sqrt(var_bi)*rand("Normal");
 ci_x = ci + sqrt(var_ci)*rand("Normal");
 
 do x = x_min to x_max by x_int;
  res = sqrt(var_eu)*rand("Normal");
  
   y = (a + ai) * (1 - (b + bi)*exp(-(c + ci) * x)) + res;
   output;
 end;
end;

drop a b c ai bi ci var_ai var_bi var_ci var_eu;

/* proc export data=simulated
    outfile='data.csv'
    dbms=csv
    replace;
run; */

goptions reset=all;
symbol i=join;
proc gplot data=simulated;
 plot y*x=id;
run;

data simulated;
  input x y mu;
datalines;
0 10.48 10
1 9.61 10.5
2 10.93 11
3 11.66 11.5
4 11.91 12
5 11.79 12.5
6 13.55 13
7 14 13.5
8 13.36 14
9 14.44 14.3
10 14.54 14.6
11 14.44 14.9
12 15.25 15
13 14.62 15
13 15.46 15
14 14.81 15
15 14.84 15
16 14.82 15
17 15.36 15
18 14.89 15
19 15.53 15
20 14.86 15
21 15.18 15
22 13.98 15
23 15.3 15
24 15.08 15
25 14.89 15
26 15.51 15
27 14.49 15
28 15.48 15
29 15.16 15
30 15.62 15
;
/* Figure 1 */
symbol1 value=' ' i=join color=blue;
symbol2 value=diamond i=none color=black;
proc gplot data=simulated;
 plot (mu y)*x/overlay;
 title1 "linear plataeu example";
 title2 "BLUE='truth'; BLACK diamonds = the observed data";
run;


proc nlin data=simulated list;
 parms a=0 b=0 c=0 ai=0 bi=0 ci=0;
 model y=(a + ai) * (1 - (b + bi)*exp(-(c + ci) * x));
run;

ods results off;

proc nlmixed data=simulated;
  parms a=50 b=0 to 1 c=0 logva=0 logvb=0 logvc=0.01 logeuv=0;
  mu=(a + ai) * (1 - (b + bi)*exp(-(c + ci) * x));
  model y ~ normal(mu, exp(logeuv));
  
  random ai bi ci ~ normal([0,0,0], [exp(logva),0,exp(logvb),0,0,exp(logvc)]) subject=id;

  predict mu out=nlm_pred;
  ods output ParameterEstimates = nlm_varest;
run;

proc export data=nlm_varest
outfile="nlm_varest.csv"
DBMS=DLM REPLACE;
DELIMITER=",";
run;

proc export data=nlm_pred
outfile="nlm_pred.csv"
DBMS=DLM REPLACE;
DELIMITER=",";
run;

proc iml;
    use simulated;
    read all;

    * Initial Values *;
    a = 10;
    b = 0.5;
    c = 0.05;

    ai_start = 0.1;
    bi_start = 0.5;
    ci_start = 0.9;

    sigma_residual = 0.5;
    sigma_random = {1, 0.1, 0.05};

    * Begin Program *;
    nobs = nrow(y);
    design = design(id);
    n_effects = ncol(design);
    n_sub = nobs/n_effects;
    crit = 1;
    niter = 0;
    print nobs n_effects n_sub;

    ai = j(n_effects,1,ai_start);
    bi = j(n_effects,1,bi_start);
    ci = j(n_effects,1,ci_start);

    beta_fixed = a//b//c;
    beta_random = ai//bi//ci;

    ai_x = ai@j(n_sub,1,1);
    bi_x = bi@j(n_sub,1,1);
    ci_x = ci@j(n_sub,1,1);

    fa = 1 - (b + bi_x) # exp(- (c + ci_x) # x);
    fb = (a + ai_x) # exp(- (c + ci_x) # x);
    fc = (a + ai_x) - (b + bi_x) # -x # exp(- (c + ci_x) # x);
    xstar = fa||fb||fc;
    
    fai = design#fa;
    fbi = design#fb;
    fci = design#fc;
    zstar = fai||fbi||fci;

    sigma_ai = i(n_effects)*sigma_random[1];
    sigma_bi = i(n_effects)*sigma_random[2];
    sigma_ci = i(n_effects)*sigma_random[3];

    g_side = block(sigma_ai, sigma_bi, sigma_ci);
    g_inv = inv(g_side);

    r_side = i(nobs)*sigma_residual;
    r_inv = inv(r_side);

    var_fun = zstar*g_side*zstar`+r_side;
    test = det(var_fun);
    print test;
    var_inv = sweep(var_fun);
    
    do while (crit>1e-12);
        yhat = (a + ai_x) # (1 - (b + bi_x) # exp(-(c + ci_x) # x));
        ystar = y - yhat + xstar*beta_fixed + zstar*beta_random;

        rss = ystar - xstar * inv(xstar` * var_inv * xstar) * xstar` * var_inv * ystar;
        log_PL = -0.5 * (log(det(var_fun)) + det(xstar` * var_inv * xstar) + rss` * var_inv * rss);

        lhs = ((xstar`*r_inv*xstar)||(xstar`*r_inv*zstar)) //
              ((zstar`*r_inv*xstar)||(zstar`*r_inv*zstar + g_inv));
        rhs = (xstar`*r_inv*ystar)//(zstar`*r_inv*ystar);
        solution = inv(lhs)*rhs;

        beta_fixed_new = solution[1:3];
        beta_random_new = solution[4:nrow(solution)];
        beta_random_matrix = (shape(beta_random_new,3,n_effects))`;

        a = beta_fixed_new[1];
        b = beta_fixed_new[2];
        c = beta_fixed_new[3];

        ai = beta_random_matrix[,1];
        bi = beta_random_matrix[,2];
        ci = beta_random_matrix[,3];

        ai_x = ai@j(n_sub,1,1);
        bi_x = bi@j(n_sub,1,1);
        ci_x = ci@j(n_sub,1,1);

        fa = 1 - (b + bi_x) # exp(- (c + ci_x) # x);
        fb = (a + ai_x) # exp(- (c + ci_x) # x);
        fc = (a + ai_x) - (b + bi_x) # -x # exp(- (c + ci_x) # x);
        xstar = fa||fb||fc;

        fai = design#fa;
        fbi = design#fb;
        fci = design#fc;
        zstar = fai||fbi||fci;

        var_fun = zstar*g_side*zstar`+r_side;
        var_inv = sweep(var_fun);

        p = var_inv-var_inv*xstar*sweep(xstar`*var_inv*xstar)*xstar`*var_inv;

        dv_ai = fai*fai`;
        dv_bi = fbi*fbi`;
        dv_ci = fci*fci`;
        dv_e = i(nobs); 

        sca = (-0.5)#trace(p*dv_ai) + 
            (1/2)#((ystar-xstar*beta_fixed_new)`*var_inv*dv_ai*var_inv*(ystar-xstar*beta_fixed_new));
        scb = (-0.5)#trace(p*dv_bi) + 
            (1/2)#((ystar-xstar*beta_fixed_new)`*var_inv*dv_bi*var_inv*(ystar-xstar*beta_fixed_new));
        scc = (-0.5)#trace(p*dv_ci) + 
            (1/2)#((ystar-xstar*beta_fixed_new)`*var_inv*dv_ci*var_inv*(ystar-xstar*beta_fixed_new));
        sce = (-0.5)#trace(p*dv_e) + 
            (1/2)#((ystar-xstar*beta_fixed_new)`*var_inv*dv_e*var_inv*(ystar-xstar*beta_fixed_new));
        score = sca//scb//scc//sce;

        h11=0.5#trace(p*dv_ai*p*dv_ai);
        h12=0.5#trace(p*dv_ai*p*dv_bi);
        h13=0.5#trace(p*dv_ai*p*dv_ci);
        h14=0.5#trace(p*dv_ai*p*dv_e);
        h21=0.5#trace(p*dv_bi*p*dv_ai);
        h22=0.5#trace(p*dv_bi*p*dv_bi);
        h23=0.5#trace(p*dv_bi*p*dv_ci);
        h24=0.5#trace(p*dv_bi*p*dv_e);
        h31=0.5#trace(p*dv_ci*p*dv_ai);
        h32=0.5#trace(p*dv_ci*p*dv_bi);
        h33=0.5#trace(p*dv_ci*p*dv_ci);
        h34=0.5#trace(p*dv_ci*p*dv_e);
        h41=0.5#trace(p*dv_e*p*dv_ai);
        h42=0.5#trace(p*dv_e*p*dv_bi);
        h43=0.5#trace(p*dv_e*p*dv_ci);
        h44=0.5#trace(p*dv_e*p*dv_e);  
        h = (h11 || h12 || h13 || h14) // 
            (h21 || h22 || h23 || h24) // 
            (h31 || h32 || h33 || h34) // 
            (h41 || h42 || h43 || h44) ;
        
        old_sigma = sigma_random // sigma_residual;
        sigma = old_sigma + inv(h)*score;
        if sigma[1,1]<0.0001 then sigma[1,1]=0.0001;
        if sigma[2,1]<0.0001 then sigma[2,1]=0.0001;
        if sigma[3,1]<0.000000001 then sigma[3,1]=0.000000001;
        sigma_random = sigma[1:3];
        print sigma_random;
        sigma_residual = sigma[4];

        sigma_ai = i(n_effects)*sigma_random[1];
        sigma_bi = i(n_effects)*sigma_random[2];
        sigma_ci = i(n_effects)*sigma_random[3];

        g_side = block(sigma_ai, sigma_bi, sigma_ci);
        g_inv = inv(g_side);

        r_side = i(nobs)*sigma_residual;
        r_inv = inv(r_side);

        var_fun = zstar*g_side*zstar`+r_side;
        test = det(var_fun);
        print test;
        var_inv = sweep(var_fun);

        yhat = (a + ai_x) # (1 - (b + bi_x) # exp(-(c + ci_x) # x));
        ystar = y - yhat + xstar*beta_fixed + zstar*beta_random;

        rss = ystar - xstar * inv(xstar` * var_inv * xstar) * xstar` * var_inv * ystar;
        new_log_PL = -0.5 * (log(det(var_fun)) + det(xstar` * var_inv * xstar) + rss` * var_inv * rss);

        crit = abs((new_log_PL - log_PL) / log_PL);
        log_PL = new_log_PL;
        beta_fixed = beta_fixed_new;
        beta_random = beta_random_new;
        niter = niter+1;
        
        print niter crit sigma old_sigma;
    end;

    print "Converaged after" niter "iternations";

    StdErrPred = sqrt(sigma_residual/n_sub);
    StdErrPred = j(nobs,1,StdErrPred);

    lower = yhat-1.96*StdErrPred;
    upper = yhat+1.96*StdErrPred;

    iml_pred = id || x || y || ystar || yhat || StdErrPred || lower || upper;

    iml_pred_colnames = {"id", "x", "y", "ystar", "Pred", "StdErrPred", "Lower", "Upper"};

    var_vi = sigma_random[1];
    var_ki = sigma_random[2];

    iml_varest = vm || km || var_vi || var_ki || sigma_residual;

    iml_varest_colnames = {"vm", "km", "var_vi", "var_ki", "var_res"};

    create iml_varest from iml_varest [colname=iml_varest_colnames];
    append from iml_varest;
    close iml_varest;

    create iml_pred from iml_pred [colname=iml_pred_colnames];
    append from iml_pred;
    close iml_pred;
run;

proc export data=iml_varest
outfile="iml_varest.csv"
DBMS=DLM REPLACE;
DELIMITER=",";
run;

proc export data=iml_pred
outfile="iml_pred.csv"
DBMS=DLM REPLACE;
DELIMITER=",";
run;


