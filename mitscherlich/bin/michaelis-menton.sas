
data simulated;
n_effects = 12;
x_min = 1;
x_max = 70;
x_int = 5;

vm = 9.6;
km = 30;

vi = 1;
ki = 0.1;

var_vi = 0.0001;
var_ki = .8;
var_eu = .1;

do id = 1 to n_effects;
 vi_x = vi + sqrt(var_vi)*rand("Normal");
 ki_x = ki + sqrt(var_ki)*rand("Normal");
 
 do x = x_min to x_max by x_int;
  res = sqrt(var_eu)*rand("Normal");
  
   y = ((vm + vi) * x) / (km + ki + x) + res;
   output;
 end;
end;

proc export data=simulated
    outfile='data.csv'
    dbms=csv
    replace;
run;

/* 
goptions reset=all;
symbol i=join;
proc gplot data=simulated;
 plot y*x=id;
run;
 */

ods results off;

proc nlmixed data=simulated;
  parms vmax=10 km=30 sdv=1 sdk=0.1 sdresid=1;
  mu=((vmax + vi) * x)/(km + ki + x);
  model y ~ normal(mu, sdresid*sdresid);
  random vi ki ~ normal([0,0],[sdv*sdv,0,sdk*sdk]) subject=id;

  predict eta out=nlm_pred;
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
    vm = 90;
    km = 10;

    vi_start = 0.1;
    ki_start = 0.1;

    sigma_residual = 14.1;
    sigma_random = {0.1, 6.7};

    * Begin Program *;
    nobs = nrow(y);
    design = design(id);
    n_effects = ncol(design);
    n_sub = nobs/n_effects;
    crit = 1;
    niter = 0;
    print nobs n_effects n_sub;

    vi = j(n_effects,1,vi_start);
    ki = j(n_effects,1,ki_start);

    beta_fixed = vm//km;
    beta_random = vi//ki;

    vi_x = vi@j(n_sub,1,1);
    ki_x = ki@j(n_sub,1,1);

    fvm = x / (km + ki_x + x);
    fkm = -((vm + vi_x) # x) / ((km + ki_x + x) ## 2);
    xstar = fvm||fkm;

    fvi = design#fvm;
    fki = design#fkm;
    zstar = fvi||fki;

    sigma_vi = i(n_effects)*sigma_random[1];
    sigma_ki = i(n_effects)*sigma_random[2];

    g_side = block(sigma_vi, sigma_ki);
    g_inv = inv(g_side);

    r_side = i(nobs)*sigma_residual;
    r_inv = inv(r_side);

    var_fun = zstar*g_side*zstar`+r_side;
    var_inv = inv(var_fun);

    do while (crit>1e-12);
        yhat = ((vm + vi_x) # x) / (km + ki_x + x);
        ystar = y - yhat + xstar*beta_fixed + zstar*beta_random;

        rss = ystar - xstar * inv(xstar` * var_inv * xstar) * xstar` * var_inv * ystar;
        log_PL = -0.5 * (log(det(var_fun)) + det(xstar` * var_inv * xstar) + rss` * var_inv * rss);

        lhs = ((xstar`*r_inv*xstar)||(xstar`*r_inv*zstar)) //
            ((zstar`*r_inv*xstar)||(zstar`*r_inv*zstar + g_inv));
        rhs = (xstar`*r_inv*ystar)//(zstar`*r_inv*ystar);
        solution = inv(lhs)*rhs;

        beta_fixed_new = solution[1:2];
        beta_random_new = solution[3:nrow(solution)];
        beta_random_matrix = (shape(beta_random_new,2,n_effects))`;

        vm = beta_fixed_new[1];
        km = beta_fixed_new[2];

        vi = beta_random_matrix[,1];
        ki = beta_random_matrix[,2];

        vi_x = vi@j(n_sub,1,1);
        ki_x = ki@j(n_sub,1,1);

        fvm = x / (km + ki_x + x);
        fkm = -((vm + vi_x) # x) / ((km + ki_x + x) ## 2);
        xstar = fvm||fkm;

        fvi = design#fvm;
        fki = design#fkm;
        zstar = fvi||fki;

        var_fun = zstar*g_side*zstar`+r_side;
        var_inv = sweep(var_fun);

        p = var_inv-var_inv*xstar*sweep(xstar`*var_inv*xstar)*xstar`*var_inv;

        dv_vi = fvi*fvi`;
        dv_ki = fki*fki`;
        dv_e = i(nobs); 

        scv = -(1/2)#trace(p*dv_vi) + 
              (1/2)*((ystar-xstar*beta_fixed_new)`*var_inv*dv_vi*var_inv*(ystar-xstar*beta_fixed_new));
        sck = -(1/2)#trace(p*dv_ki) + 
              (1/2)*((ystar-xstar*beta_fixed_new)`*var_inv*dv_ki*var_inv*(ystar-xstar*beta_fixed_new));
        sce = -(1/2)#trace(p*dv_e) + 
              (1/2)*((ystar-xstar*beta_fixed_new)`*var_inv*dv_e*var_inv*(ystar-xstar*beta_fixed_new));
        score = scv//sck//sce;

        h11=0.5*trace(p*dv_vi*p*dv_vi);
        h12=0.5*trace(p*dv_vi*p*dv_ki);
        h13=0.5*trace(p*dv_vi*p*dv_e);
        h21=0.5*trace(p*dv_ki*p*dv_vi);
        h22=0.5*trace(p*dv_ki*p*dv_ki);
        h23=0.5*trace(p*dv_ki*p*dv_e);
        h31=0.5*trace(p*dv_e*p*dv_vi);
        h32=0.5*trace(p*dv_e*p*dv_ki);
        h33=0.5*trace(p*dv_e*p*dv_e);  
        h = (h11 || h12 || h13) // 
          (h21 || h22 || h23) // 
          (h31 || h32 || h33);

        old_sigma = sigma_random // sigma_residual;
        sigma = old_sigma + inv(h)*score;
        if sigma[1,1]<0.0001 then sigma[1,1]=0.0001;
        if sigma[2,1]<0.0001 then sigma[2,1]=0.0001;
        sigma_random = sigma[1:2];
        sigma_residual = sigma[3];

        sigma_vi = i(n_effects)*sigma_random[1];
        sigma_ki = i(n_effects)*sigma_random[2];

        g_side = block(sigma_vi, sigma_ki);
        g_inv = inv(g_side);

        r_side = i(nobs)*sigma_residual;
        r_inv = inv(r_side);

        var_fun = zstar*g_side*zstar`+r_side;
        var_inv = inv(var_fun);

        yhat = ((vm + vi_x) # x) / (km + ki_x + x);
        ystar = y - yhat + xstar*beta_fixed + zstar*beta_random;

        rss = ystar - xstar * inv(xstar` * var_inv * xstar) * xstar` * var_inv * ystar;
        test = det(var_fun);
        print test;
        new_log_PL = -0.5 * (log(det(var_fun)) + det(xstar` * var_inv * xstar) + rss` * var_inv * rss);

        crit = abs((new_log_PL - log_PL) / log_PL);
        log_PL = new_log_PL;
        beta_fixed = beta_fixed_new;
        beta_random = beta_random_new;
        niter = niter+1;
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


