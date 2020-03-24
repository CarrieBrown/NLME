proc import datafile="data.csv" 
out=simulated dbms=csv replace; 
getnames=yes; 
run;

proc iml;
    use simulated;
    read all;

    * Initial Values *;
    vm = 100;
    km = 10;

    vi_start = 0.1;
    ki_start = 0.1;

    sigma_residual = 15;
    sigma_random = {1, 0.5};

    * Begin Program *;
    start;
    nobs = nrow(y);
    design = design(id);
    n_effects = ncol(design);
    n_sub = nobs/n_effects;
    print nobs n_effects n_sub;
    crit = 1;
    niter = 0;

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

        *rss = ystar - xstar * inv(xstar` * var_inv * xstar) * xstar` * var_inv * ystar;
        *log_PL = -0.5 * (log(det(var_fun)) + det(xstar` * var_inv * xstar) + rss` * var_inv * rss);

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
        var_inv = inv(var_fun);

        p = var_inv-var_inv*xstar*inv(xstar`*var_inv*xstar)*xstar`*var_inv;

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
        if (sigma<0) then sigma[loc(sigma < 0)] = 0;
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

        *rss = ystar - xstar * inv(xstar` * var_inv * xstar) * xstar` * var_inv * ystar;
        *new_log_PL = -0.5 * (log(det(var_fun)) + det(xstar` * var_inv * xstar) + rss` * var_inv * rss);

        *crit = abs((new_log_PL - log_PL) / log_PL);
        *log_PL = new_log_PL;
        crit = max(abs((old_sigma - sigma) / old_sigma));
        beta_fixed = beta_fixed_new;
        beta_random = beta_random_new;
        
        niter = niter + 1;
        if niter > 200 then goto failed;        
    end;

    success:
    print "Converged after" niter "iternations - crit: " crit;
    
    c_11 = xstar`*r_inv*xstar;
    c_12 = xstar`*r_inv*zstar;
    c_21 = zstar`*r_inv*xstar;
    c_22 = zstar`*r_inv*zstar + g_inv;

    c_matrix = (c_11 || c_12) // (c_21 || c_22);
    c_inv = inv(c_matrix);
    
    se_fixed = sqrt(vecdiag(c_inv[1:2,1:2]));

    xi = {1,5,10,20,30,50};
    
    yhat_xi = ((vm) # xi) / (km + xi);
    dfvm = xi / (km + xi);
    dfkm = -((vm) # xi) / ((km + xi) ## 2);
    
    k = j(ncol(zstar),nrow(xi),0);
    k = dfvm || dfkm || k`;

    se_yhat_xi = sqrt(vecdiag(k*c_inv*k`));

    print xi yhat_xi se_yhat_xi;
    print sigma_random sigma_residual;

    var_vi = sigma_random[1];
    var_ki = sigma_random[2];

    iml_varest = vm || km || var_vi || var_ki || sigma_residual || se_fixed`;
    iml_varest_colnames = {"vm", "km", "var_vi", "var_ki", "var_res", "se_vm", "se_km"};

    create iml_varest from iml_varest [colname=iml_varest_colnames];
    append from iml_varest;
    close iml_varest;

    iml_varest2 = xi || yhat_xi || se_yhat_xi;
    iml_varest2_colnames = {"xi", "yhat_xi", "se_yhat_xi"};

    create iml_varest2 from iml_varest2 [colname=iml_varest2_colnames];
    append from iml_varest2;
    close iml_varest2;
    


/*    iml_pred = id || x || y || ystar || yhat;
    iml_pred_colnames = {"id", "x", "y", "ystar", "Pred"};
    create iml_pred from iml_pred [colname=iml_pred_colnames];
    append from iml_pred;
    close iml_pred; */
    
    goto results;

    failed:
    print "Failed to converge after" niter "iternations - crit: " crit;
    
    results:
    
finish;
run;
quit;


proc export data=iml_varest
outfile="iml_varest.csv"
DBMS=DLM REPLACE;
DELIMITER=",";
run;

proc export data=iml_varest2
outfile="iml_varest2.csv"
DBMS=DLM REPLACE;
DELIMITER=",";
run;

/* proc export data=iml_pred
outfile="iml_pred.csv"
DBMS=DLM REPLACE;
DELIMITER=",";
run; */
