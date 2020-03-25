proc import datafile="data.csv" 
out=simulated dbms=csv replace; 
getnames=yes; 
run;

proc iml;
    use simulated;
    read all;

    * Initial Values;
    a = 9.6;
    b = 30;
    c = -5;
    d = 0.5;

    ai_start = 0.01;
    bi_start = 0.01;
    ci_start = 0.01;
    di_start = 0.01;

    sigma_residual = 1;
    sigma_random = {.75,3,0.2,0.01};

    * Begin Program;
    start;
    nobs = nrow(y);
    design = design(id);
    n_effects = ncol(design);
    n_sub = nobs/n_effects;
    crit = 1;
    niter = 1;

    ai = j(n_effects,1,ai_start);
    bi = j(n_effects,1,bi_start);
    ci = j(n_effects,1,ci_start);
    di = j(n_effects,1,di_start); 

    beta_fixed = a//b//c//d;
    beta_random = ai//bi//ci//di;

    ai_x = ai@j(n_sub,1,1);
    bi_x = bi@j(n_sub,1,1);
    ci_x = ci@j(n_sub,1,1);
    di_x = di@j(n_sub,1,1);

    fa = j(nobs,1,1);
    fb = 1 / (1 + EXP(- (c + ci_x + (d + di_x) # x)));
    fc = - (- EXP(- (c + ci_x + (d + di_x) # x)) # (b + bi_x) / (1 + EXP(- (c + ci_x + (d + di_x) # x)))##2);
    fd = - (- x # EXP(- (c + ci_x + (d + di_x) # x)) # (b + bi_x) / (1 + EXP(- (c + ci_x + (d + di_x) # x)))##2);
    xstar = fa||fb||fc||fd;

    fai = design;
    fbi = design#fb;
    fci = design#fc;
    fdi = design#fd;
    zstar = fai||fbi||fci||fdi;

    sigma_ai = i(n_effects)#sigma_random[1];
    sigma_bi = i(n_effects)#sigma_random[2];
    sigma_ci = i(n_effects)#sigma_random[3];
    sigma_di = i(n_effects)#sigma_random[4];

    g_side = block(sigma_ai, sigma_bi, sigma_ci, sigma_di);
    g_inv = inv(g_side);

    r_side = i(nobs)#sigma_residual;
    r_inv = inv(r_side);
    
    var_fun = zstar*g_side*zstar`+r_side;
    var_inv = inv(var_fun);

    do while (crit>1e-8);

        yhat = a + ai_x + ((b + bi_x) / (1 + exp(- (c + ci_x + (d + di_x) # x))));
        ystar = y - yhat + xstar*beta_fixed + zstar*beta_random;

        * rss = ystar - xstar * inv(xstar` * var_inv * xstar) * xstar` * var_inv * ystar;
        * log_PL = -0.5 * log(det(var_fun)) + det(xstar` * var_inv * xstar) + rss` * var_inv * rss;

        lhs = ((xstar`*r_inv*xstar)||(xstar`*r_inv*zstar)) //
         ((zstar`*r_inv*xstar)||(zstar`*r_inv*zstar+g_inv));
        rhs = (xstar`*r_inv*ystar)//(zstar`*r_inv*ystar);
        solution = inv(lhs)*rhs;
        beta_fixed_new = solution[1:4];
        beta_random_new = solution[5:nrow(solution)];
        beta_random_matrix = (shape(beta_random_new,4,n_effects))`;

        a = beta_fixed_new[1];
        b = beta_fixed_new[2];
        c = beta_fixed_new[3];
        d = beta_fixed_new[4];

        ai = beta_random_matrix[,1];
        bi = beta_random_matrix[,2];
        ci = beta_random_matrix[,3];
        di = beta_random_matrix[,4];

        ai_x = ai@j(n_sub,1,1);
        bi_x = bi@j(n_sub,1,1);
        ci_x = ci@j(n_sub,1,1);
        di_x = di@j(n_sub,1,1);

        yhat = a + ai_x + ((b + bi_x) / (1 + exp(- (c + ci_x + (d + di_x) # x))));

        fa = j(nobs,1,1);
        fb = 1 / (1 + EXP(- (c+ci_x + (d + di_x) # x)));
        fc = - (- EXP(- (c + ci_x + (d + di_x) # x)) # (b + bi_x) / (1 + EXP(- (c + ci_x + (d + di_x) # x)))##2);
        fd = - (- x # EXP(- (c+ci_x + (d+di_x) # x)) # (b+bi_x) / (1 + EXP(- (c+ci_x + (d+di_x) # x)))##2);
        xstar = fa||fb||fc||fd;

        fai = design;
        fbi = design#fb;
        fci = design#fc;
        fdi = design#fd;
        zstar = fai||fbi||fci||fdi;

        var_fun = zstar*g_side*zstar`+r_side;
        var_inv = inv(var_fun);

        p = var_inv-var_inv*xstar*inv(xstar`*var_inv*xstar)*xstar`*var_inv;

        dv_ai = fai*fai`;
        dv_bi = fbi*fbi`;
        dv_ci = fci*fci`;
        dv_di = fdi*fdi`;
        dv_e = i(nobs); 

        sca = (-0.5)#trace(p*dv_ai) + 
        (1/2)#((ystar-xstar*beta_fixed_new)`*var_inv*dv_ai*var_inv*(ystar-xstar*beta_fixed_new));
        scb = (-0.5)#trace(p*dv_bi) + 
        (1/2)#((ystar-xstar*beta_fixed_new)`*var_inv*dv_bi*var_inv*(ystar-xstar*beta_fixed_new));
        scc = (-0.5)#trace(p*dv_ci) + 
        (1/2)#((ystar-xstar*beta_fixed_new)`*var_inv*dv_ci*var_inv*(ystar-xstar*beta_fixed_new));
        scd = (-0.5)#trace(p*dv_di) + 
        (1/2)#((ystar-xstar*beta_fixed_new)`*var_inv*dv_di*var_inv*(ystar-xstar*beta_fixed_new));
        sce = (-0.5)#trace(p*dv_e) + 
        (1/2)#((ystar-xstar*beta_fixed_new)`*var_inv*dv_e*var_inv*(ystar-xstar*beta_fixed_new));
        score = sca//scb//scc//scd//sce;

        h11=0.5#trace(p*dv_ai*p*dv_ai);
        h12=0.5#trace(p*dv_ai*p*dv_bi);
        h13=0.5#trace(p*dv_ai*p*dv_ci);
        h14=0.5#trace(p*dv_ai*p*dv_di);
        h15=0.5#trace(p*dv_ai*p*dv_e);
        h21=0.5#trace(p*dv_bi*p*dv_ai);
        h22=0.5#trace(p*dv_bi*p*dv_bi);
        h23=0.5#trace(p*dv_bi*p*dv_ci);
        h24=0.5#trace(p*dv_bi*p*dv_di);
        h25=0.5#trace(p*dv_bi*p*dv_e);
        h31=0.5#trace(p*dv_ci*p*dv_ai);
        h32=0.5#trace(p*dv_ci*p*dv_bi);
        h33=0.5#trace(p*dv_ci*p*dv_ci);
        h34=0.5#trace(p*dv_ci*p*dv_di);
        h35=0.5#trace(p*dv_ci*p*dv_e);
        h41=0.5#trace(p*dv_di*p*dv_ai);
        h42=0.5#trace(p*dv_di*p*dv_bi);
        h43=0.5#trace(p*dv_di*p*dv_ci);
        h44=0.5#trace(p*dv_di*p*dv_di);
        h45=0.5#trace(p*dv_di*p*dv_e);
        h51=0.5#trace(p*dv_e*p*dv_ai);
        h52=0.5#trace(p*dv_e*p*dv_bi);
        h53=0.5#trace(p*dv_e*p*dv_ci);
        h54=0.5#trace(p*dv_e*p*dv_di);
        h55=0.5#trace(p*dv_e*p*dv_e);  
        h = (h11 || h12 || h13 || h14 || h15) // 
        (h21 || h22 || h23 || h24 || h25) // 
        (h31 || h32 || h33 || h34 || h35) // 
        (h41 || h42 || h43 || h44 || h45) // 
        (h51 || h52 || h53 || h54 || h55);

        old_sigma = sigma_random // sigma_residual;
        sigma = old_sigma + inv(h)*score;
        if (sigma<0) then sigma[loc(sigma < 0)] = 0;
        sigma_random = sigma[1:4];
        sigma_residual = sigma[5];

        sigma_ai = i(n_effects)#sigma_random[1];
        sigma_bi = i(n_effects)#sigma_random[2];
        sigma_ci = i(n_effects)#sigma_random[3];
        sigma_di = i(n_effects)#sigma_random[4];

        g_side = block(sigma_ai, sigma_bi, sigma_ci, sigma_di);
        g_inv = inv(g_side);

        r_side = i(nobs)#sigma_residual;
        r_inv = inv(r_side);

        var_fun = zstar*g_side*zstar`+r_side;
        var_inv = inv(var_fun);

        yhat = a + ai_x + ((b + bi_x) / (1 + exp(- (c + ci_x + (d + di_x) # x))));
        ystar = y - yhat + xstar*beta_fixed_new + zstar*beta_random_new;

        * rss = ystar - xstar * inv(xstar` * var_inv * xstar) * xstar` * var_inv * ystar;
        * new_log_PL = -0.5 * log(det(var_fun)) + det(xstar` * var_inv * xstar) + rss` * var_inv * rss;

        * crit = abs((new_log_PL - log_PL) / log_PL);
        * log_PL = new_log_PL;
        crit=max(abs((old_sigma-sigma)/old_sigma));
        beta_fixed = beta_fixed_new;
        beta_random = beta_random_new;
        niter = niter+1;
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

    se_fixed = sqrt(vecdiag(c_inv[1:4,1:4]));
    var_random = vecdiag(c_inv[5:nrow(c_inv),5:ncol(c_inv)]);

    xi = {1,5,10,15,20};

    yhat_xi = a + (b / (1 + exp(- c - (d # xi))));
    dfa = j(nrow(xi),1,1); 
    dfb = 1 / (1 + EXP(- (c + (d # xi))));
    dfc = - (- EXP(- (c + (d # xi))) * b / (1 + EXP(- (c + (d # xi))))##2);
    dfd = - (- xi # EXP(- (c + (d # xi))) * b / (1 + EXP(- (c + (d # xi))))##2);

    k = j(ncol(zstar),nrow(xi),0);
    k = dfa || dfb || dfc || dfd || k`;

    se_yhat_xi = sqrt(vecdiag(k*c_inv*k`));

    var_ai = sigma_random[1];
    var_bi = sigma_random[2];
    var_ci = sigma_random[3];
    var_di = sigma_random[4];

    iml_varest = a || b || c || d || 
    var_ai || var_bi || var_ci || var_di || sigma_residual || se_fixed`;

    iml_varest_colnames = {"a", "b", "c", "d", 
    "var_ai", "var_bi", "var_ci", "var_di", "var_res",
    "se_a", "se_b", "se_c", "se_d"};

    create iml_varest from iml_varest [colname=iml_varest_colnames];
    append from iml_varest;
    close iml_varest;

    iml_varest2 = xi || yhat_xi || se_yhat_xi;
    iml_varest2_colnames = {"xi", "yhat_xi", "se_yhat_xi"};

    create iml_varest2 from iml_varest2 [colname=iml_varest2_colnames];
    append from iml_varest2;
    close iml_varest2;

    print xi yhat_xi se_yhat_xi;
    print sigma_random sigma_residual;

    /*iml_pred = id || x || y || ystar || yhat;
    iml_pred_colnames = {"id", "x", "y", "ystar", "Pred"};

    create iml_pred from iml_pred [colname=iml_pred_colnames];
    append from iml_pred;
    close iml_pred;*/
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

/*proc export data=iml_pred
outfile="iml_pred.csv"
DBMS=DLM REPLACE;
DELIMITER=",";
run; */


