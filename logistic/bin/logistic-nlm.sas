proc import datafile="data.csv" 
out=simulated dbms=csv replace; 
getnames=yes; 
run;

ods results off;

proc nlmixed data=simulated;
 parms alpha=10 beta=30 gamma=-5 delta=0.5 logva=-1 to 1 logvb=0 to 3 logvc=-2 to 0 logvd=-5 to -3 logeuv=0;
 eta = alpha + ai + (beta + bi)/(1 + exp(-((gamma + ci) + (delta + di) * x)));
 model y~normal(eta,exp(logeuv));
 random ai bi ci di ~ normal([0,0,0,0],[exp(logva),0, exp(logvb),0,0,exp(logvc),0,0,0,exp(logvd)]) subject=id;
 
 estimate 'y-hat at x=1' alpha + (beta/(1 + exp(-(gamma + delta))));
 estimate 'y-hat at x=5' alpha + (beta/(1 + exp(-(gamma + 5*delta))));
 estimate 'y-hat at x=10' alpha + (beta/(1 + exp(-(gamma + 10*delta))));
 estimate 'y-hat at x=15' alpha + (beta/(1 + exp(-(gamma + 15*delta))));
 estimate 'y-hat at x=20' alpha + (beta/(1 + exp(-(gamma + 20*delta))));

 estimate 'ai variance' exp(logva);
 estimate 'bi variance' exp(logvb);
 estimate 'ci variance' exp(logvc);
 estimate 'di variance' exp(logvd);
 estimate 'eu variance' exp(logeuv);

 *predict eta out=nlm_pred;
 ods output ParameterEstimates=nlm_varest;
 ods output AdditionalEstimates=nlm_varest2;
run;

proc export data=nlm_varest
outfile="nlm_varest.csv"
DBMS=DLM REPLACE;
DELIMITER=",";
run;

proc export data=nlm_varest2
outfile="nlm_varest2.csv"
DBMS=DLM REPLACE;
DELIMITER=",";
run;

/*proc export data=nlm_pred
outfile="nlm_pred.csv"zs
DBMS=DLM REPLACE;
DELIMITER=",";
run;*/