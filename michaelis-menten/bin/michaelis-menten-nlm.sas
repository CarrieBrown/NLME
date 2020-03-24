proc import datafile="data.csv" 
out=simulated dbms=csv replace; 
getnames=yes; 
run;

ods results off;

proc nlmixed data=simulated;
  parms vmax=100 km=10 logvv=-2.5 logvk=-.5 logveu=2.5;
  eta=((vmax + vi) * x)/(km + ki + x);
  model y ~ normal(eta, exp(logveu));
  random vi ki ~ normal([0,0],[exp(logvv),0,exp(logvk)]) subject=id;

  estimate 'y-hat at x=1' ((vmax) / (km + 1));
  estimate 'y-hat at x=5' ((vmax * 5) / (km + 5));
  estimate 'y-hat at x=10' ((vmax * 10) / (km + 10));
  estimate 'y-hat at x=20' ((vmax * 20) / (km + 20));
  estimate 'y-hat at x=30' ((vmax * 30) / (km + 30));
  estimate 'y-hat at x=50' ((vmax * 60) / (km + 60));

  estimate 'vmax variance' exp(logvv);
  estimate 'km variance' exp(logvk);
  estimate 'eu variance' exp(logveu);

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

/* proc export data=nlm_pred
outfile="nlm_pred.csv"
DBMS=DLM REPLACE;
DELIMITER=",";
run; */