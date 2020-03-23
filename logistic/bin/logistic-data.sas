data simulated;
n_effects = 10;
n_sub = 20;

a = 10;
b = 30;
c = -7;
d = .75;

var_ai = .5;
var_bi = 2.5;
var_ci = .05;
var_di = .005;
var_eu = .5;

do id=1 to n_effects;
  ai_x = sqrt(var_ai)*rand("Normal");
  bi_x = sqrt(var_bi)*rand("Normal");
  ci_x = sqrt(var_ci)*rand("Normal");
  di_x = sqrt(var_di)*rand("Normal");

  do x=1 to n_sub;
   res = sqrt(var_eu)*rand("Normal");

   y = a + ai_x + ((b + bi_x ) / (1 + exp(- (c + ci_x + ((d + di_x ) * x ))))) + res;
   true_y = a + (b / (1 + exp(- (c + (d * x )))));
   output;
  end;
 end;

drop n_effects n_sub a b c d ai_x bi_x ci_x di_x var_ai var_bi var_ci var_di var_eu res;

/* goptions reset=all;
symbol i=join;
proc gplot data=simulated;
 plot y*x=id;
 plot true_y*x;
run; */
proc export data=simulated
    outfile='data.csv'
    dbms=csv
    replace;
run;