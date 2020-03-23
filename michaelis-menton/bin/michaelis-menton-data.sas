data simulated;
n_effects = 30;
x_min = 0;
x_max = 50;
x_int = 5;

vm = 100;
km = 10;

var_vi = 5;
var_ki = 0.5;
var_eu = 15;

do id = 1 to n_effects;
 vi_x = sqrt(var_vi)*rand("Normal");
 ki_x = sqrt(var_ki)*rand("Normal");
 
 do x = x_min to x_max by x_int;
  res = sqrt(var_eu)*rand("Normal");
  
   y = ((vm + vi_x) * x) / (km + ki_x + x) + res;
   true_y = (vm * x) / (km + x);
   output;
 end;
end;

drop n_effects x_min x_max x_int vm km vi ki var_vi var_ki var_eu;

proc export data=simulated
    outfile='data.csv'
    dbms=csv
    replace;
run;
