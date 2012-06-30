
proc fcmp 
            library = sashelp.svrtdist 
            outlib  = sasuser.StateFarm.SAS_Pareto;

   function PARETO_QUANTILE(cdf, Theta, Alpha);
      return(Theta * ((1- cdf)**(-1 / Alpha) - 1));
   endsub;

   function PARETO_ROll(seed, Theta, Alpha);
      p = ranuni(seed);
      return(PARETO_QUANTILE(p, Theta, Alpha));
   endsub;

run;
quit;

option append = (cmplib = sasuser.statefarm);
proc severity 
               data   = autoclaims
               outest = Pareto
               ;
    loss loss;
    dist Pareto;
    dist sf_pareto;
run;


proc fcmp 
            library = sashelp.svrtdist 
            outlib  = sasuser.StateFarm.SAS_Pareto;

   function PARETO_QUANTILE(cdf, Theta, Alpha);
      return(Theta * ((1- cdf)**(-1 / Alpha) - 1));
   endsub;

   function PARETO_ROll(seed, Theta, Alpha);
      p = ranuni(seed);
      return(PARETO_QUANTILE(p, Theta, Alpha));
   endsub;

run;
quit;

option append = (cmplib = sasuser.statefarm);
proc severity 
               data   = autoclaims
               outest = Pareto
               ;
    loss loss;
    dist Pareto;
    dist sf_pareto;
run;
