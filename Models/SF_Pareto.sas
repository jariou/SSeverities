/*--------------------------------------------------------------
        Name: SF_Pareto.sas
       Title: SF_PARETO distribution
 Description: Model definition for SF_PARETO distribution
     Product: SAS/ETS Software
        Keys: fitting continuous distributions
        PROC: SEVERITY
--------------------------------------------------------------*/
%macro install_sf_pareto(
                         inLib      , 
                         outLib     , 
                         outDataset , 
                         outPackage
                         );
   %let libStatement = %str();
   
   %if &inLib NE %str() %then %do;
      %let libStatement = unquote("library = &inlib.");
   %end;
                         
   proc fcmp 
            &libStatement.
            outlib  = &library..&dataset..&package
            ;

      function SF_PARETO_DESCRIPTION() $32;
         length model $32;
         model = "StateFarm Pareto Distribution";
         return(model);
      endsub;

      function SF_PARETO_PARMCOUNT();
         return(2);
      endsub;

      function SF_PARETO_PDF(x, Theta, Alpha);
         /* Theta : Scale */
         /* Alpha : Shape */
         if (x < Theta) then do;
            return (0);
         end;
         else do;
            logpdf = log(Alpha / x) + Alpha * log(Theta / x);

            if (logpdf < 174) then do;
               return ( exp(logpdf) );
            end;
            else do;
               return ( . );
            end;
         end;
      endsub;

      function SF_PARETO_CDF(x, Theta, Alpha);
         /* Theta : Scale */
         /* Alpha : Shape */
         if (x < Theta) then do;
            return (0);
         end;
         else do;
            return (1 - (Theta / x)**Alpha);
         end;
      endsub;

      subroutine SF_PARETO_PARMINIT(dim, x[*], nx[*], F[*], FType, Theta, Alpha);
         outargs Theta, Alpha;
         array m[2] / nosymbols;

         put Version  = "New";

         Theta = 10000000000;

         do i = 1 to dim;
            if not missing(x[i]) then do;
               Theta = min(x[i], Theta);
            end;
         end;

         meanLog = 0;
         countLog = 0;

         do i = 1 to dim;
            if not missing(x[i]) then do;
               countLog = countLog + 1;      
               meanLog = meanLog + log(x[i] / Theta);
            end;
         end;

         Alpha = countLog / meanLog;
      endsub;

      subroutine SF_PARETO_PARMINIT_OLD(dim, x[*], nx[*], F[*], FType, Theta, Alpha);
         outargs Theta, Alpha;
         array m[2] / nosymbols;

         /* Use Method of Moments */
         call svrtutil_rawmoments(dim, x, nx, 2, m);
         if (missing(m[1])) then do;
            Theta = .;
            Alpha = .;
         end;
         else do;
            t1 = m[2] - m[1]**2;
            t2 = 2*t1 - m[2];
            eps = constant('MACEPS');

            if (t1 < eps or t2 < eps) then do;
               Theta = m[1];
               Alpha = 2;
            end;
            else do;
               Theta = (m[1] * m[2]) / t2;
               Alpha = (2 * t1) / t2;
            end;
         end;
      endsub;

      function SF_PARETO_QUANTILE(cdf, Theta, Alpha);
         return(Theta * (1 - cdf)**(-1 / Alpha));
      endsub;

      function SF_PARETO_ROll(seed, Theta, Alpha);
         p = ranuni(seed);
         return(SF_PARETO_QUANTILE(p, Theta, Alpha));
      endsub;

   run;
   quit;
%mend install_sf_pareto;

