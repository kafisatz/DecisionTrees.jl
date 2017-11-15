%macro save_julia_results(lib=rout);
	/*runmodel*/
	%if %sysfunc(exist(runmodel)) %then %do;
		data &lib..run_&jname.;set runmodel;run;
	%end; %else %do; %PUT WARNING: (bk) Data Set Runmodel does not exist!;%end;
	/*settings*/
	%if %sysfunc(exist(settings)) %then %do;
		data &lib..Set_&jname.;set settings;run;
	%end; %else %do; %PUT WARNING: (bk) Data Set Settings does not exist!;%end;
	/*res*/
	%if %sysfunc(exist(res)) %then %do;
	data &lib..res_&jname.;set res;run;
	%end; %else %do; %PUT WARNING: (bk) Data Set res does not exist!;%end;
	/*final est*/
	%if %sysfunc(exist(final_est)) %then %do;
		data &lib..err_&jname.;set final_est;ruN;
	%end; %else %do; %PUT WARNING: (bk) Data Set final_est does not exist!;%end;
	/*save result for boosting*/
	%if %sysfunc(exist(iterations_res)) %then %do;
		data &lib..eri_&jname.;set iterations_res;ruN;
	%end;
	/*itm*/	
	%if %sysfunc(exist(itm)) %then %do;
		data &lib..itm_&jname.;set itm;ruN;
	%end;
	/*save Result for Rankoptimization*/
	%if %symexist(bROSASProduceRkOptStats) %then %do;
		%if &bROSASProduceRkOptStats. eq 1 %then %do;
			%if %sysfunc(exist(rankopt)) %then %do;
				data &lib..rko_&jname.;set rankopt;ruN;
			%end;
			%else %do; %PUT WARNING: (bk) Data Set rankopt does not exist!;%end;
		%end;
	%end;

/*copy CSV files*/options noxwait;
%let statsfile=&outfname_of_julia_macro..stats.csv;
	data _null_; length xx $1000.; xx=cat("copy ",'"',"&datafolder_of_julia_macro.\&statsfile.",'"'," ",'"',"%sysfunc(pathname(rout))\&statsfile.",'"'); 
	fname="temp";rc=filename(fname,"&datafolder_of_julia_macro.\&statsfile.");if rc = 0 and fexist(fname) then rc=system(xx);run;
/*copy XLSX file */
%let statsfile=&outfname_of_julia_macro..xlsx;
	data _null_; length xx $1000.; xx=cat("copy ",'"',"&datafolder_of_julia_macro.\&statsfile.",'"'," ",'"',"%sysfunc(pathname(rout))\&statsfile.",'"'); 
	fname="temp";rc=filename(fname,"&datafolder_of_julia_macro.\&statsfile.");if rc = 0 and fexist(fname) then rc=system(xx);run;
/*%let rkoptstatsinfile=&outfname_of_julia_macro..roptstats.csv;*/
/*	data _null_; length xx $1000.; xx=cat("copy ",'"',"&datafolder_of_julia_macro.\&rkoptstatsinfile.",'"'," ",'"',"%sysfunc(pathname(rout))\&rkoptstatsinfile.",'"'); */
/*	fname="temp";rc=filename(fname,"&datafolder_of_julia_macro.\&rkoptstatsinfile.");if rc = 0 and fexist(fname) then rc=system(xx);run;*/
%mend;
