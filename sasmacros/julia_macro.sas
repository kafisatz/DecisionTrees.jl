%macro julia_macro(juliadatafolder=default,version=prod,indatafile=julia.csv,insettingsfile=settings.julia.csv,outfilename=out_julia,jmnotes=0,put_fullstring=1);
/*This runs Julia (run.jl)*/

/*
%let juliadatafolder=r:\temp\;
%let outfilename=out_&jname.;
%let matrix=&juliadatafolder.\&outfilename._iteration_matrix.csv;
%let matrix=&juliadatafolder.out_julia_iteration_matrix.csv;
%let outcsv=out_julia.csv;
%let out=res;
*/

%global algorithm_type outfname_of_julia_macro datafolder_of_julia_macro;
%local stop_julia_macro var_order len end out mdfvecsize this_mdf in this_iter i moderationvector juliaexecutable juliabatch juliarootfolder juliaprogfolder indatafile outcsv;
%let stop_julia_macro=0;

%if %length(&jname.)>28 %then %do;
	%let stop_julia_macro=1;
	%put ERROR: (bk) Length of Jname must be less than 28 (as SAS dataset names are limited to 32 characters);
%end;
%else %do;

%if not %SYMEXIST(nthreads) %then %do; %global nthreads;%let nthreads=0;%end;

%if not %SYMEXIST(smoothEstimates) %then %do; %global smoothEstimates;%let smoothEstimates=1;%end;
%if &model_type. = build_tree %then %let smoothEstimates=0;
%if not %SYMEXIST(indata) %then %do; %global indata;%let indata=undefined;%end;
%if not %length(&indata.) %then %let indata=unused;
%if not %SYMEXIST(var_dep) %then %do;%global var_dep; %let var_dep=undefined;%end;
%if not %SYMEXIST(initial_sample_size) %then %do; %global initial_sample_size;%let initial_sample_size=max;%end;
%if not %SYMEXIST(validation_sample_size) %then %do; %global validation_sample_size;%let validation_sample_size=max;%end;
%if not %SYMEXIST(write_iteration_matrix) %then %do; %global write_iteration_matrix;%let write_iteration_matrix=1;%end;
%if not %SYMEXIST(preppedJLDFileExists) %then %do; %global preppedJLDFileExists;%let preppedJLDFileExists=0;%end;

/*todo/tbd improve these settings.... sometimes they do not work as intended*/
goptions DEVICE=SVG;
ods html close;
ods html5 path='c:\temp\' options(svg_mode="inline");
ods graphics /outputfmt=svg;

goptions DEVICE=SVG;

options NOQUOTELENMAX;

data vtrn/view=vtrn;set market(obs=&initial_sample_size. where=(training_validation_bool eq 1));*if _N_ eq 1 then put "WARNING: (bk) Make sure dates are formatted as integer values!";run;
data vval/view=vval;set market(obs=&validation_sample_size. where=(training_validation_bool eq 0));run;/*%distinct_values(data=in,out=levels);*/

%if not %SYMEXIST(standard_plot) %then %do; %global standard_plot;%let standard_plot=1;%end;

%timestamp;


%let juliaprogfolder=&algorithmsFolder.\Julia\Code\&version.\src;
/*Set default Julia Settings for the parameters which were not provided*/

%if not %symexist(boolRankOptimization) %then %do;%global boolRankOptimization;%let boolRankOptimization=0;%end;
%if not %symexist(bROSASProduceRkOptStats) %then %do;%global bROSASProduceRkOptStats;%let bROSASProduceRkOptStats=0;%end;
%if &boolRankOptimization. %then %let bROSASProduceRkOptStats=1;

%if not %length(&juliadatafolder.) %then %let juliadatafolder=&global_sas_mvar_tmpfolder.;
%if &juliadatafolder.=default %then %let juliadatafolder=&global_sas_mvar_tmpfolder.;
%let infileWithPath=&juliadatafolder.\&indatafile.;
%let settingsfileWithPath=&juliadatafolder.\&insettingsfile.;
/*%put &=settingsfileWithPath;*/
*remove \ at the end;%if %substr(%sysfunc(reverse(&juliadatafolder.)),1,1)=\ %then %let juliadatafolder=%substr(&juliadatafolder.,1,%length(&juliadatafolder.)-1);
%let out=res;
%let outcsv=&outfilename..csv;
%let juliaexecutable=&global_sas_mvar_juliafolder.julia.exe;

*determine whether a moderationvector was supplied;
/*%if %sysfunc(countw(&mf.," "))>1 %then %do;*/
/*	%let moderationvector=&mf.;*/
/*	%let mf=0.0;*/
/*%end;*/
%let mdfvecsize=%sysfunc(countw(&moderationvector.," "));
%let model_type=%sysfunc(lowcase(&model_type.));
%let algorithm_type=j_&model_type.;
%let matrix=&juliadatafolder.\&outfilename._iteration_matrix.csv;
%let statsfile=&juliadatafolder.\&outfilename..stats.csv;

%if not &preppedJLDFileExists. %then %do;
	proc datasets nolist;delete runmodel ;quit;
%end;

proc datasets nolist;delete julia_leafnrs itm rankopt final_est &out.;quit;
proc datasets nolist memtype=view;delete itm rankopt final_est &out.;quit;

*I am unsure whether the sorting by key is really necessary, but I would keep it as best practice;
/*proc sort data=vtrn;by irkey;quit;*/
/*proc sort data=val;by irkey;quit;*/

/*Numeric independent columns need to come first, then charater*/
proc contents noprint data=vtrn out=tmp_cont(keep=NAME type);quit;
data tmp_cont2;set tmp_cont;if lowcase(name) in ("numerator" "denominator" "irkey" "weight" "training_validation_bool") then delete;run;
%let char_vars=;%let num_vars=;
proc sql noprint; 
	select upcase(name) into: num_vars separated by " " from tmp_cont2 where type ne 2;
	select upcase(name) into: char_vars separated by " " from tmp_cont2 where type eq 2;
quit;
%let number_of_num_features=0;%if %length(&num_vars.) %then %let number_of_num_features=%sysfunc(countw(&num_vars.));

%let numberOfRows=-1;
data _null_;
	set vtrn(keep=irkey) vval(keep=irkey);
	call symputx("numberOfRows",_N_);
run;

/*Create Settings data set*/
/*NOTE: this macro variable needs to be consistent with the field names of the type ModelSettings*/
%let possible_settings_columns=model_type minw randomw crit max_splitting_points_num niter mf nscores adaptiveLearningRate number_of_tariffs prem_buffer BoolStartAtMean bool_write_tree number_of_num_features parallel_tree_construction using_local_variables parallel_level_threshold parallel_weight_threshold nthreads dataIdentifier algorithmsFolder startime indata var_dep indepcount spawnsmaller recursivespawning pminweightfactor pminrelpctsize pflipspawnsmalllargedepth juliaprogfolder boolMDFPerLeaf ranks_lost_pct variable_mdf_pct boolRankOptimization bROSASProduceRkOptStats boolRandomizeOnlySplitAtTopNode subsampling_prop subsampling_features_prop version preppedJLDFileExists catSortByThreshold catSortBy scorebandsstartingpoints showTimeUsedByEachIteration smoothEstimates     roptForcedPremIncr     premStep write_sas_code write_iteration_matrix write_result write_statistics boolCreateZipFile write_csharp_code boolTariffEstStats bINTERNALignoreNegRelsBoosting statsbyvariables statsRandomByVariable boolSaveResultAsJLDFile print_details;
data settings;
%do iii=1 %to %sysfunc(countw(&possible_settings_columns.));
	%let this_iii=%scan(&possible_settings_columns.,&iii.);
	%if %SYMEXIST(&this_iii.) %then &this_iii.="&&&this_iii";;
%end;
run;

%if not &preppedJLDFileExists. %then %do; /*runmodel already exists*/
	%let var_order=%nrstr(format irkey;format numerator denominator weight best. training_validation_bool;format )&num_vars. &char_vars.;
	data runmodel;
		&var_order.;
		set vtrn vval;
		if training_validation_bool ne 1 then training_validation_bool=0;
	run;
	%let len=%length(&infileWithPath.);
	%let end=%upcase(%sysfunc(lowcase(%substr(&infileWithPath.,&len.-2,3))));
	%if &end. eq SQL %then %do;
		/*
		The SQL DB will have varchar(x) columns where x is the length of the sas
		variable. As this length can be far to long, we reduce the length of the variables
		in the SAS dataset. This should yield smaller SQL tables
		*/
		%minimize_length_of_char_vars(runmodel);
		/*Unfortunately the order of the variables is now changed....*/
		data runmodel;
			&var_order.;
			set runmodel;
		ruN;
	%end;
%end;


%if %vartype(data=runmodel,var=irkey) ne C %then %do;
	%let stop_julia_macro=1;
	%put ERROR: (bk) The Data set runmodel needs a variable irkey which must be of type character. Abort;
%end;
%if not %sysfunc(exist(runmodel))  %then %do;
	%let stop_julia_macro=1;
	%put ERROR: (bk) Data set runmodel does not exist. Abort;
%end;
%else %do;
%if not %getattrn(data=runmodel) %then %do;
	%let stop_julia_macro=1;
	%put ERROR: (bk) Data set runmodel is empty. Abort;
%end;
%else %do;

%if &stop_julia_macro. %then %do;
	%put ERROR: (bk) &=stop_julia_macro. See above. Abort;
%end;
%else %do;
options xwait;
data _null_;
/*Delete the input and outputfile*/
    fname="temp";rc=filename(fname,"&matrix.");if rc = 0 and fexist(fname) then rc=fdelete(fname);rc=filename(fname);
    fname="temp";rc=filename(fname,"&juliadatafolder.\&outcsv.");if rc = 0 and fexist(fname) then rc=fdelete(fname);rc=filename(fname);
	fname="temp";rc=filename(fname,"&juliadatafolder.\&outfilename..sas");if rc = 0 and fexist(fname) then rc=fdelete(fname);rc=filename(fname);
	fname="temp";rc=filename(fname,"&juliadatafolder.\&outfilename..txt");if rc = 0 and fexist(fname) then rc=fdelete(fname);rc=filename(fname);	
	fname="temp";rc=filename(fname,"&juliadatafolder.\&outfilename..csv");if rc = 0 and fexist(fname) then rc=fdelete(fname);rc=filename(fname);	
	fname="temp";rc=filename(fname,"&juliadatafolder.\&outfilename..stats.csv");if rc = 0 and fexist(fname) then rc=fdelete(fname);rc=filename(fname);	
	fname="temp";rc=filename(fname,"&juliadatafolder.\&outfilename.roptstats.csv");if rc = 0 and fexist(fname) then rc=fdelete(fname);rc=filename(fname);	
	fname="temp";rc=filename(fname,"&juliadatafolder.\&outfilename._iteration_matrix.csv");if rc = 0 and fexist(fname) then rc=fdelete(fname);rc=filename(fname);	
	fname="temp";rc=filename(fname,"&settingsfileWithPath.");if rc = 0 and fexist(fname) then rc=fdelete(fname);rc=filename(fname);	
	if lowcase(substr(reverse("&infileWithPath."),1,4))="vsc." then do;
		fname="temp";rc=filename(fname,"&infileWithPath.");if rc = 0 and fexist(fname) then rc=fdelete(fname);rc=filename(fname);	
	end;
run;

/*Export Data*/
/*export settings column*/		
/*%delete_windows_file(c:\temp\testme.xlsb)*/
	%delete_windows_file(&settingsfileWithPath.);
	%export_csv_comma_utf8(data=settings,file=&settingsfileWithPath.);
%if &preppedJLDFileExists. eq 0 %then %do;
	%put NOTE: (bk) Exporting...;
		%if &jmnotes. %then %do;option nonotes nosource;		%end;
			%delete_windows_file(&infileWithPath.);
/*			%export_csv_comma_utf8(data=runmodel,file=&infileWithPath.);*/
			%export_for_julia(data=runmodel,file=&infileWithPath.);			
		%if &jmnotes. %then %do;option notes source;%end;
%end;
%else %do;
	%put NOTE: (bk) Exporting only settings csv file: &juliadatafolder.\&insettingsfile.;
%end;

options xmin; *x command will start minimized;

%if &nthreads.>0 %then %do;%let _addstring_=-p &nthreads. -F;%end;%else %do; %let _addstring_=-F;%end;
/*Run Julia*/
data _null_;
	*the first command (before &&) sets the buffer (heigth/width) of the window;
	fullstring='mode con:cols=116 lines=5000 &&'||quote("&juliaexecutable.")||" &_addstring_. "||quote("&juliaprogfolder.\run.jl")
		||" "||quote("&settingsfileWithPath.")||" "||quote("&infileWithPath.")||" "||quote("&outfilename.")||"%NRSTR( && IF ERRORLEVEL 0 EXIT)";	*the last term is to keep the command window open after an error ocurred;
	*Replace single backslash by double backslash;
	fullstring=TRANSTRN(fullstring,"\\","\");
	fullstring=TRANSTRN(fullstring,"\\","\");
	fullstring=TRANSTRN(fullstring,"\\\","\");
	%if &put_fullstring. %then %do;put fullstring=;%end;
	/*put fullstring=;*/
	rc=system(fullstring);
	/*	if rc eq 0 then rc2=system('exit');*/	
	/*Old CODE: a=quote("&juliabatch.");call system("cd "||quote("&juliaprogfolder."));call system(a);*/
run;

data _null_;put "NOTE: (bk) Command Line Window has closed has finished.";run;

*Timestamp;
%timestamp;%let time_julia_m=&elapsed_m_.;
*Check if Julia produced any output?;

proc datasets nolist;delete itm iterations_res &out. inn inn2 final_est out tm: &out.: &out inn2 out &out. res: final_est resiterations tmp: tmp_mat: tmpiter;quit;

%let leafnrfile=&juliadatafolder.\&outfilename..leafnumbers.csv;
%if &write_iteration_matrix. %then %do;
	%if %FileExist(&leafnrfile.) %then %do;	 		
		PROC IMPORT OUT=julia_leafnrs DATAFILE= "&leafnrfile." DBMS=DLM REPLACE;DELIMITER='2C'x; guessingrows=10001;GETNAMES=NO;DATAROW=1;RUN;
		%put WARNING: (bk) there is an issue here if the first n rows of the leafnr file have a short irkey, then the length of that variable will not be long enough and thus cut off.;
		proc datasets nolist; modify julia_leafnrs;rename var1=irkey
		%do i=1 %to &niter.;
			var%eval(&i.+1)=leaf_iter&i. 
		%end;
		;quit;
	%end;


%if %FileExist(&matrix.) %then %do;
/*Import Boosting Iterations*/
%if &jmnotes. %then %do;option nonotes nosource;%end;
PROC IMPORT OUT=tmp_mat DATAFILE= "&matrix." DBMS=DLM REPLACE;DELIMITER='2C'x; guessingrows=10001;GETNAMES=NO;DATAROW=1;RUN;

proc datasets nolist; modify tmp_mat;rename var1=irkey
%do i=0 %to &niter.;
	var%eval(&i.+2)=est_iter&i. 
%end;
;quit;

%getattrn_special(runmodel);
/*%put &=this_special_count.;*/
%if %getattrn(data=tmp_mat) ne &this_special_count. %then %do;
	%put ERROR: (bk) Iteration Matrix is incomplete (observations are missing). It is possible that Julia crashed (not enough RAM?) while 
	writing the file &matrix.;
%end;
%else %do;
proc sql noprint;
	create table tmp_mat2 as select i.*,j.numerator,j.denominator,j.numerator/j.denominator as observed,j.training_validation_bool as trn from tmp_mat as i
	left join runmodel as j on (i.irkey=j.irkey)
	order by j.training_validation_bool desc,irkey;quit;

data tmpiter;set tmp_mat2(/*keep=..*/);	
	%do this_iter=1 %to &niter.;	
		abs_difference_numerator&this_iter.=abs(numerator- est_iter&this_iter./denominator);
	%end;
run;

data itm;set tmpiter(drop=abs_:);run;

/*%put WARNING: (bk) this summary and all that follows is not meaningful yet (for ratios such as LR);*/
proc summary data=tmpiter;class trn;var abs_: numerator denominator est_iter: observed;output out=tmpout(drop=_type_ rename=_freq_=observations where=(trn ne .)) sum=;quit;

data iterations_res(drop=abs_:);
	set tmpout;
	%do this_iter=1 %to &niter.;
		error_numerator&this_iter.=abs_difference_numerator&this_iter./numerator;
	%end;
	%do this_iter=1 %to &niter.;
		error&this_iter.=observed/est_iter&this_iter.-1;
	%end;
run;

%if &jmnotes. %then %do;option notes source;%end;

	%if &standard_plot. %then %do;
	/*Plot it*/
		proc transpose data=iterations_res(keep=trn error:) out=for_graph;id trn;quit;
		data for_graph;set for_graph(rename=(_0=val _1=trn));iteration=_N_;run;

		title "&jname.";
		proc sgplot data=for_graph;
			series y=trn x=iteration / MARKERS LINEATTRS=(color=blue PATTERN=2) MARKERATTRS=(color=blue SYMBOL=circle);
			yaxis label="Error";
			xaxis label="Iteration";		
			series y=val x=iteration /MARKERS LINEATTRS=(color=blue PATTERN=35) MARKERATTRS=(color=blue SYMBOL=plus);
			legend;
		run;quit;
	%end;
%end;
%end;
%else %do;
	%if &model_type. eq boosted_tree %then %do;
		%put WARNING: (bk) File: &matrix. was not found!;
	%end;
%end;
%end; /*write_iteration_matrix=0*/

%if %FileExist(&juliadatafolder.\&outcsv.) %then %do;
	%if &jmnotes. %then %do;option nonotes nosource;%end;
	PROC IMPORT OUT=tmpres DATAFILE= "&juliadatafolder.\&outcsv." DBMS=DLM REPLACE;guessingrows=10001;DELIMITER='2C'x; GETNAMES=NO;DATAROW=1;RUN;	

	%if &model_type. eq boosted_tree %then %do;
		proc datasets nolist; modify tmpres;rename var1=irkey var2=score var3=unsmoothed_estimate var4=smoothed_estimate var5=rawEstimatePerObs;quit;
	%end;%else %do;
		proc datasets nolist; modify tmpres;rename var1=irkey var2=unsmoothed_estimate;quit; /*var3=smoothed_estimate;quit;*/
	%end;

	data tmp_order;set runmodel(keep=irkey);order_variable=_N_;run;
	proc sql noprint;
		create table &out.(drop=order_variable irkey2 irkey3) as select i.*,o.*,j.*,j.numerator/j.denominator as observed 
			from tmpres(%if &smoothEstimates.=1 %then %do;rename=(smoothed_estimate=est_from_trn)%end;%else %do;rename=(unsmoothed_estimate=est_from_trn)%end;) as i 
		left join runmodel(keep=irkey weight training_validation_bool numerator denominator rename=irkey=irkey2) as j on (j.irkey2=i.irkey)
		left join tmp_order(rename=irkey=irkey3) as o on (o.irkey3=i.irkey)
	order by order_variable;quit;

	data &out._oos(drop=training_validation_bool);set &out.(where=(training_validation_bool eq 0) keep=observed est_from_trn training_validation_bool irkey numerator denominator );rename est_from_trn=estimate;run;

	%if &jmnotes. %then %do;option notes source;%end;

	%if &bROSASProduceRkOptStats. %then %do;
		%rankoptimization_stats(&statsfile.,&jname.,&niter.,plot=1);	
		data rankopt;set Rankopt_&jname.;run;
	%end;
/*	%if &jmnotes. %then %do;option notes source;%end;*/
/*		proc datasets memtype=view nolist;*delete vtmp: vt:;quit;*/
/*		proc datasets nolist;*delete _cont_ tmp:;quit;*/
	*Produce Statistics for Testbench;
/*		%let var_dep=observed;*/
/*		%prepare_oos_error_output(in=res_oos);*/
/*	%if &jmnotes. %then %do;option notes;%end;*/

/*Calculate Error*/
	%PUT WARNING: (bk) these error statistics need to adjusted for ratios! do not consider the output for now!;
data inn2;set &out.(rename=(irkey=policyno est_from_trn=estimate training_validation_bool=trn) keep=irkey est_from_trn trainin: observed numerator denominator);	
	numerator_est=estimate/denominator;
	abs_difference_numerator=abs(numerator-numerator_est);
	if numerator not in (0 .) then abs_rel=abs_difference_numerator/numerator;
run;
proc summary data=inn2;class trn;var abs_difference_numerator numerator numerator_est denominator;output out=tmp_sum_out(drop=_type_ rename=_freq_=observations where=(trn ne .)) sum=;quit;
data final_est;
	set tmp_sum_out;
	estimate=numerator_est/denominator;
	observed=numerator/denominator;
	error=observed/estimate-1;
	error_numerator=abs_difference_numerator/numerator;
	abs_error=abs(error);
run;

%end;
%else %do;
	%put ERROR: (bk) Julia produced no output!;
%end;

%let outfname_of_julia_macro=&outfilename.;/*this is used by the save_julia_results macro*/
%let datafolder_of_julia_macro=&juliadatafolder.;

proc datasets nolist;*delete inn2 tmp:;quit;
ODS _ALL_ CLOSE;ODS HTML;ODS RESULTS ON;
%end;/*length of jname is ok*/
%end;/*runmodel is not empty*/
%end;/*runmodel exists*/
%end;/*stop_julia_macro ne 0*/
%mend julia_macro;
