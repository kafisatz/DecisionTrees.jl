%macro define_known_char_varlist(data=runmodel,var=&var_indep.,suffix=_VALS,file=c:\temp\tmp_varlist.txt,upcase=1);
	%PUT NOTE: (bk) Variable names must be shorter than %eval(32-%length(&suffix.)) otherwise this macro will fail.;
	%PUT NOTE: (bk) You can choose a shorter suffix value if you have longer variable names;
	%local i v vglobal nv vars charvars;

	%let charvars=;
	*Select charvars;
		proc sql noprint;
			select upcase(name) into: charvars separated by " " from  sashelp.vcolumn where memname="%sysfunc(upcase(&data.))" and libname="WORK" and type='char';
		quit;
		%let charvars=%upcase(&charvars.);
		%let vars=;
		%do i=1 %to %sysfunc(countw(&var.));
			%let v=%upcase(%scan(&var.,&i.));
			%if %index(%upcase(&charvars.),%upcase(&v.)) %then %let vars=&vars. &v.;	
		%end;
		%let vars=&vars.;
		%put &=vars.;

	%let nv=%sysfunc(countw(&vars.));
	%do i=1 %to &nv.;%let v=%sysfunc(upcase(%scan(&vars.,&i.)));%let vglobal=&v.&suffix.;
	/*	%put &=v.;	%put &=vglobal.; %put %symexist(&vglobal.);*/
		%if not %symexist(&vglobal.) %then %global &vglobal.;
		%let &vglobal.=;
	%end;
	sasfile &data. open;
	proc sql noprint;
		%do i=1 %to &nv.;%let v=%sysfunc(upcase(%scan(&vars.,&i.)));%let vglobal=&v.&suffix.;
		%if &upcase. %then %do;
			select distinct upcase(quote(trim(left(&v.))))   into: &vglobal. separated by " " from &data.;
		%end;
		%else %do;
			select distinct quote(trim(left(&v.)))   into: &vglobal. separated by " " from &data.;
		%end;			
		%end;
	quit;
	sasfile &data. close;

/*	%do i=1 %to &nv.;%let v=%sysfunc(upcase(%scan(&vars.,&i.)));%let vglobal=&v.&suffix.;*/
/*		%put &&&vglobal.;*/
/*		%let &vglobal.=%str(')&&&vglobal.%str(');*/
/*	%end;*/

/*	%do i=1 %to &nv.;%let v=%sysfunc(upcase(%scan(&vars.,&i.)));%let vglobal=&v.&suffix.;*/
/*		%put %nrstr(%let )&vglobal.=&&&vglobal.%str(;);*/
/*	%end;*/

	data _tmp_delme_;
		do i=1 to &nv.;
			output;
		end;
	run;
	filename tmp229x "&file.";
	data _tmp_delme_2;
/*	   format varname $32.;*/
		format msgline myquote $6000.;
		file tmp229x;
	   set _tmp_delme_;
	   %do i=1 %to &nv.;
	        %let v=%sysfunc(upcase(%scan(&vars.,&i.)));%let vglobal=&v.&suffix.;			
			%if %length(&&&vglobal.)<5900 %then %do;;
				if _N_=&i. then do;
					msgline=%unquote(%str(%'let &vglobal.=&&&vglobal.;%'));	
	/*"'%let '"||vname||val=%str(;);%str(');*/
				end;				
				%end;
				%else %do;
					%PUT NOTE: (bk) &v. results in a list that is longer than 5900 chars. No output created;
				%end;
		%end;
		myquote='%'||msgline;
	   if msgline ne "" then do;
	   		put myquote;
		end;
		else do;
		%do i=1 %to &nv.;
	        %let v=%sysfunc(upcase(%scan(&vars.,&i.)));%let vglobal=&v.&suffix.;			
			%if not %length(&&&vglobal.)<5900 %then %do;;
				if _N_=&i. then do;
					myquote="Remove condition for variable &vglobal.";	
					put myquote;
				end;
				%end;
		%end;
		
/*			myquote=*/
		end;
/*	   msgline='write this line';*/
/*	   put msgline;*/
	run;
	
%mend define_known_char_varlist;

/*
usage:
%let ZulassungsKanton_vals=32;
%put &=ZulassungsKanton_vals.;

%define_known_char_varlist(data=runmodel,var=&var_indep.);

*/
