/*
Sets nonfrequent values to "~" for all character variables in the data set.
The parameter maxvalues determines how many disticnt values are kept.
The macro will keep the most frequent values.

The "~" is hardcoded below. It could be changed if necessary. 
Note that issues with variable length (truncating) may appear if the ~ is replaced with something else.

Usage:

%define_other_for_rare_values(data=sashelp.cars,out=test,maxvalues=6);

*/

%macro define_other_for_rare_values(data=tmp,out=tmp_out,maxvalues=7);
	%let maxvalues=%eval(&maxvalues.-1);/* the additional level for "~" or "OTHER" will be added*/
	%local i varlistc c v;

	sasfile &data. open;
	%let varlistc=;
	proc freq data=&data. nlevels ;table _character_ / noprint;ods output nlevels=_tmp_(rename=(nlevels=distinct TableVar=_name_));run;
	proc sql noprint; select _name_ into : varlistc separated by ' ' from _tmp_ where distinct gt &maxvalues.;quit;

	%let varlistc=&varlistc.;
	%if %length(&varlistc.) %then %do;
		/*otherwise there are no character variables*/
		proc datasets nolist;delete tmpxx23_:;quit;
		%let c=%sysfunc(countw(&varlistc.));
		%if &c. gt 0 %then %do;	
			proc sql noprint outobs=&maxvalues. nowarn;
				create table tmp_distinctx as select "distinct_values" as distinct_values
				%do i=1 %to &c.;%let v=%scan(&varlistc.,&i.);
					,count(distinct(&v.)) as &v. 
				%end;				
				from &data.;
				%do i=1 %to &c.;%let v=%scan(&varlistc.,&i.);
					create table tmpxx23_&i. as select "&v." as name format=$32. length=32, &v. as value format=$50. length=600,count(*) as freq from &data. group by &v. order by freq desc;
					select value into :tmpvarlist&i. separated by '" "' notrim from tmpxx23_&i.;
				%end;
			;quit;
		%end;
		data &out.;
			set &data.;
			%do i=1 %to &c.;%let v=%scan(&varlistc.,&i.);
				if not (&v. in ("&&tmpvarlist&i.")) then &v.="~"; *OTHER or ___OTHER___ would be an alternative here however, the length/format of the variable might be too short; 
			%end;
		ruN;
		proc  transpose data=tmp_distinctx out=tmp_distinct;quit;
		proc sort data=tmp_distinct(rename=(col1=distinct_values));by descending distinct_values;quit;
		proc datasets nolist;*delete _tmp_ tmp_distinctx tmpxx23_:;quit;
	%end;
	%else %do;
		data &out.;
			set &data.;
		run;
	%end;
	sasfile &data. close;
%mend;

