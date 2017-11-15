%macro export_for_julia(data=runmodel,file=);
	%local len end;
	%if not %length(&file.) or not %length(data.) %then %do;
		%PUT ERROR: (bk) Export failed. No file or data parameter was found in macro export_for_julia;
	%end;
	%else %do;
		%let len=%length(&file.);
		%let end=%upcase(%sysfunc(lowcase(%substr(&file.,&len.-2,3))));
/*		%put &=end.;	*/
		%if &end. ne .DB and &end. ne SQL %then %do;
			%PUT NOTE: (bk) Exporting to CSV: &file.;
			%export_csv_comma_utf8(data=&data.,file=&file.);
		%end;
		%else %do;
			%if &end. eq SQL %then %do;
				%PUT NOTE: (bk) Exporting to SQL Database;				
				%writeSQLDB_smallspaceSQL_PT(data=&data.);
			%end;
			%else %do;
				%PUT NOTE: (bk) Exporting to SQLiteDB: &file.;
				%writetosqlite(data=&data.,file=&file.);
			%end;
		%end;
	%end;
%mend;

/*%export_for_julia(data=runmodel,file=c:\temp\mira.Dx);*/
