pro ahk_plot3csvprofiles,stem1,stem2,stem3, slitidmin=slitidmin, slitidmax=slitidmax

;+
; NAME:
;   ahk_plot3csvprofiles
;
; PURPOSE:
;   Searches the current directory for csv files beginning with strings stem1, stem2, stem3. These are expected to be the output from ahk_objextract_singlegrating
;
; CALLING SEQUENCE:
;	ahk_plot3csvprofiles, 'yfitnorm_blue2082_slitid', 'yfitnorm_blue2083_slitid', 'yfitnorm_blue2086_slitid', slitidmin=2,slitidmax=2
;	ahk_plot3csvprofiles, 'yfitnorm_lred2068_slitid', 'yfitnorm_lred2069_slitid', 'yfitnorm_lred2070_slitid', slitidmin=1,slitidmax=21
;
; INPUTS:
;
; OPTIONAL INPUTS:
;                
; OUTPUTS:
;
;
; OPTIONAL OUTPUTS:
;   
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;	
;
; REVISION HISTORY:
;   18-November-2014 -- Written by Alice Koning
;-
FOR slitid =slitidmin,slitidmax DO BEGIN

	IF FILE_TEST(stem1+STRTRIM(slitid,2)+'.csv') THEN p1 = read_csv(stem1+STRTRIM(slitid,2)+'.csv') ELSE p1={field1:[0,0],field2:[0,0]}
	IF FILE_TEST(stem2+STRTRIM(slitid,2)+'.csv') THEN p2 = read_csv(stem2+STRTRIM(slitid,2)+'.csv') ELSE p2={field1:[0,0],field2:[0,0]}
	IF FILE_TEST(stem3+STRTRIM(slitid,2)+'.csv') THEN p3 = read_csv(stem3+STRTRIM(slitid,2)+'.csv') ELSE p3={field1:[0,0],field2:[0,0]}


	plt = plot(p1.field1, p1.field2, 'b', title='slitid ' +STRTRIM(slitid,2))
	plt = plot(p2.field1, p2.field2, 'r', /OVERPLOT)
	plt = plot(p3.field1, p3.field2, /OVERPLOT)


	;pltdif = plot(p1.field1, p1.field2-p2.field2, title='1-2: slitid ' +STRTRIM(slitid,2))
	;pltdif = plot(p1.field1, p1.field2-p3.field2, title='1-3: slitid ' +STRTRIM(slitid,2))
	;pltdif = plot(p1.field1, p2.field2-p3.field2, title='2-3: slitid ' +STRTRIM(slitid,2))


ENDFOR


end
