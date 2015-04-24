pro ahk_plot3spectra,stem1,stem2,stem3, slitidmin=slitidmin, slitidmax=slitidmax

;+
; NAME:
;   ahk_plot3spectra
;
; PURPOSE:
;   Searches the current directory for fits files beginning with strings stem1, stem2, stem3. These are expected to contain a structure in the 5th fits extension
;   (with tags flux_opt and wave_opt in the structure)
;   NB: for now, only uses filenames that end in obj1, e.g. 'blue2082_slit1_obj1.fits' 
;
; CALLING SEQUENCE:
;	ahk_plot3spectra, 'blue2082_slit', 'blue2083_slit', 'blue2086_slit', slitidmin=i,slitidmax=i
;	ahk_plot3spectra, 'lred2068_slit', 'lred2069_slit', 'lred2070_slit', slitidmin=1,slitidmax=21
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

	IF FILE_TEST(stem1+STRTRIM(slitid,2)+'_obj1.fits') THEN p1 = mrdfits(stem1+STRTRIM(slitid,2)+'_obj1.fits',5,hd) ELSE p1={wave_opt:[0,0],flux_opt:[0,0]}
	IF FILE_TEST(stem2+STRTRIM(slitid,2)+'_obj1.fits') THEN p2 = mrdfits(stem2+STRTRIM(slitid,2)+'_obj1.fits',5,hd) ELSE p2={wave_opt:[0,0],flux_opt:[0,0]}
	IF FILE_TEST(stem3+STRTRIM(slitid,2)+'_obj1.fits') THEN p3 = mrdfits(stem3+STRTRIM(slitid,2)+'_obj1.fits',5,hd) ELSE p3={wave_opt:[0,0],flux_opt:[0,0]}


	plt = plot(p1.wave_opt, p1.flux_opt, 'b', title='slitid ' +STRTRIM(slitid,2))
	plt = plot(p2.wave_opt, p2.flux_opt, 'r', /OVERPLOT)
	plt = plot(p3.wave_opt, p3.flux_opt, /OVERPLOT)


	;pltdif = plot(p1.wave_opt, p1.flux_opt-p2.flux_opt, title='1-2: slitid ' +STRTRIM(slitid,2))
	;pltdif = plot(p1.wave_opt, p1.flux_opt-p3.flux_opt, title='1-3: slitid ' +STRTRIM(slitid,2))
	;pltdif = plot(p1.wave_opt, p2.flux_opt-p3.flux_opt, title='2-3: slitid ' +STRTRIM(slitid,2))


ENDFOR


end
