pro ahk_Gauss_profile_singlegrating,scifile,slitfile=slitfile,wavefile=wavefile

;+
; NAME:
;   ahk_Gauss_profile_singlegrating
;
; PURPOSE:
;   From output of long_reduce, fit Gaussian model across each slit profile then extract objects accordingly.
;   General outline -- All steps within loop over slits:
;	1) Restore saved objstr.FLUXMODEL (1-d array of summed fluxes in our lines of interest)
;	2) Fit FLUXMODEL with sum of Gaussians using ahk_profile.pro. Get out 4096x4096 profile.
;	3) Feed profile into long_extract_optimal
;
; CALLING SEQUENCE:
;	ahk_Gauss_profile_singlegrating,'Science/sci-lblue2082.fits.gz',slitfile='slits-lblue2039.fits',wavefile='wave-lblue2204.fits'
;
; INPUTS:
;
; OPTIONAL INPUTS:
;                
; OUTPUTS:
;	Saves csv for each slit with slitindex, yfitnorm as field1, field2, respectively. This is what gets used by ahk_comboprofile
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
;	ahk_profile_singleslit
;
; REVISION HISTORY:
;   3-November-2014 -- Written by Alice Koning
;-

;; Restore saved objstr.FLUXMODEL
RESTORE, repstr(scifile,'fits.gz','sav')

;; Get various required extensions out of science file
ivar = xmrdfits(scifile,1)
skyimage = xmrdfits(scifile,2)
objimage = xmrdfits(scifile,3)
outmask = xmrdfits(scifile,4)
waveimg = xmrdfits(wavefile, silent = (keyword_set(verbose) EQ 0), 0)
sciimg = xmrdfits(scifile)
imgminsky = sciimg - skyimage

;; Size of science image
nx = (size(sciimg))[1]
ny = (size(sciimg))[2]
xarr = findgen(nx)## replicate(1.0, ny)
yarr = findgen(ny)## replicate(1.0, nx)

;; Number of slits in image
slitmask = mrdfits(slitfile)
nslit = max(slitmask,/NaN)

;; Find best fit Gaussian curve for each slit.
FOR ii =1,nslit DO BEGIN
	;;Find first object in final_struct with slitid equal to ii (avoids duplicates for slits with more than one object)
	structid = min(WHERE(final_struct.slitid EQ (ii)))
	
	;;Start next iteration if ii is not a slitid in final_struct
	IF (structid EQ -1) THEN CONTINUE

	slitid = final_struct[structid].SLITID
	print, '***** SLITID: ', slitid, ' *****'
	slitindex = [] ;; Array to hold values of index along slit. May not need in final code.
	yfitnorm = [] ;; Array to hold normalized profile of each slit. May not need in final code.
	yfitparams = [] ;; Array to hold parameters that describe yfitnorm of each slit. May not need in final code.

	;; Get object profile along slit
	yfitnotnorm = ahk_Gauss_profile_singleslit(*final_struct[structid].FLUXMODEL,slitindex,yfitnorm,yfitparams,scifile=scifile,slitid=slitid)
        
        
	IF (finite(yfitnorm(1))) THEN BEGIN
		;; Save normalized profile and slit position to file (to be used for calculating average normalized profile across all 3 gratings)
		WRITE_CSV, 'yfitnorm_'+strmid(scifile,15,8,/reverse_offset)+'_slitid'+StrTrim(slitid,2)+'.csv', slitindex, yfitnorm
	ENDIF
ENDFOR ;; end loop over slits


end


