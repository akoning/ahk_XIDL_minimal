pro ahk_profile_singlegrating,scifile,slitfile=slitfile,wavefile=wavefile

;+
; NAME:
;   ahk_profile_singlegrating
;
; PURPOSE:
;  
;
; CALLING SEQUENCE:
;	ahk_profile_singlegrating,'../Science_blue600/sci-lblue2082.fits.gz',slitfile='../slits-lblue2039.fits',wavefile='../wave-lblue2204.fits'
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
	print, '***** Working on slit: ', slitid, ' *****'

        profile = *final_struct[structid].FLUXMODEL
        slitindex = findgen(N_ELEMENTS(profile))

	;; Normalize object profile along slit
	area = int_tabulated(slitindex,profile)
	normprofile = profile/abs(area)

	p = PLOT(slitindex, profile, yrange=[0.9*min(profile),1.1*max(profile)], title=slitid)
	p2 = PLOT(slitindex, normprofile, 'b', title=slitid)
        
	;;Save plot for profileCompareTests
	p.Save, 'profile_'+strmid(scifile,15,8,/reverse_offset)+'_slitid'+StrTrim(slitid,2)+'.png'
	p2.Save, 'normprofile_'+strmid(scifile,15,8,/reverse_offset)+'_slitid'+StrTrim(slitid,2)+'.png'

	IF (finite(normprofile(1))) THEN BEGIN
		;; Save normalized profile and slit position to file (to be used for calculating average normalized profile across all 3 gratings)
		WRITE_CSV, 'profile_'+strmid(scifile,15,8,/reverse_offset)+'_slitid'+StrTrim(slitid,2)+'.csv', slitindex, profile
		WRITE_CSV, 'normprofile_'+strmid(scifile,15,8,/reverse_offset)+'_slitid'+StrTrim(slitid,2)+'.csv', slitindex, normprofile
	ENDIF
ENDFOR ;; end loop over slits


end


