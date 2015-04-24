function ahk_comboprofile,slitid,slitindex_blue600,yfitnorm_blue600,slitindex_red900,yfitnorm_red900,slitindex_red400,yfitnorm_red400, ccdside=ccdside

;+
; NAME:
;   ahk_comboprofile
;
; PURPOSE:
;   **Unlike ahk_comboprofile_nopad, this code will pad slits with zeros if off by more than 10% from what is expected from plate scale conversion.
;   **Requires additional input: ccdside = 0 if on left side of chip gap (i.e. pix LT 2048), = 1 if on right side of chip gap.
;   This will determine where the slits get anchored with respect to one another.
;   
;   For a particular slit, this code will:
;   1) Calculate normalized slitindex, since blue slitindex does not equal red. This will eventually need to change to be RA/Dec.
;   2) Multiply values from the n grisms for which values are available, at a particular normalized slitindex (RA/Dec)
;   3) Renormalize
;   Code allows for cases when not all 3 profiles are defined.
;
; CALLING SEQUENCE:
;  ahk_comboprofile(slitid,slitindex_blue600,yfitnorm_blue600,slitindex_red900,yfitnorm_red900,slitindex_red400,yfitnorm_red400, ccdside=0)
;
; INPUTS:
;
; OPTIONAL INPUTS:
;                
; OUTPUTS:
; normalized slit profile.
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
; REVISION HISTORY:
;   12-Oct-2014 -- Written by Alice Koning
;-

redplatescale = 0.21746 ;;[arcsec/pix]
blueplatescale = 0.135 ;;[arcsec/pix]

CASE 1 OF
	;; None of the three profiles are defined
	(N_ELEMENTS(yfitnorm_blue600) EQ 0 AND N_ELEMENTS(yfitnorm_red900) EQ 0 AND N_ELEMENTS(yfitnorm_red400) EQ 0): BEGIN
		print, 'None of the three profiles are defined'
		comboyfitnorm = MAKE_ARRAY(10,1, VALUE = 0)
	END

	;; Only Blue600 defined
	(N_ELEMENTS(yfitnorm_red900) EQ 0 AND N_ELEMENTS(yfitnorm_red400) EQ 0): BEGIN
		print, 'Only Blue600 defined'
		comboyfitnorm = yfitnorm_blue600
	END

	;; Only Red900 defined
	(N_ELEMENTS(yfitnorm_blue600) EQ 0 AND N_ELEMENTS(yfitnorm_red400) EQ 0): BEGIN
		print, 'Only Red900 defined'
		comboyfitnorm = yfitnorm_red900
	END

	;; Only Red400 defined
	(N_ELEMENTS(yfitnorm_blue600) EQ 0 AND N_ELEMENTS(yfitnorm_red900) EQ 0): BEGIN
		print, 'Only Red400 defined'
		comboyfitnorm = yfitnorm_red400
	END

	;; Blue600 not defined
	(N_ELEMENTS(yfitnorm_blue600) EQ 0): BEGIN
		print, 'Blue600 not defined'
		; Find the normalized x position along each slit
		; Since these are both red, first check to make sure they are the same length to within 10%. If not, raise the alarm!
		red900length = redplatescale*MAX(slitindex_red900)
		red400length = redplatescale*MAX(slitindex_red400)

		IF (red900length-red400length)/red900length GT 0.1 THEN BEGIN
			print, 'WARNING!! Red900 slit length longer than Red400 by more than 10 percent for slit '+StrTrim(slitid,2)
		ENDIF ELSE IF (red400length-red900length)/red400length GT 0.1 THEN BEGIN
			print, 'WARNING!! Red400 slit length longer than Red900 by more than 10 percent for slit '+StrTrim(slitid,2)			
		ENDIF

		xnorm_red900 = slitindex_red900/MAX(slitindex_red900)
		xnorm_red400 = slitindex_red400/MAX(slitindex_red400)

		; Make red900 slitindex the reference
		; so interpolate red400 y values, at xnorm_red900 slitindices.
		yinterp_red400 = INTERPOL(yfitnorm_red400,xnorm_red400,xnorm_red900)

		; Now compute the combined profile by multiplying all the normalized profiles together (at same xnorm_red900 slitindices)
		; Divide by n-th root of the number of profiles used (2 in this case) to keep correct width
		comboyfit=(yfitnorm_red900*yinterp_red400)^(1./2)

		; Renormalize
		area = int_tabulated(slitindex_red900,comboyfit) ;; Careful here! will need to change slitindex_blue600 to RA/Dec, once known.
		comboyfitnorm = comboyfit/abs(area)
	END

	;; Red900 not defined
	(N_ELEMENTS(yfitnorm_red900) EQ 0): BEGIN
		print, 'Red900 not defined'
		; Find the normalized x position along each slit
		; First check to make sure they are the same length to within 10% (after plate scale conversion). If not, pad with zeros.
		blue600length = blueplatescale*MAX(slitindex_blue600)
		red400length = redplatescale*MAX(slitindex_red400)

		IF (red400length-blue600length)/red400length GT 0.1 THEN BEGIN
			print, 'Red400 slit length longer than Blue600 by more than 10 percent for slit '+StrTrim(slitid,2)
			blue600reqpix=FIX(red400length/blueplatescale) ; Required length of blue600 in pixels
			yfitnorm_blue600_padded=MAKE_ARRAY(blue600reqpix, VALUE=0.0)
			diff = (blue600reqpix-MAX(slitindex_blue600))>1.0
			slitindex_blue600=findgen(blue600reqpix)
			IF ccdside THEN BEGIN
				;;On right side of CCD, so anchor blue and red at 0 then pad blue
				FOR i=0,(blue600reqpix-diff) DO yfitnorm_blue600_padded[i]=yfitnorm_blue600[i]
				yfitnorm_blue600 = yfitnorm_blue600_padded
			ENDIF ELSE BEGIN
				;;On left side of CCD, so anchor blue and red at 1.0 then pad blue
				FOR i=diff,blue600reqpix DO yfitnorm_blue600_padded[i-1]=yfitnorm_blue600[i-diff]
				yfitnorm_blue600 = yfitnorm_blue600_padded
			ENDELSE
		ENDIF ELSE IF (blue600length-red400length)/blue600length GT 0.1 THEN BEGIN
			print, 'Blue600 slit length longer than Red400 by more than 10 percent for slit '+StrTrim(slitid,2)
			red400reqpix=FIX(blue600length/redplatescale) ; Required length of red400 in pixels
			yfitnorm_red400_padded=MAKE_ARRAY(red400reqpix, VALUE=0.0)
			diff = (red400reqpix-MAX(slitindex_red400))>1.0
			slitindex_red400=findgen(red400reqpix)
			IF ccdside THEN BEGIN
				;;On right side of CCD, so anchor blue and red at 0 then pad red
				FOR i=0,(red400reqpix-diff) DO yfitnorm_red400_padded[i]=yfitnorm_red400[i]
				yfitnorm_red400 = yfitnorm_red400_padded
			ENDIF ELSE BEGIN
				;;On left side of CCD, so anchor blue and red at 1.0 then pad red
				FOR i=diff,red400reqpix DO yfitnorm_red400_padded[i-1]=yfitnorm_red400[i-diff]
				yfitnorm_red400 = yfitnorm_red400_padded
			ENDELSE
		ENDIF

		xnorm_blue600 = slitindex_blue600/MAX(slitindex_blue600)
		xnorm_red400 = slitindex_red400/MAX(slitindex_red400)

		; Interpolate red400 y values, at xnorm_blue600 slitindices.
		yinterp_red400 = INTERPOL(yfitnorm_red400,xnorm_red400,xnorm_blue600)

		; Now compute the combined profile by multiplying all the normalized profiles together (at same xnorm_blue600 slitindices)
		; Divide by n-th root of the number of profiles used (2 in this case) to keep correct width
		comboyfit=(yfitnorm_blue600*yinterp_red400)^(1./2)

		; Renormalize
		area = int_tabulated(slitindex_blue600,comboyfit) ;; Careful here! will need to change slitindex_blue600 to RA/Dec, once known.
		comboyfitnorm = comboyfit/abs(area)
	END

	;; Red400 not defined
	(N_ELEMENTS(yfitnorm_red400) EQ 0): BEGIN
		print, 'Red400 not defined'
		; Find the normalized x position along each slit
		; First check to make sure they are the same length to within 10% (after plate scale conversion). If not, pad with zeros.
		blue600length = blueplatescale*MAX(slitindex_blue600)
		red900length = redplatescale*MAX(slitindex_red900)

		IF (red900length-blue600length)/red900length GT 0.1 THEN BEGIN
			print, 'Red900 slit length longer than Blue600 by more than 10 percent for slit '+StrTrim(slitid,2)
			blue600reqpix=FIX(red900length/blueplatescale) ; Required length of blue600 in pixels
			yfitnorm_blue600_padded=MAKE_ARRAY(blue600reqpix, VALUE=0.0)
			diff = (blue600reqpix-MAX(slitindex_blue600))>1.0 ; Difference between required and old length of blue600 in pixels
			slitindex_blue600=findgen(blue600reqpix) ; Replace slitindex vector with new values spanning the required length
			IF ccdside THEN BEGIN
				;;On right side of CCD, so anchor blue and red at 0 then pad blue
				FOR i=0,(blue600reqpix-diff) DO yfitnorm_blue600_padded[i]=yfitnorm_blue600[i]
				yfitnorm_blue600 = yfitnorm_blue600_padded
			ENDIF ELSE BEGIN
				;;On left side of CCD, so anchor blue and red at 1.0 then pad blue
				FOR i=diff,blue600reqpix DO yfitnorm_blue600_padded[i-1]=yfitnorm_blue600[i-diff]
				yfitnorm_blue600 = yfitnorm_blue600_padded
			ENDELSE
		ENDIF ELSE IF (blue600length-red900length)/blue600length GT 0.1 THEN BEGIN
			print, 'Blue600 slit length longer than Red900 by more than 10 percent for slit '+StrTrim(slitid,2)
			red900reqpix=FIX(blue600length/redplatescale) ; Required length of red900 in pixels
			yfitnorm_red900_padded=MAKE_ARRAY(red900reqpix, VALUE=0.0)
			diff = (red900reqpix-MAX(slitindex_red900))>1.0 ; Difference between required and old length of red900 in pixels
			slitindex_red900=findgen(red900reqpix) ; Replace slitindex vector with new values spanning the required length
			IF ccdside THEN BEGIN
				;;On right side of CCD, so anchor blue and red at 0 then pad red
				FOR i=0,(red900reqpix-diff) DO yfitnorm_red900_padded[i]=yfitnorm_red900[i]
				yfitnorm_red900 = yfitnorm_red900_padded
			ENDIF ELSE BEGIN
				;;On left side of CCD, so anchor blue and red at 1.0 then pad red
				FOR i=diff,red900reqpix DO yfitnorm_red900_padded[i-1]=yfitnorm_red900[i-diff]
				yfitnorm_red900 = yfitnorm_red900_padded
			ENDELSE
		ENDIF

		xnorm_blue600 = slitindex_blue600/MAX(slitindex_blue600)
		xnorm_red900 = slitindex_red900/MAX(slitindex_red900)

		; Interpolate red900 y values, at xnorm_blue600 slitindices.
		yinterp_red900 = INTERPOL(yfitnorm_red900,xnorm_red900,xnorm_blue600)

		; Now compute the combined profile by multiplying all the normalized profiles together (at same xnorm_blue600 slitindices)
		; Divide by n-th root of the number of profiles used (2 in this case) to keep correct width
		comboyfit=(yfitnorm_blue600*yinterp_red900)^(1./2)
		
		; Renormalize
		area = int_tabulated(slitindex_blue600,comboyfit) ;; Careful here! will need to change slitindex_blue600 to RA/Dec, once known.
		comboyfitnorm = comboyfit/abs(area)

	END

	;; All profiles defined
	ELSE: BEGIN
		print, 'All profiles defined'
		; Find the normalized x position along each slit
		; First check to make sure they are the same length to within 10% (after plate scale conversion). If not, pad with zeros.
		blue600length = blueplatescale*MAX(slitindex_blue600)
		red900length = redplatescale*MAX(slitindex_red900)
		red400length = redplatescale*MAX(slitindex_red400)

		IF (red900length-red400length)/red900length GT 0.1 THEN BEGIN
			print, 'WARNING!! Red900 slit length longer than Red400 by more than 10 percent for slit '+StrTrim(slitid,2)
		ENDIF ELSE IF (red400length-red900length)/red400length GT 0.1 THEN BEGIN
			print, 'WARNING!! Red400 slit length longer than Red900 by more than 10 percent for slit '+StrTrim(slitid,2)			
		ENDIF

		IF (red400length-blue600length)/red400length GT 0.1 THEN BEGIN
			print, 'Red400 slit length longer than Blue600 by more than 10 percent for slit '+StrTrim(slitid,2)
			blue600reqpix=FIX(red400length/blueplatescale) ; Required length of blue600 in pixels
			yfitnorm_blue600_padded=MAKE_ARRAY(blue600reqpix, VALUE=0.0)
			diff = (blue600reqpix-MAX(slitindex_blue600))>1.0
			slitindex_blue600=findgen(blue600reqpix)
			IF ccdside THEN BEGIN
                           ;;On right side of CCD, so anchor blue and red at 0 then pad blue
                           FOR i=diff,blue600reqpix DO yfitnorm_blue600_padded[i-1]=yfitnorm_blue600[i-diff]
                           yfitnorm_blue600 = yfitnorm_blue600_padded
                           ;FOR i=0,(blue600reqpix-diff) DO yfitnorm_blue600_padded[i]=yfitnorm_blue600[i]
                           ;yfitnorm_blue600 = yfitnorm_blue600_padded
			ENDIF ELSE BEGIN
                           ;;On left side of CCD, so anchor blue and
                           ;;red at 1.0 then pad blue
                           FOR i=0,(blue600reqpix-diff) DO yfitnorm_blue600_padded[i]=yfitnorm_blue600[i]
                           yfitnorm_blue600 = yfitnorm_blue600_padded
                           ;FOR i=diff,blue600reqpix DO yfitnorm_blue600_padded[i-1]=yfitnorm_blue600[i-diff]
                           ;yfitnorm_blue600 = yfitnorm_blue600_padded
			ENDELSE
		ENDIF ELSE IF (blue600length-red400length)/blue600length GT 0.1 THEN BEGIN
			print, 'Blue600 slit length longer than Red400 by more than 10 percent for slit '+StrTrim(slitid,2)
			red400reqpix=FIX(blue600length/redplatescale) ; Required length of red400 in pixels
			yfitnorm_red400_padded=MAKE_ARRAY(red400reqpix, VALUE=0.0)
			diff = (red400reqpix-MAX(slitindex_red400))>1.0
			slitindex_red400=findgen(red400reqpix)
			IF ccdside THEN BEGIN
				;;On right side of CCD, so anchor blue and red at 0 then pad red
				FOR i=0,(red400reqpix-diff) DO yfitnorm_red400_padded[i]=yfitnorm_red400[i]
				yfitnorm_red400 = yfitnorm_red400_padded
			ENDIF ELSE BEGIN
				;;On left side of CCD, so anchor blue and red at 1.0 then pad red
				FOR i=diff,red400reqpix DO yfitnorm_red400_padded[i-1]=yfitnorm_red400[i-diff]
				yfitnorm_red400 = yfitnorm_red400_padded
			ENDELSE
		ENDIF

		IF (red900length-blue600length)/red900length GT 0.1 THEN BEGIN
			print, 'Red900 slit length longer than Blue600 by more than 10 percent for slit '+StrTrim(slitid,2)
			blue600reqpix=FIX(red900length/blueplatescale) ; Required length of blue600 in pixels
			yfitnorm_blue600_padded=MAKE_ARRAY(blue600reqpix, VALUE=0.0)
			diff = (blue600reqpix-MAX(slitindex_blue600))>1.0
			slitindex_blue600=findgen(blue600reqpix)
			IF ccdside THEN BEGIN
                           ;;On right side of CCD, so anchor blue and red at 0 then pad blue
                           FOR i=diff,blue600reqpix DO yfitnorm_blue600_padded[i-1]=yfitnorm_blue600[i-diff]
                           yfitnorm_blue600 = yfitnorm_blue600_padded
                           ;FOR i=0,(blue600reqpix-diff) DO yfitnorm_blue600_padded[i]=yfitnorm_blue600[i]
                           ;yfitnorm_blue600 = yfitnorm_blue600_padded
			ENDIF ELSE BEGIN
                           ;;On left side of CCD, so anchor blue and
                           ;;red at 1.0 then pad blue
                           FOR i=0,(blue600reqpix-diff) DO yfitnorm_blue600_padded[i]=yfitnorm_blue600[i]
                           yfitnorm_blue600 = yfitnorm_blue600_padded
                           ;FOR i=diff,blue600reqpix DO yfitnorm_blue600_padded[i-1]=yfitnorm_blue600[i-diff]
                           ;yfitnorm_blue600 = yfitnorm_blue600_padded
			ENDELSE
		ENDIF ELSE IF (blue600length-red900length)/blue600length GT 0.1 THEN BEGIN
			print, 'Blue600 slit length longer than Red900 by more than 10 percent for slit '+StrTrim(slitid,2)
			red900reqpix=FIX(blue600length/redplatescale) ; Required length of red900 in pixels
			yfitnorm_red900_padded=MAKE_ARRAY(red900reqpix, VALUE=0.0)
			diff = (red900reqpix-MAX(slitindex_red900))>1.0
			slitindex_red900=findgen(red900reqpix)
			IF ccdside THEN BEGIN
				;;On right side of CCD, so anchor blue and red at 0 then pad red
				FOR i=0,(red900reqpix-diff) DO yfitnorm_red900_padded[i]=yfitnorm_red900[i]
				yfitnorm_red900 = yfitnorm_red900_padded
			ENDIF ELSE BEGIN
				;;On left side of CCD, so anchor blue and red at 1.0 then pad red
				FOR i=diff,red900reqpix DO yfitnorm_red900_padded[i-1]=yfitnorm_red900[i-diff]
				yfitnorm_red900 = yfitnorm_red900_padded
			ENDELSE
		ENDIF

                xnorm_blue600 = slitindex_blue600/MAX(slitindex_blue600)
		xnorm_red900 = slitindex_red900/MAX(slitindex_red900)
		xnorm_red400 = slitindex_red400/MAX(slitindex_red400)
                
		; Make blue slitindex the reference (because expect red slits to be cut-off at edges, not blue)
		; so interpolate red900 and red400 y values, at xnorm_blue600 slitindices.
		yinterp_red900 = INTERPOL(yfitnorm_red900,xnorm_red900,xnorm_blue600)
		yinterp_red400 = INTERPOL(yfitnorm_red400,xnorm_red400,xnorm_blue600)

		pblue = PLOT(yfitnorm_blue600, 'b', title=slitid)
                pred9 =  PLOT(yinterp_red900, 'r', /OVERPLOT)
                pred4 =  PLOT(yinterp_red400, /OVERPLOT)
                pblue.Save, 'allProfiles_switchpadside_slitid'+StrTrim(slitid,2)+'.png'

                ; Force profiles that will be used to computer comboprofile to be strictly positive
                yfitnorm_blue600 = yfitnorm_blue600 > 0.0
                yinterp_red900 = yinterp_red900 > 0.0
                yinterp_red400 = yinterp_red400 > 0.0                

		; Now compute the combined profile by multiplying all the normalized profiles together (at same xnorm_blue600 slitindices)
		; Divide by n-th root of the number of profiles used (3 in this case) to keep correct width		
                comboyfit=(yfitnorm_blue600*yinterp_red900*yinterp_red400)^(1./3)

		; Renormalize
		area = int_tabulated(slitindex_blue600,comboyfit) ;; Careful here! will need to change slitindex_blue600 to RA/Dec, once known.
		comboyfitnorm = comboyfit/abs(area)
	END
ENDCASE



p2 = PLOT(comboyfitnorm, 'b', thick=5, title=slitid)
p2.Save, 'comboyfitnorm_slitid'+StrTrim(slitid,2)+'_exp3.png'

return, comboyfitnorm

end
