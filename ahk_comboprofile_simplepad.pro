function ahk_comboprofile_simplepad,slitid,slitindex_blue600,yfitnorm_blue600,slitindex_red900,yfitnorm_red900,slitindex_red400,yfitnorm_red400, ccdflagblue=ccdflagblue ,ccdflagred=ccdflagred

;+
; NAME:
;   ahk_comboprofile
;
; PURPOSE:
;   **Similar to ahk_comboprofile, but uses additional input of
;   ccdflagblue and ccdflagred instead of just ccdside. Each ccdflag
;   will be 1 if padding should go on the right of the slit, 2 if the
;   padding should go on the left of the slit and 0 if no padding is to be applied.
;   **This code is also simplified in that it assumes all the slits
;   have a defined fluxmodel.
;
;   For a particular slit, this code will:
;   1) Calculate normalized slitindex, since blue slitindex does not equal red.
;   2) Multiply values from the n grisms for which values are available, at a particular normalized slitindex (RA/Dec)
;   3) Renormalize
;
; CALLING SEQUENCE:
;  ahk_comboprofile(slitid,slitindex_blue600,yfitnorm_blue600,slitindex_red900,yfitnorm_red900,slitindex_red400,yfitnorm_red400,ccdflagblue=ccdflagblue[ii-1],ccdflagred=ccdflagred[ii-1])
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

blue600length = blueplatescale*MAX(slitindex_blue600)
red900length = redplatescale*MAX(slitindex_red900)
red400length = redplatescale*MAX(slitindex_red400)

IF (red900length-red400length)/red900length GT 0.1 OR (red400length-red900length)/red400length GT 0.1 THEN BEGIN
   print, 'WARNING!! Red900 slit length different from Red400 by more than 10 percent for slit '+StrTrim(slitid,2)
ENDIF

; Will have one of three cases, determined from ccdflagblue and ccdflagred:
; 1) Blue slit needs padding (ccdflagblue NE 0)
; 2) Red slits need padding (ccdflagred NE 0)
; 3) None of the slits need padding

CASE 1 OF
	;Blue slit needs padding
	(ccdflagblue NE 0 AND ccdflagred EQ 0): BEGIN
		print, 'Padding blue slit '+StrTrim(slitid,2)

                ;Pad blue to average length of the two red slits
                avgredlength = (red900length + red400length)/2
                bluereqpix=FIX(avgredlength/blueplatescale)
                yfitnorm_blue600_padded=MAKE_ARRAY(bluereqpix, VALUE=0.0)
		diff = (bluereqpix-MAX(slitindex_blue600))>1.0
		slitindex_blue600=findgen(bluereqpix)

                IF ccdflagblue EQ 2 THEN BEGIN
                   ;;Need to pad on left side, so anchor blue and red at 1
                   FOR i=diff,bluereqpix DO yfitnorm_blue600_padded[i-1]=yfitnorm_blue600[i-diff]
                   yfitnorm_blue600 = yfitnorm_blue600_padded
		ENDIF ELSE IF ccdflagblue EQ 1 THEN BEGIN
                    ;;Need to pad on right side, so anchor blue and red at 0
                   FOR i=0,(bluereqpix-diff) DO yfitnorm_blue600_padded[i]=yfitnorm_blue600[i]
                   yfitnorm_blue600 = yfitnorm_blue600_padded
                ENDIF ELSE BEGIN
                   print, 'Do not recognize ccdflagblue. Exiting.'
                   return, 0
                ENDELSE

	END
	;Red slits need padding
	(ccdflagblue EQ 0 AND ccdflagred NE 0): BEGIN
		print, 'Padding red slit '+StrTrim(slitid,2)

                ;Pad both red slits to length of the blue slit
                redreqpix=FIX(blue600length/redplatescale)
                yfitnorm_red900_padded=MAKE_ARRAY(redreqpix, VALUE=0.0)
                yfitnorm_red400_padded=MAKE_ARRAY(redreqpix, VALUE=0.0)
		diff900 = (redreqpix-MAX(slitindex_red900))>1.0
		diff400 = (redreqpix-MAX(slitindex_red400))>1.0
		slitindex_red900=findgen(redreqpix)
		slitindex_red400=findgen(redreqpix)

                IF ccdflagred EQ 2 THEN BEGIN
                   ;;Need to pad on left side, so anchor blue and red at 1
                   FOR i=diff900,redreqpix DO yfitnorm_red900_padded[i-1]=yfitnorm_red900[i-diff900]
                   yfitnorm_red900 = yfitnorm_red900_padded
                   FOR i=diff400,redreqpix DO yfitnorm_red400_padded[i-1]=yfitnorm_red400[i-diff400]
                   yfitnorm_red400 = yfitnorm_red400_padded
		ENDIF ELSE IF ccdflagred EQ 1 THEN BEGIN
                    ;;Need to pad on right side, so anchor blue and red at 0
                   FOR i=0,(redreqpix-diff900) DO yfitnorm_red900_padded[i]=yfitnorm_red900[i]
                   yfitnorm_red900 = yfitnorm_red900_padded
                   FOR i=0,(redreqpix-diff400) DO yfitnorm_red400_padded[i]=yfitnorm_red400[i]
                   yfitnorm_red400 = yfitnorm_red400_padded
                ENDIF ELSE BEGIN
                   print, 'Do not recognize ccdflagred. Exiting.'
                   return, 0
                ENDELSE

	END
	;None of the slits need padding
	(ccdflagblue EQ 0  AND ccdflagred EQ 0): BEGIN
		print, 'No padding for slit '+StrTrim(slitid,2)
        END

	ELSE: BEGIN
		print, 'Do not recognize ccdflags for slit '+StrTrim(slitid,2)
        END
ENDCASE

xnorm_blue600 = slitindex_blue600/MAX(slitindex_blue600)
xnorm_red900 = slitindex_red900/MAX(slitindex_red900)
xnorm_red400 = slitindex_red400/MAX(slitindex_red400)
                
; Interpolate red900 and red400 y values at xnorm_blue600 slitindices.
yinterp_red900 = INTERPOL(yfitnorm_red900,xnorm_red900,xnorm_blue600)
yinterp_red400 = INTERPOL(yfitnorm_red400,xnorm_red400,xnorm_blue600)

pall1 = PLOT(yfitnorm_blue600, 'b', title=slitid)
pall2 =  PLOT(yinterp_red900, 'r', /OVERPLOT)
pall3 =  PLOT(yinterp_red400, /OVERPLOT)
pall1.Save, 'allProfiles_simplepad_slitid'+StrTrim(slitid,2)+'_exp1.png'

; Force profiles that will be used to compute comboprofile to be
; strictly positive
print, 'Min criterium:'
print,  (1./MAX(slitindex_blue600))
yfitnorm_blue600 = yfitnorm_blue600 > (1./MAX(slitindex_blue600))
yinterp_red900 = yinterp_red900 > (1./MAX(slitindex_blue600))
yinterp_red400 = yinterp_red400 > (1./MAX(slitindex_blue600))                

; Now compute the combined profile by multiplying all the normalized profiles together (at same xnorm_blue600 slitindices)
; Divide by n-th root of the number of profiles used to keep correct width		
comboyfit=(yfitnorm_blue600*yinterp_red900*yinterp_red400)^(1./3)

; Subtract the minimum so we have a baseline of zero.
comboyfit = comboyfit - min(comboyfit)

; Smooth it!
boxcar = 2.0 ;[arcsec] this is approx 2x the seeing
boxcarpix = boxcar/blueplatescale ;[pix]
comboyfitsmooth = SMOOTH(comboyfit, boxcarpix, /edge_truncate)

; Renormalize
area = int_tabulated(slitindex_blue600,comboyfitsmooth)
comboyfitsmoothnorm = comboyfitsmooth/abs(area)

; Make the plot and save
pcombosmooth = PLOT(comboyfitsmoothnorm, 'r', thick=5, title=slitid)
pcombosmooth.Save, 'comboyfitnorm_simplepad_smooth_nelementMask_slitid'+StrTrim(slitid,2)+'_exp1.png'

return, comboyfitsmoothnorm

end
