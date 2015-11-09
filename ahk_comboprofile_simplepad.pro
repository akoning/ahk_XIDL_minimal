function ahk_comboprofile_simplepad,slitid,slitindex_blue600,yfitnorm_blue600, $
slitindex_red900,yfitnorm_red900,slitindex_red400,yfitnorm_red400, $
ccdflagblue=ccdflagblue ,ccdflagred=ccdflagred, igrating=igrating

;+
; NAME:
;   ahk_comboprofile
;
; PURPOSE:
;   **Similar to ahk_comboprofile, but uses additional input of
;   ccdflagblue and ccdflagred instead of just ccdside. Each ccdflag
;   will be 1 if padding should go on the right of the slit, 2 if the
;   padding should go on the left of the slit and 0 if no padding is to be applied.
;
;   Also have additional input of igrating so that the appropriate de-padded comboyfitsmoothnorm
;   profile can be output for the desired grating (igrating 0 for blue600, 1 for red900, 2 for red400)
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
;   24-Jun-2015 -- Edited to handle any number of defined profiles
;-

redplatescale = 0.21746 ;;[arcsec/pix]
blueplatescale = 0.135 ;;[arcsec/pix]


IF N_ELEMENTS(yfitnorm_blue600) NE 0 THEN blue600length = blueplatescale*MAX(slitindex_blue600)
IF N_ELEMENTS(yfitnorm_red900) NE 0 THEN red900length = redplatescale*MAX(slitindex_red900)
IF N_ELEMENTS(yfitnorm_red400) NE 0 THEN red400length = redplatescale*MAX(slitindex_red400)

;;First, take care of the easy cases, when none of the profiles are
;;defined, or only one is defined
CASE 1 OF
	;; None of the three profiles are defined
	(N_ELEMENTS(yfitnorm_blue600) EQ 0 AND N_ELEMENTS(yfitnorm_red900) EQ 0 AND N_ELEMENTS(yfitnorm_red400) EQ 0): BEGIN
		print, 'None of the three profiles are defined'
		comboyfitsmoothnorm = MAKE_ARRAY(10,1, VALUE = 0)
                return, comboyfitsmoothnorm
	END

	;; Only Blue600 defined
	(N_ELEMENTS(yfitnorm_red900) EQ 0 AND N_ELEMENTS(yfitnorm_red400) EQ 0): BEGIN
		print, 'Only Blue600 defined'
		comboyfit = yfitnorm_blue600
		; Subtract the minimum so we have a baseline of zero.
                comboyfit = comboyfit - min(comboyfit)

		; Smooth it!
                boxcar = 2.0      ;[arcsec] this is approx 2x the seeing
                boxcarpix = boxcar/blueplatescale ;[pix]
                comboyfitsmooth = SMOOTH(comboyfit, boxcarpix, /edge_truncate)

		; Renormalize
                area = int_tabulated(slitindex_blue600,comboyfitsmooth)
                comboyfitsmoothnorm = comboyfitsmooth/abs(area)

		; Make the plot and save
                pcombosmooth = PLOT(comboyfitsmoothnorm, 'b', thick=5, title=slitid)
                pcombosmooth.Save, 'comboyfitnorm_simplepad_smooth_nelementMask_slitid'+StrTrim(slitid,2)+'_grating'+StrTrim(igrating,2)+'.png'

                return, comboyfitsmoothnorm
        END
	;; Only Red900 defined
	(N_ELEMENTS(yfitnorm_blue600) EQ 0 AND N_ELEMENTS(yfitnorm_red400) EQ 0): BEGIN
		print, 'Only Red900 defined'
		comboyfit = yfitnorm_red900
		; Subtract the minimum so we have a baseline of zero.
                comboyfit = comboyfit - min(comboyfit)

		; Smooth it!
                boxcar = 2.0      ;[arcsec] this is approx 2x the seeing
                boxcarpix = boxcar/blueplatescale ;[pix]
                comboyfitsmooth = SMOOTH(comboyfit, boxcarpix, /edge_truncate)

		; Renormalize
                area = int_tabulated(slitindex_red900,comboyfitsmooth)
                comboyfitsmoothnorm = comboyfitsmooth/abs(area)

		; Make the plot and save
                pcombosmooth = PLOT(comboyfitsmoothnorm, 'r', thick=5, title=slitid)
                pcombosmooth.Save, 'comboyfitnorm_simplepad_smooth_nelementMask_slitid'+StrTrim(slitid,2)+'_grating'+StrTrim(igrating,2)+'.png'

                return, comboyfitsmoothnorm
        END
	;; Only Red400 defined
	(N_ELEMENTS(yfitnorm_blue600) EQ 0 AND N_ELEMENTS(yfitnorm_red900) EQ 0): BEGIN
		print, 'Only Red400 defined'
		comboyfit = yfitnorm_red400
		; Subtract the minimum so we have a baseline of zero.
                comboyfit = comboyfit - min(comboyfit)

		; Smooth it!
                boxcar = 2.0      ;[arcsec] this is approx 2x the seeing
                boxcarpix = boxcar/blueplatescale ;[pix]
                comboyfitsmooth = SMOOTH(comboyfit, boxcarpix, /edge_truncate)

		; Renormalize
                area = int_tabulated(slitindex_red400,comboyfitsmooth)
                comboyfitsmoothnorm = comboyfitsmooth/abs(area)

		; Make the plot and save
                pcombosmooth = PLOT(comboyfitsmoothnorm, 'k', thick=5, title=slitid)
                pcombosmooth.Save, 'comboyfitnorm_simplepad_smooth_nelementMask_slitid'+StrTrim(slitid,2)+'_grating'+StrTrim(igrating,2)+'.png'

                return, comboyfitsmoothnorm
        END
	ELSE: BEGIN
		print, 'We are dealing with a more complicated case for slit '+StrTrim(slitid,2)+'. Will need to add some padding.'
        END
     ENDCASE

;;Next, deal with case when more than one profile is defined (and padding may be required)
slitindices = ptrarr(3, /allocate_heap)
yfitnorms = ptrarr(3, /allocate_heap)
yfitnorms_unpadded = ptrarr(3, /allocate_heap)
lengths = ptrarr(3, /allocate_heap)

result = [ISA(slitindex_blue600),ISA(slitindex_red900),ISA(slitindex_red400)]
igood = where(result EQ 1)
ngood = TOTAL(result)

FOREACH element, result, key DO BEGIN
   IF element EQ 0 THEN BEGIN
        *slitindices[key] = []
        *yfitnorms[key] = []
        *yfitnorms_unpadded[key] = []
	*lengths[key] = []
	CONTINUE
   ENDIF
   IF key EQ 0 THEN BEGIN
	*slitindices[key] = slitindex_blue600
        *yfitnorms[key] = yfitnorm_blue600
	*yfitnorms_unpadded[key] = yfitnorm_blue600
	*lengths[key] = blue600length
   ENDIF
   IF key EQ 1 THEN BEGIN
        *slitindices[key] = slitindex_red900
        *yfitnorms[key] = yfitnorm_red900
        *yfitnorms_unpadded[key] = yfitnorm_red900
        *lengths[key] = red900length
   ENDIF
   IF key EQ 2 THEN BEGIN
        *slitindices[key] = slitindex_red400
        *yfitnorms[key] = yfitnorm_red400
        *yfitnorms_unpadded[key] = yfitnorm_red400
        *lengths[key] = red400length
   ENDIF
ENDFOREACH

flags = [ccdflagblue, ccdflagred, ccdflagred]
scales= [blueplatescale, redplatescale, redplatescale]         

;;If both red profiles are defined, then we need to check they are close to the same length.
;;If not, the shorter one will be padded with zeroes to exactly match the length of the longer.
;;Problem with this plan is that I don't know which side to pad it on!
;;Will leave for now and come back to later if causing problems. (Will need to depad too!!)
IF where(igood EQ 1) NE -1 AND where(igood EQ 2) NE -1 THEN BEGIN
	IF (*lengths[1]-*lengths[2])/*lengths[1] GT 0.1 OR (*lengths[2]-*lengths[1])/*lengths[2] GT 0.1 THEN BEGIN
   		print, 'WARNING!! Red900 slit length different from Red400 by more than 10 percent for slit '+StrTrim(slitid,2)
	ENDIF
ENDIF

;;Apply the appropriate padding for the various possible scenarios of defined profiles.
;;Save the padding appropriate for the given igrating, so we can de-pad it at the end.
depad = [] ;;Amount of padding
depad_flag = [] ;;1 if padding on right side, 2 if padding on left side

IF ngood EQ 2 THEN BEGIN
	print, 'Two profiles are defined. Will pad each according to ccdflagblue/ccdflagred, unless it is the two red profiles defined.'
	
	IF where(igood EQ 0) EQ -1 THEN BEGIN
		print, 'No blue profile defined. No padding needed.'
	ENDIF ELSE BEGIN
		print, 'One blue and one red profile defined.' ;;igood[0] is the blue profile's index, igood[-1] is the red profile's index.
		CASE 1 OF
		        ;Blue slit needs padding
        		(ccdflagblue NE 0 AND ccdflagred EQ 0): BEGIN
                		print, 'Padding blue slit '+StrTrim(slitid,2)

                		;Pad blue to length of the red slit
                		bluereqpix=FIX(*lengths[igood[-1]]/blueplatescale)
                		yfitnorm_blue600_padded=MAKE_ARRAY(bluereqpix, VALUE=0.0)
                		diff = (bluereqpix-size(*slitindices[igood[0]],/N_ELE))>0.0
                		*slitindices[igood[0]]=findgen(bluereqpix)
				IF igrating EQ 0 THEN depad = diff

                		IF ccdflagblue EQ 2 THEN BEGIN
                   			;;Need to pad on left side, so anchor blue and red at 1
                   			FOR i=diff,(bluereqpix-1) DO yfitnorm_blue600_padded[i]=(*yfitnorms_unpadded[igood[0]])[i-diff]
                   			*yfitnorms[igood[0]] = yfitnorm_blue600_padded
					IF igrating EQ 0 THEN depad_flag = 2
                		ENDIF ELSE IF ccdflagblue EQ 1 THEN BEGIN
                    			;;Need to pad on right side, so anchor blue and red at 0
                   			FOR i=0,(bluereqpix-diff-1) DO yfitnorm_blue600_padded[i]=(*yfitnorms_unpadded[igood[0]])[i]
                   			*yfitnorms[igood[0]] = yfitnorm_blue600_padded
                                        IF igrating EQ 0 THEN depad_flag = 1
                		ENDIF ELSE BEGIN
                  			print, 'Do not recognize ccdflagblue. Exiting.'
                   			return, 0
                		ENDELSE

        		END
        		;Red slit needs padding
        		(ccdflagblue EQ 0 AND ccdflagred NE 0): BEGIN
                		print, 'Padding red slit '+StrTrim(slitid,2)
		
                		;Pad red slit to length of the blue slit
                		redreqpix=FIX(*lengths[igood[0]]/redplatescale)
                		yfitnorm_red_padded=MAKE_ARRAY(redreqpix, VALUE=0.0)
                		diff = (redreqpix-size(*slitindices[igood[-1]],/N_ELE))>0.0
                		*slitindices[igood[-1]]=findgen(redreqpix)
                                IF igrating EQ igood[-1] THEN depad = diff

                		IF ccdflagred EQ 2 THEN BEGIN
                   			;;Need to pad on left side, so anchor blue and red at 1
                   			FOR i=diff,(redreqpix-1) DO yfitnorm_red_padded[i]=(*yfitnorms_unpadded[igood[-1]])[i-diff]
                   			*yfitnorms[igood[-1]] = yfitnorm_red_padded
                                        IF igrating EQ igood[-1] THEN depad_flag = 2
                		ENDIF ELSE IF ccdflagred EQ 1 THEN BEGIN
                    			;;Need to pad on right side, so anchor blue and red at 0
                   			FOR i=0,(redreqpix-diff-1) DO yfitnorm_red_padded[i]=(*yfitnorms_unpadded[igood[-1]])[i]
                   			*yfitnorms[igood[-1]] = yfitnorm_red_padded
                                        IF igrating EQ igood[-1] THEN depad_flag = 1
                		ENDIF ELSE BEGIN
                   			print, 'Do not recognize ccdflagred. Exiting.'
                			return, 0
                		ENDELSE

        		END

                        ;Neither of the slits need padding
                        (ccdflagblue EQ 0  AND ccdflagred EQ 0): BEGIN
                                print, 'No padding for slit '+StrTrim(slitid,2)
                        END

        		ELSE: BEGIN
                		print, 'Do not recognize ccdflags for slit '+StrTrim(slitid,2)
        		END
		ENDCASE
	ENDELSE

ENDIF ELSE IF ngood EQ 3 THEN BEGIN
        print, 'All three profiles are defined. Will pad each according to ccdflagblue/ccdflagred.' ;Blue is index 0, red900 is index 1, red400 is index 2.

                CASE 1 OF
                        ;Blue slit needs padding
                        (ccdflagblue NE 0 AND ccdflagred EQ 0): BEGIN
                                print, 'Padding blue slit '+ StrTrim(slitid,2)

                                ;Pad blue to the average length of the red slits
                                bluereqpix=FIX(((*lengths[1]+*lengths[2])/2)/blueplatescale)
                                yfitnorm_blue600_padded=MAKE_ARRAY(bluereqpix, VALUE=0.0)
                                diff = (bluereqpix-size(*slitindices[0],/N_ELE))>0.0
                                *slitindices[0]=findgen(bluereqpix)
                                IF igrating EQ 0 THEN depad = diff

                                IF ccdflagblue EQ 2 THEN BEGIN
                                        ;;Need to pad on left side, so anchor blue and red at 1
                                        FOR i=diff,(bluereqpix-1) DO yfitnorm_blue600_padded[i]=(*yfitnorms_unpadded[0])[i-diff]
                                        *yfitnorms[0] = yfitnorm_blue600_padded
                                        IF igrating EQ 0 THEN depad_flag = 2
                                ENDIF ELSE IF ccdflagblue EQ 1 THEN BEGIN
                                        ;;Need to pad on right side, so anchor blue and red at 0
                                        FOR i=0,(bluereqpix-diff-1) DO yfitnorm_blue600_padded[i]=(*yfitnorms_unpadded[0])[i]
                                        *yfitnorms[0] = yfitnorm_blue600_padded
                                        IF igrating EQ 0 THEN depad_flag = 1
                                ENDIF ELSE BEGIN
                                        print, 'Do not recognize ccdflagblue. Exiting.'
                                        return, 0
                                ENDELSE

                        END
                        ;Red slits need padding
                        (ccdflagblue EQ 0 AND ccdflagred NE 0): BEGIN
                                print, 'Padding red slit '+StrTrim(slitid,2)

                                ;Pad red slit to length of the blue slit
                                redreqpix=FIX(*lengths[igood[0]]/redplatescale)
                                yfitnorm_red900_padded=MAKE_ARRAY(redreqpix, VALUE=0.0)
                                yfitnorm_red400_padded=MAKE_ARRAY(redreqpix, VALUE=0.0)
                                diffred900 = (redreqpix-size(*slitindices[1],/N_ELE))>0.0
                                diffred400 = (redreqpix-size(*slitindices[2],/N_ELE))>0.0
                                *slitindices[1]=findgen(redreqpix)
                                *slitindices[2]=findgen(redreqpix)
                                IF igrating EQ 1 THEN depad = diffred900
                                IF igrating EQ 2 THEN depad = diffred400


                                IF ccdflagred EQ 2 THEN BEGIN
                                        ;;Need to pad on left side, so anchor blue and red at 1
                                        FOR i=diffred900,(redreqpix-1) DO yfitnorm_red900_padded[i]=(*yfitnorms_unpadded[1])[i-diffred900]
                                        *yfitnorms[1] = yfitnorm_red900_padded
                                        IF igrating EQ 1 THEN depad_flag = 2
                                        FOR i=diffred400,(redreqpix-1) DO yfitnorm_red400_padded[i]=(*yfitnorms_unpadded[2])[i-diffred400]
                                        *yfitnorms[2] = yfitnorm_red400_padded
                                        IF igrating EQ 2 THEN depad_flag = 2
                                ENDIF ELSE IF ccdflagred EQ 1 THEN BEGIN
                                        ;;Need to pad on right side, so anchor blue and red at 0
                                        FOR i=0,(redreqpix-diffred900-1) DO yfitnorm_red900_padded[i]=(*yfitnorms_unpadded[1])[i]
                                        *yfitnorms[1] = yfitnorm_red900_padded
                                        IF igrating EQ 1 THEN depad_flag = 1
                                        FOR i=0,(redreqpix-diffred400-1) DO yfitnorm_red400_padded[i]=(*yfitnorms_unpadded[2])[i]
                                        *yfitnorms[2] = yfitnorm_red400_padded
                                        IF igrating EQ 2 THEN depad_flag = 1
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
ENDIF ELSE print, 'Something went wrong.'


;;Now we're ready to compute the comboprofile and renormalize!

xnorm = ptrarr(3, /allocate_heap)
yinterp = ptrarr(3, /allocate_heap)

FOREACH i, igood DO BEGIN
	*xnorm[i] = *slitindices[i]/MAX(*slitindices[i])
ENDFOREACH

;;If blue profile is defined: Interpolate red900 and red400 y values at blue600 xnorm slit indices.
;;Otherwise interpolate red400 y values at red900 xnorm slit indices
IF where(igood EQ 0) NE -1 THEN BEGIN
	FOREACH i, igood DO *yinterp[i] = INTERPOL(*yfitnorms[i],*xnorm[i],*xnorm[0])
	mincrit = 1./MAX(*slitindices[0])
	interp_flag = 1
ENDIF ELSE BEGIN
        FOREACH i, igood DO *yinterp[i] = INTERPOL(*yfitnorms[i],*xnorm[i],*xnorm[1])
        mincrit = 1./MAX(*slitindices[1])
        interp_flag = 2
ENDELSE

IF *yinterp[0] NE !NULL THEN BEGIN
   pall1 = PLOT(*yinterp[0], 'b', title=slitid)
   IF *yinterp[1] NE !NULL THEN pall2 =  PLOT(*yinterp[1], 'r', /OVERPLOT)
   IF *yinterp[2] NE !NULL THEN pall3 =  PLOT(*yinterp[2], 'k', /OVERPLOT)
   pall1.Save, 'allProfiles_simplepad_slitid'+StrTrim(slitid,2)+'_grating'+StrTrim(igrating,2)+'.png'
ENDIF ELSE IF *yinterp[1] NE !NULL THEN BEGIN
   pall2 =  PLOT(*yinterp[1], 'r', title=slitid)
   IF *yinterp[2] NE !NULL THEN pall3 =  PLOT(*yinterp[2], 'k', /OVERPLOT)
   pall2.Save, 'allProfiles_simplepad_slitid'+StrTrim(slitid,2)+'_grating'+StrTrim(igrating,2)+'.png'
ENDIF ELSE IF *yinterp[2] NE !NULL THEN BEGIN
    pall3 =  PLOT(*yinterp[2], 'k', title=slitid)
    pall3.Save, 'allProfiles_simplepad_slitid'+StrTrim(slitid,2)+'_grating'+StrTrim(igrating,2)+'.png'
ENDIF ELSE print, 'No plot available for slit '+StrTrim(slitid,2)+' grating '+StrTrim(igrating,2)

;;Force profiles that will be used to compute comboprofile to be
;;strictly positive using the previously defined mincrit
;;(no physical meaning behind mincrit being defined from the max slitindex,
;;just a convenient way to get a very small number)
print, 'Min criterium:'+strtrim(mincrit,2)
FOREACH i, igood DO BEGIN
	*yinterp[i] = *yinterp[i] > mincrit
ENDFOREACH

; Now compute the combined profile by multiplying all the normalized profiles together (at same xnorm_blue600 slitindices)
; Divide by n-th root of the number of profiles used to keep correct width
comboyfit = []
FOREACH i, igood DO BEGIN
	IF comboyfit EQ !NULL THEN comboyfit = *yinterp[i] $
	ELSE comboyfit *= *yinterp[i] 
ENDFOREACH  
comboyfit=comboyfit^(1./ngood)

; Subtract the minimum so we have a baseline of zero.
comboyfit = comboyfit - min(comboyfit)

; Smooth it!
boxcar = 2.0 ;[arcsec] this is approx 2x the seeing
boxcarpix = boxcar/blueplatescale ;[pix]
comboyfitsmooth = SMOOTH(comboyfit, boxcarpix, /edge_truncate)

; De-pad if necessary. First need to re-interpolate at original xnorm positions for red400 and sometimes red900.
comboyfitsmooth_interp = comboyfitsmooth
IF interp_flag EQ 1 THEN BEGIN
	comboyfitsmooth = INTERPOL(comboyfitsmooth_interp,*xnorm[0],*xnorm[igrating])
ENDIF ELSE IF interp_flag EQ 2 THEN BEGIN
        comboyfitsmooth = INTERPOL(comboyfitsmooth_interp,*xnorm[1],*xnorm[igrating])
ENDIF

IF depad_flag EQ !NULL THEN comboyfitsmooth_depad = comboyfitsmooth ELSE BEGIN
	IF depad_flag EQ 1 THEN comboyfitsmooth_depad = comboyfitsmooth[0:(-depad-1)] $
	ELSE IF depad_flag EQ 2 THEN comboyfitsmooth_depad = comboyfitsmooth[depad:-1] $
	ELSE print, 'depad_flag unknown.' 
ENDELSE

; Renormalize
slitindex_depad = findgen(size(comboyfitsmooth_depad,/N_ELE))
area = int_tabulated(slitindex_depad,comboyfitsmooth_depad)
comboyfitsmoothnorm_depad = comboyfitsmooth_depad/abs(area)

; Make the plot and save
pcombosmooth = PLOT(comboyfitsmoothnorm_depad, 'g', thick=5, title=slitid)
pcombosmooth.Save, 'comboyfitsmoothnorm_depad_simplepad_slitid'+StrTrim(slitid,2)+'_grating'+StrTrim(igrating,2)+'.png'

return, comboyfitsmoothnorm_depad

end


;Check that no one profile is much different than the other two
;Not working yet! Difficult to choose criteria that will always work.
;Also, requires I make an assumption about which profile is "right" and
;which is "wrong". Can't reasonably do this as long as sky subtraction
;is as bad as it is.

;Turn into keyboard entry after plotting all 3?? This would make this
;step a TON more work when reducing entire data set :(

;IF ((max(yfitnorm_blue600) GT 1.8*max(yfitnorm_red900) AND 1.8*max(yfitnorm_red400)) $
;    OR (max(yfitnorm_blue600) LT 0.4*max(yfitnorm_red900) AND 0.4*max(yfitnorm_red400))) THEN BEGIN
;   print, 'Blue600 way different than others. Setting to 1.0'
;   replicate_inplace, yfitnorm_blue600, 1.0
;ENDIF
;IF ((max(yfitnorm_red900) GT 1.8*max(yfitnorm_blue600) AND 1.8*max(yfitnorm_red400)) $
;    OR (max(yfitnorm_red900) LT 0.4*max(yfitnorm_blue600) AND  0.4*max(yfitnorm_red400))) THEN BEGIN
;   print, 'Red900 way different than others. Setting to 1.0'
;   replicate_inplace, yfitnorm_red900, 1.0
;ENDIF
;IF ((max(yfitnorm_red400) GT 1.8*max(yfitnorm_blue600) AND 1.8*max(yfitnorm_red900)) $
;    OR (max(yfitnorm_red400) LT 0.4*max(yfitnorm_blue600) AND 0.4*max(yfitnorm_red900))) THEN BEGIN
;   print, 'Red400 way different than others. Setting to 1.0'
;   replicate_inplace, yfitnorm_red400, 1.0
;ENDIF

