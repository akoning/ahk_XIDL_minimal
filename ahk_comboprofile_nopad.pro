function ahk_comboprofile_nopad,slitid,slitindex_blue600,yfitnorm_blue600,slitindex_red900,yfitnorm_red900,slitindex_red400,yfitnorm_red400

;+
; NAME:
;   ahk_comboprofile_nopad
;
; PURPOSE:
;   For a particular slit, this code will:
;   1) Calculate normalized slitindex, since blue slitindex does not equal red. This will eventually need to change to be RA/Dec.
;   2) Multiply values from the n grisms for which values are available, at a particular normalized slitindex (RA/Dec)
;   3) Renormalize
;   Code allows for cases when not all 3 profiles are defined.
;
; CALLING SEQUENCE:
;  ahk_comboprofile_nopad(slitid,slitindex_blue600,yfitnorm_blue600,slitindex_red900,yfitnorm_red900,slitindex_red400,yfitnorm_red400)
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
		; Don't expect this to be perfect for all slits until this is changed to RA/Dec, especially for slits on edges that get cut-off
		xnorm_red900 = slitindex_red900/MAX(slitindex_red900)
		xnorm_red400 = slitindex_red400/MAX(slitindex_red400)

		; Make red900 slitindex the reference
		; so interpolate red400 y values, at xnorm_red900 slitindices.
		yinterp_red400 = INTERPOL(yfitnorm_red400,xnorm_red400,xnorm_red900)

		; Now compute the combined profile by multiplying all the normalized profiles together (at same xnorm_blue600 slitindices)
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
		; Don't expect this to be perfect for all slits until this is changed to RA/Dec, especially for slits on edges that get cut-off
		xnorm_blue600 = slitindex_blue600/MAX(slitindex_blue600)
		xnorm_red400 = slitindex_red400/MAX(slitindex_red400)

		; Make blue slitindex the reference (because expect red slits to be cut-off at edges, not blue)
		; so interpolate red400 y values, at xnorm_blue600 slitindices.
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
		; Don't expect this to be perfect for all slits until this is changed to RA/Dec, especially for slits on edges that get cut-off
		xnorm_blue600 = slitindex_blue600/MAX(slitindex_blue600)
		xnorm_red900 = slitindex_red900/MAX(slitindex_red900)

		; Make blue slitindex the reference (because expect red slits to be cut-off at edges, not blue)
		; so interpolate red900 y values, at xnorm_blue600 slitindices.
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
		; Don't expect this to be perfect for all slits until this is changed to RA/Dec, especially for slits on edges that get cut-off
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
                pblue.Save, 'allProfiles_nopad_slitid'+StrTrim(slitid,2)+'_exp2.png'

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

p = PLOT(comboyfitnorm, 'b', thick=5, title=slitid)
p.Save, 'comboyfitnorm_nopad__slitid'+StrTrim(slitid,2)+'_exp2.png'

return, comboyfitnorm

end
