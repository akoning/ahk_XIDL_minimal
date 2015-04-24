function ahk_profile2imarray,slitindex,yfitnorm,slitid=slitid,slitfile=slitfile

;+
; NAME:
;   ahk_profile2imarray
;
; PURPOSE:
;   Take as input a normalized profile distribution (as a function of position on slit; slitindex) and evaluate at normalized x position to produce 4096x4096 profile.
;
; CALLING SEQUENCE:
;  ahk_profile2imarray(yfitnorm,slitid=1,slitfile='slits-lblue2039.fits')
;
; INPUTS:
;
; OPTIONAL INPUTS:
;                
; OUTPUTS:
; 4096x4096 array containing zeros everywhere outside slit, and normalized amplitudes running along slit (following long_slits2x x pos) of fitted Gaussian(s).
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

;;Find normalized x position of slit on detector and apply to yfitnorm
tset_slits = xmrdfits(slitfile,1,silent=(keyword_set(verbose)EQ 0))
ximg = long_slits2x(tset_slits, slitid=slitid)
sximg = size(ximg)

slitprofile = MAKE_ARRAY(sximg(1),sximg(2)) ;;Initialize output profile array to same size as ximg
xfitnorm = slitindex/MAX(slitindex)

FOR row=0,sximg(2)-1 DO BEGIN
	yinterp = INTERPOL(yfitnorm, xfitnorm, ximg[*,row])
	ximgmask = (ximg[*,row] GT 0)
	yinterpfinal = yinterp*ximgmask
	slitprofile[*,row]=yinterpfinal
ENDFOR

return, slitprofile

end
