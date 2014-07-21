function ahk_profile,linesflux,yfitnorm,result,slitid=slitid,slitfile=slitfile

;+
; NAME:
;   ahk_profile
;
; PURPOSE:
;   Fit gaussian(s) + linear baseline across input vector. Normalize to unit area and evaluate as a function of normalized x position to produce 4096x4096 profile.
;
; CALLING SEQUENCE:
;  ahk_profile(linesflux,slitid=1,slitfile='slits-lblue2039.fits')
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
;   26-June-2014 -- Written by Alice Koning
;-


;;Fit Gaussian using mpfitpeak and subtract. Repeat until only noise remains. Redo fit using all found mpfitpeak as starting point for mpfitexpr.
params = [] ;; Declare empty array to put params found from mpfitpeak
linesflux_rem = linesflux ;; Copy linesflux into variable which will have each new fit subtracted from it
slitindex = findgen(N_ELEMENTS(linesflux))
slitwidth = max(slitindex)
print, 'slitwidth: ', slitwidth
;print, 'linesflux: ', linesflux
;print, 'index: ', slitindex

counter=0
mpnterms=5
REPEAT BEGIN
	counter = counter+1
   	yfit = mpfitpeak(slitindex, linesflux_rem, fitresult, NTERMS=mpnterms, /NAN, /POSITIVE)
	;; Do not include extremely  broad Gaussians
	IF counter EQ 1 THEN params = [params,fitresult] ELSE BEGIN
		IF ((abs(fitresult[2]) LE (slitwidth*0.5)) AND (fitresult(0)*fitresult(2) GE 100)) THEN params = [params,fitresult]
	ENDELSE
	;testplot1 = plot(slitindex, linesflux_rem)
	;testplot2 = plot(slitindex, yfit, /OVERPLOT)
	linesflux_rem=linesflux_rem-yfit
	print, 'mpfitpeak fitresult: ', fitresult
	print, 'Approx area for mpfitpeak', fitresult(0)*fitresult(2)
ENDREP UNTIL fitresult(0)*fitresult(2) LT 50 ;; or could use: fitresult(0) LT 6*abs(fitresult(3))

CASE 1 OF
	(N_ELEMENTS(params)/mpnterms EQ 1): BEGIN
		print, 'One local maximum found in sky subtracted line profile.'
		expr = 'P[0] + P[1]*X + GAUSS1(X, P[2:4])'
		start = [params(3),params(4), params(1), params(2) , params(0)*params(2)]
		result = MPFITEXPR(expr, slitindex, linesflux, 0, start, /WEIGHTS, /QUIET)
		print, 'Gaussian result: ', result

		yfitfinal = result(0)+result(1)*slitindex+gauss1(slitindex, result(2:4))
		;;Output needs to be object profiles without baseline or peaks that are closer than 2*sigma to slitedge.
		;;Profile normalized to have unit area
		yfitnobase = 0
		IF (result(2) GE abs(2*result(3))) AND (result(2) LE (slitwidth-abs(2*result(3)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(2:4))

		area = int_tabulated(slitindex,yfitnobase)
		yfitnorm = yfitnobase/abs(area)

		p = PLOT(slitindex, linesflux, yrange=[min(linesflux)-5,max(linesflux)+5])
		pgauss = PLOT(slitindex, yfitfinal, 'r', thick=5, /OVERPLOT)
		p2 = PLOT(slitindex, yfitnorm, 'b', thick=5)

	END
	(N_ELEMENTS(params)/mpnterms EQ 2): BEGIN
		print, 'Two local maxima found in sky subtracted line profile.'

		expr = 'P[0] + P[1]*X + GAUSS1(X, P[2:4]) + GAUSS1(X, P[5:7])' 
		start = [params(3),params(4), params(1), params(2) , params(0)*params(2), params(6), params(7) , params(5)*params(7)]
		result = MPFITEXPR(expr, slitindex, linesflux, 0, start, /WEIGHTS, /QUIET)
		print, 'Gaussian result: ', result

		yfitfinal = result(0)+result(1)*slitindex+gauss1(slitindex, result(2:4)) + gauss1(slitindex, result(5:7))

		;;Output needs to be object profiles without baseline or peaks that are closer than 2*sigma to slitedge.
		;;Profile normalized to have unit area
		yfitnobase = 0
		IF (result(2) GE abs(2*result(3))) AND (result(2) LE (slitwidth-abs(2*result(3)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(2:4))
		IF (result(5) GE abs(2*result(6))) AND (result(5) LE (slitwidth-abs(2*result(6)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(5:7))

		area = int_tabulated(slitindex,yfitnobase)
		yfitnorm = yfitnobase/abs(area)

		p = PLOT(slitindex, linesflux, yrange=[min(linesflux)-5,max(linesflux)+5])
		pgauss = PLOT(slitindex, yfitfinal, 'r', thick=5, /OVERPLOT)
		p2 = PLOT(slitindex, yfitnorm, 'b', thick=5)

	END
	(N_ELEMENTS(params)/mpnterms EQ 3): BEGIN
		print, 'Three local maxima found in sky subtracted line profile.'

		expr = 'P[0] + P[1]*X + GAUSS1(X, P[2:4]) + GAUSS1(X, P[5:7]) + GAUSS1(X, P[8:10])'
		start = [params(3),params(4), params(1), params(2) , params(0)*params(2), params(6), params(7) , params(5)*params(7), $
			 params(11), params(12) , params(10)*params(12)]
		result = MPFITEXPR(expr, slitindex, linesflux, 0, start, /WEIGHTS, /QUIET)
		print, 'Gaussian result: ', result

		yfitfinal = result(0)+result(1)*slitindex+gauss1(slitindex, result(2:4)) + gauss1(slitindex, result(5:7)) $
			+ gauss1(slitindex, result(8:10))

		;;Output needs to be object profiles without baseline or peaks that are closer than 2*sigma to slitedge.
		;;Profile normalized to have unit area
		yfitnobase = 0
		IF (result(2) GE abs(2*result(3))) AND (result(2) LE (slitwidth-abs(2*result(3)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(2:4))
		IF (result(5) GE abs(2*result(6))) AND (result(5) LE (slitwidth-abs(2*result(6)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(5:7))
		IF (result(8) GE abs(2*result(9))) AND (result(8) LE (slitwidth-abs(2*result(9)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(8:10))

		area = int_tabulated(slitindex,yfitnobase)
		yfitnorm = yfitnobase/abs(area)

		p = PLOT(slitindex, linesflux, yrange=[min(linesflux)-5,max(linesflux)+5])
		pgauss = PLOT(slitindex, yfitfinal, 'r', thick=5, /OVERPLOT)
		p2 = PLOT(slitindex, yfitnorm, 'b', thick=5)

	END
	(N_ELEMENTS(params)/mpnterms EQ 4): BEGIN
		print, 'Four local maxima found in sky subtracted line profile.'

		expr = 'P[0] + P[1]*X + GAUSS1(X, P[2:4]) + GAUSS1(X, P[5:7]) + GAUSS1(X, P[8:10]) + GAUSS1(X, P[11:13])'
		start = [params(3),params(4), params(1), params(2) , params(0)*params(2), params(6), params(7) , params(5)*params(7), $
			 params(11), params(12) , params(10)*params(12), params(16), params(17) , params(15)*params(17)]
		result = MPFITEXPR(expr, slitindex, linesflux, 0, start, /WEIGHTS, /QUIET)
		print, 'Gaussian result: ', result

		yfitfinal = result(0)+result(1)*slitindex+gauss1(slitindex, result(2:4)) + gauss1(slitindex, result(5:7)) $
			+ gauss1(slitindex, result(8:10)) + gauss1(slitindex, result(11:13))

		;;Output needs to be object profiles without baseline or peaks that are closer than 2*sigma to slitedge.
		;;Profile normalized to have unit area
		yfitnobase = 0
		IF (result(2) GE abs(2*result(3))) AND (result(2) LE (slitwidth-abs(2*result(3)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(2:4))
		IF (result(5) GE abs(2*result(6))) AND (result(5) LE (slitwidth-abs(2*result(6)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(5:7))
		IF (result(8) GE abs(2*result(9))) AND (result(8) LE (slitwidth-abs(2*result(9)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(8:10))
		IF (result(11) GE abs(2*result(12))) AND (result(11) LE (slitwidth-abs(2*result(12)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(11:13))

		area = int_tabulated(slitindex,yfitnobase)
		yfitnorm = yfitnobase/abs(area)


		p = PLOT(slitindex, linesflux, yrange=[min(linesflux)-5,max(linesflux)+5])
		pgauss = PLOT(slitindex, yfitfinal, 'r', thick=5, /OVERPLOT)
		p2 = PLOT(slitindex, yfitnorm-0.01, yrange=[0,max(yfitnorm)])

	END
	(N_ELEMENTS(params)/mpnterms EQ 5): BEGIN
		print, 'Five local maxima found in sky subtracted line profile.'

		expr = 'P[0] + P[1]*X + GAUSS1(X, P[2:4]) + GAUSS1(X, P[5:7]) + GAUSS1(X, P[8:10]) + GAUSS1(X, P[11:13]) + GAUSS1(X, P[14:16])'
		start = [params(3),params(4), params(1), params(2) , params(0)*params(2), params(6), params(7) , params(5)*params(7), $
			 params(11), params(12) , params(10)*params(12), params(16), params(17) , params(15)*params(17), $
			params(21), params(22) , params(20)*params(22)]
		result = MPFITEXPR(expr, slitindex, linesflux, 0, start, /WEIGHTS, /QUIET)
		print, 'Gaussian result: ', result

		yfitfinal = result(0)+result(1)*slitindex+gauss1(slitindex, result(2:4)) + gauss1(slitindex, result(5:7)) $
			+ gauss1(slitindex, result(8:10)) + gauss1(slitindex, result(11:13)) + gauss1(slitindex, result(14:16))

		;;Output needs to be object profiles without baseline or peaks that are closer than 2*sigma to slitedge.
		;;Profile normalized to have unit area
		yfitnobase = 0
		IF (result(2) GE abs(2*result(3))) AND (result(2) LE (slitwidth-abs(2*result(3)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(2:4))
		IF (result(5) GE abs(2*result(6))) AND (result(5) LE (slitwidth-abs(2*result(6)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(5:7))
		IF (result(8) GE abs(2*result(9))) AND (result(8) LE (slitwidth-abs(2*result(9)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(8:10))
		IF (result(11) GE abs(2*result(12))) AND (result(11) LE (slitwidth-abs(2*result(12)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(11:13))
		IF (result(14) GE abs(2*result(15))) AND (result(14) LE (slitwidth-abs(2*result(15)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(14:16))

		area = int_tabulated(slitindex,yfitnobase)
		yfitnorm = yfitnobase/abs(area)
	
		p = PLOT(slitindex, linesflux, yrange=[min(linesflux)-5,max(linesflux)+5])
		pgauss = PLOT(slitindex, yfitfinal, 'r', thick=5, /OVERPLOT)
		p2 = PLOT(slitindex, yfitnorm, 'b', thick=5)

	END
	(N_ELEMENTS(params)/mpnterms EQ 6): BEGIN
		print, 'Six local maxima found in sky subtracted line profile.'

		expr = 'P[0] + P[1]*X + GAUSS1(X, P[2:4]) + GAUSS1(X, P[5:7]) + GAUSS1(X, P[8:10]) + GAUSS1(X, P[11:13]) + GAUSS1(X, P[14:16]) + GAUSS1(X, P[17:19])'
		start = [params(3),params(4), params(1), params(2) , params(0)*params(2), params(6), params(7) , params(5)*params(7), $
			 params(11), params(12) , params(10)*params(12), params(16), params(17) , params(15)*params(17), $
			params(21), params(22) , params(20)*params(22), params(26), params(27) , params(25)*params(27)]
		result = MPFITEXPR(expr, slitindex, linesflux, 0, start, /WEIGHTS, /QUIET)
		print, 'Gaussian result: ', result

		yfitfinal = result(0)+result(1)*slitindex+gauss1(slitindex, result(2:4)) + gauss1(slitindex, result(5:7)) $
			+ gauss1(slitindex, result(8:10)) + gauss1(slitindex, result(11:13)) + gauss1(slitindex, result(14:16)) $
			+ gauss1(slitindex, result(17:19))

		;;Output needs to be object profiles without baseline or peaks that are closer than 2*sigma to slitedge.
		;;Profile normalized to have unit area
		yfitnobase = 0
		IF (result(2) GE abs(2*result(3))) AND (result(2) LE (slitwidth-abs(2*result(3)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(2:4))
		IF (result(5) GE abs(2*result(6))) AND (result(5) LE (slitwidth-abs(2*result(6)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(5:7))
		IF (result(8) GE abs(2*result(9))) AND (result(8) LE (slitwidth-abs(2*result(9)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(8:10))
		IF (result(11) GE abs(2*result(12))) AND (result(11) LE (slitwidth-abs(2*result(12)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(11:13))
		IF (result(14) GE abs(2*result(15))) AND (result(14) LE (slitwidth-abs(2*result(15)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(14:16))
		IF (result(17) GE abs(2*result(18))) AND (result(17) LE (slitwidth-abs(2*result(18)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(17:19))

		area = int_tabulated(slitindex,yfitnobase)
		yfitnorm = yfitnobase/abs(area)

		p = PLOT(slitindex, linesflux, yrange=[min(linesflux)-5,max(linesflux)+5])
		pgauss = PLOT(slitindex, yfitfinal, 'r', thick=5, /OVERPLOT)
		p2 = PLOT(slitindex, yfitnorm, 'b', thick=5)

	END
	(N_ELEMENTS(params)/mpnterms EQ 7): BEGIN
		print, 'Seven local maxima found in sky subtracted line profile.'

		expr = 'P[0] + P[1]*X + GAUSS1(X, P[2:4]) + GAUSS1(X, P[5:7]) + GAUSS1(X, P[8:10]) + GAUSS1(X, P[11:13]) + GAUSS1(X, P[14:16]) + GAUSS1(X, P[17:19]) + GAUSS1(X, P[20:22])'
		start = [params(3),params(4), params(1), params(2) , params(0)*params(2), params(6), params(7) , params(5)*params(7), $
			 params(11), params(12) , params(10)*params(12), params(16), params(17) , params(15)*params(17), $
			params(21), params(22) , params(20)*params(22), params(26), params(27) , params(25)*params(27), $
			params(31), params(32) , params(30)*params(32)]
		result = MPFITEXPR(expr, slitindex, linesflux, 0, start, /WEIGHTS, /QUIET)
		print, 'Gaussian result: ', result

		yfitfinal = result(0)+result(1)*slitindex+gauss1(slitindex, result(2:4)) + gauss1(slitindex, result(5:7)) $
			+ gauss1(slitindex, result(8:10)) + gauss1(slitindex, result(11:13)) + gauss1(slitindex, result(14:16)) $
			+ gauss1(slitindex, result(17:19)) + gauss1(slitindex, result(20:22))

		;;Output needs to be object profiles without baseline or peaks that are closer than 2*sigma to slitedge.
		;;Profile normalized to have unit area
		yfitnobase = 0
		IF (result(2) GE abs(2*result(3))) AND (result(2) LE (slitwidth-abs(2*result(3)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(2:4))
		IF (result(5) GE abs(2*result(6))) AND (result(5) LE (slitwidth-abs(2*result(6)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(5:7))
		IF (result(8) GE abs(2*result(9))) AND (result(8) LE (slitwidth-abs(2*result(9)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(8:10))
		IF (result(11) GE abs(2*result(12))) AND (result(11) LE (slitwidth-abs(2*result(12)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(11:13))
		IF (result(14) GE abs(2*result(15))) AND (result(14) LE (slitwidth-abs(2*result(15)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(14:16))
		IF (result(17) GE abs(2*result(18))) AND (result(17) LE (slitwidth-abs(2*result(18)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(17:19))
		IF (result(20) GE abs(2*result(21))) AND (result(20) LE (slitwidth-abs(2*result(21)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(20:22))

		area = int_tabulated(slitindex,yfitnobase)
		yfitnorm = yfitnobase/abs(area)

		p = PLOT(slitindex, linesflux, yrange=[min(linesflux)-5,max(linesflux)+5])
		pgauss = PLOT(slitindex, yfitfinal, 'r', thick=5, /OVERPLOT)
		p2 = PLOT(slitindex, yfitnorm, 'b', thick=5)

	END
ELSE: BEGIN
	print, 'Error identifying peaks.'
	result = 0.00
	yfitfinal = 0.00
	yfitnobase = 0.00
	yfitnorm = 0.00
	p = PLOT(slitindex, linesflux, yrange=[min(linesflux)-5,max(linesflux)+5])
	p2 = PLOT(slitindex, yfitnorm, 'b', thick=5)
END
ENDCASE

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
