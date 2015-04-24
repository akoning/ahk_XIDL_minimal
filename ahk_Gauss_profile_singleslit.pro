function ahk_Gauss_profile_singleslit,linesflux,slitindex,yfitnorm,result,scifile=scifile,slitid=slitid

;+
; NAME:
;   ahk_Gauss_profile_singleslit
;
; PURPOSE:
;   Fit gaussian(s) + linear baseline across input vector. Normalize to unit area.
;
; CALLING SEQUENCE:
;  ahk_Gauss_profile_singleslit(linesflux,slitindex,yfitnorm,result,slitid=1)
;
; INPUTS:
;
; OPTIONAL INPUTS:
;                
; OUTPUTS:
; yfitfinal (unnormalized)
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
;   3-November-2014 -- Written by Alice Koning
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
ENDREP UNTIL fitresult(0)*fitresult(2) LT 65 ;; or could use: fitresult(0) LT 6*abs(fitresult(3))

CASE 1 OF
	(N_ELEMENTS(params)/mpnterms EQ 1): BEGIN
		print, 'One local maximum found in sky subtracted line profile.'
		expr = 'P[0] + P[1]*X + GAUSS1(X, P[2:4])'
		parinfo = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},5)
		parinfo(0).limited(0) = 1
		parinfo(0).limits(0) = -75
		parinfo(0).limited(1) = 1
		parinfo(0).limits(1) = 75
		start = [params(3), params(4), params(1), params(2) , params(0)*params(2)]
		result = MPFITEXPR(expr, slitindex, linesflux, 0, start, PARINFO=parinfo, WEIGHTS=1D, /QUIET)
		print, 'Gaussian result: ', result

		yfitfinal = result(0)+result(1)*slitindex+gauss1(slitindex, result(2:4))
		;;Output needs to be object profiles without baseline or peaks that are closer than 2*sigma to slitedge.
		;;Profile normalized to have unit area
		yfitnobase = 0.0*slitindex
		IF (result(2) GE abs(result(3))) AND (result(2) LE (slitwidth-abs(result(3)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(2:4))

		;IF (N_ELEMENTS(yfitnobase) GT 0) THEN BEGIN
			;area = int_tabulated(slitindex,yfitnobase)
			;yfitnorm = yfitnobase/abs(area)

			;p = PLOT(slitindex, linesflux, yrange=[0.9*min(linesflux),1.1*max(linesflux)], title=slitid)
			;pgauss = PLOT(slitindex, yfitfinal, 'r', thick=5, /OVERPLOT)
			;p2 = PLOT(slitindex, yfitnorm, 'b', thick=5, title=slitid)
		;ENDIF

	END
	(N_ELEMENTS(params)/mpnterms EQ 2): BEGIN
		print, 'Two local maxima found in sky subtracted line profile.'

		expr = 'P[0] + P[1]*X + GAUSS1(X, P[2:4]) + GAUSS1(X, P[5:7])' 
		parinfo = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},8)
		parinfo(0).limited(0) = 1
		parinfo(0).limits(0) = -75
		parinfo(0).limited(1) = 1
		parinfo(0).limits(1) = 75
		start = [params(8), params(9), params(1), params(2) , params(0)*params(2), params(6), params(7) , params(5)*params(7)]
		result = MPFITEXPR(expr, slitindex, linesflux, 0, start, PARINFO=parinfo, WEIGHTS=1D, /QUIET)
		print, 'Gaussian result: ', result

		yfitfinal = result(0)+result(1)*slitindex+gauss1(slitindex, result(2:4)) + gauss1(slitindex, result(5:7))

		;;Output needs to be object profiles without baseline or peaks that are closer than 2*sigma to slitedge.
		;;Profile normalized to have unit area
		yfitnobase = 0.0*slitindex
		IF (result(2) GE abs(result(3))) AND (result(2) LE (slitwidth-abs(result(3)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(2:4))
		IF (result(5) GE abs(result(6))) AND (result(5) LE (slitwidth-abs(result(6)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(5:7))

		;IF (N_ELEMENTS(yfitnobase) GT 0) THEN BEGIN
			;area = int_tabulated(slitindex,yfitnobase)
			;yfitnorm = yfitnobase/abs(area)

			;p = PLOT(slitindex, linesflux, yrange=[0.9*min(linesflux),1.1*max(linesflux)], title=slitid)
			;pgauss = PLOT(slitindex, yfitfinal, 'r', thick=5, /OVERPLOT)
			;p2 = PLOT(slitindex, yfitnorm, 'b', thick=5, title=slitid)
		;ENDIF

	END
	(N_ELEMENTS(params)/mpnterms EQ 3): BEGIN
		print, 'Three local maxima found in sky subtracted line profile.'

		expr = 'P[0] + P[1]*X + GAUSS1(X, P[2:4]) + GAUSS1(X, P[5:7]) + GAUSS1(X, P[8:10])'
		parinfo = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},11)
		parinfo(0).limited(0) = 1
		parinfo(0).limits(0) = -75
		parinfo(0).limited(1) = 1
		parinfo(0).limits(1) = 75
		start = [params(8), params(9), params(1), params(2) , params(0)*params(2), params(6), params(7) , params(5)*params(7), $
			 params(11), params(12) , params(10)*params(12)]
		result = MPFITEXPR(expr, slitindex, linesflux, 0, start, PARINFO=parinfo, WEIGHTS=1D, /QUIET)
		print, 'Gaussian result: ', result

		yfitfinal = result(0)+result(1)*slitindex+gauss1(slitindex, result(2:4)) + gauss1(slitindex, result(5:7)) $
			+ gauss1(slitindex, result(8:10))

		;;Output needs to be object profiles without baseline or peaks that are closer than 2*sigma to slitedge.
		;;Profile normalized to have unit area
		yfitnobase = 0.0*slitindex
		IF (result(2) GE abs(result(3))) AND (result(2) LE (slitwidth-abs(result(3)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(2:4))
		IF (result(5) GE abs(result(6))) AND (result(5) LE (slitwidth-abs(result(6)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(5:7))
		IF (result(8) GE abs(result(9))) AND (result(8) LE (slitwidth-abs(result(9)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(8:10))

		;IF (N_ELEMENTS(yfitnobase) GT 0) THEN BEGIN
			;area = int_tabulated(slitindex,yfitnobase)
			;yfitnorm = yfitnobase/abs(area)

			;p = PLOT(slitindex, linesflux, yrange=[0.9*min(linesflux),1.1*max(linesflux)], title=slitid)
			;pgauss = PLOT(slitindex, yfitfinal, 'r', thick=5, /OVERPLOT)
			;p2 = PLOT(slitindex, yfitnorm, 'b', thick=5, title=slitid)

			;;Save plot for profileCompareTests
			;p.Save, 'linesflux_slitid'+StrTrim(slitid,2)+'.ps'
			;p2.Save, 'yfitnorm_slitid'+StrTrim(slitid,2)+'.ps'
		;ENDIF

	END
	(N_ELEMENTS(params)/mpnterms EQ 4): BEGIN
		print, 'Four local maxima found in sky subtracted line profile.'

		expr = 'P[0] + P[1]*X + GAUSS1(X, P[2:4]) + GAUSS1(X, P[5:7]) + GAUSS1(X, P[8:10]) + GAUSS1(X, P[11:13])'
		parinfo = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},14)
		parinfo(0).limited(0) = 1
		parinfo(0).limits(0) = -75
		parinfo(0).limited(1) = 1
		parinfo(0).limits(1) = 75
		start = [params(8), params(9), params(1), params(2) , params(0)*params(2), params(6), params(7) , params(5)*params(7), $
			 params(11), params(12) , params(10)*params(12), params(16), params(17) , params(15)*params(17)]
		result = MPFITEXPR(expr, slitindex, linesflux, 0, start, PARINFO=parinfo, WEIGHTS=1D, /QUIET)
		print, 'Gaussian result: ', result

		yfitfinal = result(0)+result(1)*slitindex+gauss1(slitindex, result(2:4)) + gauss1(slitindex, result(5:7)) $
			+ gauss1(slitindex, result(8:10)) + gauss1(slitindex, result(11:13))

		;;Output needs to be object profiles without baseline or peaks that are closer than 2*sigma to slitedge.
		;;Profile normalized to have unit area
		yfitnobase = 0.0*slitindex
		IF (result(2) GE abs(result(3))) AND (result(2) LE (slitwidth-abs(result(3)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(2:4))
		IF (result(5) GE abs(result(6))) AND (result(5) LE (slitwidth-abs(result(6)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(5:7))
		IF (result(8) GE abs(result(9))) AND (result(8) LE (slitwidth-abs(result(9)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(8:10))
		IF (result(11) GE abs(result(12))) AND (result(11) LE (slitwidth-abs(result(12)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(11:13))

		;IF (N_ELEMENTS(yfitnobase) GT 0) THEN BEGIN
			;area = int_tabulated(slitindex,yfitnobase)
			;yfitnorm = yfitnobase/abs(area)

			;p = PLOT(slitindex, linesflux, yrange=[0.9*min(linesflux),1.1*max(linesflux)], title=slitid)
			;pgauss = PLOT(slitindex, yfitfinal, 'r', thick=5, /OVERPLOT)
			;p2 = PLOT(slitindex, yfitnorm, 'b', thick=5, title=slitid)

			;;Save plot for profileCompareTests
			;p.Save, 'linesflux_slitid'+StrTrim(slitid,2)+'.ps'
			;p2.Save, 'yfitnorm_slitid'+StrTrim(slitid,2)+'.ps'
		;ENDIF
	END
	(N_ELEMENTS(params)/mpnterms EQ 5): BEGIN
		print, 'Five local maxima found in sky subtracted line profile.'

		expr = 'P[0] + P[1]*X + GAUSS1(X, P[2:4]) + GAUSS1(X, P[5:7]) + GAUSS1(X, P[8:10]) + GAUSS1(X, P[11:13]) + GAUSS1(X, P[14:16])'
		parinfo = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},17)
		parinfo(0).limited(0) = 1
		parinfo(0).limits(0) = -75
		parinfo(0).limited(1) = 1
		parinfo(0).limits(1) = 75
		start = [params(8), params(9), params(1), params(2) , params(0)*params(2), params(6), params(7) , params(5)*params(7), $
			 params(11), params(12) , params(10)*params(12), params(16), params(17) , params(15)*params(17), $
			params(21), params(22) , params(20)*params(22)]
		result = MPFITEXPR(expr, slitindex, linesflux, 0, start, PARINFO=parinfo, WEIGHTS=1D, /QUIET)
		print, 'Gaussian result: ', result

		yfitfinal = result(0)+result(1)*slitindex+gauss1(slitindex, result(2:4)) + gauss1(slitindex, result(5:7)) $
			+ gauss1(slitindex, result(8:10)) + gauss1(slitindex, result(11:13)) + gauss1(slitindex, result(14:16))

		;;Output needs to be object profiles without baseline or peaks that are closer than 2*sigma to slitedge.
		;;Profile normalized to have unit area
		yfitnobase = 0.0*slitindex
		IF (result(2) GE abs(result(3))) AND (result(2) LE (slitwidth-abs(result(3)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(2:4))
		IF (result(5) GE abs(result(6))) AND (result(5) LE (slitwidth-abs(result(6)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(5:7))
		IF (result(8) GE abs(result(9))) AND (result(8) LE (slitwidth-abs(result(9)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(8:10))
		IF (result(11) GE abs(result(12))) AND (result(11) LE (slitwidth-abs(result(12)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(11:13))
		IF (result(14) GE abs(result(15))) AND (result(14) LE (slitwidth-abs(result(15)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(14:16))

		;IF (N_ELEMENTS(yfitnobase) GT 0) THEN BEGIN
			;area = int_tabulated(slitindex,yfitnobase)
			;yfitnorm = yfitnobase/abs(area)
	
			;p = PLOT(slitindex, linesflux, yrange=[0.9*min(linesflux),1.1*max(linesflux)], title=slitid)
			;pgauss = PLOT(slitindex, yfitfinal, 'r', thick=5, /OVERPLOT)
			;p2 = PLOT(slitindex, yfitnorm, 'b', thick=5, title=slitid)

			;;Save plot for profileCompareTests
			;p.Save, 'linesflux_slitid'+StrTrim(slitid,2)+'.ps'
			;p2.Save, 'yfitnorm_slitid'+StrTrim(slitid,2)+'.ps'
		;ENDIF
	END
	(N_ELEMENTS(params)/mpnterms EQ 6): BEGIN
		print, 'Six local maxima found in sky subtracted line profile.'

		expr = 'P[0] + P[1]*X + GAUSS1(X, P[2:4]) + GAUSS1(X, P[5:7]) + GAUSS1(X, P[8:10]) + GAUSS1(X, P[11:13]) + GAUSS1(X, P[14:16]) + GAUSS1(X, P[17:19])'
		parinfo = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},20)
		parinfo(0).limited(0) = 1
		parinfo(0).limits(0) = -75
		parinfo(0).limited(1) = 1
		parinfo(0).limits(1) = 75
		start = [params(8), params(9), params(1), params(2) , params(0)*params(2), params(6), params(7) , params(5)*params(7), $
			 params(11), params(12) , params(10)*params(12), params(16), params(17) , params(15)*params(17), $
			params(21), params(22) , params(20)*params(22), params(26), params(27) , params(25)*params(27)]
		result = MPFITEXPR(expr, slitindex, linesflux, 0, start, PARINFO=parinfo, WEIGHTS=1D, /QUIET)
		print, 'Gaussian result: ', result

		yfitfinal = result(0)+result(1)*slitindex+gauss1(slitindex, result(2:4)) + gauss1(slitindex, result(5:7)) $
			+ gauss1(slitindex, result(8:10)) + gauss1(slitindex, result(11:13)) + gauss1(slitindex, result(14:16)) $
			+ gauss1(slitindex, result(17:19))

		;;Output needs to be object profiles without baseline or peaks that are closer than 2*sigma to slitedge.
		;;Profile normalized to have unit area
		yfitnobase = 0.0*slitindex
		IF (result(2) GE abs(result(3))) AND (result(2) LE (slitwidth-abs(result(3)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(2:4))
		IF (result(5) GE abs(result(6))) AND (result(5) LE (slitwidth-abs(result(6)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(5:7))
		IF (result(8) GE abs(result(9))) AND (result(8) LE (slitwidth-abs(result(9)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(8:10))
		IF (result(11) GE abs(result(12))) AND (result(11) LE (slitwidth-abs(result(12)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(11:13))
		IF (result(14) GE abs(result(15))) AND (result(14) LE (slitwidth-abs(result(15)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(14:16))
		IF (result(17) GE abs(result(18))) AND (result(17) LE (slitwidth-abs(result(18)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(17:19))

		;IF (N_ELEMENTS(yfitnobase) GT 0) THEN BEGIN
			;area = int_tabulated(slitindex,yfitnobase)
			;yfitnorm = yfitnobase/abs(area)

			;p = PLOT(slitindex, linesflux, yrange=[0.9*min(linesflux),1.1*max(linesflux)], title=slitid)
			;pgauss = PLOT(slitindex, yfitfinal, 'r', thick=5, /OVERPLOT)
			;p2 = PLOT(slitindex, yfitnorm, 'b', thick=5, title=slitid)

			;;Save plot for profileCompareTests
			;p.Save, 'linesflux_slitid'+StrTrim(slitid,2)+'.ps'
			;p2.Save, 'yfitnorm_slitid'+StrTrim(slitid,2)+'.ps'
		;ENDIF
	END
	(N_ELEMENTS(params)/mpnterms EQ 7): BEGIN
		print, 'Seven local maxima found in sky subtracted line profile.'

		expr = 'P[0] + P[1]*X + GAUSS1(X, P[2:4]) + GAUSS1(X, P[5:7]) + GAUSS1(X, P[8:10]) + GAUSS1(X, P[11:13]) + GAUSS1(X, P[14:16]) + GAUSS1(X, P[17:19]) + GAUSS1(X, P[20:22])'
		parinfo = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},23)
		parinfo(0).limited(0) = 1
		parinfo(0).limits(0) = -75
		parinfo(0).limited(1) = 1
		parinfo(0).limits(1) = 75
		start = [params(8), params(9), params(1), params(2) , params(0)*params(2), params(6), params(7) , params(5)*params(7), $
			 params(11), params(12) , params(10)*params(12), params(16), params(17) , params(15)*params(17), $
			params(21), params(22) , params(20)*params(22), params(26), params(27) , params(25)*params(27), $
			params(31), params(32) , params(30)*params(32)]
		result = MPFITEXPR(expr, slitindex, linesflux, 0, start, PARINFO=parinfo, WEIGHTS=1D, /QUIET)
		print, 'Gaussian result: ', result

		yfitfinal = result(0)+result(1)*slitindex+gauss1(slitindex, result(2:4)) + gauss1(slitindex, result(5:7)) $
			+ gauss1(slitindex, result(8:10)) + gauss1(slitindex, result(11:13)) + gauss1(slitindex, result(14:16)) $
			+ gauss1(slitindex, result(17:19)) + gauss1(slitindex, result(20:22))

		;;Output needs to be object profiles without baseline or peaks that are closer than 2*sigma to slitedge.
		;;Profile normalized to have unit area
		yfitnobase = 0.0*slitindex
		IF (result(2) GE abs(result(3))) AND (result(2) LE (slitwidth-abs(result(3)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(2:4))
		IF (result(5) GE abs(result(6))) AND (result(5) LE (slitwidth-abs(result(6)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(5:7))
		IF (result(8) GE abs(result(9))) AND (result(8) LE (slitwidth-abs(result(9)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(8:10))
		IF (result(11) GE abs(result(12))) AND (result(11) LE (slitwidth-abs(result(12)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(11:13))
		IF (result(14) GE abs(result(15))) AND (result(14) LE (slitwidth-abs(result(15)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(14:16))
		IF (result(17) GE abs(result(18))) AND (result(17) LE (slitwidth-abs(result(18)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(17:19))
		IF (result(20) GE abs(result(21))) AND (result(20) LE (slitwidth-abs(result(21)))) THEN yfitnobase = yfitnobase + gauss1(slitindex, result(20:22))

		;IF (N_ELEMENTS(yfitnobase) GT 0) THEN BEGIN
			;area = int_tabulated(slitindex,yfitnobase)
			;yfitnorm = yfitnobase/abs(area)

			;p = PLOT(slitindex, linesflux, yrange=[0.9*min(linesflux),1.1*max(linesflux)], title=slitid)
			;pgauss = PLOT(slitindex, yfitfinal, 'r', thick=5, /OVERPLOT)
			;p2 = PLOT(slitindex, yfitnorm, 'b', thick=5, title=slitid)

			;;Save plot for profileCompareTests
			;p.Save, 'linesflux_slitid'+StrTrim(slitid,2)+'.ps'
			;p2.Save, 'yfitnorm_slitid'+StrTrim(slitid,2)+'.ps'
		;ENDIF
	END
ELSE: BEGIN
	print, 'Error identifying peaks.'
	result = 1.00
	yfitfinal = 1.00
	;;p = PLOT(slitindex, linesflux, yrange=[0.9*min(linesflux),1.1*max(linesflux)], title=slitid)
	;;p2 = PLOT(slitindex, yfitfinal, 'b', thick=5, title=slitid)
END
ENDCASE

IF (N_ELEMENTS(yfitnobase) GT 0) THEN BEGIN
	area = int_tabulated(slitindex,yfitnobase)
	yfitnorm = yfitnobase/abs(area)

	p = PLOT(slitindex, linesflux, yrange=[0.9*min(linesflux),1.1*max(linesflux)], title=slitid)
	pgauss = PLOT(slitindex, yfitfinal, 'r', thick=5, /OVERPLOT)
	p2 = PLOT(slitindex, yfitnorm, 'b', thick=5, title=slitid)

	;;Save plot for profileCompareTests
	p.Save, 'linesflux_'+strmid(scifile,15,8,/reverse_offset)+'_slitid'+StrTrim(slitid,2)+'.png'
	p2.Save, 'yfitnorm_'+strmid(scifile,15,8,/reverse_offset)+'_slitid'+StrTrim(slitid,2)+'.png'

ENDIF ELSE BEGIN
	undefine, yfitnorm
ENDELSE

return, yfitfinal

end
