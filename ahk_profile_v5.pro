pro ahk_profile_v5,scifile,wavefile=wavefile,slitfile=slitfile

;+
; NAME:
;   ahk_profile_v5
;
; PURPOSE:
;   Loop over slits and determine spatial profile of each using mpfitpeak with nterms=5.
;
; CALLING SEQUENCE:
;
; INPUTS:
;    scifile  -- science image (multi extension FITS, processed 2D image of data is 0th extension)
;    wavefile -- FITS image of wavelength solution (wave-xxx.fits output by long_reduce)
;    slitfile -- binary FITS table (or FITS image??) containing the parameters which describe the slit edges (slits-xxx.fits output by long_reduce)
;
; OPTIONAL INPUTS:
;                
; OUTPUTS:
;  Returns the line profile in counts along the slit. 
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
;   12-June-2014 -- Written by Alice Koning
;-

;; Check for each input file
if  KEYWORD_SET(scifile) then print,'Science image file: ' + scifile $
else begin
	print, 'No science file found. Leaving ahk_profile.'
	return
endelse

if  KEYWORD_SET(wavefile) then print,'Wavelength solution file: ' + wavefile $
else begin
	print, 'No wavelength solution file found. Leaving ahk_profile.'
	return
endelse

if  KEYWORD_SET(slitfile) then print,'Slit edges file: ' + slitfile $
else begin
	print, 'No slit file found. Leaving ahk_profile.'
	return
endelse

;; Line wavelength we are interested in
linewave = 5006.7

;; Determine if science frame is on red or blue side of LRIS (will be used later to choose which line to search for in slit)
if STRMATCH(scifile, '*lblue*') then colour = 0 $
else colour = 1

;; Extract i-th slit from scifile using slit edges given by slitfile
sciarray = mrdfits(scifile)
wavearray = mrdfits(wavefile)
slitarray = mrdfits(slitfile)

;; Find how many slits are in science image
nslit = max(slitarray,/NaN)
print,'Number of slits detected: ' + STRTRIM(nslit, 2)

;; Create variable to hold sky location info of each slit
skylocstruct = {slitid: 0, index:Ptr_New()}
skylocarray = REPLICATE(skylocstruct, nslit)
skylocarray.index = PtrArr(nslit, /ALLOCATE_HEAP)



;; Start loop over slits to find profile in each
;FOR slitid=1,nslit DO BEGIN
FOR slitid=11,11 DO BEGIN
	splog, '**** SLIT ID: ', slitid, ' ****'

	thismask = (slitarray EQ slitid) ;; Array of same size as slitarray, contains 1's (when in slit) and 0's (anywhere else)

	;; Get mask parameters
	masksize = size(thismask)
	ncolmask = masksize(1)
	maskbeginindex = min(WHERE(thismask))
	maskendindex = max(WHERE(thismask))
	maskxbegin = maskbeginindex MOD ncolmask
	maskxend = maskendindex MOD ncolmask

	;; Trim pixels on either side (along slit direction, not dispersion direction)
	thismask[0:maskxbegin+5,*] = 0.
	thismask[maskxend-8:-5,*] = 0.
	

        thisslit = thismask*sciarray ;; Do not change sci image values in slit, make all others zero
	wave_thisslit = thismask*wavearray

	;; Find array indices which correspond to wavelengths near thisline.
	if colour EQ 0 then begin
		print, 'BLUE!!'
		thisline = ((wave_thisslit LT (linewave+2.)) and (wave_thisslit GT (linewave-2.))) ;; Array of 1's (when near emission line wavelength) and 0's (anywhere else)
		sci_profileregion = thisline*thisslit
	endif else begin
		print, 'RED!!'
		thisline = ((wave_thisslit LT (linewave+2.)) and (wave_thisslit GT (linewave-2.)))
		sci_profileregion = thisline*thisslit
	endelse

	;; Look +/- 25 pixels away from sci_profileregion (along wavelength direction) and build up a sky brightness profile.
	;; Interpolate across thisline, then subtract from science values in thisline (sci_profileregion).

	;; Start by getting thisline parameters
	linesize = size(thisline)
	ncolline = linesize(1)
	nrowline = linesize(2)
	linebeginindex = min(WHERE(thisline))
	lineendindex = max(WHERE(thisline))
	linexbegin = linebeginindex MOD ncolline
	linexend = lineendindex MOD ncolline

	;; Next find the sky flux and wavelengths on either side of the line profile region
	lineybegin = linebeginindex / ncolline
	lineyend = lineendindex / ncolline
	;print, 'liney begin and end', lineybegin, lineyend
	sky1 = thisslit[*,lineybegin+10:lineybegin+15]
	skywave1 = wavearray[*,lineybegin+25:lineybegin+30]
	;print, 'sky1', skywave1[2348:2355,*]
	sky2 = thisslit[*,lineyend-15:lineyend-10]
	skywave2 = wavearray[*,lineyend-30:lineyend-25]
	;print, 'sky2', skywave2[2348:2355,*]
	sky = [[sky2],[sky1]]
	skywave = [[skywave2],[skywave1]]
	;print, 'sky', sky[2348:2355,*]
	;print, 'size of sky', size(sky)

	;; Define array with wavelengths in sci_profileregion, then do the linear interpretation to find sky values at these wavelengths
	wave_profileregion = thisline * wave_thisslit
	sky_profileregion = MAKE_ARRAY(ncolline,nrowline)
	FOR i = 0,ncolline-1 DO sky_profileregion(i,*)=INTERPOL(sky(i,*),skywave(i,*),wave_profileregion(i,*), /NAN)

	;; Replace -NaN values output by INTERPOL with zeros.
	;; Either of the following two lines will accomplish this. Can choose best option once rest of code is finalized.
	sky_profileregion = sky_profileregion*thisline
	;FOREACH element, sky_profileregion DO sky_profileregion[WHERE(element NE element)] = 0

	;; Check output
	;print, size(sky_profileregion)
	;print, 'sky_profileregion', sky_profileregion[2348:2355,1756:1760]

	;; Subtract sky from science when sky GT 0 (avoids NaN's)
	;print, 'sci_profileregion before', sci_profileregion[2348:2355,1756:1760]
	goodindex = WHERE(sky_profileregion GT 0)
	sci_profileregion[goodindex] = sci_profileregion[goodindex] - sky_profileregion[goodindex]
	;print, 'sci_profileregion3 after', sci_profileregion[2348:2355,1759]	
	;plot, sci_profileregion[2340:2460,1759]

	;; Find coordinate in sci_profileregion (ie near thisline) which corresponds to the maximum pixel value after sky subtraction.
	;; Should find thisline, regardless of whether or not the wavelength calibration is slightly off or not.
	;; NOTE: Does not consider possibility of finding cosmic ray instead of thisline line at the moment.
	scipeak = max(sci_profileregion,peaklocation,/NaN)
	peakindex = ARRAY_INDICES(sci_profileregion, peaklocation)
	print, 'peak index at (x,y) = (', peakindex[0], ',', peakindex[1], ')'


	;; Extract profile across row containing peak value
	slitprofile = sci_profileregion[*,peakindex[1]]
	slitsize = size(slitprofile)
  	nsamp = slitsize(1)
	slitindex = FINDGEN(nsamp)
	slitindex = slitindex + 1

	;;Only want to fit to non-zero values in slitprofile
	slitprofile2 = slitprofile[WHERE(slitprofile)]
	slitindex2 = slitindex[WHERE(slitprofile)]

	;;Fit Gaussian using mpfitpeak and subtract. Repeat until only noise remains. Redo fit using all found mpfitpeak as starting point for mpfitexpr.
	params = [] ;; Declare empty array to put params found from mpfitpeak
	slitprofile_rem = slitprofile2 ;; Copy slitprofile2 into variable which will have each new fit subtracted from it
	bkgd = min(slitprofile[WHERE(slitprofile GT 0)])


	REPEAT BEGIN
		mpnterms=5
   		yfit = mpfitpeak(slitindex2, slitprofile_rem, fitresult, NTERMS=mpnterms, /NAN, /POSITIVE)
		params = [params,fitresult]
		;testplot1 = plot(slitindex2, slitprofile_rem)
		;testplot2 = plot(slitindex2, yfit, /OVERPLOT)
		slitprofile_rem=slitprofile_rem-yfit
		print, fitresult
	ENDREP UNTIL fitresult(0)*fitresult(2) LT 400 ;; or use: fitresult(0) LT 6*abs(fitresult(3))
	
	CASE 1 OF
		(N_ELEMENTS(params)/mpnterms EQ 1): BEGIN
			print, 'One local maximum found in sky subtracted line profile.'

			expr = 'P[0] + P[1]*X + GAUSS1(X, P[2:4])'
			start = [params(3),params(4), params(1), params(2) , params(0)*params(2)]
			result = MPFITEXPR(expr, slitindex2, slitprofile2, 0, start, /WEIGHTS, /QUIET)
			print, 'Gaussian result: ', result

			yfitfinal = result(0)+result(1)*slitindex2+gauss1(slitindex2, result(2:4))

			;; Plot results
			p = PLOT(slitindex2, slitprofile2, yrange=[0,max(slitprofile2)+10], Title='MPFITPEAK Gaussian Fit - Slit '+STRTRIM(slitid,2))
			pgauss = PLOT(slitindex2, yfitfinal, 'r', thick=5, /OVERPLOT)
		END
		(N_ELEMENTS(params)/mpnterms EQ 2): BEGIN
			print, 'Two local maxima found in sky subtracted line profile.'

			expr = 'P[0] + P[1]*X + GAUSS1(X, P[2:4]) + GAUSS1(X, P[5:7])' 
			start = [params(3),params(4), params(1), params(2) , params(0)*params(2), params(6), params(7) , params(5)*params(7)]
			result = MPFITEXPR(expr, slitindex2, slitprofile2, 0, start, /WEIGHTS, /QUIET)
			print, 'Gaussian result: ', result

			yfitfinal = result(0)+result(1)*slitindex2+gauss1(slitindex2, result(2:4)) + gauss1(slitindex2, result(5:7))
	
			p = PLOT(slitindex2, slitprofile2, yrange=[0,max(slitprofile2)+10], Title='Gaussian Fit - Slit '+STRTRIM(slitid,2))
			pgauss = PLOT(slitindex2, yfitfinal, 'r', thick=5, /OVERPLOT)

		END
		(N_ELEMENTS(params)/mpnterms EQ 3): BEGIN
			print, 'Three local maxima found in sky subtracted line profile.'

			expr = 'P[0] + P[1]*X + GAUSS1(X, P[2:4]) + GAUSS1(X, P[5:7]) + GAUSS1(X, P[8:10])'
			start = [params(3),params(4), params(1), params(2) , params(0)*params(2), params(6), params(7) , params(5)*params(7), $
				 params(11), params(12) , params(10)*params(12)]
			result = MPFITEXPR(expr, slitindex2, slitprofile2, 0, start, /WEIGHTS, /QUIET)
			print, 'Gaussian result: ', result

			yfitfinal = result(0)+result(1)*slitindex2+gauss1(slitindex2, result(2:4)) + gauss1(slitindex2, result(5:7)) $
				+ gauss1(slitindex2, result(8:10))
	
			p = PLOT(slitindex2, slitprofile2, yrange=[0,max(slitprofile2)+10], Title='Gaussian Fit - Slit '+STRTRIM(slitid,2))
			pgauss = PLOT(slitindex2, yfitfinal, 'r', thick=5, /OVERPLOT)

		END
		(N_ELEMENTS(params)/mpnterms EQ 4): BEGIN
			print, 'Four local maxima found in sky subtracted line profile.'

			expr = 'P[0] + P[1]*X + GAUSS1(X, P[2:4]) + GAUSS1(X, P[5:7]) + GAUSS1(X, P[8:10]) + GAUSS1(X, P[11:13])'
			start = [params(3),params(4), params(1), params(2) , params(0)*params(2), params(6), params(7) , params(5)*params(7), $
				 params(11), params(12) , params(10)*params(12), params(16), params(17) , params(15)*params(17)]
			result = MPFITEXPR(expr, slitindex2, slitprofile2, 0, start, /WEIGHTS, /QUIET)
			print, 'Gaussian result: ', result

			yfitfinal = result(0)+result(1)*slitindex2+gauss1(slitindex2, result(2:4)) + gauss1(slitindex2, result(5:7)) $
				+ gauss1(slitindex2, result(8:10)) + gauss1(slitindex2, result(11:13))
	
			p = PLOT(slitindex2, slitprofile2, yrange=[0,max(slitprofile2)+10], Title='Gaussian Fit - Slit '+STRTRIM(slitid,2))
			pgauss = PLOT(slitindex2, yfitfinal, 'r', thick=5, /OVERPLOT)

		END
		(N_ELEMENTS(params)/mpnterms EQ 5): BEGIN
			print, 'Five local maxima found in sky subtracted line profile.'

			expr = 'P[0] + P[1]*X + GAUSS1(X, P[2:4]) + GAUSS1(X, P[5:7]) + GAUSS1(X, P[8:10]) + GAUSS1(X, P[11:13]) + GAUSS1(X, P[14:16])'
			start = [params(3),params(4), params(1), params(2) , params(0)*params(2), params(6), params(7) , params(5)*params(7), $
				 params(11), params(12) , params(10)*params(12), params(16), params(17) , params(15)*params(17), $
				params(21), params(22) , params(20)*params(22)]
			result = MPFITEXPR(expr, slitindex2, slitprofile2, 0, start, /WEIGHTS, /QUIET)
			print, 'Gaussian result: ', result

			yfitfinal = result(0)+result(1)*slitindex2+gauss1(slitindex2, result(2:4)) + gauss1(slitindex2, result(5:7)) $
				+ gauss1(slitindex2, result(8:10)) + gauss1(slitindex2, result(11:13)) + gauss1(slitindex2, result(14:16))
	
			p = PLOT(slitindex2, slitprofile2, yrange=[0,max(slitprofile2)+10], Title='Gaussian Fit - Slit '+STRTRIM(slitid,2))
			pgauss = PLOT(slitindex2, yfitfinal, 'r', thick=5, /OVERPLOT)

		END
		(N_ELEMENTS(params)/mpnterms EQ 6): BEGIN
			print, 'Six local maxima found in sky subtracted line profile.'

			expr = 'P[0] + P[1]*X + GAUSS1(X, P[2:4]) + GAUSS1(X, P[5:7]) + GAUSS1(X, P[8:10]) + GAUSS1(X, P[11:13]) + GAUSS1(X, P[14:16]) + GAUSS1(X, P[17:19])'
			start = [params(3),params(4), params(1), params(2) , params(0)*params(2), params(6), params(7) , params(5)*params(7), $
				 params(11), params(12) , params(10)*params(12), params(16), params(17) , params(15)*params(17), $
				params(21), params(22) , params(20)*params(22), params(26), params(27) , params(25)*params(27)]
			result = MPFITEXPR(expr, slitindex2, slitprofile2, 0, start, /WEIGHTS, /QUIET)
			print, 'Gaussian result: ', result

			yfitfinal = result(0)+result(1)*slitindex2+gauss1(slitindex2, result(2:4)) + gauss1(slitindex2, result(5:7)) $
				+ gauss1(slitindex2, result(8:10)) + gauss1(slitindex2, result(11:13)) + gauss1(slitindex2, result(14:16)) $
				+ gauss1(slitindex2, result(17:19))
	
			p = PLOT(slitindex2, slitprofile2, yrange=[0,max(slitprofile2)+10], Title='Gaussian Fit - Slit '+STRTRIM(slitid,2))
			pgauss = PLOT(slitindex2, yfitfinal, 'r', thick=5, /OVERPLOT)

		END
	ELSE: BEGIN
		print, 'Error identifying peaks.'
		result = 0.00
		yfitfinal = 0.00
		p = PLOT(slitindex2, slitprofile2, yrange=[0,max(slitprofile2)+10], Title='Gaussian Fit - Slit '+STRTRIM(slitid,2))
	END
	ENDCASE

	;; Append final mean, fwhm, area under curve (i.e. mpfitexpr results) to data file
	openw, 1, 'ahk_profile_test.dat', /append
	printf, 1, slitid, result
	close,1

	;; Extract sky location to pass back to long_reduce
	;; sky location defined as when fitted gaussian(s) less than some fraction of peak value. Change this?
	;skylocindex = WHERE(yfitfinal LE (1.05*result(0)))
	skylocarray[slitid-1].slitid = slitid
	*skylocarray[slitid-1].index = WHERE(yfitfinal LE (1.05*result(0))) ;;NB: this gives index starting at 0 for each new slit.
	print, *skylocarray[slitid-1].index


ENDFOR ;End for loop over all slits


end
