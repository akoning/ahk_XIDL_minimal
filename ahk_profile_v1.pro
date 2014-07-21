pro ahk_profile_v1,scifile,wavefile=wavefile,slitfile=slitfile

;+
; NAME:
;   ahk_profile_v1
;
; PURPOSE:
;   Loop over slits and determine spatial profile of each based on a single line (Halpha 6562.82 in red, Hbeta 4861.33 or OIII 5006.84 in blue).
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

  ;; Start loop over slits to find profile in each
  ;FOR slitid=1,nslit DO BEGIN
  FOR slitid=9,9 DO BEGIN
	slit_location = WHERE(slitarray eq i,count)

	;; Extract portion of science frame and wavelength solution which correspond to slit location.
	wave_slitlocation = wavearray[slit_location]
	sci_slitlocation = sciarray[slit_location]
	print, size(sci_slitlocation)

	;; Find array indices which correspond to wavelengths near Hbeta (blue) or Halpha (red).
	if colour EQ 0 then begin
		print, 'BLUE!!'
		Hbeta = WHERE((wave_slitlocation LT 4866) and (wave_slitlocation GT 4856), count)
		wave_profileregion = wave_slitlocation[Hbeta]
		print, 'wavelengths'
		print, wave_profileregion[1:5], wave_profileregion[500:505]
		sci_profileregion = sci_slitlocation[Hbeta]
	endif else begin
		print, 'RED!!'
		Halpha = WHERE((wave_slitlocation LT 6568) and (wave_slitlocation GT 6558), count)
		wave_profileregion = wave_slitlocation[Halpha]
		sci_profileregion = sci_slitlocation[Halpha]
	endelse

	;; Find coordinate in sci_profileregion which corresponds to the maximum pixel value (should tell us exactly where Hbeta or Halpha is, regardless of whether or not the wavelength calibration is slightly off or not) **Does not consider possibility of finding cosmic ray instead of Hbeta or Halpha line at the moment.
	scimax = max(sci_profileregion,maxindex,/NaN)
	print, size(sci_profileregion)
	print, i, scimax, maxindex
	print, size(sci_slitlocation)






;;;;!!!!!!!!!!
;; Want to find slitprofile = all values in same row as max, x centroid = sum x(n)val(n)/sum val(n), FWHM = literally take full width at half max (using linear interpolation?)

;; First, extract profile across row containing maximum pixel
  s = SIZE(sci_profileregion)
  ncol = s(1)
  row = maxindex / ncol
  print, 'ROW: ', row
  slitprofile = sci_slitlocation[*,row] ;This isnt doing what I want yet

;; Second, find the centroid
  

;; Third, calculate full-wdith at half maximum.Modelled after section in jds_objextract
  yhalf = 0.5 * scimax
  x0 = ROUND(xcen)
  nsamp = CEIL(max(xsize))
  if (x0 LT nsamp-1) then begin
  	i2 = (WHERE(profile[x0:nsamp-1] LT yhalf))[0]
  	xright = INTERPOL(x0+[i2-1,i2], profile[x0+i2-1:x0+i2], yhalf)
  endif else xright = 0
  if (x0 GT 0) then begin
  	i1 = (reverse(where(profile[0:x0] LT yhalf)))[0]
  	xleft = interpol([i1,i1+1], profile[i1:i1+1], yhalf)
  endif else xleft = 0
  fwhmmeas = 0.
  if (xleft NE 0 AND xright NE 0) then fwhmmeas = (xright - xleft) $
  else if (xleft NE 0) then fwhmmeas = 2*(xcen[ipeak] - xleft) $
  else if (xright NE 0) then fwhmmeas = 2*(xright - xcen[ipeak]

;; Put slit profile, FWHM, and centroid in output structure
  slitstructure = { profile: slitprofile, fwhm: fwhmmeas, centroid: xcen}
  plot, slitstructure.profile
  print, slitstructure.fwhm
  print, slitstructure.centroid

  ENDFOR ;End for loop over all slits

end
