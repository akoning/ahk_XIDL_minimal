;+
; NAME:
;   ewr_skymask
;
; PURPOSE:
;   Find a place where we like to fit the sky.  
;
; CALLING SEQUENCE:
;   objstruct = ewr_skymask( image, tset_slits=, $
;    [fwhm=, nperslit=, peakthresh=, skymask=, objmask= ] )
;
; INPUTS:
;   image      - Image for finding the slits, which would typically
;                be a flat-field image, an arc image, or a sum of those
;   tset_slits - Trace sets with slit start/end positions
;
; OPTIONAL INPUTS:
;   fwhm       - FWHM for convolving flux along the slit before peak-finding;
;                default to 3.0 pix.
;                (Also, do not allow peaks closer to each other than FWHM pix.)
;   nperslit   - Do not find more than this many objects per slit; default to 10
;   peakthresh - Flux threshhold for finding objects; the flux must be
;                at least this fraction of the brightest object in each slit;
;                default to 0.
;   absthresh  - Absolute flux threshold for finding objects; the peakflux must 
;                be at least htis large; default to 0. If both peakthresh and 
;                absthresh are set, absthresh overrides peakthresh. 
;   SIG_THRESH - Sigma threshold for objects  [default=5.0]
;   OBJTHRESH  - threshold for object masking
; OUTPUTS:
;   objstruct  - Structure with object parameters
;
; OPTIONAL OUTPUTS:
;   skymask    - Image of the cleanest sky pixels, =1 for good sky pixels
;
; COMMENTS:
;   Each slit is smashed to a single vector in the spatial direction,
;   taking care to align each wavelength in the spatial direction first.
;   On this smashed flux vector, we then search for significant peaks.
;   We then compute the expected position of each object at each wavelength,
;   by assuming that each object always falls in the same position within
;   the slit.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   extract_asymbox2()
;   find_npeaks()
;   gaussian()
;   long_slits2mask()
;   long_slits2x()
;   splog
;   traceset2xy
;
; INTERNAL SUPPORT ROUTINES:
;   long_obj_create()
;
; REVISION HISTORY:
;   10-Mar-2005  Written by D. Schlegel, LBL
;   30-Jun-2014  Hacked over to a skymasker by Erik Rosolowsky
;-
;------------------------------------------------------------------------------

function ewr_skymask, image, tset_slits=tset_slits $
                      , fwhm = fwhm1, nperslit = nperslit1 $
                      , ncoeff = ncoeff1 $
                      , peakthresh = peakthresh1, ABSTHRESH = ABSTHRESH1 $
                      , skymask = skymask_new, objmask = objmask $
                      , OBJTHRESH = OBJTHRESH, CRUDE = CRUDE $
                      , sky_fwhm = sky_fwhm, sigma_sky = sigma_sky $
                      , SILENT = SILENT, POS_SET = POS_SET $
                      , PEAK_SMTH=PEAK_SMTH, SIG_THRESH=sig_thresh1 $
                      , HAND_X = HAND_X $
                      , HAND_MINX = HAND_MINX1, HAND_MAXX = HAND_MAXX1 $
                      , HAND_Y = HAND_Y $
                      , HAND_FWHM = HAND_FWHM1, HAND_SUB = HAND_SUB1 $
                      , INVVAR = invvar $
                      , VERBOSE = VERBOSE, STDTRACE = STDTRACE, ISLIT = ISLIT, $
                      WAVEIMG = WAVEIMG, WAVEMASK = WAVEMASK, $
                      skywavemask = skywavemask, nudgelam = nudgelam, donudge = donudge

  ;;----------
                                ; Test inputs

  if (keyword_set(image) EQ 0) then $
     message, 'Must specify IMAGE'

  if NOT keyword_set(invvar) then  invvar = 1.0/(abs(image) + 100.)
  
                                ;----------
                                ; Set defaults

  dims = size(image, /dimens)
  nx = dims[0]
  if nx LE 2 then begin
     splog, 'Not enough columns in IMAGE', nx
     return, 0
  endif
  ;; Default values for TRACE_CRUDE
  if (NOT keyword_set(ksize)) then ksize = 5
  if (ksize LT 1 OR ksize GE nx/2-1) then $
     message, 'Invalid kernel size KSIZE'
  if (NOT keyword_set(nave)) then nave = 3
  if (NOT keyword_set(maxshifte)) then maxshifte = 0.1
  if (NOT keyword_set(maxshift0)) then maxshift0 = 1.0
  IF NOT KEYWORD_SET(PEAK_SMTH) THEN PEAK_SMTH=5.0
  IF NOT KEYWORD_SET(CRUDE_TOL) THEN CRUDE_TOL = 30.0D
  IF NOT KEYWORD_SET(NCOEFF1) THEN NCOEFF = 5 $
  ELSE NCOEFF = NCOEFF1
  func_crude = tset_slits[0].func
  dim_crude = size(tset_slits.COEFF, /dim)
  ncoeff_crude = dim_crude[0]

  ny = dims[1]
  ymid = ny / 2
  if (keyword_set(fwhm1)) then fwhm = fwhm1 else fwhm = 3.0 ; in pixels
;   if (NOT keyword_set(fwhm)) then fwhm = 3.0 ; in pixels
  if (keyword_set(peakthresh1)) then peakthresh = peakthresh1 $
  else peakthresh = 0.
  if (keyword_set(SIG_THRESH1)) then SIG_THRESH = SIG_THRESH1 $
  else SIG_THRESH=5.0
  if (keyword_set(absthresh1)) then absthresh = absthresh1 $
  else absthresh = 0.
  if (keyword_set(nperslit1)) then nperslit = nperslit1 $
  else nperslit = 10
  IF NOT KEYWORD_SET(OBJTHRESH) THEN OBJTHRESH = 0.5D

  objstruct = 0
  traceset2xy, tset_slits[0], yy1, xx1
  traceset2xy, tset_slits[1], yy2, xx2

  slitmask = long_slits2mask(tset_slits, nslit = nslit)
  skymask = (slitmask GT 0)
  objmask = slitmask*0

; Make gaussian kernels
  ngpix = (2 * ceil(1.5*fwhm) + 1) < (nx-1)
  width = 3
  gkern = gaussian(findgen(ngpix) - (ngpix-1)/2, [1.,0.,fwhm/2.305])
  gkern = gkern / total(gkern)

;----------
; Loop over each slit
  
  IF KEYWORD_SET(ISLIT) THEN BEGIN
     IF islit GT nslit THEN message $
        , 'ERROR: islit not found. islit cannot be larger than nslit'
     nreduce = 1
     slit_vec = [islit]
  ENDIF ELSE BEGIN
     nreduce = nslit
     slit_vec = lindgen(nslit) + 1L
  ENDELSE

; Calculate the wavelength offset for this mask w.r.t. expected
; values.

  minlam = min(waveimg[where(waveimg gt 0)])
  maxlam = max(waveimg[where(waveimg gt 0)])

  if keyword_set(donudge) then begin 
     for ii = 0,n_elements(skywavemask)-1 do begin 
        if (skywavemask[ii] lt minlam) or (skywavemask[ii] gt maxlam) then continue
        dlam = (waveimg-skywavemask[ii])
        idx = where(abs(dlam) lt 7,ct)
        if ct gt 0 then begin
           offset = (n_elements(offset) eq 0) ? [dlam[idx]] : [offset,dlam[idx]]
           vals = (n_elements(vals) eq 0) ? [image[idx]] : [vals,image[idx]] 
        endif
     endfor
     idx = where(vals gt 10*mad(vals)) ; Find the bright wing of the distribution
     nudgelam = total(offset[idx]*vals[idx])/total(vals[idx])
  endif else nudgelam = 0.0

  mask = (invvar GT 0.0)
  image_size=size(invvar)
  linemask = bytarr(image_size[1],image_size[2])
  offmask = bytarr(image_size[1],image_size[2])
  wskymask = bytarr(image_size[1],image_size[2])
  win = 3.5                     ; Angstroms
  for ii = 0,n_elements(wavemask)-1 do linemask = linemask or abs(waveimg-wavemask[ii]-nudgelam) lt win 
  for ii = 0,n_elements(wavemask)-1 do offmask = offmask or (abs(waveimg-wavemask[ii]-nudgelam) gt 1.5*win $
                                                             and abs(waveimg-wavemask[ii]-nudgelam) lt 2.5*win)
  for ii = 0,n_elements(skywavemask)-1 do wskymask = wskymask or (abs(waveimg-skywavemask[ii]-nudgelam)) lt win*2
  contmask = 1b-(linemask or wskymask)
  offmask = offmask*contmask
  skymask_out = bytarr(image_size[1],image_size[2])
  
  for jj = 0L, nreduce-1L DO BEGIN
     slitid = slit_vec[jj]
     thisimg = image * (slitmask EQ slitid)*mask
     thisline = linemask * (slitmask eq slitid)*mask
     thisoff =  offmask * (slitmask eq slitid)*mask
     thiscont =  contmask * (slitmask eq slitid)*mask
     if n_elements(waveimg) eq n_elements(image) then $  
        thiswave = waveimg * (slitmask EQ slitid)*mask

     ximg = long_slits2x(tset_slits, slitid = slitid)
     ;; Smash the image (for this slit) into a single flux vector.
     ;; How many pixels wide is the slit at each Y?
     xsize = xx2[*,slitid-1] - xx1[*,slitid-1]
     ;; How many pixels to sample this at?
     nsamp = ceil(djs_median(xsize))
;      nsamp = ceil(max(xsize)) ; changed by JFH aug 2007. Not sure if will
;;     cause problems but using the maximum nsamp doesn't seem to make sense
     
;     Mask skypixels with 2 FWHM of edge
;
     IF NOT KEYWORD_SET(SKY_FWHM) THEN SKY_FWHM = FWHM
     border = (2*sky_fwhm/nsamp) < 0.03
     skymask = skymask * (1B-((slitmask EQ slitid) AND $
                              ((ximg LT border) OR (ximg GT 1 - border))))
     
     left_asym = rebin(xx1[*,slitid-1], ny, nsamp) + $
                 (xsize/nsamp) # (findgen(nsamp))
     right_asym = left_asym + (xsize/nsamp) # replicate(1,nsamp) 

                                ; Generate rectified spectra and equivalent masks for this slit
     flux_spec = extract_asymbox2(thisimg, left_asym, right_asym)
     line_mask_spec = extract_asymbox2(thisline, left_asym, right_asym)
     off_mask_spec = extract_asymbox2(thisoff, left_asym, right_asym)
     cont_mask_spec = extract_asymbox2(thiscont, left_asym, right_asym)

     if total(linemask) eq 0 then begin
        fluxvec = djs_avsigclip(flux_spec, 1, sigrej = 25) 
        fluxsub = fluxvec - median(fluxvec)
     endif else begin

        onvec = total(flux_spec*line_mask_spec,1)/total(line_mask_spec,1)
        offvec = total(flux_spec*off_mask_spec,1)/total(off_mask_spec,1)
        contvec = total(flux_spec*cont_mask_spec,1)/total(cont_mask_spec,1)

                                ; This is the findback algorithm from JAC
                                ; Apply first the continuum region
        window = 5  ; This is the size of find-back algorithm feature. 
        back = contvec
        for j = -window,window do back = back < shift(back,j)
        for j = -window,window do back = back > shift(back,j)
        contback = smooth(back,window,/edge_trun,/nan)

        window = 5
        fluxsub = (onvec-offvec)>0
        back = fluxsub
        for j = -window,window do back = back < shift(back,j)
        for j = -window,window do back = back > shift(back,j)
        onoffback = smooth(back,window,/edge_trun,/nan)

        if total(fluxsub gt 0) lt 10 then begin
           fluxvec = djs_avsigclip(flux_spec, 1, sigrej = 25) 
           fluxsub = fluxvec - median(fluxvec)
           back = fluxsub
           for j = -window,window do back = back < shift(back,j)
           for j = -window,window do back = back > shift(back,j)
           onoffback = smooth(back,window,/edge_trun,/nan)
        endif
     endelse

                                ; Find regions that are large w.r.t. their local maximum
     skymaskvec = contvec/contback gt 1.25
     skymaskvec = morph_close(skymaskvec,fltarr(window)+1)
     l = label_region(skymaskvec)
     nelts = n_elements(skymaskvec)
     for kk = 1,max(l) do begin
        ind = where(l eq kk)
        minx = float(min(ind))/nElts
        maxx = float(max(ind))/nElts
        skymask_out = skymask_out + (ximg gt minx)*(ximg le maxx)
     endfor


;     linemaskvec = (fluxsub/onoffback) gt 1.1
;     linemaskvec = (onvec-offvec)/offvec gt 0.25
;     sind = sort(fluxsub)
;     cutoff = (fluxsub)[sind[0.15*n_elements(sind)]] 
;     linemaskvec2 = fluxsub gt cutoff
;     linemaskvec = linemaskvec < linemaskvec2

     pcts = cgPercentiles(onvec-offvec,percentiles=[0.02275,0.158])
     linemaskvec = (onvec-offvec) gt (2*pcts[1]-1*pcts[0])
     linemaskvec = morph_open(linemaskvec,fltarr(3)+1)

;     pcts = cgPercentiles(onvec-offvec,percentiles=[0.2])
;     linemaskvec = (onvec-offvec gt pcts[0])
     l = label_region(linemaskvec)
     nelts = n_elements(linemaskvec)
     for kk = 1,max(l) do begin
        ind = where(l eq kk)
        minx = float(min(ind))/nElts
        maxx = float(max(ind))/nElts
        skymask_out = skymask_out + (ximg gt minx)*(ximg le maxx);*linemask
     endfor
  endfor
  skymask_out = skymask_out < 1
  skymask_out = (1b-skymask_out)*skymask
  return, skymask_out
end

;------------------------------------------------------------------------------


;     fluxconv = convol(fluxsub, gkern, /center, /edge_truncate)   
;; ;     null = ahk_profile(fluxconv,ymodel)
;; ; Set sky as below 35%ile of the brightness.

;;      sind = sort(fluxconv)
;;      firstcut = fluxconv[sind[0.35*n_elements(fluxconv)]]

;;      skymaskvec = fluxconv le firstcut

     ;; djs_iterstat, fluxconv[where(skymaskvec)], sigma = skythresh, sigrej = 1.5, median = skymedian

     ;; finalcut = (firstcut > (2*skythresh+skymedian)) > 0.05*max(fluxconv)

     ;; if max(fluxconv) lt 5*finalcut then begin
     ;;       fluxvec = djs_avsigclip(flux_spec, 1, sigrej = 25) 
     ;;       djs_iterstat, fluxvec, sigma = skythresh, sigrej = 1.5, median = skymedian
     ;;       fluxsub = fluxvec - skymedian
     ;;       fluxconv = convol(fluxsub, gkern, /center, /edge_truncate)
     ;; endif

;     skymaskvec = fluxconv lt finalcut
