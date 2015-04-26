;+
; NAME:
;   long_proc
;
; PURPOSE:
;   Overscan subtract, bias subtract, and flat field a CCD image. The routine 
;   calls long_oscan and flat fields.
;
; CALLING SEQUENCE:
;   long_proc, filename, [ flux, invvar, hdr=, $
;    biasfile=, pixflatfile=, adderr=, /verbose ]
;
; INPUTS:
;   filename   -  name of a image file
;
; OPTIONAL INPUTS:
;   biasfile   - File with average bias
;   pixflatfile- File with pixel flat
;   adderr     - Additional error to add to the formal errors, as a fraction
;                of the flux; default to 0.01 (1 per cent).
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   flux       - Image 
;   invvar     - Inverse variance image
;   hdr        - FITS header
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   Need to add cosmic ray zapping.
;
; PROCEDURES CALLED:
;   divideflat
;   headfits()
;   mrdfits()
;   long_oscan
;
; REVISION HISTORY:
;   10-Mar-2005  Written by J. Hennawi (UCB), D. Schlegel (LBL)
;-
     ;; The following values are good for Jun2012, from
     ;; mods1r20120623.0024  -- JXP 17 Aug 2012
;------------------------------------------------------------------------------
pro long_proc, filename, flux, invvar, hdr = hdr $
               , biasfile = biasfile, pixflatfile = pixflatfile $
               , adderr = adderr, verbose = verbose, bin= bin $
               , rnoise = rnoise, gain = gain, CCDONLY = CCDONLY $
               , TRANSFORM = TRANSFORM

  if  N_params() LT 2  then begin 
     print,'Syntax - ' + $
           'long_proc, fil, flux [invvar], HDR=, BIASFILE=, GAIN= [v1.0]' 
     return
  endif 

  minval = 0.5
  if (n_elements(adderr) EQ 0) then adderr = 0.01
  hdr = xheadfits(filename)

  if (size(hdr, /tname) NE 'STRING') then begin
     splog, 'Invalid FITS header for file ', filename
     flux = 0
     invvar = 0
     return
  endif

  IF sxpar(hdr, 'COADD_2D') THEN BEGIN
     flux = xmrdfits(filename, 0, hdr, /silen)
     invvar = xmrdfits(filename, 1, /silen)
     long_oscan, filename, img_dum, ivar_dum, bin = bin
     RETURN
  ENDIF

;; XGMOS
  if sxpar(hdr,'XGMOS') then begin
     flux = xmrdfits(filename,0,hdr,/silen)
     invvar = xmrdfits(filename,1,/silen)
     cbin = strtrim(strsplit(sxpar(hdr,'CCDSUM'),/extract),2)
     bin = reverse(long(cbin))
     return
  endif

;  If FLUX and INVVAR are not requested, then return with just the header
  if (arg_present(flux) EQ 0 AND arg_present(invvar) EQ 0) then return

  if (arg_present(invvar)) then $
     long_oscan, filename, flux, invvar, hdr = hdr, bin=bin, verbose=verbose $
                 , rnoise = rnoise, gain = gain, CCDONLY = CCDONLY, TRANSFORM = TRANSFORM $
  else $
     long_oscan, filename, flux, hdr = hdr, verbose = verbose, bin = bin $
                 , rnoise = rnoise, gain = gain, CCDONLY = CCDONLY, TRANSFORM = TRANSFORM

  if (keyword_set(invvar) AND adderr GT 0) then begin
     gmask = invvar GT 0        ; =1 for good points
     invvar = gmask / (1.0 / (invvar + (1-gmask)) + $
                       adderr^2 * (abs(flux))^2)
  endif
  
  ndim = size(flux, /n_dimen)
  dims = size(flux, /dimens)

  if (keyword_set(biasfile)) then begin
     biasimg = xmrdfits(biasfile, 0, biashdr, $
                        silent = (keyword_set(verbose) EQ 0))
     if (size(biasimg, /n_dimen) NE ndim $
         OR total(size(biasimg, /dimens) NE dims) NE 0) then $
            message, 'Dimensions of image and bias image do not agree'
     
     flux = flux - biasimg
     ;; Redo INVVAR for Detector without overscan
     if strtrim(sxpar(hdr[*, 0], 'TELESCOP'), 2) EQ 'CA-2.2' then begin
        ;; This assumes that windowing was done.  Might not be true!
        invvar = 1.0/(abs(flux - sqrt(2.0)*rnoise) +rnoise^2) 
     endif
  endif

  if (keyword_set(pixflatfile)) then begin
     flatimg = xmrdfits(pixflatfile, 0, biashdr, $
                        silent = (keyword_set(verbose) EQ 0))
     if (size(flatimg, /n_dimen) NE ndim $
         OR total(size(flatimg, /dimens) NE dims) NE 0) then $
            message, 'Dimensions of image and flat image do not agree'
     
     divideflat, flux, flatimg, invvar = invvar, minval = minval, $
                 quiet = (keyword_set(verbose) EQ 0)
  endif

   if stregex(sxpar(reform(hdr,n_elements(hdr)),'INSTRUME'),'LRISBLUE',/bool)*$
      (1b-stregex(sxpar(reform(hdr,n_elements(hdr)),'SLITNAME'),'long',/bool))*$
      ((1b-sxpar(reform(hdr,n_elements(hdr)),'LAMPS') EQ '0,0,0,0,0,0') OR (1b-stregex(sxpar(reform(hdr,n_elements(hdr)),'TRAPDOOR'),'closed',/bool))) then begin
                                ; Blanks the line on the chip gap.
                                ; Except for longslit and biases
      flux[2047,*]=0.0
   endif

;; ; the LRIS bar occupies 7.0 mm of space at the position of the slit
;; ; mask.  Going into the AUTOSLIT software suggests that the
;; ; astrometric offset between these should be 7mm.  
;; ; 0.724262821 mm / arcsec
;; ; 0.135" / pixel 
;; ; So a 7 mm gap goes to 9.665" goes to 71.56925 pixels

;;      sz = size(flux)
;;      shiftpix = round(71.56925)
;;      flux_new = fltarr(sz[1]+shiftpix,sz[2])
;;      flux_new[0:2047,*] = flux[0:2047,*]
;;      flux_new[(2047+shiftpix+1):*,*] = flux[2048:*,*]
;;      if n_elements(invvar) eq n_elements(flux) then begin  
;;         invvar_new  = fltarr(sz[1]+shiftpix,sz[2])
;;         invvar_new[0:2047,*] = invvar[0:2047,*]
;;         invvar_new[(2047+shiftpix+1):*,*] = invvar[2048:*,*]
;;         invvar = invvar_new
;;      endif
;;      flux = flux_new
;;      sxaddpar,hdr,'NAXIS1',sz[1]+shiftpix

;;   endif
   bad_variances = where(invvar gt 1.0,bad_count)
   if bad_count gt 0 then begin
      invvar[bad_variances]=0.0
      flux[bad_variances]=0.0
   endif
   
  return
end

;------------------------------------------------------------------------------
