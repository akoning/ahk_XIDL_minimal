;+
; NAME:
;   long_reduce_work
;
; PURPOSE:
;
;   Main program for the Low-redux pipeline.  This set of algorithms
;   runs mainly as a black box.
;
; CALLING SEQUENCE:
;  long_reduce, planfile, /clobber, /NOZAP, /NOFLEX, /NOHELIO 
;
; INPUTS:
;  planfile  -- File created by long_plan which guides the reduction
;               process
;
; OPTIONAL INPUTS:
; /NOFLEX  -- Do not apply flexure correction [necessary if your setup
;             has not been calibrated.  Contact JH or JXP for help if
;             this is the case.]
;  HAND_FWHM -- Set the FWHM of the object profile to this value (in
;  PROF_NSIGMA= -- Extend the region to fit a profile by hand 
; /NOHELIO -- Do not correct to heliocentric velocities
; /NOZAP   -- Do not flag CRs
; wavemaskType -- Added by AHK. Specifies which wavelengths to use
;                 when making FLUXMODEL. Can be 'Full', 'Balmer',
;                 'Metal', 'HiIon', 'LoIon'
;
; OUTPUTS:
;  (1) Various calibration files
;  (2) One multi-extension FITS file in Science per exposure containing
;  the extracted data and processed images
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
;   11-Mar-2005  Written by JH + SB
;-  
;-----------------------------------------------------------------------------
; BUGS:
;   CR zapping should be improved???
;
;-----------------------------------------------------------------------------
; The code that follows is the science frame reduction code

PRO long_reduce_work, filename, scifile, ISLIT = ISLIT $
                      , slitfile = slitfile, wavefile = wavefile $
                      , biasfile = biasfile, pixflatfile = pixflatfile $
                      , illumflatfile = illumflatfile $
                      , verbose = verbose, box_rad = box_rad1 $
                      , maxobj = maxobj1 $
                      , reduxthresh = reduxthresh, SIG_THRESH=sig_thresh $
                      , NOFLEX = noflex1, NOSHIFT = NOSHIFT1 $
                      , NOHELIO = NOHELIO , MAXFLEX=maxflex1  $
                      , MAXGOOD = maxgood $
                      , NOZAP = NOZAP, HAND_X = HAND_X, HAND_Y = HAND_Y $
                      , HAND_FWHM = HAND_FWHM, STD = STD $
                      , FILESTD = FILESTD1 $
                      , CHK = CHK, TRCCHK = TRCCHK, SKYTRACE = SKYTRACE1 $
                      , SKYSAMPLE = SKYSAMPLE1 $
                      , PROF_NSIGMA=prof_nsigma1 $
                      , ALLPROF_NSIGMA = ALLPROF_NSIGMA $
                      , NOLOCAL = nolocal $
                      , NSIGMA_CR = NSIGMA_CR1 $
                      , PSF_CR_FUDGE = PSF_CR_FUDGE1 $
                      , NITER=NITER, NOVAC=NOVAC, SN_GAUSS = SN_GAUSS1 $
                      , wavemaskType=wavemaskType
t0 = systime(1)

;if (NOT keyword_set(profile_filename)) then profile_filename = 0
if (NOT keyword_set(maxobj1)) then maxobj = 0 $ 
ELSE maxobj = maxobj1

;;----------
;; Read the raw science image
long_proc, filename, sciimg, sciivar, hdr = scihdr $
           , biasfile = biasfile, pixflatfile = pixflatfile, bin = bin
;;  What telescope are we using?
telescope = strcompress(sxpar(scihdr[*, 0], 'TELESCOP'), /rem)
;; Is this a 2-d coadd being re-fed through the pipeline?
COADD_2D = sxpar(scihdr[*, 0], 'COADD_2D')
;;  Determine detector specific parameters for reduction
long_reduce_params, scihdr, bin, skyfile = skyfile, anamorph = anamorph $
                    , bsp = bsp, SN_GAUSS = SN_GAUSS, SKYTRACE = SKYTRACE $
                    , SKYSAMPLE = SKYSAMPLE $
                    , NOSHIFT = NOSHIFT, NCCD = NCCD, PEAK_SMTH = PEAK_SMTH $
                    , FILESTD = FILESTD, FWHM = FWHM, BOX_RAD = BOX_RAD $
                    , flg_skyfile=flg_skyfile
;; override param default if standard is set
IF KEYWORD_SET(FILESTD1) THEN FILESTD = FILESTD1
IF KEYWORD_SET(FILESTD) THEN BEGIN
    splog, 'Using standard star trace as crutch from ' + filestd
    print, 'Using standard star trace as crutch from ' + filestd
    stdstruct = xmrdfits(filestd, 5, /silent)
    stdmax = max(stdstruct.PEAKFLUX, stdind)
    stdtrace = stdstruct[stdind].XPOS
ENDIF ELSE stdtrace = 0

IF KEYWORD_SET(STD) THEN BEGIN
   NOFLEX = 1
   SKYTRACE = 0
;   box_rad = 20  -- Best for slitless
;   NOLOCAL= 0
   box_rad = 50
   NOLOCAL= 1

   NOSHIFT = 1
   NITER = 1
;   PROF_NSIGMA = 20.0
;;   maxobj = 1
   sigrej = 25.0
   NOZAP = 1 ;; Do not zap cosmics since they are not masked in boxcar
;;   NSIGMA_CR = 10.0d
;;   PSF_CR_FUDGE = 0.5d ;; make the PSF sharper to avoid rejection for STDs
;;   PROF_NSIGMA = 20.0
;;  Over-ridden by JXP
;   reduxthresh=0.001 ;; makes sure we reduce other objects on the slit. 
     PEAK_SMTH = 20  ; EWR Added to deal with bright standards.
ENDIF ELSE BEGIN
   ;; If any of these were passed in, overwrite long_reduce_params values
   IF KEYWORD_SET(PSF_CR_FUDGE1) GT 0 THEN PSF_CR_FUDGE = PSF_CR_FUDGE1 $
   ELSE PSF_CR_FUDGE = 1.0d
   IF KEYWORD_SET(NSIGMA_CR1) THEN NSIGMA_CR = NSIGMA_CR1
   IF KEYWORD_SET(skytrace1) THEN SKYTRACE = SKYTRACE1
   IF n_elements(skysample1) THEN SKYSAMPLE = SKYSAMPLE1
   IF n_elements(SN_GAUSS1) THEN SN_GAUSS = SN_GAUSS1
   IF keyword_set(noshift1) THEN NOSHIFT = NOSHIFT1
   IF KEYWORD_SET(noflex1) THEN NOFLEX = NOFLEX1
   IF KEYWORD_SET(maxflex1) THEN MAXFLEX = MAXFLEX1
   IF KEYWORD_SET(box_rad1)  THEN BOX_RAD = BOX_RAD1
   IF KEYWORD_SET(PROF_NSIGMA1) THEN PROF_NSIGMA=PROF_NSIGMA1
   IF KEYWORD_SET(ALLPROF_NSIGMA1) THEN ALLPROF_NSIGMA = ALLPROF_NSIGMA1
   NOZAP = 0
ENDELSE

dims = size(sciimg, /dimens)
nx = dims[0]
ny = dims[1]

;; Read in slitmask structure
tset_slits = xmrdfits(slitfile, 1, silent = (keyword_set(verbose) EQ 0))
;; FLEXURE SHIFT SLITS FOR LRIS
IF NOT KEYWORD_SET(NOSHIFT) THEN $
   xshift = long_xcorr_slits(sciimg, tset_slits, /shift) 
;;   Generate slit position, wavelengths, and slit illumination 
ximg = long_slits2x(tset_slits, edgmask = edgmask, nslit = nslit)
slitmask = long_slits2mask(tset_slits) * (ximg GT 0. AND ximg LT 1.)
IF strcmp(telescope, 'Gemini-North') THEN BEGIN
   IF KEYWORD_SET(illumflatfile) THEN $
      slit_illum = xmrdfits(illumflatfile, silent = (keyword_set(verbose) EQ 0))
   waveimg = xmrdfits(wavefile, silent = (keyword_set(verbose) EQ 0), 0)
   piximg = xmrdfits(wavefile, silent = (keyword_set(verbose) EQ 0), 0)
ENDIF ELSE BEGIN
   IF KEYWORD_SET(illumflatfile) THEN $
      slit_illum = long_slitillum(illumflatfile, slitmask, ximg, edgmask)
   ;;   Reconstruct PIXIMG WAVEIMG and using coefficient sets
   pixset  = xmrdfits(wavefile, silent = (keyword_set(verbose) EQ 0), 1)
   wavesvfile =  repstr(wavefile, '.fits', '.sav')
   restore, wavesvfile
   piximg = long_wpix2image(pixset, tset_slits, XFIT = xfit $
                            , waveimg = waveimg, ISLIT = ISLIT)
ENDELSE
;   Read in wavelength solution structures
;fwhmset = xmrdfits(wavefile, silent = (keyword_set(verbose) EQ 0), 2)

;; If ISLIT and the OUTPUT image exists, then read in stuff to save
if keyword_set(ISLIT) then begin
    thisfile = findfile(scifile+'*', count = nct)
    if nct GT 0 then begin
        print, 'long_reduce_work:  Will save other slits as previously reduced.'
        modelivar = xmrdfits(scifile,1)
        skyimage = xmrdfits(scifile,2)
        objimage = xmrdfits(scifile,3)
        outmask = xmrdfits(scifile,4)
        final_struct = xmrdfits(scifile,5)
    endif
    islitmask = where(slitmask EQ ISLIT)
endif else nct = 0
        
;--------- 
;;  If there is no slit illumination function assume that it is unity
;;  If this  is a co-add, slit-illum already applied, and don't
;;  reapply. 
IF NOT KEYWORD_SET(slit_illum) OR COADD_2D THEN slit_illum = 0.0*sciimg +1.0 
; Allocate images which we will need later
objimaget = fltarr(nx, ny)
splog, 'Applying slit illumination'
;; Don't apply illumination corrections larger than 30%
gdpix = WHERE(slit_illum GT 0.6D AND slit_illum LT 1.4)
sciimg[gdpix]  = sciimg[gdpix]/slit_illum[gdpix]
sciivar[gdpix] = sciivar[gdpix]*slit_illum[gdpix]^2 
;; Trace sky lines for sky subtraction?
;;trcchk = 1
IF KEYWORD_SET(SKYTRACE) THEN BEGIN
    splog, 'Tracing sky lines'
    wavesvfile =  repstr(wavefile, '.fits', '.sav')
    restore, wavesvfile
    wstruct = long_wstruct(scihdr)
    pixset = long_wavepix(sciimg, tset_slits, fwhm = fwhmset.MEDIAN $
                          , box_radius = wstruct.radius $
                          , sig_thresh = wstruct.sig_wpix $
                          , pkwdth = wstruct.pkwdth $
                          , TOLER = wstruct.TOLER, CHK = TRCCHK $
                          , piximg_in = piximg, ISLIT = ISLIT)
    piximg = long_wpix2image(pixset, tset_slits, XFIT = xfit $
                             , waveimg = waveimg, ISLIT = ISLIT)
ENDIF

; EWR added this because the pixel fits seem to not be good enough.
;   waveimg = xmrdfits(wavefile, silent = (keyword_set(verbose) EQ 0), 0)
; Now removed because of better fits; the fits clearly track the slit
; edges better as would be expected.  occasionally fits do still go wonky

splog, 'Zapping cosmic rays using qzap.'
;----------
; AHK's magic qzap parameters
qzap,sciimg,sciimg_new,skyfiltsize=40,nsigma=11,fluxratio=0.15
sciimg = sciimg_new


splog, 'Finding objects on the slits: First pass'
;----------
;  Find objects on the slits

skywavemask = [  5224.137,   5238.747,   5256.083,   5577.338,   5889.95 ,$
                 5895.92 ,   5915.301,   5932.862,   5953.42 ,   5977.077,$
                 6003.911,   6300.304,   6363.78 ,   6533.044,   6553.617,$
                 6577.285,   6604.135,   6634.229,   6900.833,   6923.22 ,$
                 6949.044,   6978.413,   7011.41 ,   7316.282,   7340.885,$
                 7369.365,   7401.858,   7438.473,   7479.31 ,   7794.111,$
                 7821.503,   7853.36 ,   7889.851,   7931.124,   7964.65 ,$
                 7993.333,   8025.81 ,   8062.178,   8103.   ,   8399.17 ,$
                 8430.174,   8465.358,   8504.841,   8548.709,   8885.85 ,$
                 8919.635,   8958.084,   9001.346,   9049.539,   9102.77 ,$
                 9439.65 ,   9476.83 ,   9519.354,   9567.339,   9620.965,$
                 9872.137,   9914.673,   9961.932,  10013.986,  10082.468,$
                 10124.008,  10171.512,  10171.512,  10171.888,  10225.783]

IF (wavemaskType EQ 'Full') THEN BEGIN
   wavemask = [3711.97,3721.94,3726.03,3728.82,3734.37,3750.15,3770.63,$
               3797.90,3835.38,3868.75,3888.65,3970.07,4026.21,4068.60,$
               4076.35,4101.74,4340.47,4363.21,4471.50,4713.17,4861.33,$
               4921.93,4958.91,5006.84,5015.68,5047.74,5197.90,5270.40,$
               5517.71,5537.88,5754.64,5875.66,6312.10,6363.78,6548.10,$
               6562.77,6583.50,6678.16,6716.44,6730.82,7065.25,7135.80,$
               7235.00,7281.35,7319.45,7330.20,7751.43,8268.00,8545.38,$
               8598.39,8665.02,8727.12,8869.00,9014.91,9068.60]
ENDIF ELSE IF (wavemaskType EQ 'Balmer') THEN BEGIN
   wavemask = [3711.97,3734.37,3750.15,3770.63,3797.90,3835.38,3970.07,$
               4101.74,4340.47,4861.33,6562.77,8268.00,8545.38,8598.39,$
               8665.02,9014.91]
ENDIF ELSE IF (wavemaskType EQ 'Metal') THEN BEGIN
   wavemask = [3721.94,3726.03,3728.82,3868.75,3888.65,4026.21,4068.60,$
               4076.35,4363.21,4471.50,4713.17,4921.93,4958.91,5006.84,$
               5015.68,5047.74,5197.90,5270.40,5517.71,5537.88,5754.64,$
               5875.66,6312.10,6363.78,6548.10,6583.50,6678.16,6716.44,$
               6730.82,7065.25,7135.80,7235.00,7281.35,7319.45,7330.20,$
               7751.43,8727.12,8869.00,9068.60]
ENDIF ELSE IF (wavemaskType EQ 'HiIon') THEN BEGIN
                                ; Include high ionization lines OIII,
                                ; ArIII, NeIII, SIII, ClIII (only
                                ; bright lines, i.e. those also in
                                ; Full wavemask)
   wavemask = [3721.94,3868.75,4363.21,4958.91,5006.84,5517.71,5537.88,$
               6312.10,7135.80,7751.43,9068.60]
ENDIF ELSE IF (wavemaskType EQ 'LoIon') THEN BEGIN
                                ; Include low ionization lines Balmer,
                                ; OII, NII, SII (again, only those
                                ; also in Full wavemask)
   wavemask = [3711.97,3726.03,3728.82,3734.37,3750.15,3770.63,3797.90,$
               3835.38,3970.07,4068.60,4076.35,4101.74,4340.47,4861.33,$
               5754.64,6548.10,6562.77,6583.50,6716.44,6730.82,7319.45,$
               7330.20,8268.00,8545.38,8598.39,8665.02,9014.91]
ENDIF

;wavemask = [3727.00,4861.33,$
;            4958.91,5006.84,6363.78]

; Teh Dopplerz
vsys = -190 ; km/s for M33
deltalam = wavemask*vsys/3e5
wavemask = wavemask + deltalam


; Let's make this into a skymasking routine?

;; objstruct1 = ewr_objfind(sciimg, tset_slits = tset_slits, invvar = sciivar $
;;                           , skymask = skymask1, objmask = objmask1 $
;;                           , nperslit = maxobj, peakthresh = reduxthresh $
;;                           , fwhm = FWHM, PEAK_SMTH = PEAK_SMTH $
;;                           , SIG_THRESH = SIG_THRESH $
;;                           , HAND_X = HAND_X, HAND_Y = HAND_Y $
;;                           , HAND_FWHM = HAND_FWHM, STDTRACE = STDTRACE $
;;                           , ISLIT = ISLIT, WAVEIMG = waveimg, wavemask = wavemask)

skymask1 = ewr_skymask(sciimg, tset_slits = tset_slits, invvar = sciivar $
                       , skymask = skymask1, objmask = objmask1 $
                       , nperslit = maxobj, peakthresh = reduxthresh $
                       , fwhm = FWHM, PEAK_SMTH = PEAK_SMTH $
                       , SIG_THRESH = SIG_THRESH $
                       , HAND_X = HAND_X, HAND_Y = HAND_Y $
                       , HAND_FWHM = HAND_FWHM, STDTRACE = STDTRACE $
                       , ISLIT = ISLIT, WAVEIMG = waveimg, wavemask = wavemask, $
                       skywavemask = skywavemask,nudgelam=nudgelam)

splog, 'Aperture masked sky subtraction'
skyimaget = ewr_skysub(sciimg, sciivar, piximg, slitmask, skymask1 $
                        , edgmask, bsp = bsp, ISLIT = ISLIT, CHK = chk,$
                       waveimg=waveimg,wavemask=wavemask,nudgelam=nudgelam)
;stop
if keyword_set(ISLIT) and nct GT 0 then $
  skyimage[islitmask] = skyimaget[islitmask] $
  else skyimage=skyimaget
;stop
IF NOT KEYWORD_SET(NOZAP) THEN BEGIN
    splog, 'Comods1r.20120130.0041.fits.gzsmic ray rejection'
    IF KEYWORD_SET(FWHMSET) THEN sigma_psf = $
      djs_median(fwhmset.median)/2.35482D $
    ELSE sigma_psf = 3.0D/2.35482D
    ;; For standards be stricter about CR REJECTION
    sigma_psf = sigma_psf*PSF_CR_FUDGE 
    ;;   Take the PSF width to be that of the spectral direction.  
    ;;   This prevents the routine from rejecting sky lines
    ;;   This is a description of the 3x3 core of the 2D PSF for reject_cr.pro
    ;;    
    ;;    PSFVALS[1]  PSFVALS[0]   PSFVALS[1] 
    ;;    PSFVALS[0]          1.   PSFVALS[0]
    ;;    PSFVALS[1]  PSFVALS[0]   PSFVALS[1]
    psfvals = [exp(-1.0/(2*sigma_psf^2)), exp(-2.0/(2*sigma_psf^2))]
    crmask  = psf_reject_cr(sciimg-skyimage, sciivar, psfvals $
                            , satmask = (sciimg-skyimage) GT 8d4 $
                            , nsigma = nsigma_CR)
    sciivar = sciivar*(crmask EQ 0)
    ;; Do we need to sky-subtract before CR rejection???
 ENDIF

splog, 'Finding objects in sky-subtracted image: Second pass'
; Redo object finding on sky subtracted image
;IF KEYWORD_SET(OBJSTRUCT1) THEN FWHM = djs_median(objstruct1.FWHM)
objstruct = long_objfind(sciimg-skyimage, tset_slits = tset_slits $
                         , invvar = sciivar, skymask = skymask $
                         , objmask = objmask, nperslit = maxobj $
                         , peakthresh = reduxthresh $
                         , SIG_THRESH = SIG_THRESH $
                         , fwhm = FWHM, PEAK_SMTH = PEAK_SMTH $
                         , HAND_X = HAND_X, HAND_Y = HAND_Y $
                         , HAND_FWHM = HAND_FWHM, STDTRACE = STDTRACE $
                         , ISLIT = ISLIT, wavemask = wavemask,waveimg=waveimg)
; EWR -- I like the first sky subtraction just fine!
;; splog, 'Redoing global sky subtraction'
;; skyimaget = long_skysub(sciimg, sciivar, piximg, slitmask, skymask, edgmask $
;;                         , bsp = bsp, ISLIT = ISLIT, CHK = CHK)
;; if keyword_set(ISLIT) and nct GT 0 then $
;;   skyimage[islitmask] = skyimaget[islitmask] $
;;   else skyimage=skyimaget
IF NOT keyword_set(objstruct) and nct eq 0 THEN BEGIN
    splog, 'WARNING: no objects found in this frame'
    ;;----------
    ;; Write sky subtracted image to file
    splog, 'Writing FITS file ', scifile
    mwrfits, float(sciimg), scifile, scihdr[*, 0], /create
    mwrfits, float(sciivar)*float(slitmask GT 0), scifile
    mwrfits, float(skyimage)*float(slitmask GT 0), scifile
    mwrfits, float(0*sciimg)*float(slitmask GT 0), scifile
    mwrfits, float(0*sciimg)*float(slitmask GT 0), scifile
    mwrfits, {junk:0.0}, scifile
    splog, 'Compressing ', scifile
    spawn, 'gzip -f '+scifile
    RETURN
 ENDIF 
IF NOT keyword_set(objstruct) and nct gt 0 THEN BEGIN
    splog, 'No objects found in this slit'
    splog, 'Writing FITS file ', scifile
    mwrfits, float(sciimg), scifile, scihdr[*, 0], /create
    mwrfits, float(modelivar)*float(slitmask GT 0), scifile
    mwrfits, float(skyimage)*float(slitmask GT 0), scifile
    mwrfits, float(objimage)*float(slitmask GT 0), scifile
    mwrfits, float(outmask)*float(slitmask GT 0), scifile
    mwrfits, final_struct, scifile

    splog, 'Compressing ', scifile
    spawn, 'gzip -f '+scifile
    long_plotsci, scifile, hard_ps = repstr(scifile, '.fits', '.ps')

    RETURN
 ENDIF 
;; Keep spatial flexure shift in objstruct
IF NOT KEYWORD_SET(NOSHIFT) AND KEYWORD_SET(xshift) $
  THEN objstruct.FLX_SHFT_SPA = xshift
tfinal_struct = 0
;;----------
;; Loop over each slit
orig_sky = skyimage
;; expand slit edges
traceset2xy, tset_slits[0], yy1, xx1
traceset2xy, tset_slits[1], yy2, xx2
modelivart = sciivar
;;outmaskt = fltarr(nx, ny)
outmaskt = (sciivar GT 0 AND slitmask GT 0)


IF KEYWORD_SET(islit) and KEYWORD_SET(objstruct) THEN BEGIN
    IF islit GT nslit THEN message $
      , 'ERROR: islit not found. islit cannot be larger than nslit'
    nreduce = 1
    slit_vec = [islit]
ENDIF ELSE BEGIN
    nreduce = nslit
    slit_vec = lindgen(nslit) + 1L
ENDELSE

IF KEYWORD_SET(PROF_NSIGMA) THEN BEGIN
   prof_nsigma_vec=fltarr(n_elements(objstruct))
   IF prof_nsigma LE 0.0 THEN  prof_nsigma_vec[*] = abs(PROF_NSIGMA) $
   ELSE BEGIN
      ;; Only used a dilated nsigma profile for the brightest object
      fbri = max(objstruct.PEAKFLUX, ibri)
      prof_nsigma_vec[ibri] = PROF_NSIGMA
   ENDELSE
ENDIF

FOR jj = 0L, nreduce-1L DO BEGIN
    slitid = slit_vec[jj]
    ;;  Perhaps box_rad should use FWHM estimate from objstruct?
    ii = where(objstruct.slitid EQ slitid, nobj)
    thismask = (slitmask EQ slitid)
    if nobj EQ 0 then begin
        splog, 'No objects, Skipping slit #', slitid
        continue
     endif
    fwhm_init = objstruct.FWHM
    extract_struct = long_localskysub(sciimg, sciivar, skyimage $
                                      , piximg, waveimg, ximg $
                                      , objstruct, thismask $
                                      , xx1[*, slitid-1L], xx2[*, slitid-1L] $
                                      , edgmask, bsp = bsp $
                                      , objimage = objimaget $
                                      , niter = niter $
                                      , modelivar = modelivart $
                                      , outmask = outmaskt $
                                      , indx = ii, nccd = nccd $
                                      , prof_nsigma = prof_nsigma_vec $
                                      , box_rad = box_rad $
                                      , sigrej = sigrej $
                                      , STD = STD $
                                      , scihdr = scihdr $
                                      , SN_GAUSS = SN_GAUSS $
                                      , CHK = CHK, NOLOCAL=nolocal $
                                      , skysample = skysample $
                                      , COADD_2D = COADD_2D)
    ;; Did the FWHM's change significantly? If they did the
    ;; extraction regions were wrong, so do one more iteration
    fwhm_diff = abs(fwhm_init - objstruct.FWHM)/fwhm_init
    ichange = WHERE(fwhm_diff GT 0.70d, nchange)
    IF nchange GT 0 THEN BEGIN
       splog, 'FWHM changed by more than 70%, redoing extraction with updated regions'
       ;; Start over again with original sky
       skyimage = orig_sky
       extract_struct = long_localskysub(sciimg, sciivar, skyimage $
                                      , piximg, waveimg, ximg $
                                      , objstruct, thismask $
                                      , xx1[*, slitid-1L], xx2[*, slitid-1L] $
                                      , edgmask, bsp = bsp $
                                      , objimage = objimaget $
                                      , niter = niter $
                                      , modelivar = modelivart $
                                      , outmask = outmaskt $
                                      , indx = ii, nccd = nccd $
                                      , prof_nsigma = prof_nsigma_vec $
                                      , box_rad = box_rad $
                                      , sigrej = sigrej $
                                      , STD = STD $
                                      , scihdr = scihdr $
                                      , SN_GAUSS = SN_GAUSS $
                                      , CHK = CHK, NOLOCAL=nolocal $
                                      , skysample = skysample $ 
                                      , COADD_2D = COADD_2D)

    ENDIF
    tfinal_struct = struct_append(tfinal_struct, extract_struct)
ENDFOR
stop ;check spectra for dips (i=16)

; EWR this replaces their per-slit model with the orignal model. 
skyimage = orig_sky
;; Save
if keyword_set(ISLIT) and nct GT 0 then begin 
    objimage[islitmask] = objimaget[islitmask] 
    outmask[islitmask] = outmaskt[islitmask] 
    modelivar[islitmask] = modelivart[islitmask]
endif else begin 
    objimage=objimaget
    outmask=outmaskt
    modelivar=modelivart
endelse
    
xpix = findgen(ny)
arc_struct = replicate(create_struct('ARC_FWHM_FIT', fltarr(ny) $
                                     , 'ARC_FWHM_MED', 0.0D $
                                     , 'PIX_RES',  0.0D $
                                     , 'BINNING', bin ) $
                       , n_elements(tfinal_struct))
tfinal_struct = struct_addtags(tfinal_struct, arc_struct)
stop ;check spectra
IF KEYWORD_SET(fwhmset) THEN BEGIN
    FOR slitid = 1L, nslit DO BEGIN
        slit_inds = WHERE(tfinal_struct.SLITID EQ slitid, n_obj)
        IF n_obj GT 0 AND size(fwhmset, /tname) EQ 'STRUCT' THEN BEGIN
            traceset2xy, fwhmset[slitid-1L], xpix, fwhmfit
            tfinal_struct[slit_inds].ARC_FWHM_FIT = fwhmfit
            tfinal_struct[slit_inds].ARC_FWHM_MED = fwhmset[slitid-1L].MEDIAN
        ENDIF
    ENDFOR
    long_calc_resln, tfinal_struct, anamorph = anamorph
ENDIF
stop ;check spectra
;; Convert wavelengths to vacuum
nstr = n_elements(tfinal_struct)
FOR ii = 0, nstr-1 DO BEGIN
    wvo = tfinal_struct[ii].WAVE_OPT
    wvb = tfinal_struct[ii].WAVE_BOX
    airtovac, wvo
    airtovac, wvb
    tfinal_struct[ii].WAVE_OPT = wvo
    tfinal_struct[ii].WAVE_BOX = wvb
ENDFOR

;; Tweak vacuum wavelengths to remove flexure
IF KEYWORD_SET(SKYFILE) and not keyword_set(NOFLEX) THEN begin 
    QAFILE = repstr(scifile, '.fits', '-flex.ps')
    if keyword_set(MAXFLEX) then splog, 'long_flexure: MAXFLEX = ', maxflex
    long_flexure, tfinal_struct, skyfile, flg_skyfile $
                  , QAFILE = QAFILE, MAXFLEX = maxflex $
                  , MAXGOOD=MAXGOOD
 ENDIF

;; Compute heliocentric correction 
IF NOT KEYWORD_SET(NOHELIO) THEN long_helio, scihdr, tfinal_struct

if keyword_set(ISLIT) and nct GT 0 then begin 
    msk = replicate(1B, n_elements(final_struct))
    ;; Zero out old extraction
    idx = where(final_struct.slitid EQ ISLIT, nmtch)
    if nmtch GT 0 then msk[idx] = 0B
    ;; Save good ones
    final_struct = final_struct[where(msk)]
    ;; Update ObjID
    n_newobj = n_elements(tfinal_struct)
    tfinal_struct.objid = 999L + lindgen(n_newobj)
    final_struct = struct_append(final_struct, tfinal_struct)
    srt = sort(final_struct.slitid)
    final_struct = final_struct[srt]
    final_struct.objid = lindgen(n_elements(final_struct)) + 1
endif else final_struct = tfinal_struct

;----------
; Write output file
save,file='tmpdata.sav',sciimg,modelivar,skyimage,slitmask,objimage,outmask,final_struct
for k = 0,n_elements(final_struct)-1 do $
   *final_struct[k].fluxmodel = double(*final_struct[k].fluxmodel)

splog, 'Writing FITS file ', scifile
mwrfits, float(sciimg), scifile, scihdr, /create
mwrfits, float(modelivar)*float(slitmask GT 0), scifile
mwrfits, float(skyimage)*float(slitmask GT 0), scifile
mwrfits, float(objimage)*float(slitmask GT 0), scifile
mwrfits, float(outmask)*float(slitmask GT 0), scifile
mwrfits, final_struct, scifile


splog, 'Compressing ', scifile
spawn, 'gzip -f '+scifile
save,file=repstr(scifile,'fits','sav'),final_struct

;long_plotsci, scifile, hard_ps = repstr(scifile, '.fits', '.ps')
splog, 'Elapsed time = ', systime(1)-t0, ' sec'

RETURN
END
