pro ahk_fluxandcoadd_SlitsInDiffFiles, stem1,stem2,stem3, stem4, stem5, slitidmin=slitidmin, slitidmax=slitidmax, sensfuncfile=sensfuncfile, gratingstem=gratingstem

;+
; NAME:
;   ahk_fluxandcoadd
;
; PURPOSE:
;   Loops through all exposures and slits for a single grating (e.g. blue600, red900, or red400) to get fluxed and coadded spectra with errors.
;
; PURPOSE:
;   Searches the current directory for fits files beginning with strings stem1, stem2, stem3. These are expected to be the output from ahk_objextract_multigrating
;
; CALLING SEQUENCE:
;	ahk_fluxandcoadd, '../profiles/blue2082_slit', '../profiles/blue2083_slit', '../profiles/blue2086_slit', 'fake', 'fake', slitidmin=1,slitidmax=21, sensfuncfile='../../../fluxCalibration/blue600_bsens_all.fits', gratingstem='blue600'
;	ahk_fluxandcoadd, '../profiles/lred2068_slit', '../profiles/lred2069_slit', '../profiles/lred2070_slit', 'fake', 'fake', slitidmin=1,slitidmax=21, sensfuncfile='../../../fluxCalibration/red900_bsens_all.fits', gratingstem='red900'
;	ahk_fluxandcoadd, '../profiles/lred2073_slit', 'fake', 'fake', 'fake', 'fake', slitidmin=1,slitidmax=21, sensfuncfile='../../../fluxCalibration/red400_bsens_lred2064.fits', gratingstem='red400'

;
; INPUTS:
;
; OPTIONAL INPUTS:
;                
; OUTPUTS:
;
;
; OPTIONAL OUTPUTS:
;   
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;    Must choose frm_sci more wisely! currently just uses frm_sci=1,
;    but often (at least blue side) needs this to be 2.
;
; PROCEDURES CALLED:
;	
;
; REVISION HISTORY:
;   2-December-2014 -- Written by Alice Koning
;-

FOR slitid =slitidmin,slitidmax DO BEGIN
   print, 'AHK_FLUXANDCOADD: Working on slit '+STRTRIM(slitid,2)
	FOR objid = 1, 5 DO BEGIN
		;Need a way to keep track of which files are defined. Use the following counters
		c1=0
		c2=0
		c3=0		
                c4=0
                c5=0

		; Check which fits files are defined and flux calibrate with the given sensitivity function (sensfuncfile)
		IF FILE_TEST(stem1+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'.fits.gz') THEN BEGIN
			c1=1
			LONG_FLUXCAL, stem1+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'.fits.gz',SENSFUNCFILE=sensfuncfile,OUTFIL=stem1+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits', frm_sci=1
		ENDIF
		IF FILE_TEST(stem2+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'.fits.gz') THEN BEGIN
			c2=1
			LONG_FLUXCAL, stem2+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'.fits.gz',SENSFUNCFILE=sensfuncfile,OUTFIL=stem2+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits', frm_sci=1
		ENDIF
		IF FILE_TEST(stem3+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'.fits.gz') THEN BEGIN
			c3=1
			LONG_FLUXCAL, stem3+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'.fits.gz',SENSFUNCFILE=sensfuncfile,OUTFIL=stem3+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits', frm_sci=1
		ENDIF
		IF FILE_TEST(stem4+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'.fits.gz') THEN BEGIN
			c4=1
			LONG_FLUXCAL, stem4+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'.fits.gz',SENSFUNCFILE=sensfuncfile,OUTFIL=stem4+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits', frm_sci=1
		ENDIF
		IF FILE_TEST(stem5+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'.fits.gz') THEN BEGIN
			c5=1
			LONG_FLUXCAL, stem5+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'.fits.gz',SENSFUNCFILE=sensfuncfile,OUTFIL=stem5+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits', frm_sci=1
		ENDIF

		IF c1 AND c2 AND c3 AND c4 AND NOT c5 THEN BEGIN
			print, 'Object '+STRTRIM(objid,2)+' defined in first 4 exposures.'

			; Combine spectra with long_combspec
			; Note OUTFIL from LONG_FLUXCAL contain flux, error, wavelength as first three fits extensions.
			flux1=mrdfits(stem1+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',0,hd)
			err1=mrdfits(stem1+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',1,hd)
			wave1=mrdfits(stem1+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',2,hd)
			ivar1=1./err1^2

			flux2=mrdfits(stem2+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',0,hd)
			err2=mrdfits(stem2+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',1,hd)
			wave2=mrdfits(stem2+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',2,hd)
			ivar2=1./err2^2

			flux3=mrdfits(stem3+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',0,hd)
			err3=mrdfits(stem3+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',1,hd)
			wave3=mrdfits(stem3+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',2,hd)
			ivar3=1./err3^2

			flux4=mrdfits(stem4+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',0,hd)
			err4=mrdfits(stem4+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',1,hd)
			wave4=mrdfits(stem4+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',2,hd)
			ivar4=1./err4^2
	
			influx = [[flux1],[flux2],[flux3],[flux4]]
			inivar = [[ivar1],[ivar2],[ivar3],[ivar4]]
			lam = [[wave1],[wave2],[wave3],[wave4]]
			loglam = alog(lam)/alog(10)
			LONG_COMBSPEC,  influx, inivar, loglam, newloglam=newloglam, newflux=newflux, newivar=newivar

			newlam = 10^newloglam
			newerr = 1/newivar^(0.5)

			; Save flux, error, and wave as first three extensions in new fits structure
			fileout = gratingstem+'_slit'+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux_coadd.fits'
			mwrfits,newflux,fileout,/create
			mwrfits,newerr,fileout 
			mwrfits,newlam,fileout 		
                        CONTINUE
		ENDIF
		IF c1 AND c2 AND c3 AND c4 AND c5 THEN BEGIN
			print, 'Object '+STRTRIM(objid,2)+' defined in all 5 exposures.'

			; Combine spectra with long_combspec
			; Note OUTFIL from LONG_FLUXCAL contain flux, error, wavelength as first three fits extensions.
			flux1=mrdfits(stem1+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',0,hd)
			err1=mrdfits(stem1+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',1,hd)
			wave1=mrdfits(stem1+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',2,hd)
			ivar1=1./err1^2

			flux2=mrdfits(stem2+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',0,hd)
			err2=mrdfits(stem2+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',1,hd)
			wave2=mrdfits(stem2+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',2,hd)
			ivar2=1./err2^2

			flux3=mrdfits(stem3+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',0,hd)
			err3=mrdfits(stem3+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',1,hd)
			wave3=mrdfits(stem3+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',2,hd)
			ivar3=1./err3^2

			flux4=mrdfits(stem4+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',0,hd)
			err4=mrdfits(stem4+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',1,hd)
			wave4=mrdfits(stem4+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',2,hd)
			ivar4=1./err4^2

			flux5=mrdfits(stem5+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',0,hd)
			err5=mrdfits(stem5+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',1,hd)
			wave5=mrdfits(stem5+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',2,hd)
			ivar5=1./err5^2
	
			influx = [[flux1],[flux2],[flux3],[flux4],[flux5]]
			inivar = [[ivar1],[ivar2],[ivar3],[ivar4],[ivar5]]
			lam = [[wave1],[wave2],[wave3],[wave4],[wave5]]
			loglam = alog(lam)/alog(10)
			LONG_COMBSPEC,  influx, inivar, loglam, newloglam=newloglam, newflux=newflux, newivar=newivar

			newlam = 10^newloglam
			newerr = 1/newivar^(0.5)

			; Save flux, error, and wave as first three extensions in new fits structure
			fileout = gratingstem+'_slit'+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux_coadd.fits'
			mwrfits,newflux,fileout,/create
			mwrfits,newerr,fileout 
			mwrfits,newlam,fileout 		
                        CONTINUE
		ENDIF


		IF c1 AND c2 AND c3 AND NOT c4 AND NOT c5 THEN BEGIN
			print, 'Object '+STRTRIM(objid,2)+' defined in first 3 exposures.'
                        
			; Combine spectra with long_combspec
			; Note OUTFIL from LONG_FLUXCAL contain flux, error, wavelength as first three fits extensions.
			flux1=mrdfits(stem1+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',0,hd)
			err1=mrdfits(stem1+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',1,hd)
			wave1=mrdfits(stem1+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',2,hd)
			ivar1=1./err1^2

			flux2=mrdfits(stem2+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',0,hd)
			err2=mrdfits(stem2+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',1,hd)
			wave2=mrdfits(stem2+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',2,hd)
			ivar2=1./err2^2

			flux3=mrdfits(stem3+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',0,hd)
			err3=mrdfits(stem3+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',1,hd)
			wave3=mrdfits(stem3+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',2,hd)
			ivar3=1./err3^2

			influx = [[flux1],[flux2],[flux3]]
			inivar = [[ivar1],[ivar2],[ivar3]]
			lam = [[wave1],[wave2],[wave3]]
			loglam = alog(lam)/alog(10)
			LONG_COMBSPEC,  influx, inivar, loglam, newloglam=newloglam, newflux=newflux, newivar=newivar, /NOSHARP

			newlam = 10^newloglam
			newerr = 1/newivar^(0.5)

			; Save flux, error, and wave as first three extensions in new fits structure
			fileout = gratingstem+'_slit'+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux_coadd.fits'
			mwrfits,newflux,fileout,/create
			mwrfits,newerr,fileout 
			mwrfits,newlam,fileout 		
                        CONTINUE
		ENDIF

		IF c1 AND c2 AND NOT c3 AND NOT c4 AND NOT c5 THEN BEGIN
			print, 'Object '+STRTRIM(objid,2)+' defined in 2 of 4 exposures.'

			; Combine spectra with long_combspec
			; Note OUTFIL from LONG_FLUXCAL contain flux, error, wavelength as first three fits extensions.
			flux1=mrdfits(stem1+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',0,hd)
			err1=mrdfits(stem1+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',1,hd)
			wave1=mrdfits(stem1+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',2,hd)
			ivar1=1./err1^2

			flux2=mrdfits(stem2+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',0,hd)
			err2=mrdfits(stem2+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',1,hd)
			wave2=mrdfits(stem2+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',2,hd)
			ivar2=1./err2^2
	
			influx = [[flux1],[flux2]]
			inivar = [[ivar1],[ivar2]]
			lam = [[wave1],[wave2]]
			loglam = alog(lam)/alog(10)
			LONG_COMBSPEC,  influx, inivar, loglam, newloglam=newloglam, newflux=newflux, newivar=newivar

			newlam = 10^newloglam
			newerr = 1/newivar^(0.5)

			; Save flux, error, and wave as first three extensions in new fits structure
			fileout = gratingstem+'_slit'+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux_coadd.fits'
			mwrfits,newflux,fileout,/create
			mwrfits,newerr,fileout 
			mwrfits,newlam,fileout 
                        CONTINUE
		ENDIF

		IF c1 AND c3 AND NOT c2 AND NOT c4 AND NOT c5 THEN BEGIN
			print, 'Object '+STRTRIM(objid,2)+' defined in 2 of 4 exposures.'

			; Combine spectra with long_combspec
			; Note OUTFIL from LONG_FLUXCAL contain flux, error, wavelength as first three fits extensions.
			flux1=mrdfits(stem1+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',0,hd)
			err1=mrdfits(stem1+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',1,hd)
			wave1=mrdfits(stem1+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',2,hd)
			ivar1=1./err1^2

			flux3=mrdfits(stem3+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',0,hd)
			err3=mrdfits(stem3+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',1,hd)
			wave3=mrdfits(stem3+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',2,hd)
			ivar3=1./err3^2
	
			influx = [[flux1],[flux3]]
			inivar = [[ivar1],[ivar3]]
			lam = [[wave1],[wave3]]
			loglam = alog(lam)/alog(10)
			LONG_COMBSPEC,  influx, inivar, loglam, newloglam=newloglam, newflux=newflux, newivar=newivar

			newlam = 10^newloglam
			newerr = 1/newivar^(0.5)

			; Save flux, error, and wave as first three extensions in new fits structure
			fileout = gratingstem+'_slit'+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux_coadd.fits'
			mwrfits,newflux,fileout,/create
			mwrfits,newerr,fileout 
			mwrfits,newlam,fileout 
                        CONTINUE
		ENDIF

		IF c2 AND c3 AND NOT c1 AND NOT c4 AND NOT c5 THEN BEGIN
			print, 'Object '+STRTRIM(objid,2)+' defined in 2 of 4 exposures.'

			; Combine spectra with long_combspec
			; Note OUTFIL from LONG_FLUXCAL contain flux, error, wavelength as first three fits extensions.
			flux2=mrdfits(stem2+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',0,hd)
			err2=mrdfits(stem2+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',1,hd)
			wave2=mrdfits(stem2+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',2,hd)
			ivar2=1./err2^2

			flux3=mrdfits(stem3+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',0,hd)
			err3=mrdfits(stem3+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',1,hd)
			wave3=mrdfits(stem3+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',2,hd)
			ivar3=1./err3^2
	
			influx = [[flux2],[flux3]]
			inivar = [[ivar2],[ivar3]]
			lam = [[wave2],[wave3]]
			loglam = alog(lam)/alog(10)
			LONG_COMBSPEC,  influx, inivar, loglam, newloglam=newloglam, newflux=newflux, newivar=newivar

			newlam = 10^newloglam
			newerr = 1/newivar^(0.5)

			; Save flux, error, and wave as first three extensions in new fits structure
			fileout = gratingstem+'_slit'+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux_coadd.fits'
			mwrfits,newflux,fileout,/create
			mwrfits,newerr,fileout 
			mwrfits,newlam,fileout
                        CONTINUE
		ENDIF

		IF c1 AND NOT c2 AND NOT c3 AND NOT c4 THEN BEGIN
			print, 'Object '+STRTRIM(objid,2)+' defined in 1 of 4 exposures. Skipping to next object.'
			CONTINUE
		ENDIF

		IF c2 AND NOT c1 AND NOT c3 AND NOT c4 THEN BEGIN
			print, 'Object '+STRTRIM(objid,2)+' defined in 1 of 4 exposures. Skipping to next object.'
			CONTINUE
		ENDIF

		IF c3 AND NOT c1 AND NOT c2 AND NOT c4 THEN BEGIN
			print, 'Object '+STRTRIM(objid,2)+' defined in 1 of 4 exposures. Skipping to next object.'
			CONTINUE
		ENDIF
		
		IF NOT c1 AND NOT c2 AND NOT c3 AND NOT c4 THEN CONTINUE

	ENDFOR
ENDFOR
end
