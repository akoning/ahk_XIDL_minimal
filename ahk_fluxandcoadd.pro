pro ahk_fluxandcoadd, listOfFiles, sensfuncfile=sensfuncfile, gratingstem=gratingstem
;pro ahk_fluxandcoadd, stem1,stem2,stem3, stem4, stem5, slitidmin=slitidmin, slitidmax=slitidmax, sensfuncfile=sensfuncfile, gratingstem=gratingstem

;+
; NAME:
;   ahk_fluxandcoadd
;
; PURPOSE:
;   Loops through all exposures and slits for a single grating (e.g. blue600, red900, or red400) to get fluxed and coadded spectra with errors.
;
;   Searches the current directory for fits files beginning with
;   strings stem1, stem2, stem3. These are expected to be the output
;   from ahk_objextract_multigrating
;
;   Changed from ahk_fluxandcoadd_SlitsInDiffFiles.pro to handle
;   output from ahk_objextract_multigrating_multiexposure.pro where
;   all the slits of a given exposure are in the same output file
;   (array of structures in last fits extension).
;
; CALLING SEQUENCE:
;       list = ['../blue6099_allSlitsReduced.fits.gz', '../blue6100_allSlitsReduced.fits.gz', '../blue6103_allSlitsReduced.fits.gz']
;	ahk_fluxandcoadd, list, sensfuncfile='../../../fluxCalibration/blue600_bsens_all.fits', gratingstem='blue600'
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
;   
;
; PROCEDURES CALLED:
;	
;
; REVISION HISTORY:
;   23-June-2015 -- Written by Alice Koning
;-
nFiles=size(listOfFiles,/N_ELE)
nslits=0
nobjs=0

structs = ptrarr(nFiles, /allocate_heap)
FOREACH element, listOfFiles, key DO BEGIN
   *structs[key] = mrdfits(element,5,hd)
   nslits = nslits > max((*structs[key]).slit_id)
   nobjs = nobjs > max((*structs[key]).obj_id)
ENDFOREACH

FOR slitid =1,nslits DO BEGIN
   print, 'AHK_FLUXANDCOADD: Working on slit '+STRTRIM(slitid,2)
	FOR objid = 1, nobjs DO BEGIN
           ngood = 0
           tempflux = []
           tempivar = []
           templam = []
           FOREACH element, listOfFiles, key DO BEGIN
              ;Find where each structure has this slitid and objid
              index = where((*structs[key]).slit_id EQ slitid and (*structs[key]).obj_id eq objid, count)
              ngood += count
              IF count LT 1 THEN CONTINUE

              ;Do the flux calibration
              LONG_FLUXCAL_SLITINDEX, element, slitindex=index, SENSFUNCFILE=sensfuncfile, OUTFIL=strmid(element,0,8)+'_slit'+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits', frm_sci=1

              ;Set things up so we can coadd easily after all good exposures are flux calibrated
              flux1=mrdfits(strmid(element,0,8)+'_slit'+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',0,hd)
              err1=mrdfits(strmid(element,0,8)+'_slit'+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',1,hd)
              wave1=mrdfits(strmid(element,0,8)+'_slit'+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits',2,hd)
              ivar1=1./err1^2

              tempflux = [tempflux, flux1]
              tempivar = [tempivar, ivar1]
              templam = [templam, wave1]
              nspec = size(flux1,/dim)
           ENDFOREACH

          IF ngood LT 1 THEN CONTINUE
           IF ngood EQ 1 THEN BEGIN
              print, 'Only one good exposure for slit '+STRTRIM(slitid,2)+' Object '+STRTRIM(objid,2)+'. No need to coadd. Continuing.'
              spawn, 'cp '+ strmid(element,0,8)+'_slit'+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux.fits '+gratingstem+'_slit'+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux_coadd.fits'
              CONTINUE
           ENDIF

           print, STRTRIM(ngood,2)+' good exposures for slit '+STRTRIM(slitid,2)+' Object '+STRTRIM(objid,2)

           ;Reshape flux, ivar, lam to have dimensions of [nspec,nimgs]
           flux = reform(tempflux, [nspec,ngood])
           ivar =  reform(tempivar, [nspec,ngood])
           lam = reform(templam, [nspec,ngood])

           loglam = alog(lam)/alog(10)
           
           LONG_COMBSPEC,  flux, ivar, loglam, newloglam=newloglam, newflux=newflux, newivar=newivar, /NOSHARP

           newlam = 10^newloglam
           newerr = 1/newivar^(0.5)

                                ; Save flux, error, and wave as first three extensions in new fits structure
           fileout = gratingstem+'_slit'+STRTRIM(slitid,2)+'_obj'+STRTRIM(objid,2)+'_flux_coadd.fits'
           mwrfits,newflux,fileout,/create
           mwrfits,newerr,fileout 
           mwrfits,newlam,fileout 		

       ENDFOR
ENDFOR
end
