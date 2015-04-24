; NAME:
;   long_skysub
;
; PURPOSE:
;
; CALLING SEQUENCE:
; skyimage = long_skysub( sciimg, sciivar, piximg, slitmask, skymask) 
;    
; INPUTS:
;  sciimg -- Science image
;  sciivar -- Inverse variance
;  piximg -- 
;  slitmask -- 
;  skymask -- 
;
; OUTPUTS: 
;  skyimage -- Sky model
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
;   27-May-2005 Written by J. Hennawi (UCB)
;-
;------------------------------------------------------------------------------
function ewr_skysub, sciimg, sciivar, piximg, slitmask, skymask, edgmask $
                      , subsample = subsample, npoly = npoly $
                      , nbkpts = nbkpts $
                      , bsp = bsp, islit = islit, CHK = CHK, $
                     waveimg = waveimg,wavemask=wavemask,$
                     ximg=ximg

   IF NOT KEYWORD_SET(BSP) THEN BSP = 0.6D
   IF NOT KEYWORD_SET(SIGREJ) THEN SIGREJ = 3.0
   nx = (size(sciimg))[1] 
   ny = (size(sciimg))[2] 
   nslit = max(slitmask)
   sky_image = sciimg * 0.
   if NOT keyword_set(skymask) then skymask = slitmask*0+1
   
   sky_slitmask = slitmask*(skymask*(sciivar GT 0) AND EDGMASK EQ 0)


   IF KEYWORD_SET(islit) THEN BEGIN
       IF islit GT nslit THEN message $
         , 'ERROR: islit not found. islit cannot be larger than nslit'
       nreduce = 1
       slit_vec = islit
   ENDIF ELSE BEGIN
       nreduce = nslit
       slit_vec = lindgen(nslit) + 1L
   ENDELSE
   for jj = 0L, nreduce-1L DO BEGIN
      slitid = slit_vec[jj]
      ;; Select only the pixels on this slit
      all = where((slitmask EQ slitid)*(waveimg ne 0), nall)
      isky = where((sky_slitmask EQ slitid)*(waveimg ne 0), nsky)
      if (nsky LT 10) then begin
         splog, 'Not enough sky pixels found in slit ', slitid, nsky, nall
         continue
      endif
      
      isky = isky[sort(waveimg[isky])]
      wsky = waveimg[isky]
      sky = sciimg[isky]
      sky_ivar = sciivar[isky]
      
;      if keyword_set(nbkpts) then everyn =  1.0*nsky / (nbkpts+1) $
;      else everyn = 0.6 * nsky / ny
      pos_sky = where(sky GT 1.0 AND sky_ivar GT 0., npos)
      if npos GT ny then begin
         lsky = alog(sky[pos_sky])
         lsky_ivar = lsky*0.+0.1
         
         skybkpt = bspline_bkpts(wsky[pos_sky], nord = 4, bkspace = bsp $
                                 , /silent)
         lskyset = bspline_longslit(wsky[pos_sky], lsky, lsky_ivar, $
                                    pos_sky*0.+1, fullbkpt = skybkpt $
                                    , upper = sigrej, lower = sigrej $
                                    , /silent, yfit = lsky_fit $
                                    , /groupbadpix)
         res = (sky[pos_sky] - exp(lsky_fit))*sqrt(sky_ivar[pos_sky])
         lmask = (res LT 5.0 AND res GT -4.0)
         sky_ivar[pos_sky] = sky_ivar[pos_sky] * lmask
      endif
      fullbkpt = bspline_bkpts(wsky, nord = 4, bkspace = bsp, /silent)

; Run a 30%-ile rolling filter across the data
      sind = sort(wsky)
      nElts = n_elements(sind) 
      sky_out = fltarr(nElts)
      sky_mad = fltarr(nElts)
      sky_min = fltarr(nElts)
;      window= ((nElts/ny)/6.) > 16
      window = nElts/ny
      madsky = mad(sky)
; Pass 1
      for k = 0,nElts-1 do begin
         subset_inds = sind[((k-window)>0):(k+window)<(nElts-1)]
         subset = sky[subset_inds]
;         sind2 = sort(subset)
;         nElts2 = n_elements(sind2)
;         sky_out[k] = subset[sind2[0.3*nElts2]]
         sky_mad[k] = mad(subset)
         sky_min[k] = min(subset)
      endfor 
; Pass 2

      skymad_im = bspline_iterfit(wsky,sky_mad, $
                                  everyn=255L*30,/groupbadpix,maxrej=10,$
                                  lower=sigrej,upper=sigrej,yfit=madsky)

      for k = 0,nElts-1 do begin
         subset_inds = sind[((k-window)>0):(k+window)<(nElts-1)]
         subset = sky[subset_inds]
         if min(abs(wsky[k]-wavemask)) lt 5.5 then begin
            sind2 = where(subset lt sky_min[k]+3*madsky[k])
            sky_out[k] = median([subset[sind2]])
         endif else sky_out[k] = median([subset])
      endfor
stop
      skyset = bspline_longslit(wsky[sind], sky_out, sky_ivar[sind], isky*0.+1. $
                                , /groupbadpix, maxrej = 10 $
                                , fullbkpt = fullbkpt, upper = sigrej $
                                , lower = sigrej, /silent, yfit=yfit)

      ;;;;;;;;;;;;;;;;;;;
      ;; JXP -- Have had to kludge this when using a kludged Arc frame
;      skyset = bspline_iterfit(wsky, sky $ ;nvvar=sky_ivar  $
;                               , /groupbadpix, maxrej = 10 $
;                               , everyn=255L, upper = sigrej $
;                                , lower = sigrej, /silent, yfit=yfit)
      ;;;;;;;;;;;;;;;;;;;
      sky_image[all] = (bspline_valu(waveimg[all], skyset) )>0
      IF KEYWORD_SET(CHK) THEN $
         x_splot, wsky, sky, psym1 = 3, xtwo = wsky, ytwo = yfit, /block

   endfor

   return, sky_image
end
