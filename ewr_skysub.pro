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
                      , bsp = bsp_in, islit = islit, CHK = CHK, $
                     waveimg = waveimg,wavemask=wavemask,$
                     ximg=ximg, nudgelam = nudgelam
   IF NOT KEYWORD_SET(BSP_in) THEN BSP_in = 0.6D
   bsp = bsp_in
   IF NOT KEYWORD_SET(SIGREJ) THEN SIGREJ = 3.0
   if n_elements(nudgelam) eq 0 then nudgelam = 0.0 
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

      sind = sort(wsky)
      nElts = n_elements(sind) 
      sky_out = fltarr(nElts)
      sky_est1 = fltarr(nElts)
      sky_est2 = fltarr(nElts)
      sky_mad = fltarr(nElts)
      sky_min = fltarr(nElts)
      sky_2sig = fltarr(nElts)
      sky_1sig = fltarr(nElts)
      sky_50 = fltarr(nElts)
      sky_75 = fltarr(nElts)
      sky_95 = fltarr(nElts)
;      window= ((nElts/ny)/6.) > 16
      window = (nElts/ny)*2
      madsky = mad(sky)
      mask = bytarr(nElts)+1
      nPass = 10
      for i = 0,nPass-1 do begin
         idx = where(mask)
         sky_mad = median(abs(sky[sind[idx]]-shift(sky[sind[idx]],1)),window*2)/0.6745
         sky_med = median(sky[sind[idx]],window)
         newmask = sky[sind[idx]] lt sky_mad*2+sky_med
         mask[idx] = newmask
      endfor
      idx = where(mask)
      sky_out = median(sky[sind[idx]],window/2)
; Pass 1
;;       for k = 0,nElts-1 do begin
;;          subset_inds = sind[((k-window)>0):(k+window)<(nElts-1)]
;;          subset = sky[subset_inds]
;; ;         sind2 = sort(subset)
;; ;         nElts2 = n_elements(sind2)
;; ;         sky_out[k] = subset[sind2[0.3*nElts2]]
;;          sky_min[k] = min(subset)
;;          pcts = cgPercentiles(subset,percentiles=[0.02275,0.158,0.5,0.75,0.95])
;;          sky_2sig[k] = pcts[0]
;;          sky_1sig[k] = pcts[1]
;;          sky_50[k] = pcts[2]
;;          sky_75[k] = pcts[3]
;;          sky_95[k] = pcts[4]
;;       endfor 

;; ; Pass 2
;; ;      skymad_im = bspline_iterfit(wsky,sky_mad, $
;; ;                                  everyn=255L*30,/groupbadpix,maxrej=10,$
;; ;                                  lower=sigrej,upper=sigrej,yfit=madsky)

;;       for k = 0,nElts-1 do begin
;;          subset_inds = sind[((k-window)>0):(k+window)<(nElts-1)]
;;          subset = sky[subset_inds]
;;          sky_est_regular = median(subset)
;; ;         if min(abs(wsky[sind[k]]-wavemask-nudgelam)) lt 6.0 then begin
;;          sind2 = where((subset lt (madsky*2+sky_1sig[k])) and $
;;                        (subset gt (sky_2sig[k])))
;;          sky_est_lower = median(subset[sind2])
;; ;         endif else sky_est_lower = sky_est_regular
;;          sky_est1[k] = sky_est_regular
;;          sky_est2[k] = sky_est_lower
;;          sky_out[k] = sky_est_regular < sky_est_lower
;;       endfor
      fullbkpt = bspline_bkpts(wsky[sind[idx]], nord = 4, everyn=window*0.2, /silent)

      skyset = bspline_longslit(wsky[sind[idx]], sky_out, sky_ivar[sind[idx]], isky[sind[idx]]*0.+1. $
                                , /groupbadpix, maxrej = 10 $
                                , fullbkpt = fullbkpt, upper = sigrej $
                                , lower = sigrej, /silent, yfit=yfit,everyn=window*0.2)

      ;;;;;;;;;;;;;;;;;;;
      ;; JXP -- Have had to kludge this when using a kludged Arc frame
;      skyset = bspline_iterfit(wsky, sky $ ;nvvar=sky_ivar  $
;                               , /groupbadpix, maxrej = 10 $
;                               , everyn=255L, upper = sigrej $
;                                , lower = sigrej, /silent, yfit=yfit)
      ;;;;;;;;;;;;;;;;;;;
      sky_image[all] = (bspline_valu(waveimg[all], skyset) )>0

;stop
      IF KEYWORD_SET(CHK) THEN $
         x_splot, wsky, sky, psym1 = 3, xtwo = wsky, ytwo = yfit, /block

   endfor

   return, sky_image
end
