;+
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
;
;       Wed Jul 22 10:29:14 2015, <erosolo@noise.siglab.ok.ubc.ca>
;
;		Additional revisions to deal with diffuse emission
;
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
   sky_image2 = sky_image
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

; Generate an all-mask sky model.

   nbd = where((slitmask gt 0)*(waveimg ne 0), nnbd)
   wsky_nbd = waveimg[nbd]
   sky_nbd = sciimg[nbd]
   sky_ivar_nbd = sciivar[nbd]

   sindn = sort(wsky_nbd)
   nElts = n_elements(sindn) 
   madsky = mad(sky_nbd)
   mask = bytarr(nElts)+1
   nPass = 10
   nwindow = (nElts/ny)*0.1 
   for i = 0,nPass-1 do begin
      idxn = where(mask)
      sky_mad = median(abs(sky_nbd[sindn[idxn]]-shift(sky_nbd[sindn[idxn]],1)),nwindow*2)/0.6745
      sky_med = median(sky_nbd[sindn[idxn]],nwindow)
      newmask = sky_nbd[sindn[idxn]] lt sky_mad*2+sky_med
      mask[idxn] = newmask
   endfor
   idxn = where(mask)
   sky_out_nbd = median(sky_nbd[sindn[idxn]],nwindow/2)
   nbd_fullbkpt = bspline_bkpts(wsky_nbd[sindn[idxn]], $
                                nord = 4, everyn=nwindow*0.2, /silent)
   
   nbd_skyset = bspline_longslit(wsky_nbd[sindn[idxn]], sky_out_nbd, $
                                 sky_ivar_nbd[sindn[idxn]], nbd[sindn[idxn]]*0.+1. $
                                 , /groupbadpix, maxrej = 10 $
                                 , fullbkpt = nbd_fullbkpt, upper = sigrej $
                                 , lower = sigrej, /silent, $
                                 yfit=yfit,everyn=nwindow*0.2)
   sky_image2[nbd] = (bspline_valu(waveimg[nbd], nbd_skyset) )>0
   
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



      fullbkpt = bspline_bkpts(wsky[sind[idx]], nord = 4, $
                               everyn=window*0.2, /silent)

      skyset = bspline_longslit(wsky[sind[idx]], sky_out, $
                                sky_ivar[sind[idx]], isky[sind[idx]]*0.+1. $
                                , /groupbadpix, maxrej = 10 $
                                , fullbkpt = fullbkpt, upper = sigrej $
                                , lower = sigrej, /silent, $
                                yfit=yfit,everyn=window*0.2)



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
   sky_image = (sky_image < sky_image2)
   return, sky_image
end
