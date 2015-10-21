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
function long_skysub, sciimg, sciivar, piximg, slitmask, skymask, edgmask $
                      , subsample = subsample, npoly = npoly $
                      , nbkpts = nbkpts $
                      , bsp = bsp, islit = islit, CHK = CHK

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
      all = where(slitmask EQ slitid, nall)
      isky = where(sky_slitmask EQ slitid, nsky)
      if (nsky LT 10) then begin
         splog, 'Not enough sky pixels found in slit ', slitid, nsky, nall
         continue
      endif
      
      isky = isky[sort(piximg[isky])]
      wsky = piximg[isky]
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

;      bsp = 1.5
      sigreg=1.5
      sind = sort(wsky)
      mask = 1b 
      nPass = 10

      window = n_elements(wsky)/nx 
      mask = bytarr(n_elements(wsky))+1b 
      for i = 0,nPass-1 do begin
         idx = where(mask)
         sky_mad = median(abs(sky[sind[idx]]-$
                              shift(sky[sind[idx]],1)),window*2)/0.6745
         sky_med = median(sky[sind[idx]],window)
         newmask = sky[sind[idx]] lt sky_mad*2+sky_med
         mask[idx] = newmask
      endfor

      sigrej = 1.0
      sigma = interpol(sky_mad,wsky[idx],wsky)
      spacepts = (n_elements(wsky)/ny)*1.0
      fullbkpt = bspline_bkpts(wsky, nord = 4, $
                               /silent, $
                               everyn=spacepts)
      ;; fullbkpt = bspline_bkpts(wsky, nord = 4, $
      ;;                          /silent, $
      ;;                          bkspace=0.6)


      ivar = sky_ivar < 1/sigma^2
      numiter=10

      skyset = bspline_longslit(wsky, sky, ivar, ivar*0.+1. $
                                , /groupbadpix, maxrej = 3 $
                                , fullbkpt = fullbkpt, upper = sigrej $
                                , lower = sigrej, /silent, yfit=yfit $
                                , maxiter=50, numiter=numiter,everyn=spacepts)
      ;; skyset = bspline_longslit(wsky, sky, ivar, ivar*0.+1. $
      ;;                           , /groupbadpix, maxrej = 3 $
      ;;                           , fullbkpt = fullbkpt, upper = sigrej $
      ;;                           , lower = sigrej, /silent, yfit=yfit $
      ;;                           , maxiter=50, numiter=numiter,bkspace=bsp)
     



      ;;;;;;;;;;;;;;;;;;;
      ;; JXP -- Have had to kludge this when using a kludged Arc frame
;      skyset = bspline_iterfit(wsky, sky $ ;nvvar=sky_ivar  $
;                               , /groupbadpix, maxrej = 10 $
;                               , everyn=255L, upper = sigrej $
;                                , lower = sigrej, /silent, yfit=yfit)
      ;;;;;;;;;;;;;;;;;;;
      sky_image[all] = bspline_valu(piximg[all], skyset) 
      IF KEYWORD_SET(CHK) THEN $
         x_splot, wsky, sky, psym1 = 3, xtwo = wsky, ytwo = yfit, /block
   endfor
;   stop
   return, sky_image
end
