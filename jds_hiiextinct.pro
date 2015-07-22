;+
; NAME:
;   jds_hiiextinct
;
; PURPOSE:
;   Given a set of emission line fluxes from an HII region, determine 
;   the equivalent width of the stellar Balmer absorption lines and the
;   extinction
;
; CALLING SEQUENCE:
;   jds_hiiextinct
;
; INPUTS:
;   frame           - CCD exposure number
;   slit            - slit number on the mask
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   chb             - Best-fit value for C(Hbeta)
;   dchb            - uncertainty on C(Hbeta)
;   eqw             - Best-fit equivalent width for stellar Balmer absorption
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;                   - FIND A WAY TO MAKE SURE THAT THE FIT CONVERGES!!!!!
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   24-Aug-2006  Converted from jds_hiiextinct_old by J. Simon, Caltech
;
;------------------------------------------------------------------------------
pro jds_hiiextinct,w,spec,dspec,linename,linewave,flux,fluxerr,chb,d_chb,eqw, $
                   d_eqw,d_chbmc,d_eqwmc,LRISBLUE=LRISBLUE,FITRANGE=fitrange, $
                   SILENT=SILENT,PN=PN,PLOT=PLOT,VARIABLE_EW=VARIABLE_EW


resolve_routine,'extinctfunc',/is_func


;GET BALMER LINE POSITIONS
hb = where(linename eq '         Hb',got_hbeta)
hg = where(linename eq '         Hg',got_hgamma)
hd = where(linename eq '         H6',got_hdelta)
h9 = where(linename eq '         H9',got_h9)
h10 = where(linename eq '         H10',got_h10)
h11 = where(linename eq '         H11',got_h11)
h12 = where(linename eq '         H12',got_h12)

;IF HBETA/HGAMMA AREN'T THERE, WE CAN'T DERIVE EXTINCTION
if got_hbeta eq 0 then stop,'No flux for Hbeta!'
if got_hgamma eq 0 then stop,'No flux for Hgamma!'

;GET BALMER LINE FLUXES
f_hb = flux[hb]
df_hb = fluxerr[hb]
f_hg = flux[hg]
df_hg = fluxerr[hg]
lambda = [linewave[hb],linewave[hg]]

;NOW SEE WHAT OTHER LINES WE CAN USE
if got_hdelta eq 1 then begin
    if flux[hd] ge 50.*fluxerr[hd] then begin
        f_hd = flux[hd]
        df_hd = fluxerr[hd]
        lambda = [lambda,linewave[hd]]
        print,'Strong detection of Hdelta . . .'
    endif else got_hdelta = 0
endif
if got_h9 eq 1 then begin
     if flux[h9] ge 20.*fluxerr[h9] then begin
        f_h9 = flux[h9]
        df_h9 = fluxerr[h9]
        lambda = [lambda,linewave[h9]]
        print,'Strong detection of H9 . . .'
    endif else got_h9 = 0
endif
if got_h10 eq 1 then begin
     if flux[h10] ge 20.*fluxerr[h10] then begin
        f_h10 = flux[h10]
        df_h10 = fluxerr[h10]
        lambda = [lambda,linewave[h10]]
        print,'Strong detection of H10 . . .'
    endif else got_h10 = 0
endif
if got_h11 eq 1 then begin
     if flux[h11] ge 20.*fluxerr[h11] then begin
        f_h11 = flux[h11]
        df_h11 = fluxerr[h11]
        lambda = [lambda,linewave[h11]]
        print,'Strong detection of H11 . . .'
    endif else got_h11 = 0
endif
if got_h12 eq 1 then begin
     if flux[h12] ge 20.*fluxerr[h12] then begin
        f_h12 = flux[h12]
        df_h12 = fluxerr[h12]
        lambda = [lambda,linewave[h12]]
        print,'Strong detection of H12 . . .'
    endif else got_h12 = 0
endif


;FIND OUT HOW MANY LINES WE CAN USE TO ESTIMATE THE EXTINCTION
nlines = got_hgamma + got_hdelta + got_h9 + got_h10 + got_h11 + got_h12


;INTRINSIC BALMER LINE FLUX RATIOS FOR Te = 10000 K, ne = 100 cm^-3 (HII)
;OR Te = 12500 K, ne = 1000 cm^-3 (PN)
;FROM HUMMER & STOREY 1987
;BETTER Te ESTIMATES ARE 9500 K FOR HII REGIONS AND 12000 K FOR PNE
;hii_region_fluxratios = [0.468,0.259,0.0731,0.0530,0.0397,0.0305]
;pne_fluxratios = [0.471,0.261,0.0737,0.0534,0.0400,0.0307]



;CALCULATE CONTINUUM LEVELS NEAR EACH LINE
junk = min(abs(w-4790),p4790)
junk = min(abs(w-4840),p4840)
junk = min(abs(w-4270),p4270)
junk = min(abs(w-4320),p4320)
junk = min(abs(w-4120),p4120)
junk = min(abs(w-4170),p4170)
junk = min(abs(w-3842),p3842)
junk = min(abs(w-3860),p3860)

continuum_hbeta = median(spec[p4840:p4790])
dcontinuum_hbeta = 1.25*sqrt(total(dspec[p4840:p4790]^2))/(p4790-p4840+1.)
continuum_hgamma = median(spec[p4320:p4270])
dcontinuum_hgamma = 1.25*sqrt(total(dspec[p4320:p4270]^2))/(p4270-p4320+1.)
continuum_hdelta = median(spec[p4170:p4120])
dcontinuum_hdelta = 1.25*sqrt(total(dspec[p4170:p4120]^2))/(p4120-p4170+1.)
if got_h9 then begin
    continuum_h9 = median(spec[p3860:p3842])
    dcontinuum_h9 = 1.25*sqrt(total(dspec[p3860:p3842]^2))/(p3842-p3860+1.)
    continuum_h10 = continuum_h9
    dcontinuum_h10 = dcontinuum_h9
    continuum_h11 = continuum_h9
    dcontinuum_h11 = dcontinuum_h9
    continuum_h12 = continuum_h9
    dcontinuum_h10 = dcontinuum_h9
endif



;SET INITIAL GUESSES OF BALMER ABSORPTION EW = 1 A AND C(HB) = 0.2
initial_guesses = dblarr(2)
initial_guesses[0] = 1.0d
initial_guesses[1] = 0.2d

;MAKE SURE NEITHER PARAMETER IS ALLOWED TO GO NEGATIVE
;AND PREVENT RIDICULOUSLY HIGH POSITIVE VALUES TOO
param_control = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D], $
                          step:0}, 2)
param_control[0].limited[0] = 1
param_control[0].limits[0] = 0
param_control[0].limited[1] = 1
param_control[0].limits[1] = 10
;param_control[0].step = 0
param_control[1].limited[0] = 1
param_control[1].limits[0] = 0
param_control[1].limited[1] = 1
param_control[1].limits[1] = 3


x=[f_hg,f_hb,continuum_hgamma,continuum_hbeta]
dx=[df_hg,df_hb,dcontinuum_hgamma,dcontinuum_hbeta]

;HAVE NOT IMPLEMENTED ABILITY TO USE H10 OR HIGHER LINES YET . . .
if nlines ge 3 then begin
    x = [x,f_hd,continuum_hdelta,f_h9,continuum_h9]
    dx = [dx,df_hd,dcontinuum_hdelta,df_h9,dcontinuum_h9]
endif else if nlines eq 2 then begin
    x = [x,f_hd,continuum_hdelta]
    dx = [dx,df_hd,dcontinuum_hdelta]
endif

variable_ew=1
functargs = {x: x, $
             dx: dx, $
             nlines: nlines, $
             variable_ew: keyword_set(variable_ew)}

try_again:
outparams = tnmin('extinctfunc', $
                  initial_guesses,functargs=functargs, /autoderiv, $
                  parinfo=param_control,bestmin=chisq1, $
                  quiet=keyword_set(silent), $
                  STATUS=status, ERRMSG=errmsg)

if status LE 0 then $
  if errmsg eq 'ERROR: Line search failed to converge' then begin
    print,errmsg
    goto,skip_error 
    endif else message, errmsg
skip_error:

;FIND UNCERTAINTIES ON PARAMETERS BY CHANGING THEM UNTIL DELTA CHI^2 = 1
for i = 0,1 do begin
    p = outparams 
    chisqtest=0
    while chisqtest lt (1 + chisq1) AND abs(p[i]) lt 1000 do begin 
        xr1 = (x[0]/x[1])*( (1 + p[0]*x[2]/x[0])/(1 + p[0]*x[3]/x[1]) ) * $
          10^(0.1565*p[1]) 
        if nlines ge 3 then begin

            if keyword_set(variable_ew) then begin
                hd_ew_multiplier = 1.2 
                h9_ew_multiplier = 0.9
            endif else begin
                hd_ew_multiplier = 1.0
                h9_ew_multiplier = 1.0
            endelse
            xr2 = (x[4]/x[1])*( (1 + hd_ew_multiplier*p[0]*x[5]/x[4])/(1 + p[0]*x[3]/x[1]) ) * $
              10^(0.2295*p[1])
            xr3 = (x[6]/x[1])*( (1 + h9_ew_multiplier*p[0]*x[7]/x[6])/(1 + p[0]*x[3]/x[1]) ) * $
              10^(0.2989*p[1])

        endif else begin

            if nlines eq 2 then begin
                if keyword_set(variable_ew) then hd_ew_multiplier = 1.2 else $
                  hd_ew_multiplier = 1.0
                xr2 = (x[4]/x[1])*( (1 + hd_ew_multiplier*p[0]*x[5]/x[4])/(1 + p[0]*x[3]/x[1]) ) * $
                  10^(0.2295*p[1]) 
            endif    

        endelse
        
        sig_xr1 = sqrt( ((10^(0.1565*p[1]))/(x[1] + p[0]*x[3]) * dx[0])^2 + $
                        ((x[0] + p[0]*x[2])*(10^(0.1565*p[1]))/(x[1] + p[0]*x[3])^2 * dx[1])^2 + $
                        (p[0]*(10^(0.1565*p[1]))/(x[1] + p[0]*x[3]) * dx[2])^2 + $
                        (p[0]*(x[0] + p[0]*x[2])*(10^(0.1565*p[1]))/(x[1] + p[0]*x[3])^2 * dx[3])^2 ) 
        if nlines ge 3 then begin
            sig_xr2 = sqrt( ((10^(0.2295*p[1]))/(x[1] + p[0]*x[3]) * dx[4])^2 + $
                            ((x[4] + hd_ew_multiplier*p[0]*x[5])*(10^(0.2295*p[1]))/(x[1] + p[0]*x[3])^2 * dx[1])^2 + $
                            (hd_ew_multiplier*p[0]*(10^(0.2295*p[1]))/(x[1] + p[0]*x[3]) * dx[5])^2 + $
                            (p[0]*(x[4] + hd_ew_multiplier*p[0]*x[5])*(10^(0.2295*p[1]))/(x[1] + p[0]*x[3])^2 * dx[3])^2 )
            sig_xr3 = sqrt( ((10^(0.2989*p[1]))/(x[1] + p[0]*x[3]) * dx[6])^2 + $
                            ((x[6] + h9_ew_multiplier*p[0]*x[7])*(10^(0.2989*p[1]))/(x[1] + p[0]*x[3])^2 * dx[1])^2 + $
                            (h9_ew_multiplier*p[0]*(10^(0.2989*p[1]))/(x[1] + p[0]*x[3]) * dx[7])^2 + $
                            (p[0]*(x[6] + h9_ew_multiplier*p[0]*x[7])*(10^(0.2989*p[1]))/(x[1] + p[0]*x[3])^2 * dx[3])^2 )
        endif else if nlines eq 2 then $
            sig_xr2 = sqrt( ((10^(0.2295*p[1]))/(x[1] + p[0]*x[3]) * dx[4])^2 + $
                            ((x[4] + hd_ew_multiplier*p[0]*x[5])*(10^(0.2295*p[1]))/(x[1] + p[0]*x[3])^2 * dx[1])^2 + $
                            (hd_ew_multiplier*p[0]*(10^(0.2295*p[1]))/(x[1] + p[0]*x[3]) * dx[5])^2 + $
                            (p[0]*(x[4] + hd_ew_multiplier*p[0]*x[5])*(10^(0.2295*p[1]))/(x[1] + p[0]*x[3])^2 * dx[3])^2 )

        chisqtest = (xr1 - 0.468)^2/sig_xr1^2 
        if nlines ge 2 then chisqtest = chisqtest + (xr2 - 0.259)^2/sig_xr2^2 
        if nlines ge 3 then  chisqtest = chisqtest + (xr3 - 0.0731)^2/sig_xr3^2  
        
        p[i] = p[i] + 0.001 
    endwhile

    diff1 = p[i] - outparams[i]

    p = outparams 
    chisqtest=0
    while chisqtest lt (1 + chisq1) AND abs(p[i]) lt 1000 do begin 
        xr1 = (x[0]/x[1])*( (1 + p[0]*x[2]/x[0])/(1 + p[0]*x[3]/x[1]) ) * $
          10^(0.1565*p[1]) 
        if nlines ge 3 then begin
            xr2 = (x[4]/x[1])*( (1 + hd_ew_multiplier*p[0]*x[5]/x[4])/(1 + p[0]*x[3]/x[1]) ) * $
              10^(0.2295*p[1])
            xr3 = (x[6]/x[1])*( (1 + h9_ew_multiplier*p[0]*x[7]/x[6])/(1 + p[0]*x[3]/x[1]) ) * $
              10^(0.2989*p[1])
        endif else if nlines eq 2 then $
            xr2 = (x[4]/x[1])*( (1 + hd_ew_multiplier*p[0]*x[5]/x[4])/(1 + p[0]*x[3]/x[1]) ) * $
              10^(0.2295*p[1]) 
        
        sig_xr1 = sqrt( ((10^(0.1565*p[1]))/(x[1] + p[0]*x[3]) * dx[0])^2 + $
                        ((x[0] + p[0]*x[2])*(10^(0.1565*p[1]))/(x[1] + p[0]*x[3])^2 * dx[1])^2 + $
                        (p[0]*(10^(0.1565*p[1]))/(x[1] + p[0]*x[3]) * dx[2])^2 + $
                        (p[0]*(x[0] + p[0]*x[2])*(10^(0.1565*p[1]))/(x[1] + p[0]*x[3])^2 * dx[3])^2 ) 
        if nlines ge 3 then begin
            sig_xr2 = sqrt( ((10^(0.2295*p[1]))/(x[1] + p[0]*x[3]) * dx[4])^2 + $
                            ((x[4] + hd_ew_multiplier*p[0]*x[5])*(10^(0.2295*p[1]))/(x[1] + p[0]*x[3])^2 * dx[1])^2 + $
                            (hd_ew_multiplier*p[0]*(10^(0.2295*p[1]))/(x[1] + p[0]*x[3]) * dx[5])^2 + $
                            (p[0]*(x[4] + hd_ew_multiplier*p[0]*x[5])*(10^(0.2295*p[1]))/(x[1] + p[0]*x[3])^2 * dx[3])^2 )
            sig_xr3 = sqrt( ((10^(0.2989*p[1]))/(x[1] + p[0]*x[3]) * dx[6])^2 + $
                            ((x[6] + h9_ew_multiplier*p[0]*x[7])*(10^(0.2989*p[1]))/(x[1] + p[0]*x[3])^2 * dx[1])^2 + $
                            (h9_ew_multiplier*p[0]*(10^(0.2989*p[1]))/(x[1] + p[0]*x[3]) * dx[7])^2 + $
                            (p[0]*(x[6] + h9_ew_multiplier*p[0]*x[7])*(10^(0.2989*p[1]))/(x[1] + p[0]*x[3])^2 * dx[3])^2 )
        endif else if nlines eq 2 then $
            sig_xr2 = sqrt( ((10^(0.2295*p[1]))/(x[1] + p[0]*x[3]) * dx[4])^2 + $
                            ((x[4] + hd_ew_multiplier*p[0]*x[5])*(10^(0.2295*p[1]))/(x[1] + p[0]*x[3])^2 * dx[1])^2 + $
                            (hd_ew_multiplier*p[0]*(10^(0.2295*p[1]))/(x[1] + p[0]*x[3]) * dx[5])^2 + $
                            (p[0]*(x[4] + hd_ew_multiplier*p[0]*x[5])*(10^(0.2295*p[1]))/(x[1] + p[0]*x[3])^2 * dx[3])^2 )

        chisqtest = (xr1 - 0.468)^2/sig_xr1^2 
        if nlines ge 2 then chisqtest = chisqtest + (xr2 - 0.259)^2/sig_xr2^2 
        if nlines ge 3 then  chisqtest = chisqtest + (xr3 - 0.0731)^2/sig_xr3^2  
        
        p[i] = p[i] - 0.001 
    endwhile

    diff2 = outparams[i] - p[i]

    if i eq 0 then d_eqw = (diff1 + diff2)/2. else d_chb = (diff1 + diff2)/2.

endfor


eqw = outparams[0]
chb = outparams[1]    

print,'Chi^2 minimization results:'
print,'Chi^2 = ',chisq1, format='(a8,f6.2)'
print,'C(Hbeta) = ',chb,' +/- ',d_chb, format='(a11,f6.3,a5,f5.3)'
print,'Equivalent width = ',eqw,' +/- ',d_eqw,' A', format='(a19,f6.3,a5,f5.3,a2)'
print,''

extinction_montecarlo,x,dx,nlines,mc_chb,mc_eqw,plot=keyword_set(plot)
d_chbmc = stdev(mc_chb)
d_eqwmc = stdev(mc_eqw)
print,'Monte Carlo results:'
print,'C(Hbeta) = ',mean(mc_chb),' +/- ',d_chbmc, format='(a11,f6.3,a5,f5.3)'
print,'Equivalent width = ',mean(mc_eqw),' +/- ',d_eqwmc,' A',$
  format='(a19,f6.3,a5,f5.3,a2)'

return
end


