;+
; NAME:
;   extinction_montecarlo
;
; PURPOSE:
;   Run a Monte Carlo simulation to estimate uncertainties on equivalent
;   width and extinction values derived by jds_hiiextinct.pro (following 
;   Olive & Skillman 2001)
; 
; CALLING SEQUENCE:
;   extinction_montecarlo,x,dx,nlines,mc_chb,mc_eqw
;
; INPUTS:
;   xorig           - Array of line fluxes to feed into extinction fitter
;   dxorig          - Array of line flux errors to feed into extinction fitter
;   nlines          - Number of Balmer lines being used in the fitting
;                     (not counting Hbeta)
;
; OPTIONAL INPUTS:
;   n               - Number of iterations to run the Monte Carlo (default 1000)
;
; KEYWORDS:
;   verbose         - Print the output of tnmin
;   plot            - Plot the Monte Carlo results
;
; OUTPUTS:
;   mc_chb          - Monte Carlo results for c(Hb)
;   mc_eqw          - Monte Carlo results for Balmer EW
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
;   tnmin (Markwardt)
;   extinctfunc (JDS)   
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   26-Aug-2006  J. Simon, Caltech
;
;------------------------------------------------------------------------------
pro extinction_montecarlo,xorig,dxorig,nlines,mc_chb,mc_eqw, $
                          N=N, PLOT=PLOT, VERBOSE=verbose


;SET KEYWORDS
if not keyword_set(n) then n=1000
if keyword_set(verbose) then silent = 0 else silent = 1

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
;param_control[0].step = 0.0
param_control[1].limited[0] = 1
param_control[1].limits[0] = 0
param_control[1].limited[1] = 1
param_control[1].limits[1] = 3

mc_chb = fltarr(n)
mc_eqw = fltarr(n)
dx = dxorig
mcseed = long(systime(1) - 1.155d+09)
;RUN MONTE CARLO
for i = 0,n-1 do begin
  x = 0*xorig  
  x[0] = xorig[0] + randomn(mcseed)*dx[0]
  x[1] = xorig[1] + randomn(mcseed)*dx[1]
  x[2] = xorig[2] + randomn(mcseed)*dx[2]
  x[3] = xorig[3] + randomn(mcseed)*dx[3]

  if nlines ge 2 then begin
      x[4] = xorig[4] + randomn(mcseed)*dx[4]
      x[5] = xorig[5] + randomn(mcseed)*dx[5]
  endif

  if nlines ge 3 then begin
      x[6] = xorig[6] + randomn(mcseed)*dx[6]
      x[7] = xorig[7] + randomn(mcseed)*dx[7]
  endif

  functargs = {x: x, $
               dx: dx, $
               nlines: nlines}

  outparams = tnmin('extinctfunc', $
                    initial_guesses,functargs=functargs, /autoderiv, $
                    parinfo=param_control,bestmin=chisq1, $
                    quiet=keyword_set(silent), $
                    STATUS=status, ERRMSG=errmsg)

  mc_eqw[i] = outparams[0]
  mc_chb[i] = outparams[1]
endfor

if keyword_set(plot) then begin
    wset,0
    plot,mc_chb,mc_eqw,ps=3
    window,1
    oldpmulti= !p.multi
    !p.multi=[0,2,1]
    junk = jdshist(mc_eqw,min=0,max=1.1*max(mc_eqw), $
                   bin=0.05*max(mc_eqw) > 0.01,xax,/plot)
    junk = jdshist(mc_chb,min=0,max=1.1*max(mc_chb), $
                   bin=0.05*max(mc_chb) > 0.01,xax,/plot)
    !p.multi = oldpmulti
endif

end
