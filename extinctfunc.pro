function extinctfunc, p, x=x, dx=dx, nlines=nlines, variable_ew = variable_ew
  
;p[0] - EW
;p[1] - c(HBeta)

;x[0] = Hgamma flux
;x[1] = Hbeta flux
;x[2] = Hgamma continuum level
;x[3] = Hbeta continuum level
;x[4] = Hdelta flux
;x[5] = Hdelta continuum level
;x[6] = H9 flux
;x[7] = H9 continuum level
if nlines eq 1 then begin

  xr1 = (x[0]/x[1])*( (1 + p[0]*x[2]/x[0])/(1 + p[0]*x[3]/x[1]) ) * $
    10^(0.1565*p[1])

  sig_xr1 = sqrt( ((10^(0.1565*p[1]))/(x[1] + p[0]*x[3]) * dx[0])^2 + $
    ((x[0] + p[0]*x[2])*(10^(0.1565*p[1]))/(x[1] + p[0]*x[3])^2 * dx[1])^2 + $
    (p[0]*(10^(0.1565*p[1]))/(x[1] + p[0]*x[3]) * dx[2])^2 + $
    (p[0]*(x[0] + p[0]*x[2])*(10^(0.1565*p[1]))/(x[1] + p[0]*x[3])^2 * dx[3])^2 )

  return, (xr1 - 0.468)^2/sig_xr1^2
endif

if nlines eq 2 then begin

if keyword_set(variable_ew) then hd_ew_multiplier = 1.2 else $
  hd_ew_multiplier = 1.0

  xr1 = (x[0]/x[1])*( (1 + p[0]*x[2]/x[0])/(1 + p[0]*x[3]/x[1]) ) * $
    10^(0.1565*p[1])
  xr2 = (x[4]/x[1])*( (1 + hd_ew_multiplier*p[0]*x[5]/x[4])/(1 + p[0]*x[3]/x[1]) ) * $
    10^(0.2295*p[1])

  sig_xr1 = sqrt( ((10^(0.1565*p[1]))/(x[1] + p[0]*x[3]) * dx[0])^2 + $
    ((x[0] + p[0]*x[2])*(10^(0.1565*p[1]))/(x[1] + p[0]*x[3])^2 * dx[1])^2 + $
    (p[0]*(10^(0.1565*p[1]))/(x[1] + p[0]*x[3]) * dx[2])^2 + $
    (p[0]*(x[0] + p[0]*x[2])*(10^(0.1565*p[1]))/(x[1] + p[0]*x[3])^2 * dx[3])^2 )


  sig_xr2 = sqrt( ((10^(0.2295*p[1]))/(x[1] + p[0]*x[3]) * dx[4])^2 + $
    ((x[4] + hd_ew_multiplier*p[0]*x[5])*(10^(0.2295*p[1]))/(x[1] + p[0]*x[3])^2 * dx[1])^2 + $
    (hd_ew_multiplier*p[0]*(10^(0.2295*p[1]))/(x[1] + p[0]*x[3]) * dx[5])^2 + $
    (p[0]*(x[4] + hd_ew_multiplier*p[0]*x[5])*(10^(0.2295*p[1]))/(x[1] + p[0]*x[3])^2 * dx[3])^2 )


  return, (xr1 - 0.468)^2/sig_xr1^2 + (xr2 - 0.259)^2/sig_xr2^2 
endif

if nlines eq 3 then begin

if keyword_set(variable_ew) then begin
  hd_ew_multiplier = 1.2 
  h9_ew_multiplier = 0.9
endif else begin
  hd_ew_multiplier = 1.0
  h9_ew_multiplier = 1.0
endelse

  xr1 = (x[0]/x[1])*( (1 + p[0]*x[2]/x[0])/(1 + p[0]*x[3]/x[1]) ) * $
    10^(0.1565*p[1])
  xr2 = (x[4]/x[1])*( (1 + hd_ew_multiplier*p[0]*x[5]/x[4])/(1 + p[0]*x[3]/x[1]) ) * $
    10^(0.2295*p[1])
  xr3 = (x[6]/x[1])*( (1 + h9_ew_multiplier*p[0]*x[7]/x[6])/(1 + p[0]*x[3]/x[1]) ) * $
    10^(0.2989*p[1])

  sig_xr1 = sqrt( ((10^(0.1565*p[1]))/(x[1] + p[0]*x[3]) * dx[0])^2 + $
    ((x[0] + p[0]*x[2])*(10^(0.1565*p[1]))/(x[1] + p[0]*x[3])^2 * dx[1])^2 + $
    (p[0]*(10^(0.1565*p[1]))/(x[1] + p[0]*x[3]) * dx[2])^2 + $
    (p[0]*(x[0] + p[0]*x[2])*(10^(0.1565*p[1]))/(x[1] + p[0]*x[3])^2 * dx[3])^2 )
  sig_xr2 = sqrt( ((10^(0.2295*p[1]))/(x[1] + p[0]*x[3]) * dx[4])^2 + $
    ((x[4] + hd_ew_multiplier*p[0]*x[5])*(10^(0.2295*p[1]))/(x[1] + p[0]*x[3])^2 * dx[1])^2 + $
    (hd_ew_multiplier*p[0]*(10^(0.2295*p[1]))/(x[1] + p[0]*x[3]) * dx[5])^2 + $
    (p[0]*(x[4] + hd_ew_multiplier*p[0]*x[5])*(10^(0.2295*p[1]))/(x[1] + p[0]*x[3])^2 * dx[3])^2 )
  sig_xr3 = sqrt( ((10^(0.2989*p[1]))/(x[1] + p[0]*x[3]) * dx[6])^2 + $
    ((x[6] + h9_ew_multiplier*p[0]*x[7])*(10^(0.2989*p[1]))/(x[1] + p[0]*x[3])^2 * dx[1])^2 + $
    (h9_ew_multiplier*p[0]*(10^(0.2989*p[1]))/(x[1] + p[0]*x[3]) * dx[7])^2 + $
    (p[0]*(x[6] + h9_ew_multiplier*p[0]*x[7])*(10^(0.2989*p[1]))/(x[1] + p[0]*x[3])^2 * dx[3])^2 )


  return, (xr1 - 0.468)^2/sig_xr1^2 + (xr2 - 0.259)^2/sig_xr2^2 + $
    (xr3 - 0.0731)^2/sig_xr3^2 
endif


end
