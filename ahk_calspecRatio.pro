pro ahk_calspecRatio,obsSpecFile,calSpecFile

;+
; NAME:
;   ahk_calspecRatio
;
; PURPOSE:
;  Given the actual observed spectrum of a standard star and the
;  expected spectrum, plot their ratio as a function of wavelength.
;
; CALLING SEQUENCE:
;  .compile ahk_calspecRatio
;  ahk_calspecRatio,obsSpecFile,calSpecFile
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
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   28-Jan-2015 -- Written by Alice Koning
;-

;;Read in observed spectrum wavelength and flux
obsWave=mrdfits(obsSpecFile,2,hd)
obsFlux1=mrdfits(obsSpecFile,0,hd)
obsFlux=1d-17*obsFlux1

;;Read in calibration spectrum wavelength and flux
calSpec=mrdfits(calSpecFile,1,hd)
calWave=calSpec.wavelength
calFlux=calSpec.flux

;;Interpolate calSpec Flux values at obsSpec wavelengths
calFlux_interp = INTERPOL(calFlux,calWave,obsWave)

;;Calculate ratio
ratio=obsFlux/calFlux_interp

;pobs = PLOT(obsWave, obsFlux, 'b', thick=5, title='Observed Spec (black) and Calibration Spec (red)')
;pcal = PLOT(calWave, calFlux, 'r', thick=5, /OVERPLOT)

p = PLOT(obsWave, ratio, 'b', thick=2, title='Ratio: Observed/Calibration Spectrum')


end
