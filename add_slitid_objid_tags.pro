function add_slitid_objid_tags, struct, slitid=slit, objid=obj

; This code will return a new struct that contains all the same info
; as the input structure, plus the added tags of SLIT_ID and OBJ_ID
; with values specified by the input line


nt = size(struct.wave_opt, /dimen)
box_rad = struct.box_rad

new_struct = create_struct('WAVE_OPT', dblarr(nt) $      ; optimal wavelengths
                       , 'FLUX_OPT', fltarr(nt) $ ; optimal flux
                       , 'SIVAR_OPT', fltarr(nt) $ ; optimal inverse var
                       , 'IVAR_OPT', fltarr(nt) $ ; model optimal inverse var
                       , 'SKY_OPT', fltarr(nt) $ ; optimally exttracted sky
                       , 'RN_OPT', fltarr(nt) $ ; sigma from RN in combined 
                       , 'NIVAR_OPT', fltarr(nt) $ ; optimal sky + RN ivar
                       , 'MASK_OPT', bytarr(nt) $ ; optimal mask
                       , 'FRAC_USE', fltarr(nt) $ ; frac of profile pixels used
                       , 'CHI2', fltarr(nt) $      ; chi2 of optimal model fit
                       , 'WAVE_BOX', dblarr(nt) $  ; boxcar wavelengths
                       , 'FLUX_BOX', fltarr(nt) $  ; boxcar flux
                       , 'SIVAR_BOX', fltarr(nt) $  ; boxcar inverse var
                       , 'IVAR_BOX', fltarr(nt) $ ; boxcar model inverse var
                       , 'NIVAR_BOX', fltarr(nt) $ ; boxcar sky + RN ivar
                       , 'SKY_BOX', fltarr(nt) $   ; boxcar sky
                       , 'RN_BOX', fltarr(nt) $    ; sigma from RN in combined
                       , 'MASK_BOX', bytarr(nt) $    ; optimal mask
                       , 'MINCOL', 0L, 'MAXCOL', 0L $ ; minmax of profile fits
                       , 'BOX_RAD', box_rad $          ; boxcar radius
                       , 'SLIT_ID', 0 $          ; Slit ID
                       , 'OBJ_ID', 0 )          ; Object ID

struct_assign, struct, new_struct
new_struct.SLIT_ID=slit
new_struct.OBJ_ID=obj

return, new_struct
end
