function ahk_fluxmodelextract_multigrating,scifileblue,scifilered9,scifilered4,slitfileblue=slitfileblue,wavefileblue=wavefileblue,sensfuncfileblue=sensfuncfileblue,slitfilered9=slitfilered9,wavefilered9=wavefilered9,sensfuncfilered9=sensfuncfilered9,slitfilered4=slitfilered4,wavefilered4=wavefilered4,sensfuncfilered4=sensfuncfilered4

;+
; NAME:
;   ahk_fluxmodelextract_multigrating
;
; PURPOSE:
;  call long_extract_optimal after using objstruct.fluxmodel to
;  determine where objects are in the slits
;
;
; CALLING SEQUENCE:
;	result =  ahk_objextractandflux_multigrating('../../../lblue_600_4000/Science/sci-lblue2082.fits.gz','../../../lred_900_5500/Science/sci-lred2068.fits.gz','../../../lred_400_8500/Science/sci-lred2073.fits.gz',slitfileblue='../../../lblue_600_4000/slits-lblue2039.fits',wavefileblue='../../../lblue_600_4000/wave-lblue2027.fits',sensfuncfileblue='feige34_bsens_lblue2166.fits',slitfilered9='../../../lred_900_5500/slits-lred2031.fits',wavefilered9='../../../lred_900_5500/wave-lred2158.fits',sensfuncfilered9='feige34_bsens_lred2119.fits',slitfilered4='../../../lred_400_8500/slits-lred2122.fits',wavefilered4='../../../lred_400_8500/wave-lred2019.fits',sensfuncfilered4='feige34_bsens_lred2117.fits')

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
;	ahk_profile_singlegrating
;	ahk_profile_singleslit
;	ahk_comboprofile
;	ahk_profile2imarray
;	long_extract_optimal
;
; REVISION HISTORY:
;   3-November-2014 -- Written by Alice Koning
;-

;; Read in one of the slitfiles to be used later to determine which side of the ccd each slit is on
slitim=xmrdfits(slitfileblue,0)

;; Number of slits in image
slitmaskblue = mrdfits(slitfileblue)
nslitblue = max(slitmaskblue,/NaN)
slitmaskred9 = mrdfits(slitfilered9)
nslitred9 = max(slitmaskred9,/NaN)
slitmaskred4 = mrdfits(slitfilered4)
nslitred4 = max(slitmaskred4,/NaN)

IF (nslitblue NE nslitred9) OR (nslitblue NE nslitred4) OR (nslitred9 NE nslitred4) THEN BEGIN
	print, 'Number of slits does not match in all 3 slitfiles'
	return, 0
     ENDIF

;; This will be used to find the master profile in each slit
;; Define comboyfitnorm as an array of pointers so each row can have a different length (where each row gives comboyfitnorm for a given slit)
comboyfitnorm = replicate(Ptr_New(),nslitblue)
comboyfitnorm = PtrArr(nslitblue, /Allocate_heap)

scifileNames = [scifileblue,scifilered9,scifilered4] ;;maybe this isn't the best plan of attack afterall...

FOREACH scifile, scifileNames DO BEGIN

   print, 'Science filename is ', scifile

   ;; Restore saved objstr.FLUXMODEL
   RESTORE, repstr(scifile,'fits.gz','sav')
   stop ;Check everything is in order here.

   ;; Get various required extensions out of science file
   ivar = xmrdfits(scifile,1)
   skyimage = xmrdfits(scifile,2)
   objimage = xmrdfits(scifile,3)
   outmask = xmrdfits(scifile,4)
   waveimg = xmrdfits(wavefile, silent = (keyword_set(verbose) EQ 0), 0)
   sciimg = xmrdfits(scifile)
   imgminsky = sciimg - skyimage
 
   FOR ii =1,nslitblue DO BEGIN
	;;Find first object in final_struct with slitid equal to ii (avoids duplicates for slits with more than one object)
	structid = min(WHERE(final_struct.slitid EQ (ii)))
	
	;;Start next iteration if ii is not a slitid in final_struct
	IF (structid EQ -1) THEN BEGIN
           print, '***** SLITID ', slitid, ' is not in final_struct *****'
           CONTINUE
        ENDIF

	print, '***** SLITID: ', slitid, ' *****'
	slitid = final_struct[structid].SLITID
        fluxmodel = final_struct[structid].fluxmodel

	ccdloc = MAX(WHERE(slitim EQ ii))
	ccdind = ARRAY_INDICES(slitim, ccdloc)
	IF ccdind[0] LE 2048 THEN ccdside=0 ELSE ccdside=1

	undefine, csvblue, slitindexblue, yfitnormblue
	search = file_search('yfitnorm_'+strmid(scifileblue,15,8,/reverse_offset)+'_slitid'+StrTrim(ii,2)+'.csv',count=c1)
	IF (c1 GT 0) THEN BEGIN
		csvblue = read_csv('yfitnorm_'+strmid(scifileblue,15,8,/reverse_offset)+'_slitid'+StrTrim(ii,2)+'.csv')
		slitindexblue = csvblue.field1
		yfitnormblue = csvblue.field2
	ENDIF

	undefine, csvred9, slitindexred9, yfitnormred9
	search = file_search('yfitnorm_'+strmid(scifilered9,15,8,/reverse_offset)+'_slitid'+StrTrim(ii,2)+'.csv',count=c2)
	IF (c2 GT 0) THEN BEGIN
		csvred9 = read_csv('yfitnorm_'+strmid(scifilered9,15,8,/reverse_offset)+'_slitid'+StrTrim(ii,2)+'.csv')
		slitindexred9 = csvred9.field1
		yfitnormred9 = csvred9.field2
	ENDIF

	undefine, csvred4, slitindexred4, yfitnormred4
	search = file_search('yfitnorm_'+strmid(scifilered4,15,8,/reverse_offset)+'_slitid'+StrTrim(ii,2)+'.csv',count=c3)
	IF (c3 GT 0) THEN BEGIN
		csvred4 = read_csv('yfitnorm_'+strmid(scifilered4,15,8,/reverse_offset)+'_slitid'+StrTrim(ii,2)+'.csv')
		slitindexred4 = csvred4.field1
		yfitnormred4 = csvred4.field2
	ENDIF

	;*comboyfitnorm[ii-1]=ahk_comboprofile_nopad(ii,slitindexblue,yfitnormblue,slitindexred9,yfitnormred9,slitindexred4,yfitnormred4)
	*comboyfitnorm[ii-1]=ahk_comboprofile(ii,slitindexblue,yfitnormblue,slitindexred9,yfitnormred9,slitindexred4,yfitnormred4,ccdside=ccdside)

   ENDFOR ;; end loop over slits to find comboyfitnorm in each


;; Loop over slits in each grating individually now. Convert comboyfitnorm to an array, then feed to long_extract_optimal.
;; Blue grating:
FOR ii =1,nslitblue DO BEGIN
	search = file_search('yfitnorm_'+strmid(scifileblue,15,8,/reverse_offset)+'_slitid'+StrTrim(ii,2)+'.csv',count=c1)
	;; If CSV file exists, and therefore there is flux in that slit, then continue with the analysis. Otherwise skip to next slit.
	IF (c1 LT 1) THEN CONTINUE
		print, '***** Working on BLUE SLITID: ', ii, ' *****'

		;; Get various required extensions out of science file
		ivar = xmrdfits(scifileblue,1)
		skyimage = xmrdfits(scifileblue,2)
		objimage = xmrdfits(scifileblue,3)
		outmask = xmrdfits(scifileblue,4)
		waveimg = xmrdfits(wavefileblue, silent = (keyword_set(verbose) EQ 0), 0)
		sciimg = xmrdfits(scifileblue)
		imgminsky = sciimg - skyimage

		;; Size of science image
		nx = (size(sciimg))[1]
		ny = (size(sciimg))[2]
		xarr = findgen(nx)## replicate(1.0, ny)
		yarr = findgen(ny)## replicate(1.0, nx)

		;; Must determine the normalized slitindex. Use our knowledge that comboyfitnorm abscissa values run from 0 to 1.
		slitindex = findgen(N_ELEMENTS(*comboyfitnorm[ii-1]))
		normslitindex = slitindex/MAX(slitindex)

		;; Get object profile along slit (i.e. convert this slit's comboyfitnorm to a 2-d array with non-zero values running along slit length)
		slitprofile = ahk_profile2imarray(normslitindex,*comboyfitnorm[ii-1],slitid=ii,slitfile=slitfileblue)

		;; Number of objects in slit. If more than one, must get 1D spectrum for each.
		ymask = (slitprofile GT 0.01) ;; Assumes yfit is normalized.
		label = label_region(ymask)
		nobj = max(label)
		print, 'Number of objects in slit: ', nobj

		IF nobj EQ 0 THEN BEGIN
			print, 'No objects in profile GT 0.01'
		ENDIF ELSE BEGIN	
			FOR objid=1,nobj DO BEGIN
				trace=MAKE_ARRAY(ny)
				objlength=MAKE_ARRAY(ny)
				;; Need to loop over rows to find trace (xcen) and average number of included pixels on either side to be box_rad
				FOR row=1,ny-2 DO BEGIN ;; First and last rows of label are all zeros, so skip here
					objindx = []
					objindx = WHERE(label[*,row] EQ objid)
					;print, 'ROW', row, 'OBJINDX', objindx
					trace[row] = MEDIAN(objindx)
					objlength[row] = N_ELEMENTS(objindx)/2.			
				ENDFOR

				box_rad = MEAN(objlength)

				objmask = MAKE_ARRAY(nx,ny)
				FOR row=1,ny-2 DO BEGIN
					objmask[*,row] = ((xarr[row,*] GT (trace[row] - 1.5*box_rad)) AND (xarr[row,*] LE (trace[row] + 1.5*box_rad)))
				ENDFOR

				spec = long_extract_optimal(waveimg, imgminsky, ivar, slitprofile, outmask*objmask, skyimage, trace, BOX_RAD=box_rad)
	
				;; Output all images etc as they were, but now last science extension contains spec
				rdfits_struct,scifileblue,scistr ;; read original scifile into structure
				imgname = STRMID(scifileblue, 15, 8, /REVERSE_OFFSET)
				fileout = imgname+'_slit'+STRTRIM(ii,2)+'_obj'+STRTRIM(objid,2)+'.fits'
				mwrfits,scistr.im0,fileout,scistr.hdr0,/create
				mwrfits,scistr.im1,fileout,scistr.hdr1 
				mwrfits,scistr.im2,fileout,scistr.hdr2 
				mwrfits,scistr.im3,fileout,scistr.hdr3 
				mwrfits,scistr.im4,fileout,scistr.hdr4 
				mwrfits,spec,fileout,hdr_struct
			
				;; Flux calibration
				fluxedfileout = imgname+'_slit'+STRTRIM(ii,2)+'_obj'+STRTRIM(objid,2)+'_fluxed.fits'
				long_fluxcal, fileout, SENSFUNCFILE=sensfuncfileblue, OUTFIL=fluxedfileout, frm_sci=1
			
			ENDFOR ;; end loop over objects

		ENDELSE

ENDFOR ;; end loop over blue slits

;; red9 grating:
FOR ii =1,nslitred9 DO BEGIN
	search = file_search('yfitnorm_'+strmid(scifilered9,15,8,/reverse_offset)+'_slitid'+StrTrim(ii,2)+'.csv',count=c1)
	;; If CSV file exists, and therefore there is flux in that slit, then continue with the analysis. Otherwise skip to next slit.
	IF (c1 LT 1) THEN CONTINUE
		print, '***** Working on RED9 SLITID: ', ii, ' *****'

		;; Get various required extensions out of science file
		ivar = xmrdfits(scifilered9,1)
		skyimage = xmrdfits(scifilered9,2)
		objimage = xmrdfits(scifilered9,3)
		outmask = xmrdfits(scifilered9,4)
		waveimg = xmrdfits(wavefilered9, silent = (keyword_set(verbose) EQ 0), 0)
		sciimg = xmrdfits(scifilered9)
		imgminsky = sciimg - skyimage

		;; Size of science image
		nx = (size(sciimg))[1]
		ny = (size(sciimg))[2]
		xarr = findgen(nx)## replicate(1.0, ny)
		yarr = findgen(ny)## replicate(1.0, nx)

		;; Must determine the normalized slitindex. Use our knowledge that comboyfitnorm abscissa values run from 0 to 1.
		slitindex = findgen(N_ELEMENTS(*comboyfitnorm[ii-1]))
		normslitindex = slitindex/MAX(slitindex)

		;; Get object profile along slit (i.e. convert this slit's comboyfitnorm to a 2-d array with non-zero values running along slit length)
		slitprofile = ahk_profile2imarray(normslitindex,*comboyfitnorm[ii-1],slitid=ii,slitfile=slitfilered9)

		;; Number of objects in slit. If more than one, must get 1D spectrum for each.
		ymask = (slitprofile GT 0.01) ;; Assumes yfit is normalized.
		label = label_region(ymask)
		nobj = max(label)
		print, 'Number of objects in slit: ', nobj

		IF nobj EQ 0 THEN BEGIN
			print, 'No objects in profile GT 0.01'
		ENDIF ELSE BEGIN	
			FOR objid=1,nobj DO BEGIN
				trace=MAKE_ARRAY(ny)
				objlength=MAKE_ARRAY(ny)
				;; Need to loop over rows to find trace (xcen) and average number of included pixels on either side to be box_rad
				FOR row=1,ny-2 DO BEGIN ;; First and last rows of label are all zeros, so skip here
					objindx = []
					objindx = WHERE(label[*,row] EQ objid)
					;print, 'ROW', row, 'OBJINDX', objindx
					trace[row] = MEDIAN(objindx)
					objlength[row] = N_ELEMENTS(objindx)/2.			
				ENDFOR

				box_rad = MEAN(objlength)

				objmask = MAKE_ARRAY(nx,ny)
				FOR row=1,ny-2 DO BEGIN
					objmask[*,row] = ((xarr[row,*] GT (trace[row] - 1.5*box_rad)) AND (xarr[row,*] LE (trace[row] + 1.5*box_rad)))
				ENDFOR

				spec = long_extract_optimal(waveimg, imgminsky, ivar, slitprofile, outmask*objmask, skyimage, trace, BOX_RAD=box_rad)
	
				;; Output all images etc as they were, but now last science extension contains spec
				rdfits_struct,scifilered9,scistr ;; read original scifile into structure
				imgname = STRMID(scifilered9, 15, 8, /REVERSE_OFFSET)
				fileout = imgname+'_slit'+STRTRIM(ii,2)+'_obj'+STRTRIM(objid,2)+'.fits'
				mwrfits,scistr.im0,fileout,scistr.hdr0,/create
				mwrfits,scistr.im1,fileout,scistr.hdr1 
				mwrfits,scistr.im2,fileout,scistr.hdr2 
				mwrfits,scistr.im3,fileout,scistr.hdr3 
				mwrfits,scistr.im4,fileout,scistr.hdr4 
				mwrfits,spec,fileout,hdr_struct
			
				;; Flux calibration
				fluxedfileout = imgname+'_slit'+STRTRIM(ii,2)+'_obj'+STRTRIM(objid,2)+'_fluxed.fits'
				long_fluxcal, fileout, SENSFUNCFILE=sensfuncfilered9, OUTFIL=fluxedfileout, frm_sci=1
			
			ENDFOR ;; end loop over objects

		ENDELSE

ENDFOR ;; end loop over red9 slits


;; red4 grating:
FOR ii =1,nslitred4 DO BEGIN
	search = file_search('yfitnorm_'+strmid(scifilered4,15,8,/reverse_offset)+'_slitid'+StrTrim(ii,2)+'.csv',count=c1)
	;; If CSV file exists, and therefore there is flux in that slit, then continue with the analysis. Otherwise skip to next slit.
	IF (c1 LT 1) THEN CONTINUE
		print, '***** Working on RED4 SLITID: ', ii, ' *****'

		;; Get various required extensions out of science file
		ivar = xmrdfits(scifilered4,1)
		skyimage = xmrdfits(scifilered4,2)
		objimage = xmrdfits(scifilered4,3)
		outmask = xmrdfits(scifilered4,4)
		waveimg = xmrdfits(wavefilered4, silent = (keyword_set(verbose) EQ 0), 0)
		sciimg = xmrdfits(scifilered4)
		imgminsky = sciimg - skyimage

		;; Size of science image
		nx = (size(sciimg))[1]
		ny = (size(sciimg))[2]
		xarr = findgen(nx)## replicate(1.0, ny)
		yarr = findgen(ny)## replicate(1.0, nx)

		;; Must determine the normalized slitindex. Use our knowledge that comboyfitnorm abscissa values run from 0 to 1.
		slitindex = findgen(N_ELEMENTS(*comboyfitnorm[ii-1]))
		normslitindex = slitindex/MAX(slitindex)

		;; Get object profile along slit (i.e. convert this slit's comboyfitnorm to a 2-d array with non-zero values running along slit length)
		slitprofile = ahk_profile2imarray(normslitindex,*comboyfitnorm[ii-1],slitid=ii,slitfile=slitfilered4)

		;; Number of objects in slit. If more than one, must get 1D spectrum for each.
		ymask = (slitprofile GT 0.01) ;; Assumes yfit is normalized.
		label = label_region(ymask)
		nobj = max(label)
		print, 'Number of objects in slit: ', nobj

		IF nobj EQ 0 THEN BEGIN
			print, 'No objects in profile GT 0.01'
		ENDIF ELSE BEGIN	
			FOR objid=1,nobj DO BEGIN
				trace=MAKE_ARRAY(ny)
				objlength=MAKE_ARRAY(ny)
				;; Need to loop over rows to find trace (xcen) and average number of included pixels on either side to be box_rad
				FOR row=1,ny-2 DO BEGIN ;; First and last rows of label are all zeros, so skip here
					objindx = []
					objindx = WHERE(label[*,row] EQ objid)
					;print, 'ROW', row, 'OBJINDX', objindx
					trace[row] = MEDIAN(objindx)
					objlength[row] = N_ELEMENTS(objindx)/2.			
				ENDFOR

				box_rad = MEAN(objlength)

				objmask = MAKE_ARRAY(nx,ny)
				FOR row=1,ny-2 DO BEGIN
					objmask[*,row] = ((xarr[row,*] GT (trace[row] - 1.5*box_rad)) AND (xarr[row,*] LE (trace[row] + 1.5*box_rad)))
				ENDFOR

				spec = long_extract_optimal(waveimg, imgminsky, ivar, slitprofile, outmask*objmask, skyimage, trace, BOX_RAD=box_rad)
	
				;; Output all images etc as they were, but now last science extension contains spec
				rdfits_struct,scifilered4,scistr ;; read original scifile into structure
				imgname = STRMID(scifilered4, 15, 8, /REVERSE_OFFSET)
				fileout = imgname+'_slit'+STRTRIM(ii,2)+'_obj'+STRTRIM(objid,2)+'.fits'
				mwrfits,scistr.im0,fileout,scistr.hdr0,/create
				mwrfits,scistr.im1,fileout,scistr.hdr1 
				mwrfits,scistr.im2,fileout,scistr.hdr2 
				mwrfits,scistr.im3,fileout,scistr.hdr3 
				mwrfits,scistr.im4,fileout,scistr.hdr4 
				mwrfits,spec,fileout,hdr_struct

			ENDFOR ;; end loop over objects

		ENDELSE

ENDFOR ;; end loop over red4 slits

return, spec

end


