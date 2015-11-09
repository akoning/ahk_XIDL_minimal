function ahk_objextractandflux_multigrating,scifileblue,scifilered9,scifilered4,slitfileblue=slitfileblue,wavefileblue=wavefileblue,sensfuncfileblue=sensfuncfileblue,slitfilered9=slitfilered9,wavefilered9=wavefilered9,sensfuncfilered9=sensfuncfilered9,slitfilered4=slitfilered4,wavefilered4=wavefilered4,sensfuncfilered4=sensfuncfilered4

;+
; NAME:
;   ahk_objextractandflux_multigrating
;
; PURPOSE:
;   **Same as ahk_objextract_multigrating, with the additional input sensitivity function for each grating, and additional final step to call long_fluxcal for flux calibration.
;   **Param to watch out for: frm_sci=1 for every call to long_fluxcal. Need to double check if this works correctly/as expected.
;
;   Uses output of long_reduce. First restore fluxmodel of each slit profile, in each grating. Combine these profiles to get a master profile.
;   Create a 2D array of this profile running along the slit, then extract objects accordingly with long_extract_optimal.
;   General outline
;	1) Determine number of slits in each grating from
;	slitfileblue, slitfilered9, slitfilered4. If not equal and red
;	has one slit less than blue, then ask for user input of whether the red side is missing a slit
;	at the beginning or end and  make the appropriate
;	adjustements. Otherwise, exit.
;	-- Steps which are also repeated for each grating:
;	2) Restore saved objstr.FLUXMODEL (1-d array of summed fluxes in our lines of interest) -- done by ahk_profile_singlegrating
;	3) save csv of  normalized FLUXMODEL -- done by ahk_profile_singleslit.pro within ahk_profile_singlegrating.
;	-- Steps in loop over slits:
;	4) Combine all DEFINED profiles (i.e. csv file exists) to get master profile for each slit using ahk_comboprofile.pro
;	-- Steps in new loop over slits, repeated for each grating:
;	5) Create 4096x4096 profile -- done by ahk_profile2imarray.pro
;	6) Feed profile into long_extract_optimal
;	7) Replace last fits extension in science image with output from long_extract_optimal
;	8) Call to long_fluxcal
;
; CALLING SEQUENCE:
;	result =  ahk_objextractandflux_multigrating('../Science_blue600/sci-lblue2082.fits.gz','../Science_red900/sci-lred2068.fits.gz','../Science_red400/sci-lred2073.fits.gz',slitfileblue='../slits-lblue2039.fits',wavefileblue='../wave-lblue2027.fits',sensfuncfileblue='../../../fluxCalibration/blue600_bsens_all.fits',slitfilered9='../slits-lred2031.fits',wavefilered9='../wave-lred2158.fits',sensfuncfilered9='../../../fluxCalibration/red900_bsens_all.fits',slitfilered4='../slits-lred2122.fits',wavefilered4='../wave-lred2019.fits',sensfuncfilered4='../../../fluxCalibration/red400_bsens_lred2064.fits')
;	result =  ahk_objextractandflux_multigrating('../Science_blue600/sci-lblue2092.fits.gz','../Science_red900/sci-lred2074.fits.gz','../Science_red400/sci-lred2079.fits.gz',slitfileblue='../slits-lblue2044.fits',wavefileblue='../wave-lblue2029.fits',sensfuncfileblue='feige34_bsens_lblue2166.fits',slitfilered9='../slits-lred2035.fits',wavefilered9='../wave-lred2012.fits',sensfuncfilered9='feige34_bsens_lred2119.fits',slitfilered4='../slits-lred2078.fits',wavefilered4='../wave-lred2169.fits',sensfuncfilered4='feige34_bsens_lred2117.fits')
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

;;Search for case when red gratings have one fewer slits than
;;blue. Ask user to specify which slit is missing in the red and
;;adjust the slitmasks accordingly **THIS HAS NOT BEEN FULLY TESTED
missSlit=0   
IF (nslitred9 EQ nslitred4) AND (nslitblue EQ nslitred9+1) THEN BEGIN
        read, 'Enter 1 if first slit in red missing, enter 2 if last slit in red missing: ',missSlit
        CASE missSlit OF
           1: BEGIN
              slitmaskred4(where(slitmaskred4 NE 0)) = slitmaskred4(where(slitmaskred4 NE 0)) + 1.0
              slitmaskred9(where(slitmaskred9 NE 0)) = slitmaskred9(where(slitmaskred9 NE 0)) + 1.0
              print, 'First slit missing. OK.'
           END
           2: print, 'Last slit missing. OK'
           ELSE: print, 'Do not understand input.'
        ENDCASE
ENDIF ELSE IF (nslitblue NE nslitred9) OR (nslitblue NE nslitred4) OR (nslitred9 NE nslitred4) THEN BEGIN
	print, 'Number of slits does not match in all 3 slitfiles. Does not look like an easy fix. Exiting.'
	return, 0
ENDIF

;; Use ahk_profile_singlegrating to find the profiles in all the slits in each grating. This procedure saves yfitnorm to a csv file.
;ahk_profile_singlegrating,scifileblue,slitfile=slitfileblue,wavefile=wavefileblue
;ahk_profile_singlegrating,scifilered9,slitfile=slitfilered9,wavefile=wavefilered9
;ahk_profile_singlegrating,scifilered4,slitfile=slitfilered4,wavefile=wavefilered4


;; Now it's time to find the master profile in each slit
;; Define combonormprofile as an array of pointers so each row can have a different length (where each row gives comboyfitnorm for a given slit)
combonormprofile = replicate(Ptr_New(),nslitblue)
combonormprofile = PtrArr(nslitblue, /Allocate_heap)

;This flag will be 1 if padding should go on the right of the slit, 2
;if the padding should go on the left of the slit, and 0 if no padding
;is to be applied.
ccdflagblue = MAKE_ARRAY(nslitblue, /INTEGER, VALUE=0)
ccdflagred = MAKE_ARRAY(nslitblue, /INTEGER, VALUE=0) ;This is the
;same length as nslitblue on purpose, to account for times when a red
;slit is cutoff/missing.

ccdside = MAKE_ARRAY(nslitblue, /INTEGER, VALUE=0)
FOR ii =1,nslitblue DO BEGIN
	ccdloc = MAX(WHERE(slitim EQ ii))
	ccdind = ARRAY_INDICES(slitim, ccdloc)
	IF ccdind[0] LE 2048 THEN ccdside[ii-1]=0 ELSE ccdside[ii-1]=1
ENDFOR

;Blue slit number immediately to left of centre gets padded on the right
ccdflagblue[array_indices(ccdside,max(where(ccdside eq 0)))] = 1.0

;Blue lit number immediately to right of centre gets padded on the left
ccdflagblue[array_indices(ccdside,min(where(ccdside eq 1)))] = 2.0

;First red slit is padded on the left unless the first red slit was
;missing
IF missSlit NE 1 THEN ccdflagred[0] = 2.0

;Last red slit is padded on the right unless the last red slit was
;missing
IF missSlit NE 2 THEN ccdflagred[-1] = 1.0

FOR ii =1,nslitblue DO BEGIN
	undefine, csvblue, slitindexblue, normprofileblue
	search = file_search('normprofile_'+strmid(scifileblue,15,8,/reverse_offset)+'_slitid'+StrTrim(ii,2)+'.csv',count=c1)
	IF (c1 GT 0) THEN BEGIN
		csvblue = read_csv('normprofile_'+strmid(scifileblue,15,8,/reverse_offset)+'_slitid'+StrTrim(ii,2)+'.csv')
		slitindexblue = csvblue.field1
		normprofileblue = csvblue.field2
	ENDIF

	undefine, csvred9, slitindexred9, normprofilered9
	search = file_search('normprofile_'+strmid(scifilered9,15,8,/reverse_offset)+'_slitid'+StrTrim(ii,2)+'.csv',count=c2)
	IF (c2 GT 0) THEN BEGIN
		csvred9 = read_csv('normprofile_'+strmid(scifilered9,15,8,/reverse_offset)+'_slitid'+StrTrim(ii,2)+'.csv')
		slitindexred9 = csvred9.field1
		normprofilered9 = csvred9.field2
	ENDIF

	undefine, csvred4, slitindexred4, normprofilered4
	search = file_search('normprofile_'+strmid(scifilered4,15,8,/reverse_offset)+'_slitid'+StrTrim(ii,2)+'.csv',count=c3)
	IF (c3 GT 0) THEN BEGIN
		csvred4 = read_csv('normprofile_'+strmid(scifilered4,15,8,/reverse_offset)+'_slitid'+StrTrim(ii,2)+'.csv')
		slitindexred4 = csvred4.field1
		normprofilered4 = csvred4.field2
	ENDIF

;*combonormprofile[ii-1]=ahk_comboprofile_nopad(ii,slitindexblue,normprofileblue,slitindexred9,normprofilered9,slitindexred4,normprofilered4)
;*combonormprofile[ii-1]=ahk_comboprofile(ii,slitindexblue,normprofileblue,slitindexred9,normprofilered9,slitindexred4,normprofilered4,ccdside=ccdside[ii-1])
	*combonormprofile[ii-1]=ahk_comboprofile_simplepad(ii,slitindexblue,normprofileblue,slitindexred9,normprofilered9,slitindexred4,normprofilered4,ccdflagblue=ccdflagblue[ii-1],ccdflagred=ccdflagred[ii-1])

ENDFOR ;; end loop over slits to find combonormprofile in each
stop
;; Loop over slits in each grating individually now. Convert combonormprofile to an array, then feed to long_extract_optimal.
;; Blue grating:
FOR ii =1,nslitblue DO BEGIN
	search = file_search('normprofile_'+strmid(scifileblue,15,8,/reverse_offset)+'_slitid'+StrTrim(ii,2)+'.csv',count=c1)
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

		;; Must determine the normalized slitindex. Use our knowledge that combonormprofile abscissa values run from 0 to 1.
		slitindex = findgen(N_ELEMENTS(*combonormprofile[ii-1]))
		normslitindex = slitindex/MAX(slitindex)

		;; Get object profile along slit (i.e. convert this slit's combonormprofile to a 2-d array with non-zero values running along slit length)
		slitprofile = ahk_profile2imarray(normslitindex,*combonormprofile[ii-1],slitid=ii,slitfile=slitfileblue)

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
				;fluxedfileout = imgname+'_slit'+STRTRIM(ii,2)+'_obj'+STRTRIM(objid,2)+'_fluxed.fits'
				;long_fluxcal, fileout, SENSFUNCFILE=sensfuncfileblue, OUTFIL=fluxedfileout, frm_sci=1
			
			ENDFOR ;; end loop over objects

		ENDELSE

ENDFOR ;; end loop over blue slits

;; red9 grating:
FOR ii =1,nslitred9 DO BEGIN
	search = file_search('normprofile_'+strmid(scifilered9,15,8,/reverse_offset)+'_slitid'+StrTrim(ii,2)+'.csv',count=c1)
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

		;; Must determine the normalized slitindex. Use our knowledge that combonormprofile abscissa values run from 0 to 1.
		slitindex = findgen(N_ELEMENTS(*combonormprofile[ii-1]))
		normslitindex = slitindex/MAX(slitindex)

		;; Get object profile along slit (i.e. convert this slit's combonormprofile to a 2-d array with non-zero values running along slit length)
		slitprofile = ahk_profile2imarray(normslitindex,*combonormprofile[ii-1],slitid=ii,slitfile=slitfilered9)

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
				;fluxedfileout = imgname+'_slit'+STRTRIM(ii,2)+'_obj'+STRTRIM(objid,2)+'_fluxed.fits'
				;long_fluxcal, fileout, SENSFUNCFILE=sensfuncfilered9, OUTFIL=fluxedfileout, frm_sci=1
			
			ENDFOR ;; end loop over objects

		ENDELSE

ENDFOR ;; end loop over red9 slits


;; red4 grating:
FOR ii =1,nslitred4 DO BEGIN
	search = file_search('normprofile_'+strmid(scifilered4,15,8,/reverse_offset)+'_slitid'+StrTrim(ii,2)+'.csv',count=c1)
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

		;; Must determine the normalized slitindex. Use our knowledge that combonormprofile abscissa values run from 0 to 1.
		slitindex = findgen(N_ELEMENTS(*combonormprofile[ii-1]))
		normslitindex = slitindex/MAX(slitindex)

		;; Get object profile along slit (i.e. convert this slit's combonormprofile to a 2-d array with non-zero values running along slit length)
		slitprofile = ahk_profile2imarray(normslitindex,*combonormprofile[ii-1],slitid=ii,slitfile=slitfilered4)

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

				;; FUTURE STEP: coadd before flux calibration. This might mean we have to put coadding and fluxing together in a separate script from objextraction.
			
				;; Flux calibration
				;fluxedfileout = imgname+'_slit'+STRTRIM(ii,2)+'_obj'+STRTRIM(objid,2)+'_fluxed.fits'
				;long_fluxcal, fileout, SENSFUNCFILE=sensfuncfilered4, OUTFIL=fluxedfileout, frm_sci=1

			ENDFOR ;; end loop over objects

		ENDELSE

ENDFOR ;; end loop over red4 slits

return, spec

end


