function ahk_objextract_multigrating,sciblue_1,scired9_1,scired4_1,sciblue_2=sciblue_2,sciblue_3=sciblue_3,slitfileblue=slitfileblue,wavefileblue=wavefileblue,scired9_2=scired9_2,scired9_3=scired9_3,slitfilered9=slitfilered9,wavefilered9=wavefilered9,slitfilered4=slitfilered4,wavefilered4=wavefilered4

;+
; NAME:
;   ahk_objextract_multigrating
;
; PURPOSE:
;   Uses output of long_reduce. First fit Gaussian model across each slit profile, in each grating. Combine these profiles to get a master profile.
;   Create a 2D array of this profile running along the slit, then extract objects accordingly with long_extract_optimal.
;   General outline
;	1) Determine number of slits in each grating from slitfileblue, slitfilered9, slitfilered4. If not equal print error message and exit
;	-- Steps which are also repeated for each grating:
;	2) Restore saved objstr.FLUXMODEL (1-d array of summed fluxes in our lines of interest) -- done by ahk_profile_singlegrating
;	3) Fit FLUXMODEL with sum of Gaussians and save csv of best normalized fit -- done by ahk_profile_singleslit.pro within ahk_profile_singlegrating.
;	-- Steps in loop over slits:
;	4) Combine all DEFINED profiles (i.e. csv file exists) to get master profile for each slit using ahk_comboprofile.pro
;	-- Steps in new loop over slits, repeated for each grating:
;	5) Create 4096x4096 profile -- done by ahk_profile2imarray.pro (CAREFUL: should this step go in loop over just slits??)
;	6) Feed profile into long_extract_optimal
;	7) Replace last fits extension in science image with output from long_extract_optimal
;
; CALLING SEQUENCE:
;	result =  ahk_objextract_multigrating('Science_blue600/sci-lblue2082.fits.gz','Science_red900/sci-lred2068.fits.gz','Science_red400/sci-lred2073.fits.gz',slitfileblue='slits-lblue2039.fits',wavefileblue='wave-lblue2027.fits',slitfilered9='slits-lred2031.fits',wavefilered9='wave-lred2158.fits',slitfilered4='slits-lred2122.fits',wavefilered4='wave-lred2019.fits')

;	result =  ahk_objextract_multigrating('Science_blue600/sci-lblue2082.fits.gz','Science_red900/sci-lred2068.fits.gz','Science_red400/sci-lred2073.fits.gz',sciblue_2='Science_blue600/sci-lblue2083.fits.gz',sciblue_3='Science_blue600/sci-lblue2086.fits.gz',slitfileblue='slits-lblue2039.fits',wavefileblue='wave-lblue2027.fits',scired9_2='Science_red900/sci-lred2069.fits.gz',scired9_3='Science_red900/sci-lred2070.fits.gz',slitfilered9='slits-lred2031.fits',wavefilered9='wave-lred2158.fits',slitfilered4='slits-lred2122blankit.fits',wavefilered4='wave-lred2019.fits')

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
; Assumes there are 3 exposures in blue600, 3 in red900, and 1 in red400. Uses first exposure in each to calculate comboyfit, then reuses that for other 2 exposures in blue600 and red900. If there are different numbers of exposures in each, will have to adjust code. If we want to calculate comboyfit from all (e.g.) 7 exposures at once, will have to adjust code.
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

;; Use ahk_profile_singlegrating to find the profiles in all the slits in each grating. This procedure saves yfitnorm to a csv file.
ahk_profile_singlegrating,sciblue_1,slitfile=slitfileblue,wavefile=wavefileblue
ahk_profile_singlegrating,scired9_1,slitfile=slitfilered9,wavefile=wavefilered9
ahk_profile_singlegrating,scired4_1,slitfile=slitfilered4,wavefile=wavefilered4


;; Now it's time to find the master profile in each slit
;; Define comboyfitnorm as an array of pointers so each row can have a different length (where each row gives comboyfitnorm for a given slit)
comboyfitnorm = replicate(Ptr_New(),nslitblue)
comboyfitnorm = PtrArr(nslitblue, /Allocate_heap)

FOR ii =1,nslitblue DO BEGIN
	ccdloc = MAX(WHERE(slitim EQ ii))
	ccdind = ARRAY_INDICES(slitim, ccdloc)
	IF ccdind[0] LE 2048 THEN ccdside=0 ELSE ccdside=1
	undefine, csvblue, slitindexblue, yfitnormblue
	search = file_search('yfitnorm_'+strmid(sciblue_1,15,8,/reverse_offset)+'_slitid'+StrTrim(ii,2)+'.csv',count=c1)
	IF (c1 GT 0) THEN BEGIN
		csvblue = read_csv('yfitnorm_'+strmid(sciblue_1,15,8,/reverse_offset)+'_slitid'+StrTrim(ii,2)+'.csv')
		slitindexblue = csvblue.field1
		yfitnormblue = csvblue.field2
	ENDIF

	undefine, csvred9, slitindexred9, yfitnormred9
	search = file_search('yfitnorm_'+strmid(scired9_1,15,8,/reverse_offset)+'_slitid'+StrTrim(ii,2)+'.csv',count=c2)
	IF (c2 GT 0) THEN BEGIN
		csvred9 = read_csv('yfitnorm_'+strmid(scired9_1,15,8,/reverse_offset)+'_slitid'+StrTrim(ii,2)+'.csv')
		slitindexred9 = csvred9.field1
		yfitnormred9 = csvred9.field2
	ENDIF

	undefine, csvred4, slitindexred4, yfitnormred4
	search = file_search('yfitnorm_'+strmid(scired4_1,15,8,/reverse_offset)+'_slitid'+StrTrim(ii,2)+'.csv',count=c3)
	IF (c3 GT 0) THEN BEGIN
		csvred4 = read_csv('yfitnorm_'+strmid(scired4_1,15,8,/reverse_offset)+'_slitid'+StrTrim(ii,2)+'.csv')
		slitindexred4 = csvred4.field1
		yfitnormred4 = csvred4.field2
	ENDIF

	;*comboyfitnorm[ii-1]=ahk_comboprofile_nopad(ii,slitindexblue,yfitnormblue,slitindexred9,yfitnormred9,slitindexred4,yfitnormred4)
	*comboyfitnorm[ii-1]=ahk_comboprofile(ii,slitindexblue,yfitnormblue,slitindexred9,yfitnormred9,slitindexred4,yfitnormred4,ccdside=ccdside)

ENDFOR ;; end loop over slits to find comboyfitnorm in each

;; Loop over slits in each grating individually now. Convert comboyfitnorm to an array, then feed to long_extract_optimal.
;; Blue grating:
FOR ii =1,nslitblue DO BEGIN
	FOREACH scifile, [sciblue_1,sciblue_2,sciblue_3] DO BEGIN
		search = file_search('yfitnorm_'+strmid(sciblue_1,15,8,/reverse_offset)+'_slitid'+StrTrim(ii,2)+'.csv',count=c1)
		;; If CSV file exists, and therefore there is flux in that slit, then continue with the analysis. Otherwise skip to next slit.
		IF (c1 LT 1) THEN CONTINUE
		print, '***** Working on BLUE SLITID: ', ii, ' *****'

		;; Get various required extensions out of science file
		ivar = xmrdfits(scifile,1)
		skyimage = xmrdfits(scifile,2)
		objimage = xmrdfits(scifile,3)
		outmask = xmrdfits(scifile,4)
		waveimg = xmrdfits(wavefileblue, silent = (keyword_set(verbose) EQ 0), 0)
		sciimg = xmrdfits(scifile)
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
				rdfits_struct,scifile,scistr ;; read original scifile into structure
				imgname = STRMID(scifile, 15, 8, /REVERSE_OFFSET)
				fileout = imgname+'_slit'+STRTRIM(ii,2)+'_obj'+STRTRIM(objid,2)+'.fits'
				mwrfits,scistr.im0,fileout,scistr.hdr0,/create
				mwrfits,scistr.im1,fileout,scistr.hdr1 
				mwrfits,scistr.im2,fileout,scistr.hdr2 
				mwrfits,scistr.im3,fileout,scistr.hdr3 
				mwrfits,scistr.im4,fileout,scistr.hdr4 
				mwrfits,spec,fileout,hdr_struct
			
			ENDFOR ;; end loop over objects
		ENDELSE
	ENDFOREACH
ENDFOR ;; end loop over blue slits

;; red9 grating:
FOR ii =1,nslitred9 DO BEGIN
	FOREACH scifile, [scired9_1,scired9_2,scired9_3] DO BEGIN
		search = file_search('yfitnorm_'+strmid(scired9_1,15,8,/reverse_offset)+'_slitid'+StrTrim(ii,2)+'.csv',count=c1)
		;; If CSV file exists, and therefore there is flux in that slit, then continue with the analysis. Otherwise skip to next slit.
		IF (c1 LT 1) THEN CONTINUE
		print, '***** Working on RED9 SLITID: ', ii, ' *****'

		;; Get various required extensions out of science file
		ivar = xmrdfits(scifile,1)
		skyimage = xmrdfits(scifile,2)
		objimage = xmrdfits(scifile,3)
		outmask = xmrdfits(scifile,4)
		waveimg = xmrdfits(wavefilered9, silent = (keyword_set(verbose) EQ 0), 0)
		sciimg = xmrdfits(scifile)
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
				rdfits_struct,scifile,scistr ;; read original scifile into structure
				imgname = STRMID(scifile, 15, 8, /REVERSE_OFFSET)
				fileout = imgname+'_slit'+STRTRIM(ii,2)+'_obj'+STRTRIM(objid,2)+'.fits'
				mwrfits,scistr.im0,fileout,scistr.hdr0,/create
				mwrfits,scistr.im1,fileout,scistr.hdr1 
				mwrfits,scistr.im2,fileout,scistr.hdr2 
				mwrfits,scistr.im3,fileout,scistr.hdr3 
				mwrfits,scistr.im4,fileout,scistr.hdr4 
				mwrfits,spec,fileout,hdr_struct
			
			ENDFOR ;; end loop over objects

		ENDELSE
	ENDFOREACH

ENDFOR ;; end loop over red9 slits

;; red4 grating:
FOR ii =1,nslitred4 DO BEGIN
	search = file_search('yfitnorm_'+strmid(scired4_1,15,8,/reverse_offset)+'_slitid'+StrTrim(ii,2)+'.csv',count=c1)
	;; If CSV file exists, and therefore there is flux in that slit, then continue with the analysis. Otherwise skip to next slit.
	IF (c1 LT 1) THEN CONTINUE
		print, '***** Working on RED4 SLITID: ', ii, ' *****'

		;; Get various required extensions out of science file
		ivar = xmrdfits(scired4_1,1)
		skyimage = xmrdfits(scired4_1,2)
		objimage = xmrdfits(scired4_1,3)
		outmask = xmrdfits(scired4_1,4)
		waveimg = xmrdfits(wavefilered4, silent = (keyword_set(verbose) EQ 0), 0)
		sciimg = xmrdfits(scired4_1)
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
				rdfits_struct,scired4_1,scistr ;; read original scifile into structure
				imgname = STRMID(scired4_1, 15, 8, /REVERSE_OFFSET)
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


