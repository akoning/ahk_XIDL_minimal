function ahk_objextract_multigrating_multiexposure, scifileblue, scifilered9, scifilered4 $
   , slitfileblue=slitfileblue, wavefileblue=wavefileblue $
   , slitfilered9=slitfilered9, wavefilered9=wavefilered9 $
   , slitfilered4=slitfilered4, wavefilered4=wavefilered4


;+
; NAME:
;   ahk_objextract_multigrating_multiexposure
;
; PURPOSE:
;   Intended to be used after long_reduce (or long_reduce_work) and
;   before ahk_fluxandcoadd. Uses first exposure of each grating to
;   create a master slit profile, then uses this master profile to
;   rerun long_extract_optimal on all exposures of all gratings.
;   Changed from
;   ahk_objextract_multigrating_multiexposure_SlitsInDiffFiles to now
;   save the output for all slits of a given mask to one output
;   file. This file's first 5 extensions are the same as
;   before, the last extension is now an array of structures, where
;   each structure is for a different slit+object.
;
;   **Param to watch out for: frm_sci=1 for every call to long_fluxcal. Need to double check if this works correctly/as expected.
;
; INPUTS:
;   scifile****     - file containing object structure output from
;                     long_reduce for all **** exposures. Should
;                     be a list, but can be any length (including 1).
;
;   slitfile****    - slits-****.fits file output from long_reduce.
;
;   wavefile****    - wave-****.fits file output from long_reduce.
;;
; CALLING SEQUENCE:
;
;	bluefiles =
;	['../Science_blue600/sci-lblue2082.fits.gz', '../Science_blue600/sci-lblue2083.fits.gz', '../Science_blue600/sci-lblue2086.fits.gz']
;	red9files =
;	['../Science_red900/sci-lred2068.fits.gz', '../Science_red900/sci-lred2069.fits.gz', '../Science_red900/sci-lred2070.fits.gz']
;	red4files = ['../Science_red400/sci-lred2073.fits.gz']
;
;	result =  ahk_objextract_multigrating_multiexposure(bluefiles,red9files, red4files, slitfileblue='../slits-lblue2039.fits',wavefileblue='../wave-lblue2027.fits',slitfilered9='../slits-lred2031.fits', wavefilered9='../wave-lred2019.fits',slitfilered4='../slits-lred2122.fits',wavefilered4='../wave-lred2158.fits')
;
;       For 2008nov06/M33r3_7:
;IDL> bluefiles=['../Science_blue600/sci-lblue6099.fits.gz','../Science_blue600/sci-lblue6100.fits.gz','../Science_blue600/sci-lblue6103.fits.gz']
;IDL> red9files=['../Science_red900/sci-lred6090.fits.gz','../Science_red900/sci-lred6091.fits.gz']
;IDL> red4files=['../Science_red400/sci-lred6094.fits.gz']
;IDL> result =  ahk_objextract_multigrating_multiexposure(bluefiles,red9files, red4files, slitfileblue='../slits-lblue6053.fits',wavefileblue='../wave-lblue6035.fits',slitfilered9='../slits-lred6046.fits', wavefilered9='../wave-lred6026.fits',slitfilered4='../slits-lred6066.fits',wavefilered4='../wave-lred6036.fits')
;
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
IF KEYWORD_SET(slitfileblue) THEN slitim=xmrdfits(slitfileblue,0) $
ELSE IF KEYWORD_SET(slitfilered9) THEN slitim=xmrdfits(slitfilered9,0) $
ELSE IF KEYWORD_SET(slitfilered4) THEN slitim=xmrdfits(slitfilered4,0) $
ELSE print, 'WARNING! No slitfiles are defined'

;; Number of slits in image
IF KEYWORD_SET(slitfileblue) THEN BEGIN
   slitmaskblue = mrdfits(slitfileblue)
   nslitblue = max(slitmaskblue,/NaN)
ENDIF ELSE nslitblue = 0
IF KEYWORD_SET(slitfilered9) THEN BEGIN
   slitmaskred9 = mrdfits(slitfilered9)
   nslitred9 = max(slitmaskred9,/NaN)
ENDIF ELSE nslitred9 = 0
IF KEYWORD_SET(slitfilered4) THEN BEGIN
   slitmaskred4 = mrdfits(slitfilered4)
   nslitred4 = max(slitmaskred4,/NaN)
ENDIF ELSE nslitred4 = 0

IF nslitblue NE 0 THEN nslits = nslitblue $
ELSE IF nslitred9 NE 0 THEN nslits = nslitred9 $
ELSE IF nslitred4 NE 0 THEN nslits = nslitred4

;;Search for case when red gratings have different number of slits
;;than blue. Ask user to specify which slit is missing in the red and
;;adjust the slitmasks accordingly.
missSlit=0   
IF (nslitred9 EQ nslitred4) AND (nslitblue NE nslitred9) THEN BEGIN
   read, 'Enter 1 if first slit in red missing, enter 2 if last slit in red missing, enter 3 if both the first and last slits are missing, enter 4 if first TWO slits are missing: ',missSlit
        CASE missSlit OF
           1: BEGIN
              slitmaskred4(where(slitmaskred4 NE 0)) = slitmaskred4(where(slitmaskred4 NE 0)) + 1.0
              slitmaskred9(where(slitmaskred9 NE 0)) = slitmaskred9(where(slitmaskred9 NE 0)) + 1.0
              ;nslitred4 += 1
              ;nslitred9 += 1
              print, 'First slit missing. OK.'
           END
           2: BEGIN
              ;nslitred4 += 1
              ;nslitred9 += 1
              print, 'Last slit missing. OK'
           END
           3: BEGIN
              slitmaskred4(where(slitmaskred4 NE 0)) = slitmaskred4(where(slitmaskred4 NE 0)) + 1.0
              slitmaskred9(where(slitmaskred9 NE 0)) = slitmaskred9(where(slitmaskred9 NE 0)) + 1.0
              ;nslitred4 += 2
              ;nslitred9 += 2
              print, 'First and last slits missing. OK.'
           END
           4: BEGIN
              slitmaskred4(where(slitmaskred4 NE 0)) = slitmaskred4(where(slitmaskred4 NE 0)) + 2.0
              slitmaskred9(where(slitmaskred9 NE 0)) = slitmaskred9(where(slitmaskred9 NE 0)) + 2.0
              ;nslitred4 += 2
              ;nslitred9 += 2
              print, 'First two slits missing. OK.'
           END
           ELSE: print, 'Do not understand input.'
        ENDCASE
ENDIF ELSE IF (nslitred9 NE nslitred4) AND (nslitred9 EQ 0) THEN BEGIN
   print, 'Red900 not defined.'
   IF (nslitblue NE nslitred4) THEN BEGIN
      read, 'Enter 1 if first slit in red missing, enter 2 if last slit in red missing, enter 3 if both the first and last slits are missing, enter 4 if first TWO slits are missing: ',missSlit
      CASE missSlit OF
         1: BEGIN
            slitmaskred4(where(slitmaskred4 NE 0)) = slitmaskred4(where(slitmaskred4 NE 0)) + 1.0
            print, 'First slit missing. OK.'
         END
         2: BEGIN
            print, 'Last slit missing. OK'
         END
         3: BEGIN
            slitmaskred4(where(slitmaskred4 NE 0)) = slitmaskred4(where(slitmaskred4 NE 0)) + 1.0
            print, 'First and last slits missing. OK.'
         END
         4: BEGIN
            slitmaskred4(where(slitmaskred4 NE 0)) = slitmaskred4(where(slitmaskred4 NE 0)) + 2.0
            print, 'First two slits missing. OK.'
         END
         ELSE: print, 'Do not understand input.'
      ENDCASE
   ENDIF
ENDIF ELSE IF (nslitred9 NE nslitred4) AND (nslitred4 EQ 0) THEN BEGIN
   print, 'Red400 not defined.'
   IF (nslitblue NE nslitred9) THEN BEGIN
      read, 'Enter 1 if first slit in red missing, enter 2 if last slit in red missing, enter 3 if both the first and last slits are missing, enter 4 if first TWO slits are missing: ',missSlit
      CASE missSlit OF
         1: BEGIN
            slitmaskred9(where(slitmaskred9 NE 0)) = slitmaskred9(where(slitmaskred9 NE 0)) + 1.0
            print, 'First slit missing. OK.'
         END
         2: BEGIN
            print, 'Last slit missing. OK'
         END
         3: BEGIN
            slitmaskred9(where(slitmaskred9 NE 0)) = slitmaskred9(where(slitmaskred9 NE 0)) + 1.0
            print, 'First and last slits missing. OK.'
         END
         4: BEGIN
            slitmaskred9(where(slitmaskred9 NE 0)) = slitmaskred9(where(slitmaskred9 NE 0)) + 2.0
            print, 'First two slits missing. OK.'
         END
         ELSE: print, 'Do not understand input.'
      ENDCASE
   ENDIF

ENDIF ELSE IF (nslitred9 NE nslitred4) THEN BEGIN
   print, 'Number of slits does not match in the two red slitfiles. Requires some human help.'
   stop
ENDIF

;; Use ahk_profile_singlegrating to find the profiles in all the slits
;; in the first exposure of each grating. This procedure saves
;; yfitnorm to a csv file.
IF KEYWORD_SET(slitfileblue) THEN ahk_profile_singlegrating,scifileblue[0],slitmask=slitmaskblue,wavefile=wavefileblue,missSlit=0
IF KEYWORD_SET(slitfilered9) THEN ahk_profile_singlegrating,scifilered9[0],slitmask=slitmaskred9,wavefile=wavefilered9,missSlit=missSlit
IF KEYWORD_SET(slitfilered4) THEN ahk_profile_singlegrating,scifilered4[0],slitmask=slitmaskred4,wavefile=wavefilered4,missSlit=missSlit


;; Now it's time to find the master profile in each slit
;; Define combonormprofile as an array of pointers with a different vector of comboprofile yvalues for each grating AND each slit
;; Blue600 = row index 0, Red900 = row index 1, Red400 = row index 2
combonormprofile = replicate(Ptr_New(),nslits)
combonormprofile = PtrArr(3, nslits, /Allocate_heap)

;This flag will be 1 if padding should go on the right of the slit, 2
;if the padding should go on the left of the slit, and 0 if no padding
;is to be applied.
ccdflagblue = MAKE_ARRAY(nslits, /INTEGER, VALUE=0)
ccdflagred = MAKE_ARRAY(nslits, /INTEGER, VALUE=0)

ccdside = MAKE_ARRAY(nslits, /INTEGER, VALUE=0)
FOR ii =1,nslits DO BEGIN
	ccdloc = MAX(WHERE(slitim EQ ii))
	ccdind = ARRAY_INDICES(slitim, ccdloc)
	IF ccdind[0] LE ccdind[1]/2 THEN ccdside[ii-1]=0 ELSE ccdside[ii-1]=1
ENDFOR

;Blue slit number immediately to left of centre gets padded on the right
ccdflagblue[array_indices(ccdside,max(where(ccdside eq 0)))] = 1.0

;Blue slit number immediately to right of centre gets padded on the left
ccdflagblue[array_indices(ccdside,min(where(ccdside eq 1)))] = 2.0

;First red slit that is not missing is padded on the left
IF missSlit EQ 1 OR missSlit EQ 3 THEN ccdflagred[1] = 2.0 $
ELSE IF missSlit EQ 4 THEN ccdflagred[2] = 2.0 $
ELSE ccdflagred[0] = 2.0

;Last red slit that is not missing is padded on the right
IF missSlit EQ 2 OR missSlit EQ 3 THEN ccdflagred[-2] = 1.0 $
ELSE ccdflagred[-1] = 1.0


; Loop over slits to find combonormprofile in each 
FOR ii =1,nslits DO BEGIN
	undefine, csvblue, slitindexblue, normprofileblue
	search = file_search('normprofile_'+strmid(scifileblue[0],15,8,/reverse_offset)+'_slitid'+StrTrim(ii,2)+'.csv',count=c1)
	IF (c1 GT 0) THEN BEGIN
		csvblue = read_csv('normprofile_'+strmid(scifileblue[0],15,8,/reverse_offset)+'_slitid'+StrTrim(ii,2)+'.csv')
		slitindexblue = csvblue.field1
		normprofileblue = csvblue.field2
	ENDIF

	undefine, csvred9, slitindexred9, normprofilered9
	search = file_search('normprofile_'+strmid(scifilered9[0],15,8,/reverse_offset)+'_slitid'+StrTrim(ii,2)+'.csv',count=c2)
	IF (c2 GT 0) THEN BEGIN
		csvred9 = read_csv('normprofile_'+strmid(scifilered9[0],15,8,/reverse_offset)+'_slitid'+StrTrim(ii,2)+'.csv')
		slitindexred9 = csvred9.field1
		normprofilered9 = csvred9.field2
	ENDIF

	undefine, csvred4, slitindexred4, normprofilered4
	search = file_search('normprofile_'+strmid(scifilered4[0],15,8,/reverse_offset)+'_slitid'+StrTrim(ii,2)+'.csv',count=c3)
	IF (c3 GT 0) THEN BEGIN
		csvred4 = read_csv('normprofile_'+strmid(scifilered4[0],15,8,/reverse_offset)+'_slitid'+StrTrim(ii,2)+'.csv')
		slitindexred4 = csvred4.field1
		normprofilered4 = csvred4.field2
	ENDIF
	
        goodfiles=[ISA(slitindexblue),ISA(slitindexred9),ISA(slitindexred4)]
        igood = where(goodfiles NE 0)

	FOREACH igrating, igood DO $
	*combonormprofile[igrating,ii-1]=ahk_comboprofile_simplepad(ii,slitindexblue,normprofileblue,slitindexred9,normprofilered9,slitindexred4,normprofilered4, $
	ccdflagblue=ccdflagblue[ii-1],ccdflagred=ccdflagred[ii-1], igrating=igrating)

ENDFOR ;; end loop over slits to find combonormprofile in each
save, combonormprofile, filename='combonormprofile.sav'
save, /all, filename='investigatingCrap.sav' ;and set nslitblue=0 so it will skip to red900
;; Loop over slits in each exposure in each grating individually now. Convert combonormprofile to an array, then feed to long_extract_optimal.
;; Blue grating:
IF nslitblue EQ 0 THEN GOTO, Red900
nexpblue = (size(scifileblue))[1]
igrating = 0

FOR jj = 0 , nexpblue-1 DO BEGIN 
   ;; Get various required extensions out of science file
   ivar = xmrdfits(scifileblue[jj],1)
   skyimage = xmrdfits(scifileblue[jj],2)
   objimage = xmrdfits(scifileblue[jj],3)
   outmask = xmrdfits(scifileblue[jj],4)
   waveimg = xmrdfits(wavefileblue, silent = (keyword_set(verbose) EQ 0), 0)
   sciimg = xmrdfits(scifileblue[jj])
   imgminsky = sciimg - skyimage

   ;Search to make sure this hasn't already been done, if it has, do not overwrite
   imgname = STRMID(scifileblue[jj], 15, 8, /REVERSE_OFFSET)
   fileout = imgname+'_allSlitsReduced.fits'
   search = file_search(fileout,count=c2)
   IF c2 GT 0 THEN BEGIN
      print, 'Do not overwrite existing file '+fileout
      CONTINUE
   ENDIF

   ;; Output all fits extensions as they were, except last science
   ;; extension b/c that will contain an array of structures
   rdfits_struct,scifileblue[jj],scistr ;; read original scifile into structure
   mwrfits,scistr.im0,fileout,scistr.hdr0,/create
   mwrfits,scistr.im1,fileout,scistr.hdr1 
   mwrfits,scistr.im2,fileout,scistr.hdr2 
   mwrfits,scistr.im3,fileout,scistr.hdr3 
   mwrfits,scistr.im4,fileout,scistr.hdr4 


   FOR ii =1,nslitblue DO BEGIN
      search = file_search('normprofile_'+strmid(scifileblue[0],15,8,/reverse_offset)+'_slitid'+StrTrim(ii,2)+'.csv',count=c1)
      ;; If CSV file exists, and therefore there is flux in that slit, then continue with the analysis. Otherwise skip to next slit.
      IF (c1 LT 1) THEN CONTINUE
      print, '***** Working on BLUE EXPOSURE ', jj, ', SLITID ', ii, ' *****'


      ;; Size of science image
      nx = (size(sciimg))[1]
      ny = (size(sciimg))[2]
      xarr = findgen(nx)## replicate(1.0, ny)
      yarr = findgen(ny)## replicate(1.0, nx)

      ;; Must determine the normalized slitindex. Use our knowledge that combonormprofile abscissa values run from 0 to 1.
      slitindex = findgen(N_ELEMENTS(*combonormprofile[igrating,ii-1]))
      normslitindex = slitindex/MAX(slitindex)

      ;; Get object profile along slit (i.e. convert this slit's combonormprofile to a 2-d array with non-zero values running along slit length)
      slitprofile = ahk_profile2imarray(normslitindex,*combonormprofile[igrating,ii-1],slitid=ii,slitfile=slitfileblue)

      ;; Number of objects in slit. If more than one, must get 1D spectrum for each.
      ymask = (slitprofile GT 0.01) ;; Assumes yfit is normalized.
      label = label_region(ymask)
      nobj = max(label)
      print, 'Number of objects in slit: ', nobj

      IF nobj EQ 0 OR nobj GT 5 THEN BEGIN
         print, 'No objects in profile GT 0.01, or too many objects (something went wrong).'
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
            ;; Add slitid and objid tags to spec
            ;; structure
            newspec = add_slitid_objid_tags(spec, slitid=ii, objid=objid)
            ;; Concatenate spec (a structure) to
            ;; the array of structures that will
            ;; be output in fileout.
            IF ii EQ 1 and objid EQ 1 THEN specarray = newspec ELSE specarray = [specarray, newspec]
            
         ENDFOR ;; end loop over objects

      ENDELSE
   ENDFOR ;; end loop over blue slits
   mwrfits,specarray,fileout,hdr_struct
   spawn, 'gzip -f '+fileout
ENDFOR ;; end loop over blue exposures

;; red9 grating:
Red900: print, 'Starting Red900 reduction now'
IF nslitred9 EQ 0 THEN GOTO, Red400
nexpred9 = (size(scifilered9))[1]
igrating = 1

FOR jj = 0 , nexpred9-1 DO BEGIN 
   ;; Get various required extensions out of science file
   ivar = xmrdfits(scifilered9[jj],1)
   skyimage = xmrdfits(scifilered9[jj],2)
   objimage = xmrdfits(scifilered9[jj],3)
   outmask = xmrdfits(scifilered9[jj],4)
   waveimg = xmrdfits(wavefilered9, silent = (keyword_set(verbose) EQ 0), 0)
   sciimg = xmrdfits(scifilered9[jj])
   imgminsky = sciimg - skyimage
   stop ;AHK check for negative numbers

   ;Search to make sure this hasn't already been done, if it has, do not overwrite
   imgname = STRMID(scifilered9[jj], 15, 8, /REVERSE_OFFSET)
   fileout = imgname+'_allSlitsReduced.fits'
   search = file_search(fileout,count=c2)
   IF c2 GT 0 THEN BEGIN
      print, 'Do not overwrite existing file '+fileout
      CONTINUE
   ENDIF

   ;; Output all fits extensions as they were, except last science
   ;; extension b/c that will contain an array of structures
   rdfits_struct,scifilered9[jj],scistr ;; read original scifile into structure
   mwrfits,scistr.im0,fileout,scistr.hdr0,/create
   mwrfits,scistr.im1,fileout,scistr.hdr1 
   mwrfits,scistr.im2,fileout,scistr.hdr2 
   mwrfits,scistr.im3,fileout,scistr.hdr3 
   mwrfits,scistr.im4,fileout,scistr.hdr4

   FOR ii =1,nslitred9 DO BEGIN
        slitid = ii
        IF missSlit EQ 1 OR missSlit EQ 3 THEN slitid = ii + 1
        IF missSlit EQ 4 THEN slitid = ii + 2

	search = file_search('normprofile_'+strmid(scifilered9[0],15,8,/reverse_offset)+'_slitid'+StrTrim(slitid,2)+'.csv',count=c1)
	;; If CSV file exists, and therefore there is flux in that slit, then continue with the analysis. Otherwise skip to next slit.
	IF (c1 LT 1) THEN CONTINUE
        print, '***** Working on RED9 EXPOSURE ', jj, ', SLITID ', slitid, ', ii ', ii, ' *****'

        ;; Size of science image
        nx = (size(sciimg))[1]
        ny = (size(sciimg))[2]
        xarr = findgen(nx)## replicate(1.0, ny)
        yarr = findgen(ny)## replicate(1.0, nx)

        ;; Must determine the normalized slitindex. Use our knowledge that combonormprofile abscissa values run from 0 to 1.
        slitindex = findgen(N_ELEMENTS(*combonormprofile[igrating,slitid-1]))
        normslitindex = slitindex/MAX(slitindex)

        ;; Get object profile along slit (i.e. convert this slit's combonormprofile to a 2-d array with non-zero values running along slit length)
        slitprofile = ahk_profile2imarray(normslitindex,*combonormprofile[igrating,slitid-1],slitid=ii,slitfile=slitfilered9)
        stop ;check profile
        ;; Number of objects in slit. If more than one, must get 1D spectrum for each.
        ymask = (slitprofile GT 0.01) ;; Assumes yfit is normalized.
        label = label_region(ymask)
        nobj = max(label)
        print, 'Number of objects in slit: ', nobj

        IF nobj EQ 0 OR nobj GT 5 THEN BEGIN
           print, 'No objects in profile GT 0.01, or too many objects (something went wrong).'
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
              stop ;check stuff before long extract optimal
              spec = long_extract_optimal(waveimg, imgminsky, ivar, slitprofile, outmask*objmask, skyimage, trace, BOX_RAD=box_rad)
              stop ;check spec
              ;; Add slitid and objid tags to spec
              ;; structure
              newspec = add_slitid_objid_tags(spec, slitid=ii, objid=objid)
              ;; Concatenate spec (a structure) to
              ;; the array of structures that will
              ;; be output in fileout.
              IF ii EQ 1 and objid EQ 1 THEN specarray = newspec ELSE specarray = [specarray, newspec]

           ENDFOR ;; end loop over objects

        ENDELSE

     ENDFOR ;; end loop over red9 slits
   mwrfits,specarray,fileout,hdr_struct			
   spawn, 'gzip -f '+fileout
ENDFOR ;; end loop over red9 exposures

;; red4 grating:
Red400: print, 'Starting Red400 reduction now'
IF nslitred4 EQ 0 THEN GOTO, FINISH
nexpred4 = (size(scifilered4))[1]
igrating = 2

FOR jj = 0 , nexpred4-1 DO BEGIN 
   ;; Get various required extensions out of science file
   ivar = xmrdfits(scifilered4[jj],1)
   skyimage = xmrdfits(scifilered4[jj],2)
   objimage = xmrdfits(scifilered4[jj],3)
   outmask = xmrdfits(scifilered4[jj],4)
   waveimg = xmrdfits(wavefilered4, silent = (keyword_set(verbose) EQ 0), 0)
   sciimg = xmrdfits(scifilered4[jj])
   imgminsky = sciimg - skyimage

   ;Search to make sure this hasn't already been done, if it has, do not overwrite
   imgname = STRMID(scifilered4[jj], 15, 8, /REVERSE_OFFSET)
   fileout = imgname+'_allSlitsReduced.fits'
   search = file_search(fileout,count=c2)
   IF c2 GT 0 THEN BEGIN
      print, 'Do not overwrite existing file '+fileout
      CONTINUE
   ENDIF

   ;; Output all images etc as they were, but now last science extension contains spec
   rdfits_struct,scifilered4[jj],scistr ;; read original scifile into structure
   mwrfits,scistr.im0,fileout,scistr.hdr0,/create
   mwrfits,scistr.im1,fileout,scistr.hdr1 
   mwrfits,scistr.im2,fileout,scistr.hdr2 
   mwrfits,scistr.im3,fileout,scistr.hdr3 
   mwrfits,scistr.im4,fileout,scistr.hdr4 

   FOR ii =1,nslitred4 DO BEGIN
        slitid = ii
        IF missSlit EQ 1 OR missSlit EQ 3 THEN slitid = ii + 1
        IF missSlit EQ 4 THEN slitid = ii + 2
        
	search = file_search('normprofile_'+strmid(scifilered4[0],15,8,/reverse_offset)+'_slitid'+StrTrim(slitid,2)+'.csv',count=c1)
	;; If CSV file exists, and therefore there is flux in that slit, then continue with the analysis. Otherwise skip to next slit.
        IF (c1 LT 1) THEN CONTINUE
        print, '***** Working on RED4 EXPOSURE ', jj, ', SLITID ', slitid, ', ii ', ii, ' *****'
        

        ;; Size of science image
        nx = (size(sciimg))[1]
        ny = (size(sciimg))[2]
        xarr = findgen(nx)## replicate(1.0, ny)
        yarr = findgen(ny)## replicate(1.0, nx)

        ;; Must determine the normalized slitindex. Use our knowledge that combonormprofile abscissa values run from 0 to 1.
        slitindex = findgen(N_ELEMENTS(*combonormprofile[igrating,slitid-1]))
        normslitindex = slitindex/MAX(slitindex)

        ;; Get object profile along slit (i.e. convert this slit's combonormprofile to a 2-d array with non-zero values running along slit length)
        slitprofile = ahk_profile2imarray(normslitindex,*combonormprofile[igrating,slitid-1],slitid=ii,slitfile=slitfilered4)

        ;; Number of objects in slit. If more than one, must get 1D spectrum for each.
        ymask = (slitprofile GT 0.01) ;; Assumes yfit is normalized.
        label = label_region(ymask)
        nobj = max(label)
        print, 'Number of objects in slit: ', nobj

        IF nobj EQ 0 OR nobj GT 5 THEN BEGIN
           print, 'No objects in profile GT 0.01, or too many objects (something went wrong).'
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
              ;; Add slitid and objid tags to spec
              ;; structure
              newspec = add_slitid_objid_tags(spec, slitid=ii, objid=objid)
              ;; Concatenate spec (a structure) to
              ;; the array of structures that will
              ;; be output in fileout.
              IF ii EQ 1 and objid EQ 1 THEN specarray = newspec ELSE specarray = [specarray, newspec]
              

           ENDFOR ;; end loop over objects

        ENDELSE

     ENDFOR ;; end loop over red4 slits
   mwrfits,specarray,fileout,hdr_struct
   spawn, 'gzip -f '+fileout
ENDFOR ;; end loop over red4 exposures
FINISH: print, 'Finished'
return, 'done'

end


