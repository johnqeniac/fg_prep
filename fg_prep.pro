
;-----------------------------------------------------------------------------
PRO fg_prep, input1, input2, index_out, data_out,                      $
             darkdir=darkdir,	          user_dark=user_dark,	       $
             user_dindex=user_dindex,     dark_image=dark_image,       $
             dark_index=dark_index,                                    $
             flatdir=flatdir, 	          user_flat=user_flat,         $
             user_findex=user_findex,     flat_image=flat_image,       $
             flat_index=flat_index,                                    $
             no_shiftpix=no_shiftpix,     shiftscale=shiftscale,       $
             no_roll_scale=no_roll_scale, do_roll_scale=do_roll_scale, $
             no_do_roll=no_do_roll,                                    $
             no_badpix=no_badpix,         no_darksub=no_darksub,       $
             no_flatfield=no_flatfield,   no_pointing=no_pointing,     $
             tf_deripple=tf_deripple,                                  $
             redspot=redspot,             doppler=doppler,             $
             polarcal=polarcal,                                        $
             despike=despike,             nofloat=nofloat,             $
             x0=x0, x1=x1,                y0=y0, y1=y1,                $
             subimgx=subimgx,             subimgy=subimgy,             $
             center=center,                                            $
             nodata=nodata,               no_calib=no_calib,	       $
             original=original,	                                       $
             outdir=outdir,	          outflatfits=outflatfits,     $
             outfiletemplate=outfiletemplate,  prefix=prefix,          $
             display=display,             run_time=run_time,           $
             topdir_genx=topdir_genx,     corr_db_good=corr_db_good,   $
             n_corr_db_good=n_corr_db_good,                            $
             qstop=qstop, quiet=quiet, verbose=verbose,                $
             _extra=_extra, version=progver, name=prognam

;+
; NAME:
;	FG_PREP
;
; PURPOSE:
;	Process an SOT FG dataset: BFI or NFI filtergram, magnetogram, dopplergram, or Stokes set.
;
;       The steps performed are:
;		1. Read in FITS file(s) from a filelist
;		     or read in a datacube and structure array.
;               For each data file:
;               a. Identify type of data: BFI or NFI image or NFI
;               image set (e.g. Stokes IV [m,n,2] array).

;               b. Correct camera readout defects (e.g. missing rows at center).

;               c. Subtract ADC offset and dark current from the frame.
;               d. Muliply by gain image (1/flatfield image).
;               e. Optional: remove radiation-belt/cosmic-ray spikes and 
; 		     streaks.
;               f. Optional: remove BFI 6684 red spot. NYI.
;               g. Optional NFI: Correct for NFI tunable filter
;                  I-ripple. NYI.
;		h. Optional NFI: Calibrate velocities for dopplergram
;		   modes. NYI.
;               i. Optional NFI: Apply polarization calibration to
;                  magnetogram modes. NYI. 
;		j. Optional: shift/scale to reference frame coordinates
;               k. Optional: extract a subimage from each image
;	        l. Output the corrected image(s) in an updated structure
;		     and data cube, and optionally output a FITS file
;		     with 1 or more images plus a binary extension or as
;		     2D flat FITS files (1 per image)
;
; CATEGORY:
;	FITS processing
;
; SAMPLE CALLING SEQUENCE:
;       FG_PREP, index, data, index_out, data_out
;            [darkdir=darkdir],	[user_dark=user_dark],
;            [user_dindex=user_dindex], [dark_image=dark_image],       
;            [dark_index=dark_index],                                    
;            [flatdir=flatdir], [user_flat=user_flat],         
;            [user_findex=user_findex], [flat_image=flat_image],       
;            [flat_index=flat_index],                                    
;            [/no_shiftpix], [/shiftscale],   
;            [/no_badpix], [/no_darksub],       
;            [/no_flatfield], [/tf_deripple],     
;            [/redspot], [/doppler], [/polarcal],                                        
;            [/despike], [/nofloat],            
;            [x0=x0], [x1=x1], [y0=y0], [y1=y1],                
;            [subimgx=subimgx], [subimgy=subimgy],             
;            [center=center],                                            
;            [nodata=nodata], [no_calib=no_calib],	       
;            [original=original],       
;            [outdir=outdir], [/outflatfits],     
;            [/qstop], [/quiet], [/verbose],	       
;            [run_time=run_time],           
;            [version=progver], [name=prognam]
;
;       1. Simplest use with filenames as input:
;       fg_prep, filename, -1, index_out, image_out
;	fg_prep, file_list(indir,infile), ss(0:8), index_out, data_out  
;		    (with ss=lindgen(total))
;
;       Example for the beginner:
;	fg_prep, '/net/kokuten/archive/prelaunch/sci/2006/01/22/FG/H0130/FG2006D022_013040.502.fits', -1, $
;                 index, image ,/despike ,/outflatfits ,outdir=outdir
;
;       Corrected image is returned in IMAGE. FITS header structure is in INDEX.
;
;       2. Simplest use with INDEX and DATA already read into memory:
;       fg_prep, index, data, index_out, data_out
;
;       3. Recommended useage for large datasets. No data loaded into memory:
;       sot_cat,'29-aug-2007T10:45:00','31-aug-2007T10:00:00',cat,files,/level0, $
;                search_array=['wave=G*band*4305']    
;       ; Results in 3391 G-band filenames.
;       fg_prep,files,-1,index,/outflatfits,outfiletemplate='/data1/SOT/FG/gband/sot_20070829_gband.#####.fits'
;       ;Note that no images are returned from this call, only header structures.
;
; INPUTS:
;		 ** There are 2 methods for calling FG_PREP, same as TRACE_PREP **
;  Case	1. input1 - The input SOT FITS filename(s)
;          input2 - The list number(s) to extract and process. If this is -1, then all files are read. 
;
;  Case	2. input1 - The index structure for each of the input images
;		       (e.g., index from output of read_sot with option /image)
;	   input2 - The input data array (cube)
;
; 
; OUTPUTS (OPTIONAL):
;
;	INDEX_OUT - The updated index structure of the input images.
;
;	IMAGE_OUT - Processed output SOT images (data cube).  
;			Default data type is I*2
;
; OPTIONAL INPUT KEYWORD PARAMETERS:
;	
;	DARKDIR   - Directory for ~monthly dark current processed fits images,
;			default =  "concat_dir('$tdb','tdc_darks')",  which 
;			should work at most sites (including vestige
;			and EOF).
;
;	USER_DARK - If set, use the user-supplied dark current image.
;                       This image is assumed to be sized, binned, 
;			and pixel aligned already.
;
;	FLATDIR   - Directory for ~quarterly flat field processed fits images,
;			default =  "concat_dir('$tdb','tdc_darks')",  which 
;			should work at most sites (including vestige
;			and EOF).
;
;       USER_FLAT  - If set, use  the user-supplied flatfield image.
;                       This image is assumed to be sized, binned, 
;			and pixel aligned already.
;
;       NO_SHIFTPIX- Set this keyword to skip the CCD readout pixel corrections.
;                    [default: perform it]
;
;       NO_BADPIX  - Set this keyword to skip the bad pixel correction.
;                    [default: perform it]
;
;	NO_DARKSUB - Set this keyword to skip dark subtraction.
;                    [default: perform it]
;  
; 	NO_FLATFIELD - Set this keyword to skip flat fielding correction.
;                      [default: perform it]
;
;       NO_POINTING - Set this keyword to suppress correction to
;                      XCEN and YCEN in the image header for
;                      SOT-to-spacecraft pointing offset. The default
;                      is to correct for the offset. Note that
;                      the images are not changed, only the header
;                      values of XCEN and YCEN are altered.
;
;       TF_DERIPPLE - Set this keyword to perform removal of Tunable filter I-ripple.
;
;       REDSPOT - Set this keyword to perform removal of the red spot in BFI 668.4nm images.
;
;       DOPPLER   - Set this keyword to perform the doppler velocity calibration of FGDG images.
;
;       POLARCAL  - Set this keyword to perform the polarization calibration of images.
;
;	DESPIKE   - Set to remove radiation-belt/
;			cosmic-ray spikes. This
;			method uses convolution and thresholding to remove 
;		        spikes and may remove small real features.
;
;	NOFLOAT   - Set if you want to return I*2 images for filtergrams. [Default is Float].
;
;	X0        - Set to lower left x-position of subimage to be extracted.
;       X1        - Set to upper right x-position of subimage to be
;                   extracted. Overrides SUBIMG size, if that is also
;                   set.  
;       Y0        - Set to lower left y-position of subimage to be extracted.
;       Y1        - Set to upper right y-position of subimage to be
;                   extracted. Overrides SUBIMGY size, if that is also
;                   set. 
;	SUBIMGX   - Set to x size of subimage for extraction at LL
;                   corner (X0,Y0).
;
;	SUBIMGY   - Set to y size of subimage for extraction at LL
;                   corner (X0,Y0).
;       CENTER    - Set with SUBIMGX and SUBIMGY specified to do a
;                   centered extraction. Overrides X0 and Y0 settings.
;
;       SHIFTSCALE- Set to shift and scale images to a reference image
;                   position and scale. Default (/SHIFTSCALE) is
;                   G-band reference. SHIFTSCALE='red' shifts and scales to
;                   Red continuum reference. Strings allowed:
;                   'cn','ca','gband','blue','green','red'
;
;       NODATA    - Set this if you don't want output array (auto set if
;                   lt 3 params).
;
;	NO_CALIB  - Set to skip all calibration steps but not various output
;                   options such as reference scaling or sub-frame extraction. 
;
;	ORIGINAL  - Set if input files need more keywords for dark
;                   subtraction.
;
;       OUTDIR    - Destination directory for prepped FITS files.
;
;       OUTFLATFITS - Set if want to output a 2D flat FITS file (1 per
;                     data file).
;
;       OUTFILETEMPLATE- filename template to use in writing out FITS
;                   files. Default filename is "sotYYYYMMDD_hhmmss.[index].fits". 
;                   Use "#" signs in template to designate numerical
;                   series. Example: "sot_20070816_gband.####.fits". 
;                   Output files will thus be named
;                   "sot_20070816_gband.0001.fits", etc. 
;
;       PREFIX    - Set to a string value that will preceed the
;                   default output filename. E.g. 'sot' will result in
;                   filenames such as 'sot20070829_104648.9.fits'
;
;	QSTOP	  - Set to stop in this procedure for debuging
;                   purposes.
;
;	QUIET     - Set for fewer messages, default is loud.
;
;	VERBOSE   - Set for lots of messages and intermediate data listings;
;			suppresses quiet
;
; OPTIONAL OUTPUT KEYWORD PARAMETERS:
;	run_time  - The run time in seconds for FG_PREP
;
; COMMON BLOCKS: none.
;
; RESTRICTIONS:
;
;       Flatfield images do not yet exist for Mg Ib 517.2nm or Fe I
;       525.0 nm and 557.6nm NFI modes. 
;
;       NFI Shutterless Stokes mode is not properly corrected yet. 
;
;       Images are returned as floating point. Set /NOFLOAT to force
;       return of integer data. 
;
; PROCEDURE:
;       Each raw image of the FITS filelist, or of the index and data array, 
;	are processed one at a time. Missing pixels are handled with a TBD algorithm.
;       Then a dark image appropriate to the exposure and summing/binning of the image is
;       subtracted. Negative values are set to 0. This is followed by multiplication
;       by a gain image g = 1./(F - Df) where F-Df is the flatfield image appropriate 
;       for the image wavelength, exposure, and summing/binning, corrected for dark noise.
;       NFI images are corrected for the tunable filter intensity ripple.
;       Optional procedures such as radiation spike removal, image cropping, and 
;       alignment to a reference image or coordinate system can be applied.
;
;	In the future, corrections to the image may also be made for response
;	normalization to physical units.
;
;	The images are returned as an index structure and a datacube, and 
;	optionally as a SOT binary extension FITS file with an updated header,
;	1 or more images per file, and a filename of "TBD", or as
;	a 2D flat FITS file (1 image per file) with a filename of 
;	"TBD", where # = sai integer (0, 1, or 2).
;
;
; MODIFICATION HISTORY:
;v0.1    Outline created 14-Sep-2006 by T. Berger based on trace_prep.pro
;v0.2    Draft release 13-Oct-2006. TEB.
;v0.3    Expanded applicability of nfi_pcalx to include simple filtergrams and doppler.
;        Switched polarization calibration to be before doppler calibration.
;        Changed 'fgmg' switch to 'polcal'. TEB 13-Oct-2006.
;v0.4    Draft release suitable for bad row/col correction and dark subtraction
;        of FG images. Handles 2x2 full FOV images. Not tested on 4x4 summing or
;        2x2 sub-field readout images (because they don't exist in the database).
;v0.5    Replaced fg_bad_pix with fg_shift_pix in correction chain. Fg_bad_pix now
;        after fg_flatfield.
;v0.6    Changed calls to subroutines to include Katsukawa's masks. TEB 25-Apr-2007.
;v1.0    Refined extraction method. Added X1 & Y1 keywords. TEB 26-Apr-2007.
;v1.1    Made redspot removal, tf_deripple, doppler calibration, and polarization 
;        calibration optional instead of anti-optional. Removed /update keyword from
;        shift_pix, flatfield, and bad_pix. TEB 14-May-2007.
;v1.2    Changed despike routine to sot_nospike.pro and added this routine to SSW.
;        TEB 14-May-2007.
;v1.3    Added shift/alignment capability. TEB 10-Oct-2007.
;v1.3.1   Bug fix. case -> switch on parameter check. TEB 24-Oct-2007.
;v1.4    Major modifications to handle multi-image modes. Header
;        updates moved to subroutines. TEB, LMSAL 25-Oct-2007.
;v1.4.1   Extraction code redesigned. TEB LMSAL 26-Oct-2007.
;v1.4.2   Revised required tags: drop exptime to handle
;        shutterless. TEB 2-Nov-2007.
;v1.4.3   Revised input parameter handling and data_out structure. TEB
;        5-Nov-2007. 
;v1.4.4   Default output is now floating point. Option to change to I*2
;        for BFI filtergrams. Fixed extraction bugs. TEB 7-Nov-2007. 
;v1.4.5   Minor bug fixes. YK 8-Nov-2007.
;v1.4.6   Spelling of "TEMPORARY" corrected. TEB 20-Nov-07.
;v1.4.7   Minor bug fix (flat_findex -> user_findex). TEB 8-Jan-08.
;v1.4.8   Minor typo/bug fix on line 476. TEB 30-jan-08.
;v1.4.9   Changed polarcal switch to include gen_id 17. TEB 18-Apr-2008.
;v1.5    Changed the dark sub and flatfield tests to a switch statement. TEB 18-Apr-2008.
;v1.5.1   Minor modification in handling the LHS vignetting of NFI. YK 3-Feb-2009.
;v1.5.2   Added NO_POINTING keyword and Greg Slater's pointing correction call.
;        Also moved the file writing segment into the main loop to
;        allow file output one at a time. Updated sub-frame extraction code
;        also. TEB 5-Mar-2009.
;v1.5.3   Moved call to pointing update outside of main loop.  Determine corrected values
;        for all index records there and store in special reduced tag 'index_pnt'.
;        Then during looping update each index0 with values from
;        index_pnt. G.Slater 16-Mar-2009 LMSAL.
;v1.5.3.1 Removed "data_out =0" for /outflatfits case. Can now write
;        data and get it returned in a variable. T.Berger 23-June-2009 LMSAL.
;v1.5.4   Corrected error in implementation of sot boresight offset correction to
;         XCEN/YCEN tags
; TO DO:
;         - Check for, and flag or delete, bad images (lots of 0s in an
;           image).
;
;         - Move sub-frame extraction to start of loop and only
;           calibrate small sub-frame.
;-
;-----------------------------------------------------------------------------
;;Name, version, and timing information
prognam = 'FG_PREP.PRO'
progver = 'V1.5.3.1'
t0 = SYSTIME(1)
t1 = t0					; Keep track of running time

;Constants
xccdmax = 4095   ;maximum x-value on FG CCD array
yccdmax = 2047   ;maximum y-value
xccd0 = 0        ;default LL corner of image
yccd0 = 0        ;default LL corner of image

;Keyword processing
loud = 1 - KEYWORD_SET(quiet)
verbose = KEYWORD_SET(verbose)
if (verbose eq 1) then loud = 1
if (loud) then PRINT, 'Running ', prognam, ' ', progver
no_shiftpix = KEYWORD_SET(no_shiftpix)
no_badpix = KEYWORD_SET(no_badpix)
no_darksub = KEYWORD_SET(no_darksub)
no_flatfield = KEYWORD_SET(no_flatfield)
tf_deripple = KEYWORD_SET(tf_deripple)
redspot = KEYWORD_SET(redspot)
doppler = KEYWORD_SET(no_doppler)
polarcal = KEYWORD_SET(polarcal)
original = KEYWORD_SET(original)
if KEYWORD_SET(outflatfits) then if (N_ELEMENTS(outdir) le 0) then begin
   MESSAGE,/info,'Warning: no output directory specified with /OUTFLATFITS keyword. Files will be written to current directory.'
   outdir = curdir()
end
if KEYWORD_SET(shiftscale) then begin
   ssflag = 1
   if shiftscale eq 1 then sswave='blue' else sswave=shiftscale
end else ssflag = 0
if not KEYWORD_SET(outdir) then outdir=curdir()
if KEYWORD_SET(outfiletemplate) then begin
  ;delimiter: ASCII 35B = '#'
   dl = 35B
   ts = BYTE(outfiletemplate)
   nel = N_ELEMENTS(ts)
   nps = 0
   for i=0,nel-1 do if ts(i) eq dl then nps = nps + 1
   if nps eq 0 then begin
      MESSAGE,'Invalid filename template: no #-signs'
      STOP ; RETURN
   end
end

;Parse the input arguments
nfiles = 0
read_file = 0
indexin = 0
switch N_PARAMS() of
   0:
   1: begin
      MESSAGE,'At least 2 input parameters are required.'
      STOP ; RETURN
   end
   2: PRINT, 'WARNING: no output header index or data array specified - nothing will be returned.'
   3: PRINT, 'WARNING: no output data array specified - only header index structures will be returned.'
   4: begin
      case datatype(input1) of
         'STR': read_file = 1               ;Case 1: filename or list of filenames.
         'STC': indexin = 1                 ;Case 2: index structure.
         else: begin
                 MESSAGE,'Positional parameter 1: illegal type. Filename string array or SSW Index structure array required. Returning...'
                 STOP ; RETURN
               end
      endcase
     
      ;check input
      nfiles = N_ELEMENTS(input1) ;preliminary number of data elements to process.
      sz2 = SIZE(input2)        ;either -1, a positive integer, an array of positive ints, 
                                ; a 2d image, an array of 2d images,
                                ; a 3d image, or an arry of 3d images.

      d2 = sz2[0]               ;number of dimensions for input2
      ;Case 1: input1 is a filename or STRARR of filenames. 
      if (read_file) then begin      
         case d2 of 
           ; -1 for "read all files in list" or a single file in the list
           0: if (input2 eq -1) then ifile = INDGEN(nfiles) else begin
              nfiles = 1
              ifile = input2 
           end        
           ;list of indices of file list (supplied as INPUT1) to read
           1: if (sz2[1] gt nfiles) then begin    
              MESSAGE,'Number of indices ('+STRTRIM(sz2[1],2)+') is greater than number of filenames in list ('+STRTRIM(nfiles,2)+'). Returning...'
              STOP ; RETURN
           end else begin
              nfiles = sz2[1]
              ifile = input2
           end

           ;Bad input parameter:
           else: begin
              MESSAGE,'Positional parameter 2: illegal type. -1 or list of integers required. Returning...'
              STOP ; RETURN
           end 
        endcase
      end 
      
      ;Case 2: input1 is an SSW Index structure or an array thereof. 
      if (indexin) then begin 
         case d2 of 
            ;a single 2D image
            2: begin
                 nimages = 1
                 ndim = 2
                 if nfiles gt 1 then begin ;single 2D image
                    MESSAGE,'Number of Index structures ('+STRTRIM(nfiles,2)+') exceeds number of images (1). Returning...'
                    STOP ; RETURN
                 end 
              end 
            ;3D array. Could be a list of BFI 2D images 
            ; or it could be a single NFI Stokes set. 
            3: begin
               if nfiles eq 1 then begin  ;a single index structure with a single NFI data element.
                  nimages = 1 
                  ndim = 3
               end else begin
                  nimages = sz2[3] ;array of 2-D images
                  ndim = 2
                  if nimages ne nfiles then begin
                     MESSAGE,'Number of index structures ('+STRTRIM(nfiles,2)+') exceeds number of images ('+STRTRIM(nimages,2)+'). Returning...'
                     STOP ; RETURN
                  end
               end
            end
            ;4D array. A list of 3D NFI Stokes data (e.g. IV image pairs).
            4: begin
               if sz2[4] ne nfiles then begin
                  MESSAGE,'Number of index structures ('+STRTRIM(nfiles,2)+') exceeds number of data files ('+STRTRIM(sz2[4],2)+'). Returning...'
                  STOP ; RETURN
               end else begin
                  nimages = sz2[4]
                  ndim = 3
               end
            end

            else: begin
               MESSAGE,'Positional parameter2: illegal type. 2D image or array of 2D images required. Returning...'
               STOP ; RETURN
            end
         endcase ;d2
      end ;indexin
      BREAK
   end  ;case 4

   else: begin
      MESSAGE,'Useage: fg_prep, input1, input2, [output1, output2], /KEYWORDS...'
      STOP ; RETURN
   end
endswitch

;If requested, and in the case INPUT1/INPUT2 are index/data then create arrays of CRVAL1 and CRVAL2
;  corrected for SC offset:
if KEYWORD_SET(no_pointing) then begin
   if (loud) then MESSAGE,/info,'Pointing corrections not applied.'
end else begin
   if (loud) then box_message,'Establishing pointing correction for all headers in dataset.'
   if (indexin) then sot_offsets = get_shimizu(input1.date_obs)
end

;Read the index structures for sizing the output.
if (read_file) then begin
   if (loud) then begin
      MESSAGE, /INFO, 'Reading first FITS file in list.'
      MESSAGE, /INFO, 'The total number of files to be read = '+STRTRIM(nfiles,2)
   end
   read_sot, input1[ifile[0]], index, data
   sz = SIZE(data)
   ndim = sz[0]                 ;dimensions of each data element in list
   data = 0
end 

if (indexin) then begin
   index = input1
   sz = SIZE(input2[*,*,*,0])  ;subscripting with 0 works for any 2, 3, or 4D array
end
;Nominal size of the images in the data array:
nx0 = sz[1]
ny0 = sz[2]
nx1 = nx0
ny1 = ny0

;Sub-image extraction checks. 
extract = 0
if KEYWORD_SET(x0) then begin
   x0 = LONG(x0)
   if (x0 lt 0) or (x0 gt (nx0-1)) then begin
      MESSAGE,'Subimage X0 parameter < 0 or greater than image size. Returning...'
      STOP ; RETURN
   end
   case 1 of 
      KEYWORD_SET(x1): x1 = LONG(x1)
      KEYWORD_SET(subimgx): x1 = x0 + LONG(subimgx)
      else: x1 = nx0-1
   endcase
   if (x1 gt nx0) then begin
      MESSAGE,/INFO,'X1 parameter larger than image size. Setting to XMAX...'
      x1 = nx0-1
   end
   extract = 1
   nx0 = x1-x0+1
end

if KEYWORD_SET(x1) and not KEYWORD_SET(x0) then begin
   x1 = LONG(x1)
   if (x1 lt 0) or (x1 gt (nx0-1)) then begin
      MESSAGE,'Subimage X1 parameter < 0 or greater than image size. Returning...'
      STOP ; RETURN
   end 
   case 1 of 
      KEYWORD_SET(x0): x0 = LONG(x0)
      KEYWORD_SET(subimgx): x0 = x1 - LONG(subimgx)
      else: x0 = 0
   endcase
   if (x0 lt 0) then begin
      MESSAGE,/INFO,'X0 parameter less than 0. Setting to 0...'
      x0 = 0
   end
   extract = 1
   nx0 = x1-x0+1
end

if KEYWORD_SET(y0) then begin
   y0 = LONG(y0)
   if (y0 lt 0) or (y0 gt (ny0-1)) then begin
      MESSAGE,'Subimage Y0 parameter < 0 or greater than image size. Returning...'
      STOP ; RETURN
   end
   case 1 of 
      KEYWORD_SET(y1): y1 = LONG(y1)
      KEYWORD_SET(subimgy): y1 = y0 + LONG(subimgy)
      else: y1 = ny0-1
   endcase
   if (y1 gt ny0) then begin
      MESSAGE,/INFO,'Y1 parameter larger than image size. Setting to YMAX...'
      y1 = ny0-1
   end
   extract = 1
   ny0 = y1-y0+1
end
if KEYWORD_SET(y1) and not KEYWORD_SET(y0) then begin
   y1 = LONG(y1)
   if (y1 lt 0) or (y1 gt (ny0-1)) then begin
      MESSAGE,'Subimage Y1 parameter < 0 or greater than image size. Returning...'
      STOP ; RETURN
   end 
   case 1 of 
      KEYWORD_SET(y0): y0 = LONG(y0)
      KEYWORD_SET(subimgy): y0 = y1 - LONG(subimgy)
      else: y0 = 0
   endcase
   if (y0 lt 0) then begin
      MESSAGE,/INFO,'Y0 parameter less than 0. Setting to 0...'
      y0 = 0
   end
   extract = 1
   ny0 = y1-y0+1
end

if KEYWORD_SET(subimgx) then begin
   if KEYWORD_SET(x0) and KEYWORD_SET(x1) then begin
      MESSAGE,'SUBIMGX set along with X0 and X1: too many extraction keywords are set. Returning...'
      STOP ; RETURN
   end
   if subimgx gt nx0 then begin
      MESSAGE,'Subimage size in X larger than image size. Returning...'
      STOP ; RETURN
   end
   case 1 of 
      KEYWORD_SET(x0): begin
         x0 = LONG(x0)
         if (x0 lt 0) or (x0 gt nx0) then begin
            MESSAGE,'Subimage X0 parameter < 0 or greater than image size. Returning...'
            STOP ; RETURN
         end
         if (x0 + subimgx) gt nx0 then begin
            MESSAGE,'Subimage size not compatible with specified value of X0. Returning...'
            STOP ; RETURN
         end
      end
      KEYWORD_SET(center): begin
         x0 = (nx0-subimgx)/2
         if x0 lt 0 then begin
            MESSAGE,'Subimage size in X not possible with CENTER keyword set. Returning...'
            STOP ; RETURN
         end
      end
      else: begin
         MESSAGE,/INFO,'Subimage size in X specified without specifying XO. Assuming X0 = 0...'
         x0 =0
      end  
   endcase
   nx0 = subimgx
   x1 = x0 + nx0 - 1
   extract = 1
end

if KEYWORD_SET(subimgy) then begin
   if KEYWORD_SET(y0) and KEYWORD_SET(y1) then begin
      MESSAGE,'SUBIMGY set along with Y0 and Y1: too many extraction keywords are set. Returning...'
      STOP ; RETURN
   end
   if subimgy gt ny0 then begin
      MESSAGE,'Subimage size in Y larger than image size. Returning...'
      STOP ; RETURN
   end
   case 1 of 
      KEYWORD_SET(y0): begin
         y0 = LONG(y0)
         if (y0 lt 0) or (y0 gt ny0) then begin
            MESSAGE,'Subimage Y0 parameter < 0 or greater than image size. Returning...'
            STOP ; RETURN
         end
         if (y0 + subimgy) gt ny0 then begin
            MESSAGE,'Subimage size not compatible with specified value of Y0. Returning...'
            STOP ; RETURN
         end
      end
      KEYWORD_SET(center): begin
         y0 = (ny0-subimgy)/2
         if y0 lt 0 then begin
            MESSAGE,'Subimage size in Y not possible with CENTER keyword set. Returning...'
            STOP ; RETURN
         end
      end
      else: begin
         MESSAGE,/INFO,'Subimage size in Y specified without specifying YO. Assuming Y0 = 0...'
         y0 =0
      end  
   endcase
   ny0 = subimgy
   y1 = y0 + ny0 - 1
   extract = 1
end

;Allocate output data storage
case 1 of 
   KEYWORD_SET(nodata): if (loud) then MESSAGE, /INFO, 'Warning: No output data will be created.'
   KEYWORD_SET(outflatfits): begin
      if (loud) then MESSAGE,/INFO,'Writing output files to '+STRTRIM(outdir,2)
      ;data_out=0
   end
   else: begin  ;allocate output data array
      case ndim of 
         2: begin
            if nfiles gt 1 then data_out = FLTARR(nx0,ny0,nfiles) else data_out = FLTARR(nx0,ny0)
            if (index[0].waveid ge 0) and (index[0].waveid le 6) then begin ;BFI simple FG
               if KEYWORD_SET(nofloat) then data_out = FIX(data_out)
            end
         end
         3: if nfiles gt 1 then data_out = FLTARR(nx0,ny0,sz[3],nfiles) else data_out = FLTARR(nx0,ny0,sz[3])
      endcase
   end
endcase
;release memory
index = 0

if (loud) then MESSAGE,/INFO,'Number of data files to process = '+STRTRIM(nfiles,2)
if (KEYWORD_SET(no_calib) and loud) then MESSAGE,/INFO, 'WARNING: no calibation will be applied to images.'
   
;Main loop.
;Process one file at a time.
ti = SYSTIME(1)
for i=0, nfiles-1  do begin

   file_no = STRTRIM(i,2)    ;was image_no
   if (loud) then begin
      PRINT,'***********************************'
      MESSAGE,/info,'Working on file '+file_no+'....'
      PRINT,'***********************************'
   end

   if (read_file) then begin    ;Read FITS image.
      filename = input1[ifile[i]]
      if (loud) then MESSAGE,/info,'Reading '+filename
      read_sot, filename, index0, data0       
   end 
   if (indexin) then begin
      index0 = input1[i]
      case ndim of 
         2: data0 = input2[*,*,i]
         3: data0 = input2[*,*,*,i]
      end
   end

;Check data file properties
   if (N_ELEMENTS(data0) lt 1) then begin
      MESSAGE,'No images in data file, returning...'           
      STOP ; RETURN
   end
   if not required_tags(index0,'date_obs,waveid,gen_id,obs_id,naxis1,naxis2,camgain,camamp,campsum,camssum,t_fgccd,t_fgceb,xscale,yscale,fgbinx,fgbiny,fgccdix0,fgccdiy0',$
   missing=missing) then begin
      box_message,'The following missing tags are required: '+missing
      STOP ; RETURN
   end
   waveid = gt_tagval(index0,/waveid)
   if (loud) then MESSAGE,/INFO,'Wavelength of data file '+file_no+' = '+STRTRIM(fg_waveid(waveid),2)
   gen_id = gt_tagval(index0,/gen_id)
   obs_id = gt_tagval(index0,/obs_id)
   if (loud) then MESSAGE,/INFO,'Obs_id of data file '+file_no+' = '+STRTRIM(obs_id,2)
   ;number of effective exposures. >1 for shutterless modes.

   if KEYWORD_SET(qstop) then STOP

   nexp = fg_num_eff_exp(index0, quiet=(1-loud))
   image_no = N_ELEMENTS(nexp)  ;the number of images per file

;Now have a valid data file in DATA0 and its index in INDEX0. Proceed with
;appropriate calibration.

; CALIBRATION: 
; All observable types undergo the basic steps of bad pixel removal and flat fielding. Filtergrams and shutterless IQUV
; images undergo dark subtraction. Note that original headers are overwritten in each step. 

   if (KEYWORD_SET(no_calib)) then MESSAGE,/INFO,'WARNING: no_calib keyword set: no calibration steps applied.' $
   else begin

;0.  Determine image type and set flags appropriately
      ; Only way to tell if it's NFI is by wavelength since 'FG simple' occurs for BFI and NFI both.
      if (index0.waveid ge 0) and (index0.waveid le 6) then begin  ;BFI simple FG
         nfi=0 
         polcal = 0
         fgdg = 0
         fgmgdg = 0
      end else begin
         nfi = 1
         ;dopplergram 
         if (gen_id ge 9) and (gen_id le 11) then fgdg = 1 else fgdg = 0
         ;polarization calibration
         if (gen_id ge 9) and (gen_id lt 17) then polcal = 0
         if (not KEYWORD_SET(polarcal)) then begin
            if (loud) then MESSAGE,/INFO,'Warning: no polarization calibration will be applied.'
            polcal = 0
         end else polcal=1
         ;magnetogram/dopplergram mode (TBD)
         fgmgdg=0
      end


;1.  Correct camera readout defects.
;    check .HISTORY (only apply the correction once per image).
      if (no_shiftpix) then begin
         MESSAGE,/INFO, 'No readout corrections applied.'
         index1 = TEMPORARY(index0)
         data1 = TEMPORARY(data0)
      end else begin
         previous = get_history(index0,caller='FG_SHIFT_PIX',found=found)
         if (not found) then begin 
            fg_shift_pix, index0, data0, index1, data1, mask1 $
              , no_vignet=no_flatfield $
              , version=version, verbose=verbose, /update
            index0 = 0
            data0 = 0
         end else begin
            if (loud) then box_message,'Readout corrections already applied. Skippping...'
            index1 = TEMPORARY(index0)
            data1 = TEMPORARY(data0)
         end
      end

;2. Subtract dark pedestal (adc offset) and dark current (small) 
;    check .HISTORY (only apply the correction one time) 
     switch 1 of 
        no_darksub:
        fg_image_type(index1,0) eq 'VI':
        STREGEX(fg_image_type(index1,0),'^D',/bool): begin
           data2 = FLOAT(TEMPORARY(data1))
           index2 = TEMPORARY(index1)
           mask2 = TEMPORARY(mask1)
           if (loud) then MESSAGE,/INFO,'No dark correction applied.' 
           BREAK
        end
        else: begin   
           previous = get_history(index1,caller='FG_DARK_SUB',found=found)
           if (not found) then begin   
              fg_dark_sub, index1, data1, index2, data2, mask2, user_mask=mask1, $
                           user_dark=user_dark, user_dindex=user_dindex, $
                           dark_image=dark_image, dark_index=dark_index, $
                           retrn=retrn, verbose=verbose, /update
              if (retrn ne 0) then begin
                 MESSAGE,'Dark subtraction failed. Returning.'
                 STOP ; RETURN
              end
              data1 = 0
              index1 = 0
           end else begin
              if (loud) then box_message,'Dark correction already applied, skipping...'
              data2 = FLOAT(TEMPORARY(data1))
              index2 = TEMPORARY(index1)
           end
        end
     endswitch

;3. Apply flatfield correction  (unless no_flatfield is set or if images are V/I or DG/I images)
     switch 1 of
        no_flatfield:
        fg_image_type(index2,0) eq 'VI':
        STREGEX(fg_image_type(index2,0),'^D',/bool): begin
           data3 = TEMPORARY(data2)
           index3 = TEMPORARY(index2)
           mask3 = TEMPORARY(mask2)
           if (loud) then MESSAGE,/INFO,'No flatfield correction applied.'
           BREAK
        end
        else: begin
           previous = get_history(index2,caller='FG_FLATFIELD',found=found)
           if (not found ) then begin 
              fg_flatfield, index2, data2, index3, data3, mask3, user_mask=mask2, $
			    flatdir=flatdir, user_file=user_file, $
                            user_flat=user_flat, user_index=user_findex, $
                            flat_image=flat_image, flat_index=flat_index, version=version, $
                            retrn=retrn, loud=loud, verbose=verbose, /update

              if (retrn ne 0) then begin
                 MESSAGE,'Flatfield correction failed. Returning.'
                 STOP ; RETURN
              end
              data2 = 0
              index2 = 0
           end else  begin
              if (loud) then box_message,'FG_FLAT correction already applied, skipping...'
              index3 = TEMPORARY(index2)
              data3 = TEMPORARY(data2)
           end
        end
     endswitch

;4. Handle bad pixels:
      if (no_badpix) then begin
         MESSAGE,/INFO, 'No bad pixel correction applied.'
         index4 = TEMPORARY(index3)
         data4 = TEMPORARY(data3)
         mask4 = TEMPORARY(mask3)
      end else begin
         previous = get_history(index3,caller='FG_BAD_PIX',found=found)
         if (not found) then begin 
            fg_bad_pix, index3, data3, index4, data4, mask4 $
              , user_mask=mask3, nbadpix=nbadpix, no_vignet=no_flatfield $
              ,  version=version, verbose=verbose, /update
            if (loud) then for j=0,image_no-1 do begin
               jbadpix = nbadpix[j]
               if (jbadpix gt 0) then MESSAGE, /INFO,' Corrected '+STRTRIM(jbadpix,2)+' bad pixels in image'+STRTRIM(j,2)
            end
            index3 = 0
            data3 = 0
         end else begin
            if (loud) then box_message,'Pixel corrections already applied. Skippping...'
            index4 = TEMPORARY(index3)
            data4 = TEMPORARY(data3)
         end
      end

;Nominal index and data arrays are INDEX4 and DATA4.

;Correct pointing values. No changes to data, just to index.
      if KEYWORD_SET(no_pointing) then begin    ; do nothing but warn user...
         if (loud) then MESSAGE, /INFO, "NOTE: SOT Pointing update not done ('NO_POINTING' keyword set)"
      end else begin

         previous = get_history(index4, caller='sot_fix_pointing', found=found)
         if (not found) then begin
            index4a = sot_fix_pointing(index4, use_shimizu=use_shimizu, $
                                       corr_db_good=corr_db_good, n_corr_db_good=n_corr_db_good, $
                                       topdir_genx=topdir_genx, loud=loud, verbose=verbose, _extra=extra)
if size(index4a, /type) eq 8 then begin
            index4 = TEMPORARY(index4a)
            history_rec = 'pointing updated'
            update_history, index4, history_rec, caller='sot_fix_pointing', version=version
            if (loud) then MESSAGE, /INFO,' HISTORY record updated for data file '+file_no+':  ' + history_rec
endif else if loud then box_message, 'Could not update pointing using SOT SDO correlation database.  No appropriate matches found.'
         end else if (loud) then box_message,'pointing update already applied. Skipping...'

      end

;Correct roll and plate scale values. No changes to data, just to index.
      if not KEYWORD_SET(do_roll_scale) then begin    ; do nothing but warn user...
         if (loud) then MESSAGE, /INFO, "NOTE: SOT roll/scale update not done ('NO_ROLL_SCALE' keyword set)"
      end else begin

         previous = get_history(index4, caller='sot_fix_roll_scale', found=found)
         if (not found) then begin
            index4a = sot_fix_roll_scale(index4, interp=interp, _extra=extra)
            index4 = TEMPORARY(index4a)
;           history_rec = 'roll and plate scale updated'
;           update_history, index4, history_rec, caller='sot_fix_roll_scale', version=version
;           if (loud) then MESSAGE, /INFO,' HISTORY record updated for data file '+file_no+':  '+history_rec
         end else if (loud) then box_message,'roll and plate scale update already applied. Skipping...'

      end

;Optional: Radiation despike if /despike is set:
     if KEYWORD_SET(despike) then begin
        sot_nospike, index4, data4, index5, data5, nspikes=nspikes,/update
        index4 = TEMPORARY(index5)
        data4 = TEMPORARY(data5)
     end 

;BFI 6684 red spot removal.
     if (index4.waveid eq 6) then $
        if (redspot) then begin
           previous = get_history(index4, caller='FG_REDSPOT_REMOVE', found=found)
           if (not found) then begin
              fg_redspot_remove, index4, data4, index5, data5, _extra=extra
              index4 = TEMPORARY(index5)
              data4 = TEMPORARY(data5)
              tagval5 = 'Image corrected for 668.4nm red spot'
              update_history, index4, tagval5, caller='FG_REDSPOT_REMOVE', version=version
              if (loud) then MESSAGE, /INFO,' HISTORY record updated for data file '+file_no+':  '+tagval5
           end else if (loud) then box_message,'BFI 668.4nm red spot correction already applied. Skipping...'
        end else if (loud) then MESSAGE,/INFO,'Red spot removal not performed.'

;NFI-only steps:
     if (nfi) then begin

        if (tf_deripple) then $
        begin
           previous = get_history(index4, caller='FG_TF_DERIPPLE', found=found)
           if (not found) then begin
              fg_tf_deripple, index4, data4, index5, data5, _extra=extra
              index4 = TEMPORARY(index5)
              data4 = TEMPORARY(data5)
              tagval4 = 'Image corrected for TF I-ripple'
              update_history,index4, tagval4, caller='FG_TF_DERIPPLE', version=version
              if (loud) then MESSAGE, /INFO,' HISTORY record updated for data file '+file_no+':  '+tagval4
           end else if (loud) then MESSAGE,/INFO,'TF Ripple correction already applied. Skipping...'
        end else if (loud) then MESSAGE,/INFO,'No TF I-ripple correction applied.'
     
;Apple polarization calibration. If simple filtergram, just return the Q,U,V->I crosstalk.
        if (polarcal) then begin
           previous = get_history(index4, caller='FG_POLAR_CAL', found=found)
           if (not found) then begin 
              fg_polar_cal, index4, data4, index5, data5, id_table=id_table, calver=calver, version=version
              index4 = TEMPORARY(index5)
              data4 = TEMPORARY(data5)
              tagval8 = STRING('Polarization calibration applied. Version = ',version, format='(a,a)')
              update_history,index4,tagval8, caller='FG_POLAR_CAL', version=version
              if (loud) then MESSAGE, /INFO,'HISTORY record updated for data file '+file_no,':  '+tagval8
           end else if (loud) then STOP, 'Polarization calibration already applied. Skipping...'
        end

;If dopplergram, apply velocity calibration:
        if (fgdg) then $
           if (doppler) then begin 
              previous = get_history(index4, caller='FG_DG_CAL', found=found)
              if (not found) then begin 
                 fg_dg_cal,index4, data4, index5, data5, vtable=vtable, version=version
                 index4 = TEMPORARY(index5)
                 data4 = TEMPORARY(data5)
                 tagval7 = STRING('Doppler velocity calibrated. Vtable = ',vtable,'. Version = ',version, format='(a,a,a,a)')
                 update_history,index4, tagval7, caller='FG-DG_CAL',version=version
                 if (loud) then MESSAGE, /INFO,'HISTORY record updated for data file '+file_no+':  '+tagval7
              end else if (loud) then MESSAGE,/INFO,'Doppler velocity calibration alread applied. Skipping...'
           end else if (loud) then MESSAGE, /INFO, 'No doppler velocity calibration applied.'
      
;If doppler/magnetogram, apply different calibration. Note that Polarization cal has been applied above.
        if (fgmgdg) then $
           if (doppler) then begin 
              previous = get_history(index4, caller='FG_MGDG_CAL', found=found)
              if (not found) then begin 
                 fg_mgdg_cal,index4, data4, index5, data5, vtable=vtable, version=version
                 index4 = TEMPORARY(index5)
                 data4 = TEMPORARY(data5)
                 tagval7 = STRING('Magnetogram/Doppler velocity calibrated. Vtable = ',vtable,'. Version = ',version, format='(a,a,a,a)')
                 update_history,index4, tagval7, caller='FG_MGDG_CAL',version=version
                 if (loud) then MESSAGE, /INFO,'HISTORY record updated for data file '+file_no+':  '+tagval7
              end else if (loud) then MESSAGE,/INFO,'Magnetogram/Doppler velocity calibration alread applied. Skipping...'
           end else  MESSAGE, /INFO, 'No magnetogram/doppler velocity calibration applied.'

     
     end ;of NFI-only steps.

   end                            ; end of all calibration steps.

;Optional: scale/align with reference image:
   if (ssflag) then begin
      ssstring = STRTRIM(sswave,2)
      if (loud) then MESSAGE,/INFO, 'Shifting and scaling image to reference image: '+ssstring
      previous = get_history(index4,caller='FG_REG_WAVE',found=found)
      if (not found) then begin 
         if (loud) then fg_reg_wave, index4, data4, index5, data5, WAVE_REF=sswave else $
            fg_reg_wave, index4, data4, index5, data5, WAVE_REF=sswave, /QUIET
         index4 = TEMPORARY(index5)
         data4 = TEMPORARY(data5)
         tagval6 = STRING('Aligned and scaled to reference image:'+ssstring,format='(A)')
         update_history, index4, tagval6, caller='FG_REG_WAVE',version=version
         if (loud) then MESSAGE,/INFO, 'HISTORY record updated for data file'+file_no+':  '+tagval6
      end else if (loud) then box_message,'Image already scaled and aligned: '+previous
   end

;Optional: Extract a subimage if /extract is set to new x axis:
;     note: .HISTORY is checked and updated within function 
   if (extract) then begin
      fg_extract, index4, data4, index5, data5, x0=x0, y0=y0, x1=x1, y1=y1, $
                  nx1=subimgx, ny1=subimgy, $
                  version=version, verbose=verbose
      nx  = subimgx
      ny  = subimgy
      index4 = TEMPORARY(index5)
      data4 = TEMPORARY(data5)
      tagval5a = STRING('Subimage extracted with lower left corner at ',  $
                        x0,y0,' for size ', nx,' x ', ny, format='(a,i5,",",i5,a,i5,a,i5)')
      update_history,index4,tagval5a, caller='FG_EXTRACT',version=version
      if (loud) then MESSAGE, /INFO, 'HISTORY record updated for data file'+file_no+':  '+tagval5a
   end

; DATA4 is floating point or double precision at this point. For BFI
; filtergrams only, transform to I*2 if /NOFLOAT is set:
   if not (nfi) and KEYWORD_SET(nofloat) then data4 = FIX(ROUND(data4))

;Load output arrays if needed:
; Concatenate the headers

   if N_PARAMS() gt 2 then $
   if (i eq 0) then index_out = index4 else begin
      out_str = index_out[0]
      index_out = [str_copy_tags(out_str,index_out), str_copy_tags(out_str,index4)]
   end

   if N_PARAMS() gt 3 then $
      case ndim of 
      2: data_out[*,*,i] = data4
      3: data_out[*,*,*,i] = data4
   end

;done calibrating one image   
   systm = SYSTIME(1)
   loop_time = systm - ti
   if (loud) then MESSAGE, /INFO, STRING('****************Single data file took ', loop_time, ' seconds',  $
                                         format='(a,f5.1,a)')
   ti = systm

;Write out images one at a time:
   if (KEYWORD_SET(outflatfits)) then begin
      if (KEYWORD_SET(outfiletemplate)) then begin
         snum = STRTRIM(STRING(i),2)
         while STRLEN(snum) lt nps do snum = '0' + snum
         temp2 = outfiletemplate
         STRPUT,temp2,snum,STRPOS(temp2,'#')
         if (loud) then write_sot,index4,data4,outfile=temp2,/flat_fits else $
            write_sot,index4,data4,outfile=temp2,/flat_fits,/quiet
      end else begin
         indir = file_break(filename, /path)
         infil = file_break(filename)
         outfil = infil
         ss_und = strpos(filename,'_')
         prefix = strmid(filename,0,ss_und)
;        outdir = str_replace(indir, '/oberon/archive1/hinode/sot/level0', '~/data/sot/level1')
         outdir = str_replace(indir, 'level0', 'level1')
         mk_dir, outdir
         outfile = concat_dir(outdir, outfil)
;        if (loud) then write_sot,index4,data4,outdir=outdir,prefix=prefix,/flat_fits else $
;           write_sot,index4,data4,outdir=outdir,prefix=prefix,/flat_fits,/quiet
;;       if (loud) then write_sot,index4,data4,outfile=outfile,/flat_fits else $
;;          write_sot,index4,data4,outfile=outfile,/flat_fits,/quiet
         write_sot,index4,data4,outfile=outfile,/flat_fits
      end
   endif

; Write a FITS binary file
   if (KEYWORD_set(fitsbinary)) then begin
      if (loud) then write_sot,index4,data4,outdir=outdir,prefix='sot' else $
         write_sot,index4,data4,outdir=outdir,prefix='sot',/quiet
   endif 

end  ;Main FOR loop end.

runtime = SYSTIME(1)-t0
if (KEYWORD_SET(run_time)) then run_time=runtime
if (loud) then MESSAGE, /info, string('Total processing time of ', runtime,' seconds for ', $  
         nfiles, ' files' ,format='(a,f8.1,a,i5,a)')

if (KEYWORD_SET(qstop)) then STOP

RETURN
END


 
