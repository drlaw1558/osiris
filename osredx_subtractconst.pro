pro osredx_subtractconst, zaphot=zaphot, MASK=mask, maskbox=maskbox
;
; Program loops through a specified set of images which have had
; first-pass sky subtraction and subtracts off the median value at
; each spectral channel.  Uses a bad pixel map made from the first
; channel of the first image.
;
; Version 2 written 5/22/07
;
; David R. Law <drlaw@astro.caltech.edu>
; Written 10/11/2006
; Last modified: 8/19/08
;
; Loosely based on code by Randy Campbell
; rdc july 2006
;
; Specifying /zaphot will zap hot values in the final cube.
; This zeroes out values greater than zapval

zapval = 0.15

; Edited 8/19/08 to add MASK keyword
; This is a same-dimension image that lets you specify
; regions to exlude when determining the median at each
; wavelength.  It will still keep them in the final image though.

; get list of files and loop on list
;
file_name = dialog_pickfile( title='Read Files to Process', $
                                     filter='*0.fits', $
                                     get_path=new_path, $
                                     /MUST_EXIST        , $
                                     /MULTIPLE_FILES)

; Read in first file to define constants
dat = readfits(file_name[0],head,/SILENT)
dat=transpose(dat,[2,1,0])

; If no mask specified, include all pixels
if NOT(keyword_set(MASK)) then begin
  mask=fltarr((size(dat))[1],(size(dat))[2])
endif

nspec = sxpar(head,'NAXIS1');Number of spectral pixels
nx = sxpar(head,'NAXIS3');Number of pixels along X axis
ny = sxpar(head,'NAXIS2');Number of pixels along Y axis
nframes = (size(file_name))[1];Number of frames
filter = sxpar(head,'SFILTER');Filter name
filter = STRCOMPRESS(filter,/REMOVE_ALL); Remove whitespace from filter name
sscale = float(sxpar(head,'SSCALE')); Lenslet scale

; Define a bad pixel map
; Anything that has a zero value in 1st spectral slice is called a bad
; pixel (this works with old IFU data with bad channel too).
; Values of 0 are good, values of 1 are bad.
badpix=fltarr(nx,ny)
temp=dat[*,*,0]
badpix[where(temp eq 0)] = 1
; Modify mask file to mask bad pixels out of the median calculation too
mask[where(badpix eq 1)]=1

; Kludge out a certain area
if (keyword_set(maskbox)) then begin
  mask[maskbox[0]:maskbox[1],maskbox[2]:maskbox[3]]=1
  print,'kludged mask!'
endif

for i = 0, nframes-1 do begin
  dat = readfits(file_name[i],head,/SILENT)
  noise = readfits(file_name[i], exten_no = 1,/SILENT)
  qual = readfits(file_name[i], exten_no = 2,/SILENT)

  dat=transpose(dat,[2,1,0])
; find median value of image at each spectral slice and
; subtract it off
; need to leave already zero values alone
; this will only work for cases where both positive and negative
; image of source fall in region used for median, and when
; source is faint compared to the background
  medval=1.

  for k = 0, (nspec-1) do begin
    temp=dat[*,*,k]
    ; Calculate the median
    medval=median(temp[where(mask eq 0)])
    ; And subtract the median from these good pixel values
;    print,'WARNING: MODDED!'
;    medval=meanclipdrl(temp[where(badpix eq 0)])
    temp[where(badpix eq 0)] = temp[where(badpix eq 0)] - medval

;    if (k eq 334) then print,medval

    ; Zero out any really bad pixels
    temp[where(badpix eq 1)] = 0.
    dat[*,*,k]=temp
  endfor

  ; Zap hot pixels??
  if keyword_set(zaphot) then begin
    for l=0,nx-1 do begin
      for m=0,ny-1 do begin
        for k=0,nspec-1 do begin
          if abs(dat[l,m,k]) gt zapval then dat[l,m,k]=0.
        endfor
      endfor
    endfor
  endif

  dat=transpose(dat,[2,1,0])

  t_filename = file_name[i]+'_rs.fits'
;  t_filename = dialog_pickfile( title='Enter name of file to write', $
;                                     filter='*.fits', $
;                                     /write )
  if t_filename ne "" then begin
    writefits, t_filename, dat, head
    writefits, t_filename, noise, /append
    writefits, t_filename, byte(qual), /append
  endif
endfor

return
end
