; This code makes a collapsed line map from an OSIRIS
; data cube.
;
;For example,
;osredx_makeim,'q0100-bx172-o3.fits',gfwhm=1.5,zstart=319,zstop=330

pro osredx_makeim,file,gfwhm=gfwhm,lfwhm=lfwhm,mincov=mincov,zstart=zstart,zstop=zstop,outdir=outdir

a=mrdfits(file,0,hdr0)
covmap=mrdfits(file,3)
b=transpose(a,[2,1,0]) 
b_cov=transpose(covmap,[2,1,0])
maxcov=max(b_cov)
nz=fxpar(hdr0,'NAXIS1')

; Output subdirectory
if (~keyword_set(outdir)) then outdir='./images/'

; Default parameters
; Minimum coverage threshhold
if (~keyword_set(mincov)) then mincov=0.7*maxcov
; Value (in spaxels) to spatially smooth by
if (~keyword_set(gfwhm)) then gfwhm=1.5
; Start/stop range for the collapse of cube
if (~keyword_set(zstart)) then zstart=0
if (~keyword_set(zstop)) then zstop=nz-1

; Create output directory
if file_test(outdir,/directory) eq 0 then spawn, '\mkdir -p '+outdir

if (keyword_set(lfwhm)) then begin
   new=b
   for i=0,(size(b))[1]-1 do begin
      for j=0,(size(b))[2]-1 do begin
        new[i,j,*]=smooth(b[i,j,*],lfwhm)
      endfor
   endfor
   b=new
endif
   

c=b
for i=0,nz-1 do begin
   newim=filter_image(b[*,*,i],fwhm_gaussian=gfwhm)
   temp=b_cov[*,*,i]
   indx=where(temp lt mincov)
   newim[indx]=0.
   c[*,*,i]=newim
endfor

if (zstart eq zstop) then im=c[*,*,zstart]
if (zstart ne zstop) then begin
  test=c[*,*,zstart:zstop]
  im=median(test,dimension=3)
endif
d=transpose(c,[2,1,0])
; Write out the smoothed/masked data cube
;outname=ml_strreplace(file,'.fits','_smooth.fits')
;writefits,outname,d,hdr0

; Make an image header and copy over WCS
mkhdr,hdrim,im
fxaddpar,hdrim,'CTYPE1',fxpar(hdr0,'CTYPE3')
fxaddpar,hdrim,'CTYPE2',fxpar(hdr0,'CTYPE2')
fxaddpar,hdrim,'CUNIT1',fxpar(hdr0,'CUNIT3')
fxaddpar,hdrim,'CUNIT2',fxpar(hdr0,'CUNIT2')
fxaddpar,hdrim,'CRVAL1',fxpar(hdr0,'CRVAL3')
fxaddpar,hdrim,'CRVAL2',fxpar(hdr0,'CRVAL2')
fxaddpar,hdrim,'CRPIX1',fxpar(hdr0,'CRPIX3')
fxaddpar,hdrim,'CRPIX2',fxpar(hdr0,'CRPIX2')
fxaddpar,hdrim,'CDELT1',fxpar(hdr0,'CDELT3')
fxaddpar,hdrim,'CDELT2',fxpar(hdr0,'CDELT2')
fxaddpar,hdrim,'PC1_1',fxpar(hdr0,'PC3_3')
fxaddpar,hdrim,'PC1_2',fxpar(hdr0,'PC3_2')
fxaddpar,hdrim,'PC2_2',fxpar(hdr0,'PC2_2')
fxaddpar,hdrim,'PC2_1',fxpar(hdr0,'PC2_3')

; Write out the image
outname=ml_strreplace(file,'.fits','_image.fits')
writefits,concat_dir(outdir,outname),im,hdrim

return
end
