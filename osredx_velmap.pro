; This is the code to make velocity maps and do kinematic fitting.
; Mostly it wraps MPFIT
; Example: osredx_velmap,'q0142-bx195-o3.fits',2.3807,5008.24,'o3double',gfwhm=2.5,minsig=30.




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Two-component Gaussian function tying the
; two [O III] lines
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function o3_twocomponent, wavevec, p

; Vacuum wavelengths of the O III lines
waveo3_1=4960.30
waveo3_2=5008.24

; p[0]=strong o3 wave, p[1]=width, p[2]=strong o3 strength, p[3]=constant
wave1=p[0]*waveo3_1/waveo3_2
wave2=p[0]
sig1=p[1]
sig2=p[1]
flux1=p[2]/2.985; Ratio from Osterbrock Table 3.8
flux2=p[2]
; Our gaussian function actually needs amplitudes though;
; where flux=amplitude*sigma*sqrt(2pi)
amp1=flux1/(sig1*sqrt(2*!PI))
amp2=flux2/(sig2*sqrt(2*!PI))

; Construct using the 'gaussian' function
temp1=gaussian(wavevec,[amp1,wave1,sig1])
temp2=gaussian(wavevec,[amp2,wave2,sig2])

model=temp1+temp2+p[3]

return,model
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Generic single-component Gaussian
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function onecomponent, wavevec, p

; p[0]=wavelength, p[1]=width, p[2]=strength, p[3]=constant
wave1=p[0]
sig1=p[1]

flux1=p[2]
; Our gaussian function actually needs amplitudes though;
; where flux=amplitude*sigma*sqrt(2pi)
amp1=flux1/(sig1*sqrt(2*!PI))

; Construct using the 'gaussian' function
temp1=gaussian(wavevec,[amp1,wave1,sig1])

model=temp1+p[3]

return,model
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Farm fitting of sky spectra to get the LSF
; to a subroutine.  Return spectral resolution
; R=lambda/FWHM
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function fit_skyspec, wavevec, skyspec
  ; List of fairly isolated skylines with which to measure the LSF
  ohlines=[15332.402,15597.631,16030.831,16840.481,16903.679,17008.757,17123.659,20412.680]

  ; Set up basic boundaries for the skyline fits
  ; 4 parameters; wavelength, sigma, flux, baseline
  parinfo = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0]}, 4)
  parinfo[*].value = [0., 3., 0.01, 0.]; starting values
  ; Initial guess
  parms0=[0., 3., 0.1, 0.]
  ;Bound the velocities
  parinfo[0].limited[0] = 1
  parinfo[0].limited[1] = 1
  ; Bound the sigma
  parinfo[1].limited[0] = 1
  parinfo[1].limited[1] = 1
  parinfo[1].limits[0]  = 0.
  parinfo[1].limits[1]  = 150.
  ; Require positive flux
  parinfo[2].limited[0] = 1
  parinfo[2].limits[0]  = 0.D
  mockerr=replicate(1e-3,n_elements(wavevec))

  ; Loop over good skylines in the wavelength range
  goodsky=where((ohlines gt min(wavevec+100.))and(ohlines lt max(wavevec-100.)),ngood)
  allR=fltarr(ngood)
  for i=0,ngood-1 do begin
    parinfo[*].value = [ohlines[goodsky[i]], 3., 0.01, 0.]; starting values
    parms0=[ohlines[goodsky[i]], 3., 0.1, 0.]
    parinfo[0].limits[0]  = ohlines[goodsky[i]]-10.
    parinfo[0].limits[1]  = ohlines[goodsky[i]]+10.

    dlam=where(abs(ohlines[goodsky[i]]-wavevec) gt 10.)
    thissky=skyspec
    thissky[dlam]=0.

    parms=mpfitfun('onecomponent',wavevec,thissky,mockerr,parms0,parinfo=parinfo,yfit=yfit,/quiet)
    allR[i]=ohlines[goodsky[i]]/(parms[1]*2.35); Spectral Resolution=wave/fwhm
  endfor

  ; Use the median to clip out any bad skylines that were blended
  finalR=median(allR)

return,finalR
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Wrapper code to the analysis
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro osredx_velmap,file,zsys,lrest,ltype,gfwhm=gfwhm,mincov=mincov,tag=tag,outdir=outdir,maxvel=maxvel,maxsig=maxsig,minsig=minsig,verbose=verbose,skyfile=skyfile,sncut=sncut

; zrest is systemic z
; lrest is rest wavelength of line in Angstroms
; ltype is line type; 'single' or 'o3double'

device,decomposed=0
loadct,39

flux=mrdfits(file,0,hdr0)
err=mrdfits(file,1,hdr1)
covmap=mrdfits(file,3,hdr3)
flux=transpose(flux,[2,1,0])
err=transpose(err,[2,1,0])
covmap=transpose(covmap,[2,1,0])

; Optional skycube
if (keyword_set(skyfile)) then begin
  skycube=mrdfits(skyfile,0)
  skycube=transpose(skycube,[2,1,0])
endif

; Coverage threshhold
maxcov=max(covmap)


; Default parameters
; Minimum coverage threshhold
if (~keyword_set(mincov)) then mincov=0.7*maxcov
; Value (in spaxels) to spatially smooth by
if (~keyword_set(gfwhm)) then gfwhm=1.5
; Descriptive tag for output filenames
if (~keyword_set(tag)) then tag='image'
; Output subdirectory
if (~keyword_set(outdir)) then outdir='./images/'

; Default fitting parameters
if (~keyword_set(maxvel)) then maxvel=1000. ; +- 1000 km/s
if (~keyword_set(maxsig)) then maxsig=1000. ; 1000 km/s
if (~keyword_set(minsig)) then minsig=-1 ; <0 means use the LSF as a lower bound
if (~keyword_set(maxsig)) then maxsig=1000. ; 1000 km/s

; QA parameters
if (~keyword_set(sncut)) then sncut=6

; Create output directory
if file_test(outdir,/directory) eq 0 then spawn, '\mkdir -p '+outdir



nx=fxpar(hdr0,'NAXIS3')
ny=fxpar(hdr0,'NAXIS2')
nz=fxpar(hdr0,'NAXIS1')


; Get the wavelength vector (vacuum)
wavevec=fxpar(hdr0,'CRVAL1')+findgen(nz)*fxpar(hdr0,'CDELT1')
; Heliocentric correction using header information
; Rough since headers don't account for midpoint time of ALL
; exposures, but good enough since dominated by seasonal variation
; instead of few-hour variation
ra=fxpar(hdr0,'RA')
dec=fxpar(hdr0,'DEC')
jd=fxpar(hdr0,'MJD-OBS')+2400000.5d
long=360.-155.4747d
lat=19.8260
alt=4205.
helio=heliocentric(ra,dec,jd=jd,longitude=long,latitude=lat,altitude=alt)
; Correct wave vector and convert to Angstroms
wavevec = wavevec / (1 + helio/299792.458) * 10

; Average the error vector over entire field to get skylines
errspec=fltarr(nz)
; Error spectrum seems to have weird normalization, apply an empirical
; factor such that data-model spectrum divided by the noise vector is 
; gaussian with sigma=1
; This is only used for accurate weights in the model fitting
fudge=15.
for i=0,nz-1 do begin
  errspec[i]=median(err[*,*,i])/fudge
endfor


; Spatial smoothing and zero out stuff where the coverage is too low
newflux=flux
covmask=median(covmap,dimension=3)
goodindx=where(covmask ge mincov,complement=badindx,ngood)
; Convert to x,y good locations
ind=array_indices(covmask,goodindx)
goodx=ind[0,*]
goody=ind[1,*]
print,'ngood=',ngood

for i=0,nz-1 do begin
  newim=filter_image(flux[*,*,i],fwhm_gaussian=gfwhm)
  ; zero out bad index locations
  newim[badindx]=0.
  newflux[*,*,i]=newim

  ; Smooth the sky cube if provided
  if (keyword_set(skyfile)) then begin
    skyim=filter_image(skycube[*,*,i],fwhm_gaussian=gfwhm)
    skycube[*,*,i]=skyim
  endif
endfor
; write out smoothed cube
temp=transpose(newflux,[2,1,0])
writefits,concat_dir(outdir,'cube_smoothed.fits'),temp,hdr0




; 4 parameters; wavelength, sigma, flux, baseline
parinfo = replicate({value:0.D, fixed:0, limited:[0,0], $
                       limits:[0.D,0]}, 4)

c=3e5; speed of light in km/s
lambda_guess=(1.+zsys)*lrest
lambda_ul=(1.d + maxvel/c)*lambda_guess; Upper limit on wavelength in Angstrom
lambda_ll=(1.d - maxvel/c)*lambda_guess; Lower limit on wavelength in Angstrom
sigma_ll=(minsig/c)*lambda_guess; Lower limit on sigma in Angstrom
sigma_ul=(maxsig/c)*lambda_guess; Upper limit on sigma in Angstrom

; Parameters are [wavelength, sigma, flux, baseline]
parinfo[*].value = [lambda_guess, 5., 0.01, 0.]; starting values
; Initial guess
parms0=[lambda_guess, 5., 0.1, 0.]
;Bound the velocities
parinfo[0].limited[0] = 1
parinfo[0].limited[1] = 1
parinfo[0].limits[0]  = lambda_ll
parinfo[0].limits[1]  = lambda_ul
; Bound the sigma
parinfo[1].limited[0] = 1
parinfo[1].limited[1] = 1
parinfo[1].limits[0]  = sigma_ll
parinfo[1].limits[1]  = sigma_ul
; Require positive flux
parinfo[2].limited[0] = 1
parinfo[2].limits[0]  = 0.D

velmap=replicate(-999.,nx,ny)
sigmap=replicate(-999.,nx,ny)
fluxmap=fltarr(nx,ny)
snmap=fltarr(nx,ny)

; Create a map of the LSF (sigma in km/s).  By default it is constant
; and set to typical OSIRIS properties
lsfmap=replicate(37.5,nx,ny); Sigma in km/s

for i=0,ngood-1 do begin
  ;i=471
  spectrum=newflux[goodx[i],goody[i],*]

  ; Optional sky spectrum fitting to get the LSF
  if (keyword_set(skyfile)) then begin
    skyspec=skycube[goodx[i],goody[i],*]
    skyR=fit_skyspec(wavevec, skyspec); spectral resolution
    thislsf=c/skyR/2.35
    lsfmap[goodx[i],goody[i]]=thislsf; lsf map of sigma in km/s
  endif

  ; If minsig parameter was less than zero, update
  ; the bounded range for the fitter so that sigma cannot
  ; be below instrumental resolution in Angstroms
  thisR=c/lsfmap[goodx[i],goody[i]]/2.35
  if (minsig lt 0) then $
    parinfo[1].limits[0]=lambda_guess/skyR/2.35

  parms=mpfitfun('o3_twocomponent',wavevec,spectrum,errspec,parms0,parinfo=parinfo,yfit=yfit,perror=perr,/quiet)
  diff=spectrum-yfit
  rms=sqrt((moment(diff/errspec))[1])
  if (keyword_set(verbose)) then $
    print,i,' Noise scaling rms: ',rms
  plot,wavevec,spectrum
  oplot,wavevec,yfit,color=250

  thevel=(parms[0]-lambda_guess)/lambda_guess*3e5
  thesig=(parms[1])/parms[0]*3e5
  theflux=(parms[2])
  thesnr=parms[2]/perr[2]

  ; Subtract instrumental line width in quadrature
  thislsf=thislsf*lambda_guess/parms[0]; Correct for difference in guess vs real line center
  if (keyword_set(verbose)) then begin
    print,'Fitted snr: ',thesnr
    print,'Fitted sigma: ',thesig,' km/s'
    print,'Instrumental sigma: ',thislsf,' km/s'
  endif
  thesig=sqrt(thesig*thesig-thislsf*thislsf)
  if (~finite(thesig)) then thesig=0.
  if (keyword_set(verbose)) then $
    print,'Corrected sigma: ',thesig,' km/s'

  ; Remove any bad fits with zero or negative flux
  if (theflux eq 0.) then begin
    thevel=-999.
    thesig=-999.
    theflux=0.
    thesnr=0.
  endif

  ; Remove bad fits with below threshhold snr
  if ((thesnr lt sncut)or(~(finite(thesnr)))) then begin
    thevel=-999.
    thesig=-999.
    theflux=0.
    thesnr=0.
  endif 



  velmap[goodx[i],goody[i]]=thevel
  sigmap[goodx[i],goody[i]]=thesig
  fluxmap[goodx[i],goody[i]]=theflux
  snmap[goodx[i],goody[i]]=thesnr

  if ((keyword_set(verbose))and(snmap[goodx[i],goody[i]] gt sncut)) then wait, 1.0
endfor

; Make headers for the new images
mkhdr,hdrim,fluxmap
; Add WCS
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


writefits,concat_dir(outdir,'fluxmap.fits'),fluxmap,hdrim
writefits,concat_dir(outdir,'velmap.fits'),velmap,hdrim
writefits,concat_dir(outdir,'sigmap.fits'),sigmap,hdrim
writefits,concat_dir(outdir,'snmap.fits'),snmap,hdrim
if (keyword_set(skyfile)) then $
  writefits,concat_dir(outdir,'lsfmap.fits'),lsfmap,hdrim
stop


return
end
