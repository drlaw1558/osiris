;+
; NAME:
;   osredx_velmap
;
; PURPOSE:
;   This is the code to make velocity maps and do kinematic fitting.
;   Mostly it wraps MPFIT.
;
; CALLING SEQUENCE:
;   
;
; INPUTS:
;   
;
; OPTIONAL INPUTS:
;   
;
; PROCEDURES CALLED:
;
; EXAMPLE:
;   osredx_velmap,'q0142-bx195-o3.fits',2.3807,5008.24,'o3double',gfwhm=2.5,minsig=-1,skyfile='./q0142-bx195-o3sky.fits',sncut=6.,outdir='velfit_o3'
;
; REVISION HISTORY: 
;   29-May-2017 Written (D. Law; dlaw@stsci.edu)
;-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Two-component Gaussian function tying the
; two [O III] lines
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function o3_twocomponent, wavevec, p

; Vacuum wavelengths of the O III lines
waveo3_1=4960.30
waveo3_2=5008.24

; p[0]=strong o3 wave, p[1]=width, p[2]=strong o3 amplitude, p[3]=constant
wave1=p[0]*waveo3_1/waveo3_2
wave2=p[0]
sig1=p[1]
sig2=p[1]
amp1=p[2]/2.985; Ratio from Osterbrock Table 3.8
amp2=p[2]

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

; p[0]=wavelength, p[1]=width, p[2]=amplitude, p[3]=constant
wave1=p[0]
sig1=p[1]

amp1=p[2]

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
  ohlines=[15332.402,15597.631,15972.596,$
           16030.831,16079.753,16128.608,16502.365,16840.481,16903.679,$
           17008.757,17123.659,$
           17449.967,17880.298,17993.962,$
           20412.680,$
           21802.312,21873.518,21955.637]

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
  goodsky=where((ohlines gt min(wavevec+50.))and(ohlines lt max(wavevec-50.)),ngood)
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

pro osredx_velmap,file,zsys,lrest,ltype,gfwhm=gfwhm,mincov=mincov,tag=tag,outdir=outdir,maxvel=maxvel,maxsig=maxsig,minsig=minsig,verbose=verbose,skyfile=skyfile,sncut=sncut,writesb=writesb

; Output subdirectory
if (~keyword_set(outdir)) then outdir='./velfit/'
; Create output directory
if file_test(outdir,/directory) eq 0 then spawn, '\mkdir -p '+outdir

; Memory tracking
thismem = memory()
maxmem = 0
;Open and start a log file
logfile=concat_dir(outdir,'velfit.log')
;cpbackup, logfile
splog, filename=logfile
stime0 = systime(1)
splog, 'Log file ' + logfile + ' opened ' + systime()

  
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

; QA parameters
if (~keyword_set(sncut)) then sncut=6

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

for i=0,nz-1 do begin
  errspec[i]=median(err[*,*,i])
endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Apply spatial smoothing per wavelength channel and
; mask the low-coverage regions.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

newflux=flux
covmask=median(covmap,dimension=3)
goodindx=where(covmask ge mincov,complement=badindx,ngood)
; Convert to x,y good locations
ind=array_indices(covmask,goodindx)
goodx=ind[0,*]
goody=ind[1,*]
splog,'ngood=',ngood

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


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Initialize fitting defaults
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; default of the median sky spectrum is the error spectrum
; (indicative of skylines ONLY, not right in any numeric sense,
; will generally replace this later)
medianskyspec=errspec

; 4 parameters; wavelength, sigma, amplitude, baseline
parinfo = replicate({value:0.D, fixed:0, limited:[0,0], $
                       limits:[0.D,0]}, 4)

; Default fitting parameters
if (~keyword_set(maxvel)) then maxvel=1000. ; +- 1000 km/s
if (~keyword_set(maxsig)) then maxsig=1000. ; 1000 km/s
if (~keyword_set(minsig)) then minsig=-1 ; <0 means use the LSF as a lower bound
if (~keyword_set(maxsig)) then maxsig=1000. ; 1000 km/s

c=3e5; speed of light in km/s
lambda_guess=(1.+zsys)*lrest
lambda_ul=(1.d + maxvel/c)*lambda_guess; Upper limit on wavelength in Angstrom
lambda_ll=(1.d - maxvel/c)*lambda_guess; Lower limit on wavelength in Angstrom
sigma_ll=(minsig/c)*lambda_guess; Lower limit on sigma in Angstrom
sigma_ul=(maxsig/c)*lambda_guess; Upper limit on sigma in Angstrom

; Parameters are [wavelength, sigma, amplitude, baseline]
parinfo[*].value = [lambda_guess, 6., 0.01, 0.]; starting values
; Initial guess
parms0=[lambda_guess, 6., 0.1, 0.]
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
; Require positive amplitude
parinfo[2].limited[0] = 1
parinfo[2].limits[0]  = 0.D

velmap=replicate(-999.,nx,ny)
sigmap=replicate(-999.,nx,ny)
velmap_err=replicate(-999.,nx,ny)
sigmap_err=replicate(-999.,nx,ny)
fluxmap=fltarr(nx,ny)
snmap=fltarr(nx,ny)

; Create a map of the LSF (sigma in km/s).  By default it is constant
; and set to typical OSIRIS properties
lsfmap=replicate(37.5,nx,ny); Sigma in km/s

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; First loop will do an initial fit to determine the typical
; covariance correction to error spectrum due to our smoothing.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Loop over at most 100 spectra, anything more is overkill
covar=fltarr(ngood<100)
for i=0,(ngood<100)-1 do begin
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
    parinfo[1].limits[0]=lambda_guess/thisR/2.35

  ; First fit uses the input error spectrum
  if (ltype eq 'o3double') then $
    parms=mpfitfun('o3_twocomponent',wavevec,spectrum,errspec,parms0,parinfo=parinfo,yfit=yfit,perror=perr,/quiet)
  if (ltype eq 'single') then $
    parms=mpfitfun('onecomponent',wavevec,spectrum,errspec,parms0,parinfo=parinfo,yfit=yfit,perror=perr,/quiet)
  
  ; Compute model-subtracted spectrum; compare residuals against 
  ; nominal error spectrum to figure out a covariance correction factor
  ; that accounts for our spatial smoothing
  diff=spectrum-yfit
  rms=sqrt((moment(diff/errspec))[1])
  covar[i]=rms*sqrt(fxpar(hdr0,'CDELT1')*10.); Error spectrum should have 1/sqrt(pixsize) standard dev
endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Second pass loop using the error spectrum corrected for typical
; covariance derived from the initial fit
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Apply the typical covariance correction to the error spectrum
errspec=errspec*median(covar)

for i=0,ngood-1 do begin
  spectrum=newflux[goodx[i],goody[i],*]

  ; Optional sky spectrum fitting to get the LSF
  if (keyword_set(skyfile)) then begin
    skyspec=skycube[goodx[i],goody[i],*]
    skyR=fit_skyspec(wavevec, skyspec); spectral resolution
    thislsf=c/skyR/2.35
    lsfmap[goodx[i],goody[i]]=thislsf; lsf map of sigma in km/s
  endif else begin
    thislsf=lsfmap[goodx[i],goody[i]]
  endelse

  ; If minsig parameter was less than zero, update
  ; the bounded range for the fitter so that sigma cannot
  ; be below instrumental resolution in Angstroms
  thisR=c/lsfmap[goodx[i],goody[i]]/2.35
  if (minsig lt 0) then $
    parinfo[1].limits[0]=lambda_guess/thisR/2.35
  
  ; Re-fit using covariance-corrected error spectrum
  if (ltype eq 'o3double') then $
    parms=mpfitfun('o3_twocomponent',wavevec,spectrum,errspec,parms0,parinfo=parinfo,yfit=yfit,perror=perr,/quiet,bestnorm=bestnorm,dof=dof)
  if (ltype eq 'single') then $
    parms=mpfitfun('onecomponent',wavevec,spectrum,errspec,parms0,parinfo=parinfo,yfit=yfit,perror=perr,/quiet,bestnorm=bestnorm,dof=dof)

  ; Flux is amplitude x sigma x sqrt(2pi)
  theflux=parms[2]*parms[1]*sqrt(2*!PI)
  ; Error propagation
  dflux=theflux*sqrt((perr[2]/parms[2])^2 + (perr[1]/parms[1])^2)*sqrt(bestnorm/dof)
  thesnr=theflux/dflux
  thevel=(parms[0]-lambda_guess)/lambda_guess*3e5
  thevel_err=(perr[0])/lambda_guess*3e5*sqrt(bestnorm/dof)
  thesig=(parms[1])/parms[0]*3e5
  thesig_err=perr[1]/parms[0]*3e5*sqrt(bestnorm/dof)
  
  ; Subtract instrumental line width in quadrature
  thislsf=thislsf*lambda_guess/parms[0]; Correct for difference in guess vs real line center
  if (keyword_set(verbose)) then begin
    splog,'Fitted snr: ',thesnr
    splog,'Fitted sigma: ',thesig,' km/s'
    splog,'Instrumental sigma: ',thislsf,' km/s'
  endif
  thesig=sqrt(thesig*thesig-thislsf*thislsf)
  if (~finite(thesig)) then thesig=0.
  if (keyword_set(verbose)) then $
    splog,'Corrected sigma: ',thesig,' km/s'

    ; Remove any bad fits with zero or negative flux
    if (theflux eq 0.) then begin
      thevel=-999.
      thesig=-999.
      thevel_err=-999.
      thesig_err=-999.
      theflux=0.
      thesnr=0.
    endif

    ; Remove bad fits with below threshhold snr
    if ((thesnr lt sncut)or(~(finite(thesnr)))) then begin
      thevel=-999.
      thesig=-999.
      thevel_err=-999.
      thesig_err=-999.
      theflux=0.
      thesnr=0.
    endif
    
    ; Remove bad fits within 5% of upper sigma bound
    if (thesig gt (0.95 * maxsig)) then begin
      thevel=-999.
      thesig=-999.
      thevel_err=-999.
      thesig_err=-999.
      theflux=0.
      thesnr=0.
    endif
   
  velmap[goodx[i],goody[i]]=thevel
  sigmap[goodx[i],goody[i]]=thesig
  velmap_err[goodx[i],goody[i]]=thevel_err
  sigmap_err[goodx[i],goody[i]]=thesig_err
  fluxmap[goodx[i],goody[i]]=theflux
  snmap[goodx[i],goody[i]]=thesnr

  if ((keyword_set(verbose))and(snmap[goodx[i],goody[i]] gt sncut)) then begin
    plot,wavevec,spectrum
    oplot,wavevec,yfit,color=250
    thestring=strcompress('SNR = '+string(snmap[goodx[i],goody[i]]))
    xyouts,median(wavevec),median(spectrum)+2/3.*(max(spectrum)-median(spectrum)),thestring,color=250,charsize=2
    wait, 1.0
  endif
endfor
   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Make composite spectra from all of the spaxels that passed the
; quality cuts
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Composite spectrum of all regions
indx_all=where(snmap gt sncut,nhighsn_all)
ind_all=array_indices(snmap,indx_all)
goodx_all=ind_all[0,*]
goody_all=ind_all[1,*]
sumspec_all=fltarr(nz)
for i=0,nhighsn_all-1 do $
  sumspec_all += newflux[goodx_all[i],goody_all[i],*]
if (keyword_set(verbose)) then splog,'nhighsn_all=',nhighsn_all

; Composite spectrum of low-dispersion regions
indx_low=where((snmap gt sncut)and(sigmap lt 300.),nhighsn_low)
if (nhighsn_low gt 0) then begin
  ind_low=array_indices(snmap,indx_low)
  goodx_low=ind_low[0,*]
  goody_low=ind_low[1,*]
  sumspec_low=fltarr(nz)
  for i=0,nhighsn_low-1 do $
    sumspec_low += newflux[goodx_low[i],goody_low[i],*]
  if (keyword_set(verbose)) then splog,'nhighsn_low=',nhighsn_low
endif

; Composite spectrum of high-dispersion regions
indx_high=where((snmap gt sncut)and(sigmap ge 300.),nhighsn_high)
if (nhighsn_high gt 0) then begin
  ind_high=array_indices(snmap,indx_high)
  goodx_high=ind_high[0,*]
  goody_high=ind_high[1,*]
  sumspec_high=fltarr(nz)
  for i=0,nhighsn_high-1 do $
    sumspec_high += newflux[goodx_high[i],goody_high[i],*]
  if (keyword_set(verbose)) then splog,'nhighsn_high=',nhighsn_high
endif

; Median sky spectrum across all high snr regions
if (keyword_set(skyfile)) then begin
  medianskyspec[*]=0.
  for i=0,nz-1 do begin
    temp=skycube[*,*,i]
    medianskyspec[i]=median(temp[goodx_all,goody_all]) 
  endfor
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Fit the composite spectra
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Need to use round() on covar; otherwise get some odd numerical artifacts that
; are not meaningful

; Fit the composite spectrum of all regions
; First fit uses incorrect error spectrum (not scaled)
if (ltype eq 'o3double') then $
  parms_all=mpfitfun('o3_twocomponent',wavevec,sumspec_all,errspec,parms0,parinfo=parinfo,yfit=yfit_all,perror=perr_all,/quiet)
if (ltype eq 'single') then $
  parms_all=mpfitfun('onecomponent',wavevec,sumspec_all,errspec,parms0,parinfo=parinfo,yfit=yfit_all,perror=perr_all,/quiet)
; Derive error correction term to account for coadded covariant spectra
diff=sumspec_all-yfit_all
rms=sqrt((moment(diff/errspec))[1])
covar=round(rms*sqrt(fxpar(hdr0,'CDELT1')*10.))
; Second fit uses corrected errors
if (ltype eq 'o3double') then $
  parms_all=mpfitfun('o3_twocomponent',wavevec,sumspec_all,errspec*covar,parms0,parinfo=parinfo,yfit=yfit_all,perror=perr_all,/quiet,bestnorm=bestnorm,dof=dof)
if (ltype eq 'single') then $
  parms_all=mpfitfun('onecomponent',wavevec,sumspec_all,errspec*covar,parms0,parinfo=parinfo,yfit=yfit_all,perror=perr_all,/quiet,bestnorm=bestnorm,dof=dof)
theflux=parms_all[2]*parms_all[1]*sqrt(2*!PI)
dflux=theflux*sqrt((perr_all[2]/parms_all[2])^2 + (perr_all[1]/parms_all[1])^2)*sqrt(bestnorm/dof)
thevel=(parms_all[0]-lambda_guess)/lambda_guess*3e5
thevel_err=(perr_all[0])/lambda_guess*3e5*sqrt(bestnorm/dof)
thesig=(parms_all[1])/parms_all[0]*3e5
thesig_err=perr_all[1]/parms_all[0]*3e5*sqrt(bestnorm/dof)
splog,'Sum_all flux: ',theflux/10.,' +- ',dflux/10.,' e-16 erg/s/cm2'
splog,strcompress('Sum_all wave: '+string(parms_all[0])+' um, (z='+string(parms_all[0]/lrest-1.)+')')
splog,'Sum_all velocity: ',thevel,' +- ',thevel_err,' km/s'
splog,'Sum_all sigma: ',thesig,' +- ',thesig_err,' km/s'
splog,'SNR: ',theflux/dflux

; Fit the composite spectrum of low sigma regions
; First fit uses incorrect error spectrum (not scaled)
if (nhighsn_low gt 0) then begin
if (ltype eq 'o3double') then $
  parms_low=mpfitfun('o3_twocomponent',wavevec,sumspec_low,errspec,parms0,parinfo=parinfo,yfit=yfit_low,perror=perr_low,/quiet)
if (ltype eq 'single') then $
  parms_low=mpfitfun('onecomponent',wavevec,sumspec_low,errspec,parms0,parinfo=parinfo,yfit=yfit_low,perror=perr_low,/quiet)
; Derive error correction term to account for coadded covariant spectra
diff=sumspec_low-yfit_low
rms=sqrt((moment(diff/errspec))[1])
covar=round(rms*sqrt(fxpar(hdr0,'CDELT1')*10.))
; Second fit uses corrected errors
if (ltype eq 'o3double') then $
  parms_low=mpfitfun('o3_twocomponent',wavevec,sumspec_low,errspec*covar,parms0,parinfo=parinfo,yfit=yfit_low,perror=perr_low,/quiet,bestnorm=bestnorm,dof=dof)
if (ltype eq 'single') then $
  parms_low=mpfitfun('onecomponent',wavevec,sumspec_low,errspec*covar,parms0,parinfo=parinfo,yfit=yfit_low,perror=perr_low,/quiet,bestnorm=bestnorm,dof=dof)
theflux=parms_low[2]*parms_low[1]*sqrt(2*!PI)
dflux=theflux*sqrt((perr_low[2]/parms_low[2])^2 + (perr_low[1]/parms_low[1])^2)*sqrt(bestnorm/dof)
thevel=(parms_low[0]-lambda_guess)/lambda_guess*3e5
thevel_err=(perr_low[0])/lambda_guess*3e5*sqrt(bestnorm/dof)
thesig=(parms_low[1])/parms_low[0]*3e5
thesig_err=perr_low[1]/parms_low[0]*3e5*sqrt(bestnorm/dof)
splog,'Sum_low flux: ',theflux/10.,' +- ',dflux/10.,' e-16 erg/s/cm2'
splog,strcompress('Sum_low wave: '+string(parms_low[0])+' um, (z='+string(parms_low[0]/lrest-1.)+')')
splog,'Sum_low velocity: ',thevel,' +- ',thevel_err,' km/s'
splog,'Sum_low sigma: ',thesig,' +- ',thesig_err,' km/s'
splog,'SNR: ',theflux/dflux
endif

; Fit the composite spectrum of high sigma regions
; First fit uses incorrect error spectrum (not scaled)
if (nhighsn_high gt 0) then begin
if (ltype eq 'o3double') then $
  parms_high=mpfitfun('o3_twocomponent',wavevec,sumspec_high,errspec,parms0,parinfo=parinfo,yfit=yfit_high,perror=perr_high,/quiet)
if (ltype eq 'single') then $
  parms_high=mpfitfun('onecomponent',wavevec,sumspec_high,errspec,parms0,parinfo=parinfo,yfit=yfit_high,perror=perr_high,/quiet)
; Derive error correction term to account for coadded covariant spectra
diff=sumspec_high-yfit_high
rms=sqrt((moment(diff/errspec))[1])
covar=round(rms*sqrt(fxpar(hdr0,'CDELT1')*10.))
; Second fit uses corrected errors
if (ltype eq 'o3double') then $
  parms_high=mpfitfun('o3_twocomponent',wavevec,sumspec_high,errspec*covar,parms0,parinfo=parinfo,yfit=yfit_high,perror=perr_high,/quiet,bestnorm=bestnorm,dof=dof)
if (ltype eq 'single') then $
  parms_high=mpfitfun('onecomponent',wavevec,sumspec_high,errspec*covar,parms0,parinfo=parinfo,yfit=yfit_high,perror=perr_high,/quiet,bestnorm=bestnorm,dof=dof)
theflux=parms_high[2]*parms_high[1]*sqrt(2*!PI)
dflux=theflux*sqrt((perr_high[2]/parms_high[2])^2 + (perr_high[1]/parms_high[1])^2)*sqrt(bestnorm/dof)
thevel=(parms_high[0]-lambda_guess)/lambda_guess*3e5
thevel_err=(perr_high[0])/lambda_guess*3e5*sqrt(bestnorm/dof)
thesig=(parms_high[1])/parms_high[0]*3e5
thesig_err=perr_high[1]/parms_high[0]*3e5*sqrt(bestnorm/dof)
splog,'Sum_high flux: ',theflux/10.,' +- ',dflux/10.,' e-16 erg/s/cm2'
splog,strcompress('Sum_high wave: '+string(parms_high[0])+' um, (z='+string(parms_high[0]/lrest-1.)+')')
splog,'Sum_high velocity: ',thevel,' +- ',thevel_err,' km/s'
splog,'Sum_high sigma: ',thesig,' +- ',thesig_err,' km/s'
splog,'SNR: ',theflux/dflux
endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Convert flux map to surf-b map
; in units of /arcsec^2
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
dx=fxpar(hdr0,'CDELT2')*3600.   ;arcsec/spaxel
surfb=fluxmap/dx/dx/10.; 1e-16 erg/s/cm/arcsec2

  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Write out the results
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

writefits,concat_dir(outdir,'wavevec.fits'),wavevec
writefits,concat_dir(outdir,'sumspec_all.fits'),sumspec_all
writefits,concat_dir(outdir,'skyspec.fits'),medianskyspec
writefits,concat_dir(outdir,'yfit_all.fits'),yfit_all
if (nhighsn_low gt 0) then begin
  writefits,concat_dir(outdir,'sumspec_low.fits'),sumspec_low
  writefits,concat_dir(outdir,'yfit_low.fits'),yfit_low
endif
if (nhighsn_high gt 0) then begin
  writefits,concat_dir(outdir,'sumspec_high.fits'),sumspec_high
  writefits,concat_dir(outdir,'yfit_high.fits'),yfit_high
endif

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
writefits,concat_dir(outdir,'velmap_err.fits'),velmap_err,hdrim
writefits,concat_dir(outdir,'sigmap_err.fits'),sigmap_err,hdrim
writefits,concat_dir(outdir,'snmap.fits'),snmap,hdrim
if (keyword_set(skyfile)) then $
  writefits,concat_dir(outdir,'lsfmap.fits'),lsfmap,hdrim
if (keyword_set(writesb)) then $
  writefits,concat_dir(outdir,'surfbmap.fits'),surfb,hdrim

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Clean up
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Track memory usage
thismem = memory()
maxmem = maxmem > thismem[3]
splog, 'Max memory usage = ', string(maxmem/1e6,format='(f7.1)'), ' MB'

splog, 'Total time for osredx_velmap = ', systime(1)-stime0, ' seconds', format='(a,f6.0,a)'
splog, 'Successful completion of osredx_velmap at ' + systime()

; Close logfile
if (keyword_set(logfile)) then splog, /close

return
end
