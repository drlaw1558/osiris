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

pro osredx_velmap

device,decomposed=0
loadct,39

file='q0142-bx195-o3.fits'
flux=mrdfits(file,0,hdr0)
err=mrdfits(file,1,hdr1)
covmap=mrdfits(file,3,hdr3)

flux=transpose(flux,[2,1,0])
err=transpose(err,[2,1,0])
covmap=transpose(covmap,[2,1,0])

nx=fxpar(hdr0,'NAXIS3')
ny=fxpar(hdr0,'NAXIS2')
nz=fxpar(hdr0,'NAXIS1')

; Average the error vector over entire field to get skylines
errspec=fltarr(nz)
; Error spectrum seems weird normalization, apply an empirical
; factor such that pure-noise spectra of individual spaxels
; divided by the noise vector is gaussian with sigma=1
fudge=15.
for i=0,nz-1 do begin
  errspec[i]=median(err[*,*,i])/fudge
endfor


; Get the wavelength vector (vacuum)
wavevec=fxpar(hdr0,'CRVAL1')+findgen(nz)*fxpar(hdr0,'CDELT1')
; Heliocentric correction
; Kinda rough, but good enough
; (Works on assumption SAME for all exposures)
ra=fxpar(hdr0,'RA')
dec=fxpar(hdr0,'DEC')
jd=fxpar(hdr0,'MJD-OBS')+2400000.5d
long=360.-155.4747d
lat=19.8260
alt=4205.
helio=heliocentric(ra,dec,jd=jd,longitude=long,latitude=lat,altitude=alt)
; Correct wave vector and convert to Angstroms
wavevec = wavevec / (1 + helio/299792.458) * 10

; Coverage threshhold
maxcov=max(covmap)
mincov=0.7*maxcov

; Spatial smoothing
newflux=flux
for i=0,nz-1 do begin
  newim=filter_image(flux[*,*,i],fwhm_gaussian=2.5)
  temp=covmap[*,*,i]
  indx=where(temp lt mincov)
  newim[indx]=0.
  newflux[*,*,i]=newim
endfor
; write out smoothed cube
temp=transpose(newflux,[2,1,0])
writefits,'temp.fits',temp,hdr0


; Median combine to get an overall mask for analysis
themed=median(newflux,dimension=3)
goodval=where(themed ne 0,ngood)
; Convert to x,y good locations
ind=array_indices(themed,goodval)
goodx=ind[0,*]
goody=ind[1,*]
print,'ngood=',ngood

; 4 parameters; wavelength, sigma, flux, baseline
parinfo = replicate({value:0.D, fixed:0, limited:[0,0], $
                       limits:[0.D,0]}, 4)
parinfo[*].value = [16918., 20., 0.01, 0.]; starting values
; Initial guess
parms0=[16918., 20., 0.1, 0.]
;Bound the velocities
parinfo[0].limited[0] = 1
parinfo[0].limited[1] = 1
parinfo[0].limits[0]  = 16844.D
parinfo[0].limits[1]  = 16986.D
; Bound the sigma
parinfo[1].limited[0] = 1
parinfo[1].limited[1] = 1
parinfo[1].limits[0]  = 5.D
parinfo[1].limits[1]  = 100.D
; Require positive flux
parinfo[2].limited[0] = 1
parinfo[2].limits[0]  = 0.D

wsys=3.3807*5008.24
velmap=replicate(-999.,nx,ny)
sigmap=replicate(-999.,nx,ny)
fluxmap=fltarr(nx,ny)
snmap=fltarr(nx,ny)
sncut=6
for i=0,ngood-1 do begin
  ;i=471
  spectrum=newflux[goodx[i],goody[i],*]
  parms=mpfitfun('o3_twocomponent',wavevec,spectrum,errspec,parms0,parinfo=parinfo,yfit=yfit,perror=perr,/quiet)
  diff=spectrum-yfit
  rms=sqrt((moment(diff/errspec))[1])
  print,i,' Noise scaling rms: ',rms
  plot,wavevec,spectrum
  oplot,wavevec,yfit,color=250

  thevel=(parms[0]-wsys)/wsys*3e5
  thesig=(parms[1])/wsys*3e5
  theflux=(parms[2])
  thesnr=parms[2]/perr[2]

  if (theflux eq 0.) then begin
    thevel=-999.
    thesig=-999.
    theflux=0.
    thesnr=0.
  endif

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

  ;if (snmap[goodx[i],goody[i]] gt 5.) then wait, 1.0
endfor

writefits,'fluxmap.fits',fluxmap
writefits,'velmap.fits',velmap  
writefits,'sigmap.fits',sigmap  
writefits,'snmap.fits',snmap 

stop


return
end
