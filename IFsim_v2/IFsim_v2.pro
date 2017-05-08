pro IFsim_v2
; Program IFsim_v2
; This is the IDL port of the OSIRIS simulation code
; This version 2.2 for release.
;
; Author:
;   David R. Law (Caltech/UCLA)
;   drlaw@astro.ucla.edu
;   First written: 10/2/06
;   Last modified: 10/8/09
;
; Modified on 9/5/08 to update to IDL 7.0
; IDL no longer accepts more than 8 parameters in a function call
; therefore had to condense them into a vector of parameters
; This regrettably makes the code much harder to read.
; Modified 9/22/08 to include all dependent subroutines within this file

; Modified 10/7/08 to allow me to read in and simulate TSG UVLGs
;
; Modified 9/24/09 to cut out some algorithms not needed for public
; release, add in some more general features.
;
; Modified 10/8/09 to add in wavelength info to fits header

; Set default values for basic globals
  SKYFILE=''
  TRANSFILE=''
  PARMVEC=dblarr(18)
  SetParms,PARMVEC,SKYFILE,TRANSFILE

  print,'Tel size: ',PARMVEC[0], ' cm'

; Find out what to do
  MainMenu,Task

; And do it
  IF(Task EQ 1) THEN BEGIN
    Task1,PARMVEC,SKYFILE,TRANSFILE
  ENDIF
  IF(Task EQ 2) THEN BEGIN
    Task2,PARMVEC,SKYFILE,TRANSFILE
  ENDIF
  IF(Task EQ 3) THEN BEGIN
    return
  ENDIF
end

; Pre-define some later functions
Forward_Function PixelScale, LineFluxToEFlux, MakeSNmap
Forward_Function MakeModelCube


;--------------------------------------------------------------------------------------------

; Basic bleeding due to finite PSF
; Bleeds a data cube

function BasicBleed,InCube,gfwhm_pix
  OutCube=InCube
  for k=0,(size(InCube))[3]-1 do OutCube[*,*,k]=filter_image(InCube[*,*,k],fwhm_gaussian=gfwhm_pix)
  
  return,OutCube
end

;--------------------------------------------------------------------------------------------
; Return a 2d moffat function array
; Give it an array size, function center
; alpha, and beta
; Alpha controls how fat it is, beta how fast the tail falls.
; Be warned that values beta<=1 may give bizarre results

function moffatIFS,xsize,ysize,xcen,ycen,alpha,beta,totalflux

img=fltarr(xsize,ysize)

for i=0,xsize-1 do begin
  for j=0,ysize-1 do begin
    radius=sqrt((i-xcen)^2+(j-ycen)^2)
    value=(1+(radius/alpha)^2)^(-beta)
    img[i,j]=value
  endfor
endfor

; Scale to desired total flux
proftot=3.14159*alpha*alpha/(beta-1)
img=img*totalflux/proftot

return,img
end

;--------------------------------------------------------------------------------------------
; Similarly to how filter_image smooths an image by a Gaussian
; kernel, moffatblur.pro smooths an image by a Moffat profile
; Input is a 2-d flux table, as is output
; Use alpha = 92 mas (converted to pixel units), and beta=1.6??

function moffatblurIFS,InpImg,alpha,beta

xsize=(size(InpImg))[1]
ysize=(size(InpImg))[2]
OutpImg=fltarr(xsize,ysize)

for i=0,xsize-1 do begin
  for j=0,ysize-1 do begin
    if (InpImg[i,j] ne 0.) then begin
      temp2=moffatIFS(xsize,ysize,i,j,alpha,beta,InpImg[i,j])
      OutpImg=OutpImg+temp2
    endif
  endfor
endfor

return,OutpImg
end

;--------------------------------------------------------------------------------------------

; Main menu subroutine
; Output main menu to terminal, get back the choice of which
; task to perform
pro MainMenu,Task
  print,'**************************************'
  print,'OSIRIS Simulation Interface (IDL port): v2.1'
  print,'David R. Law <drlaw@astro.ucla.edu>'
  print,'Last modified September 24, 2009'
  print,'Choose one of the following actions to perform:'
  print,'(1) Generate a simulated data cube and S/N ratio map'
  print,'(2) Combine a flux, velocity, sigma map and observe.'
  print,'(3) Quit.'
  read,Task,prompt='Choice: '
  return
  end

;--------------------------------------------------------------------------------------------
  
; Subroutine to read the map in from a data file, rotate it, and rebin
; to desired lenslet scale
pro InputData,InpFile,InpPixScale,InpRot,Xobject,Yobject,PltScl,Filter,Data
  ; Get FOV from subroutine
  GetFOV,Filter,PltScl,xfov,yfov
  ; Figure out what size region in pixels to crop out
  xcrop=xfov/2./InpPixScale
  ycrop=yfov/2./InpPixScale

  ; First read in large enough chunk that can do arbitrary rotation
  ; without clipping the image.  Pick 2x largest dimension
  clipsize=xcrop
  if (ycrop gt xcrop) then clipsize=ycrop
  fxread,InpFile,data0,header,LONG(Xobject-2*clipsize),LONG(Xobject+2*clipsize),LONG(Yobject-2*clipsize),LONG(Yobject+2*clipsize)

  ; Rotate it to desired orientation.  Pivot about target coordinates
  ; so that coordinates stay the same in rotated image.
  ; Note that the 1.0 is for magnification
  data1=rot(data0,InpRot,1.0,2*clipsize,2*clipsize,/pivot,/interp)

  ; Crop out the section of the rotated image that we want
  data2=data1[2*clipsize-xcrop:2*clipsize+xcrop,2*clipsize-ycrop:2*clipsize+ycrop]

  ; Rebin to desired lenslet scale
  xelem=xfov/PltScl
  yelem=yfov/PltScl
  IF (PltScl NE InpPixScale) THEN BEGIN
    Data=frebin(data2,xelem,yelem,/total)
  ENDIF ELSE BEGIN
    Data=data2
  ENDELSE

  return
  end

;--------------------------------------------------------------------------------------------

; Returns the FWHM of the Gaussian seeing psf in arcseconds at a given lambda.
; NOTE- Looks like program happier if 'function's come above 'pro's in this code
function FWHMSeeing,lambda,PARMVEC
  return,PARMVEC[12]*((lambda/5500.)^(-0.2))
  end

;--------------------------------------------------------------------------------------------

; Returns the FWHM of the (modified Gaussian scaled by 2.4/2.95) diffraction limited 
; psf in arcseconds at a given lambda
function FWHMDiffraction,lambda,PARMVEC
; 1.22 lambda/d, wavelength to cm, divide for sigma, convert to arcseconds from radians
  return,1.22*lambda*1e-8/PARMVEC[0]/2.95*2.4*(180.*3600./3.1415926536)
  end

;--------------------------------------------------------------------------------------------

; Returns the strehl of AO correction at a given lambda using exponential formula with wavefront error
function GetStrehl,lambda,PARMVEC
  return, 2.718281828^(-((2.*3.1415926536*PARMVEC[13]*10./lambda)^2.))
  end

;--------------------------------------------------------------------------------------------

; Calculate the pixel scale (Angstroms/pixel) at a given lambda
; Mimics OSIRIS resampled values for reduced data cubes
function PixelScale,lambda
  value=-1.
  IF ((lambda GT 10190.)AND(lambda LT 11740)) THEN value=1.2 ELSE $
  IF ((lambda GT 11740.)AND(lambda LT 13890)) THEN value=1.5 ELSE $
  IF ((lambda GT 14660.)AND(lambda LT 18080)) THEN value=2.0 ELSE $
  IF ((lambda GT 19550.)AND(lambda LT 24080)) THEN value=2.5
  return,value
  end
  
;--------------------------------------------------------------------------------------------

; Returns the average throughput for the narrowband filter used at the
; given wavelength.  Note that this doesn't account for *exact*
; filter response curves, but works from average throughput values.
; This function is of course very specific to OSIRIS
function FilterTput,lambda
  tput=0.
  ; z bands
  IF((lambda GT 10190.)AND(lambda LT 10700.)) THEN tput=0.7 ELSE $
  IF ((lambda GT 10600.)AND(lambda LT 11130.)) THEN tput=0.72 ELSE $
  IF ((lambda GT 11030.)AND(lambda LT 11580.)) THEN tput=0.78 ELSE $
  IF ((lambda GT 11480.)AND(lambda LT 11740.)) THEN tput=0.73   ELSE $
  ;Overlap of filters: jump to J1 at 1174 nm
  ; J bands
  IF ((lambda GT 11740.)AND(lambda LT 12320.)) THEN tput=0.78 ELSE $
  IF ((lambda GT 12280.)AND(lambda LT 12890.)) THEN tput=0.78 ELSE $
  IF ((lambda GT 12750.)AND(lambda LT 13390.)) THEN tput=0.84 ELSE $
  IF ((lambda GT 13230.)AND(lambda LT 13890.)) THEN tput=0.82 ELSE $
  ; H bands
  IF ((lambda GT 14660.)AND(lambda LT 15410.)) THEN tput=0.81 ELSE $
  IF ((lambda GT 15320.)AND(lambda LT 16100.)) THEN tput=0.75 ELSE $
  IF ((lambda GT 15940.)AND(lambda LT 16760.)) THEN tput=0.80 ELSE $
  IF ((lambda GT 16520.)AND(lambda LT 17370.)) THEN tput=0.83 ELSE $
  IF ((lambda GT 17210.)AND(lambda LT 18080.)) THEN tput=0.75 ELSE $
  ; K bands
  IF ((lambda GT 19550.)AND(lambda LT 20550.)) THEN tput=0.85 ELSE $
  IF ((lambda GT 20360.)AND(lambda LT 21410.)) THEN tput=0.83 ELSE $
  IF ((lambda GT 21210.)AND(lambda LT 22290.)) THEN tput=0.84 ELSE $
  IF ((lambda GT 22080.)AND(lambda LT 23200.)) THEN tput=0.75 ELSE $
  IF ((lambda GT 22920.)AND(lambda LT 24080.)) THEN tput=0.8

  return,tput
  end

;--------------------------------------------------------------------------------------------

; Returns the grating throughput at a given wavelength (assumes ideal
; grating), modelled as sinc^2 function.
; This function is of course quite specific to OSIRIS
function BlazeTput,lambda
  PI=3.1415926536
  m=0.; Order
  lc0=6.5*1e4; Blazed at 6.5 microns (in Angstroms)

  IF ((lambda GT 10190.)AND(lambda LT 11740.)) THEN m=6 ELSE $
  IF ((lambda GT 11740.)AND(lambda LT 13890.)) THEN m=5 ELSE $
  IF ((lambda GT 14660.)AND(lambda LT 18080.)) THEN m=4 ELSE $
  IF ((lambda GT 19550.)AND(lambda LT 24080.)) THEN m=3

  lc=lc0/m; Blaze wavelength in an order
  X=PI*m*(lambda/lc-1.); Argument of blaze function
  Answer=0.
  IF (abs(X) GT 1e-8) THEN BEGIN; Avoid dividing by zero
    Answer=(sin(X)/X)^2
  ENDIF
  
 ; Grating is not perfect, based on Larkin's throughput results then
 ; Grating should be taken into account in ptrans parameter
 ; Answer=0.5*Answer

  return, Answer
  end

;--------------------------------------------------------------------------------------------

; Returns the atmospheric transmission at a certain wavelength
function AtmosTput,lambda,LineWidth,TRANSFILE
  ; Read the background sky flux values for all lambda from file
  READCOL,TRANSFILE,translam,transval,SKIPLINE=2,/silent
; Multiply wavelengths by 1e4 to convert from microns to Angstroms
  translam=translam*1e4

; Compute the average sky transmission value for wavelength covered by the
; width of your line.
  avgskytrans=0.
  navg=0.
  FOR i=0,(size(translam))[1]-1 DO BEGIN
    IF((translam[i] GE (lambda-LineWidth/2.))AND(translam[i] LE (lambda+LineWidth/2.))) THEN BEGIN
      avgskytrans=avgskytrans+transval[i]
      navg=navg+1
    ENDIF
  ENDFOR
  avgskytrans=avgskytrans/navg

  IF (navg EQ 0) THEN print,'WARNING- Sky transmission not defined!'

  return,avgskytrans
  end

;--------------------------------------------------------------------------------------------

; Returns the sky flux at a certain wavelength
function SkyFlux,lambda,LineWidth,SKYFILE
  ; Read the background sky flux values for all lambda from file (#/s/nm/m^2/as^2)
  READCOL,SKYFILE,skylam,skyval,SKIPLINE=10,/SILENT
; Multiply wavelengths by ten to convert from nm to Angstroms
  skylam=skylam*10.
; Multiply fluxes input from file by 1e3 to convert to #/s/cm^2/cm/as^2
  skyval=skyval*1e3

; Compute the average sky flux value for wavelength covered by the
; width of your line.  This is average sky background in #/s/cm^2/cm/as^2
  avgskyflux=0.
  navg=0.
  FOR i=0,(size(skylam))[1]-1 DO BEGIN
    IF((skylam[i] GE (lambda-LineWidth/2.))AND(skylam[i] LE (lambda+LineWidth/2.))) THEN BEGIN
      avgskyflux=avgskyflux+skyval[i]
      navg=navg+1
    ENDIF
  ENDFOR
  avgskyflux=avgskyflux/navg

  return,avgskyflux
  end

;--------------------------------------------------------------------------------------------

; Returns the background flux at a certain wavelength in #/s/lenslet/pixel
function BackgroundFlux,PARMVEC,lambda,LineWidth,PltScl,type,SKYFILE,TRANSFILE
  PI=3.1415926536
  radpas=4.8481e-6; Radians per arcsecond
  c=2.9979e10; Speed of light (cm/s)
  hc=1.986e-16; Planck constant times speed of light (erg cm)
  k=1.381e-16; Boltzmann constant (erg/K)
  A=PI*PARMVEC[0]*PARMVEC[0]/4.; Area of primary (cm^2)
  PS=PixelScale(lambda); Pixel scale (Angstroms/pixel)

; Calculate the background sky flux
  avgskyflux=SkyFlux(lambda,LineWidth,SKYFILE)
; Figure out the flux of thermal photons from telescope and AO system
; in #/s/cm^2/cm/as^2
  bgAO=PARMVEC[15]*radpas*radpas*2.*c/((lambda*1e-8)^4)/(exp(hc/(lambda*1e-8*k*PARMVEC[7]))-1.)
  bgTEL=PARMVEC[16]*radpas*radpas*2.*c/((lambda*1e-8)^4)/(exp(hc/(lambda*1e-8*k*PARMVEC[8]))-1.)

; Overload this function: if type=0 then return background in
; #/s/lenslet/pixel, otherwise in #/s/nm/m^2/arcsec^2
  IF(type EQ 0) THEN BEGIN
    retvalue=(avgskyflux*PARMVEC[11]*PARMVEC[9]+bgTEL*PARMVEC[9]+bgAO)*PARMVEC[10]*BlazeTput(lambda)*FilterTput(lambda)*A*PltScl*PltScl*PS*1e-8/PARMVEC[2]
  ENDIF ELSE BEGIN
    retvalue=(avgskyflux+bgTEL/PARMVEC[11]+bgAO/PARMVEC[11]/PARMVEC[9])*1e4/1e7
  ENDELSE

  return, retvalue
  end

;--------------------------------------------------------------------------------------------

; Given some flux incident on the atmosphere from the source in some
; spectral line (erg/s/cm^2/lenslet), find the resulting count rate it
; would produce on the detector in e-/s/lenslet/pixel
; Assumes line flux has box profile between FWHM points, ie gives
; AVERAGE counts/pixel.  This routine neglects strehl losses (dealt
; with elsewhere), but accounts for all throughput sources
function LineFluxToEFlux,PARMVEC,lambda,LineWidth,TRANSFILE
  PI=3.1415926536
  hc=1.986e-16; Planck constant times speed of light (erg cm)
  A=PI*PARMVEC[0]*PARMVEC[0]/4.; Area of primary (cm^2)
  PS=PixelScale(lambda); Pixel scale (Angstroms/pixel)

  Through=AtmosTput(lambda,LineWidth,TRANSFILE)*PARMVEC[11]*PARMVEC[9]*PARMVEC[10]*BlazeTput(lambda)*FilterTput(lambda)

  return,1.*lambda*1e-8/hc*A/LineWidth*PS/PARMVEC[2]*Through
  end

;--------------------------------------------------------------------------------------------

;This function combines a 2-d Strehl-bled flux map with background count values to produce a 
; 2-d S/N ratio map.  Input flux map will be propagated through all necessary throughputs, 
; background counts have already had this done.
; BledSignalData is table of source counts contributing to signal
; BledNoiseData is table of source counts contributing to noise
function MakeSNmap,PARMVEC,BledData,lambda,LineWidth,PltScl,SKYFILE,TRANSFILE,xtime
  ; Determine background flux
  Rbg=BackgroundFlux(PARMVEC,lambda,LineWidth,PltScl,0,SKYFILE,TRANSFILE)

  CalcStoN,PARMVEC,BledData,SNmap,Rbg,lambda,LineWidth,xtime
  ; Zero out negative values as unphysical
  zvals=where(SNmap lt 0.)
  if ((size(zvals))[0] ne 0) then begin
    SNmap[zvals]=0.
  endif
  return,SNmap
  end

;--------------------------------------------------------------------------------------------

; Given some count rate averaged over some number of spectral pixels,
; calculate StoN ratio in a given time xtime.  Requires input counts,
; non-signal background counts from source, background counts from sky
; and emissivity, wavelength, time, effective spectral resolution for
; number of spectral pixels summed
pro CalcStoN,PARMVEC,BledData,SNmap,Rbg,lambda,LineWidth,xtime
  PS=PixelScale(lambda); Pixel scale (Angstroms/pixel)
  Nspec=LineWidth/PS; Number of pixels summed in spectral direction

;  SNmap=BledData*N1*Nspec*xtime/sqrt((BledData+Rbg+RDK)*N1*Nspec*xtime+RN*RN*xtime/T0*N1*Nspec)
  SNmap=BledData*Nspec*PARMVEC[2]*xtime/sqrt((BledData+Rbg+PARMVEC[4])*Nspec*PARMVEC[2]*xtime+PARMVEC[5]*PARMVEC[5]*xtime/PARMVEC[6]*PARMVEC[2]*Nspec)

  return
  end

;--------------------------------------------------------------------------------------------

; This function modifies a 2-d flux table (x by y) to take into account losses to the seeing halo 
; in AO correction.  I.e., it "bleeds" flux between pixels
; Note this does not account for any kind of atmospheric/AO
; throughput, which is done elsewhere.
; This version uses a 2-component Gaussian model
pro StrehlBleedT,PARMVEC,Data,BledData,lambda,PltScl
  FWHM_SEE=FWHMSeeing(lambda,PARMVEC)
  FWHM_DIFF=FWHMDiffraction(lambda,PARMVEC)
  STREHL=GetStrehl(lambda,PARMVEC)

; Do the Gaussian blurring
; Assume light in diffraction core is recovered signal, light in seeing halo
; is noise
  BledNoiseData=(1-STREHL)*filter_image(Data[*,*],fwhm_gaussian=FWHM_SEE/PltScl)
  BledSignalData=STREHL*filter_image(Data[*,*],fwhm_gaussian=FWHM_DIFF/PltScl)

  BledData=BledNoiseData+BledSignalData
  return
  end

;--------------------------------------------------------------------------------------------

; This function modifies a 2-d flux table (x by y) to take into account losses to the seeing halo 
; in AO correction.  I.e., it "bleeds" flux between pixels
; Note this does not account for any kind of atmospheric/AO
; throughput, which is done elsewhere.
; This version uses an empirical Moffat model
pro StrehlBleedM,PARMVEC,Data,BledData,PltScl
  alpha=PARMVEC[17]/1000.; arcseconds
  beta=1.6
  BledData=moffatblurIFS(Data,alpha/PltScl,beta)
  return
  end

;--------------------------------------------------------------------------------------------

; This function modifies a 3-d data cube (x by y by z) to take into account losses to the seeing halo 
; in AO correction.  I.e., it "bleeds" flux between pixels
; Note this does not account for any kind of atmospheric/AO
; throughput, which is done elsewhere.
; This version uses an empirical Moffat model
pro StrehlBleedM3d,PARMVEC,Cube,BledCube,PltScl
  alpha=PARMVEC[17]/1000.; arcseconds
  beta=1.6
  BledCube=Cube
  for i=0,(size(Cube))[3]-1 do begin
    BledCube[*,*,i]=moffatblurIFS(Cube[*,*,i],alpha/PltScl,beta)
  endfor
  return
  end

;--------------------------------------------------------------------------------------------

; Rescale data numbers from FITS map to line fluxes
pro MakeFluxMap,Data,LineStrength
; Calculate 1-sigma threshhold above the background
; note this doesn't do any fancy clipping of the source
; so it will only be an approximation that gets much
; worse for very large/bright sources.
; In this case, replace the runningsum where values
; are greater than threshhold with the simple total

  thresh=sqrt((moment(Data))[1])+(moment(Data))[0]
;  runningsum=total(Data)
  runningsum=total(Data[where(Data gt thresh)])

  Data=LineStrength*1e-17*Data/runningsum
  return
  end
;--------------------------------------------------------------------------------------------

; Rescale data numbers from FITS map to line fluxes
pro MakeFluxMap2,Data,LineStrength
; Assume no noise in specified map, so total is simply the total
;  thresh=sqrt((moment(Data))[1])+(moment(Data))[0]
  runningsum=total(Data)
;  runningsum=total(Data[where(Data gt thresh)])

  Data=LineStrength*1e-17*Data/runningsum
  return
  end
;--------------------------------------------------------------------------------------------

; Figure out roughly what the OSIRIS field of view is (in arcsec)
pro GetFOV,Filter,PltScl,xfov,yfov
  IF(Filter EQ 'BB') THEN BEGIN
    xfov=16.*PltScl
    yfov=64.*PltScl
  ENDIF ELSE BEGIN
    xfov=48.*PltScl
    yfov=64.*PltScl
  ENDELSE
  RETURN
  END

;--------------------------------------------------------------------------------------------

; Global parameter subroutine- reads basic globals from parameter file
pro SetParms,PARMVEC,SKYFILE,TRANSFILE
  openr, lun, '/Users/drlaw/GenCode/IDL/IFsim/IFsim_v2/IFparam_Keck', /GET_LUN
  readf, lun, temp;DIAM_TEL
  PARMVEC[0]=temp
  readf, lun, temp;RES_SPEC
  PARMVEC[1]=temp
  readf, lun, temp;N1
  PARMVEC[2]=temp
  readf, lun, temp;N2
  PARMVEC[3]=temp
  readf, lun, temp;RDK
  PARMVEC[4]=temp
  readf, lun, temp;RN
  PARMVEC[5]=temp
  readf, lun, temp;T0
  PARMVEC[6]=temp
  readf, lun, temp;TEMP_AO
  PARMVEC[7]=temp
  readf, lun, temp;TEMP_TEL
  PARMVEC[8]=temp
  readf, lun, temp;ETA_AO
  PARMVEC[9]=temp
  readf, lun, temp;ETA_PSYS
  PARMVEC[10]=temp
  readf, lun, temp;ETA_TEL
  PARMVEC[11]=temp
  readf, lun, temp;SEE_V
  PARMVEC[12]=temp
  readf, lun, temp;WFE
  PARMVEC[13]=temp
  readf, lun, SKYFILE
  readf, lun, TRANSFILE
  readf, lun, temp;ETA_FILT
  PARMVEC[14]=temp
  readf, lun, temp;EPS_AO
  PARMVEC[15]=temp
  readf, lun, temp;EPS_TEL
  PARMVEC[16]=temp
  readf, lun, temp;MALPHA
  PARMVEC[17]=temp
  close,lun
  free_lun,lun
  return
end

;--------------------------------------------------------------------------------------------


function MakeModelCube,PARMVEC,VelModel,BledData,CubeParms,SKYFILE,TRANSFILE

  c=299790.;speed of light (km/s)
  xsize=(size(VelModel))[1]
  ysize=(size(VelModel))[2]
  zsize=400; Fix this for simplicity
  Lam0=CubeParms[3]-zsize/2.*CubeParms[0]
;  zsize=FIX(2*(LamCenter-Lam0)/PixScl);Ensure that integer

  PS=PixelScale(CubeParms[3]); Pixel scale (Angstroms/pixel)
  Nspec=CubeParms[2]/PS; Number of pixels summed in spectral direction

  ldisppix=CubeParms[2]/CubeParms[0]/2.35; Sigma of the line in spec channels

  MaxVelRange=(CubeParms[3]-Lam0)/Lam0*c

  Cube=fltarr(xsize,ysize,zsize)

  for i=0,xsize-1 do begin
    for j=0,ysize-1 do begin
      if ((BledData[i,j] ne 0) and (VelModel[i,j] gt -MaxVelRange)) then begin
        zloc=(VelModel[i,j]/c*CubeParms[3]+CubeParms[3]-Lam0)/CubeParms[0]
        ; Make a normalized gaussian at the requisite location
        spec=psf_gaussian([1.,zloc,ldisppix],npixel=zsize,/NORMAL)

;        spec=psf_gaussian(NDIMEN=1,NPIXEL=zsize,CENTROID=zloc,FWHM=ldisppix,/NORMAL)
        ; Scale the values so that it is
        ; normalized to the input flux for the pixel
        ; Note that need to xply BledData by Nspec to get proper units
        for k=0,zsize-1 do begin
          Cube[i,j,k]=spec[k]*BledData[i,j]*Nspec
        endfor
      endif
    endfor
  endfor

  Rmap=fltarr(xsize,ysize)
  FluxMap=fltarr(xsize,ysize)
  Ntotmap=fltarr(xsize,ysize)
  Nskymap=fltarr(xsize,ysize)
  Sigmap=fltarr(xsize,ysize)
  Specmap=fltarr(xsize,ysize)
  SNmap=fltarr(xsize,ysize)
;  systime_SEED=-1.

  FinalCube=fltarr(xsize,ysize,zsize)

  ; Now make it realistic by converting to data numbers and adding noise
  for k=0,zsize-1 do begin
    lamtemp=Lam0+k*CubeParms[0]
;    convfactor=LineFluxToEFlux(lamtemp,LineWidth,DIAM_TEL,N1,N2,ETA_AO,ETA_TEL,ETA_PSYS,TRANSFILE)
    FluxMap=Cube[*,*,k]
    Rmap=FluxMap
    Rbg=BackgroundFlux(PARMVEC,lamtemp,CubeParms[2],CubeParms[1],0,SKYFILE,TRANSFILE)
    Ntotmap=PARMVEC[2]*CubeParms[4]*(Rmap+Rbg+PARMVEC[4]) + PARMVEC[5]*(CubeParms[4]/PARMVEC[6])*PARMVEC[2]

    Nskymap=PARMVEC[2]*CubeParms[4]*(Rbg+PARMVEC[4]) + PARMVEC[5]*(CubeParms[4]/PARMVEC[6])*PARMVEC[2]
    Sigmap=sqrt((Rmap+Rbg+PARMVEC[4])*PARMVEC[2]*CubeParms[4]+PARMVEC[5]*PARMVEC[5]*(CubeParms[4]/PARMVEC[6])*PARMVEC[2])

    Specmap=RANDOMN(systime_SEED,xsize,ysize)*Sigmap+Ntotmap-Nskymap
    FinalCube[*,*,k]=Specmap/CubeParms[4]
  endfor

  return,FinalCube
end

;--------------------------------------------------------------------------------------------


function MakeModelCube2,PARMVEC,VelModel,SigModel,Data,CubeParms,SKYFILE,TRANSFILE

  c=299790.;speed of light (km/s)
  xsize=(size(VelModel))[1]
  ysize=(size(VelModel))[2]
  zsize=400; Fix this for simplicity
  Lam0=CubeParms[3]-zsize/2.*CubeParms[0]
;  zsize=FIX(2*(LamCenter-Lam0)/PixScl);Ensure that integer

  PS=PixelScale(CubeParms[3]); Pixel scale (Angstroms/pixel)
  Nspec=CubeParms[2]/PS; Number of pixels summed in spectral direction

  MaxVelRange=(CubeParms[3]-Lam0)/Lam0*c

  Cube=fltarr(xsize,ysize,zsize)

; Make gaussian emission lines for each input pixel
  for i=0,xsize-1 do begin
    for j=0,ysize-1 do begin
      if ((Data[i,j] ne 0) and (VelModel[i,j] gt -MaxVelRange) and (SigModel[i,j] gt 0.)) then begin
        zloc=(VelModel[i,j]/c*CubeParms[3]+CubeParms[3]-Lam0)/CubeParms[0]
        ldisppix=SigModel[i,j]/CubeParms[0]/2.35; Sigma of the line in spec channels
        ; Make a normalized gaussian at the requisite location
        spec=psf_gaussian([1.,zloc,ldisppix],npixel=zsize,/NORMAL)

;        spec=psf_gaussian(NDIMEN=1,NPIXEL=zsize,CENTROID=zloc,FWHM=ldisppix,/NORMAL)
        ; Scale the values so that it is
        ; normalized to the input flux for the pixel
        ; Note that need to xply BledData by Nspec to get proper units
        ; This simply takes out a factor introduced earlier, it's
        ; ok that it isn't exactly the right value
        for k=0,zsize-1 do begin
          Cube[i,j,k]=spec[k]*Data[i,j]*Nspec
        endfor
      endif
    endfor
  endfor

; Now rebin to real size
  ; these sizes are the rebinned sizes
  x1size=xsize*CubeParms[5]/CubeParms[1]
  y1size=ysize*CubeParms[5]/CubeParms[1]
  CubeRebin=fltarr(x1size,y1size,zsize)
  print,'cuberebin: ',size(CubeRebin)
  for i=0,zsize-1 do begin
    CubeRebin[*,*,i]=frebin(Cube[*,*,i],x1size,y1size,/total)
  endfor
 
; Steps to crop out only the OSIRIS FOV were coded below, but
; comment them out for the general release code.
; Replace with simple copy-over.
  x2size=x1size
  y2size=y1size
  CubeRebinFull=CubeRebin

; Get FOV (arcsec) from subroutine
;  Filter='Kn3'
;  GetFOV,Filter,CubeParms[0],xfov,yfov
; and make the actual FOV this size
;  x2size=xfov/CubeParms[0]
;  y2size=yfov/CubeParms[0]
;  CubeRebinFull=fltarr(x2size,y2size,zsize)
;  xmid=x2size/2
;  ymid=y2size/2
;  for i=0,x1size-1 do begin
;    for j=0,y1size-1 do begin
;       for k=0,zsize-1 do begin
;          CubeRebinFull[xmid+i,ymid+j,k]=CubeRebin[i,j,k]
;       endfor
;    endfor
; endfor

; Do the strehl bleeding due to seeing blurring plus AO corrected core   
  StrehlBleedM3d,PARMVEC,CubeRebinFull,BledCube,CubeParms[1]
;  BledCube=Cube

  Rmap=fltarr(x2size,y2size)
  FluxMap=fltarr(x2size,y2size)
  Ntotmap=fltarr(x2size,y2size)
  Nskymap=fltarr(x2size,y2size)
  Sigmap=fltarr(x2size,y2size)
  Specmap=fltarr(x2size,y2size)
  SNmap=fltarr(x2size,y2size)
;  systime_SEED=-1.

print,'xsize =',xsize
print,'xsize2=',x2size

  FinalCube=fltarr(x2size,y2size,zsize)

  ; Now make it realistic by converting to data numbers and adding noise
  for k=0,zsize-1 do begin
    lamtemp=Lam0+k*CubeParms[0]
;    convfactor=LineFluxToEFlux(lamtemp,LineWidth,DIAM_TEL,N1,N2,ETA_AO,ETA_TEL,ETA_PSYS,TRANSFILE)
    FluxMap=BledCube[*,*,k]
    Rmap=FluxMap
    Rbg=BackgroundFlux(PARMVEC,lamtemp,CubeParms[2],CubeParms[1],0,SKYFILE,TRANSFILE)
    Ntotmap=PARMVEC[2]*CubeParms[4]*(Rmap+Rbg+PARMVEC[4]) + PARMVEC[5]*(CubeParms[4]/PARMVEC[6])*PARMVEC[2]

    Nskymap=PARMVEC[2]*CubeParms[4]*(Rbg+PARMVEC[4]) + PARMVEC[5]*(CubeParms[4]/PARMVEC[6])*PARMVEC[2]
    Sigmap=sqrt((Rmap+Rbg+PARMVEC[4])*PARMVEC[2]*CubeParms[4]+PARMVEC[5]*PARMVEC[5]*(CubeParms[4]/PARMVEC[6])*PARMVEC[2])

    Specmap=RANDOMN(systime_SEED,x2size,y2size)*Sigmap+Ntotmap-Nskymap

    FinalCube[*,*,k]=Specmap/CubeParms[4]
  endfor

  return,FinalCube
  end

;--------------------------------------------------------------------------------------------

; Task1 subroutine
; Do all of the stuff associated with Task1=Making S/N map and cube
pro Task1,PARMVEC,SKYFILE,TRANSFILE
  print,''
  print,'This routine requires an input FITS map'
  print,'of the source, which will be cropped, scaled by an'
  print,'input line flux and artificially observed.'
  print,''
; Specify an input FITS file

  InpFile=''
  SNFile=''
  CubeFile=''
  InpFile = dialog_pickfile( title='Input FITS image', $
                                     filter='*.fits', $
                                     get_path=new_path, $
                                     /MUST_EXIST)

  print,'Selected input file: ', InpFile
;  InpFile = '/home/varda/drlaw/Research/HiZ/GOODSN/ACS/Stacked/acs_sum2.fits'
; What rotation angle to use (clockwise rotation)
  read,InpRot,prompt='Input rotation angle: '
; And where the object is
  read,Xobject,prompt='Target X coordinate: '
  read,Yobject,prompt='Target Y coordinate: '
  read,InpPixScale,prompt='Input pixel scale (arcsec/pixel): '

; Specify OSIRIS filter and lenslet size
  Filter=''
  read,Filter,prompt='Broadband (BB) or narrow (NB): '
  read,PltScl,prompt='Lenslet scale (arcsec): '
; Get FOV (arcsec) from subroutine
  GetFOV,Filter,PltScl,xfov,yfov

; Read in data from file, rotating and rebinning
  InputData,InpFile,InpPixScale,InpRot,Xobject,Yobject,PltScl,Filter,Data
;  writefits,'InpDataTest.fits',Data

; Get spectral line information and observing information
  read,RestWave,prompt='Spectral line rest wavelength (Angstroms): '
  print,'Rest frame line FWHM (km/s):'
  read,LineWidth,prompt='Enter -1 for unresolved: '
  read,Redshift,prompt='Redshift: '
  lambda=RestWave*(1.+Redshift)
  PixScl=PixelScale(lambda)
  ; Convert line width to Angstroms
  LineWidth=LineWidth/300000.*lambda;
  IF (LineWidth LE 0) THEN BEGIN
    LineWidth = lambda/PARMVEC[1]
  ENDIF
  print,'Line Width: ',LineWidth, ' Angstroms'
  read,LineStrength,prompt='Line Strength (units of e-17 erg/s/cm^2): '
  read,xtime,prompt='Exposure time (s): '
  print,'2-component (T) or Moffat (M) PSF?'
  psftype=''
  read,psftype,prompt='(M more accurate for high-strehl cases): '

; Rescale the data numbers from FITS map to line fluxes
  MakeFluxMap,Data,LineStrength

; Do the strehl bleeding due to seeing blurring plus AO corrected core
  if (psftype ne 'M') then begin
    print,'Using 2-component PSF model'
    StrehlBleedT,PARMVEC,Data,BledData,lambda,PltScl
  endif else begin
    print,'Using Moffat PSF model'
    StrehlBleedM,PARMVEC,Data,BledData,PltScl
  endelse

  ; Calculate scale for energy flux to e- flux, going through all throughputs
  ; Note that e- flux is per pixel (both spectrally and spatially)
  EnScale=LineFluxToEFlux(PARMVEC,lambda,LineWidth,TRANSFILE)
  ; Scale by this value
  BledData=BledData*EnScale

; Make the S/N map
  SNmap=MakeSNmap(PARMVEC,BledData,lambda,LineWidth,PltScl,SKYFILE,TRANSFILE,xtime)

; Make the data cube
  VelModel=Data
  VelModel[*,*]=0.
  CubeParms=dblarr(5)
  CubeParms[0]=PixScl
  CubeParms[1]=PltScl
  CubeParms[2]=LineWidth
  CubeParms[3]=lambda
  CubeParms[4]=xtime
  Cube=MakeModelCube(PARMVEC,VelModel,BledData,CubeParms,SKYFILE,TRANSFILE)
  ; Make a smoothed data cube analagous to my OSIRIS reduction method
  SmthCube=Cube
  for k=0,(size(Cube))[3]-1 do begin
    SmthCube[*,*,k]=filter_image(Cube[*,*,k],fwhm_gaussian=1.6)
  endfor

; Define output files
  SNFile = dialog_pickfile( title='Output S/N ratio map', $
                                     filter='*.fits', $
                                     /write)
  CubeFile = dialog_pickfile( title='Output data cube', $
                                     filter='*.fits', $
                                     /write)
  CubeSmoothFile = dialog_pickfile( title='Output smoothed data cube', $
                                     filter='*.fits', $
                                     /write)

  print,'Output SN file: ', SNFile
  print,'Output Data cube: ',CubeFile
  print,'Output Smoothed Data cube: ',CubeSmoothFile

; Transpose output cube to be lambda,y,x like for OSIRIS data cubes
  Cube=transpose(Cube,[2,1,0])
  SmthCube=transpose(SmthCube,[2,1,0])

; Make a minimal header
  mkhdr, head, Cube
; add wavelength info
  fxaddpar, head, 'CRVAL1', (lambda-200.*PixScl)/10.
  fxaddpar, head, 'CUNIT1', 'nm'
  fxaddpar, head, 'CDELT1', PixScl/10.
  fxaddpar, head, 'CRPIX1', 1

; Write the output files
  if ((SNFile ne "") and (CubeFile ne "") and (CubeSmoothFile ne "")) then begin
    writefits,SNFile,SNmap
    writefits,CubeFile,Cube,head
    writefits,CubeSmoothFile,SmthCube,head
  endif
  return
  end
  
;--------------------------------------------------------------------------------------------

;Task2 subroutine for release
; Do all the stuff for simulating a flux map + velocity field
; superimposed to make an observed data cube
; Note that things in this routine are done is a very different order
; than for Task1!!!!!!
pro Task2,PARMVEC,SKYFILE,TRANSFILE
; Require 3 input FITS files.  They *MUST* match exactly in dimension
; InpFlux: contains pixel flux values
; InpVel: contains relative velocity shift values
; InpSig: contains velocity dispersion values
; These are for the intrinsic source, just binned into a finite grid.
; Neglect cosmology: relative velocity, dispersion are in the source frame.

  InpFluxFile=''
  InpVelFile=''
  InpSigFile=''
  InpFluxFile = dialog_pickfile( title='Input FITS Flux image', $
                                     filter='*.fits', $
                                     get_path=new_path, $
                                     /MUST_EXIST)
  InpVelFile = dialog_pickfile( title='Input FITS Velocity image', $
                                     filter='*.fits', $
                                     get_path=new_path, $
                                     /MUST_EXIST)
  InpSigFile = dialog_pickfile( title='Input FITS Sigma image', $
                                     filter='*.fits', $
                                     get_path=new_path, $
                                     /MUST_EXIST)

; There are too many possible uses to code a general-purpose rotating
; and rebinning algorithm that works correctly in all cases.
; Therefore, the input files will not be rotated, and are assumed
; to actually fit in the OSIRIS FOV.  They will be rebinned by
; pixel scale, but not cropped.
; Responsibility is on the user to check whether it would actually fit.

  read,InpPixScale,prompt='Input pixel scale (arcsec/pixel): '

  ; read in data
  InpSN0=readfits(InpFluxFile)
  InpVel0=readfits(InpVelFile)
  InpSig0=readfits(InpSigFile)

; Calculate what a 'typical' sigma for the map is for
; noise calculation purposes.
  LineWidth=median(InpSig0[where(InpSig0 ne 0.)])


; Specify OSIRIS filter and lenslet size
  Filter=''
  read,Filter,prompt='Broadband (BB) or narrow (NB): '
  read,PltScl,prompt='Lenslet scale (arcsec): '

; Get spectral line information and observing information
  read,RestWave,prompt='Spectral line rest wavelength (Angstroms): '
  read,Redshift,prompt='Redshift: '
  lambda=RestWave*(1.+Redshift)
  PixScl=PixelScale(lambda)

  ; Convert line width to Angstroms in observed frame
  ; and convolve with the instrumental line width
  LineWidth=LineWidth/300000.*lambda
  LineWidth=sqrt(LineWidth*LineWidth + lambda*lambda/PARMVEC[1]/PARMVEC[1])
  InpSig0=sqrt(InpSig0*InpSig0 +  lambda*lambda/PARMVEC[1]/PARMVEC[1])

  ; Read scaling factors: line strength, exposure time, PSF type
  read,LineStrength,prompt='Line Strength (units of e-17 erg/s/cm^2): '
  read,xtime,prompt='Exposure time (s): '
  print,'2-component (T) or Moffat (M) PSF?'
  psftype=''
  read,psftype,prompt='(M more accurate for high-strehl cases): '

  SNFile=''
  CubeFile=''

;;;;;;;;;
  
; Rescale the data numbers from FITS map to line fluxes
  MakeFluxMap2,InpSN0,LineStrength
;  writefits,'temp4.fits',InpSN0
  ; Calculate scale for energy flux to e- flux, going through all throughputs
  ; Note that e- flux is per pixel (both spectrally and spatially)
  EnScale=LineFluxToEFlux(PARMVEC,lambda,LineWidth,TRANSFILE)
  ; Scale by this value
  InpSN0=InpSN0*EnScale

; Make the data cube
  CubeParms=dblarr(6)
  CubeParms[0]=PixScl
  CubeParms[1]=PltScl
  CubeParms[2]=LineWidth
  CubeParms[3]=lambda
  CubeParms[4]=xtime
  CubeParms[5]=InpPixScale
  ; Unlike first routine, the strehl bleeding in this case is done
  ; INSIDE the MakeCube routine so that it is wavelength dependent
  Cube=MakeModelCube2(PARMVEC,InpVel0,InpSig0,InpSN0,CubeParms,SKYFILE,TRANSFILE)


  ; Make a smoothed data cube analagous to my OSIRIS reduction method
  SmthCube=Cube
  for k=0,(size(Cube))[3]-1 do begin
    SmthCube[*,*,k]=filter_image(Cube[*,*,k],fwhm_gaussian=1.6)
  endfor


  CubeFile = dialog_pickfile( title='Output data cube', $
                                     filter='*.fits', $
                                     /write)
  CubeSmoothFile = dialog_pickfile( title='Output smoothed data cube', $
                                     filter='*.fits', $
                                     /write)

; Transpose output cube to be lambda,y,x like for OSIRIS data cubes
  Cube=transpose(Cube,[2,1,0])
  SmthCube=transpose(SmthCube,[2,1,0])

; Make a minimal header
  mkhdr, head, Cube
; add wavelength info
  fxaddpar, head, 'CRVAL1', (lambda-200.*PixScl)/10.
  fxaddpar, head, 'CUNIT1', 'nm'
  fxaddpar, head, 'CDELT1', PixScl/10.
  fxaddpar, head, 'CRPIX1', 1

; Write the output files
  if ((CubeFile ne "") and (CubeSmoothFile ne "")) then begin
    writefits,CubeFile,Cube,head
    writefits,CubeSmoothFile,SmthCube,head
  endif

  return
  end
;;;;;;;;;;;;;


;--------------------------------------------------------------------------------------------
