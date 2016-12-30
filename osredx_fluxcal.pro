; This program flux calibrates an OSIRIS data cube
; and does telluric corrections at the same time.
; You can specify a telluric file, or leave it blank
; to ignore this.
; Requires a Vega spectrum too, default is included.
; OSIRIS input data units are ADU/s
;
; Imported to github generic version Dec 2016.
;
; calmag is REQUIRED Vega magnitude in relevant band
; calbox is a 4-ple of calibration star location
; bgbox is an optional background location
; tellbox is a 4ple of telluric location
; /notell will not try to apply a telluric correction
; /nowrite will just report the calibration values and not apply them.
; /verbose writes intermediate files
;
; David R. Law <dlaw@stsci.edu>
; Written 6/21/2009
; Last modified: 12/26/16

Forward_Function ORFC_lam, ORFC_dlam, ORFC_ABvega

; Return lambda for a given channel and filter
function ORFC_lam,Channel,filter
  temp=0.
  case Filter of
    'Hn4':begin
       temp=16520+Channel*2.0
    end
    'Hn5':begin
       temp=17210+Channel*2.0
    end
    'Kn1':begin
       temp=19550+Channel*2.5
    end
    'Kn2':begin
       temp=20360+Channel*2.5
    end
    'Kn3':begin
       temp=21210+Channel*2.5
    end
    'Kn4':begin
       temp=22080+Channel*2.5
    end
    'Kn5':begin
       temp=22920+Channel*2.5
    end
    else:begin
       print,'Error: invalid filter'
    end
 endcase

return,temp
end

; Return wavelength step for a given filter (2.0 or 2.5 ang/pix)
function ORFC_dlam,filter
  temp=0.
  case Filter of
    'Hn4':begin
      temp=2.0
    end
    'Hn5':begin
      temp=2.0
    end
    'Kn1':begin
      temp=2.5
    end
    'Kn2':begin
      temp=2.5
    end
    'Kn3':begin
      temp=2.5
    end
    'Kn4':begin
      temp=2.5
    end
    'Kn5':begin
      temp=2.5
    end
    else:begin
       print,'Error: invalid filter'
    end
 endcase

return,temp
end

; Return correction from AB to Vega magnitude
function ORFC_ABvega,filter
  temp=0.
  case Filter of
    'Hn4':begin
      temp=1.37
    end
    'Hn5':begin
      temp=1.37
    end
    'Kn1':begin
      temp=1.88
    end
    'Kn2':begin
      temp=1.88
    end
    'Kn3':begin
      temp=1.88
    end
    'Kn4':begin
      temp=1.88
    end
    'Kn5':begin
      temp=1.88
    end
    else:begin
       print,'Error: invalid filter'
    end
 endcase

return,temp
end

pro osredx_fluxcal,calmag,nowrite=nowrite,notell=notell,tellbox=tellbox,tellbgbox=tellbgbox,calbox=calbox,calbgbox=calbgbox,verbose=verbose

; Select input science files
files = dialog_pickfile( title='Read Files to Process', $
                                     filter='*rs.fits', $
                                     get_path=new_path, $
                                     /MUST_EXIST        , $
                                     /MULTIPLE_FILES)
outfiles=ml_strreplace(files,'rs.fits','rsfc.fits')
nfiles=n_elements(files)

; Input telluric files
if (~keyword_set(notell)) then $
  tellfile = dialog_pickfile( title='Telluric Star File', $
                                     filter='*.fits', $
                                     get_path=new_path, $
                                     /MUST_EXIST)
if (~keyword_set(tellbox)) then tellbox=[28,32,28,32] ;[xstart,stop,ystart,ystop]

; Input flux calibration files (either TT or telluric)
calibfile = dialog_pickfile( title='Calibration Star File', $
                                     filter='*.fits', $
                                     get_path=new_path, $
                                     /MUST_EXIST)
if (~keyword_set(calbox)) then calbox=[28,32,28,32] ;[xstart,stop,ystart,ystop]

; Read in first file to define constants
dat = readfits(files[0],head,/SILENT)
dat=transpose(dat,[2,1,0])
nspec=(size(dat))[3]
filter=sxpar(head,'SFILTER')
filter = STRCOMPRESS(filter,/REMOVE_ALL); Remove whitespace from filter name
fnu0=double(3.631e-20) ; erg/s/cm^2/Hz
c=double(2.9979e10); cm/s
hc=double(1.986e-16); Plank times speed of light (erg cm)
gain=fxpar(head,'GAIN')
if (gain eq 0.) then gain=fxpar(head,'SYSGAIN')
if (gain eq 0.) then gain=0.23 ;e-/DN
PI=!DPI
TelDiam=1000.; cm
A=PI*TelDiam*TelDiam/4. ; cm^2 Telescope area
lam0=fxpar(head,'CRVAL1')*10.; Starting lambda in angstroms
if (lam0 eq 0.) then lam0=ORFC_lam(0,filter)
dlam=fxpar(head,'CDELT1')*10.; Lambda step in angstroms
if (dlam eq 0.) then dlam=ORFC_dlam(filter)
lambda=findgen(nspec)*dlam+lam0; Wavelength vector

; Define Vega reference spectrum
vegadir=concat_dir(ml_getenv('OSIRISTOOLS_DIR'),'vegaspec/')
if ((filter eq 'Hn1')or(filter eq 'Hn2')or(filter eq 'Hn3')or(filter eq 'Hn4')or(filter eq 'Hn5')) then $
  vegaspec=readfits(vegadir+'vega_2d_h.fits')
if ((filter eq 'Kn1')or(filter eq 'Kn2')or(filter eq 'Kn3')or(filter eq 'Kn4')or(filter eq 'Kn5')) then $
  vegaspec=readfits(vegadir+'vega_2d_k.fits')

; Default telluric correction vector is 1.0 everywhere
tellspec_div=replicate(1.0,nspec)
; If doing so, derive telluric correction
if (~keyword_set(NOTELL)) then begin
  tellcube=readfits(tellfile,/silent)
  tellcube=transpose(tellcube,[2,1,0])
  tellspec=fltarr(nspec)
  for i=0,nspec-1 do tellspec[i]=total(tellcube[tellbox[0]:tellbox[1],tellbox[2]:tellbox[3],i])
  ; If passed a background region, subtract it
  if (keyword_set(tellbgbox)) then begin
    for i=0,nspec-1 do tellspec[i]=tellspec[i]-total(tellcube[tellbgbox[0]:tellbgbox[1],tellbgbox[2]:tellbgbox[3],i])
  endif
  if (keyword_set(verbose)) then writefits,'tellspec.fits',tellspec

  ; Vega spectrum
  ; Dimension 0 is lambda, 1 is flux value
  ; sampled every 0.5 Angstroms
  vegadlam=vegaspec[0,1]-vegaspec[0,0]

  ; Divide the telluric spectrum by the vega spectrum
  ; Find out which pixel j in the vega spectrum best matches i=0 in
  ; the telluric spectrum
  vegalam=vegaspec[0,*]

  vegaresamp=interpol(vegaspec[1,*],vegalam,lambda)
  tellspec_div=tellspec/vegaresamp
  tellspec_div=tellspec_div/median(tellspec_div)
  if (keyword_set(verbose)) then writefits,'tellspec_div.fits',tellspec_div
endif


; Read the calibration star
calcube=readfits(calibfile,/silent)
calcube=transpose(calcube,[2,1,0])

; Extract the calibration star spectrum
calspec=fltarr(nspec)
for i=0,nspec-1 do calspec[i]=total(calcube[calbox[0]:calbox[1],calbox[2]:calbox[3],i])
; If passed a background region, subtract it
if (keyword_set(bgbox)) then begin
  for i=0,nspec-1 do calspec[i]=calspec[i]-total(calcube[calbgbox[0]:calbgbox[1],calbgbox[2]:calbgbox[3],i])
endif
; Apply telluric correction to calibration star
calspec_orig=calspec
calspec=calspec/tellspec_div
if (keyword_set(verbose)) then writefits,'calspec.fits',calspec

; Work out correction factor needed to convert from ADU/s/channel
; to erg/s/cm2/Angstrom
corvec=gain*hc/(lambda*1e-8)/A/dlam; erg/s/cm2/Angstrom
; and apply it
calspec_cor=calspec*corvec

; Now work out what the provided magnitude is
calmag_hz=10^(-(calmag+ORFC_ABvega(filter)+48.6)/2.5);erg/s/cm2/Hz
; provided magnitude in erg/s/cm2/Ang
calmag_ang=calmag_hz*c/(lambda*1e-8)/(lambda*1e-8)/1e8

; My calibration factor is the median of the ratio between observed
; spectrum and the nominal spectrum
calfac=median(calmag_ang/calspec_cor)
print,'calfac=',calfac

; Now apply the calibration factors to the science data
if (~keyword_set(nowrite)) then begin
  for i=0,nfiles-1 do begin
    dat = readfits(files[i],head,/SILENT)
    noise = readfits(files[i], exten_no = 1,/SILENT)
    qual = readfits(files[i], exten_no = 2,/SILENT)

    dat=transpose(dat,[2,1,0])
    noise=transpose(noise,[2,1,0])
    for j=0,nspec-1 do begin
      dat[*,*,j]=dat[*,*,j]*calfac*corvec[j]
      noise[*,*,j]=noise[*,*,j]*calfac*corvec[j]
    endfor
    dat=transpose(dat,[2,1,0])*1d17; Convert to 1e-17 cgs
    noise=transpose(noise,[2,1,0])*1d17
    fxaddpar,head,'BUNIT','1e-17 erg/s/cm2/Ang'
    writefits, outfiles[i], dat,head
    writefits, outfiles[i], noise, /append
    writefits, outfiles[i], byte(qual), /append
  endfor
endif

return
end
