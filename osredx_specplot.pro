; Assumes that input is wavelengths in nm, flux in 1e-18 erg/s/cm2/Ang
; filename must end in .ps

pro osredx_specplot,filename,wave,flux,fit,sky,xrange=xrange,yrange_sci=yrange_sci,yrange_sky=yrange_sky,remove=remove

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Main science plot
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if (~keyword_set(xrange)) then xrange=[min(wave),max(wave)] 
if (~keyword_set(yrange_sci)) then yrange_sci=[min(flux),max(flux)] 

set_plot,'ps'
loadct,39
origp=!p
origx=!x
origy=!y
!p.multi=[0,1,2]
!p.thick=0 & !x.thick=0 & !y.thick=0 & !p.charthick=0
dfpsplot, filename, /color,ysize=8.,xsize=8.

plot,wave,flux,yrange=yrange_sci,xrange=xrange,/xstyle,/ystyle,xthick=4,ythick=4,thick=4,charthick=4,ytitle=textoidl('Flux (1e-18 erg s cm^{-2} \AA^{-1})'),XTICKFORMAT="(A1)",charsize=1.3,psym=10
oplot,wave,fit,color=250,thick=4

;plot,xbin_hb,hist_hb,xrange=[-50,50],/xstyle,xtitle='Delta Velocity (km/s)',ytitle='Number',charsize=1.5,xthick=6,ythick=6,thick=6,/nodata,charthick=6,yrange=[0,230],/ystyle,title='H-Beta'
;oplot,xbin_hb,hist_hb,thick=6,psym=10
;oplot,xbin_hb,yfit_hb,thick=6,color=250,psym=10
;xyouts,20,max(yfit_hb),strcompress(string(gparm_hb[1])+'+-' + string(gparm_hb[2]))

dfpsclose
; Convert to PDF and remove old PS file
spawn, strcompress('ps2pdf '+filename+' '+ml_strreplace(filename,'.ps','.pdf'))
if (keyword_set(remove)) then spawn, strcompress('rm -f '+filename)
!p.multi=0
!p=origp
!x=origx
!y=origy
cleanplot

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Sky plot
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

set_plot,'ps'
loadct,39
origp=!p
origx=!x
origy=!y
!p.multi=[0,1,2]
!p.thick=0 & !x.thick=0 & !y.thick=0 & !p.charthick=0
filename_sky=ml_strreplace(filename,'.ps','_sky.ps')
dfpsplot, filename_sky, /color,ysize=3.,xsize=8.

if (~keyword_set(yrange_sky)) then yrange_sky=[0,max(sky)]

plot,wave,sky,xrange=xrange,/xstyle,xtitle=textoidl('Wavelength (nm)'),yrange=yrange_sky,/ystyle,thick=4,xthick=4,ythick=4,charthick=4,yticks = 1, ytickname = [' ', ' '],xticklen=0.1,charsize=1.3,psym=10
;plot,xbin_hb,hist_hb,xrange=[-50,50],/xstyle,xtitle='Delta Velocity (km/s)',ytitle='Number',charsize=1.5,xthick=6,ythick=6,thick=6,/nodata,charthick=6,yrange=[0,230],/ystyle,title='H-Beta'
;oplot,xbin_hb,hist_hb,thick=6,psym=10
;oplot,xbin_hb,yfit_hb,thick=6,color=250,psym=10
;xyouts,20,max(yfit_hb),strcompress(string(gparm_hb[1])+'+-' + string(gparm_hb[2]))

dfpsclose
; Convert to PDF and remove old PS file
spawn, strcompress('ps2pdf '+filename_sky+' '+ml_strreplace(filename_sky,'.ps','.pdf'))
if (keyword_set(remove)) then spawn, strcompress('rm -f '+filename_sky)
!p.multi=0
!p=origp
!x=origx
!y=origy
cleanplot

return
end
