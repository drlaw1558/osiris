; Break a plot into segments by not plotting certain parts
; onoff is the vector saying whether or not to plot,
; 0=no, 1=yes
pro breakplot,xvec,yvec,onoff,_EXTRA=KeywordsForPlot

if ((n_elements(xvec) ne n_elements(yvec))or(n_elements(xvec) ne n_elements(onoff))or(n_elements(yvec) ne n_elements(onoff))) then begin
  print,'Wrong size inputs!'
  return
endif
  
; This will say what segment a point is in
; Values of 0 are not to be plotted, 1+ are
segment=lonarr(n_elements(xvec))

thissegment=0
started=0
for i=0,n_elements(xvec)-1 do begin
  ; Case where we're in the middle of the vector
  if ((started eq 1)) then begin
    ; Both this element and last are on
    if ((onoff[i] eq 1)and(onoff[i-1] eq 1)) then segment[i]=thissegment

    ; This element is on, last is off
    if ((onoff[i] eq 1)and(onoff[i-1] eq 0)) then begin
     thissegment+=1 ; increment the counting
     segment[i]=thissegment
   endif
    
    ; This element is off (do nothing)

 endif

  ; Case where we're finding the data start point
  if ((started eq 0)and(onoff[i] eq 1)) then begin
    started=1
    thissegment+=1
    segment[i]=thissegment
  endif
endfor

; Now loop over the segments calling plot routine
nseg=max(segment)
for i=1,nseg do begin
  indx=where(segment eq i)
  oplot,xvec[indx],yvec[indx],_EXTRA=KeywordsForPlot
endfor

return
end
