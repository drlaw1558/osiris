; Very simplistic, for EAM only
; Grow the mask by 1 pixel in all directions

function ml_growmask,mask,radius

newmask=mask

index=where(mask ne 0,nindex)

for i=0L,nindex-1 do begin
  here=array_indices(mask,index[i])
  thisx=here[0]
  thisy=here[1]

  xmin=(thisx-radius) > 0
  xmax=(thisx+radius) < ((size(mask))[1]-1)
  ymin=(thisy-radius) > 0
  ymax=(thisy+radius) < ((size(mask))[2]-1)

  newmask[xmin:xmax,ymin:ymax]=1
endfor

return,newmask
end
