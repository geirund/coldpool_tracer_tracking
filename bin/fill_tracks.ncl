begin
 path=path
 fio = addfile(path+"irt_objects_mask.nc","r") 
 fit = addfile(path+"irt_tracks_mask.nc","r")

 objects  = fio->var1
 tracksin = fit->var1

 fo = addfile(path+"irt_tracks_mask_fulltime.nc","c")
       fAtt               = True            ; assign file attributes
       fAtt@title         = "precipitaion track mask without time gaps"
       fAtt@source_file   =  "irt_tracks_mask.nc"
       fAtt@Conventions   = "None"
       fAtt@creation_date = systemfunc ("date")
     filedimdef(fo,"time",-1,True)

; tt!2 ="xt"
; tt!3 ="yt"
; tt!1 ="zt"
; tt!0="time"
 dims = dimsizes(objects)
 dimstrack = dimsizes(tracksin)
 trackout = where(objects .eq. 0, 0, -1)
 nt = dims(0)
 nttrack = dimstrack(0)
 counter = 0

print(dims)
print(dimsizes(tracksin))
 do it = 0,nt-1
  if (counter .lt. nttrack)
   if any(objects(it,:,:,:) .gt. 0)
     objectmask :=where(objects(it,:,:,:) .eq. 0, 0, 1)
     trackmask :=where(tracksin(counter,:,:,:) .eq. 0, 0, 1)
     if all(objectmask .eq. trackmask) then
       print(it + " "+ counter)
       print(typeof(trackout))
       print(typeof(tracksin))
       trackout(it,:,:,:) = tointeger(tracksin(counter,:,:,:))
       counter = counter +1
     end if
   end if
   end if
 end do
 fo->trackID = trackout

end
