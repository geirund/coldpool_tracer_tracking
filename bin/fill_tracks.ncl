begin
 path=path
 fio = addfile(path+"irt_objects_mask.nc","r") 
 fit = addfile(path+"irt_tracks_mask.nc","r")

 objects  = fio->var1
 tracksin = fit->var1
 printVarSummary(objects)
 dsizes = dimsizes(objects)
 exists = isfilepresent(path+"irt_tracks_mask_fulltime.nc")
 if(exists) then
   print("file there")
   del = systemfunc("rm " + path + "irt_tracks_mask_fulltime.nc")
   print(del)
 end if

 fo = addfile(path+"irt_tracks_mask_fulltime.nc","c")
       fAtt               = True            ; assign file attributes
       fAtt@title         = "precipitaion track mask without time gaps"
       fAtt@source_file   =  "irt_tracks_mask.nc"
       fAtt@Conventions   = "None"
       fAtt@creation_date = systemfunc ("date")
       dim_names = (/"time","lev","y","x" /)
       dim_sizes = (/-1,dsizes(1),dsizes(2),dsizes(3)/)
       dimUnlim = (/ True, False, False, False/)
       filedimdef(fo,dim_names,dim_sizes,dimUnlim)

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
