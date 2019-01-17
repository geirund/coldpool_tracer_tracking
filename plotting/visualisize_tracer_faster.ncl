begin
; 
startzeit = 45 
;--------------------------------
;--- GET DATA
;------------------------------
; get tracer heade to identify self imerged events
 header =readAsciiTable("$results/coldpool/lindp2K/input/cp/irt_tracks_header_sorted.txt",5,"float",0)
 IDheader     = header(:,0)   ; track ID given in header 
 startLabel   = header(:,3)   ; how does tarck starts ?
 indself = ind(startLabel .eq. 0) ; where does from nothing
 track0 = IDheader(indself) ; trackIDs for tracks starting from 0
; get rain file to get original cogs
 body=readAsciiTable("$results/coldpool/lindp2K/output/raincell/irt_tracks_output_pure.txt",12,"float",0)
;--- tracer data
; data_all   = readAsciiTable("cp_history.txt",10,"double",1)
 data_all   = readAsciiTable("$results/coldpool/lindp2K/output/cp/coldpool_tracer_out.txt",12,"float",0) 
 timei=0
 agei=1
 traceIDi= 2
 cpIDi = 3
 Xpos = 4
 Ypos = 5
 Xgp = 7
 Ygp = 8
 data = data_all(ind(data_all(:,cpIDi) .eq. any(track0)),:)
 print("number of tracer: " +max(data(:,traceIDi)))
 print("number of cps: "+count_unique_values(data_all(:,cpIDi)))
 print("number of cps: "+count_unique_values(data(:,cpIDi)))
; print(get_unique_values(data(:,cpIDi)))
;-------------------
;---modelo output
; filein = addfile("$modelo/test1plus4K/level1/test1plus4K.out.vol.t.nc","r")
; plotvar = filein->t(134,:,:,0)
 filein = addfile("$modelo/../UCLA_lind/lind_p2K/lind_p2K.out.vol.u.nc","r")
 dims = getfilevardimsizes(filein,"u")
; u = filein->u(startzeit-1::,:,:,0)
 filein = addfile("$modelo/../UCLA_lind/lind_p2K//lind_p2K.out.vol.v.nc","r")
; v = filein->v(startzeit-1::,:,:,0)
; filein=addfile("irt_tracks_mask.nc","r")
; plotvar = filein->var1

  resv =True
  resv@gsnFrame     = False
  resv@gsnDraw      = False
;  resv@gsnMaximize  = True     ; Maximize plot in frame.
;  resv@trGridType = "TriangularMesh"  ;makes it faster
;  resv@vcVectorDrawOrder       = "PostDraw"        ; draw vectors last
;  resv@vpClipOn     = False    ; vector over plot border
  resv@vcMinDistanceF             = 0.01
  resv@vcRefLengthF               = 0.05
  resv@vcRefMagnitudeF         = 5.0              ; make vectors larger
  resv@vcGlyphStyle            = "CurlyVector"    ; turn on curly vectors;   
  resv@vcPositionMode ="ArrowCenter"; "ArrowHead";Default: ArrowCenter (/"ArrowCenter","ArrowTail","ArrowHead"/)
  resv@tiMainString = "Vector plot with default clipping"
  resv@vcRefAnnoJust = 0
  resv@vcRefAnnoSide ="Right"
  resv@vcRefAnnoParallelPosF = 1-0.6
  resv@vcRefAnnoOrthogonalPosF = -0.3
  resv@trXMaxF =1024 ;320 ; 150 ;170 ;150; 200;170; ;200
  resv@trYMaxF =1024; 320; 150 ;240 ;210;  250;160 ;250
  resv@trXMinF =0 ;150 ;70  ;130;  100;100;; 100
  resv@trYMinF =0; 150 ;140 ;190;  150;90 ; 150

; res = True
; res@gsnDraw = False
; res@gsnFrame = False
; res@cnFillPalette        = "MPL_Pastel2"
; res@cnFillOn             = True       ; turn on color fill
; res@cnLinesOn           = False
; res@cnLineLabelsOn       = False      ; turns off contour line labels
; res@cnInfoLabelOn        = False     
; res@cnFillMode          = "RasterFill"
; res@cnRasterSmoothingOn = True
; res@cnRasterSampleFactorF = 0.9 ; speed 
; pl = gsn_csm_contour(wks,plotvar,res)     ; create plotnd
 

;--- select CP Tracer
;find first CPid:
 cpIDs = get_unique_values(data(:,cpIDi))      ; individual CP numbers 
 noCPs = count_unique_values(data(:,cpIDi))    ; number od CPs 
 notracer = count_unique_values(data(:,traceIDi))

 colors= (/"orange","darkgreen","magenta","blue","red","seagreen"/)
 colortables=(/"MPL_Blues","MPL_Greens","MPL_Purples","MPL_Oranges","MPL_Reds","MPL_Greys","MPL_PuRd"/)
 marker=new((/noCPs/),graphic)
 markers=new((/noCPs/),graphic)
 text = new((/noCPs/),graphic)
  nmarker = 36
  xleg = new(nmarker,float)
  xleg(0::4) = fspan(0.07,0.87,9) 
  xleg(1::4) = fspan(0.07,0.87,9)    
  xleg(2::4) = fspan(0.07,0.87,9)
  xleg(3::4) = fspan(0.07,0.87,9)

  xtxt = xleg+ 0.02 
  yleg = new(nmarker,float)
  yleg(0::4)=0.12
  yleg(1::4)=0.09
  yleg(2::4)=0.06
  yleg(3::4)=0.03

  ytxt=yleg 

 legmark = new(nmarker,graphic)
 legtext = new(nmarker,graphic)
 not =new(nmarker,integer)
 lres = True
 lres@gsLineThicknessF = 3.
 mres = True
 mres@gsEdgeColor  = "Black"
 cmap =read_colormap_file("default");uniform");BlAqGrYeOrReVi200")
 ii = new(nmarker,integer)
print("start main loop")
do t = 0,dim(0) ;sizes(u(:,0,0))-2
  
  u := filein->u(startzeit-1+t::,:,:,0)
  v := filein->v(startzeit-1+t::,:,:,0)
  wks = gsn_open_wks("png","plots/testtracer_new_wholedomain_stop"+sprinti("%0.4i",t) )

  pl = gsn_csm_vector(wks, u(:,:),v(:,:), resv)

  if (.not. all(ismissing(ind(data(:,0) .eq. t+startzeit))))
  datatime = data(ind(data(:,0) .eq. t+startzeit),:) 
  c = 0
  do i = 0,noCPs-1 ; loop troughy cps 
;   print(i)
   indcpID = ind(datatime(:,cpIDi) .eq. cpIDs(i))
   if (.not. all( ismissing(indcpID)))
     if (c .lt. nmarker-1)
       ii(c) = i
       not(c) = count_unique_values(datatime(indcpID,traceIDi))
       c = c+1
     end if
     mres@gsMarkerColor =  cmap(mod(i,30),:)
     xpt =  datatime(indcpID,Xpos)
     ypt=  datatime(indcpID,Ypos)
     mres@gsMarkerIndex = 16;j
     mres@gsMarkerSizeF = 4
     marker(i) = gsn_add_polymarker(wks,pl,xpt,ypt,mres) 
     delete([/xpt,ypt/])
     bodycpindex = min(ind(body(:,0) .eq. cpIDs(i)))
     cogx = body(bodycpindex,8)
     cogy = body(bodycpindex,9)

     mkres               = True
     mkres@gsMarkerIndex = 16
     mres@gsMarkerSizeF = 10 
     mkres@gsMarkerColor = cmap(mod(i,30),:)
   ;  markers(i) = gsn_add_polymarker(wks,pl,cogx,cogy,mkres)
     mkres@txFontColor = cmap(mod(i,30),:)
     mkres@txFontHeightF =0.02 
     ;mkres@txFontThicknessF = 0.5
     markers(i) = gsn_add_text(wks,pl,cpIDs(i),cogx,cogy,mkres)
   end if
   delete([/indcpID/])
;print("5")
  end do
  delete(datatime)
  draw(pl)

  txres               = True
  txres@txJust        = "CenterLeft"
  txres@txFontHeightF = 0.01
;     txres@txFontColor = cmap(mod(i,30),:)
;     text(i) = gsn_add_text(wks,pl,"  " + cpIDs(i),0.3,0.1,txres)
  gsres = True
  gsres@gsMarkerIndex = 16
  gsres@gsMarkerThicknessF = 0.7;*(i+1)
  
  do i = 0,nmarker-1
print(i)
print(ii(i))
     if (.not. ismissing(ii(i))) then
       gsres@gsMarkerColor      = cmap(mod(ii(i),30),:) 
       gsn_polymarker_ndc(wks,          xleg(i),yleg(i),gsres)
       gsn_text_ndc(wks,cpIDs(ii(i))+" ("+not(i)+")",xtxt(i),ytxt(i),txres)
     end if
  end do
  
  frame(wks)
  end if
end do
end
