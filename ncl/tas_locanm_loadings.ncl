; map TraCE loadings

mapname = "loadings_tas_locanm_02"

ncpath = "../data/SVD/"
ncfile = "TraCE_monthly_RSpectra_tas_locanm_loadings.nc"

label1a = "SVD 2 m Air Temp. Local Anom."
masknum = 2151  ;

begin
; open file and read in data

; Bart's paths
maskdata = "/Users/bartlein/Dropbox/WorkCurrent/TraCE_diag/nc_files/mask/"\
    +"Trace21_landmask_ice-5g_v2.nc"

f00 = addfile(ncpath + ncfile + ".nc", "r")

m1 = addfile(maskdata,"r")

; landmask
lat   = m1->lat
lon   = m1->lon   ; 0 ~ 360
t     = m1->time
landmask = new ((/48,96/), integer)
landmask = m1->landmask2(masknum,:,:) ; 0 ~ 360
printVarSummary(landmask)
;landmask = lonFlip(LANDMASK) ; -180 ~ +180
printVarSummary(landmask)
landmask!0 = "lat"
landmask!1 = "lon"
landmask&lat@units = "degrees_north"
landmask&lon@units = "degrees_east"

do j = 0,95
    do k = 0,47
        if (landmask(k,j) .ge. 1) then
            landmask(k,j) = 1
        end if
    end do
end do

paleo_outline(landmask,lat,lon,1,"land")

; data
; open files and read data
f1 = addfile(ncpath + ncfile, "r")
loadings = f1->tas_locanm_loadings(:,:,:)
lon = f1->lon
lat = f1->lat
comp = f1->comp
loadings&lon       = lon
loadings&lat       = lat
loadings&comp      = comp
; loadings = lonFlip(loadings)
printVarSummary(loadings)

; loadings(1,:,:) = -1.0 * (loadings(1,:,:))
; loadings(3,:,:) = -1.0 * (loadings(3,:,:))

; loadings
; level0 = new((/21/),float);
; level0  = (/-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,\
;            0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0/) 
level0 = new((/11/),float);
level0  = (/-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0/) 
level0begclr = 73 ; 178
level0endclr = 60 ;195

; create plot
; wks_type = "newpdf"
wks_type = "png"
wks_type@wkWidth = 2048
wks_type@wkHeight = 2048
wks_type@wkOrientation = "portrait"; portrait ;landscape

wks  = gsn_open_wks (wks_type, mapname)
plot = new(4,graphic)

; colormap
gsn_define_colormap(wks,"ClimA5")

; drawNDCGrid(wks)               ; Draw an NDC grid for reference.

res = True      


; gsn high-level interfaces ---------------------------------------------------
res@gsnDraw                = False              ; don't draw
res@gsnFrame               = False              ; don't advance frame
res@gsnSpreadColors        = True               ; use total colormap
 
; contour ----------------------------------------------------------------------
;res@gsnScalarContour      = True               ; contours desired
res@mpGridLineColor        = "transparent"      ; trick ncl into drawing perimeter
res@cnInfoLabelOn          = False              ; turn off cn info label
res@cnFillOn               = True               ; color fill  
res@cnLinesOn              = False              ; no contour lines
res@cnLineLabelsOn         = False              ; no contour labels
res@cnInfoLabelOn          = False              ; no contour info label
res@cnFillMode             = "CellFill"; "RasterFill" or "CellFill"
res@cnLevelSelectionMode   = "ExplicitLevels" ;
res@cnLabelBarEndStyle     = "ExcludeOuterBoxes"

; maps -------------------------------------------------------------------------
res@mpOutlineOn            = True;     ; continental outline
res@mpFillOn               = False       ; turn off map fill
res@mpGridAndLimbOn        = True      ; Turn on lat/lon grid
res@mpPerimOn              = False      ; turn off perimeter  
res@mpProjection           = "Robinson"       ; choose map projection

res@mpGeophysicalLineColor = "black"    ; make outlines thicker
res@mpGeophysicalLineThicknessF = 4    ; make outlines thicker
res@mpDataBaseVersion      = "Ncarg4_1"     ; choose new outline database
res@mpDataSetName          = "./land" ; data base we just created
res@mpLimbLineThicknessF    = 5
res@mpGridLatSpacingF       = 30
res@mpGridLonSpacingF       = 30 
res@mpGridLineColor         = "black"
res@mpGridLineThicknessF    = 2

res@mpLimitMode             = "LatLon" ; use lat/lon coordinates to limit area
res@mpMinLatF               =  -90. 
res@mpMaxLatF               =   90. 
res@mpMinLonF               =  -180. 
res@mpMaxLonF               =   180.; common setting for all panels

; ; gsn high-level interfaces 
; res@gsnDraw                = False              ; don't draw
; res@gsnFrame               = False              ; don't advance frame
; res@gsnSpreadColors        = True               ; use total colormap

res@gsnLeftStringFontHeightF   = 0.0175
res@gsnRightStringFontHeightF  = 0.0

res@mpPerimLineThicknessF = 2.0
  
; res@mpMinLatF              =  -00. ;-60.
; res@mpMaxLatF              =  90. ;80.  
; res@mpMinLonF              =  -190 ;-60.
; res@mpMaxLonF              =   -40. ;80.

res@pmLabelBarOrthogonalPosF = -.05
res@gsnLeftStringOrthogonalPosF = -.005

; label bar *********************************
res@lbLabelBarOn        = True             ; turn off individual cb's
res@lbLabelFontHeightF  = 0.015
;res@lbLabelAutoStride   = True    ; Nice spacing for labelbar labels.
;res@lbLabelAutoStride   = True    ; Nice spacing for labelbar labels.
res@txFontHeightF = 0.025
res@txFont       = 21   

res@lbBoxLineThicknessF = 4

res@gsnSpreadColorStart     = level0begclr
res@gsnSpreadColorEnd       = level0endclr 
res@cnLevels  = level0 


; plot loadings
do i=0,3                                
  res@gsnLeftString = "~F21~ Component" + " " + tostring(i+1)
  res@gsnrightString = " "
  plot(i) = gsn_csm_contour_map(wks,loadings(i,:,:),res)  
end do

; lables
txres1     = True
txres1@txFontHeightF = 0.015
txres1@txFont        = 21

label1b = ""

txres1@txJust        = "topLeft" 
gsn_text_ndc(wks,label1a,0.05,0.82,txres1)


; Panel Plot 
resP      = True
resP@gsnPaperOrientation = "portrait"   
resP@gsnPanelLeft  = 0.05
resP@gsnPanelRight = 0.95 
resP@gsnPanelTop   = 0.9
resP@gsnPanelBottom = 0.1  
resP@txFont             = 21
resP@gsnPanelYWhiteSpacePercent = 0
resP@lbTitleOn        = False                  ; turn on title
resP@txFontHeightF = 0.0

gsn_panel(wks,plot,(/2,2/),resP)

end
exit()