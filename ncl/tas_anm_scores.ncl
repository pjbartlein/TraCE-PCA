; time series plots of SVD scores

begin

plotfile = "scores_tas_anm"

datapath = "../data/SVD/"
scoresfile   = "TraCE_monthly_RSpectra_tas_anm_scores.csv"

; colors    
linecols1 = (/ "darkgreen" /)
; linecols2 = (/ "black", "red" /)

; labels
label01 = "SVD Scores 2 m Air Temp. Anomalies"
labpos01 = (/0.10, 0.94/)
; label02 = "TraCE-21ka"
; labpos02 = (/0.10, 0.955/)

; Y-axis info: single-height drawn axis
ymid = 0.0
yhalf = 4.0

ytickmin = -4.0
ytickmax = 4.0
ytickinc = 2.0
ytickminor = 1

; viewport dimensions and locations
vph   = 0.10        ; y-axis height as drawn
vphx2 = 2.0 * vph   ; double-size plot height
vpw   = 0.30        ; x-axis width
vpx   = 0.10        ; top-left origin for plots
vpy   = 0.95        ;
; delta = 0.0675      ; offset between plots
delta = 0.08    ; offset between plots

; number of plots
nplots = 4

; plot labels
monlab = (/"Component 1", "Component 2", "Component 3", "Component 4" /)

; plot device 
; wks_type = "newpdf"
; wks_type@wkOrientation = "portrait"; portrait ;landscape
wks_type = "png"
wks_type@wkWidth = 2000
wks_type@wkHeight = 2000

wks = gsn_open_wks (wks_type, plotfile)

; drawNDCGrid(wks)               ; Draw an NDC grid for reference.

; no changes below here....
;-------------------------------------------------------------------------------

; y-axis info:  double-height plotting axis
ymin = ymid - yhalf
ymax = ymid + yhalf
yminx2 = ymid - 2.0 * yhalf ; double-size plot
ymaxx2 = ymid + 2.0 * yhalf ; double-size plot

; X-axis info
begyr = 22.0
endyr = -0.039
xmax = 22.25
xmin = -0.25

xtickmin = 0
xtickmax = 22
xtickinc = 2
xtickminor = 1

; read csv files
scores = readAsciiTable(datapath + scoresfile, 11, "float", 1)
printVarSummary(scores)

; flip signs ==========================================
scores(:,1) = -1.0 * scores(:,1)
; scores(:,3) = -1.0 * scores(:,3)

; get ageBP
ageBP = -1.0 * scores(:,0)
printVarSummary(ageBP)
printMinMax(ageBP, 0)

plotlines1 = new(4, graphic)
; plotlines2 = new(12, graphic)
plotaxes = new(4, graphic)

; set resources common to all 2x plots
;-------------------------------------------------------------------------------
resPlt                          = True
resPlt@gsnFrame                 = False

resPlt@trYMinF                  = yminx2
resPlt@trYMaxF                  = ymaxx2
resPlt@trXMinF                  = xmin
resPlt@trXMaxF                  = xmax

resPlt@vpWidthF                 = vpw        
resPlt@vpHeightF                = vphx2      ; double-size plot
resPlt@vpXF                     = vpx        ; x-origin      

resPlt@xyLineThicknesses        = (/ 0.5, 0.5 /)
resPlt@xyLineColors             = linecols1
resPlt@xyDashPattern            = 0          ; Make curves all solid
resPlt@xyMarkLineModes          = "Lines"
resPlt@xyCurveDrawOrder         = "PostDraw"

resPlt@tmBorderThicknessF       = 1.0

; X-axis
resPlt@tiXAxisFont              = "Helvetica"
resPlt@tiXAxisFontColor         = "black"
resPlt@tiXAxisFontHeightF       = 0.01
resPlt@tiXAxisFontThicknessF    = 1.5

resPlt@trXReverse               = True         ; reverse the X-axis
resPlt@tmXBLabelDeltaF          = -0.5         ; decrease space b/w ticks and labels
resPlt@tmXBLabelFontHeightF     = 0.009
resPlt@tmXBMajorThicknessF      = 1.0
resPlt@tmXTBorderOn             = False
resPlt@tmXTOn                   = False

; Y-axes
resPlt@tiYAxisFont              = "Helvetica"
resPlt@tiYAxisFontColor         = "black"
resPlt@tiYAxisFontHeightF       = 0.01
resPlt@tiYAxisFontThicknessF    = 1.5

resPlt@tmYLLabelDeltaF          = -0.5         ; decrease space b/w ticks and labels
resPlt@tmYLLabelFontHeightF     = 0.009
resPlt@tmYLMajorThicknessF      = 1.0

resPlt@tmYRLabelDeltaF          = -0.5         ; decrease space b/w ticks and labels
resPlt@tmYRLabelFontHeightF     = 0.009
resPlt@tmYRMajorThicknessF      = 1.0

; turn off axes
resPlt@tmYLBorderOn             = False
resPlt@tmYLOn                   = False
resPlt@tmYRBorderOn             = False
resPlt@tmYROn                   = False

resPlt@tmXBBorderOn             = False
resPlt@tmXBOn                   = False

; plotting resources for axes
;-------------------------------------------------------------------------------
resPltAxes = resPlt

delete(resPltAxes@xyLineThicknesses)       
delete(resPltAxes@xyLineColors)            
delete(resPltAxes@xyDashPattern)           
delete(resPltAxes@xyMarkLineModes)         
delete(resPltAxes@xyCurveDrawOrder)      

resPltAxes@tmXBMode              = "Manual"
resPltAxes@tmXBFormat            = "f"
resPltAxes@tmXBTickStartF        = xtickmin
resPltAxes@tmXBTickEndF          = xtickmax
resPltAxes@tmXBTickSpacingF      = xtickinc
resPltAxes@tmXBMinorPerMajor     = xtickminor

resPlt@gsnYRefLine              = 0.0             ; create a reference line   
resPlt@gsnYRefLineThicknessF    = 0.5
resPlt@gsnYRefLineColor         = "black"

resPltAxes@tmYLMode              = "Manual"
resPltAxes@tmYLFormat            = "f"
resPltAxes@tmYLTickStartF        = ytickmin
resPltAxes@tmYLTickEndF          = ytickmax
resPltAxes@tmYLTickSpacingF      = ytickinc
resPltAxes@tmYLMinorPerMajor     = ytickminor

resPltAxes@tmYRLabelsOn          = True
resPltAxes@tmYRMode              = "Manual"
resPltAxes@tmYRFormat            = "f"
resPltAxes@tmYRTickStartF        = ytickmin
resPltAxes@tmYRTickEndF          = ytickmax
resPltAxes@tmYRTickSpacingF      = ytickinc
resPltAxes@tmYRMinorPerMajor     = ytickminor

resPltAxes@vpHeightF            = vph
resPltAxes@trYMinF              = ymin
resPltAxes@trYMaxF              = ymax

; ; resources for second set of plots
; resPlt2 = resPlt
; resPlt2@xyLineThicknesses        = (/ 1.0, 1.0 /)
; resPlt2@xyLineColors             = linecols2

; plotting resources for text and labels
;-------------------------------------------------------------------------------
resTxt1 = True
resTxt1@txFontHeightF   = 0.013
resTxt1@txFont          = 21
resTxt1@txJust          = "topleft"

resTxt2 = True
resTxt2@txFontHeightF   = 0.012
resTxt2@txFont          = 21
resTxt2@txJust          = "topleft"

resTxt3 = True
resTxt3@txFontHeightF   = 0.010
resTxt3@txFont          = 21
resTxt3@txJust          = "topleft"

; plots
;-------------------------------------------------------------------------------

; main title
gsn_text_ndc(wks,label01,labpos01(0),labpos01(1),resTxt1)
; ; subtitle
; gsn_text_ndc(wks,label02,labpos02(0),labpos02(1),resTxt2)

; ; adjusted/not adjusted label
; resTxt3@txFontColor = linecols2(1)
; gsn_text_ndc(wks,label03,labpos03(0),labpos03(1),resTxt3)
; resTxt3@txFontColor = linecols2(0)
; gsn_text_ndc(wks,label04,labpos04(0),labpos04(1),resTxt3)

; plots from bottom (November) to top (December)
do i = 0, nplots - 1
    plotnum = nplots - 1 - i

    resPlt@vpYF = vpy - (plotnum * delta)
    ; plotData(0,:) = scores(:,plotnum+1)
    ; plotData(1,:) = tsanmadj(:,plotnum+1)

    plotlines1(plotnum) = gsn_csm_xy(wks, ageBP, scores(:,plotnum+1), resPlt)

    resPltAxes@vpYF = vpy - (plotnum * delta) - (0.5 * vph)

    if (mod(plotnum, 2) .eq. 0) then
        ; left-hand axis
        resPltAxes@tiYAxisSide          = "Left"
        resPltAxes@tiYAxisString        = monlab(plotnum)
        resPltAxes@tiYAxisAngleF        = 90

        resPltAxes@tmYRBorderOn         = False
        resPltAxes@tmYROn               = False
        resPltAxes@tmYLBorderOn         = True
        resPltAxes@tmYLOn               = True
    else    
        ; right-hand axis
        resPltAxes@tiYAxisSide          = "Right"
        resPltAxes@tiYAxisString        = monlab(plotnum)
        resPltAxes@tiYAxisAngleF        = -90

        resPltAxes@tmYRBorderOn         = True
        resPltAxes@tmYROn               = True
        resPltAxes@tmYLBorderOn         = False
        resPltAxes@tmYLOn               = False
    end if

    ; plot x-axis?
    if (plotnum .eq. nplots - 1)
        resPltAxes@tiXAxisString        = "Age (ka)"
        resPltAxes@tmXBBorderOn         = True
        resPltAxes@tmXBOn               = True
    else
        resPltAxes@tiXAxisString        = ""
        resPltAxes@tmXBBorderOn         = False
        resPltAxes@tmXBOn               = False
    end if

    plotaxes(plotnum) = gsn_csm_blank_plot(wks, resPltAxes)

    ; resPlt2@vpYF = vpy - (plotnum * delta)
    ; plotData(0,:) = tsnoadj(:,plotnum+1)
    ; plotData(1,:) = scores(:,plotnum+1)
    ; plotlines2(plotnum) = gsn_csm_xy(wks, ageBP, plotData, resPlt2)

end do

frame(wks)

end
exit()