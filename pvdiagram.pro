@integrate_map
;+
; NAME:
;       PVDIAGRAM
; PURPOSE:
;       A GUI program for making PV diagram using FITS cube 
; EXPLANATION:
;       Use an image header to compute X and Y positions, given the
;       RA and Dec (or longitude, latitude) in decimal degrees.  
;
; CALLING SEQUENCE:
;       PVDIAGRAM,fits=fits,[vel_range=[min,max],scale=scale]
;
; INPUTS:
;       FITS   - FITS cube exported from CASA or GILDAS
;
; OPTIONAL INPUTS:
;       VEL_RANGE - velocity range for the PV diagram.
;       SCALE     - scaling factor for integrated map       
;       XYSTART   - positiion of begining point for PV analysis
;       XYEND     - positiion of ending point for PV analysis
;       PUNIT     - position unit of X-axis, 0: offset, 1: offset in RA, 2: offset in DEC
;
; OPTIONAL OUTPUT:
;
; OPERATIONAL NOTES:
;
;
; REVISION HISTORY:
;       Chi-Hung Yan            ASIAA     June, 2013      
;-




PRO GETPVDIAGRAM,CUBE=cube, header=header,XPOS=xpos, YPOS=ypos, VEL_RANGE=vel_range,pvwin_id=pvwin_id, PUNIT=punit
      ;pa=90.0
;  fitsname = fits         ;fits input file
;  
;  if keyword_set(vel_range) then begin
;     integrate_map,fits=fitsname,$
;       vel_range=vel_range,eta=1.0,outfits='/tmp/mom0.fits'
;  endif else begin
;     integrate_map,fits=fitsname,$
;       eta=1.0,outfits='/tmp/mom0.fits'
;  
;  endelse
;  
;  
;  mom0=readfits('/tmp/mom0.fits')
;  iImage,rot(mom0,pa,CUBIC=-0.5), Margin=0.01, /Axes, CTIndex=4, /Brewer, /Reverse,/fit_to_view
  pa=0
  ;sigmat = 0.2                 ;Jy/beam one sigma value of image ;1.12
  levels = 3                         ;contour levels to draw


 ;color tables and intensity ranges:
  ; For IDL version earlier than 8.0, use color table less than 40
  if float(!VERSION.release) le 8.0 then begin
     n_colTable_cont_map = 15
  endif else begin
      n_colTable_cont_map = 54 
  endelse
  ;n_colTable_cont_map = 54 ;default: 15 or 16
  ;n_colTable_cont_map = 15 ;default: 15 or 16
  intensity_min = -1       ;default: -9999 (automatic) otherwise in Jy/beam (do not use 0.0)
  intensity_max = -9999    ;default: -9999 (automatic) otherwise in Jy/beam 


  ;contour_col  = [128,128,128] 
  ;contour_col  = [0,0,0]

 ;=====================================================
  ;image=REFORM(readfits(fitsname,header))
  image=reform(cube)

  GETVAL,header,'NAXIS1',xnpix    ;reference pixel
  GETVAL,header,'NAXIS2',ynpix    ;reference pixel
  GETVAL,header,'CRPIX1',ix    ;reference pixel
  GETVAL,header,'CRPIX2',iy    ;reference pixel
  GETVAL,header,'CRPIX3',iz    ;reference pixel

  GETVAL,header,'CDELT1',dx    ;deg/pixel
  GETVAL,header,'CDELT2',dy    ;deg/pixel
  GETVAL,header,'CDELT3',dz    ;m/s/pixel
  dx_arcsec = ABS(dx*3600.0)   ;arcsec/pixel
  dy_arcsec = ABS(dy*3600.0)   ;arcsec/pixel

  GETVAL,header,'CRVAL1',ra    ;deg
  GETVAL,header,'CRVAL2',dec   ;deg
  GETVAL,header,'CRVAL3',vel0  ;m/s


  GETVAL,header,'OBSRA',ra_px
  GETVAL,header,'OBSDEC',dec_px

  GETVAL,header,'BMAJ',beam1   ;deg
  GETVAL,header,'BMIN',beam2   ;deg
  GETVAL,header,'BPA',bpa      ;deg
  beam1 = beam1*3600.0         ;arcsec
  beam2 = beam2*3600.0         ;arcsec

  GETVAL,header,'RESTFREQ',fc  ;Hz
  GETVAL,header,'VOBS',vobs    ;m/s
  funit = SXPAR(header,'BUNIT')

  ;get correct velocity axis:
  ;iz   ... reference pixel
  ;dz   ... velocity interval
  ;vel0 ... velocity at reference pixel
  ;vobs ... velocity of observatory at time of observation
  ;=> velocity at first channel:
  vel_start = vel0 - dz*(iz-1)

  PRINT, '=================================='
  ;ADXY,header,ra,dec,xc,yc,/PRINT
  
  
  xc=xnpix/2
  yc=ynpix/2

  nx=(SIZE(image))[1]  &  x=(-FINDGEN(nx)+xc)*dx_arcsec            ;arcsec
  ny=(SIZE(image))[2]  &  y=(+FINDGEN(ny)-yc)*dy_arcsec            ;arcsec
  nz=(SIZE(image))[3]  &  z=((+FINDGEN(nz))*dz+vel_start) / 1000.0 ;km/s

  ;=====================================================

  PRINT, '=================================='
  PRINT, 'GETIM: OBSRA  = ', ra,  ' deg'
  PRINT, 'GETIM: OBSDEC = ', dec, ' deg'
  PRINT, 'GETIM: BMAJ   = ', beam1, ' arcsec'
  PRINT, 'GETIM: BMIN   = ', beam2, ' arcsec'
  PRINT, 'GETIM: BPA    = ', bpa, ' deg'
  PRINT, 'GETIM: CDELT1 = ', dx_arcsec, ' arcsec/pixel'
  PRINT, 'GETIM: CDELT2 = ', dy_arcsec, ' arcsec/pixel'
  PRINT, 'GETIM: CDELT3 = ', dz, ' m/s/pixel'
  PRINT, 'GETIM: CRVAL1 = ', ra,  ' deg'
  PRINT, 'GETIM: CRVAL2 = ', dec, ' deg'
  PRINT, 'GETIM: CRVAL3 = ', vel0,' m/s'
  PRINT, 'GETIM: CRPIX1 = ', ix, '(ref pixel)'
  PRINT, 'GETIM: CRPIX2 = ', iy, '(ref pixel)'
  PRINT, 'GETIM: CRPIX3 = ', iz, '(ref pixel)'
  PRINT, 'GETIM: RESTFQ = ', fc,     ' Hz'
  PRINT, 'GETIM:        = ', fc/1e9, ' GHz'
  PRINT, 'GETIM: BUNIT  = ', funit

  IF nx NE ny THEN PRINT, 'ERROR: x and y dimension are not the same'
  IF dx_arcsec NE dy_arcsec THEN STOP
  ;IF nx NE ny THEN STOP

  ;rotate cube first (channel by channel)
  image_rot = image
  z_vel = z

  ;rotate image in clockwise direction (pos_angle is in deg)
  pos_angle=pa
  FOR i=0, nz-1 DO image_rot[*,*,i] = ROT(image[*,*,i], pos_angle, CUBIC=-0.5)
  

  ;cut along from [0,+y] to [0,-y] (basically from N to S of the image is not rotated)
  ;if the image is rotated by 10 deg then the cut corresponds to a cut along a PA of 10 deg!
  ny=n_elements(xpos)
  image_pv = MAKE_ARRAY(ny,nz,/DOUBLE)
  
  ;FOR i=0, nz-1 DO BEGIN
  ;  FOR j=0, ny-1 DO BEGIN
  ;      image_pv[j,i] = image_rot[nx/2,j,i]
  ;  ENDFOR
  ;ENDFOR
  
  FOR i=0, nz-1 DO BEGIN
    FOR j=0, n_elements(xpos)-1 DO BEGIN
        ;print,xpos[j],ypos[j]
        ;image_pv[j,i] = image_rot[fix(xpos[j]),fix(ypos[j]),i]
        image_pv[j,i] = image[floor(xpos[j]),floor(ypos[j]),i]
    ENDFOR
  ENDFOR
  
  ; Calculate the distance between two coordinate 
  xdist=(max(xpos)-min(xpos))*dx_arcsec
  ydist=(max(ypos)-min(ypos))*dy_arcsec
  dist=sqrt(xdist^2+ydist^2)
  
  ; Calculate the RA & DEC of the center
  xcenter=min(xpos)+(max(xpos)-min(xpos))*0.5
  ycenter=min(ypos)+(max(ypos)-min(ypos))*0.5
  xyad,header,xcenter,ycenter, racenter, deccenter
  radec, ra, dec, ihr, imin, xsec, ideg, imn, xsc
  
  if keyword_set(punit) then begin
    case punit of
    0: d=dist
    1: d=xdist
    2: d=ydist
    3: d=xdist
    4: d=ydist
    end
  endif else begin
    d=dist
  endelse
  y=findgen(n_elements(xpos))
  y=((y/max(y))-0.5)*d
  
  ; find out the angle of the velocity mark.
  xmax=xpos[where(xpos eq max(xpos))]
  xmin=xpos[where(xpos eq min(xpos))]
  ymax=ypos[where(xpos eq max(xpos))]
  ymin=ypos[where(xpos eq min(xpos))]
  
  rect_coord = [xmax-xmin, ymax-ymin]   
  polar_coord = CV_COORD(FROM_RECT=rect_coord, /TO_POLAR,/degrees) 
  pos_angle=90-polar_coord[0]
  if pos_angle lt 0 then pos_angle=180+pos_angle
  IF pos_angle GE 0   AND pos_angle LT  10 THEN string_pa = '00'+STRCOMPRESS(STRING(ROUND(pos_angle)),/REMOVE_ALL)
  IF pos_angle GE 10  AND pos_angle LT 100 THEN string_pa = '0' +STRCOMPRESS(STRING(ROUND(pos_angle)),/REMOVE_ALL)
  IF pos_angle GE 100 AND pos_angle LE 360 THEN string_pa =     +STRCOMPRESS(STRING(ROUND(pos_angle)),/REMOVE_ALL)


  ;=====================================================
  ind=WHERE(FINITE(image, /NAN),COMPLEMENT=inx)
  rotind=WHERE(FINITE(image_rot, /NAN),COMPLEMENT=rotinx)
  pvind=WHERE(FINITE(image_pv, /NAN),COMPLEMENT=pvinx)
  PRINT, '=================================='
;  PRINT, 'GETIM: image all max = ', MAX(image),' Jy/beam'
;  PRINT, 'GETIM: image all min = ', MIN(image),' Jy/beam'
;  PRINT, 'GETIM: image rot max = ', MAX(image_rot),' Jy/beam'
;  PRINT, 'GETIM: image rot min = ', MIN(image_rot),' Jy/beam'
;  PRINT, 'GETIM: image pv max  = ', MAX(image_pv),' Jy/beam'
;  PRINT, 'GETIM: image pv min  = ', MIN(image_pv),' Jy/beam'
  PRINT, 'GETIM: image all max = ', MAX(image[inx]),' '+funit
  PRINT, 'GETIM: image all min = ', MIN(image[inx]),' '+funit
  PRINT, 'GETIM: image all std = ', STDDEV(image[inx]),' '+funit
  PRINT, 'GETIM: image rot max = ', MAX(image_rot[rotinx]),' '+funit
  PRINT, 'GETIM: image rot min = ', MIN(image_rot[rotinx]),' '+funit
  PRINT, 'GETIM: image pv max  = ', MAX(image_pv[pvinx]),' '+funit
  PRINT, 'GETIM: image pv min  = ', MIN(image_pv[pvinx]),' '+funit
  PRINT, 'GETIM: image pv std  = ', STDDEV(image_pv[pvinx]),' '+funit
  
  if keyword_set(sigmat) then begin
    sigmat=sigmat
  endif else begin
    sigmat= STDDEV(image[inx])
  endelse
  
  
  ;setup for plotting pv contour map
  image_plot = image_pv
  ;image_plot[pvind]=-999.0
  
  yaxis_val  = z_vel
  xaxis_val  = y
  ytitle = 'Velocity (km/s)'
  xtitle = 'Position (arcsec)'
  
  if keyword_set(punit) then begin
      case punit of
        0: xtitle = 'Position (arcsec)'
        1: begin
           xrange=[max(xaxis_val),min(xaxis_val)]
           xtitle = 'Position offset in Right Ascension (arcsec)'
        end   
        2: xtitle = 'Position offset in Declination (arcsec)'
        3: begin
           xaxis_val=(xaxis_val/3600.0)+racenter
           xrange=[max(xaxis_val),min(xaxis_val)]
           xtitle = 'Right Ascension (J2000)'
        end
        4: begin
           xaxis_val=(xaxis_val/3600.0)+deccenter
           xtitle = 'Declination (J2000)'
        end    
      end
  endif else begin
      xtitle = 'Position (arcsec)'
  endelse
  ;print,y
  title = '(PV-diagram at position angle '+string_pa+'!9'+String("260B)+'!X'+')'

  IF ROUND(intensity_max) EQ -9999 THEN aZ = MAX(image[inx]) ELSE aZ = intensity_max
  IF ROUND(intensity_min) EQ -9999 THEN iZ = MIN(image[inx]) ELSE iZ = intensity_min

  FOR i=0, N_ELEMENTS(image_plot)-1 DO BEGIN
     IF image_plot[i] LT iZ THEN image_plot[i] = iZ
     IF image_plot[i] GT aZ THEN image_plot[i] = aZ
  ENDFOR

  c_number = CEIL((aZ-iZ)/(sigmat*levels))+1         ;number of levels
  dZ = c_number * (sigmat*levels)                    ;number of levels + sigma
  iZ_min = FLOOR(iZ/(sigmat*levels)) * sigmat*levels ;min level (makes sure that it goes through 0.0)
  iZ_min=0
  
  
  c_value = MAKE_ARRAY(c_number,/DOUBLE)
  rgb_indicies = MAKE_ARRAY(c_number,/INTEGER)
  
  FOR g=0, c_number-1 DO c_value[g] = iZ_min + g * 1.0/(c_number) * dZ
  IF n_colTable_cont_map EQ 15 THEN FOR g=0, c_number-1 DO rgb_indicies[g] = ROUND(80.0+g*(240-80.0)/(c_number-1.0))
  IF n_colTable_cont_map NE 15 THEN FOR g=0, c_number-1 DO rgb_indicies[g] = ROUND((g*255.0)/(c_number-1.0))


  PRINT, '=================================='
  PRINT, 'GETIM: RGB indices and contour levels '+funit
  PRINT, 'GETIM: (derived from original image)'
  FOR i=0, N_ELEMENTS(c_value)-1 DO BEGIN
    PRINT, 'GETIM: ', rgb_indicies[i], c_value[i]
  ENDFOR
  PRINT, 'GETIM: n contours = ', c_number

  ;=====================================================
  ;make pv contour map or contour map of one channel according to cplot_kind:
  ;(1 or 3 all three contour maps, 4 only the best contour map)

  resolution = 600
  
  if keyword_set(vel_range) then begin
    yrange=vel_range
  endif else begin
    yrange=[MAX(yaxis_val),MIN(yaxis_val)]
  endelse
  
  ; Find out where the level is negtive.
  ;   
  ind=where(c_value le 0, complement=inx)

  ;A: contour plot with filled contours + contour lines
  dist_c = MAX(c_value) - MIN(c_value) ;=255, for n_colTable_cont_map15: =240-80=160 of 255
        iCONTOUR, image_plot, xaxis_val, yaxis_val, XRANGE=xrange,YRANGE=yrange,$
                 ANISOTROPY=[1,1,1],XTICKFONT_SIZE=18,YTICKFONT_SIZE=18,XTICKFONT_STYLE=1,YTICKFONT_STYLE=1,XTITLE=xtitle,YTITLE=ytitle,$
                 TITLE=title,RGB_TABLE=n_colTable_cont_map,TRANSPARENCY=0,RGB_INDICES=rgb_indicies,C_VALUE=c_value,$
                 /NO_SAVEPROMPT,LOCATION=[100,0],DIMENSIONS=[1000,850],NAME=name,VIEW_ZOOM=1.0,SHADING=1,IDENTIFIER=idTool_2,$
                 xstyle=1,ystyle=1,/Fill
       
       iCONTOUR, image_plot, xaxis_val, yaxis_val, XRANGE=xrange,YRANGE=yrange,$
                 C_THICK=1.0,TRANSPARENCY=70,RGB_INDICES=rgb_indicies,C_VALUE=c_value[inx],$
                 /NO_SAVEPROMPT,SHADING=1,/OVERPLOT,c_label_show=0,RGB_TABLE=0

       ;iCONTOUR, image_plot, xaxis_val, yaxis_val, XRANGE=xrange,YRANGE=yrange,$
       ;          C_THICK=1.0,TRANSPARENCY=70,RGB_INDICES=rgb_indicies,C_VALUE=c_value[ind],$
       ;          /NO_SAVEPROMPT,SHADING=1,/OVERPLOT,c_label_show=0,C_LINESTYLE=2,RGB_TABLE=0
  IF n_colTable_cont_map NE 15 THEN iPLOT, [MAX(xaxis_val),MAX(xaxis_val)],[MIN(c_value),MAX(c_value)],/OVERPLOT,RGB_TABLE=n_colTable_cont_map,INSERT_COLORBAR=[-0.5,-0.9]
  IF n_colTable_cont_map EQ 15 THEN iPLOT, [MAX(xaxis_val),MAX(xaxis_val)],[MIN(c_value)-80./255.*dist_c,MAX(c_value)+15./255.*dist_c],/OVERPLOT,RGB_TABLE=n_colTable_cont_map,INSERT_COLORBAR=[-0.5,-0.9]
     ;  file_write = 'images/pv_'+file_out+'_plotA_pa'+string_pa+file_exten
   
  if float(!VERSION.release) le 8.0 then begin
    
  endif else begin
      ISAVE, '~/Desktop/test.tif', RESOLUTION=resolution, TARGET_IDENTIFIER=idTool_2
      spawn,'convert ~/Desktop/test.tif ~/Desktop/test.pdf'
  endelse
  
  
 END





PRO DROPLIST_EVENTS, event
    Widget_Control, event.top, Get_UValue=info, /No_Copy
    
    ; Which button caused this event?
    
    Widget_Control, event.id, Get_Value=thisButtonValue
    
    info.punit=event.index

    Widget_Control, event.top, Set_UValue=info, /No_Copy

end


PRO PVDIAGRAM_Event, event

    ; All program button events handled here.
    
    Widget_Control, event.top, Get_UValue=info, /No_Copy
    
       ; Which button caused this event?
    
    Widget_Control, event.id, Get_Value=thisButtonValue

    CASE thisButtonValue OF
    
       'Quit': BEGIN
          Widget_Control, event.top, /Destroy
          RETURN
          ENDCASE
    
       'Erase Lines': BEGIN
          WSet, info.wid
          cgimage,info.Image,ct=16
          info.x0=-1
          info.x1=-1
          info.y0=-1
          info.y1=-1
          ENDCASE
    
       'Slice': BEGIN
          if strmatch('NONE',info.itoolsid) ne 1 then idelete,info.itoolsid
          x0=info.x0/info.scale
          x1=info.x1/info.scale
          y0=info.y0/info.scale
          y1=info.y1/info.scale
          if (abs(info.x0-info.x1) ge abs(info.y0-info.y1)) then begin
              xx=getseries(x0,x1,1)
              yy=interpol([y0,y1],[x0,x1],xx)
              
          endif else begin
              yy=getseries(y0,y1,1)
              xx=interpol([x0,x1],[y0,y1],yy)
          
          endelse
          PlotS, xx*info.scale, yy*info.scale, $
               /Device, Color=info.drawColor,thick=2
          
          ; Sorting Dec array so that 
          newyy=yy(sort(yy))
          newxx=xx(sort(yy))
          getpvdiagram,cube=info.cube,header=info.hd,xpos=newxx,ypos=newyy,vel_range=info.vel_range,pvwin_id=test,punit=info.punit
          
          ; Printing the XY coordinate as a reference 
          print,'-----------------------------'
          print,'X0 = '+string(x0,format='(F7.2)')+' Y0 = '+string(y0,format='(F7.2)')
          print,'X1 = '+string(x1,format='(F7.2)')+' Y1 = '+string(y1,format='(F7.2)')
          print,'-----------------------------'
          
          ENDCASE
    
    ENDCASE
    Widget_Control, event.top, Set_UValue=info, /No_Copy
END
;------------------------------------------------------------------



PRO PVDIAGRAM_Draw_Events, event

    ; All program draw widget events handled here.
    
       ; Deal only with DOWN, UP, and MOTION events.
    
    IF event.type GT 2 THEN RETURN
    
       ; Get the info structure.
    
    Widget_Control, event.top, Get_UValue=info, /No_Copy
    
       ; What kind of event is this?
    
    eventTypes = ['DOWN', 'UP', 'MOTION']
    thisEvent = eventTypes[event.type]
    
    whichButton = ['NONE', 'LEFT', 'MIDDLE', 'NONE', 'RIGHT']
    
    CASE thisEvent OF
    
       'DOWN': BEGIN
    
             ; Which button was used? LEFT or RIGHT?
    
          info.buttonUsed = whichButton(event.press)
    
             ; Turn motion events on for the draw widget.
    
          Widget_Control, info.drawID, Draw_Motion_Events=1
    
             ; Create a pixmap. Store its ID. Copy window contents into it.
    
          Window, /Free, /Pixmap, XSize=info.xsize, YSize=info.ysize
          info.pixID = !D.Window
          Device, Copy=[0, 0, info.xsize, info.ysize, 0, 0, info.wid]
    
             ; Initialize the starting coordinates of the line.
    
          IF info.buttonUsed EQ 'LEFT' THEN BEGIN
             info.xstart = event.x
             info.ystart = event.y
             info.x0=event.x
             info.y0=event.y             
             
          ENDIF ELSE BEGIN
             info.xvalues = Ptr_New([event.x])
             info.yvalues = Ptr_New([event.y])
          ENDELSE
    
          ENDCASE
    
       'UP': BEGIN
    
             ; Erase the last line drawn. Destroy the pixmap.
    
          WSet, info.wid
          Device, Copy=[0, 0, info.xsize, info.ysize, 0, 0, info.pixID]
          WDelete, info.pixID
    
             ; Turn draw motion events off. Clear events queued for widget.
    
          Widget_Control, info.drawID, Draw_Motion_Events=0, Clear_Events=1
    
             ; Draw the final line.
    
          IF info.buttonUsed EQ 'LEFT' THEN BEGIN
    
             PlotS, [info.xstart, event.x], [info.ystart, event.y], $
               /Device, Color=info.drawColor
    
             ; Reinitialize the line starting coordinates.
    
             info.xstart = -1
             info.ystart = -1
             
             info.x1=event.x
             info.y1=event.y    
          ENDIF ELSE BEGIN
    
             PlotS, *info.xvalues, *info.yvalues, /Device, $
                Color=info.drawColor
    
             ; Reinitialize the pointers.
    
             Ptr_Free, info.xvalues
             Ptr_Free, info.yvalues
             info.xvalues = Ptr_New()
             info.yvalues = Ptr_New()
    
          ENDELSE
          ENDCASE
    
       'MOTION': BEGIN
    
             ; Here is where the actual line is drawn and erased.
             ; First, erase the last line.
    
          WSet, info.wid
          Device, Copy=[0, 0, info.xsize, info.ysize, 0, 0, info.pixID]
    
          IF info.buttonUsed EQ 'LEFT' THEN BEGIN
    
             ; Draw a straight line.
    
          PlotS, [info.xstart, event.x], [info.ystart, event.y], /Device, $
                Color=info.drawColor
    
          ENDIF ELSE BEGIN
    
             ; Get the points of the new free-hand line and draw it.
    
          *info.xvalues = [*info.xvalues, event.x]
          *info.yvalues = [*info.yvalues, event.y]
          PlotS, *info.xvalues, *info.yvalues, /Device, Color=info.drawColor
    
          ENDELSE
    
          ENDCASE
    
    ENDCASE
    
       ; Store the info structure.
    
    Widget_Control, event.top, Set_UValue=info, /No_Copy
END
;------------------------------------------------------------------



PRO PVDIAGRAM, FITS=fits, VEL_RANGE=vel_range, SCALE=scale, XYSTART=xystart, XYEND=xyend, PUNIT=punit 

   if keyword_set(vel_range) then begin
      integrate_map,fits=fits,$
         vel_range=vel_range,eta=1.0,outfits='/tmp/mom0.fits'
   endif else begin
      integrate_map,fits=fits,$
         eta=1.0,outfits='/tmp/mom0.fits'
   endelse

   if keyword_set(punit) then begin
      punit=punit
   endif else begin
      punit=0
   endelse
   
   if keyword_set(scale) then begin
      scale=scale
   endif else begin
      scale=1
   endelse
  
   if keyword_set(vel_range) then begin
     vel_range=vel_range
   endif else begin
     vel_range=[-100,100]
   endelse
   
   
   cube=readfits(fits,hd)
   mom0=readfits('/tmp/mom0.fits')


   ; Open an image data set.
   data=readfits('/tmp/mom0.fits')
   im_size=SIZE(data)
   image = CONGRID(data, scale*im_size[1], scale*im_size[2], /INTERP)
   
   xsize = (Size(image))[1]
   ysize = (Size(image))[2]

   ; Create the TLB.
   tlb = Widget_Base(Title='PV Diagram Widget Program', $
      /column)

 
   ; Create the draw widget graphics window. Turn button events ON.

   drawID = Widget_Draw(tlb, XSize=xsize, YSize=ysize, Button_Events=1, $
   Event_Pro='PVDIAGRAM_Draw_Events')


   button_base = widget_base(tlb, column=3)

   point1_button = widget_button(button_base, $
                              value = 'Slice', $
                              uvalue = 'Slice',/align_center,xoffset=40,xsize=70)
   point2_button = widget_button(button_base, $
                              value = 'Erase Lines', $
                              uvalue = 'Erase Lines',/align_center,xoffset=40,xsize=90)
   exit_button = widget_button(button_base, $
                              value = 'Quit', $
                              uvalue = 'Quit',/align_center,xsize=70)
   
   drop_base = widget_base(tlb, column=3)
   defaultValues = ['Angular Distance', 'Offset in RA', 'Offset in DEC', 'RA (J2000)', 'DEC (J2000)']
   combo =  widget_droplist(drop_base, VALUE=defaultValues, title='Position Unit',$
                               SCR_XSIZE=200, UVALUE='defaultValues',  Event_Pro='DROPLIST_EVENTS')                           

   ; Realize widgets and make draw widget the current window.
   Widget_Control, combo,SET_DROPLIST_SELECT=punit
   Widget_Control, tlb, /Realize
   Widget_Control, drawID, Get_Value=wid
   WSet, wid

   ; Load drawing color and display the image.

drawColor = !D.N_Colors-1
cgimage,Image,ct=16

if (keyword_set(xystart) and keyword_set(xyend)) then begin
  PlotS, scale*[xystart[0], xyend[0]], scale*[xystart[1], xyend[1]], /Device, $
                Color=drawColor,thick=2
  ; Create an "info" structure with information to run the program.
  info = { image:image, $        ; The image data.
           scale:scale, $
           cube:cube, $
           vel_range:vel_range, $
           hd:hd, $
           itoolsid:'NONE', $
           wid:wid, $            ; The window index number.
           drawID:drawID, $      ; The draw widget identifier.
           pixID:-1, $           ; The pixmap identifier.
           xsize:xsize, $        ; The X size of the graphics window.
           ysize:ysize, $        ; The Y size of the graphics window.
           xstart:-1, $          ; The starting X coordinate of the line.
           ystart:-1, $          ; The starting Y coordinate of the line.
           xvalues:Ptr_New(), $  ; The X coordinates of the free-hand line.
           yvalues:Ptr_New(), $  ; The Y coordinates of the free-hand line.
           buttonUsed:'NONE', $  ; A flag to indicate which button is used.
           drawColor:drawColor, $; The rubberband box color.
           punit:punit, $        ; Passing the position unit 
           x0:xystart[0]*scale, $
           x1:xyend[0]*scale, $
           y0:xystart[1]*scale, $
           y1:xyend[1]*scale }
                
                
endif else begin
  ; Create an "info" structure with information to run the program.
  info = { image:image, $        ; The image data.
           scale:scale, $
           cube:cube, $
           vel_range:vel_range, $
           hd:hd, $
           itoolsid:'NONE', $
           wid:wid, $            ; The window index number.
           drawID:drawID, $      ; The draw widget identifier.
           pixID:-1, $           ; The pixmap identifier.
           xsize:xsize, $        ; The X size of the graphics window.
           ysize:ysize, $        ; The Y size of the graphics window.
           xstart:-1, $          ; The starting X coordinate of the line.
           ystart:-1, $          ; The starting Y coordinate of the line.
           xvalues:Ptr_New(), $  ; The X coordinates of the free-hand line.
           yvalues:Ptr_New(), $  ; The Y coordinates of the free-hand line.
           buttonUsed:'NONE', $  ; A flag to indicate which button is used.
           drawColor:drawColor, $; The rubberband box color.
           punit:punit, $        ; Passing the position unit
           x0:-1.0, $
           x1:-1.0, $
           y0:-1.0, $
           y1:-1.0 }


endelse

  
  info.scale=scale


   ; Store the info structure.

Widget_Control, tlb, Set_UValue=info, /No_Copy

   ; Start the program going.

XManager, 'pvdiagram', tlb, /No_Block
END