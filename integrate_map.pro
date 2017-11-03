PRO INTEGRATE_MAP,fits=fits,vel_range=vel_range,ex_range=ex_range,outfits=outfits,eta=eta,cont=cont
    
  file=fits
  
  im=readfits(file,hd)
  
  ; Making sure the unit is T
  bunit=sxpar(hd,'BUNIT')
  if strmatch(strtrim(bunit),'JY*',/fold_case) eq 1 then begin
    print,'Unit is JY/BEAM. Convert this unit to T before calculating aboudance.'
  endif
  
  xpix=sxpar(hd,'NAXIS1')
  ypix=sxpar(hd,'NAXIS2')
  
  
  del_v=sxpar(hd,'CDELT3')/1000.0
  ref_v=sxpar(hd,'CRVAL3')/1000.0
  ref_c=sxpar(hd,'CRPIX3')
  
  v0=(ref_v-ref_c*del_v)
  
  channel=sxpar(hd,'NAXIS3')
  
  array_vel=findgen(channel)*del_v+(ref_v-ref_c*del_v)
  ; Convert velocity to km s-1
  ch_v=findgen(channel)*del_v+(ref_v-ref_c*del_v)
  
  
  if keyword_set(vel_range) then begin
    vel_range=vel_range
  endif else begin
    vel_range=[max(ch_v),min(ch_v)]
  endelse
  ; Convert from velocity range to channel range
  chl_range=round((vel_range-v0)/del_v)-1
  
  ; If the exclusion velocity range is set, calculate the correspoding channel ranges
  if keyword_set(ex_range) then chex_range=round((ex_range-v0)/del_v)-1
  
  iim=fltarr(xpix,ypix)
  
  ; If the velocity range is inverse sequence, the step should be -1!
  if chl_range[1] ge chl_range[0] then step=1 else step=-1
  
  ; Integrate the velocity range
  for i=chl_range[0],chl_range[1],step do begin
    ;if keyword_set(cont) then begin
    iim[*,*]=iim[*,*]+im[*,*,i]
    ;endif else begin
    ;    iim[*,*]=iim[*,*]+im[*,*,i]*abs(del_v)
    ;endelse
  endfor
  ; Integrate the exclusion velocity range
  if keyword_set(ex_range) then begin
    exiim=fltarr(xpix,ypix)
    if chex_range[1] ge chex_range[0] then step=1 else step=-1
    for i=chex_range[0],chex_range[1],step do begin
      exiim[*,*]=exiim[*,*]+im[*,*,i]
    endfor
    iim=iim-exiim
  endif
  
  ; Specify this keyword to normalize the integrated flux to per-
  if keyword_set(cont) then begin
    ;  cont=fltarr(xpix,ypix)
    ; for i=0,channel-1 do begin
    ;   cont[*,*]=cont[*,*]+im[*,*,i]
    ; endfor
    ;  contim=cont-iim
    nchan=(abs(chl_range[0]-chl_range[1])+1)-(abs(chex_range[0]-chex_range[1])+1)
    iim=iim/nchan
    
  endif else begin
    iim=iim*abs(del_v)
  endelse
  
  if keyword_set(eta) then begin
    print,'correct beam efficiency.'
    iim=iim/eta
  endif
  
  writefits,outfits,iim,hd
  
  
END


PRO GETVAL,header,keys,values

  keyword = strtrim(strmid(header, 0, 8), 2)
  space = strpos( header, ' ', 10) + 1
  slash = strpos( header, '/', 10)  > space
  
  N = N_elements(header)
  len = (slash -10) > 20
  len = reform(len,1,N)
  lvalue = strtrim(strmid(header, 10, len),2)
  
  nkey = n_elements(keys)
  values = dblarr(nkey)
  for i=0,nkey-1 do values[i] = lvalue[where(keyword eq keys[i])]
  if nkey eq 1 then values=values[0]
  
END


FUNCTION GETSERIES, first, last, step
  total=abs(last-first)/step+1
  
  if first ge last then begin
    return,(findgen(total))*(-step)+first
  endif else begin
    return,(findgen(total))*step+first
    
  endelse
  
  
END

