FUNCTION GETSERIES, first, last, step
  total=abs(last-first)/step+1  
  
  if first ge last then begin
    return,(findgen(total))*(-step)+first
  endif else begin
    return,(findgen(total))*step+first
  
  endelse
  
  
END
