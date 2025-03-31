function sol_los_int, cube, r_sun=r_sun
  if n_elements(r_sun) eq 0 then r_sun = 50.762134
  cube1 = cube
  sz = (size(cube))[1]
  xp = findgen(sz)-0.5*(sz-1)
  xxp = rebin(xp, sz, sz, sz)
  yyp = rebin(reform(xp, 1, sz), sz, sz, sz)
  zzp = rebin(reform(xp, 1, 1, sz), sz, sz, sz)
  inside = where((xxp^2.+zzp^2. lt r_sun^2. and yyp gt 0) or $
                 (xxp^2.+yyp^2.+zzp^2. lt r_sun^2.))
  if (size(cube))[0] eq 3 then cube1[inside] = 0. else begin
    for i=0, (size(cube))[4]-1 do begin
      dum = reform(cube1[*, *, *, i])
      dum[inside] = 0.
      cube1[*, *, *, i] = dum 
    endfor
  endelse
  
;  cube1[where(xxp^2.+zzp^2. lt r_sun^2. and yyp gt 0)] = 0.
;  cube1[where(xxp^2.+yyp^2.+zzp^2. lt r_sun^2.)] = 0.
;  stop
  return, total(cube1, 2, /nan)
end