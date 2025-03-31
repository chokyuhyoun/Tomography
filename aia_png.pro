pro aia_png, file, _extra = extra
  window, /free, xs=8d2, ys=8d2, _extra=extra
  name = change_extension(file)
  read_sdo, file[0], index
  for i=0, n_elements(file)-1 do begin
    read_sdo, file[i], index, data
    dum1 = congrid(data, 8d2, 8d2)
    tvscl, aia_intscale(dum1, exptime=index.exptime, wavelnth=index.wavelnth)
    write_png, name[i], tvrd() 
;    stop
  endfor
end