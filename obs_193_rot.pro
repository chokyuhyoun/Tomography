cd, '/data/home/chokh/pinto/sdo_data/img5'
restore, '/data/home/chokh/pinto/sdo_data/sdo_sinogram_0193.sav'

w1 = window(dim=[8d2, 8d2], background_color='black')
aia_lct, rr, gg, bb, wave=193
im1 = image_kh(sinogram[*, 0, *], /current, min=0, max=500, $
               pos=[0.05, 0.05, 0.95, 0.95], rgb_table=[[rr], [gg], [bb]])
time = t_obs+julday(8, 21, 2017, 0, 0, 0)
time_str = string(time, f='(c(CYI, CMOI02, CDI02, "   ", CHI02, ":", CMI02, ":", CSI02))')
t1 = text(0.05, 0.05, time_str[0], color='white', $
          font_name='malgun gothic', font_size=25)
stop
for i=0, n_elements(time)-1 do begin
  t1.delete
  setdata, im1, reform(sinogram[*, i, *])
  t1 = text(0.05, 0.05, time_str[i], color='white', $
            font_name='malgun gothic', font_size=25)
  im1.save, string(i, f='(i03)')+'.png', resol=96          
endfor

end