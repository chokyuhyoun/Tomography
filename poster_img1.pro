path='/data/home/chokh/pinto/sdo_data'
cd, path
sz=4d2
dum=file_search('AIA2017*.fits')
read_sdo, dum[0], index
xc=index.crpix1*sz/index.naxis1
yc=index.crpix2*sz/index.naxis2
r_sun=index.r_sun*sz/index.naxis1

f=file_search('box*.sav')
restore, f[3]
aia_lct, rr, gg, bb, wave=193, /load

xp=findgen(sz)-xc
xxp=rebin(xp, sz, sz, sz)
yyp=rebin(reform(xp, 1, sz), sz, sz, sz)
zzp=rebin(reform(xp, 1, 1, sz), sz, sz, sz)
dist_2d=sqrt(xxp^2.+zzp^2.)
dist_3d=sqrt(xxp^2.+yyp^2.+zzp^2.)
box[where(dist_3d lt r_sun)]=!values.f_nan

res_path=path+'/img'
file_mkdir, res_path
cd, res_path
num=0
min=1d3
max=5d4
w1=window(dim=[8d2, 8d2], background_color='black')
t1=text(50, 100, 'Default : 2017. 8. 21. 18:35 UT', $
        /dev, color='white', font_size=20)

if 0 then begin
for i=0, 360 do begin
  num=num+1
  t3d, /reset, rotate=[0, 0, -i]
  pos=!p.t[0:2, 0:2]##[[xxp[*]], [yyp[*]], [zzp[*]]]
  pos=temporary(pos)+0.5*(sz-1)
  result=interpolate(box, pos[*, 0], pos[*, 1], pos[*, 2], $
                     missing=!values.f_nan)
  result=reform(temporary(result), sz, sz, sz)
  t3d, /reset
  result[where(dist_2d lt r_sun and yyp gt 0.)]=0
  if n_elements(im) eq 0 then begin
    im=image(total(result, 2, /nan), pos=[0.1, 0.1, 0.9, 0.9], $
             /current, rgb_table=[[rr], [gg], [bb]])
  endif else im.setdata, total(result, 2, /nan)
  im.min=min
  im.max=max
  if n_elements(t2) then t2.delete
  t2=text(50, 50, 'Angle : ('+ $
                    '0$\deg$, '+ $
                    '0$\deg$, '+ $
                    string(i mod 360, f='(i0)')+'$\deg$)', $
          /dev, color='white', font_size=20)
;  stop                  
  w1.save, string(num, f='(i04)')+'.png', resol=96
endfor
endif
num=360
cum=0
for i=0, 280 do begin
  num=num+1
  xang=70*sin(i*2.*!dpi/(70*4))
  t3d, /reset, rotate=[xang, 0, 0]
  pos=!p.t[0:2, 0:2]##[[xxp[*]], [yyp[*]], [zzp[*]]]
  pos=temporary(pos)+0.5*(sz-1)
  result=interpolate(box, pos[*, 0], pos[*, 1], pos[*, 2], $
    missing=!values.f_nan)
  result=reform(temporary(result), sz, sz, sz)
  t3d, /reset
  result[where(dist_2d lt r_sun and yyp gt 0.)]=0
  
  if n_elements(im) eq 0 then begin
    im=image(total(result, 2, /nan), pos=[0.1, 0.1, 0.9, 0.9], $
             /current, rgb_table=[[rr], [gg], [bb]])
  endif else im.setdata, total(result, 2, /nan)
  im.min=min
  im.max=max
  if n_elements(t2) then t2.delete
  t2=text(50, 50, 'Angle : ('+$
                    string(xang, f='(i0)')+'$\deg$, '+ $
                    '0$\deg$, '+ $
                    '0$\deg$)', $
    /dev, color='white', font_size=20)
  w1.save, string(num, f='(i04)')+'.png', resol=96    
endfor

end