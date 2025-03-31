;; up (0) - right(90) - down (180) - left (270)
;; theta in degree 
;; Fourier domain calculation with ramp filter (Zero padding)
;; dc component correction !!!

function tomography, sinogram, theta, sinc=sinc, _extra=extra
  if n_elements(display) eq 0 then display=0
  sz=size(sinogram)
  result=fltarr(sz[1], sz[1])
  ntheta=n_elements(theta)
  xc=(sz[1]-1)*0.5
  yc=(sz[1]-1)*0.5
  xp=findgen(sz[1])-xc
  yp=findgen(sz[1])-yc
  xxp=rebin(xp, sz[1], sz[1])
  yyp=rebin(transpose(yp), sz[1], sz[1])
  dist = sqrt(xxp^2.+yyp^2.)
  result[where(dist gt 0.5*sz[1])] = !values.f_nan 
      
  n_work = 2^(ceil(alog(sz[1])/alog(2.))+1)
  omega = fft_freq(findgen(n_work))
; -- sinc function as the window function 
  if keyword_set(sinc) then begin
    wl = 0
    wh = 0.5  ; Niquist freq.
    xx = !dpi*(abs(omega)-wl)/(wh-wl)
    qs = sin(xx)/xx
    qs[where(abs(omega) le wl, /null)] = 1.
    qs[where(abs(omega) gt wh, /null)] = 0.
    w = abs(omega)*qs
  endif else w = abs(omega)    
  obs = fltarr(n_work)
  for i=0, ntheta-1  do begin
    obs[0:sz[1]-1] = sinogram[*, i]
    obs_f = fft(obs, -1)
    g_tt = (float(fft(obs_f*w, 1)))[0:sz[1]-1]
    dist = xxp*cos(theta[i]*!dtor)-yyp*sin(theta[i]*!dtor)
    dum = interpol(g_tt, xp, dist, _extra=extra)
    result = result+dum
  endfor
  result1=result*!dpi/ntheta
  result1 = result1-mean(result1, /nan)$
                   +mean(total(sinogram, 1))/total(finite(result))
;  stop
  return, result1
end

path='/data/home/chokh/IDL_default/idl_lib/tomography'
cd, path
restore, 'sinogram1.sav', /verbose
sav_filename = 'tomo_result1.sav'
;drange = [0.98, 1.06]
;dif_drange = 0.2*[-1, 1]
;hist_binsize = 0.01

drange = [0, 255]
dif_drange = 20.*[-1, 1]
hist_binsize = 1

res=tomography(sinogram, theta, /spline)
sz=size(sinogram)
xc=(sz[1]-1)*0.5
yc=(sz[1]-1)*0.5
xp=findgen(sz[1])-xc
yp=findgen(sz[1])-yc
xxp=rebin(xp, sz[1], sz[1])
yyp=rebin(transpose(yp), sz[1], sz[1])
if n_elements(r_sun) then begin
  data[where(sqrt(xxp^2.+yyp^2.) lt r_sun)] = !values.f_nan
  res[where(sqrt(xxp^2.+yyp^2.) lt r_sun)] = !values.f_nan
endif

diff=res-data
real=where(finite(res))
hist=histogram(diff[real], location=x1, binsize=hist_binsize, $
               min=dif_drange[0], max=dif_drange[1])
w=window(dim=[8d2, 8d2])
p=plot(x1, hist, xthick=2, ythick=2, thick=2, /current, $
  axis_style=2, xtitle='$\Delta$ value', ytitle='# of pixel', $
  title='Difference Image Histogram', $
  font_size=15, /histogram, pos=[0.16, 0.13, 0.92, 0.88], $
  xr=dif_drange)
p.title.font_size=15
p.title.font_style=1
p1=plot([0, 0], [0, p.yr[1]], /over, '--2')
pos1=[150, 520, 300, 670]
pos_del=[0, -200, 0, -200]
pos_c=[120, 0, 30, 0]
im1=image(data, /current, title='Original', min=drange[0], max=drange[1], $
  pos=pos1, /dev)
im2=image(res, /current, title='Reconstructed', min=drange[0], max=drange[1], $
  pos=pos1+pos_del, /dev)
im3=image(diff, /current, title='Difference', min=dif_drange[0], max=dif_drange[1], $
  pos=pos1+pos_del*2., /dev, rgb_table=39)
c1=colorbar(target=im1, orient=1, pos=[1.05, 0, 1.1, 1], /relative, $
  tickv=[drange[0], mean(drange), drange[1]], textpos=1)
c2=colorbar(target=im2, orient=1, pos=[1.05, 0, 1.1, 1], /relative, $
  tickv=[drange[0], mean(drange), drange[1]], textpos=1)
c3=colorbar(target=im3, orient=1, pos=[1.05, 0, 1.1, 1], /relative, $
  tickv=[dif_drange[0], mean(dif_drange), dif_drange[1]], textpos=1, minor=3)
im4=image(sinogram, /current, title='Sinogram', $
          aspect_ratio=0, pos=[0.75, 0.4, 0.85, 0.8])
;stop
;w.save, 'tomography.pdf', page_size=8.5*[1, 1], width=8.5
save, data, res, diff, r_sun, sinogram, theta, $
      filename=sav_filename
end