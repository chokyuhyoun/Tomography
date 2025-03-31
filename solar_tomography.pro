;; theta in degree (input : -180(up, past) -90 (right) 0 (down = obs.) 90 (left, future)
;; Fourier domain calculation with ramp filter (Zero padding)
;; in the case of the screening by the solar disk
;; dc component correction !!!
;; sinogram + opposite part sinogram method - > DC component
;; results are almost same with traditional tomography
;; Update : 2019. 4. 8


function solar_tomography, sinogram1, theta1, r_sun=r_sun, $
                            lat=lat, sinc=sinc, _extra=extra
  use = where(theta1 ge -180. and theta1 lt 180.)
  theta = theta1[use]
  sinogram = sinogram1[*, use]
   
;  if n_elements(lat) eq 0 then lat = 0.
;  theta = theta1+180d0    ; new angle. 180 = observer
;  theta = (theta+360d0) mod 360d0
;  if min(theta) lt 0. or max(theta) gt 360. then begin
;    dtheta = (theta-shift(theta, 1))[1:*]
;    dum1 = where(dtheta lt -350.)
;    sinogram = sinogram1[*, dum1[0]+1:dum1[1]]
;    theta = theta[dum1[0]+1:dum1[1]]  
;  endif else sinogram = sinogram1
  
  sz=size(sinogram)
  result=fltarr(sz[1], sz[1])
  ntheta=n_elements(theta)
  xc=(sz[1]-1)*0.5
  yc=(sz[1]-1)*0.5
  xp=findgen(sz[1])-xc
  yp=findgen(sz[1])-yc
  xxp=rebin(xp, sz[1], sz[1])
  yyp=rebin(transpose(yp), sz[1], sz[1])
  dist1 = sqrt(xxp^2.+yyp^2.)
  r_eff = r_sun*cos(lat)
  result[where(dist1 gt 0.5*sz[1])] = !values.f_nan
  result[where(dist1 le r_eff)] = !values.f_nan

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

;  obs = fltarr(n_work)
;  im1 = image_kh(result)    ;-------
  
  if r_eff gt 0 then begin
    op_theta = theta+180.
    op_theta[where(op_theta gt 180.)] = $
      temporary(op_theta[where(op_theta gt 180.)])-360.
    op_ind = interpol(findgen(sz[2]), theta, op_theta)
    op_sinogram = interpolate(sinogram, findgen(sz[1]), op_ind, /grid)
    xxp1 = rebin(xp, sz[1], sz[2])
    op_sinogram[where(abs(xxp1) gt r_eff)] = 0.
    dc = sinogram+reverse(op_sinogram, 1)  
  endif else dc = sinogram
  real = where(theta ge -90. and theta lt 90.)
  r_theta = -theta[real]
  real_dc = dc[*, real]                     
  
;  oppo_ang = (theta+180.) mod 360.
;  oppo_ind = arr_eq(theta, oppo_ang)
;  oppo = sinogram[*, oppo_ind]
;  xxp1 = rebin(xp, sz[1], ntheta)
;  oppo[where(abs(xxp1) gt r_eff)] = 0.
;  dc = sinogram+reverse(oppo, 1)
;  real_part = where(theta ge 90. and theta lt 270.) ;; 180 = observer direction
;  r_theta = -theta[real_part]-180.  ;return to the obs. angle (0 = observer)
;  real_dc = dc[*, real_part]
;  stop
  for i=0, n_elements(real)-1  do begin
    obs = fltarr(n_work)
    obs[0:sz[1]-1] = real_dc[*, i]    
    obs_f = fft(obs, -1)
    g_tt1 = (float(fft(obs_f*w, 1)))[0:sz[1]-1]
    dist = xxp*cos(r_theta[i]*!dtor)+yyp*sin(r_theta[i]*!dtor)
    dum = interpol(g_tt1, xp, dist, _extra=extra)
    dum1 = reform(dum, sz[1], sz[1])
    result = result+dum1
;    im1.setdata, result ;------------
;    stop
  endfor
;  stop
  result1 = result*!dpi/n_elements(real)
  result1 = result1-mean(result1, /nan)$
                   +mean(total(real_dc, 1))/total(finite(result))
  return, result1
end

path='/data/home/chokh/IDL_default/idl_lib/tomography'
cd, path
restore, 'sinogram1_fanbeam.sav', /verbose
sav_filename = 'solar_tomo_fanbeam_result.sav'
;drange = [0.98, 1.06]
;dif_drange = 0.2*[-1, 1]
;hist_binsize = 0.01

drange = [0, 255]
dif_drange = 20.*[-1, 1]
hist_binsize = 1

;stop
res=solar_tomography(sinogram, theta, r_sun=r_sun, /spline)
sz=size(sinogram)
xc=(sz[1]-1)*0.5
yc=(sz[1]-1)*0.5
xp=findgen(sz[1])-xc
yp=findgen(sz[1])-yc
xxp=rebin(xp, sz[1], sz[1])
yyp=rebin(transpose(yp), sz[1], sz[1])
data[where(sqrt(xxp^2.+yyp^2.) lt r_sun)] = 255.

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
;w.save, 'solar_tomography.pdf', page_size=8.5*[1, 1], width=8.5
;w.save, 'solar_tomography.png', resol=200
save, data, res, diff, r_sun, sinogram, theta, $
  filename=sav_filename
end
