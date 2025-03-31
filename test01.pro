cd, '/data/home/chokh/tomography-1902'
restore, 'dem_result.sav' ;; result
restore, 'tdem_result.sav' ;; emiss, tot_e, avg_te

sz = 128
xp = (findgen(sz)-0.5*(sz-1))/r_sun
tbin_l = 10.^(lgtmin+findgen(nlgt)*dlgt)
tbin_h = 10.^(lgtmin+(findgen(nlgt)+1.)*dlgt)
tbin = 0.5*(tbin_l+tbin_h)
tbin2 = rebin(reform(tbin, 1, 1, nlgt), sz, sz, nlgt) 
tbin3 = rebin(reform(tbin, 1, 1, 1, nlgt), sz, sz, sz, nlgt)


w01 = window(dim=[8d2, 8d2])
ct01 = colortable(74, /rev)
trange1 = [1d6, 2d6]
trange2 = [1d6, 2d6]
xr = 0.7*[-1, 1]
yr = [-1.25, 0.0]
cbpos = [0, -0.05, 1, 0]

;; 1. Observation -> DEM
obs_dem = reform(result[*, *, 820, *])
obs_em_t = total(obs_dem*tbin2, 3, /nan)/total(obs_dem, 3, /nan)
im01 = image_kh(obs_em_t, xp, xp, /current, layout=[2, 2, 1], $
                rgb_table=ct01, min=trange1[0], max=trange1[1], $
                xr=xr, yr=yr, xtickformat='(a1)', $
                title='EM weighted $T_e$ from obs.')
p01 = plot(cos(findgen(360)*!dtor), sin(findgen(360)*!dtor), over=im01, /data)
cb01 = colorbar(target=im01, /relative, pos=cbpos, $
                title='$T_e$', orient=0, textpos=0, $
                font_style=1, font_size=11, font_name='malgun gothic', $
                /border, ticklen=0.5, minor=1, thick=1.5, $
                tickinterval=5d5)
 
;; 2. reconstruction -> ne^2 weighted
recon_dem = sol_los_int(emiss, r_sun=r_sun)
te1 = total(recon_dem*tbin2, 3, /nan)/total(recon_dem, 3, /nan)
im02 = image_kh(te1, xp, xp, /current, layout=[2, 2, 2], $
                rgb_table=ct01, min=trange1[0], max=trange1[1], $
                xr=xr, yr=yr, xtickformat='(a1)', $
                title='EM weighted $T_e$ from Tomography')
p02 = plot(cos(findgen(360)*!dtor), sin(findgen(360)*!dtor), over=im02, /data)
cb02 = colorbar(target=im02, /relative, pos=cbpos, $
                title='$T_e$', orient=0, textpos=0, $
                font_style=1, font_size=11, font_name='malgun gothic', $
                /border, ticklen=0.5, minor=1, thick=1.5, tickinterval=5d5)

 
;; 3. reconstruction -> ne^2 > 0 weighted
te_dist = total(emiss*tbin3, 4)/total(emiss, 4)
te2 = sol_los_int(te_dist, r_sun=r_sun)
im03 = image_kh(te2, xp, xp, /current, layout=[2, 2, 3], $
                rgb_table=ct01, min=trange2[0], max=trange2[1], $
                xr=xr, yr=yr, xtickformat='(a1)', $
                title='EM>0 weighted $T_e$ from Tomography')
p03 = plot(cos(findgen(360)*!dtor), sin(findgen(360)*!dtor), over=im03, /data)
cb03 = colorbar(target=im03, /relative, pos=cbpos, $
                title='$T_e$', orient=0, textpos=0, $
                font_style=1, font_size=11, font_name='malgun gothic', $
                /border, ticklen=0.5, minor=1, thick=1.5, tickinterval=5d5)
stop
;; 4. reconstruction -> ne weighted
ne1 = sqrt(total(emiss, 4, /nan))
recon_dem3 = sol_los_int(ne1, r_sun=r_sun)
te3 = total(recon_dem3*tbin2, 3, /nan)/total(recon_dem3, 3, /nan)
im04 = image_kh(te3, xp, xp, /current, layout=[2, 2, 4], $
                rgb_table=ct01, min=trange2[0], max=trange2[1], $
                xr=xr, yr=yr, xtickformat='(a1)', $
                title='$n_e$ weighted $T_e$ from Tomography')
p04 = plot(cos(findgen(360)*!dtor), sin(findgen(360)*!dtor), over=im04, /data)
cb04 = colorbar(target=im04, /relative, pos=cbpos, $
                title='$T_e$', orient=0, textpos=0, $
                font_style=1, font_size=11, font_name='malgun gothic', $
                /border, ticklen=0.5, minor=1, thick=1.5, tickinterval=5d5)


end