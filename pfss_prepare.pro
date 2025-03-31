;; PFSS
cd, '/data/home/chokh/tomography-1902'

pfss_sav=file_search('pfss.sav', count=n)
@pfss_data_block
if n eq 0 then begin
  pfss_restore, pfss_time2file('2019-02-15_00:00:00', /ssw_cat, /url)
  save, br, bth, bph, nr, nlat, nlon, rix, $
    lat, lon, theta, phi, phiat, phibt, l0, b0, now, $
    filename='pfss.sav'
endif else restore, 'pfss.sav'

l_phi = [0 : 360 : 30]
l_theta = [30 : 150 : 30]
n_phi = n_elements(l_phi)
n_theta = n_elements(l_theta)
n_lines = n_phi*n_theta
str = replicate(1.01, n_lines)
stph = rebin(l_phi, n_phi, n_theta)*!dtor
stth = rebin(transpose(l_theta), n_phi, n_theta)*!dtor

pfss_trace_field, kind
pfss_to_spherical, pfss_data

save, ptr, ptth, ptph, kind, pfss_data, filename='pfss_result.sav'
end