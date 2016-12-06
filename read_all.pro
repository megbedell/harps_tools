PRO read_all
;outputs the RV and all activity-based information for all stars at each observational epoch


data_dir = '/Users/mbedell/Documents/Research/HARPSTwins/Data/Reduced/'
out_dir = '/Users/mbedell/Documents/Research/HARPSTwins/Results/'

;filenames = [file_search(data_dir+'*/HARPS*ccf_G2_A.fits'), file_search(data_dir+'18Sco/*/HARPS*ccf_G2_A.fits'), file_search(data_dir+'archive/*/HARPS*ccf_G2_A.fits')]
filenames = [file_search(data_dir+'*/HARPS*ccf_G2_A.fits'), file_search(data_dir+'archive/*/HARPS*ccf_G2_A.fits')]


;bad_stars = ['HIP65708','HIP83276','HIP19911','HIP6407','HIP18844','HIP67620','HIP73241','HIP103983']

;read in some external info about the sample:
restore,'/Users/mbedell/Documents/Research/HARPSTwins/Data/save.dat'
HIPname = strcompress(string(long64(c1)),/remove_all)
HDname = strcompress(string(long64(c2)),/remove_all)
bvs = c9   ; needed for log(R'HK)
Teffs = c3   ; needed for jitter prediction


;retrieve RVs for all data, sorted by date:

n_files = n_elements(filenames)
names = strarr(n_files)
t_rv = dblarr(n_files)
t_rv_ccfiter = dblarr(n_files)
t_rv_ccf = dblarr(n_files)
t_sig = dblarr(n_files)
t_date = dblarr(n_files)
t_am = dblarr(n_files)
t_snr = dblarr(n_files)
t_drift = dblarr(n_files)
t_exp = dblarr(n_files)
t_ra = dblarr(n_files)
t_dec = dblarr(n_files)
t_berv = dblarr(n_files)
for i=0,n_files-1 do begin
  header = headfits(filenames[i])
  names[i] = strcompress(sxpar(header,'OBJECT'),/remove_all)
  if (names[i] eq 'HIP-31592') then continue
  line = where(strmid(header,0,24) eq 'HIERARCH ESO DRS CCF RVC')
  t_rv[i] =double(strmid(header[line],strpos(header[line],'=')+2,(strpos(header[line],'/')-strpos(header[line],'=')-2)))
  line = where(strmid(header,0,29) eq 'HIERARCH ESO DRS DRIFT SPE RV')
  t_drift[i] = double(strmid(header[line],strpos(header[line],'=')+2,(strpos(header[line],'/')-strpos(header[line],'=')-2)))
  line = where(strmid(header,0,26) eq 'HIERARCH ESO DRS CCF NOISE')
  t_sig[i] = double(strmid(header[line],strpos(header[line],'=')+2,(strpos(header[line],'/')-strpos(header[line],'=')-2)))
  line = where(strmid(header,0,20) eq 'HIERARCH ESO DRS BJD')
  t_date[i] = double(strmid(header[line],strpos(header[line],'=')+2,(strpos(header[line],'/')-strpos(header[line],'=')-2))) 
  ;t_date[i] = double(strcompress(sxpar(header,'MJD-OBS'),/remove_all)) 
  line = where(strmid(header,0,27) eq 'HIERARCH ESO TEL AIRM START')
  t_am[i] = double(strmid(header[line],strpos(header[line],'=')+2,(strpos(header[line],'/')-strpos(header[line],'=')-2)))
  line = where(strmid(header,0,27) eq 'HIERARCH ESO DRS SPE EXT SN')  ; SNR for each order listed separately
  sns = []
  for j=0,n_elements(line)-1 do sns = [sns, double(strmid(header[line[j]],strpos(header[line[j]],'=')+2,(strpos(header[line[j]],'/')-strpos(header[line[j]],'=')-2)))]
  t_snr[i] = median(sns)  ; save median SNR

  t_exp[i] = double(strcompress(sxpar(header,'EXPTIME'),/remove_all))
  t_ra[i] = double(strcompress(sxpar(header,'RA'),/remove_all))
  t_dec[i] = double(strcompress(sxpar(header,'DEC'),/remove_all))
  

  line = (where(strmid(header,0,21) eq 'HIERARCH ESO DRS BERV'))[0]
  t_berv[i] = double(strmid(header[line],strpos(header[line],'=')+2,(strpos(header[line],'/')-strpos(header[line],'=')-2)))  

  
endfor

feb2012 = where(t_date gt 2455983.5 and t_date lt 2455987.5)
if (feb2012[0] ne -1) then t_rv[feb2012] = t_rv[feb2012] - t_drift[feb2012]/1.0d3
apr2016 = where(t_date gt 2457505.4 and t_date lt 2457506.4)
if (apr2016[0] ne -1) then t_rv[apr2016] = t_rv[apr2016] - t_drift[apr2016]/1.0d3

jul2015 = where(t_date gt 2457218.5)   ; select all measurements post-fiber upgrade
if (jul2015[0] ne -1) then t_rv[jul2015] = t_rv[jul2015] - 15.4/1.0d3
if (jul2015[0] ne -1) then t_sig[jul2015] = sqrt((t_sig[jul2015])^2 + (0.4/1.0d3)^2)


index = sort(t_date)
t_date = t_date[index]  ; sort data by date
t_rv = t_rv[index]
;t_rv_ccfiter = t_rv_ccfiter[index]
;t_rv_ccf = t_rv_ccf[index]
t_drift = t_drift[index]
t_sig = t_sig[index]
names = names[index]
filenames = filenames[index]
t_am = t_am[index]
t_snr = t_snr[index]
t_exp = t_exp[index]
t_ra = t_ra[index]
t_dec = t_dec[index]
t_berv = t_berv[index]

;rename HD & other aliases from archive files:

remchar, names, '-'

if ((where(names eq '18_Sco'))[0] ne -1) then names[where(names eq '18_Sco')] = 'HIP79672'
if ((where(names eq '18Sco'))[0] ne -1) then names[where(names eq '18Sco')] = 'HIP79672'
if ((where(names eq 'HD59711A'))[0] ne -1) then names[where(names eq 'HD59711A')] = 'HIP36512'
if ((where(names eq 'HDTE221287'))[0] ne -1) then names[where(names eq 'HDTE221287')] = 'HD221287'
if ((where(names eq 'HD157347_std'))[0] ne -1) then names[where(names eq 'HD157347_std')] = 'HIP85042'
if ((where(names eq 'HIP7505'))[0] ne -1) then names[where(names eq 'HIP7505')] = 'HIP7585'


for i=0,n_elements(HDname)-1 do begin
  if ((where(names eq 'HD'+HDname[i]))[0] ne -1) then names[where(names eq 'HD'+HDname[i])] = 'HIP'+HIPname[i]
endfor

save, names, t_date, t_rv, t_sig, t_am, t_snr, t_drift, t_exp, t_ra, t_dec, t_berv, filenames, filename=out_dir+'all_result.dat'

;loop through all stars and derive activity parameters for each:

index1 = sort(names)
index2 = uniq(names[index1])

n_unique = n_elements(index2)


for i=0, n_unique-1 do begin
 star_name = names[index1[index2[i]]]

 ;check = where(star_name eq bad_stars, count)
 ;if (count gt 0) then continue

 index3 = where(names eq star_name)
 files = filenames[index3]
 date = t_date[index3]
 rv = t_rv[index3]
 sig = t_sig[index3]
 airm = t_am[index3]
 snr = t_snr[index3]
 exp = t_exp[index3]
 ra = t_ra[index3]
 dec = t_dec[index3]
 berv = t_berv[index3]
 bv = bvs[where('HIP'+HIPname eq star_name)]
 Teff = Teffs[where('HIP'+HIPname eq star_name)]

 drift=t_drift[index3]

 good = indgen(n_elements(date))
 if (star_name eq 'HIP10303') then good = where(date lt 2456164.0 or date gt 2456166.0)  ; eliminate data from wrong star
 ;if (star_name eq 'HIP11915') then good = where(date lt 2456708.0 or date gt 2456709.0)  ; eliminate un-drift-corrected data
 if (star_name eq 'HIP29432') then good = where(date gt 2455700.0)  ; eliminate bad archive data
 if (star_name eq 'HIP64713') then good = where(date lt 2456299.0 or date gt 2456300.0)  ; eliminate bad data
 if (star_name eq 'HIP54102') then good = where(date lt 2455983.6 or date gt 2455983.8)
 if (star_name eq 'HIP64673') then good = where(date gt 2455985.0)
 if (star_name eq 'HIP74432') then good = where(date gt 2456048.0)
 if (star_name eq 'HIP10175') then good = where(date lt 2457026.0 or date gt 2457028.0)  ; eliminate wrong star
 date = date[good]
 rv = rv[good]
 sig = sig[good]
 airm = airm[good]
 snr = snr[good]
 exp = exp[good]
 ra = ra[good]
 dec = dec[good]
 
 berv = berv[good]

 drift=drift[good]

 n_pts = n_elements(good)
 logRhk = fltarr(n_pts)

 Shk = fltarr(n_pts)
 Shk_err = fltarr(n_pts)
 jitter = fltarr(n_pts)  ; predicted jitter for activity level
 bis = fltarr(n_pts)  ; bisector span, km/s
 fwhm = fltarr(n_pts) ; CCF FWHM, km/s
 for j=0,n_pts-1 do begin
  spec_file = strmid(filenames[index3[j]],0,strpos(filenames[index3[j]],'ccf'))+'s1d_A.fits'
  bis_file = strmid(filenames[index3[j]],0,strpos(filenames[index3[j]],'ccf'))+'bis_G2_A.fits'
  Rprime_HK = calc_rhk(spec_file, rv[j], bv, Teff)
  tmp = calc_rhk(spec_file, rv[j], bv, Teff, /shk)
  shk[j] = tmp[0]
  shk_err[j] = tmp[1]
  logRhk[j] = alog10(Rprime_HK)
  jitter[j] = calc_rhk(spec_file, rv[j], bv, Teff, /jitter)
  header = headfits(bis_file)
  line = where(strmid(header,0,25) eq 'HIERARCH ESO DRS BIS SPAN')
  bis[j] = double(strmid(header[line],strpos(header[line],'=')+2,(strpos(header[line],'/')-strpos(header[line],'=')-2)))
  line = where(strmid(header,0,25) eq 'HIERARCH ESO DRS CCF FWHM')
  fwhm[j] = double(strmid(header[line],strpos(header[line],'=')+2,(strpos(header[line],'/')-strpos(header[line],'=')-2)))
 endfor

  jul2015 = where(date gt 2457218.5)   ; select all measurements post-fiber upgrade
  if (jul2015[0] ne -1) then fwhm[jul2015] = fwhm[jul2015] - 15.7/1.0d3
  if (jul2015[0] ne -1) then bis[jul2015] = bis[jul2015] - 13.7/1.0d3
  

  dat_file = out_dir+star_name+'_result.dat'
  vels_file = out_dir+star_name+'_result.vels'
 
 
 save,files,date,rv,sig,airm,berv,snr,exp,ra,dec,Shk,Shk_err,logRhk,jitter,bis,fwhm,filename=dat_file

 console_gen,dat_file,vels_file
endfor

END


