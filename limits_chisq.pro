PRO limits_chisq, e_val,infile,outfile
; e_val: eccentricity
; infile: IDL savefile contining dates, resids, sig arrays
; outfile: were the results should be saved

;outfile = 'test.dat'

;Defaults & constants
;infile = 'resids.dat'                    ;file with dates and residuals 
m_star = 1.0                           ;stellar mass in M_sun
G = 6.674d-11      ;whatever G is in SI
M_sun = 1.989d30   ;kg
M_Jup = 1.898d27   ;kg
AU = 1.496d11   ;meters



;n_a = 350                              ;number of semi-major axis steps
;semi_vals = 10.^(dindgen(n_a)/100)/100. + 0.0092200144d0    ;values in AU
;n_a = 4500
;semi_vals = 10.^(dindgen(n_a)/1000)/100. + 0.00001d0    ;values in AU

n_p = 701
p_vals = 10.^(dindgen(n_p)/100)/100. + 0.1d0 ;period values in days
n_iter = 100                           ;number of random orbits to analyze
;k_vals = dindgen(400)/10. + 0.1         ;velocity semi-amplitude values to step through
k_vals = dindgen(20000)/4. + 2.0


;Restore residuals file.
restore, infile

;Generate false residuals:
;dates = interpol(dates,100)
;resids_gen = interpol(resids,100)
;resids = randomize(resids_gen)
;sig_gen = interpol(sig,100)
;sig = randomize(sig_gen)

epochs = n_elements(dates)

k_limits = fltarr(n_p)
detect = intarr(n_iter)

fit = poly_fit(dates,resids,0,measure_errors=sig,chisq=chisq,/double,status=status)
errors = sig * sqrt(chisq/(epochs-1))

print,"Reduced chisquared of flat line = ",chisq/(epochs-1)



for i=0,n_p-1 do begin
 limit = 0
 j = -1
 while not(limit) do begin
  j = j + 1
  k = -1
  detect[*] = 1
  while (k lt n_iter-1) and (n_elements(where(detect lt 1)) le 0.01*n_iter) do begin
   k = k + 1
   Tp = min(dates) + (randomu(seed) mod 1) * max(p_vals)
   w = (randomu(seed) mod 1) * 360
   rand_ind = fix(n_elements(dates)*randomu(seed,n_elements(dates)))   ;bootstrapping of residuals
   data = calc_rv(Tp, p_vals[i], k_vals[j],e_val,w,dates) + resids[rand_ind]
   fit = poly_fit(dates,data,0,measure_errors=errors[rand_ind],chisq=chisq,/double,status=status)
   if (status ne 0) then print,"Fit has gone wrong at p-value",p_vals[i]
   if (mpchitest(chisq,epochs-1,/sigma) lt 3) then detect[k] = 0 else detect[k] = 1
  endwhile
  a = where(detect eq 1)
  if (n_elements(a) ge n_iter*0.99) then limit = 1  
  if (j eq n_elements(k_vals)-1) then begin
   print,"Reached upper limit of k at p-value",p_vals[i]
   limit = 1
  endif
 endwhile    
 k_limits[i] = k_vals[j]
 print,"limit placed for period value",p_vals[i]
endfor

semi_vals = fltarr(n_p)

for i=0,n_p-1 do semi_vals[i] = calc_a(p_vals[i],k_limits[i],e_val,m_star)

m_limits = k_limits * sqrt((1-e_val^2)/G) * sqrt(semi_vals * AU * m_star * M_sun) / M_Jup  ;mass limits in M_Jup

save, p_vals, m_limits, k_limits, e_val, semi_vals, filename=outfile
cgplot,p_vals,m_limits,xtitle="Period (days)",ytitle="Mass (M_Jup)",/xlog,/ylog


END
