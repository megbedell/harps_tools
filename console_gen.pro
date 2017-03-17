; Generate text file of RVs for input to systemic console

PRO console_gen, file_in, file_out


restore, file_in


ref_date = floor(date[0])


;date = date-ref_date

;date = date - 2400000.5d0  ; MJD


rv = rv - avgsigclip(rv, 3, /med, /iter)
rv = rv * 1e3
sig = sig * 1e3

bad_sig = where(sig le 0.1d0)
if (bad_sig[0] ne -1) then sig[bad_sig] = 1.0d0

;good = where(sig le 2.5)
;date=date[good]
;rv=rv[good]
;sig=sig[good]

short = where(exp le 600.0)
if (short[0] ne -1) then sig[short] = sqrt((sig[short])^2 + 1.1d0^2) ; inflate error if exposure <= 10 min

bin, date, rv, 20./60./24., date2, rv2, sig2, error=sig  ; bin by 20 minutes
;bin, date, rv, 2.5/24., date2, rv2, sig2, error=sig  ; bin by 2.5 hours
date = date2
rv = rv2
sig = sig2

;sig[where(sig lt 1.0)] = 1.0


n = n_elements(date)

openw, u, file_out, /get_lun
for i=0,n-1 do printf,u,date[i],rv[i],sig[i],format='(f13.4,8x,f12.3,8x,f8.3)'
free_lun,u

END
