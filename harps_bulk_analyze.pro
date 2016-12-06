PRO harps_bulk_analyze

dir = '/Users/mbedell/Documents/Research/HARPSTwins/Results/'
savefiles = file_search(dir+'HIP*_result.dat')

openw, u, 'harps_bulk_log.txt', /get_lun

printf,u,'HIP','RV RMS (m/s)','# RVs','Time Baseline (yr)','Age (Gyr)','log(R_HK)','SHK corr','BIS corr','FWHM corr','Linear trend','Periodogram peaks (d)',format='(a,15x,a,10x,a,5x,a,5x,a,5x,a,5x,a,7x,a,7x,a,5x,a,5x,a)'

for i=0,n_elements(savefiles)-1 do begin
	starname = strmid(savefiles[i],56,strpos(savefiles[i],'_')-56)
	starname = long(starname)
	
	
	if (starname eq 31592) then continue
	
	restore,savefiles[i]
	
	;if (avg(logrhk) gt -4.85) then continue
	
	rv = rv - avgsigclip(rv, 3, /med, /iter)
	rv = rv * 1d3
	sig = sig * 1d3
	bis = bis * 1d3
	
	short = where(exp le 600.0)
	sig[short] = sqrt((sig[short])^2 + 1.1d0^2) ; inflate error if exposure <= 10 min
	
	; if sigshk contains one or more NaN then don't use it:
	if (total(finite(shk_err)) lt n_elements(shk_err)) then shk_err = dblarr(n_elements(shk_err)) + 0.2 
	

	; bin everything by 20 minutes
	bin, date, rv, 20./60./24., date2, rv2, sig2, error=sig  ; bin by 20 minutes
	bin, date, shk, 20./60./24., date2, shk2, sigshk2, error=shk_err  ; bin by 20 minutes
	bin, date, bis, 20./60./24., date2, bis2, sig3, error=sig  ; bin by 20 minutes
	bin, date, fwhm, 20./60./24., date2, fwhm2, sig4, error=sig  ; bin by 20 minutes
	date = date2
	rv = rv2
	sig = sig2
	shk = shk2
	shk_err = sigshk2
	bis = bis2
	fwhm = fwhm2

	sig[where(sig lt 1.0)] = 1.0
	
	; if there's a significant (non-activity) linear trend then subtract it:
	
	linflag=0
	
	if (starname eq 14501 OR starname eq 14614 OR starname eq 54582 OR starname eq 64150 OR starname eq 72043 OR starname eq 73241 OR starname eq 81746 OR starname eq 87769 OR starname eq 18844 OR starname eq 19911 OR starname eq 103983 OR starname eq 6407 OR starname eq 67620 OR starname eq 83276) then linflag = 1
	
	if (linflag eq 1) then begin
		fit = linfit(date,rv,measure_errors=sig,sigma=sigslope)
		rv = rv - (fit[0]+fit[1]*date)
	endif
	
	
	shkflag = 0
	if ((pearson_func(rv,shk))[1] lt 0.05) then begin
		fit = linfit(shk,rv,measure_errors=sig,sigma=sigslope)
		rv = rv - (fit[0]+fit[1]*shk)
		shkflag = 1
	endif
	bisflag = 0
	if ((pearson_func(rv,bis))[1] lt 0.05) then begin
		fit = linfit(bis,rv,measure_errors=sig,sigma=sigslope)
		rv = rv - (fit[0]+fit[1]*bis)
		bisflag = 1
	endif
	fwhmflag = 0
	if ((pearson_func(rv,fwhm))[1] lt 0.05) then begin
		fit = linfit(fwhm,rv,measure_errors=sig,sigma=sigslope)
		rv = rv - (fit[0]+fit[1]*fwhm)
		fwhmflag = 1
	endif
	
	;if (shkflag gt 0 OR bisflag gt 0 OR fwhmflag gt 0) then continue
	
	
	planet_periodogram, date-min(date), rv, 2, 0.1, period, power
	peak = 1.0/period[reverse(sort(power))]  ; periods ranked from highest power to lowest
	
	date2 = date[sort(date)]
	baseline = (date2[-1] - date2[0])/365.0  ; time baseline in years
	
	logrhk2 = logrhk  ; save this because ivan's file will overwrite it
	restore,'../../ivan_parameters.dat'
	age = tau[where(hip eq starname)]
	
	printf,u,starname,stdev(rv),n_elements(rv),baseline,age,avg(logrhk2),shkflag,bisflag,fwhmflag,linflag,peak[0],peak[1],peak[2],peak[3],peak[4],format="(i,8x,f8.1,8x,i,8x,f8.1,8x,f8.1,8x,f10.2,8x,i,8x,i,8x,i,8x,i,8x,f8.2,',',f8.2,',',f8.2,',',f8.2,',',f8.2)"
	

	
	
endfor

free_lun,u

END