FUNCTION calc_rhk, spec_file, rv, bv, Teff, shk = shk, jitter = jitter
;calculates the stellar activity index given a HARPS 1D spectrum file, its RV in km/s, and the star's B-V color and T_effective
;uses procedure from Lovis et al. 2011
;optional keyword shk returns S_hk index
;optional keyword jitter returns Wright 2005's median jitter prediction (OPTIMIZED FOR SOLAR-TYPE B-V)
;otherwise, R'_hk is returned

if (n_params() ne 4) then begin
 print, 'Rprime_HK = calc_Rhk(spec_file, RV [km/s], B-V, Teff, shk = shk, jitter = jitter)'
 retall
endif

;retrieve spectrum and reconstruct wavelength scale:
spec = mrdfits(spec_file,/silent)
header = headfits(spec_file)

wave = findgen(n_elements(spec))
start = sxpar(header,'CRVAL1')
slope = sxpar(header,'CDELT1')
wave = start + wave*slope

wave = wave/(1 + rv/299792.5)  ; shift to star's rest frame

; trim unneeded parts of the spectrum (this helps avoid NaN errors from weird points at edges):
keep = where(wave gt 3850.0 and wave lt 4050.0)
spec = spec[keep]
wave = wave[keep]

; mask out negative points:
neg = where(spec lt 0.0,count)
if (count gt 0) then spec[neg] = 0.0

;construct bandpass "filters":
Kcenter = 3933.664
Kfwhm = 1.09
Kfilter = 1. - (1./Kfwhm) * abs(wave - Kcenter)
Kfilter[where(Kfilter lt 0)] = 0.0

Hcenter = 3968.470
Hfwhm = 1.09
Hfilter = 1. - (1./Hfwhm) * abs(wave - Hcenter)
Hfilter[where(Hfilter lt 0)] = 0.0

Vcenter = 3901.070
Vfwhm = 20.0
Vfilter = 1. - (1./Vfwhm) * abs(wave - Vcenter)
Vfilter[where(Vfilter lt 0)] = 0.0

Rcenter = 4001.070
Rfwhm = 20.0
Rfilter = 1. - (1./Rfwhm) * abs(wave - Rcenter)
Rfilter[where(Rfilter lt 0)] = 0.0

;get flux & photon error in each bandpass:
Kflux = total(spec * Kfilter)
Kerr = sqrt(total((sqrt(spec) * Kfilter)^2))
Hflux = total(spec * Hfilter)
Herr = sqrt(total((sqrt(spec) * Hfilter)^2))
Vflux = total(spec * Vfilter)
Verr = sqrt(total((sqrt(spec) * Vfilter)^2))
Rflux = total(spec * Rfilter)
Rerr = sqrt(total((sqrt(spec) * Rfilter)^2))


S_HARPS = (Hflux/Hfwhm + Kflux/Kfwhm)/(Rflux/Rfwhm + Vflux/Vfwhm)  ;construct S_HARPS from mean fluxes
err_S_HARPS = S_HARPS * sqrt(((Kerr/Kfwhm)^2 + (Herr/Hfwhm)^2)/(Hflux/Hfwhm + Kflux/Kfwhm)^2 + ((Rerr/Rfwhm)^2 + (Verr/Vfwhm)^2)/(Rflux/Rfwhm + Vflux/Vfwhm)^2)
S_MW = 1.111*S_HARPS + 0.0153  ;conversion to Mount Wilson value
err_S_MW = 1.111*err_S_HARPS


logC_cf = 1.13*bv^3 - 3.91*bv^2 + 2.84*bv - 0.47
;logR_phot = -4.02 - 1.40*bv
logR_phot = -4.898 + 1.918*bv^2 -2.893*bv^3

Rprime_hk = (1.340 * 10.^(-4)) * 10.^(logC_cf) * S_MW - 10.^(logR_phot)

;Rutten's index & jitter prediction from Wright 2005:
logC_cf2 = 0.25*bv^3 - 1.33*bv^2 + 0.43*bv + 0.24
F_CaII = S_MW * 10.^(logC_cf2) * Teff^4 * 10.^(-14)
F_base = 1.41  ; bv ~ 0.65
dF = F_CaII - F_base
sig_Wright = 0.83 + 0.96*dF + 0.61*dF^2  ; median prediction


if keyword_set(shk) then return, [S_MW,err_S_MW]
if keyword_set(jitter) then return, sig_Wright
return, Rprime_hk



END
