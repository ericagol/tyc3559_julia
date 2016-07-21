function model_tyc3559_norv(time,param)
# Input parameters (x=param) are:
# x[1] = P  (units of day)
# x[2] = inc = inclination angle
# x[3] = p = R_p/R_* = radius of planet in units of radius of star
# x[4] = t0 = mid-point of transit
# x[5] = u1 = linear limb-darkening coefficient
# x[6] = u2 = quadratic limb-darkening coefficient
# x[7] = f0 = uneclipsed flux
# x[8] = a/R_* = semi-major axis divided by R_*
# x[9] = e = eccentricity
# x[10] = omega = longitude of pericentre
# x[11]= a1c
# x[12]= a1s
# x[13]= a2c
# x[14]= a2s
# x[15]= ratio of flux of 3rd star to total flux
# x[16]= fraction of EB flux due to 2nd star
nt = length(time)
model=zeros(nt)
dt = median(time[2:nt]-time[1:nt-1])
# Photometric data points have positive times, while
# RV data points have negative times:
period=param[1]
t0=param[4]
a1c=param[11]
a1s=param[12]
a2c=param[13]
a2s=param[14]
param2=zeros(10)
for i=1:10
  param2[i]=param[i]
end
param2[3]=1.0/param2[3]
param2[8]=param2[8]*param2[3]
param2[7]=1.0
param2[4]=param2[4]+param2[1]/2.0
param2[10]=param2[10]+180.0 # (for now we're assuming circular, but if we add eccentricity, we'll need this)
param1=zeros(10)
for i=1:10
  param1[i] = param[i]
end
param1[7] = 1.0
tdur = 2*(1+param[3])*period/(2pi*param[8])
#println("T_dur: ",tdur)
#print,'Impact parameter: ',cos(param[2]*pi/180.)*param[8]
nsub = 10
time_phot1=0.0
model_phot2 =0.0
for i=1:nt
  time_phot=time[i]
  if abs(mod(time_phot-t0+period/2,period)-period/2) < (0.51*tdur)
    model_phot1 = transit_orb(time_phot,param1,dt,nsub)
  else
    model_phot1 = 1.0
  end
  model_phot1 += a1c*cos((time_phot-t0)*2.0*pi/period)+a1s*sin((time_phot-t0)*2.0*pi/period)
  model_phot1 += a2c*cos((time_phot-t0)*4.0*pi/period)
  if abs(mod(time_phot-t0,period)-period/2) < (0.51*tdur)
    model_phot2 = transit_orb(time_phot,param2,dt,nsub)
  else
    model_phot2 = 1.0
  end
  model_phot2 += a1s*sin((time_phot-t0)*2.0*pi/period)
  model[i]=((model_phot1*(1.0-param[16])+model_phot2*param[16])* (1.0-param[15])+param[15])*param[7]
end
#if(total(finite(model) eq 0) gt 0) then begin
#   print,'Model is non-finite: ',where(finite(model) eq 0)
#   char=get_kbrd(1)
#endif
return model
end
