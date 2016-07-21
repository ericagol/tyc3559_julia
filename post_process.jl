function hist_param(param,p1,p2,nhist)
np = length(param)
param_grid = zeros(2*nhist)
for i=1:nhist
  param_grid[2*i-1] = p1 + (p2-p1) * (i-1.0)/nhist
  param_grid[2*i  ] = p1 + (p2-p1) * i/nhist
end
param_hist = zeros(2*nhist)
nmiss = 0
for i=1:np
  ibin = floor(Int64,(param[i]-p1)/(p2-p1)*nhist)+1
  if ibin >= 1 && ibin <= nhist
    param_hist[2*ibin-1] += 1
    param_hist[2*ibin  ] += 1
  else
    nmiss +=1
  end
end
return [p1;param_grid;p2],[0;param_hist;0],nmiss
end

function rho_pearson(p1,p2)
@assert(length(p1)==length(p2))
np = length(p1)
rho = 0.0
m1 = mean(p1)
m2 = mean(p2)
s1 = std(p1)
s2 = std(p2)
for i=1:np
  rho += (p1[i]-m1)*(p2[i]-m2)
end
return rho/s1/s2/(np-1)
end

function acf_par(pvec,nstep_corr)
# Computes ACF:
nstep = length(pvec)
mpar = mean(pvec)
spar = std(pvec)
# Now compute ACF:
acf = zeros(nstep_corr)
for k=1:nstep_corr
  for m= 1:nstep-k
    acf[k] += (pvec[m]-mpar)*(pvec[m+k]-mpar)
  end
  acf[k] = acf[k]/spar^2/(nstep-k)
end
return acf
end

# Make a binned light curve:
include("model_tyc3559_v2.jl")
include("transit_orb.jl")
include("occultquad.jl")
include("kepler.jl")

using JLD
using PyPlot

@load "tyc3559_mcmc_norv_fixedkernel_v03.jld"

par_mcmc[:,:,4]=par_mcmc[:,:,4].*(1.0-par_mcmc[:,:,15])
par_mcmc[:,:,13]=par_mcmc[:,:,13].*(1.0-par_mcmc[:,:,15])
pname[4] = "D*(1-f3)"
pname[13] = "a2c*(1-f3)"
# Derived parameters:

# 1). Density of star B:

using CGS
period = (par_mcmc[:,:,2]-par_mcmc[:,:,1])/nlast
tdur = par_mcmc[:,:,5]
b = par_mcmc[:,:,3]
rhoB = 3.0/pi^2 .*period./tdur.^3/(24.*3600.)^2/GRAV.*sqrt(1.0-b.^2)
mrhoB = mean(rhoB); srhoB = std(rhoB)
println("Density of star B: ",mrhoB," +- ",srhoB)

# make histogram:
nhist = 20
p1 = mrhoB-3*srhoB; p2 = mrhoB+3*srhoB
param_grid,param_hist,nmiss = hist_param(rhoB,p1,p2,nhist)

clf()
plot(param_grid,param_hist/maximum(param_hist))
plot([mrhoB-srhoB,mrhoB-srhoB],[0,1])
plot([mrhoB+srhoB,mrhoB+srhoB],[0,1])
plot([mrhoB,mrhoB],[0,1])
xlabel("Density of star B [g/cc]")
ylabel("Relative probability")
  
read(STDIN,Char)

# 2). Doppler mass of brown dwarf:

teff = 6500.0
x = PLANCK * C / (6000e-8)/BOLTZMANN/teff
# Approximate beaming with Blackbody [NEED TO REDO WITH ATMOSPHERE & KEPLER BANDPASS]:
#alpha_beam = (exp(x)*(3.0-x)-3.0)/(exp(x)-1.0)
alpha_beam = -0.898
# Compute radial velocity in km/s:
vr = par_mcmc[:,:,12]/(3.0-alpha_beam)*C/1e5
# Convert to mass of brown dwarf, assuming mass of star:
mB = 1.0
semi = (GRAV * mB * MSUN)^(1//3).* (period.*24*3600/2/pi).^(2//3)
aonr = period./tdur.*sqrt(1.0-b.^2)/pi
inc = acos(b./aonr)
vorb = 2*pi.*semi./(period.*24.*3600.)./1e5
# Mass of brown dwarf:
mBb1 = vr./vorb.*mB.*MSUN./MJUPITER./sin(inc)
mmBb1=mean(mBb1)
smBb1=std(mBb1)
println("Mass of Bb from Doppler: ",mmBb1," +- ",smBb1," M_Jupiter")

# 3). Tidal distortion mass of brown dwarf:
# Equation 2 from Faigler & Mazeh:
q1 = par_mcmc[:,:,6]; q2 = par_mcmc[:,:,7]
u1 = 2.*sqrt(q1).*q2
u2 = sqrt(q1).*(1-2.*q2)
g = 0.3+rand(size(q1))*0.7
#g = 0.65
alpha_ellip = 0.15.*(15.0+u1).*(1.0+g)./(3.0-u1)
mBb2 = mB * MSUN .* aonr.^3./alpha_ellip.* abs(par_mcmc[:,:,13])./MJUPITER./sin(inc).^2
mmBb2 = mean(mBb2)
smBb2 = std(mBb2)
println("Mass of Bb from tides: ",mmBb2," +- ",smBb2," M_Jupiter")
p1 = minimum([mmBb2-3smBb2,mmBb1-3smBb1])
p2 = maximum([mmBb2+3smBb2,mmBb1+3smBb1])
param_grid,param_hist1,nmiss = hist_param(mBb1,p1,p2,nhist)
param_grid,param_hist2,nmiss = hist_param(mBb2,p1,p2,nhist)
clf()
plot(param_grid,param_hist1)
plot(param_grid,param_hist2)
xlabel("Mass of brown dwarf [M_Jupiter]")
ylabel("Relative probability")
read(STDIN,Char)

# 4). Brightness of brown dwarf from phase-curve & secondary eclipse:
fBb_phase = 2.*abs(par_mcmc[:,:,11])*1e6
fBb_occult = abs(par_mcmc[:,:,16])*1e6
mfBb_phase = mean(fBb_phase)
sfBb_phase = std(fBb_phase)
mfBb_occult = mean(fBb_occult)
sfBb_occult = std(fBb_occult)
println("Flux ratio of Bb to B+Bb, phase:  ",mfBb_phase," +- ",sfBb_phase)
println("Flux ratio of Bb to B+Bb, occult: ",mfBb_occult," +- ",sfBb_occult)
p1 = minimum([mfBb_phase-3*sfBb_phase,mfBb_occult-3*sfBb_occult])
p2 = maximum([mfBb_phase+3*sfBb_phase,mfBb_occult+3*sfBb_occult])
param_grid,param_hist1,nmiss = hist_param(fBb_phase,p1,p2,nhist)
param_grid,param_hist2,nmiss = hist_param(fBb_occult,p1,p2,nhist)
clf()
plot(param_grid,param_hist1)
plot(param_grid,param_hist2)
xlabel("fBb/(fB+fBb) [ppm]")
ylabel("Relative probability")
read(STDIN,Char)

# 5). Radius ratio of brown dwarf to star B:
depth = par_mcmc[:,:,4]
x = 1.0-sqrt(1.0-b.^2)
r21 = sqrt(depth.*(1.-u1./3-u2./6)./(1.0-u1.*x-u2.*x.^2))
mr21 = mean(r21); sr21=std(r21)
println("R_Bb/R_B: ",mr21," +- ",sr21)

# 6). Geometric albedo of the brown dwarf from occultation:
ag = fBb_occult.*1e-6.*(aonr./r21).^2
println("Geometric albedo: ",mean(ag)," +- ",std(ag))

# 7). Radius of brown dwarf:
rBb = r21.*(mB*MSUN./rhoB.*3/4/pi).^(1//3)./RJUPITER
mrBb = mean(rBb); srBb = std(rBb)
println("R_Bb: ",mrBb," +- ",srBb," [R_Jupiter]")
read(STDIN,Char)

# # Compute autocorrelation functions:
# nstep_corr = 2000
# acorrlength = zeros(nvary,nwalkers)
# for i=1:nvary
#   clf()
#   acf = zeros(nwalkers,nstep_corr)
#   for j=1:nwalkers
# # First estimate mean & variance:
#     pvec = vec(par_mcmc[j,:,ivary[i]])
#     acf[j,:]=acf_par(pvec,nstep_corr)
#     acorrlength[i,j] = 0.0
#     k=1
#     while (acf[j,k] > 0.5) && (k < nstep_corr-1)
#       k+=1 
#     end
#     acorrlength[i,j] = k
#     plot(vec(acf[j,:]))
#   end
#   println("Autocorrelation length for parameter ",pname[ivary[i]]," is ",maximum(acorrlength[i,:]))
# #  read(STDIN,Char)
# end
# 
# println("Autocorrelation length for chains is: ",maximum(acorrlength))
# println("Effective length of MCMC: ",nwalkers*nsteps/maximum(acorrlength))
# 
# read(STDIN,Char)
clf()

# Compute inter-parameter correlation
for i=1:nvary-1
  for j=i+1:nvary
    p1=vec(par_mcmc[:,:,ivary[i]])
    p2=vec(par_mcmc[:,:,ivary[j]])
    rho = rho_pearson(p1,p2)
    if abs(rho) > 0.5
      println("Parameters ",pname[ivary[i]]," and ",pname[ivary[j]]," correlate strongly  as ",rho)
      scatter(p1,p2)
      read(STDIN,Char)
      clf()
    end
  end
end

clf()

for i=1:nvary
  mp = zeros(nsteps)
  sp = zeros(nsteps)
  for j=1:nwalkers
    plot(vec(par_mcmc[j,1:nsteps,ivary[i]]))
  end
  for j=1:nsteps
    mp[j]=mean(par_mcmc[:,j,ivary[i]])
    sp[j]=std(par_mcmc[:,j,ivary[i]])
  end
  plot(mp+sp)
  plot(mp-sp)
  xlabel("MCMC step")
  ylabel(pname[ivary[i]])
  println("Hit return to continue")
  read(STDIN,Char)
  clf()
end

log_like_best =  -1e100
ibest = 0
jbest = 0
param_best = zeros(nparam)
for i=1:nsteps
  for j=1:nwalkers
    if log_like_mcmc[j,i] > log_like_best
      ibest = i
      jbest = j
      log_like_best = log_like_mcmc[j,i]
      param_best = vec(par_mcmc[j,i,:])
    end
  end
end

model = model_tyc3559_v2(time_phot,param_best)

# Compute residuals:

y = flux - model

period = (param_best[2]-param_best[1])/nlast

t0 = param_best[1]

nbin = 128
phasebin_resid = zeros(nbin)
phasebin_resid_med = zeros(nbin)
phasebin_raw = zeros(nbin)
phasebin_raw_med = zeros(nbin)
phasebin_model = zeros(nbin)
phase0 = mod((time_phot-t0)+period/4,period)/period
is = sortperm(phase0)
j=1
for ibin=1:nbin
  i=j
  while phase0[is[j+1]] < (ibin/nbin) && (j+1) < ntime
    j=j+1
  end
  phasebin_resid[ibin] = mean(y[is[i:j]])
  phasebin_resid_med[ibin] = median(y[is[i:j]])
  phasebin_raw[ibin] = mean(flux[is[i:j]])
  phasebin_raw_med[ibin] = median(flux[is[i:j]])
  phasebin_model[ibin] = mean(model[is[i:j]])
end

phase = zeros(nbin)
for i=1:nbin
  phase[i] = (i-0.5)/nbin
end

clf()
#scatter(phase,phasebin_resid)
scatter(phase,phasebin_raw)
scatter(phase,phasebin_raw_med)
plot(phase,phasebin_model)
plot(phase,phasebin_model)
model_reflect=((1+param_best[11]*cos((phase-0.25)*2pi))*(1.0-param_best[15])+param_best[15])*param_best[8]
model_beam   =((1+param_best[12]*sin((phase-0.25)*2pi))*(1.0-param_best[15])+param_best[15])*param_best[8]
model_ellip  =((1+param_best[13]*cos((phase-0.25)*4pi))*(1.0-param_best[15])+param_best[15])*param_best[8]
#  model[i]=((model_phot1*(1.0-param[16])+model_phot2*param[16])* (1.0-param[15])+param[15])*f0

plot(phase,model_reflect)
plot(phase,model_beam)
plot(phase,model_ellip)
