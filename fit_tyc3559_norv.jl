using PyPlot
using JLD
include("model_tyc3559_norv.jl")
include("model_tyc3559_v2.jl")
include("transit_orb.jl")
include("occultquad.jl")
include("kepler.jl")
include("compute_likelihood_save.jl")
include("lorentz_likelihood_hermitian_band_init.jl")
include("lorentz_likelihood_hermitian_band_save.jl")
include("bandec.jl")
include("banbks.jl")

#function fit_tyc3559_norv()

data = readdlm("tyc3559_pdc.dat")
time_phot  = vec(data[:,1])-54900.0
fflat = vec(data[:,2])
ntphot = 65266
#ntphot = 100
#time_phot = time_phot[1:ntphot]
#fflat = fflat[1:ntphot]

#param=[pbest,87.0,.023.0,.857.0*pbest,.5,0.0,1.0+6e-5,5.0,0.0,0.0,0e-5,0e-5,-1.3e-4,0.0,2450.0,-1334.0]
param_phot = vec(readdlm("param_best_norv.dat"))
pname1 = ["Period [day]","Inc [deg]","R_2/R_1","t0 [JD-2459000]","u1","u2","f0","a/R_2","ecc","omega","a1c","a1s","a2c","a2s","f_3/(f_1+f_2)","f_2/(f_1+f_2)"]
param_phot1=zeros(16)
for i=1:16
  param_phot1[i]=param_phot[i]
end

ntime=length(time_phot)
t_first = minimum(time_phot)
t_last = maximum(time_phot)
p0 = param_phot[1]
t00 = param_phot[4]
nfirst = floor((t_first-t00)/p0)
nlast = ceil((t_last-t00)/p0)
tn = nlast*p0+t00
aonr = param_phot[8]
u1 = param_phot[5]
u2 = param_phot[6]
r2onr1 = param_phot[3]
inc = param_phot[2]*pi/180.0
b = cos(inc)*aonr
tdur = 1.0/aonr*p0/pi*sqrt(1-b^2)
x = (1-sqrt(1-b^2))
depth = (r2onr1)^2*(1.0-u1*x-u2*x^2)/(1.0-u1/3-u2/6)
q1 = (u1+u2)^2
q2 = u1/(2*(u1+u2))
f0 = param_phot[7]

# Convert to new set of variables:
pname = ["t0 [JD - 24590000]","tN [JD-2459000]","b","depth","T_dur","q1","q2","f0","ecc","omega","a1c","a1s","a2c","a2s","f_3/(f_1+f_2)","f_2/(f_1+f_2)"]
param_phot[1] = t00
param_phot[2] = tn
param_phot[3] = b
param_phot[4] = depth
param_phot[5] = tdur
param_phot[6] = q1
param_phot[7] = q2
param_phot[8] = f0

#time_phot=reform(phasebin[0,*])
dt = median(time_phot[2:ntime]-time_phot[1:ntime-1])
flux_phot=fflat[1:ntime]

scattermin = 0.00040507201
weight_phot=ones(ntime)./scattermin^2
weight=weight_phot
flux=flux_phot
model1 = model_tyc3559_norv(time_phot,param_phot1)
model2 = model_tyc3559_v2(time_phot,param_phot)
@time model = model_tyc3559_v2(time_phot,param_phot)
chisquare=sum((flux-model).^2.*weight)
println("Chi-square: ",chisquare)
println("length(model1): ",length(model1))
println("length(model2): ",length(model2))
println("max(abs(model1-model2)): ",maximum(abs(model1-model2)))
plot(mod(time_phot,p0),model1)
plot(mod(time_phot,p0),model2)
#read(STDIN,Char)

# Now compute correlated likelihood:
w = 0.03027 * 1.6467e-7
kernel_param = [w,1.0428542*1.6467e-7, -0.38361831*1.6467e-7, 0.30345984/2*1.6467e-7, 0.1229159/dt, 0.48922908/dt, 0.09086397/dt, 2pi/12.203317/dt]
# Set kernel to a delta function:
#kernel_param = [w,1.0428542*1.6467e-7, -0.38361831*1.6467e-7, 0.30345984/2*1.6467e-7, 1e5 , 1e5 , 1e5 , 0.]
param_tot = [param_phot;kernel_param]
w = param_tot[17]
alpha = [param_tot[18],param_tot[19],param_tot[20]/2,param_tot[20]/2]
p=length(alpha)
beta = [complex(param_tot[21],0.0),complex(param_tot[22],0.0),complex(param_tot[23],param_tot[24]),complex(param_tot[23],-param_tot[24])]
aex,al_small,indx,logdeta = lorentz_likelihood_hermitian_band_init(alpha,beta,w,time_phot)
#Profile.clear()
#for i=1:100
#  @profile log_like = compute_likelihood_save(time_phot,param_tot,flux,aex,al_small,indx,logdeta,p)
#end
#Profile.print()
#read(STDIN,Char)
log_like = compute_likelihood_save(time_phot,param_tot,flux,aex,al_small,indx,logdeta,p)

nparam = 24
ivary = [1,2,3,4,5,6,7,8,11,12,13,15,16]
nvary = length(ivary)

# Run a Markov chain:
#errors = [0.0001,0.2,0.001,0.01,0.001,0.1,1e-6,0.1,0.,0.,1e-7,1e-7,1e-7,1e-7,1e-3,1e-7,0.,0.,0.,0.,0.,0.,0.,0.]
errors = [0.001,0.001,0.01,0.00001,0.01,0.1,0.1,1e-6,0.,0.,1e-7,1e-7,1e-7,0.,1e-2,1e-7,0.,0.,0.,0.,0.,0.,0.,0.]
#pname = ["t0 [JD - 24590000]","tN [JD-2459000]","b","depth","T_dur","q1","q2","f0","ecc","omega","a1c","a1s","a2c","a2s","f_3/(f_1+f_2)","f_2/(f_1+f_2)"]

nwalkers = nvary * 2
nsteps = 10000
#nsteps = 100
# Set up arrays to hold the results:
par_mcmc = zeros(nwalkers,nsteps,nparam)
log_like_mcmc = zeros(nwalkers,nsteps)
#logprior_mcmc = zeros(nwalkers,nsteps)
# Initialize walkers:
par_trial = param_tot
log_like_best = log_like
for j=1:nwalkers
# Select from within uncertainties:
  log_like_trial = -1e100
# Only initiate models with reasonable chi-square values:
  while log_like_trial < log_like_best - 1000.0
    par_trial = param_tot + errors.*randn(nparam)
    aonr2 = (1./par_trial[5]*p0/pi)^2*(1-par_trial[3]^2)
    if (par_trial[6] < 0) || (par_trial[6] > 1) || (par_trial[7] < 0) || (par_trial[7] > 1) || (par_trial[4] < 0) || (par_trial[5] < 0) || (abs(par_trial[3]) > 1) || (par_trial[3]^2 > aonr2) then
      log_like_trial = -1e100
   else
      log_like_trial = compute_likelihood_save(time_phot,par_trial,flux,aex,al_small,indx,logdeta,p)
   end
#    println("log_like_trial: ",log_like_trial," params: ",par_trial)
  end
  log_like_mcmc[j,1]=log_like_trial
  par_mcmc[j,1,:]=par_trial
#  logprior_mcmc[j,1] = log_prior+log_det
  println("Success: ",par_trial,log_like_trial)
end

# Initialize scale length & acceptance counter:
ascale = 2.0
accept = 0
tic()
# Next, loop over steps in markov chain:
for i=2:nsteps
  for j=1:nwalkers
    ipartner = j
# Choose another walker to 'breed' a trial step with:
    while ipartner == j
      ipartner = ceil(Int,rand()*nwalkers)
    end
# Now, choose a new set of parameters for the trial step:
    z=(rand()*(sqrt(ascale)-1.0/sqrt(ascale))+1.0/sqrt(ascale))^2
    par_trial = vec(par_mcmc[j,i-1,:])
    for itrial=1:nvary
      par_trial[ivary[itrial]]=z*par_mcmc[j,i-1,ivary[itrial]]+(1.0-z)*par_mcmc[ipartner,i-1,ivary[itrial]]
    end
# Compute model & chi-square:
#    log_prior,log_det,log_like_trial = compute_likelihood(time_phot,par_trial,flux)
#    log_like_trial = compute_likelihood(time_phot,par_trial,flux)
    aonr2 = (1./par_trial[5]*p0/pi)^2*(1-par_trial[3]^2)
    if (par_trial[6] < 0) || (par_trial[6] > 1) || (par_trial[7] < 0) || (par_trial[7] > 1) || (par_trial[4] < 0) || (par_trial[5] < 0) || (abs(par_trial[3]) > 1) || (par_trial[3]^2 > aonr2) then
      log_like_trial = -1e100
    else
      log_like_trial = compute_likelihood_save(time_phot,par_trial,flux,aex,al_small,indx,logdeta,p)
    end
# Next, determine whether to accept this trial step:
    alp = z^(nvary-1)*exp(log_like_trial - log_like_mcmc[j,i-1])
#    if rand() < 0.0001
#      println("Step: ",i," Walker: ",j," Chi-square: ",log_like_trial," Prob: ",alp," Frac: ",accept/(mod(i-1,1000)*nwalkers+j))
#    end
    if alp >= rand()
# If step is accepted, add it to the chains!
      par_mcmc[j,i,:] = par_trial
      log_like_mcmc[j,i,:] = log_like_trial
#      logprior_mcmc[j,i,:] = log_prior + log_det
      accept = accept + 1
    else
# If step is rejected, then copy last step:
      par_mcmc[j,i,:] = par_mcmc[j,i-1,:]
      log_like_mcmc[j,i,:] = log_like_mcmc[j,i-1]
#      logprior_mcmc[j,i,:] = logprior_mcmc[j,i-1,:]
    end
  end
  if mod(i,100) == 0
    eltime = toq()
    println("Number of steps: ",i," acceptance rate: ",accept/(100*nwalkers)," elapsed time [hr]: ",eltime/3600.)
    accept = 0
    tic()
  end
end

for i=1:nvary
  for j=1:nwalkers
    plot(vec(par_mcmc[j,1:nsteps,ivary[i]]))
  end
  xlabel("MCMC step")
  ylabel(pname[ivary[i]])
  println("Hit return to continue")
  read(STDIN,Char)
  clf()
end

@save "tyc3559_mcmc_norv_fixedkernel_v02.jld"
