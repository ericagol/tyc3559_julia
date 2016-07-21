using PyPlot
using JLD

include("model_tyc3559_v2.jl")
include("transit_orb.jl")
include("occultquad.jl")
include("kepler.jl")
include("compute_likelihood_save.jl")
include("lorentz_likelihood_hermitian_band_save.jl")
include("banbks.jl")

@load "tyc3559_mcmc_norv_fixedkernel_v02.jld"

# Initialize walkers from last step of prior run:
for j=1:nwalkers
  log_like_mcmc[j,1]=log_like_mcmc[j,nsteps]
  par_mcmc[j,1,:]=par_mcmc[j,nsteps,:]
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
    aonr2 = (1./par_trial[5]*p0/pi)^2*(1-par_trial[3]^2)
    if (par_trial[6] < 0) || (par_trial[6] > 1) || (par_trial[7] < 0) || (par_trial[7] > 1) || (par_trial[4] < 0) || (par_trial[5] < 0) || (abs(par_trial[3]) > 1) || (par_trial[3]^2 > aonr2) then
      log_like_trial = -1e100
    else
      log_like_trial = compute_likelihood_save(time_phot,par_trial,flux,aex,al_small,indx,logdeta,p)
    end
# Next, determine whether to accept this trial step:
    alp = z^(nvary-1)*exp(log_like_trial - log_like_mcmc[j,i-1])
    if alp >= rand()
# If step is accepted, add it to the chains!
      par_mcmc[j,i,:] = par_trial
      log_like_mcmc[j,i,:] = log_like_trial
      accept = accept + 1
    else
# If step is rejected, then copy last step:
      par_mcmc[j,i,:] = par_mcmc[j,i-1,:]
      log_like_mcmc[j,i,:] = log_like_mcmc[j,i-1]
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

@save "tyc3559_mcmc_norv_fixedkernel_v03.jld"

