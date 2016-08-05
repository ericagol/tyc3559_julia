# Takes full dataset and thins by a factor of 100:

using JLD
@load "tyc3559_mcmc_norv_fixedkernel_v03.jld"

nthin = 100
par_mcmc_thin = zeros(nwalkers,nthin,nparam)
log_like_mcmc_thin =  zeros(nwalkers,nthin)

for i=1:nwalkers
  for j=1:nthin
    par_mcmc_thin[i,j,:]=par_mcmc[i,(j-1)*100+1,:]
    log_like_mcmc_thin[i,j,:]=log_like_mcmc[i,(j-1)*100+1]
  end
end

save("tyc3559_mcmc_norv_fixedkernel_v03_thin.jld","par_mcmc",par_mcmc_thin,"log_like_mcmc",log_like_mcmc_thin,"nwalkers",nwalkers,"nsteps",nthin,"nparam",nparam,"nvary",nvary,"ivary",ivary,"time_phot",time_phot,"flux",flux,"nlast",nlast)
