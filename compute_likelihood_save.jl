function compute_likelihood_save(time_phot,param,data,aex,al_small,indx,logdeta,p)

# First call prior:

#log_prior = compute_prior(param)
log_prior = 0.0

# Then call model:

model = model_tyc3559_v2(time_phot,param[1:16])

# Then compute log-likelihood:

# Subtract off the model from the data:
y = data - model

log_like = lorentz_likelihood_hermitian_band_save(p,y,aex,al_small,indx,logdeta)
# Apply prior on flux ratio:
log_like -= 0.5*(param[15]-0.692)^2/0.041^2

return log_like
end
