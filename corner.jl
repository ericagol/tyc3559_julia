using JLD
using PyCall
@pyimport corner
@pyimport numpy as py
@pyimport matplotlib.pyplot as plt

@load "tyc3559_mcmc_norv_fixedkernel_v03.jld"

#nsamp = 10000
nsamp = nwalkers*nsteps
npar = 8
#npar = 2
data = zeros(nsamp,npar)
ipar = [3,4,5,11,12,13,15,16]
#ipar = [3,4]
for i1=1:nwalkers
  for i2=1:nsteps
#  i1 = ceil(Int64,rand()*nwalkers)
#  i2 = ceil(Int64,rand()*nsteps)
    for j=1:npar
      data[(i1-1)*nsteps+i2,j]=par_mcmc[i1,i2,ipar[j]]
    end
  end
end

data[:,2] *= 100.
data[:,3] *= 24.
data[:,4] *= 1e6
data[:,5] *= 1e6
data[:,6] *= 1e6
data[:,7] *= 100.
data[:,8] *= 1e6

writedlm("tyc3559_mcmc_norv_fixedkernel_v03.csv",data,",")

#labname=[r"$b$",r"$D$",r"$T$",r"$a_{1c}$",r"$a_{1s}$",r"$a_{2c}$",r"$f_A/(f_A+f_B+f_{Bb})$",r"$f_{Bb}/(f_B+f_{Bb})$"]

#figure = corner.corner(data,labels=[r"\$b\$","$D$","$T$","$a_{1c}$","$a_{1s}$","$a_{2c}$","$f_A/(f_A+f_B+f_{Bb})$","$f_{Bb}/(f_B+f_{Bb})$"], quantiles=[0.1585, 0.5, 0.8415])
#figure = corner.corner(data, labels=PyAny([r"$b$",r"$D$"]), quantiles=[0.1585, 0.5, 0.8415])
#                       quantiles=[0.1585, 0.5, 0.8415],
#                       show_titles=True, title_kwargs={"fontsize": 12})

#figure = pycall(corner.corner,PyAny, data, labels=[r"$b$",r"$D$"], quantiles=[0.1585, 0.5, 0.8415])
#figure = pycall(corner.corner,PyAny, data, labels=['$b$','$D$'], quantiles=[0.1585, 0.5, 0.8415])
#plt.show(figure)
