data=readdlm("MIST_iso_1470351180.iso")

using CGS
rho_sun = MSUN/(4pi/3.*RSUN^3)
# Convert to density:
mass = vec(data[:,4])
density = mass./10.^(3.*data[:,12])*rho_sun
teff = vec(10.^data[:,11])

using PyPlot

fig,axes = subplots(2,1)
ax = axes[1,1]
ax[:plot](density,teff,"r-",alpha=0.6,label="t = 1.5 Gyr, [Fe/H]=0.0, MIST",linewidth=3)

rhoA=0.001167
srhoA=0.000052
tA = 4750.
stA = 250.
ax[:plot](rhoA+srhoA*[-1,1],tA*[1,1],"g-",linewidth=3,alpha=0.6)
ax[:plot](rhoA*[1,1],tA+stA*[-1,1],"g-",linewidth=3,alpha=0.6,label="Star A")


rhoB=0.0824
srhoB=0.0034
tB = 6500.
stB = 500.
ax[:semilogx](rhoB+srhoB*[-1,1],tB*[1,1],"b-",linewidth=3,alpha=0.6)
ax[:semilogx](rhoB*[1,1],tB+stB*[-1,1],"b-",linewidth=3,alpha=0.6,label="Star B")
ax[:axis]([1e-4,1e0,4e3,8e3])
#ax[:set_xlabel]("Density [g/cc]")
ax[:set_ylabel]("Teff [K]")
ax[:legend](loc = "upper left")

ax = axes[2,1]

ax[:semilogx](density,mass,"r-",alpha=0.6,label="t = 1.5 Gyr, [Fe/H]=0.0, MIST",linewidth=3)
mA = 1.9
smA = 0.27
ax[:semilogx](rhoA+srhoA*[-1,1],mA*[1,1],"g-",linewidth=3,alpha=0.6)
ax[:semilogx](rhoA*[1,1],mA+smA*[-1,1],"g-",linewidth=3,alpha=0.6,label="Star A")
ax[:semilogx](rhoB+srhoB*[-1,1,1,-1,-1],[0,0,10,10,0],"b-",linewidth=3,alpha=0.6,label="Star B")
ax[:axis]([1e-4,1e0,1.7,2.0])
ax[:set_xlabel]("Density [g/cc]")
ax[:set_ylabel]("Mass [Msun]")

PyPlot.savefig("isochrone.png",bbox_inches="tight")

