from matplotlib import rcParams

import corner

rcParams["font.size"] = 16
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"

ndim = 2
nsamples = 10000

import numpy as np
fname = 'tyc3559_mcmc_norv_fixedkernel_v03.csv'
data = np.loadtxt(fname,delimiter=",")

figure = corner.corner(data, labels=[r"$b$",r"$D$ [pct]",r"$T$ [hr]",r"$a_{1c}$ [ppm]",r"$a_{1s}$ [ppm]",r"$a_{2c}$ [ppm]",r"$f_A/f_0$ [pct]",r"$f_{Bb}/f_B$ [ppm]"],quantiles=[.1585,.5,.8415],show_titles=True,title_kwargs={"fontsize": 12})

figure.savefig("demo.png", dpi=300)
