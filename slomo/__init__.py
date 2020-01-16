from julia import Julia

def get_sysimage():
    return "/home/asher/work/software/Slomo.jl/build/sys.so"

def get_env():
    return "/home/asher/work/software/Slomo.jl"

jl = Julia(sysimage=get_sysimage())
jl.eval("import Pkg; Pkg.activate(\"{}\")".format(get_env()))

from julia import Slomo as slomo
from julia.Slomo import *

halos = slomo.Halos

__all__ = ["slomo", "halos", "SersicModel", "ConstantBetaModel", "RSBetaModel",
           "JeansModel", "mass", "density", "density2d", "beta", "update",
           "sigma_los", "sigma_los_parallel", "kappa_los"]

