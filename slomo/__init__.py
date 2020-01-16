from julia import Julia

def get_default_sysimage():
    return "/home/asher/.julia/packages/PackageCompiler/4yNnV/sysimg/sys.so"

def get_slomo_env():
    return "/home/asher/work/software/Slomo.jl"

jl = Julia(sysimage=get_default_sysimage())
jl.eval("import Pkg; Pkg.activate(\"{}\")".format(get_slomo_env()))

from julia import Slomo as slomo
from julia.Slomo import *

halos = slomo.Halos

__all__ = ["slomo", "halos", "SersicModel", "ConstantBetaModel", "RSBetaModel",
           "JeansModel", "mass", "density", "density2d", "beta", "update",
           "sigma_los", "sigma_los_parallel", "kappa_los"]

