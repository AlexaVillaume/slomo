from os import path
from julia import Julia

def get_sysimage():
    here = path.abspath(path.dirname(__file__))
    built_sysimage = path.join(here, "Slomo.jl", "build", "sys.so")
    if path.isfile(built_sysimage):
        return built_sysimage
    return None

def get_env():
    here = path.abspath(path.dirname(__file__))
    env = path.join(here, "Slomo.jl")
    return env

jl = Julia(sysimage=get_sysimage())
jl.eval("import Pkg; Pkg.activate(\"{}\")".format(get_env()))

def get_nthreads():
    return jl.eval("Threads.nthreads()")

from julia import Slomo as slomo
from julia.Slomo import *

halos = slomo.Halos

__all__ = ["slomo", "halos", "SersicModel", "ConstantBetaModel", "RSBetaModel",
           "JeansModel", "mass", "density", "density2d", "beta", "update",
           "sigma_los", "sigma_los_parallel", "kappa_los"]

