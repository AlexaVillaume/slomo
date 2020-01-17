# slomo

Python bindings for [Slomo.jl](https://github.com/adwasser/Slomo.jl)


## Dependencies

You'll have to install [pyjulia](https://pyjulia.readthedocs.io/en/latest/) to get Python and Julia to play nice with each other.

## Installation

1. Download and install Julia

	[Link to Julia site](https://julialang.org/)
   
2. Clone this (slomo) repository

	We need to also copy over the `Slomo.jl` submodule by using the `--recurse-submodules` option.
	
	```shell
	git clone --recurse-submodules https://github.com/adwasser/slomo.git
	```

3. Download and install Julia dependencies

	From the shell:
	
	```shell
	cd slomo/slomo/Slomo.jl
	julia
	```
	
	From the Julia REPL:
	
	```julia
	using Pkg
	Pkg.activate(".")
	Pkg.instantiate()
	Pkg.resolve()
	```
	
4. Build a Julia system image for faster load times (optional)

	This will take a *long* time.
	Go make a cup of tea and read a book while you wait.
	The benefit is that the startup time should be *much* faster if you have a precompiled Julia system image.
	
	Run the build script in the Slomo.jl directory:
	
	```shell
	cd slomo/slomo/Slomo.jl
	mkdir build
	julia build.jl
	```
	
5. Install the python package

	From the top-level package directory:
	
	```shell
	pip install .
	```
	
If all goes well, you should be able to do the following from a Python shell:

```python
import numpy as np
from slomo import *
model = JeansModel(NFWModel(), SersicModel(), ConstantBetaModel())
R = np.logspace(-1, 1)
sigma_profile = sigma_los(model, R)
```

The first call of `sigma_los` will be incredibly slow due to the JIT compilation from Julia; subsequent calls should be much faster.
