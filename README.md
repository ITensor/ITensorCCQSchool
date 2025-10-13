# ITensor CCQ School Hands-On Tutorials

Welcome to the ITensor CCQ School hands-on tutorials.

To run the tutorials:

1. Download and install the latest release (v1.12) of Julia following the official
instructions here: https://julialang.org/install/
2. Start the Julia REPL by executing the `julia` command, which should now be available on
your computer if you followed the installation instructions in step 1.:
```
$ julia
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.12.0 (2025-10-07)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org release
|__/                   |


julia>
```
3. Create a local copy of the tutorial code in the directory `ITensorCCQSchool/` in you
current directory by running:
```julia
julia> using LibGit2: LibGit2

julia> LibGit2.clone("https://github.com/ITensor/ITensorCCQSchool", ".")
```
4. Enter the `ITensorCCQSchool/` directory and install the dependencies from the Julia REPL:
```julia
julia> cd("ITensorCCQSchool/")

julia> using Pkg: Pkg

julia> Pkg.activate(".")
  Activating project at `[...]/ITensorCCQSchool`

julia> Pkg.instantiate()
    Updating registry at `~/.julia/registries/General.toml`
    Updating `[...]/ITensorCCQSchool/Project.toml`
  [0d1a4710] + ITensorMPS v0.3.22
  [9136182c] + ITensors v0.9.13
    Updating `[...]/ITensorCCQSchool/Manifest.toml`
  [7d9f7c33] + Accessors v0.1.42
  [79e6a3ab] + Adapt v4.4.0
  [dce04be8] + ArgCheck v2.5.0
  [...]
```
5. Use `include` to run the
[first tutorial](https://github.com/ITensor/ITensorCCQSchool/blob/main/day1/01-julia-intro.jl)
from the REPL:
```julia
julia> include("day1/01-julia-intro.jl");
n = 100000
pi_approx = 3.1416026534897203
error = 9.999899927226608e-6
```
The `;` at the end of the line suppresses printing the output of the script, to avoid
getting a potentially large output to your terminal. To access the objects that are
returned from the `let` block in the script so you can analyze them interactively, you
can call the script like this:
```julia
julia> res = include("day1/01-julia-intro.jl");
n = 100000
pi_approx = 3.1416026534897203
error = 9.999899927226608e-6
```
`res = ` captures the output of the script to the
[NamedTuple](https://docs.julialang.org/en/v1/base/base/#Core.NamedTuple) `res`, which you
can think of as an anonymous struct. You can access values from `res` as follows:
```julia
julia> res.pi_approx
3.1416026534897203

julia> res.error
9.999899927226608e-6

julia> (; pi_approx, error) = res

julia> pi_approx
3.1416026534897203

julia> error
9.999899927226608e-6
```
6. Edit the file to change parameters, printing, etc.:
```julia
julia> edit("day1/01-julia-intro.jl")
```
which will open your file in a text editor determined by Julia (see the
[documentation for `edit`](https://docs.julialang.org/en/v1/stdlib/InteractiveUtils/#InteractiveUtils.edit-Tuple{AbstractString,%20Integer})
for more details). Otherwise, open the file with your text editor or IDE of choice, such as
Vim, Emacs, VS Code, etc. Simply save the file and rerun
`res = include("day1/01-julia-intro.jl");` in your existing Julia session to see your changes
reflected in the output. Try increasing or decreasing the number of terms in the series for `pi`
to see the error decrease or increase. How does the error converge with the number of terms
in the series?
