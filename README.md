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
3. Clone this Github repository to a directory on your computer by running:
```julia
julia> using LibGit2: LibGit2

julia> LibGit2.clone("https://github.com/ITensor/ITensorCCQSchool", ".")
```
4. Enter the directory and install the dependencies from the Julia REPL:
```julia
julia> cd("ITensorCCQSchool")

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
5. Use `include` to run the first example from the REPL:
```julia
julia> include("day1/01-julia-intro.jl");
n = 100000
pi_approx = 3.1416026534897203
error = 9.999899927226608e-6
```
The `;` at the end of the line suppresses printing the output of the script, to avoid
getting a potentially large output to your terminal. To access the objects that are
returned from the `let` block in the script, you can call the script like this:
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
```
6. Edit the files to change parameters, measurements, etc.:
```julia
julia> edit("day1/01-julia-intro.jl")
```
which will open your file in a text editor determined by Julia (see
https://docs.julialang.org/en/v1/stdlib/InteractiveUtils/#InteractiveUtils.edit-Tuple{AbstractString,%20Integer}). Otherwise, open the file with your text editor or IDE of choice, such as
Vim, Emacs,VS Code, etc. Simply save the file and run
`res = include("day1/01-julia-intro.jl");` in the Julia REPL again to see your changes
reflected in the output.
