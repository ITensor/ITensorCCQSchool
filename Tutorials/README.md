# ITensor CCQ School Hands-On Tutorials Introduction

Welcome to the ITensor CCQ School hands-on tutorials. Here you will find instructions on how to install Julia and get started running the hands-on tutorials. The tutorials for each day can be found in the corresponding subdirectories [Day1](./Day1/), [Day2](./Day2), and [Day3](./Day3).

## Installation Instructions

To run the tutorials:

1. Download and install the latest release (v1.12) of Julia following the official
instructions here: https://julialang.org/install/

2. Start the Julia REPL by executing the `julia` command, which should now be available
on your computer if you followed the installation instructions in step 1.:
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


julia> 1 + 1
2

```
Try typing a command (such as `1 + 1` shown above) to get a feel for how it works. A number of math operations are available out-of-the-box, such as `sin`, `cos`, etc., while other functionality (such as [linear algebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/)) requires loading packages. The surface level syntax is similar to other high level interactive languages like Python and MATLAB. The Julia documentation provides a helpful guide [comparing Julia to other comparable languages](https://docs.julialang.org/en/v1/manual/noteworthy-differences/).

3. Create a local copy of the tutorial code in a new directory `ITensorCCQSchool` in your current directory by running:
```julia
julia> using LibGit2: clone

julia> clone("https://github.com/ITensor/ITensorCCQSchool", ".")
```
Here we use Julia's  [LibGit2 standard library](https://docs.julialang.org/en/v1/stdlib/LibGit2/). Alternatively you can execute `git clone https://github.com/ITensor/ITensorCCQSchool` directly from a shell.

4. Now that you have Julia installed and the tutorial code available, we will give an introduction to running the first tutorial for day 1 ([Day1/1-julia-intro.jl](./Day1/1-julia-intro.jl)). Enter the `ITensorCCQSchool/Tutorials/Day1` directory using Julia's [`cd`](https://docs.julialang.org/en/v1/base/file/#Base.Filesystem.cd-Tuple{AbstractString}) function and install the dependencies from the Julia REPL:
```julia
julia> cd("ITensorCCQSchool/Tutorials/Day1")

julia> ]

(Day1) pkg> activate .
  Activating project at `[...]/ITensorCCQSchool/Tutorials/Day1/`

(Day1) pkg> instantiate
    Updating registry at `~/.julia/registries/General.toml`
    Updating `[...]/ITensorCCQSchool/Tutorials/Day1/Project.toml`
  [0d1a4710] + ITensorMPS v0.3.22
  [9136182c] + ITensors v0.9.13
  [...]

```
Executing `]` at the REPL enables the Pkg REPL, which is more convenient for entering Pkg commands. Press delete/backspace to exit the Pkg REPL and go back to the standard Julia REPL prompt. `activate .` enables the local environment/project in [Tutorials/Day1](./Day1/), where the package dependencies for the tutorials on the first day of the school are defined. `instantiate` installs those dependencies and performs some compilation. It may take some time but it will only need to be done once for each project (so in our case, once for each day of the school).

5. Use `include` to load the [first tutorial](./Day1/1-julia-intro.jl) into the REPL. That will introduce the function `main` which you can execute to run the tutorial:
```julia
julia> include("1-julia-intro.jl")
main (generic function with 1 method)

julia> main();
nterm = 100000
pi_approx = 3.1416026534897203
error = 9.999899927226608e-6

```
The script is approximating the value of `pi` using the [Leibniz formula](https://en.wikipedia.org/wiki/Leibniz_formula_for_%CF%80) using `nterm` terms in the series. To access the values that are returned from the `main` function so you can analyze them interactively, you can call the script like this:
```julia
julia> res = main();
nterm = 100000
pi_approx = 3.1416026534897203
error = 9.999899927226608e-6

julia> res
(nterm = 100000, pi_approx = 3.1416026534897203, error = 9.999899927226608e-6)

```
`res = main()` captures the output of the script to the [NamedTuple](https://docs.julialang.org/en/v1/base/base/#Core.NamedTuple) `res`, which you can think of as an anonymous struct. You can access values from `res` as follows:
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

6. You can run the script with different parameters as follows:
```julia
julia> res = main(; nterm=10^7);
Number of terms: 10000000
Approximate pi: 3.1415927535897814
Error: 9.999998829002266e-8

julia> res.pi_approx
3.1415927535897814

julia> pi # Julia's built-in definition of pi
π = 3.1415926535897...

```
and you can learn about other parameters by printing the documentation using Julia's [help mode](https://docs.julialang.org/en/v1/stdlib/REPL/#Help-mode):
```julia
julia> ?

help?> main
search: main @main min Main Pair map join in mark asin wait sin max tan map!

  main(; kwargs...)

  Approximate pi using the Leibniz formula and compute the error.

  Keywords
  ≡≡≡≡≡≡≡≡

    •  nterm::Int = 10^5: The number of terms to use in the approximation of pi.
    •  outputlevel::Int = 1: Controls the verbosity of the output.

  Outputs
  ≡≡≡≡≡≡≡

    •  nterm::Int: The number of terms used in the approximation of pi.
    •  pi_approx::Float64: The approximate value of pi.
    •  error::Float64: The absolute error in the approximation of pi.

julia>

```

7. Note that you can analyze which directory you are in and what tutorial files are available directly from the Julia REPL using functions such as [`pwd`](https://docs.julialang.org/en/v1/base/file/#Base.Filesystem.pwd) and [`readdir`](https://docs.julialang.org/en/v1/base/file/#Base.Filesystem.readdir):
```julia
julia> pwd()
"[...]/ITensorCCQSchool/Tutorials/Day1"

julia> readdir()
[...]-element Vector{String}:
 "1-julia-intro.jl"
 "2-dmrg.jl"
 "3-dmrg-measure.jl"
[...]

```
Of course, you can also use your terminal and/or file manager as usual. Additionally, if you want to modify the script itself, one way to do that is:
```julia
julia> edit("1-julia-intro.jl")

```
which will open your file in a text editor determined by Julia (see the [documentation for `edit`](https://docs.julialang.org/en/v1/stdlib/InteractiveUtils/#InteractiveUtils.edit-Tuple{AbstractString,%20Integer}) for more details). Otherwise, open the file with your text editor or IDE of choice, such as Vim, Emacs, VS Code, etc. A convenient way to do that is by entering Julia's shell mode by executing `;` at the REPL:
```julia
julia> ;

shell> vi 1-julia-intro.jl

```
and press delete/backspace to go back to the Julia REPL. When you are done editing the tutorial script, simply save the file and `include` the file again to get a new `main` function to execute in your existing Julia session to see your changes reflected in the output:
```julia
julia> include("1-julia-intro.jl")
main (generic function with 1 method)

julia> res = main();
[...]

```
Note that if you don't call `include` again, you won't see the changes you make to the file reflected when you call the `main` function. (For advanced users, note that you can use [`Revise.includet`](https://timholy.github.io/Revise.jl/stable/cookbook/#includet-usage) as an alternative to `include` which would automatically track changes to the file and update `main` without having to call `include` each time.)

**Note:** We highly recommend keeping your Julia session open throughout each day of the tutorial, which will keep your package environment active and ensure you don't incur re-compilation of precompiled code. If at some point you close your Julia session, make sure to enter the directory corresponding to the tutorial day (i.e. `Tutorials/Day1`) and execute:
```julia
julia> ]

pkg> activate .
```
to activate the environment, which will ensure you have the correct dependencies available to run the tutorial scripts for that day.
