# Labmnivory
[![Reproduce the analysis](https://github.com/McCannLab/Labmnivory/actions/workflows/reproduce.yaml/badge.svg)](https://github.com/McCannLab/Labmnivory/actions/workflows/reproduce.yaml)


## Installation

First install a recent version of [Julia](https://julialang.org/) (note that the 
code is only tested for Julia v1.6.0) and the following packages (which can be installed using `src/install_packages.jl` in this repository, see below):

|Package              | Links                                                  |
|:--------------------|:-------------------------------------------------------|
|DifferentialEquations| [repo](https://github.com/SciML/DifferentialEquations.jl)|
|ForwardDiff          | [repo](https://github.com/JuliaDiff/ForwardDiff.jl)    |
|LinearAlgebra        | [doc](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/)|
|NLsolve              | [repo](https://github.com/JuliaNLSolvers/NLsolve.jl)   |
|Parameters           | [repo](https://github.com/mauro3/Parameters.jl)        |
|PyPlot               | [repo](https://github.com/JuliaPy/PyPlot.jl)           |
|QuadGK               | [repo](https://github.com/JuliaMath/QuadGK.jl)         |
|RecursiveArrayTool   | [repo](https://github.com/SciML/RecursiveArrayTools.jl)|
|Statistics           | [repo](https://docs.julialang.org/en/v1/stdlib/Statistics/)|


Then clone this repository using [git](https://git-scm.com/)

```
git clone https://github.com/McCannLab/Labmnivory.git
```

or downloaded the [zip file](https://github.com/McCannLab/Labmnivory/archive/refs/heads/master.zip). If you haven't yet installed the package you can actually 
use this repository. In a terminal, set your working directory at the root of the freshly cloned/downloaded repository, then use the following command line&nbsp;:

```julia
julia src/install_packages.jl
```

You should now be ready to reproduce the analysis.


## File structure 

- `src`: source code 
- `figs`: folder where figures are exported



## How to reproduce the analysis 

In a terminal, assuming the working directory is set to be the the root of the freshly cloned/downloaded repository, use the following command line&nbsp;:

```julia
julia src/pipeline.jl
```




## Style guide

This section describes the coding style rules I (Gabriel) prefer. They are heavily based
on the `JuMP` style guides, though there are exceptions.

### Formatting

Julia unfortunately does not have an autoformatting tool like
[gofmt](https://blog.golang.org/go-fmt-your-code). Until a reliable
autoformatting tool is available, we adopt the following conventions.

#### Whitespace

For conciseness, never use more than one blank line within a function, and never
begin a function with a blank line.

Bad:
```julia
function foo(x)
    y = 2 * x


    return y
end

function foo(x)

    y = 2 * x
    return y
end
```

Julia is mostly insensitive to whitespace characters within lines.
For consistency:

- Use spaces between binary operators (with some exceptions, see below)
- Use a single space after commas and semicolons
- Do not use extra spaces for unary operators, parentheses, or braces
- Indent within new blocks (except `module`) using 4 spaces

Good:
```julia
f(x, y) = [3 * dot(x, y); x']
```

Bad:
```julia
f(x,y) = [ 3*dot(x,y) ; x' ]
```

Good:
```julia
module Foo

function f(x)
    return x + 1
end

end # module Foo
```

##### Exceptions

We make an exception for the `:` operator when it is used to form a range.

Good:
```julia
x = 1:5
```

Bad:
```julia
x = 1 : 5
```

One reason is that it can be confused with Julia's conditional statement:
`cond ? x : y` which requires whitespace around the `:`.

Good:
```julia
2 * x  # This is preferred if space is not an issue.
2 * (x + 1)
```

Bad:
```julia
2(x + 1)
```

#### Return statements

To avoid situations in which it is unclear whether the author intended to return
a certain value or not, always use an explicit `return` statement to exit from a
function. If the return from a function is `nothing`, use `return` instead of
`return nothing`.

We make an exception for assignment-form one-line functions (`f(x) = 2 * x`).

Good:
```julia
foo(x) = 2 * x  # Acceptable if one line
function foo(x)
    return 2 * x
end
function foo!(x)
    x[1] += 1
    return
end
```

Bad:
```julia
function foo(x)
    2x
end
function foo!(x)
    x[1] += 1
    return nothing
end
```

#### Figure plotting

We will use the package [PyPlot](https://github.com/JuliaPy/PyPlot.jl) for plotting. Go [here](https://gist.github.com/gizmaa/7214002) for example coding with PyPlot.

To automatically open plots outside of Atom, disable the Plot Pane in the UI Options of the Juno Settings.

Whenever plotting a figure ensure that the code follows below:

```julia
    plot_name = figure()
    plot()
    plotoptions
    return plot_name
```

To ease figure creation, place plotting code in a let block.
 
