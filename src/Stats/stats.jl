# Statistics

@doc raw"""
    dispersion_index(psth; corrected=true, kwargs...) -> Real

Also known as the `Fano factor`.

```math
D = \frac{\sigma^2}{\mu}
```

## optional arguments:
- `dims`: calculate the index over dimension `dims`.

Check out more on [wikipedia](https://en.wikipedia.org/wiki/Index_of_dispersion).
"""
function dispersion_index(psth; corrected=true, kwargs...)
    var(psth; corrected, kwargs...) ./ mean(psth; kwargs...)
end