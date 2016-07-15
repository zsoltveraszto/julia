# This file is a part of Julia. License is MIT: http://julialang.org/license

function runtests(name, isolate=true)
    if isolate
        mod_name = Symbol("TestMain_", basename(name))
        m = Module(mod_name, true)
        eval(m, :(const $mod_name = $m;
                  eval(x) = Core.eval($m, x);
                  eval(m, x) = Core.eval(m, x)))
        # Make this visible from Main module.
        # Certain functions (e.g. `deserialize`, `subtypes` and `methodswith`
        # relies on this).
        eval(Main, :($mod_name = $m))
        # Non-anonymous version
        # mod_name = Symbol("TestMain_", basename(name))
        # m = eval(Main, :(module $mod_name end))
    else
        m = Main
    end
    eval(m, :(using Base.Test))
    @printf("     \033[1m*\033[0m \033[31m%-21s\033[0m", name)
    tt = @elapsed eval(m, :(include($"$name.jl")))
    if isolate
        eval(Main, :($mod_name = nothing))
    end
    rss = Sys.maxrss()
    @printf(" in %6.2f seconds, maxrss %7.2f MB\n", tt, rss / 2^20)
    rss
end

# looking in . messes things up badly
filter!(x->x!=".", LOAD_PATH)
