include("basic_omnivory_module.jl")
include("top_heavy.jl")

function thplot(fc = "foodchain", force = "unforced")
    u0 = [1.0, 0.5, 0.1]
    t_span = (0.0, 1000.0)
    par = OmnPar()
    kvals = 1.5:0.1:3
    thvals = fill(0.0, length(kvals))

    if fc == "foodchain"
        par.ω = 0.0
        par.a_RP = 0.0
    end
    if force == "unforced"
        par.A = 0.0
    end

    for (ki, kval) in enumerate(kvals)
        par.K = kval
        prob = ODEProblem(model!, u0, t_span, par)
        sol = solve(prob)
        thvals[ki] = top_heavy(sol)[3]
    end

    plot(kvals, thvals)
    xlabel("K")
    ylabel("Predator/Consumer")
    ylim(0,12)
end

let
    th = figure()
    subplots_adjust(hspace = 0.3, wspace = 0.5, top = 0.9, bottom = 0.1, right = 0.75) #wspace = 0.4
    subplot(2,2,1)
    title("Unforced", fontsize = 15)
    thplot("foodchain", "unforced")
    subplot(2,2,2)
    title("Forced", fontsize = 15)
    thplot("foodchain", "forced")
    subplot(2,2,3)
    thplot("omn", "unforced")
    subplot(2,2,4)
    thplot("omn", "forced")
    annotate("Food Chain", (350, 220), xycoords = "figure points", fontsize = 15, rotation = 90)
    annotate("Omnivory", (350, 75), xycoords = "figure points", fontsize = 15, rotation = 90)
    return th
end


#NOTE when K < 1.5 returnng with error - Warning: dt <= dtmin. Aborting. There is either an error in your model specification or the true solution is unstable.
