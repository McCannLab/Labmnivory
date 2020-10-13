# # Feral Pig Masting Events
# I want to try to model the situation where we have a relatively low
# productivity environment, but we have period high productivity pulses,
# like a oak mastings etc.
#
# Objectives are that we can see a response of omnivory in the face of this
# variable productivity situation, and test that we are measuring omnivory
# in a way that is enlightening.
first_mast = 1000.0
mast_freq = 1000.0
mast_length = 1.0
mast_strength = 2.0
t_end = first_mast + mast_freq
t_span = (0.0, t_end)
t_grid = range(0.0, t_end, length = 1000)
t_start = 75.0

dt = .01
# time after perturbation
t_aft = (first_mast + mast_length):dt:t_end
# perturbation + after perturbation
t_all = first_mast:dt:t_end

# only one event
mast_starts = [1000]
mast_ends = mast_starts .+ mast_length
mast_event_times = sort(union(mast_starts, mast_ends))

masting_event(u, t, integrator) = t ∈ mast_event_times
function forcing!(integrator)
    if integrator.t ∈ mast_starts
        integrator.p.K = mast_strength + integrator.p.K_base
    elseif integrator.t ∈ mast_ends
        integrator.p.K = integrator.p.K_base
    end
    return
end
cb = DiscreteCallback(masting_event, forcing!)
