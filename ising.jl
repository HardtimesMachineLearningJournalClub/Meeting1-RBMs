using Base.Test

function pick_trial_sites(trial_x, trial_y, L::Int, periodic_boundaries::Bool)
    if periodic_boundaries
        trial_site_down  = trial_y - 1 > 0  ? trial_y - 1 : L
        trial_site_up    = trial_y + 1 <= L ? trial_y + 1 : 1
        trial_site_right = trial_x + 1 <= L ? trial_x + 1 : 1
        trial_site_left  = trial_x - 1 > 0  ? trial_x - 1 : L
    else
        trial_site_down  = trial_y - 1 > 0  ? trial_y - 1 : nothing
        trial_site_up    = trial_y + 1 <= L ? trial_y + 1 : nothing
        trial_site_right = trial_x + 1 <= L ? trial_x + 1 : nothing
        trial_site_left  = trial_x - 1 > 0  ? trial_x - 1 : nothing
    end
    neighbors = ((trial_x, trial_site_up),
                 (trial_x, trial_site_down),
                 (trial_site_right, trial_y),
                 (trial_site_left, trial_y))
    return neighbors
end

@test pick_trial_sites(1, 1, 4, true) == ((1,2), (1,4), (2,1), (4,1))
@test pick_trial_sites(4, 4, 4, true) == ((4,1), (4,3), (1,4), (3,4))
@test pick_trial_sites(1, 2, 4, true) == ((1,3), (1,1), (2,2), (4,2))

function generate_Ising_configurations(β::Real, L::Int, sample_count::Int=1_000_000, periodic_boundaries::Bool=true)
    starting_state = Array(bitrand(L, L))
    configurations = Vector{Int}[]
    state = copy(starting_state)
    for step in 1:sample_count #do local Metropolis update
        trial_x = rand(1:L)
        trial_y = rand(1:L)
        neighbors = pick_trial_sites(trial_x, trial_y, L, periodic_boundaries)
        diff = 0
        for (neighbor_x, neighbor_y) in neighbors
            if neighbor_x != nothing && neighbor_y != nothing
                diff -= 2*(2*Int(state[trial_x, trial_y] $ state[neighbor_x, neighbor_y])-1)
            end
        end
        if diff < 0 || rand() < exp(-diff*β)
            state[trial_x, trial_y] $= true
        end
        push!(configurations, Vector{Int}(vec(state)))
    end
    return configurations
end

function magnetization(state)
    return abs(sum(2*state - 1))
end

#=for T in linspace(1., 3., 100)
    states = generate_Ising_configurations(1/T, 8)
    magnet = mean(map(x->magnetization(x), states))/64
    @show T, magnet
end=#
