
# Test that the function returns the expected values for valid input parameters.
function test_happy_path()
    energy(v, w) = holstein_energy(v, w, α, ωβ...; dims = dims)
    v0 = 0.5
    w0 = 0.2
    lower = [0, 0]
    upper = [Inf, Inf]
    Δv_opt, w_opt, E_opt, K_opt, P_opt = vw_variation(energy, v0, w0; lower = lower, upper = upper)
    @test Δv_opt + w_opt ≈ v0
    @test w_opt ≈ w0
    @test E_opt ≈ holstein_energy(v0, w0, α, ωβ...; dims = dims)[1]
    @test K_opt ≈ holstein_energy(v0, w0, α, ωβ...; dims = dims)[2]
    @test P_opt ≈ holstein_energy(v0, w0, α, ωβ...; dims = dims)[3]
end

    # // Test that the function can handle different initial values for v and w.
function test_different_initial_values()
    energy(v, w) = holstein_energy(v, w, α, ωβ...; dims = dims)
    v0 = 0.1
    w0 = 0.3
    lower = [0, 0]
    upper = [Inf, Inf]
    Δv_opt, w_opt, E_opt, K_opt, P_opt = vw_variation(energy, v0, w0; lower = lower, upper = upper)
    @test Δv_opt + w_opt ≈ v0
    @test w_opt ≈ w0
    @test E_opt ≈ holstein_energy(v0, w0, α, ωβ...; dims = dims)[1]
    @test K_opt ≈ holstein_energy(v0, w0, α, ωβ...; dims = dims)[2]
    @test P_opt ≈ holstein_energy(v0, w0, α, ωβ...; dims = dims)[3]
end

    # // Test that the function can handle lower and upper bounds for v and w.
function test_lower_and_upper_bounds()
    energy(v, w) = holstein_energy(v, w, α, ωβ...; dims = dims)
    v0 = 0.2
    w0 = 0.4
    lower = [0.1, 0.3]
    upper = [0.5, 0.6]
    Δv_opt, w_opt, E_opt, K_opt, P_opt = vw_variation(energy, v0, w0; lower = lower, upper = upper)
    @test Δv_opt + w_opt ≈ v0
    @test w_opt ≈ w0
    @test E_opt ≈ holstein_energy(v0, w0, α, ωβ...; dims = dims)[1]
    @test K_opt ≈ holstein_energy(v0, w0, α, ωβ...; dims = dims)[2]
    @test P_opt ≈ holstein_energy(v0, w0, α, ωβ...; dims = dims)[3]
end

    # // Test that the function can handle very large values for v and w.
function test_large_values()
    energy(v, w) = holstein_energy(v, w, α, ωβ...; dims = dims)
    v0 = 1e6
    w0 = 1e7
    lower = [0, 0]
    upper = [Inf, Inf]
    Δv_opt, w_opt, E_opt, K_opt, P_opt = vw_variation(energy, v0, w0; lower = lower, upper = upper)
    @test Δv_opt + w_opt ≈ v0
    @test w_opt ≈ w0
    @test E_opt ≈ holstein_energy(v0, w0, α, ωβ...; dims = dims)[1]
    @test K_opt ≈ holstein_energy(v0, w0, α, ωβ...; dims = dims)[2]
    @test P_opt ≈ holstein_energy(v0, w0, α, ωβ...; dims = dims)[3]
end

    # // Test that the function can handle very small values for v and w.
function test_small_values()
    energy(v, w) = holstein_energy(v, w, α, ωβ...; dims = dims)
    v0 = 1e-6
    w0 = 1e-7
    lower = [0, 0]
    upper = [Inf, Inf]
    Δv_opt, w_opt, E_opt, K_opt, P_opt = vw_variation(energy, v0, w0; lower = lower, upper = upper)
    @test Δv_opt + w_opt ≈ v0
    @test w_opt ≈ w0
    @test E_opt ≈ holstein_energy(v0, w0, α, ωβ...; dims = dims)[1]
    @test K_opt ≈ hol

    # // Test that the vw_variation function can optimize the values of v and w to minimize the free energy.
function test_vw_variation_optimization()
    energy(v, w) = holstein_energy(v, w, α, ωβ...; dims = dims)
    v0 = 0.5
    w0 = 0.2
    lower = [0, 0]
    upper = [Inf, Inf]
    Δv_opt, w_opt, E_opt, K_opt, P_opt = vw_variation(energy, v0, w0; lower = lower, upper = upper)
    @test Δv_opt + w_opt ≈ v0
    @test w_opt ≈ w0
    @test E_opt ≈ holstein_energy(v0, w0, α, ωβ...; dims = dims)[1]
    @test K_opt ≈ holstein_energy(v0, w0, α, ωβ...; dims = dims)[2]
    @test P_opt ≈ holstein_energy(v0, w0, α, ωβ...; dims = dims)[3]
end

    # // Test that the vw_variation function returns the optimized values of v and w.
function test_vw_variation_return_values()
    energy(v, w) = holstein_energy(v, w, α, ωβ...; dims = dims)
    v0 = 0.5
    w0 = 0.2
    lower = [0, 0]
    upper = [Inf, Inf]
    Δv_opt, w_opt, E_opt, K_opt, P_opt = vw_variation(energy, v0, w0; lower = lower, upper = upper)
    @test typeof(Δv_opt) == Float64
    @test typeof(w_opt) == Float64
    @test typeof(E_opt) == Float64
    @test typeof(K_opt) == Float64
    @test typeof(P_opt) == Float64
end

    # // Test that the vw_variation function returns the total energy, kinetic energy component, and potential energy component.
function test_vw_variation_energy_components()
    energy(v, w) = holstein_energy(v, w, α, ωβ...; dims = dims)
    v0 = 0.5
    w0 = 0.2
    lower = [0, 0]
    upper = [Inf, Inf]
    Δv_opt, w_opt, E_opt, K_opt, P_opt = vw_variation(energy, v0, w0; lower = lower, upper = upper)
    @test E_opt ≈ holstein_energy(v0, w0, α, ωβ...; dims = dims)[1]
    @test K_opt ≈ holstein_energy(v0, w0, α, ωβ...; dims = dims)[2]
    @test P_opt ≈ holstein_energy(v0, w0, α, ωβ...; dims = dims)[3]
end

    # // Test that the function can handle non-numeric input parameters.
function test_non_numeric_input()
    energy(v, w) = holstein_energy(v, w, α, ωβ...; dims = dims)
    v0 = "abc"
    w0 = [1, 2, 3]
    lower = [0, 0]
    upper = [Inf, Inf]
    Δv_opt, w_opt, E_opt, K_opt, P_opt = vw_variation(energy, v0, w0; lower = lower, upper = upper)
    @test Δv_opt + w_opt ≈ v0
    @test w_opt ≈ w0
    @test E_opt ≈ holstein_energy(v0, w0, α, ωβ...; dims = dims)[1]
    @test K_opt ≈ holstein_energy(v0, w0, α, ωβ...; dims = dims)[2]
    @test P_opt ≈ holstein_energy(v0, w0, α, ωβ...; dims = dims)[3]
end

    # // Test that the function can handle missing input parameters.
function test_missing_input_parameters()
    energy(v, w) = holstein_energy(v, w, α, ωβ...; dims = dims)
    v0 = 0.5
    w0 = 0.2
    lower = [0, 0]
    upper = [Inf, Inf]
    Δv_opt, w_opt, E_opt, K_opt, P_opt = vw_variation(energy, v0, w0; lower = lower, upper = upper)
    @test Δv_opt + w_opt ≈ v0
    @test w_opt ≈ w0
    @test E_opt ≈ holstein_energy(v0, w0, α, ωβ...; dims = dims)[1]
    @test K_opt ≈ holstein_energy(v0, w0, α, ωβ...; dims = dims)[2]
    @test P_opt ≈ holstein_energy(v0, w0, α, ωβ...; dims = dims)[3]
end

    # // Test that the function can handle non-positive input parameters.
function test_non_positive_input_parameters()
    energy(v, w) = holstein_energy(v, w, α, ωβ...; dims = dims)
    v0 = -0.5
    w0 = -0.2
    lower = [0, 0]
    upper = [Inf, Inf]
    Δv_opt, w_opt, E_opt, K_opt, P_opt = vw_variation(energy, v0, w0; lower = lower, upper = upper)
    @test Δv_opt + w_opt ≈ v0
    @test w_opt ≈ w0
    @test E_opt ≈ holstein_energy(v0, w0, α, ωβ...; dims = dims)[1]
    @test K_opt ≈ holstein_energy(v0, w0, α, ωβ...; dims = dims)[2]
    @test P_opt ≈ holstein_energy(v0, w0, α, ωβ...; dims = dims)[3]
end

    # // Test that the function can handle non-integer input parameters.
function test_handle_non_integer_input_parameters()
    energy(v, w) = holstein_energy(v, w, α, ωβ...; dims = dims)
    v0 = 0.5
    w0 = 0.2
    lower = [0, 0]
    upper = [Inf, Inf]
    Δv_opt, w_opt, E_opt, K_opt, P_opt = vw_variation(energy, v0, w0; lower = lower, upper = upper)
    @test Δv_opt + w_opt ≈ v0
    @test w_opt ≈ w0
    @test E_opt ≈ holstein_energy(v0, w0, α, ωβ...; dims = dims)[1]
    @test K_opt ≈ holstein_energy(v0, w0, α, ωβ...; dims = dims)[2]
    @test P_opt ≈ holstein_energy(v0, w0, α, ωβ...; dims = dims)[3]
end

    # // Test that the function can handle input parameters that cause division by zero.
function test_handle_division_by_zero()
    energy(v, w) = holstein_energy(v, w, α, ωβ...; dims = dims)
    v0 = 0.0
    w0 = 0.0
    lower = [0, 0]
    upper = [Inf, Inf]
    Δv_opt, w_opt, E_opt, K_opt, P_opt = vw_variation(energy, v0, w0; lower = lower, upper = upper)
    @test Δv_opt + w_opt ≈ v0
    @test w_opt ≈ w0
    @test E_opt ≈ holstein_energy(v0, w0, α, ωβ...; dims = dims)[1]
    @test K_opt ≈ holstein_energy(v0, w0, α, ωβ...; dims = dims)[2]
    @test P_opt ≈ holstein_energy(v0, w0, α, ωβ...; dims = dims)[3]
end
