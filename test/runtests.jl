
using InteractiveUtils
using Test

using Trixi

@testset "TrixiDebug.jl" begin

versioninfo(verbose=true)

@testset "3D DGMulti failures" begin
  @testset "elixir_euler_weakform.jl (Hexahedral elements)" begin
    # trixi_include(
    #   joinpath(examples_dir(), "dgmulti_3d", "elixir_euler_weakform.jl"),
    #   element_type=Hex()
    # )
    dg = DGMulti(polydeg = 3, element_type = Hex(),
                 surface_integral = SurfaceIntegralWeakForm(FluxHLL()),
                 volume_integral = VolumeIntegralWeakForm())

    equations = CompressibleEulerEquations3D(1.4)
    initial_condition = initial_condition_convergence_test
    source_terms = source_terms_convergence_test

    # example where we tag two separate boundary segments of the mesh
    top_boundary(x, tol=50*eps()) = abs(x[2] - 1) < tol
    rest_of_boundary(x, tol=50*eps()) = !top_boundary(x, tol)
    is_on_boundary = Dict(:top => top_boundary, :rest => rest_of_boundary)

    mesh = DGMultiMesh(dg, cells_per_dimension=(4, 4, 4), is_on_boundary=is_on_boundary)

    boundary_condition_convergence_test = BoundaryConditionDirichlet(initial_condition)
    boundary_conditions = (; :top => boundary_condition_convergence_test,
                            :rest => boundary_condition_convergence_test)

    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, dg,
                                        source_terms = source_terms,
                                        boundary_conditions = boundary_conditions)

    tspan = (0.0, 0.1)
    ode = semidiscretize(semi, tspan)

    u_ode = copy(ode.u0)
    du_ode = zero(u_ode)
    Trixi.rhs!(du_ode, u_ode, semi, first(tspan))
    @info "finished calling rhs!" sum(abs, du_ode)
  end

  @testset "elixir_euler_curved.jl (Hex elements, GaussSBP, flux differencing)" begin
    # trixi_include(
    #   joinpath(examples_dir(), "dgmulti_3d", "elixir_euler_curved.jl"),
    #   approximation_type=GaussSBP()
    # )
    dg = DGMulti(polydeg = 3, element_type = Hex(), approximation_type=GaussSBP(),
                surface_integral = SurfaceIntegralWeakForm(FluxHLL()),
                volume_integral = VolumeIntegralFluxDifferencing(flux_ranocha))

    equations = CompressibleEulerEquations3D(1.4)
    initial_condition = initial_condition_convergence_test
    source_terms = source_terms_convergence_test

    # example where we tag two separate boundary segments of the mesh
    top_boundary(x, tol=50*eps()) = abs(x[2] - 1) < tol
    rest_of_boundary(x, tol=50*eps()) = !top_boundary(x, tol)
    is_on_boundary = Dict(:top => top_boundary, :rest => rest_of_boundary)

    function mapping(xi, eta, zeta)
      x = xi   + 0.1 * sin(pi * xi) * sin(pi * eta)
      y = eta  + 0.1 * sin(pi * xi) * sin(pi * eta)
      z = zeta + 0.1 * sin(pi * xi) * sin(pi * eta)
      return SVector(x, y, z)
    end
    cells_per_dimension = (4, 4, 4)
    mesh = DGMultiMesh(dg, cells_per_dimension, mapping, is_on_boundary=is_on_boundary)

    boundary_condition_convergence_test = BoundaryConditionDirichlet(initial_condition)
    boundary_conditions = (; :top => boundary_condition_convergence_test,
                            :rest => boundary_condition_convergence_test)

    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, dg,
                                        source_terms = source_terms,
                                        boundary_conditions = boundary_conditions)

    tspan = (0.0, 0.1)
    ode = semidiscretize(semi, tspan)

    u_ode = copy(ode.u0)
    du_ode = zero(u_ode)
    Trixi.rhs!(du_ode, u_ode, semi, first(tspan))
    @info "finished calling rhs!" sum(abs, du_ode)
  end

  @testset "elixir_euler_weakform_periodic.jl (Hexahedral elements)" begin
    # trixi_include(
    #   joinpath(examples_dir(), "dgmulti_3d", "elixir_euler_weakform_periodic.jl"),
    #   element_type=Hex()
    # )
    dg = DGMulti(polydeg = 3, element_type = Tet(), approximation_type = Polynomial(),
                surface_integral = SurfaceIntegralWeakForm(FluxHLL()),
                volume_integral = VolumeIntegralWeakForm())

    equations = CompressibleEulerEquations3D(1.4)
    initial_condition = initial_condition_convergence_test
    source_terms = source_terms_convergence_test

    mesh = DGMultiMesh(dg, cells_per_dimension=(4, 4, 4), periodicity=true)

    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, dg,
                                        source_terms = source_terms)

    tspan = (0.0, 0.1)
    ode = semidiscretize(semi, tspan)

    u_ode = copy(ode.u0)
    du_ode = zero(u_ode)
    Trixi.rhs!(du_ode, u_ode, semi, first(tspan))
    @info "finished calling rhs!" sum(abs, du_ode)
  end
end

end # "TrixiDebug.jl"
