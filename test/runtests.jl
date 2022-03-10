
using InteractiveUtils
using Test

using Trixi

@testset "TrixiDebug.jl" begin

versioninfo(verbose=true)

@testset "3D DGMulti failures" begin
  @testset "elixir_euler_weakform.jl (Hexahedral elements)" begin
    trixi_include(
      joinpath(get_examples(), "dgmulti_3d", "elixir_euler_weakform.jl"),
      element_type=Hex()
    )
  end

  @testset "elixir_euler_curved.jl (Hex elements, GaussSBP, flux differencing)" begin
    trixi_include(
      joinpath(get_examples(), "dgmulti_3d", "elixir_euler_curved.jl"),
      approximation_type=GaussSBP()
    )
  end

  @testset "elixir_euler_weakform_periodic.jl (Hexahedral elements)" begin
    trixi_include(
      joinpath(get_examples(), "dgmulti_3d", "elixir_euler_weakform_periodic.jl"),
      element_type=Hex()
    )
  end
end

end # "TrixiDebug.jl"
