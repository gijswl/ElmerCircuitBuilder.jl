using Test
using ElmerCircuitBuilder

begin
    c = create_circuits(1)
    
    add_component!(c[1], CurrentSource("I1", (1, 2), 1.0))
    add_component!(c[1], ElmerComponentMassive("T1", (1, 2), 1, [1]))

    write_circuits(c, "circuit1.definition"; path = @__DIR__)
end


begin
    c = create_circuits(1)
    
    add_component!(c[1], CurrentSource("I1", (1, 2), 1.0 * exp(+0im * 2π/3)))
    add_component!(c[1], CurrentSource("I2", (1, 3), 1.0 * exp(+1im * 2π/3)))
    add_component!(c[1], CurrentSource("I3", (1, 4), 1.0 * exp(-1im * 2π/3)))
    add_component!(c[1], ElmerComponentMassive("T1", (1, 2), 1, [4]))
    add_component!(c[1], ElmerComponentMassive("T2", (1, 3), 2, [5]))
    add_component!(c[1], ElmerComponentMassive("T3", (1, 4), 3, [6]))

    write_circuits(c, "circuit2.definition"; path = @__DIR__)
end