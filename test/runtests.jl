using Test

using Elmer
using ElmerCircuitBuilder

data_path = joinpath(@__DIR__, "simdata\\")

begin
    c = create_circuits(1)

    insert_component!(c[1], CurrentSource("I1", (1, 2), 1.0))
    insert_component!(c[1], ElmerComponent("T1", (1, 2), 1, [1]))

    write_circuits(c, "circuit1.definition"; path=data_path)

    sim = Simulation(7, Elmer.CoordinateCartesian(), Elmer.SimulationScanning())
    sif = Elmer.SolverInformationFile("case1", sim, data_path=data_path)

    add_body!(sif, "Coil")

    components, body_forces = add_circuits!(sif, c)

    add_coil_data!(sif, 1, CoilStranded(20, 0), TerminalClosed(100e-6))

    write(sif)
end


begin
    c = create_circuits(1)

    insert_component!(c[1], CurrentSource("I1", (1, 2), 1.0 * exp(+0im * 2π / 3)))
    insert_component!(c[1], CurrentSource("I2", (1, 3), 1.0 * exp(+1im * 2π / 3)))
    insert_component!(c[1], CurrentSource("I3", (1, 4), 1.0 * exp(-1im * 2π / 3)))
    insert_component!(c[1], ElmerComponent("T1", (1, 2), 1, [1]))
    insert_component!(c[1], ElmerComponent("T2", (1, 3), 2, [2]))
    insert_component!(c[1], ElmerComponent("T3", (1, 4), 3, [3]))

    write_circuits(c, "circuit2.definition"; path=data_path)

    sim = Simulation(7, Elmer.CoordinateCartesian(), Elmer.SimulationScanning())
    sif = Elmer.SolverInformationFile("case2", sim, data_path=data_path)

    add_body!(sif, "Coil1")
    add_body!(sif, "Coil2")
    add_body!(sif, "Coil3")

    components, body_forces = add_circuits!(sif, c)

    add_coil_data!(sif, 1, CoilMassive())
    add_coil_data!(sif, 2, CoilMassive())
    add_coil_data!(sif, 3, CoilMassive())

    write(sif)
end