module ElmerCircuitBuilder

using LinearAlgebra

export Resistor, Inductor, Capacitor, VoltageSource, CurrentSource, ElmerComponentMassive, ElmerComponentStranded, ElmerComponentFoil
include("Components.jl")

export Circuit
export create_circuits, add_component!
include("Circuits.jl")

end # module ElmerCircuitBuilder
