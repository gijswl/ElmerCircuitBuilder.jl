module ElmerCircuitBuilder

import Dates: format, now

using Elmer
using LinearAlgebra: Diagonal
using OrderedCollections: OrderedDict

export Resistor, Inductor, Capacitor, VoltageSource, CurrentSource, ElmerComponent
include("Components.jl")

export Circuit
export create_circuits, insert_component!
include("Circuits.jl")

export write_circuits
include("Write.jl")

function add_circuits!(sif::SolverInformationFile, circuits::Vector{Circuit})
    components = Dict()
    body_forces = Dict()

    for c ∈ circuits
        bf_data = OrderedDict()

        for component ∈ get_components(c)
            if (typeof(component) <: ElmerComponent)
                id = add_component!(sif, component.name; master_bodies=component.master_bodies)
                @assert id == component.component_id "Circuit component numbering should be monotonically increasing."
                components[id] = component
            elseif (typeof(component) <: SourceComponent)
                entry = format_source(component.name)

                if (typeof(component.value) <: Complex)
                    bf_data["$entry re"] = "Real $(real(component.value))"
                    bf_data["$entry im"] = "Real $(imag(component.value))"
                elseif (typeof(component.value) <: Real)
                    bf_data[entry] = "Real $(component.value)"
                end
            end
        end

        body_forces[c.index] = add_body_force!(sif, "Circuit $(c.index)"; data=bf_data)
    end

    return components, body_forces
end

end # module ElmerCircuitBuilder
