"""
    abstract type AbstractComponent

Represents an electrical circuit component.
"""
abstract type AbstractComponent end

"""
    mutable struct Resistor

Represents a resistor named `name` with a `value` in [Ω]. It is connected between two circuit `nodes`.

It is used to build the Elmer circuit's B-matrix.
"""
mutable struct Resistor <: AbstractComponent
    name::String
    nodes::NTuple{2,Int}
    value::Real
end

"""
    mutable struct Inductor

Represents a inductor named `name` with a `value` in [H]. It is connected between two circuit `nodes`.

It is used to build the Elmer circuit's A-matrix.
"""
mutable struct Inductor <: AbstractComponent
    name::String
    nodes::NTuple{2,Int}
    value::Real
end

"""
    mutable struct Capacitor

Represents a capacitor named `name` with a `value` in [F]. It is connected between two circuit `nodes`.

It is used to build the Elmer circuit's A-matrix.
"""
mutable struct Capacitor <: AbstractComponent
    name::String
    nodes::NTuple{2,Int}
    value::Real
end

"""
    mutable struct VoltageSource

Represents a voltage source named `name` with a `value` in [V]. It is connected between two circuit `nodes`.

It is used to build the Elmer circuit's B-matrix.
"""
mutable struct VoltageSource <: AbstractComponent
    name::String
    nodes::NTuple{2,Int}
    value::Real
end

"""
    mutable struct CurrentSource

Represents a current source named `name` with a `value` in [A]. It is connected between two circuit `nodes`.
Current is defined as flowing _out_ of the positive terminal `nodes[1]`.

It is used to build the Elmer circuit's B-matrix.
"""
mutable struct CurrentSource <: AbstractComponent
    name::String
    nodes::NTuple{2,Int}
    value::Real
end

"""
    mutable struct ElmerComponentMassive

Represents a massive 2D or 3D coil named `name` connected between two circuit `nodes`.
The component is associated with the Component (`component_id`) section in the SIF and the physical bodies in `master_bodies`.

Sector represents the integer associated with Symmetry Coefficient under Elmer's component. 
By default the value is 1, describing the full dimention of the circuit. 
Change the value accordingly depending on the symmetry of the problem at hand. 
For example if you're modeling half of the coil, then the value of sector is 0.5.

If the 3D coil is not `closed`, two terminal `boundaries` must be defined.

It is used to build the Elmer circuit's B-matrix.
"""
mutable struct ElmerComponentMassive <: AbstractComponent
    name::String
    nodes::NTuple{2,Int}
    component_id::Int
    master_bodies::Vector{Int}
    sector::Real
    closed::Bool
    boundaries::NTuple{2,Int}

    function ElmerComponentMassive(name, nodes, component_id, master_bodies; sector=1, closed=true, boundaries=(-1, -1))
        new(name, nodes, component_id, master_bodies, sector, closed, boundaries)
    end
end

"""
    mutable struct ElmerComponentStranded

Represents a stranded 2D or 3D coil named `name` connected between two circuit `nodes`.
The stranded coil has a `number_turns` turns and a DC resistance of `resistance`.
The component is associated with the Component (`component_id`) section in the SIF and the physical bodies in `master_bodies`.

Sector represents the integer associated with Symmetry Coefficient under Elmer's component. 
By default the value is 1, describing the full dimention of the circuit. 
Change the value accordingly depending on the symmetry of the problem at hand. 
For example if you're modeling half of the coil, then the value of sector is 0.5.

If the 3D coil is not `closed`, two terminal `boundaries` must be defined.

It is used to build the Elmer circuit's B-matrix.
"""
mutable struct ElmerComponentStranded <: AbstractComponent
    name::String
    nodes::NTuple{2,Int}
    component_id::Int
    master_bodies::Vector{Int}
    sector::Real
    closed::Bool
    boundaries::NTuple{2,Int}

    number_turns::Int
    resistance::Real

    function ElmerComponentStranded(name, nodes, component_id, master_bodies, number_turns, resistance; sector=1, closed=true, boundaries=(-1, -1))
        new(name, nodes, component_id, master_bodies, sector, closed, boundaries, number_turns, resistance)
    end
end

mutable struct ElmerComponentFoil <: AbstractComponent

end

const SourceComponent = Union{VoltageSource,CurrentSource}
const ElmerComponent = Union{ElmerComponentMassive,ElmerComponentStranded,ElmerComponentFoil}