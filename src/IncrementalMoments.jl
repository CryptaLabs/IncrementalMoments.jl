@doc readstring(joinpath(dirname(@__DIR__), "README.md")) ->
module IncrementalMoments
using Combinatorics: multinomial
using ArgCheck
import StatsBase: moment

export update_moments, update_moments!, IncrementalMoment
export update!, update

macro lintpragma(s) end
@lintpragma("Ignore use of undeclared variable update!")
@lintpragma("Ignore use of undeclared variable update_moments!")

"""
Object to compute moments of a streaming distribution
"""
abstract type IncrementalMoment end

"""
Incremental moment specialized for scalar inputs.
"""
mutable struct ScalarIncrementalMoment{T <: AbstractFloat} <: IncrementalMoment
    moments::Vector{T}
    nvalues::Int64
end

"""
Incremental moments specialized for array inputs
"""
mutable struct ArrayIncrementalMoment{T <: AbstractArray} <: IncrementalMoment
    moments::T
    nvalues::Int64
end

IncrementalMoment(T::Type{<: AbstractFloat}, order=2) =
    ScalarIncrementalMoment(zeros(typeof(oneunit(T) / 2), order), 0)

IncrementalMoment(T::Type{<: AbstractFloat},
                  dims::Tuple{Integer, Vararg{Integer}},
                  order=2) =
    ArrayIncrementalMoment(zeros(typeof(oneunit(T) / 2), (dims..., order)), 0)

ArrayIncrementalMoment(array::AbstractArray) = begin
    @argcheck ndims(array) > 1 "Last dimension runs over moments"
    ArrayIncrementalMoment(zeros(array), 0)
end

IncrementalMoment(value::AbstractArray{<: Real}, order=2) = begin
    m = similar(value ./ 2, tuple(size(value)..., order))
    fill!(m, 0)
    result = ArrayIncrementalMoment(m, 0)
    update!(result, value)
    result
end

update_moments(n::Integer, value::Real, moments::AbstractVector) = begin
    update_moments!(n, value, moments, deepcopy(moments))
end
update_moments(n::Integer, value::AbstractArray, moments::AbstractArray) = begin
    update_moments!(n, value, moments, deepcopy(moments))
end

update_moments!(n::Integer, value::Real, moments::AbstractVector) = begin
    update_moments!(n, value, moments, moments)
end
update_moments!(n::Integer, value::AbstractArray, moments::AbstractArray) = begin
    update_moments!(n, value, moments, moments)
end

update_moments!(n::Integer, value::Real, moments::AbstractVector{<: AbstractFloat},
                out::AbstractVector{<: AbstractFloat}) = begin
    length(moments) == 0 && return out

    δ = value - moments[1]
    δ_div_n = δ/n

    out[1] = moments[1] + δ_div_n
    length(moments) == 1 && return out

    for p in Iterators.drop(eachindex(moments), 1)
        result = moments[p]
        factor = δ_div_n
        for k in 1:(p - 2)
            result -= multinomial(k, p - k) * factor * moments[p - k]
            factor *= δ_div_n
        end
        out[p] = result + δ^p - δ*factor
    end
    out
end

@inline slice_last(a, n) = view(a, (Base.OneTo(u) for u in Base.front(size(a)))..., n)

update_moments!{T1 <: Real, T2 <:AbstractFloat}(
                n::Integer, value::AbstractArray{T1},
                moments::AbstractArray{T2},
                out::AbstractArray{T2}) = begin

    @argcheck ndims(value) + 1 == ndims(moments)
    @argcheck size(value) == size(moments)[1:end - 1]
    @argcheck size(moments) == size(out)
    length(moments) == 0 && return out

    δ = value .- slice_last(moments, 1)
    δⁿ = δ .* δ
    δ_div_n = δ/n

    slice_last(out, 1) .= slice_last(moments, 1) .+ δ_div_n
    size(moments, ndims(moments)) == 1 && return out

    for p in 2:size(moments, ndims(moments))
        result = slice_last(moments, p)
        factor = copy(δ_div_n)
        for k in 1:(p - 2)
            result .-= multinomial(k, p - k) .* factor .* slice_last(out, p - k)
            factor .*= δ_div_n
        end
        slice_last(out, p) .= result .+ δⁿ .- δ.*factor
        δⁿ .*= δ
    end
    out
end

update!(m::ScalarIncrementalMoment, value) = begin
    if m.nvalues == 0
        m.nvalues = 1
        m.moments[1] = value
    else
        m.nvalues += 1
        update_moments!(m.nvalues, value, m.moments)
    end
    m
end


update(m::IncrementalMoment, value) = update!(deepcopy(m), value)
Base.mean(m::ScalarIncrementalMoment) = m.moments[1]
Base.var(m::ScalarIncrementalMoment; corrected=true) = begin
    m.moments[2] / (corrected ? m.nvalues - 1: m.nvalues)
end
moment(m::ScalarIncrementalMoment, n::Integer) = begin
    n == 1 ? convert(eltype(m.moments), 0): m.moments[n] / m.nvalues
end
Base.length(m::ScalarIncrementalMoment) = length(m.moments)
Base.size(m::ScalarIncrementalMoment) = length(m.moments)
Base.size(m::IncrementalMoment) = size(m.moments)
Base.size(m::IncrementalMoment, i::Integer) = size(m.moments, i)
Base.length(m::IncrementalMoment) = length(m.moments)
Base.ndims(::ScalarIncrementalMoment) = 1
Base.ndims(m::ArrayIncrementalMoment) = ndims(m.moments)
Base.mean(m::ArrayIncrementalMoment) = copy(slice_last(m.moments, 1))
Base.var(m::ArrayIncrementalMoment; corrected=true) = begin
    s = slice_last(m.moments, 2)
    similar(s) .= s ./ (corrected ? m.nvalues - 1: m.nvalues)
end
moment(m::ArrayIncrementalMoment, n::Integer) = begin
    if n == 1
        result = similar(m.moments, size(m.moments)[1:end - 1])
        fill!(result, 0)
        result
    else
        slice_last(m.moments, n) / m.nvalues
    end
end
update!(m::ArrayIncrementalMoment, value) = begin
    if m.nvalues == 0
        m.nvalues = 1
        slice_last(m.moments, 1) .= value
    else
        m.nvalues += 1
        update_moments!(m.nvalues, value, m.moments)
    end
    m
end


end # module
