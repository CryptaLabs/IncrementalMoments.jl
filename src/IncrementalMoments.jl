module IncrementalMoments
using Combinatorics: multinomial
import StatsBase: moment
export update_moments, update_moments!, IncrementalMoment
export update!, update


"""
Object to compute moments of a streaming distribution
"""
abstract IncrementalMoment
type ScalarIncrementalMoment{T <: AbstractFloat} <: IncrementalMoment
    moments::Vector{T}
    nvalues::Int64
end

(::Type{IncrementalMoment}){T <: AbstractFloat}(::Type{T}, order=2) = begin
    ScalarIncrementalMoment(zeros(typeof(one(T) / 2), order), 0)
end

update_moments(n::Integer, value::Real, moments::AbstractVector) = begin
    update_moments!(n, value, moments, deepcopy(moments))
end

update_moments!{T <: AbstractFloat}(n::Integer, value::Real,
                                    moments::AbstractVector{T}) = begin
    update_moments!(n, value, moments, moments)
end

update_moments!{T <: AbstractFloat}(n::Integer, value::Real,
                                    moments::AbstractVector{T},
                                    out::AbstractVector{T}) = begin
    length(moments) == 0 && return out

    δ = value .- moments[1]
    δ_div_n = δ/n

    out[1] = moments[1] + δ_div_n
    length(moments) == 1 && return out

    for p in drop(eachindex(moments), 1)
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

update!(m::IncrementalMoment, value) = begin
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
Base.mean(m::IncrementalMoment) = m.moments[1]
Base.var(m::IncrementalMoment; corrected=true) = begin
    m.moments[2] / (corrected ? m.nvalues - 1: m.nvalues)
end
moment(m::IncrementalMoment, n::Integer) = begin
    n == 1 ? convert(eltype(m.moments), 0): m.moments[n] / m.nvalues
end
Base.size(m::IncrementalMoment) = length(m.moments)

end # module
