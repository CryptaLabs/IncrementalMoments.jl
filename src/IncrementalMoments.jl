module IncrementalMoments
using Combinatorics: multinomial
export update_moments, update_moments!

update_moments{T}(n::Integer, value::T, moments::AbstractVector{T}) = begin
    update_moments!(n, value, moments::AbstractVector{T}, copy(moments))
end

update_moments!{T}(n::Integer, value::T, moments::AbstractVector{T}) = begin
    update_moments!(n, value, moments::AbstractVector{T}, moments)
end

update_moments!{T <: Number}(n::Integer, value::T, 
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

end # module
