using IncrementalMoments
using Base.Test
using StatsBase

@testset "Average and variance and skewness" begin
    values = Float64[5, 6, 5, 7, 8]
    moments = [values[1], 0, 0]
    for i in 2:length(values)
        update_moments!(i, values[i], moments)
        @test mean(values[1:i]) ≈ moments[1]
        @test moment(values[1:i], 2) ≈ moments[2] / i
        @test moment(values[1:i], 3) ≈ moments[3] / i
    end

    values = rand(Float64, 100)
    moments = [values[1], 0, 0]
    update_moments!(2, values[2], moments)
    for i in 3:length(values)
        update_moments!(i, values[i], moments)
        @test mean(values[1:i]) ≈ moments[1]
        @test moment(values[1:i], 2) ≈ moments[2] / i
        # @test moment(values[1:i], 3) ≈ moments[3] / i
    end
end


@testset "n=10 moments" begin
    N = 10
    values = rand(Float64, 100)
    moments = zeros(Float64, N)
    moments[1] = values[1]
    for i in 2:(N - 1)
        update_moments!(i, values[i], moments)
    end
    for i in N:N
        update_moments!(i, values[i], moments)
        @test mean(values[1:i]) ≈ moments[1]
        for n in 2:N
            @test moment(values[1:i], n) ≈ moments[n] / i
        end
    end
end


@testset "n=10 IncrementalMoment" begin
    N = 10
    values = rand(Float64, 100)
    moments = IncrementalMoment(Float64, N)
    for i in 1:(N - 1)
        update!(moments, values[i])
    end
    for i in N:length(values)
        update!(moments, values[i])
        @test mean(values[1:i]) ≈ mean(moments)
        for n in 2:N
            @test moment(values[1:i], n) ≈ moment(moments, n)
        end
    end
end

@testset "Check slicing by last dim" begin
    values = rand(Int64, (10, 10, 6))
    @test IncrementalMoments.slice_last(values, 1) == values[:, :, 1]
    @test IncrementalMoments.slice_last(values, 5) == values[:, :, 5]

    values = rand(Int64, (10, 6))
    @test IncrementalMoments.slice_last(values, 1) == values[:, 1]
    @test IncrementalMoments.slice_last(values, 5) == values[:, 5]
end

@testset "Incremental moments over arrays" begin
    values = rand(Float64, (10, 10, 100))
    moments = IncrementalMoment(Float64, size(values)[1:end - 1], 4)
    N = size(moments, ndims(moments))
    for i in 1:(N - 1)
        update!(moments, @view values[:, :, i])
    end
    for i in N:size(values, ndims(values))
        update!(moments, @view values[:, :, i])
        @test mean(values[:, :, 1:i], ndims(values)) ≈ mean(moments)
        for n in 2:N
            expected = [moment(values[u, v, 1:i], n)
                        for u = 1:size(values, 1), v = 1:size(values, 2)]
            @test expected ≈ moment(moments, n)
        end
    end
end

