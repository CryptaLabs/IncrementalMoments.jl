using IncrementalMoments
using Base.Test
using StatsBase

@testset "Average only" begin
    values = Float64[5, 6, 5, 7, 8]
    moments = [values[1]]
    for i in 2:length(values)
        @test mean(values[1:i]) ≈ update_moments!(i, values[i], moments)[1]
    end

    values = rand(Float64, 100)
    moments = [values[1]]
    for i in 2:length(values)
        @test mean(values[1:i]) ≈ update_moments!(i, values[i], moments)[1]
    end
end

@testset "Average and variance" begin
    values = Float64[5, 6, 5, 7, 8]
    moments = [values[1], 0]
    for i in 2:length(values)
        update_moments!(i, values[i], moments)
        @test mean(values[1:i]) ≈ moments[1]
        @test var(values[1:i]) ≈ moments[2] / (i - 1)
    end

    values = rand(Float64, 100)
    moments = [values[1], 0]
    for i in 2:length(values)
        update_moments!(i, values[i], moments)
        @test mean(values[1:i]) ≈ moments[1]
        @test var(values[1:i]) ≈ moments[2] / (i - 1)
    end
end

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



