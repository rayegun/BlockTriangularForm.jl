using BlockTriangularForm
using Test
using SparseArrays
include("libbtf.jl")
using .LibBTF
@testset "BlockTriangularForm.jl" begin
    n = rand(25:400)
    s = sprand(n, n, 0.4)
    Q, P, R, nmatch, nblocks = BlockTriangularForm.btf(s)

    P2 = Vector{Int64}(undef, n)
    Q2 = Vector{Int64}(undef, n)
    R2 = Vector{Int64}(undef, n + 1)
    nmatch2 = Ref{Int64}()
    maxwork = -1.0
    work = Ref{Float64}()
    Work = Vector{Int64}(undef, 5n)
    Ap = SparseArrays.getcolptr(s) .- 1
    Ai = SparseArrays.getrowval(s) .- 1
    nblocks2 = LibBTF.btf_l_order(n, Ap, Ai, maxwork, work, P2, Q2, R2, nmatch2, Work)
    @test P2 .+ 1 == P
    @test Q2 .+ 1 == Q
    # R[2] is supposed to be `n` I believe, so incrementing it is unhelpful.
    # TODO: VERIFY THIS WITH TIM
    R2 .+= 1; R2[2] -= 1
    @test R2 == R
    @test nblocks == nblocks2
    @test nmatch == nmatch2[]
end
