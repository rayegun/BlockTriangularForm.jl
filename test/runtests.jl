using BlockTriangularForm
using Test
using SparseArrays
using MatrixMarket
using SuiteSparseMatrixCollection
ssmc = ssmc_db()
include("libbtf.jl")
using .LibBTF
@testset "BlockTriangularForm.jl" begin
    matrix = ssmc[ssmc.name .== "west0479", :]
    path = fetch_ssmc(matrix, format="MM")
    A = mmread(joinpath(path[1], "$(matrix.name[1]).mtx"))
    n = size(A, 1)
    Q, P, R, nmatch, nblocks = BlockTriangularForm.order(A)

    P2 = Vector{Int64}(undef, n)
    Q2 = Vector{Int64}(undef, n)
    R2 = Vector{Int64}(undef, n + 1)
    nmatch2 = Ref{Int64}()
    maxwork = -1.0
    work = Ref{Float64}()
    Work = Vector{Int64}(undef, 5n)
    Ap = SparseArrays.getcolptr(A) .- 1
    Ai = SparseArrays.getrowval(A) .- 1
    nblocks2 = LibBTF.btf_l_order(n, Ap, Ai, maxwork, work, P2, Q2, R2, nmatch2, Work)
    @test P2 .+ 1 == P
    @test Q2 .+ 1 == Q
    println(nblocks)
    println(nblocks2)
    println(R2[1:nblocks2])
    println(R[1:nblocks])
    # R[2] is supposed to be `n` I believe, so incrementing it is unhelpful.
    # TODO: VERIFY THIS WITH TIM
    R2 .+= 1
    @test R2[1:nblocks2] == R[1:nblocks]
    @test nblocks == nblocks2
    @test nmatch == nmatch2[]
end
