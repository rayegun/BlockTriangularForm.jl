using SparseArrays
using BenchmarkTools
include("libbtf.jl")
n = 2000
s = sprand(n, n, 0.4)
Q, P, R, nmatch, nblocks = @btime BlockTriangularForm.btf(s)

P2 = Vector{Int64}(undef, n)
Q2 = Vector{Int64}(undef, n)
R2 = Vector{Int64}(undef, n + 1)
nmatch2 = Ref{Int64}()
maxwork = -1.0
work = Ref{Float64}()
Work = Vector{Int64}(undef, 5n)
Ap = SparseArrays.getcolptr(s) .- 1
Ai = SparseArrays.getrowval(s) .- 1
nblocks2 = @btime LibBTF.btf_l_order(n, Ap, Ai, maxwork, work, P2, Q2, R2, nmatch2, Work)
