module BlockTriangularForm

flip(j) = -j - 2
isflipped(j) = j < -1
unflip(j) = isflipped(j) ? flip(j) : j
const EMPTY = -1
const UNVISITED = -2
const UNASSIGNED = -1

# ==========================================================================
# === augment ==============================================================
# ==========================================================================
# Perform a depth-first-search starting at column k, to find an augmenting
# path.  An augmenting path is a sequence of row/column pairs (i1,k), (i2,j1),
# (i3,j2), ..., (i(s+1), js), such that all of the following properties hold:
#
#      * column k is not matched to any row
#      * entries in the path are nonzero
#      * the pairs (i1,j1), (i2,j2), (i3,j3) ..., (is,js) have been 
#          previously matched to each other
#      * (i(s+1), js) is nonzero, and row i(s+1) is not matched to any column
#
# Once this path is found, the matching can be changed to the set of pairs
# path.  An augmenting path is a sequence of row/column pairs
#
#      (i1,k), (i2,j1), (i3,j2), ..., (i(s+1), js)
#
# Once a row is matched with a column it remains matched with some column, but
# not necessarily the column it was first matched with.
#
# In the worst case, this function can examine every nonzero in A.  Since it
# is called n times by maxtrans, the total time of maxtrans can be as high as
# O(n*nnz(A)).  To limit this work, pass a value of maxwork > 0.  Then at
# most O((maxwork+1)*nnz(A)) work will be performed; the maximum matching might
# not be found, however.
#
# This routine is very similar to the dfs routine in klu_kernel.c, in the
# KLU sparse LU factorization package.  It is essentially identical to the
# cs_augment routine in CSparse, and its recursive version (augment function
# in cs_maxtransr_mex.c), except that this routine allows for the search to be
# terminated early if too much work is being performed.
#
# The algorithm is based on the paper "On Algorithms for obtaining a maximum
# transversal" by Iain Duff, ACM Trans. Mathematical Software, vol 7, no. 1,
# pp. 315-330, and "Algorithm 575: Permutations for a zero-free diagonal",
# same issue, pp. 387-390.  The code here is a new implementation of that
# algorithm, with different data structures and control flow.  After writing
# this code, I carefully compared my algorithm with MC21A/B (ACM Algorithm 575)
# Some of the comparisons are partial because I didn't dig deeply into all of
# the details of MC21A/B, such as how the stack is maintained.  The following
# arguments are essentially identical between this code and MC21A:
#
# maxtrans     MC21A,B
# --------     -------
# n            N           identical
# k            JORD        identical
# Ap           IP          column / row pointers
# Ai           ICN         row / column indices
# Ap[n]        LICN        length of index array (# of nonzeros in A)
# Match        IPERM       output column / row permutation
# nmatch       NUMNZ       # of nonzeros on diagonal of permuted matrix
# Flag         CV          mark a node as visited by the depth-first-search
#
# The following are different, but analogous:
#
# Cheap        ARP         indicates what part of the a column / row has
#                          already been matched.
#
# The following arguments are very different:
#
# -            LENR        # of entries in each row/column (unused in maxtrans)
# Pstack       OUT         Pstack keeps track of where we are in the depth-
#                          first-search scan of column j.  I think that OUT
#                          plays a similar role in MC21B, but I'm unsure.
# Istack       PR          keeps track of the rows in the path.  PR is a link
#                          list, though, whereas Istack is a stack.  Maxtrans
#                          does not use any link lists.
# Jstack       OUT? PR?    the stack for nodes in the path (unsure)
#
# The following control structures are roughly comparable:
#
# maxtrans                     MC21B
# --------                     -----
# for (k = 0 ; k < n ; k++)    DO 100 JORD=1,N
# while (head >= 0)            DO 70 K=1,JORD
# for (p = Cheap [j] ; ...)    DO 20 II=IN1,IN2
# for (p = head ; ...)         DO 90 K=1,JORD

function augment!(k, Ap, Ai, Match, Cheap, Flag, Istack, Jstack, Pstack, work, maxwork)
    found = false
    head = 1
    Jstack[1] = k
    quick = maxwork > 0
    while head >= 1
        j = Jstack[head]
        pend = Ap[j + 1]

        local i

        if Flag[j] != k
            # prework for node j
            Flag[j] = k
            p = Cheap[j]
            while p < pend && !found
                i = Ai[p]
                found = (Match[i] == EMPTY)
                p += 1
            end
            Cheap[j] = p
            if found
                Istack[head] = i
                break
            end
            Pstack[head] = Ap[j]
        end

        if quick && work > maxwork
            return EMPTY
        end

        pstart = Pstack[head]
        p = pstart
        for outer p ∈ pstart:(pend - 1)
            i = Ai[p]
            j2 = Match[i]
            if Flag[j2] != k
                Pstack[head] = p + 1
                Istack[head] = i
                head += 1
                Jstack[head] = j2
                break
            end
        end
        work += p - pstart + 1

        p == pend && (head -= 1)
    end

    if found
        for p ∈ head:-1:1
            j = Jstack[p]
            i = Istack[p]
            Match[i] = j
        end
    end
    return found, work
end

function maxtrans(nrow, ncol, Ap::Vector{Ti}, Ai::Vector{Ti}, maxwork) where {Ti}
    Cheap = Vector{Ti}(undef, ncol)
    Flag = Vector{Ti}(undef, ncol)
    Istack = Vector{Ti}(undef, ncol)
    Jstack = Vector{Ti}(undef, ncol)
    Pstack = Vector{Ti}(undef, ncol)
    
    Flag .= EMPTY
    Cheap .= view(Ap, 1:ncol)

    Match = fill(EMPTY, nrow)
    maxwork > 0 && (maxwork *= Ap[ncol + 1])

    work = 0

    nmatch = 0
    work_lim_reached = false
    for k ∈ 1:ncol
        result, work = augment!(k, Ap, Ai, Match, Cheap, Flag, Istack, Jstack, Pstack, work, maxwork)
        if result == 1
            nmatch += 1
        elseif result == EMPTY
            work_lim_reached = true
        end
    end

    if work_lim_reached
        work = EMPTY
    end
    return work, Match
end

function dfs!(j, Ap, Ai, Q, Time, Flag, Low, nblocks, timestamp, Cstack, Jstack, Pstack)
    chead = 1
    jhead = 1
    Jstack[1] = j
    @assert Flag[j] == UNVISITED
    while jhead > 0
        j = Jstack[jhead] # grab the node j from the top of Jstack

        # determine which column jj of the A is column j of A*Q
        jj = Q === nothing ? j : unflip(Q[j])
        pend = Ap[jj + 1]
        
        if Flag[j] == UNVISITED
            # node j is visited for the first time
            Cstack[chead += 1] = j # push j onto the stack
            timestamp += 1 # get a timestamp
            Time[j] = timestamp # give the timestamp to j
            Low[j] = timestamp
            Flag[j] = UNASSIGNED # flag j as visited

            # set Pstack [jhead] to the first entry in column j to scan
            Pstack[jhead] = Ap[jj]
        end
        # DFS rooted at node j (start or continue where left off)
        p = Pstack[jhead]
        while p < (pend)
            i = Ai[p]
            if Flag[i] == UNVISITED
                # node i has not been visited, DFS from node i
                # keep track of where we left off in scan of adj list for j
                # so we can restart j where we left off
                Pstack[jhead] = p + 1
                # push i onto the stack and break to recurse on i.
                Jstack[jhead += 1] = i
                @assert Time[i] == EMPTY
                @assert Low[i] == EMPTY
                break
            elseif Flag[i] == UNASSIGNED
                # node i is visited but assigned to a block
                # if Time[i] < Time[j] this is a backedge
                @assert Time[i] > 0
                @assert Low[i] > 0
                Low[j] = min(Low[j], Time[i])
                p += 1
            end
        end
        if p == pend
            # if all adj nodes are visited pop j and do the post work for j
            jhead -= 1

            if Low[j] == Time[j]
                while true
                    @assert chead > 0
                    i = Cstack[chead]
                    chead -= 1
                    @assert i > 0
                    @assert Flag[i] == UNASSIGNED
                    Flag[i] = nblocks
                    i == j && break
                end
                nblocks += 1
            end
            if jhead > 1
                parent = Jstack[jhead]
                Low[parent] = min(Low[parent], Low[j])
            end
        end
    end
    return nblocks, timestamp
end

function strongcomp!(n, Ap::Vector{Ti}, Ai::Vector{Ti}, Q) where {Ti}
    P = fill(Ti(EMPTY), n)
    R = fill(Ti(EMPTY), n + 1)
    Time = fill(Ti(EMPTY), n)
    Flag = fill(Ti(UNVISITED), n)
    Low = P
    Cstack = R
    Jstack = fill(Ti(EMPTY), n)
    Pstack = fill(Ti(EMPTY), n)

    timestamp = 1
    nblocks = 1
    for j ∈ 1:n
        if Flag[j] == UNVISITED
            nblocks, timestamp = dfs!(
                j, Ap, Ai, Q, Time, Flag, Low, nblocks, timestamp,
                Cstack, Jstack, Pstack
            )
        end
    end
    println(Flag)
    R[1:nblocks] .= 1
    for j ∈ 1:n
        R[Flag[j]] += 1
    end
    Time[1] = 1
    cumsum!(R, R)
    R[nblocks] = n

    for j ∈ 1:n
        P[Time[Flag[j]] += 1] = j
    end

    if Q !== nothing
        for k ∈ 1:n
            Time[k] = Q[P[k]]
        end
        Q .= Time
    end
    return nblocks, P, R
end

function order(n, Ap::Vector{Ti}, Ai::Vector{Ti}, maxwork) where {Ti}
    nmatch, Q = maxtrans(n, n, Ap, Ai, maxwork)
    if nmatch < n
        Flag = zeros(Bool, n)
        for i ∈ 1:n
            j = Q[i]
            if j != EMPTY
                Flag[j] = 1
            end
        end
        badcol = Vector{Ti}(undef, n)
        nbadcol = 0
        for j ∈ n:-1:1
            if !Flag[j]
                nbadcol += 1
                badcol[nbadcol] = j
            end
        end

        for i ∈ 1:n
            if Q[i] == EMPTY && nbadcol > 0
                j = badcol[nbadcol]
                nbadcol -= 1
                Q[i] = flip(j)
            end
        end
    end
    nblocks, P, R = strongcomp!(n, Ap, Ai, Q)
    return Q, P, R, nmatch, nblocks
end

end

using Libdl
