module LibBTF
using SuiteSparse_jll

libbtf = SuiteSparse_jll.get_libbtf_path()
struct SuiteSparse_config_struct
    malloc_func::Ptr{Cvoid}
    calloc_func::Ptr{Cvoid}
    realloc_func::Ptr{Cvoid}
    free_func::Ptr{Cvoid}
    printf_func::Ptr{Cvoid}
    hypot_func::Ptr{Cvoid}
    divcomplex_func::Ptr{Cvoid}
end

function SuiteSparse_start()
    ccall((:SuiteSparse_start, libbtf), Cvoid, ())
end

function SuiteSparse_finish()
    ccall((:SuiteSparse_finish, libbtf), Cvoid, ())
end

function SuiteSparse_malloc(nitems, size_of_item)
    ccall((:SuiteSparse_malloc, libbtf), Ptr{Cvoid}, (Csize_t, Csize_t), nitems, size_of_item)
end

function SuiteSparse_calloc(nitems, size_of_item)
    ccall((:SuiteSparse_calloc, libbtf), Ptr{Cvoid}, (Csize_t, Csize_t), nitems, size_of_item)
end

function SuiteSparse_realloc(nitems_new, nitems_old, size_of_item, p, ok)
    ccall((:SuiteSparse_realloc, libbtf), Ptr{Cvoid}, (Csize_t, Csize_t, Csize_t, Ptr{Cvoid}, Ptr{Cint}), nitems_new, nitems_old, size_of_item, p, ok)
end

function SuiteSparse_free(p)
    ccall((:SuiteSparse_free, libbtf), Ptr{Cvoid}, (Ptr{Cvoid},), p)
end

function SuiteSparse_tic(tic)
    ccall((:SuiteSparse_tic, libbtf), Cvoid, (Ptr{Cdouble},), tic)
end

function SuiteSparse_toc(tic)
    ccall((:SuiteSparse_toc, libbtf), Cdouble, (Ptr{Cdouble},), tic)
end

function SuiteSparse_time()
    ccall((:SuiteSparse_time, libbtf), Cdouble, ())
end

function SuiteSparse_hypot(x, y)
    ccall((:SuiteSparse_hypot, libbtf), Cdouble, (Cdouble, Cdouble), x, y)
end

function SuiteSparse_divcomplex(ar, ai, br, bi, cr, ci)
    ccall((:SuiteSparse_divcomplex, libbtf), Cint, (Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), ar, ai, br, bi, cr, ci)
end

function SuiteSparse_version(version)
    ccall((:SuiteSparse_version, libbtf), Cint, (Ptr{Cint},), version)
end

function btf_maxtrans(nrow, ncol, Ap, Ai, maxwork, work, Match, Work)
    ccall((:btf_maxtrans, libbtf), Cint, (Cint, Cint, Ptr{Cint}, Ptr{Cint}, Cdouble, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}), nrow, ncol, Ap, Ai, maxwork, work, Match, Work)
end

function btf_l_maxtrans(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
    ccall((:btf_l_maxtrans, libbtf), Clong, (Clong, Clong, Ptr{Clong}, Ptr{Clong}, Cdouble, Ptr{Cdouble}, Ptr{Clong}, Ptr{Clong}), arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8)
end

function btf_strongcomp(n, Ap, Ai, Q, P, R, Work)
    ccall((:btf_strongcomp, libbtf), Cint, (Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), n, Ap, Ai, Q, P, R, Work)
end

function btf_l_strongcomp(arg1, arg2, arg3, arg4, arg5, arg6, arg7)
    ccall((:btf_l_strongcomp, libbtf), Clong, (Clong, Ptr{Clong}, Ptr{Clong}, Ptr{Clong}, Ptr{Clong}, Ptr{Clong}, Ptr{Clong}), arg1, arg2, arg3, arg4, arg5, arg6, arg7)
end

function btf_order(n, Ap, Ai, maxwork, work, P, Q, R, nmatch, Work)
    ccall((:btf_order, libbtf), Cint, (Cint, Ptr{Cint}, Ptr{Cint}, Cdouble, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), n, Ap, Ai, maxwork, work, P, Q, R, nmatch, Work)
end

function btf_l_order(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10)
    ccall((:btf_l_order, libbtf), Clong, (Clong, Ptr{Clong}, Ptr{Clong}, Cdouble, Ptr{Cdouble}, Ptr{Clong}, Ptr{Clong}, Ptr{Clong}, Ptr{Clong}, Ptr{Clong}), arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10)
end

end # module
