#Symmetric and Hermitian matrices
immutable Symmetric{T,S<:AbstractMatrix} <: AbstractMatrix{T}
    data::S
    uplo::Char
end

function Symmetric(A::AbstractMatrix, uplo::Symbol=:U)
    chksquare(A)
    return Symmetric{eltype(A),typeof(A)}(A, char_uplo(uplo))
end

immutable Hermitian{T,S<:AbstractMatrix} <: AbstractMatrix{T}
    data::S
    uplo::Char
end

function Hermitian(A::AbstractMatrix, uplo::Symbol=:U)
    chksquare(A)
    return Hermitian{eltype(A),typeof(A)}(A, char_uplo(uplo))
end

typealias HermOrSym{T,S} Union(Hermitian{T,S},
                               Symmetric{T,S})

typealias RealHermSymComplexHerm{T<:Real,S} Union(Hermitian{T,S},
                                                  Symmetric{T,S},
                                                  Hermitian{Complex{T},S})

size(A::HermOrSym, args...) = size(A.data, args...)

function getindex(A::Symmetric, i::Integer, j::Integer)
    if A.uplo == 'U' && i < j
        return getindex(A.data, i, j)
    else
        return getindex(A.data, j, i)
    end
end

function getindex(A::Hermitian, i::Integer, j::Integer)
    if A.uplo == 'U' && i < j
        return getindex(A.data, i, j)
    else
        return conj(getindex(A.data, j, i))
    end
end

full(A::Symmetric) = copytri!(copy(A.data), A.uplo)
full(A::Hermitian) = copytri!(copy(A.data), A.uplo, true)

convert{T,S<:AbstractMatrix}(::Type{Symmetric{T,S}},A::Symmetric{T,S}) = A

convert{T,S<:AbstractMatrix}(::Type{Symmetric{T,S}},A::Symmetric) =
    Symmetric{T,S}(convert(S,A.data),A.uplo)

convert{T}(::Type{AbstractMatrix{T}}, A::Symmetric) =
    Symmetric(convert(AbstractMatrix{T}, A.data), symbol(A.uplo))

convert{T,S<:AbstractMatrix}(::Type{Hermitian{T,S}},A::Hermitian{T,S}) = A

convert{T,S<:AbstractMatrix}(::Type{Hermitian{T,S}},A::Hermitian) =
    Hermitian{T,S}(convert(S,A.data),A.uplo)

convert{T}(::Type{AbstractMatrix{T}}, A::Hermitian) =
    Hermitian(convert(AbstractMatrix{T}, A.data), symbol(A.uplo))

copy{T,S}(A::Symmetric{T,S}) = Symmetric{T,S}(copy(A.data),A.uplo)
copy{T,S}(A::Hermitian{T,S}) = Hermitian{T,S}(copy(A.data),A.uplo)

ishermitian(A::Hermitian) = true
ishermitian{T<:Real,S}(A::Symmetric{T,S}) = true
ishermitian{T<:Complex,S}(A::Symmetric{T,S}) = all(imag(A.data) .== 0)

issym{T<:Real,S}(A::Hermitian{T,S}) = true
issym{T<:Complex,S}(A::Hermitian{T,S}) = all(imag(A.data) .== 0)
issym(A::Symmetric) = true

transpose(A::Symmetric) = A
ctranspose(A::Hermitian) = A

## Matvec
function A_mul_B!{T<:BlasFloat,S<:AbstractMatrix}(y::StridedVector{T},
                                                  A::Symmetric{T,S},
                                                  x::StridedVector{T})
    return BLAS.symv!(A.uplo, one(T), A.data, x, zero(T), y)
end

function A_mul_B!{T<:BlasComplex,S<:AbstractMatrix}(y::StridedVector{T},
                                                    A::Hermitian{T,S},
                                                    x::StridedVector{T})
    return BLAS.hemv!(A.uplo, one(T), A.data, x, zero(T), y)
end

## Matmat
function A_mul_B!{T<:BlasFloat,S<:AbstractMatrix}(C::StridedMatrix{T},
                                                  A::Symmetric{T,S},
                                                  B::StridedMatrix{T})
    return BLAS.symm!(A.uplo, one(T), A.data, B, zero(T), C)
end

function A_mul_B!{T<:BlasComplex,S<:AbstractMatrix}(y::StridedMatrix{T},
                                                    A::Hermitian{T,S},
                                                    x::StridedMatrix{T})
    return BLAS.hemm!(A.uplo, one(T), A.data, B, zero(T), C)
end

(*)(A::HermOrSym, B::HermOrSym) = full(A) * full(B)
(*)(A::StridedMatrix, B::HermOrSym) = A * full(B)

factorize(A::HermOrSym) = bkfact(A.data, symbol(A.uplo), issym(A))

(\)(A::HermOrSym, B::StridedVecOrMat) =
    (\)(bkfact(A.data, symbol(A.uplo), issym(A)), B)

inv{T<:BlasFloat,S<:StridedMatrix}(A::Hermitian{T,S}) =
    Hermitian{T,S}(inv(bkfact(A.data, symbol(A.uplo))), A.uplo)

inv{T<:BlasFloat,S<:StridedMatrix}(A::Symmetric{T,S}) =
    Symmetric{T,S}(inv(bkfact(A.data, symbol(A.uplo), true)), A.uplo)

eigfact!{T<:BlasReal,S<:StridedMatrix}(A::RealHermSymComplexHerm{T,S}) =
    Eigen(LAPACK.syevr!('V', 'A', A.uplo, A.data, 0.0, 0.0, 0, 0, -1.0)...)

# Because of #6721 it is necessay to specify the parameters explicitly here.
function eigfact{T1<:Real,T2}(A::RealHermSymComplexHerm{T1,T2})
    T = eltype(A)
    S = promote_type(Float32, typeof(zero(T) / norm(one(T))))
    if S == T
        return eigfact!(copy(A))
    else
        return eigfact!(convert(AbstractMatrix{S}, A))
    end
end

function eigfact!{T<:BlasReal,S<:StridedMatrix}(A::RealHermSymComplexHerm{T,S},
                                                irange::UnitRange)
    return Eigen(LAPACK.syevr!('V', 'I', A.uplo, A.data,
                               0.0, 0.0, irange.start, irange.stop, -1.0)...)
end

# Because of #6721 it is necessay to specify the parameters explicitly here.
function eigfact{T1<:Real,T2}(A::RealHermSymComplexHerm{T1,T2},
                              irange::UnitRange)
    T = eltype(A)
    S = promote_type(Float32, typeof(zero(T) / norm(one(T))))
    if S == T
        return eigfact!(copy(A), irange)
    else
        return eigfact!(convert(AbstractMatrix{S}, A), irange)
    end
end

function eigfact!{T<:BlasReal,S<:StridedMatrix}(A::RealHermSymComplexHerm{T,S},
                                                vl::Real, vh::Real)
    return Eigen(LAPACK.syevr!('V', 'V', A.uplo, A.data,
                               convert(T, vl), convert(T, vh), 0, 0, -1.0)...)
end

# Because of #6721 it is necessay to specify the parameters explicitly here.
function eigfact{T1<:Real,T2}(A::RealHermSymComplexHerm{T1,T2}, vl::Real, vh::Real)
    T = eltype(A)
    S = promote_type(Float32, typeof(zero(T) / norm(one(T))))
    if S == T
        return eigfact!(copy(A), vl, vh)
    else
        return eigfact!(convert(AbstractMatrix{S}, A), vl, vh)
    end
end

eigvals!{T<:BlasReal,S<:StridedMatrix}(A::RealHermSymComplexHerm{T,S}) =
    LAPACK.syevr!('N', 'A', A.uplo, A.data, 0.0, 0.0, 0, 0, -1.0)[1]

# Because of #6721 it is necessay to specify the parameters explicitly here.
function eigvals{T1<:Real,T2}(A::RealHermSymComplexHerm{T1,T2})
    T = eltype(A)
    S = promote_type(Float32, typeof(zero(T) / norm(one(T))))
    if S == T
        return eigvals!(copy(A))
    else
        return eigvals!(convert(AbstractMatrix{S}, A))
    end
end

function eigvals!{T<:BlasReal,S<:StridedMatrix}(A::RealHermSymComplexHerm{T,S},
                                                irange::UnitRange)
    return LAPACK.syevr!('N', 'I', A.uplo, A.data,
                         0.0, 0.0, irange.start, irange.stop, -1.0)[1]
end

# Because of #6721 it is necessay to specify the parameters explicitly here.
function eigvals{T1<:Real,T2}(A::RealHermSymComplexHerm{T1,T2},
                              irange::UnitRange)
    T = eltype(A)
    S = promote_type(Float32, typeof(zero(T) / norm(one(T))))
    if S == T
        return eigvals!(copy(A), irange)
    else
        return eigvals!(convert(AbstractMatrix{S}, A), irange)
    end
end

function eigvals!{T<:BlasReal,S<:StridedMatrix}(A::RealHermSymComplexHerm{T,S},
                                                vl::Real, vh::Real)
    return LAPACK.syevr!('N', 'V', A.uplo, A.data,
                         convert(T, vl), convert(T, vh), 0, 0, -1.0)[1]
end

# Because of #6721 it is necessay to specify the parameters explicitly here.
function eigvals{T1<:Real,T2}(A::RealHermSymComplexHerm{T1,T2},
                              vl::Real, vh::Real)
    T = eltype(A)
    S = promote_type(Float32, typeof(zero(T) / norm(one(T))))
    if S == T
        return eigvals!(copy(A), vl, vh)
    else
        return eigvals!(convert(AbstractMatrix{S}, A), vl, vh)
    end
end

eigmax{T<:Real,S<:StridedMatrix}(A::RealHermSymComplexHerm{T,S}) =
    eigvals(A, size(A, 1):size(A, 1))[1]

eigmin{T<:Real,S<:StridedMatrix}(A::RealHermSymComplexHerm{T,S}) =
    eigvals(A, 1:1)[1]

function eigfact!{T<:BlasReal,S<:StridedMatrix}(A::HermOrSym{T,S},
                                                B::HermOrSym{T,S})
    if A.uplo == B.uplo
        vals, vecs, _ = LAPACK.sygvd!(1, 'V', A.uplo, A.data, B.data)
    else
        vals, vecs, _ = LAPACK.sygvd!(1, 'V', A.uplo, A.data, B.data')
    end
    return GeneralizedEigen(vals, vecs)
end

function eigfact!{T<:BlasComplex,S<:StridedMatrix}(A::Hermitian{T,S},
                                                   B::Hermitian{T,S})
    if A.uplo == B.uplo
        vals, vecs, _ = LAPACK.sygvd!(1, 'V', A.uplo, A.data, B.data)
    else
        vals, vecs, _ = LAPACK.sygvd!(1, 'V', A.uplo, A.data, B.data')
    end
    GeneralizedEigen(vals, vecs)
end

function eigvals!{T<:BlasReal,S<:StridedMatrix}(A::HermOrSym{T,S},
                                                B::HermOrSym{T,S})
    if A.uplo == B.uplo
        return LAPACK.sygvd!(1, 'N', A.uplo, A.data, B.data)[1]
    else
        return LAPACK.sygvd!(1, 'N', A.uplo, A.data, B.data')[1]
    end
end

function eigvals!{T<:BlasComplex,S<:StridedMatrix}(A::Hermitian{T,S},
                                                   B::Hermitian{T,S})
    if A.uplo == B.uplo
        return LAPACK.sygvd!(1, 'N', A.uplo, A.data, B.data)[1]
    else
        return LAPACK.sygvd!(1, 'N', A.uplo, A.data, B.data')[1]
    end
end

#Matrix-valued functions
function expm{T<:Real}(A::RealHermSymComplexHerm{T})
    F = eigfact(A)
    return F.vectors * Diagonal(exp(F.values)) * F.vectors'
end

function sqrtm{T<:Real}(A::RealHermSymComplexHerm{T})
    F = eigfact(A)
    if isposdef(F)
        return F.vectors * Diagonal(sqrt(F.values)) * F.vectors'
    else
        return F.vectors * Diagonal(sqrt(complex(F.values))) * F.vectors'
    end
end
