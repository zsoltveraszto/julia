## Create an extractor that extracts the modified original matrix, e.g.
## LD for BunchKaufman, UL for CholeskyDense, LU for LUDense and
## define size methods for Factorization types using it.

immutable BunchKaufman{T,S<:AbstractMatrix} <: Factorization{T}
    LD::S
    ipiv::Vector{BlasInt}
    uplo::Char
    symmetric::Bool

    function BunchKaufman(LD::AbstractMatrix{T},
                          ipiv::Vector{BlasInt},
                          uplo::Char,
                          symmetric::Bool)
        return new(LD, ipiv, uplo, symmetric)
    end
end

function BunchKaufman{T}(LD::AbstractMatrix{T},
                         ipiv::Vector{BlasInt},
                         uplo::Char,
                         symmetric::Bool)
    BunchKaufman{T,typeof(LD)}(LD, ipiv, uplo, symmetric)
end

function bkfact!{T<:BlasReal}(A::StridedMatrix{T},
                              uplo::Symbol=:U,
                              symmetric::Bool=issym(A))
    if !symmetric
        throw(ArgumentError(
            "Bunch-Kaufman decomposition is only valid for symmetric real matrices"))
    end
    LD, ipiv = LAPACK.sytrf!(char_uplo(uplo) , A)
    return BunchKaufman(LD, ipiv, char_uplo(uplo), symmetric)
end

function bkfact!{T<:BlasComplex}(A::StridedMatrix{T},
                                 uplo::Symbol=:U,
                                 symmetric::Bool=issym(A))
    if symmetric
        LD, ipiv = LAPACK.sytrf!(char_uplo(uplo), A)
    else
        LD, ipiv = LAPACK.hetrf!(char_uplo(uplo) , A)
    end
    return BunchKaufman(LD, ipiv, char_uplo(uplo), symmetric)
end

bkfact{T<:BlasFloat}(A::StridedMatrix{T}, uplo::Symbol=:U, symmetric::Bool=issym(A)) =
    bkfact!(copy(A), uplo, symmetric)

function bkfact{T}(A::StridedMatrix{T}, uplo::Symbol=:U, symmetric::Bool=issym(A))
    MT = promote_type(Float32, typeof(sqrt(one(T))))
    return bkfact!(convert(Matrix{MT}, A), uplo, symmetric)
end

convert{T}(::Type{BunchKaufman{T}}, B::BunchKaufman) =
    BunchKaufman(convert(Matrix{T}, B.LD), B.ipiv, B.uplo, B.symmetric)

convert{T}(::Type{Factorization{T}}, B::BunchKaufman) =
    convert(BunchKaufman{T}, B)

size(B::BunchKaufman) = size(B.LD)
size(B::BunchKaufman, d::Integer) = size(B.LD,d)

issym(B::BunchKaufman) = B.symmetric
ishermitian(B::BunchKaufman) = !B.symmetric

inv{T<:BlasReal}(B::BunchKaufman{T}) =
    copytri!(LAPACK.sytri!(B.uplo, copy(B.LD), B.ipiv), B.uplo, true)

function inv{T<:BlasComplex}(B::BunchKaufman{T})
    if issym(B)
        return copytri!(LAPACK.sytri!(B.uplo, copy(B.LD), B.ipiv), B.uplo)
    else
        return copytri!(LAPACK.hetri!(B.uplo, copy(B.LD), B.ipiv), B.uplo, true)
    end
end

A_ldiv_B!{T<:BlasReal}(B::BunchKaufman{T}, R::StridedVecOrMat{T}) =
    LAPACK.sytrs!(B.uplo, B.LD, B.ipiv, R)

function A_ldiv_B!{T<:BlasComplex}(B::BunchKaufman{T}, R::StridedVecOrMat{T})
    if issym(B)
        return LAPACK.sytrs!(B.uplo, B.LD, B.ipiv, R)
    else
        return LAPACK.hetrs!(B.uplo, B.LD, B.ipiv, R)
    end
end
