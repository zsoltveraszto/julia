## Triangular

# could be renamed to Triangular when than name has been fully deprecated
abstract AbstractTriangular{T,S<:AbstractMatrix} <: AbstractMatrix{T}

# Methods that don't need special care for upper/lower and unit diagonal

## Lower Triangular

immutable LowerTriangular{T,S<:AbstractMatrix} <: AbstractTriangular{T,S}
    data::S
end

function LowerTriangular(A::AbstractMatrix)
    chksquare(A)
    return LowerTriangular{eltype(A), typeof(A)}(A)
end

size(A::LowerTriangular, args...) = size(A.data, args...)

convert{T,S}(::Type{LowerTriangular{T}}, A::LowerTriangular{T,S}) = A

function convert{Tnew,Told,S}(::Type{LowerTriangular{Tnew}},
                              A::LowerTriangular{Told,S})
    Anew = convert(AbstractMatrix{Tnew}, A.data)
    return LowerTriangular(Anew)
end

convert{Tnew,Told,S}(::Type{AbstractMatrix{Tnew}}, A::LowerTriangular{Told,S}) =
    convert(LowerTriangular{Tnew}, A)

convert{T,S}(::Type{Matrix}, A::LowerTriangular{T,S}) =
    convert(Matrix{T}, A)

function similar{T,S,Tnew}(A::LowerTriangular{T,S}, ::Type{Tnew}, dims::Dims)
    if dims[1] != dims[2]
        throw(ArgumentError("Triangular matrix must be square"))
    end
    if length(dims) != 2
        throw(ArgumentError("Triangular matrix must have two dimensions"))
    end
    B = similar(A.data, Tnew, dims)
    return LowerTriangular(B)
end

copy{T,S}(A::LowerTriangular{T,S}) = LowerTriangular{T,S}(copy(A.data))
big(A::LowerTriangular) = LowerTriangular(big(A.data))

real{T<:Real}(A::LowerTriangular{T}) = A
function real{T<:Complex}(A::LowerTriangular{T})
    B = real(A.data)
    return LowerTriangular(B)
end

## Unit Lower Triangular

immutable UnitLowerTriangular{T,S<:AbstractMatrix} <: AbstractTriangular{T,S}
    data::S
end
function UnitLowerTriangular(A::AbstractMatrix)
    chksquare(A)
    return UnitLowerTriangular{eltype(A), typeof(A)}(A)
end

size(A::UnitLowerTriangular, args...) = size(A.data, args...)

convert{T,S}(::Type{UnitLowerTriangular{T}}, A::UnitLowerTriangular{T,S}) = A

function convert{Tnew,Told,S}(::Type{UnitLowerTriangular{Tnew}},
                              A::UnitLowerTriangular{Told,S})
    Anew = convert(AbstractMatrix{Tnew}, A.data)
    return UnitLowerTriangular(Anew)
end

convert{Tnew,Told,S}(::Type{AbstractMatrix{Tnew}}, A::UnitLowerTriangular{Told,S}) =
    convert(UnitLowerTriangular{Tnew}, A)

convert{T,S}(::Type{Matrix}, A::UnitLowerTriangular{T,S}) =
    convert(Matrix{T}, A)

function similar{T,S,Tnew}(A::UnitLowerTriangular{T,S}, ::Type{Tnew}, dims::Dims)
    if dims[1] != dims[2]
        throw(ArgumentError("Triangular matrix must be square"))
    end
    if length(dims) != 2
        throw(ArgumentError("Triangular matrix must have two dimensions"))
    end
    B = similar(A.data, Tnew, dims)
    return UnitLowerTriangular(B)
end

copy{T,S}(A::UnitLowerTriangular{T,S}) = UnitLowerTriangular{T,S}(copy(A.data))

big(A::UnitLowerTriangular) = UnitLowerTriangular(big(A.data))

real{T<:Real}(A::UnitLowerTriangular{T}) = A

function real{T<:Complex}(A::UnitLowerTriangular{T})
    B = real(A.data)
    return UnitLowerTriangular(B)
end

## Upper Triangular

immutable UpperTriangular{T,S<:AbstractMatrix} <: AbstractTriangular{T,S}
    data::S
end

function UpperTriangular(A::AbstractMatrix)
    chksquare(A)
    return UpperTriangular{eltype(A), typeof(A)}(A)
end

size(A::UpperTriangular, args...) = size(A.data, args...)
convert{T,S}(::Type{UpperTriangular{T}}, A::UpperTriangular{T,S}) = A

function convert{Tnew,Told,S}(::Type{UpperTriangular{Tnew}}, A::UpperTriangular{Told,S})
    Anew = convert(AbstractMatrix{Tnew}, A.data)
    return UpperTriangular(Anew)
end

convert{Tnew,Told,S}(::Type{AbstractMatrix{Tnew}}, A::UpperTriangular{Told,S}) =
    convert(UpperTriangular{Tnew}, A)

convert{T,S}(::Type{Matrix}, A::UpperTriangular{T,S}) =
    convert(Matrix{T}, A)

function similar{T,S,Tnew}(A::UpperTriangular{T,S}, ::Type{Tnew}, dims::Dims)
    if dims[1] != dims[2]
        throw(ArgumentError("Triangular matrix must be square"))
    end
    if length(dims) != 2
        throw(ArgumentError("Triangular matrix must have two dimensions"))
    end
    B = similar(A.data, Tnew, dims)
    return UpperTriangular(B)
end

copy{T,S}(A::UpperTriangular{T,S}) = UpperTriangular{T,S}(copy(A.data))
big(A::UpperTriangular) = UpperTriangular(big(A.data))

real{T<:Real}(A::UpperTriangular{T}) = A
real{T<:Complex}(A::UpperTriangular{T}) = (B = real(A.data); UpperTriangular(B))


## Unit Upper Triangular

immutable UnitUpperTriangular{T,S<:AbstractMatrix} <: AbstractTriangular{T,S}
    data::S
end

function UnitUpperTriangular(A::AbstractMatrix)
    chksquare(A)
    return UnitUpperTriangular{eltype(A), typeof(A)}(A)
end

size(A::UnitUpperTriangular, args...) = size(A.data, args...)

convert{T,S}(::Type{UnitUpperTriangular{T}}, A::UnitUpperTriangular{T,S}) = A

function convert{Tnew,Told,S}(::Type{UnitUpperTriangular{Tnew}},
                              A::UnitUpperTriangular{Told,S})
    Anew = convert(AbstractMatrix{Tnew}, A.data)
    return UnitUpperTriangular(Anew)
end

convert{Tnew,Told,S}(::Type{AbstractMatrix{Tnew}}, A::UnitUpperTriangular{Told,S}) =
    convert(UnitUpperTriangular{Tnew}, A)

convert{T,S}(::Type{Matrix}, A::UnitUpperTriangular{T,S}) =
    convert(Matrix{T}, A)

function similar{T,S,Tnew}(A::UnitUpperTriangular{T,S}, ::Type{Tnew}, dims::Dims)
    if dims[1] != dims[2]
        throw(ArgumentError("Triangular matrix must be square"))
    end
    if length(dims) != 2
        throw(ArgumentError("Triangular matrix must have two dimensions"))
    end
    B = similar(A.data, Tnew, dims)
    return UnitUpperTriangular(B)
end

copy{T,S}(A::UnitUpperTriangular{T,S}) = UnitUpperTriangular{T,S}(copy(A.data))

big(A::UnitUpperTriangular) = UnitUpperTriangular(big(A.data))

real{T<:Real}(A::UnitUpperTriangular{T}) = A

function real{T<:Complex}(A::UnitUpperTriangular{T})
    B = real(A.data)
    return UnitUpperTriangular(B)
end

full(A::AbstractTriangular) = convert(Matrix, A)

fill!(A::AbstractTriangular, x) = (fill!(A.data, x); A)

# then handle all methods that requires specific handling of upper/lower and unit diagonal

function convert{Tret,T,S}(::Type{Matrix{Tret}}, A::LowerTriangular{T,S})
    B = Array(Tret, size(A, 1), size(A, 1))
    copy!(B, A.data)
    tril!(B)
    return B
end

function convert{Tret,T,S}(::Type{Matrix{Tret}}, A::UnitLowerTriangular{T,S})
    B = Array(Tret, size(A, 1), size(A, 1))
    copy!(B, A.data)
    tril!(B)
    for i = 1:size(B,1)
        B[i,i] = 1
    end
    return B
end

function convert{Tret,T,S}(::Type{Matrix{Tret}}, A::UpperTriangular{T,S})
    B = Array(Tret, size(A, 1), size(A, 1))
    copy!(B, A.data)
    triu!(B)
    return B
end

function convert{Tret,T,S}(::Type{Matrix{Tret}}, A::UnitUpperTriangular{T,S})
    B = Array(Tret, size(A, 1), size(A, 1))
    copy!(B, A.data)
    triu!(B)
    for i = 1:size(B,1)
        B[i,i] = 1
    end
    return B
end

function full!{T,S}(A::LowerTriangular{T,S})
    B = A.data
    tril!(B)
    return B
end

function full!{T,S}(A::UnitLowerTriangular{T,S})
    B = A.data
    tril!(B)
    for i = 1:size(A,1)
        B[i,i] = 1
    end
    return B
end

function full!{T,S}(A::UpperTriangular{T,S})
    B = A.data
    triu!(B)
    return B
end

function full!{T,S}(A::UnitUpperTriangular{T,S})
    B = A.data
    triu!(B)
    for i = 1:size(A,1)
        B[i,i] = 1
    end
    return B
end

function getindex(A::AbstractTriangular, i::Integer)
    m, n = divrem(i - 1, size(A, 1))
    return A[n + 1, m + 1]
end

function getindex{T,S}(A::UnitLowerTriangular{T,S}, i::Integer, j::Integer)
    if i == j
        return one(T)
    elseif i > j
        return A.data[i,j]
    else
        return zero(A.data[i,j])
    end
end

function getindex{T,S}(A::LowerTriangular{T,S}, i::Integer, j::Integer)
    if i >= j
        return A.data[i,j]
    else
        return zero(A.data[i,j])
    end
end

function getindex{T,S}(A::UnitUpperTriangular{T,S}, i::Integer, j::Integer)
    if i == j
        return one(T)
    elseif i < j
        return A.data[i,j]
    else
        return zero(A.data[i,j])
    end
end

function getindex{T,S}(A::UpperTriangular{T,S}, i::Integer, j::Integer)
    if i <= j
        return A.data[i,j]
    else
        return zero(A.data[i,j])
    end
end

function setindex!(A::UpperTriangular, x, i::Integer, j::Integer)
    if i > j
        throw(BoundsError())
    end
    A.data[i,j] = x
    return A
end

function setindex!(A::UnitUpperTriangular, x, i::Integer, j::Integer)
    if i >= j
        throw(BoundsError())
    end
    A.data[i,j] = x
    return A
end

function setindex!(A::LowerTriangular, x, i::Integer, j::Integer)
    if i < j
        throw(BoundsError())
    end
    A.data[i,j] = x
    return A
end

function setindex!(A::UnitLowerTriangular, x, i::Integer, j::Integer)
    if i <= j
        throw(BoundsError())
    end
    A.data[i,j] = x
    return A
end

istril{T,S}(A::LowerTriangular{T,S}) = true
istril{T,S}(A::UnitLowerTriangular{T,S}) = true
istril{T,S}(A::UpperTriangular{T,S}) = false
istril{T,S}(A::UnitUpperTriangular{T,S}) = false
istriu{T,S}(A::LowerTriangular{T,S}) = false
istriu{T,S}(A::UnitLowerTriangular{T,S}) = false
istriu{T,S}(A::UpperTriangular{T,S}) = true
istriu{T,S}(A::UnitUpperTriangular{T,S}) = true

transpose{T,S}(A::LowerTriangular{T,S}) = UpperTriangular{T, S}(transpose(A.data))
transpose{T,S}(A::UnitLowerTriangular{T,S}) = UnitUpperTriangular{T, S}(transpose(A.data))
transpose{T,S}(A::UpperTriangular{T,S}) = LowerTriangular{T, S}(transpose(A.data))
transpose{T,S}(A::UnitUpperTriangular{T,S}) = UnitLowerTriangular{T, S}(transpose(A.data))
ctranspose{T,S}(A::LowerTriangular{T,S}) = UpperTriangular{T, S}(ctranspose(A.data))
ctranspose{T,S}(A::UnitLowerTriangular{T,S}) = UnitUpperTriangular{T, S}(ctranspose(A.data))
ctranspose{T,S}(A::UpperTriangular{T,S}) = LowerTriangular{T, S}(ctranspose(A.data))
ctranspose{T,S}(A::UnitUpperTriangular{T,S}) = UnitLowerTriangular{T, S}(ctranspose(A.data))

transpose!{T,S}(A::LowerTriangular{T,S}) = UpperTriangular{T, S}(copytri!(A.data, 'L'))
transpose!{T,S}(A::UnitLowerTriangular{T,S}) = UnitUpperTriangular{T, S}(copytri!(A.data, 'L'))
transpose!{T,S}(A::UpperTriangular{T,S}) = LowerTriangular{T, S}(copytri!(A.data, 'U'))
transpose!{T,S}(A::UnitUpperTriangular{T,S}) = UnitLowerTriangular{T, S}(copytri!(A.data, 'U'))
ctranspose!{T,S}(A::LowerTriangular{T,S}) = UpperTriangular{T, S}(copytri!(A.data, 'L' , true))
ctranspose!{T,S}(A::UnitLowerTriangular{T,S}) = UnitUpperTriangular{T, S}(copytri!(A.data, 'L' , true))
ctranspose!{T,S}(A::UpperTriangular{T,S}) = LowerTriangular{T, S}(copytri!(A.data, 'U' , true))
ctranspose!{T,S}(A::UnitUpperTriangular{T,S}) = UnitLowerTriangular{T, S}(copytri!(A.data, 'U' , true))

diag(A::LowerTriangular) = diag(A.data)
diag(A::UnitLowerTriangular) = ones(eltype(A), size(A,1))
diag(A::UpperTriangular) = diag(A.data)
diag(A::UnitUpperTriangular) = ones(eltype(A), size(A,1))

# Unary operations
(-)(A::LowerTriangular) = LowerTriangular(-A.data)
(-)(A::UpperTriangular) = UpperTriangular(-A.data)

function (-)(A::UnitLowerTriangular)
    Anew = -A.data
    for i = 1:size(A, 1)
        Anew[i, i] = -1
    end
    return LowerTriangular(Anew)
end

function (-)(A::UnitUpperTriangular)
    Anew = -A.data
    for i = 1:size(A, 1)
        Anew[i, i] = -1
    end
    return UpperTriangular(Anew)
end

# Binary operations
(+)(A::UpperTriangular, B::UpperTriangular) = UpperTriangular(A.data + B.data)
(+)(A::LowerTriangular, B::LowerTriangular) = LowerTriangular(A.data + B.data)
(+)(A::UpperTriangular, B::UnitUpperTriangular) = UpperTriangular(A.data + triu(B.data, 1) + I)
(+)(A::LowerTriangular, B::UnitLowerTriangular) = LowerTriangular(A.data + tril(B.data, -1) + I)
(+)(A::UnitUpperTriangular, B::UpperTriangular) = UpperTriangular(triu(A.data, 1) + B.data + I)
(+)(A::UnitLowerTriangular, B::LowerTriangular) = LowerTriangular(tril(A.data, -1) + B.data + I)
(+)(A::UnitUpperTriangular, B::UnitUpperTriangular) = UpperTriangular(triu(A.data, 1) + triu(B.data, 1) + 2I)
(+)(A::UnitLowerTriangular, B::UnitLowerTriangular) = LowerTriangular(tril(A.data, -1) + tril(B.data, -1) + 2I)
(+)(A::AbstractTriangular, B::AbstractTriangular) = full(A) + full(B)

(-)(A::UpperTriangular, B::UpperTriangular) = UpperTriangular(A.data - B.data)
(-)(A::LowerTriangular, B::LowerTriangular) = LowerTriangular(A.data - B.data)
(-)(A::UpperTriangular, B::UnitUpperTriangular) = UpperTriangular(A.data - triu(B.data, 1) - I)
(-)(A::LowerTriangular, B::UnitLowerTriangular) = LowerTriangular(A.data - tril(B.data, -1) - I)
(-)(A::UnitUpperTriangular, B::UpperTriangular) = UpperTriangular(triu(A.data, 1) - B.data + I)
(-)(A::UnitLowerTriangular, B::LowerTriangular) = LowerTriangular(tril(A.data, -1) - B.data + I)
(-)(A::UnitUpperTriangular, B::UnitUpperTriangular) = UpperTriangular(triu(A.data, 1) - triu(B.data, 1))
(-)(A::UnitLowerTriangular, B::UnitLowerTriangular) = LowerTriangular(tril(A.data, -1) - tril(B.data, -1))
(-)(A::AbstractTriangular, B::AbstractTriangular) = full(A) - full(B)

######################
# BlasFloat routines #
######################

A_mul_B!(A::Tridiagonal, B::AbstractTriangular) = A * full!(B)

A_mul_B!(C::AbstractVecOrMat, A::AbstractTriangular, B::AbstractVecOrMat) =
    A_mul_B!(A, copy!(C, B))

A_mul_Bc!(C::AbstractVecOrMat, A::AbstractTriangular, B::AbstractVecOrMat) =
    A_mul_Bc!(A, copy!(C, B))

# Vector multiplication
A_mul_B!{T<:BlasFloat,S<:StridedMatrix}(A::LowerTriangular{T,S}, b::StridedVector{T}) =
    BLAS.trmv!('L', 'N', 'N', A.data, b)

# Matrix multiplication
A_mul_B!{T<:BlasFloat,S<:StridedMatrix}(A::LowerTriangular{T,S}, B::StridedMatrix{T}) =
    BLAS.trmm!('L', 'L', 'N', 'N', one(T), A.data, B)

A_mul_B!{T<:BlasFloat,S<:StridedMatrix}(A::StridedMatrix{T}, B::LowerTriangular{T,S}) =
    BLAS.trmm!('R', 'L', 'N', 'N', one(T), B.data, A)

Ac_mul_B!{T<:BlasComplex,S<:StridedMatrix}(A::LowerTriangular{T,S}, B::StridedMatrix{T}) =
    BLAS.trmm!('L', 'L', 'C', 'N', one(T), A.data, B)

Ac_mul_B!{T<:BlasReal,S<:StridedMatrix}(A::LowerTriangular{T,S}, B::StridedMatrix{T}) =
    BLAS.trmm!('L', 'L', 'T', 'N', one(T), A.data, B)

A_mul_Bc!{T<:BlasComplex,S<:StridedMatrix}(A::StridedMatrix{T}, B::LowerTriangular{T,S}) =
    BLAS.trmm!('R', 'L', 'C', 'N', one(T), B.data, A)

A_mul_Bc!{T<:BlasReal,S<:StridedMatrix}(A::StridedMatrix{T}, B::LowerTriangular{T,S}) =
    BLAS.trmm!('R', 'L', 'T', 'N', one(T), B.data, A)

# Left division
A_ldiv_B!{T<:BlasFloat,S<:StridedMatrix}(A::LowerTriangular{T,S}, B::StridedVecOrMat{T}) =
    LAPACK.trtrs!('L', 'N', 'N', A.data, B)

Ac_ldiv_B!{T<:BlasReal,S<:StridedMatrix}(A::LowerTriangular{T,S}, B::StridedVecOrMat{T}) =
    LAPACK.trtrs!('L', 'T', 'N', A.data, B)

Ac_ldiv_B!{T<:BlasComplex,S<:StridedMatrix}(A::LowerTriangular{T,S}, B::StridedVecOrMat{T}) =
    LAPACK.trtrs!('L', 'C', 'N', A.data, B)

# Right division
A_rdiv_B!{T<:BlasFloat,S<:StridedMatrix}(A::StridedVecOrMat{T}, B::LowerTriangular{T,S}) =
    BLAS.trsm!('R', 'L', 'N', 'N', one(T), B.data, A)

A_rdiv_Bc!{T<:BlasReal,S<:StridedMatrix}(A::StridedMatrix{T}, B::LowerTriangular{T,S}) =
    BLAS.trsm!('R', 'L', 'T', 'N', one(T), B.data, A)

A_rdiv_Bc!{T<:BlasComplex,S<:StridedMatrix}(A::StridedMatrix{T}, B::LowerTriangular{T,S}) =
    BLAS.trsm!('R', 'L', 'C', 'N', one(T), B.data, A)

# Matrix inverse
inv{T<:BlasFloat,S<:StridedMatrix}(A::LowerTriangular{T,S}) =
    LowerTriangular{T,S}(LAPACK.trtri!('L', 'N', copy(A.data)))

# Error bounds for triangular solve
function errorbounds{T<:BlasFloat,S<:StridedMatrix}(A::LowerTriangular{T,S},
                                                    X::StridedVecOrMat{T},
                                                    B::StridedVecOrMat{T})
    return LAPACK.trrfs!('L', 'N', 'N', A.data, B, X)
end

# Condition numbers
function cond{T<:BlasFloat,S}(A::LowerTriangular{T,S}, p::Real=2)
    chksquare(A)
    if p==1
        return inv(LAPACK.trcon!('O', 'L', 'N', A.data))
    elseif p==Inf
        return inv(LAPACK.trcon!('I', 'L', 'N', A.data))
    else #use fallback
        return cond(full(A), p)
    end
end

# Vector multiplication
A_mul_B!{T<:BlasFloat,S<:StridedMatrix}(A::UnitLowerTriangular{T,S}, b::StridedVector{T}) =
    BLAS.trmv!('L', 'N', 'U', A.data, b)

# Matrix multiplication
A_mul_B!{T<:BlasFloat,S<:StridedMatrix}(A::UnitLowerTriangular{T,S}, B::StridedMatrix{T}) =
    BLAS.trmm!('L', 'L', 'N', 'U', one(T), A.data, B)

A_mul_B!{T<:BlasFloat,S<:StridedMatrix}(A::StridedMatrix{T}, B::UnitLowerTriangular{T,S}) =
    BLAS.trmm!('R', 'L', 'N', 'U', one(T), B.data, A)

Ac_mul_B!{T<:BlasComplex,S<:StridedMatrix}(A::UnitLowerTriangular{T,S}, B::StridedMatrix{T}) =
    BLAS.trmm!('L', 'L', 'C', 'U', one(T), A.data, B)

Ac_mul_B!{T<:BlasReal,S<:StridedMatrix}(A::UnitLowerTriangular{T,S}, B::StridedMatrix{T}) =
    BLAS.trmm!('L', 'L', 'T', 'U', one(T), A.data, B)

A_mul_Bc!{T<:BlasComplex,S<:StridedMatrix}(A::StridedMatrix{T}, B::UnitLowerTriangular{T,S}) =
    BLAS.trmm!('R', 'L', 'C', 'U', one(T), B.data, A)

A_mul_Bc!{T<:BlasReal,S<:StridedMatrix}(A::StridedMatrix{T}, B::UnitLowerTriangular{T,S}) =
    BLAS.trmm!('R', 'L', 'T', 'U', one(T), B.data, A)

# Left division
A_ldiv_B!{T<:BlasFloat,S<:StridedMatrix}(A::UnitLowerTriangular{T,S}, B::StridedVecOrMat{T}) =
    LAPACK.trtrs!('L', 'N', 'U', A.data, B)

Ac_ldiv_B!{T<:BlasReal,S<:StridedMatrix}(A::UnitLowerTriangular{T,S}, B::StridedVecOrMat{T}) =
    LAPACK.trtrs!('L', 'T', 'U', A.data, B)

Ac_ldiv_B!{T<:BlasComplex,S<:StridedMatrix}(A::UnitLowerTriangular{T,S}, B::StridedVecOrMat{T}) =
    LAPACK.trtrs!('L', 'C', 'U', A.data, B)

# Right division
A_rdiv_B!{T<:BlasFloat,S<:StridedMatrix}(A::StridedVecOrMat{T}, B::UnitLowerTriangular{T,S}) =
    BLAS.trsm!('R', 'L', 'N', 'U', one(T), B.data, A)

A_rdiv_Bc!{T<:BlasReal,S<:StridedMatrix}(A::StridedMatrix{T}, B::UnitLowerTriangular{T,S}) =
    BLAS.trsm!('R', 'L', 'T', 'U', one(T), B.data, A)

A_rdiv_Bc!{T<:BlasComplex,S<:StridedMatrix}(A::StridedMatrix{T}, B::UnitLowerTriangular{T,S}) =
    BLAS.trsm!('R', 'L', 'C', 'U', one(T), B.data, A)

# Matrix inverse
inv{T<:BlasFloat,S<:StridedMatrix}(A::UnitLowerTriangular{T,S}) =
    UnitLowerTriangular{T,S}(LAPACK.trtri!('L', 'U', copy(A.data)))

# Error bounds for triangular solve
function errorbounds{T<:BlasFloat,S<:StridedMatrix}(A::UnitLowerTriangular{T,S},
                                                    X::StridedVecOrMat{T},
                                                    B::StridedVecOrMat{T})
    return LAPACK.trrfs!('L', 'N', 'U', A.data, B, X)
end

# Condition numbers
function cond{T<:BlasFloat,S}(A::UnitLowerTriangular{T,S}, p::Real=2)
    chksquare(A)
    if p == 1
        return inv(LAPACK.trcon!('O', 'L', 'U', A.data))
    elseif p == Inf
        return inv(LAPACK.trcon!('I', 'L', 'U', A.data))
    else #use fallback
        return cond(full(A), p)
    end
end

# Vector multiplication
A_mul_B!{T<:BlasFloat,S<:StridedMatrix}(A::UpperTriangular{T,S}, b::StridedVector{T}) =
    BLAS.trmv!('U', 'N', 'N', A.data, b)

# Matrix multiplication
A_mul_B!{T<:BlasFloat,S<:StridedMatrix}(A::UpperTriangular{T,S}, B::StridedMatrix{T}) =
    BLAS.trmm!('L', 'U', 'N', 'N', one(T), A.data, B)

A_mul_B!{T<:BlasFloat,S<:StridedMatrix}(A::StridedMatrix{T}, B::UpperTriangular{T,S}) =
    BLAS.trmm!('R', 'U', 'N', 'N', one(T), B.data, A)

Ac_mul_B!{T<:BlasComplex,S<:StridedMatrix}(A::UpperTriangular{T,S}, B::StridedMatrix{T}) =
    BLAS.trmm!('L', 'U', 'C', 'N', one(T), A.data, B)

Ac_mul_B!{T<:BlasReal,S<:StridedMatrix}(A::UpperTriangular{T,S}, B::StridedMatrix{T}) =
    BLAS.trmm!('L', 'U', 'T', 'N', one(T), A.data, B)

A_mul_Bc!{T<:BlasComplex,S<:StridedMatrix}(A::StridedMatrix{T}, B::UpperTriangular{T,S}) =
    BLAS.trmm!('R', 'U', 'C', 'N', one(T), B.data, A)

A_mul_Bc!{T<:BlasReal,S<:StridedMatrix}(A::StridedMatrix{T}, B::UpperTriangular{T,S}) =
    BLAS.trmm!('R', 'U', 'T', 'N', one(T), B.data, A)

# Left division
A_ldiv_B!{T<:BlasFloat,S<:StridedMatrix}(A::UpperTriangular{T,S}, B::StridedVecOrMat{T}) =
    LAPACK.trtrs!('U', 'N', 'N', A.data, B)

Ac_ldiv_B!{T<:BlasReal,S<:StridedMatrix}(A::UpperTriangular{T,S}, B::StridedVecOrMat{T}) =
    LAPACK.trtrs!('U', 'T', 'N', A.data, B)

Ac_ldiv_B!{T<:BlasComplex,S<:StridedMatrix}(A::UpperTriangular{T,S}, B::StridedVecOrMat{T}) =
    LAPACK.trtrs!('U', 'C', 'N', A.data, B)

# Right division
A_rdiv_B!{T<:BlasFloat,S<:StridedMatrix}(A::StridedVecOrMat{T}, B::UpperTriangular{T,S}) =
    BLAS.trsm!('R', 'U', 'N', 'N', one(T), B.data, A)

A_rdiv_Bc!{T<:BlasReal,S<:StridedMatrix}(A::StridedMatrix{T}, B::UpperTriangular{T,S}) =
    BLAS.trsm!('R', 'U', 'T', 'N', one(T), B.data, A)

A_rdiv_Bc!{T<:BlasComplex,S<:StridedMatrix}(A::StridedMatrix{T}, B::UpperTriangular{T,S}) =
    BLAS.trsm!('R', 'U', 'C', 'N', one(T), B.data, A)

# Matrix inverse
inv{T<:BlasFloat,S<:StridedMatrix}(A::UpperTriangular{T,S}) =
    UpperTriangular{T,S}(LAPACK.trtri!('U', 'N', copy(A.data)))

# Error bounds for triangular solve
function errorbounds{T<:BlasFloat,S<:StridedMatrix}(A::UpperTriangular{T,S},
                                                    X::StridedVecOrMat{T},
                                                    B::StridedVecOrMat{T})
    return LAPACK.trrfs!('U', 'N', 'N', A.data, B, X)
end

# Condition numbers
function cond{T<:BlasFloat,S}(A::UpperTriangular{T,S}, p::Real=2)
    chksquare(A)
    if p == 1
        return inv(LAPACK.trcon!('O', 'U', 'N', A.data))
    elseif p == Inf
        return inv(LAPACK.trcon!('I', 'U', 'N', A.data))
    else #use fallback
        return cond(full(A), p)
    end
end

# Vector multiplication
A_mul_B!{T<:BlasFloat,S<:StridedMatrix}(A::UnitUpperTriangular{T,S}, b::StridedVector{T}) =
    BLAS.trmv!('U', 'N', 'U', A.data, b)

# Matrix multiplication
A_mul_B!{T<:BlasFloat,S<:StridedMatrix}(A::UnitUpperTriangular{T,S}, B::StridedMatrix{T}) =
    BLAS.trmm!('L', 'U', 'N', 'U', one(T), A.data, B)

A_mul_B!{T<:BlasFloat,S<:StridedMatrix}(A::StridedMatrix{T}, B::UnitUpperTriangular{T,S}) =
    BLAS.trmm!('R', 'U', 'N', 'U', one(T), B.data, A)

Ac_mul_B!{T<:BlasComplex,S<:StridedMatrix}(A::UnitUpperTriangular{T,S}, B::StridedMatrix{T}) =
    BLAS.trmm!('L', 'U', 'C', 'U', one(T), A.data, B)

Ac_mul_B!{T<:BlasReal,S<:StridedMatrix}(A::UnitUpperTriangular{T,S}, B::StridedMatrix{T}) =
    BLAS.trmm!('L', 'U', 'T', 'U', one(T), A.data, B)

A_mul_Bc!{T<:BlasComplex,S<:StridedMatrix}(A::StridedMatrix{T}, B::UnitUpperTriangular{T,S}) =
    BLAS.trmm!('R', 'U', 'C', 'U', one(T), B.data, A)

A_mul_Bc!{T<:BlasReal,S<:StridedMatrix}(A::StridedMatrix{T}, B::UnitUpperTriangular{T,S}) =
    BLAS.trmm!('R', 'U', 'T', 'U', one(T), B.data, A)

# Left division
A_ldiv_B!{T<:BlasFloat,S<:StridedMatrix}(A::UnitUpperTriangular{T,S}, B::StridedVecOrMat{T}) =
    LAPACK.trtrs!('U', 'N', 'U', A.data, B)

Ac_ldiv_B!{T<:BlasReal,S<:StridedMatrix}(A::UnitUpperTriangular{T,S}, B::StridedVecOrMat{T}) =
    LAPACK.trtrs!('U', 'T', 'U', A.data, B)

Ac_ldiv_B!{T<:BlasComplex,S<:StridedMatrix}(A::UnitUpperTriangular{T,S}, B::StridedVecOrMat{T}) =
    LAPACK.trtrs!('U', 'C', 'U', A.data, B)

# Right division
A_rdiv_B!{T<:BlasFloat,S<:StridedMatrix}(A::StridedVecOrMat{T}, B::UnitUpperTriangular{T,S}) =
    BLAS.trsm!('R', 'U', 'N', 'U', one(T), B.data, A)

A_rdiv_Bc!{T<:BlasReal,S<:StridedMatrix}(A::StridedMatrix{T}, B::UnitUpperTriangular{T,S}) =
    BLAS.trsm!('R', 'U', 'T', 'U', one(T), B.data, A)

A_rdiv_Bc!{T<:BlasComplex,S<:StridedMatrix}(A::StridedMatrix{T}, B::UnitUpperTriangular{T,S}) =
    BLAS.trsm!('R', 'U', 'C', 'U', one(T), B.data, A)

# Matrix inverse
inv{T<:BlasFloat,S<:StridedMatrix}(A::UnitUpperTriangular{T,S}) =
    UnitUpperTriangular{T,S}(LAPACK.trtri!('U', 'U', copy(A.data)))

# Error bounds for triangular solve
function errorbounds{T<:BlasFloat,S<:StridedMatrix}(A::UnitUpperTriangular{T,S},
                                                    X::StridedVecOrMat{T},
                                                    B::StridedVecOrMat{T})
    return LAPACK.trrfs!('U', 'N', 'U', A.data, B, X)
end

# Condition numbers
function cond{T<:BlasFloat,S}(A::UnitUpperTriangular{T,S}, p::Real=2)
    chksquare(A)
    if p==1
        return inv(LAPACK.trcon!('O', 'U', 'U', A.data))
    elseif p==Inf
        return inv(LAPACK.trcon!('I', 'U', 'U', A.data))
    else #use fallback
        return cond(full(A), p)
    end
end

function errorbounds{T<:Union(BigFloat,Complex{BigFloat}),
                     S<:StridedMatrix}(A::AbstractTriangular{T,S},
                                       X::StridedVecOrMat{T},
                                       B::StridedVecOrMat{T})
    error("not implemented yet! Please submit a pull request.")
end

function errorbounds{TA<:Number,
                     S<:StridedMatrix,
                     TX<:Number,
                     TB<:Number}(A::AbstractTriangular{TA,S},
                                 X::StridedVecOrMat{TX},
                                 B::StridedVecOrMat{TB})
    TAXB = promote_type(TA, TB, TX, Float32)
    errorbounds(convert(AbstractMatrix{TAXB}, A),
                convert(AbstractArray{TAXB}, X),
                convert(AbstractArray{TAXB}, B))
end

# Eigensystems
## Notice that trecv works for quasi-triangular matrices and therefore
# the lower sub diagonal must be zeroed before calling the subroutine
eigvecs{T<:BlasFloat,S<:StridedMatrix}(A::UpperTriangular{T,S}) =
    LAPACK.trevc!('R', 'A', BlasInt[], triu!(A.data))

#TODO This mutates A, should be eigvecs?
function eigvecs{T<:BlasFloat,S<:StridedMatrix}(A::UnitUpperTriangular{T,S})
    for i = 1:size(A, 1)
        A.data[i,i] = 1
    end
    return LAPACK.trevc!('R', 'A', BlasInt[], triu!(A.data))
end

eigvecs{T<:BlasFloat,S<:StridedMatrix}(A::LowerTriangular{T,S}) =
    LAPACK.trevc!('L', 'A', BlasInt[], tril!(A.data)')

#TODO This mutates A, should be eigvecs?
function eigvecs{T<:BlasFloat,S<:StridedMatrix}(A::UnitLowerTriangular{T,S})
    for i = 1:size(A, 1)
        A.data[i,i] = 1
    end
    return LAPACK.trevc!('L', 'A', BlasInt[], tril!(A.data)')
end

####################
# Generic routines #
####################

(*)(A::UpperTriangular, x::Number) = UpperTriangular(A.data*x)
(*)(A::LowerTriangular, x::Number) = LowerTriangular(A.data*x)

function (*)(A::UnitUpperTriangular, x::Number)
    B = A.data*x
    for i = 1:size(A, 1)
        B[i,i] = x
    end
    return UpperTriangular(B)
end

function (*)(A::UnitLowerTriangular, x::Number)
    B = A.data*x
    for i = 1:size(A, 1)
        B[i,i] = x
    end
    return LowerTriangular(B)
end

(*)(x::Number, A::UpperTriangular) = UpperTriangular(x*A.data)
(*)(x::Number, A::LowerTriangular) = LowerTriangular(x*A.data)

function (*)(x::Number, A::UnitUpperTriangular)
    B = x*A.data
    for i = 1:size(A, 1)
        B[i,i] = x
    end
    return UpperTriangular(B)
end

function (*)(x::Number, A::UnitLowerTriangular)
    B = x*A.data
    for i = 1:size(A, 1)
        B[i,i] = x
    end
    return LowerTriangular(B)
end

(/)(A::UpperTriangular, x::Number) = UpperTriangular(A.data/x)
(/)(A::LowerTriangular, x::Number) = LowerTriangular(A.data/x)

function (/)(A::UnitUpperTriangular, x::Number)
    B = A.data*x
    invx = inv(x)
    for i = 1:size(A, 1)
        B[i,i] = x
    end
    return UpperTriangular(B)
end

function (/)(A::UnitLowerTriangular, x::Number)
    B = A.data*x
    invx = inv(x)
    for i = 1:size(A, 1)
        B[i,i] = x
    end
    return LowerTriangular(B)
end

(\)(x::Number, A::UpperTriangular) = UpperTriangular(x\A.data)
(\)(x::Number, A::LowerTriangular) = LowerTriangular(x\A.data)

function (\)(x::Number, A::UnitUpperTriangular)
    B = x\A.data
    invx = inv(x)
    for i = 1:size(A, 1)
        B[i,i] = invx
    end
    return UpperTriangular(B)
end

function (\)(x::Number, A::UnitLowerTriangular)
    B = x\A.data
    invx = inv(x)
    for i = 1:size(A, 1)
        B[i,i] = invx
    end
    return LowerTriangular(B)
end

## Generic triangular multiplication
function A_mul_B!(A::UpperTriangular, B::StridedVecOrMat)
    m, n = size(B, 1), size(B, 2)
    if m != size(A, 1)
        throw(DimensionMismatch("left and right hand side does not fit"))
    end
    for j = 1:n
        for i = 1:m
            Bij = A.data[i,i]*B[i,j]
            for k = i + 1:m
                Bij += A.data[i,k]*B[k,j]
            end
            B[i,j] = Bij
        end
    end
    return B
end

function A_mul_B!(A::UnitUpperTriangular, B::StridedVecOrMat)
    m, n = size(B, 1), size(B, 2)
    if m != size(A, 1)
        throw(DimensionMismatch("left and right hand side does not fit"))
    end
    for j = 1:n
        for i = 1:m
            Bij = B[i,j]
            for k = i + 1:m
                Bij += A.data[i,k]*B[k,j]
            end
            B[i,j] = Bij
        end
    end
    return B
end

function A_mul_B!(A::LowerTriangular, B::StridedVecOrMat)
    m, n = size(B, 1), size(B, 2)
    if m != size(A, 1)
        throw(DimensionMismatch("left and right hand side does not fit"))
    end
    for j = 1:n
        for i = m:-1:1
            Bij = A.data[i,i]*B[i,j]
            for k = 1:i - 1
                Bij += A.data[i,k]*B[k,j]
            end
            B[i,j] = Bij
        end
    end
    return B
end

function A_mul_B!(A::UnitLowerTriangular, B::StridedVecOrMat)
    m, n = size(B, 1), size(B, 2)
    if m != size(A, 1)
        throw(DimensionMismatch("left and right hand side does not fit"))
    end
    for j = 1:n
        for i = m:-1:1
            Bij = B[i,j]
            for k = 1:i - 1
                Bij += A.data[i,k]*B[k,j]
            end
            B[i,j] = Bij
        end
    end
    return B
end

function Ac_mul_B!(A::UpperTriangular, B::StridedVecOrMat)
    m, n = size(B, 1), size(B, 2)
    if m != size(A, 1)
        throw(DimensionMismatch("left and right hand side does not fit"))
    end
    for j = 1:n
        for i = m:-1:1
            Bij = A.data[i,i]*B[i,j]
            for k = 1:i - 1
                Bij += A.data[k,i]'B[k,j]
            end
            B[i,j] = Bij
        end
    end
    return B
end

function Ac_mul_B!(A::UnitUpperTriangular, B::StridedVecOrMat)
    m, n = size(B, 1), size(B, 2)
    if m != size(A, 1)
        throw(DimensionMismatch("left and right hand side does not fit"))
    end
    for j = 1:n
        for i = m:-1:1
            Bij = B[i,j]
            for k = 1:i - 1
                Bij += A.data[k,i]'B[k,j]
            end
            B[i,j] = Bij
        end
    end
    return B
end

function Ac_mul_B!(A::LowerTriangular, B::StridedVecOrMat)
    m, n = size(B, 1), size(B, 2)
    if m != size(A, 1)
        throw(DimensionMismatch("left and right hand side does not fit"))
    end
    for j = 1:n
        for i = 1:m
            Bij = A.data[i,i]*B[i,j]
            for k = i + 1:m
                Bij += A.data[k,i]'B[k,j]
            end
            B[i,j] = Bij
        end
    end
    return B
end

function Ac_mul_B!(A::UnitLowerTriangular, B::StridedVecOrMat)
    m, n = size(B, 1), size(B, 2)
    if m != size(A, 1)
        throw(DimensionMismatch("left and right hand side does not fit"))
    end
    for j = 1:n
        for i = 1:m
            Bij = B[i,j]
            for k = i + 1:m
                Bij += A.data[k,i]'B[k,j]
            end
            B[i,j] = Bij
        end
    end
    return B
end

function A_mul_B!(A::StridedMatrix, B::UpperTriangular)
    m, n = size(A)
    if size(B, 1) != n
        throw(DimensionMismatch("left and right hand side does not fit"))
    end
    for i = 1:m
        for j = n:-1:1
            Aij = A[i,j]*B[j,j]
            for k = 1:j - 1
                Aij += A[i,k]*B.data[k,j]
            end
            A[i,j] = Aij
        end
    end
    return A
end

function A_mul_B!(A::StridedMatrix, B::UnitUpperTriangular)
    m, n = size(A)
    if size(B, 1) != n
        throw(DimensionMismatch("left and right hand side does not fit"))
    end
    for i = 1:m
        for j = n:-1:1
            Aij = A[i,j]
            for k = 1:j - 1
                Aij += A[i,k]*B.data[k,j]
            end
            A[i,j] = Aij
        end
    end
    return A
end

function A_mul_B!(A::StridedMatrix, B::LowerTriangular)
    m, n = size(A)
    if size(B, 1) != n
        throw(DimensionMismatch("left and right hand side does not fit"))
    end
    for i = 1:m
        for j = 1:n
            Aij = A[i,j]*B[j,j]
            for k = j + 1:n
                Aij += A[i,k]*B.data[k,j]
            end
            A[i,j] = Aij
        end
    end
    return A
end

function A_mul_B!(A::StridedMatrix, B::UnitLowerTriangular)
    m, n = size(A)
    if size(B, 1) != n
        throw(DimensionMismatch("left and right hand side does not fit"))
    end
    for i = 1:m
        for j = 1:n
            Aij = A[i,j]
            for k = j + 1:n
                Aij += A[i,k]*B.data[k,j]
            end
            A[i,j] = Aij
        end
    end
    return A
end

function A_mul_Bc!(A::StridedMatrix, B::UpperTriangular)
    m, n = size(A)
    if size(B, 1) != n
        throw(DimensionMismatch("left and right hand side does not fit"))
    end
    for i = 1:m
        for j = 1:n
            Aij = A[i,j]*B[j,j]
            for k = j + 1:n
                Aij += A[i,k]*B.data[j,k]'
            end
            A[i,j] = Aij
        end
    end
    return A
end

function A_mul_Bc!(A::StridedMatrix, B::UnitUpperTriangular)
    m, n = size(A)
    if size(B, 1) != n
        throw(DimensionMismatch("left and right hand side does not fit"))
    end
    for i = 1:m
        for j = 1:n
            Aij = A[i,j]
            for k = j + 1:n
                Aij += A[i,k]*B.data[j,k]'
            end
            A[i,j] = Aij
        end
    end
    return A
end

function A_mul_Bc!(A::StridedMatrix, B::LowerTriangular)
    m, n = size(A)
    if size(B, 1) != n
        throw(DimensionMismatch("left and right hand side does not fit"))
    end
    for i = 1:m
        for j = n:-1:1
            Aij = A[i,j]*B[j,j]
            for k = 1:j - 1
                Aij += A[i,k]*B.data[j,k]'
            end
            A[i,j] = Aij
        end
    end
    return A
end

function A_mul_Bc!(A::StridedMatrix, B::UnitLowerTriangular)
    m, n = size(A)
    if size(B, 1) != n
        throw(DimensionMismatch("left and right hand side does not fit"))
    end
    for i = 1:m
        for j = n:-1:1
            Aij = A[i,j]
            for k = 1:j - 1
                Aij += A[i,k]*B.data[j,k]'
            end
            A[i,j] = Aij
        end
    end
    return A
end

#Generic solver using naive substitution
function naivesub!(A::UpperTriangular, b::AbstractVector, x::AbstractVector=b)
    N = size(A, 2)
    if N != length(b) != length(x)
        throw(DimensionMismatch())
    end
    for j = N:-1:1
        x[j] = b[j]
        for k = j+1:1:N
            x[j] -= A[j,k] * x[k]
        end
        x[j] = A[j,j]==0 ? throw(SingularException(j)) : A[j,j]\x[j]
    end
    return x
end

function naivesub!(A::UnitUpperTriangular, b::AbstractVector, x::AbstractVector=b)
    N = size(A, 2)
    if N != length(b) != length(x)
        throw(DimensionMismatch())
    end
    for j = N:-1:1
        x[j] = b[j]
        for k = j+1:1:N
            x[j] -= A[j,k] * x[k]
        end
    end
    return x
end

function naivesub!(A::LowerTriangular, b::AbstractVector, x::AbstractVector=b)
    N = size(A, 2)
    if N != length(b) != length(x)
        throw(DimensionMismatch())
    end
    for j = 1:N
        x[j] = b[j]
        for k = 1:j-1
            x[j] -= A[j,k] * x[k]
        end
        x[j] = A[j,j]==0 ? throw(SingularException(j)) : A[j,j]\x[j]
    end
    return x
end

function naivesub!(A::UnitLowerTriangular, b::AbstractVector, x::AbstractVector=b)
    N = size(A, 2)
    if N != length(b) != length(x)
        throw(DimensionMismatch())
    end
    for j = 1:N
        x[j] = b[j]
        for k = 1:j-1
            x[j] -= A[j,k] * x[k]
        end
    end
    return x
end

function A_rdiv_B!(A::StridedMatrix, B::UpperTriangular)
    m, n = size(A)
    if size(A, 1) != n
        throw(DimensionMismatch("left and right hand side does not fit"))
    end
    for i = 1:m
        for j = 1:n
            Aij = A[i,j]
            for k = 1:j - 1
                Aij -= A[i,k]*B.data[k,j]
            end
            A[i,j] = Aij/B[j,j]
        end
    end
    return A
end

function A_rdiv_B!(A::StridedMatrix, B::UnitUpperTriangular)
    m, n = size(A)
    if size(A, 1) != n
        throw(DimensionMismatch("left and right hand side does not fit"))
    end
    for i = 1:m
        for j = 1:n
            Aij = A[i,j]
            for k = 1:j - 1
                Aij -= A[i,k]*B.data[k,j]
            end
            A[i,j] = Aij
        end
    end
    return A
end

function A_rdiv_B!(A::StridedMatrix, B::LowerTriangular)
    m, n = size(A)
    if size(A, 1) != n
        throw(DimensionMismatch("left and right hand side does not fit"))
    end
    for i = 1:m
        for j = n:-1:1
            Aij = A[i,j]
            for k = j + 1:n
                Aij -= A[i,k]*B.data[k,j]
            end
            A[i,j] = Aij/B[j,j]
        end
    end
    return A
end

function A_rdiv_B!(A::StridedMatrix, B::UnitLowerTriangular)
    m, n = size(A)
    if size(A, 1) != n
        throw(DimensionMismatch("left and right hand side does not fit"))
    end
    for i = 1:m
        for j = n:-1:1
            Aij = A[i,j]
            for k = j + 1:n
                Aij -= A[i,k]*B.data[k,j]
            end
            A[i,j] = Aij
        end
    end
    return A
end

function A_rdiv_Bc!(A::StridedMatrix, B::UpperTriangular)
    m, n = size(A)
    if size(A, 1) != n
        throw(DimensionMismatch("left and right hand side does not fit"))
    end
    for i = 1:m
        for j = n:-1:1
            Aij = A[i,j]
            for k = j + 1:n
                Aij -= A[i,k]*B.data[j,k]'
            end
            A[i,j] = Aij/B[j,j]
        end
    end
    return A
end

function A_rdiv_Bc!(A::StridedMatrix, B::UnitUpperTriangular)
    m, n = size(A)
    if size(A, 1) != n
        throw(DimensionMismatch("left and right hand side does not fit"))
    end
    for i = 1:m
        for j = n:-1:1
            Aij = A[i,j]
            for k = j + 1:n
                Aij -= A[i,k]*B.data[j,k]'
            end
            A[i,j] = Aij
        end
    end
    return A
end

function A_rdiv_Bc!(A::StridedMatrix, B::LowerTriangular)
    m, n = size(A)
    if size(A, 1) != n
        throw(DimensionMismatch("left and right hand side does not fit"))
    end
    for i = 1:m
        for j = 1:n
            Aij = A[i,j]
            for k = 1:j - 1
                Aij -= A[i,k]*B.data[j,k]'
            end
            A[i,j] = Aij/B[j,j]
        end
    end
    return A
end

function A_rdiv_Bc!(A::StridedMatrix, B::UnitLowerTriangular)
    m, n = size(A)
    if size(A, 1) != n
        throw(DimensionMismatch("left and right hand side does not fit"))
    end
    for i = 1:m
        for j = 1:n
            Aij = A[i,j]
            for k = 1:j - 1
                Aij -= A[i,k]*B.data[j,k]'
            end
            A[i,j] = Aij
        end
    end
    return A
end

# Promotion
## Promotion methods in matmul don't apply to triangular multiplication since it is inplace.
# Hence we have to make very similar definitions, but without allocation of a result array.

# For multiplication and unit diagonal division the element type doesn't have to be stable
# under division whereas that is necessary in the general triangular solve problem.

## Some Triangular-Triangular cases. We might want to write taylored methods for these cases, but I'm not sure it is worth it.
(*)(A::Tridiagonal, B::UpperTriangular) = A_mul_B!(full(A), B)
(*)(A::Tridiagonal, B::UnitUpperTriangular) = A_mul_B!(full(A), B)
(*)(A::Tridiagonal, B::LowerTriangular) = A_mul_B!(full(A), B)
(*)(A::Tridiagonal, B::UnitLowerTriangular) = A_mul_B!(full(A), B)

(*)(A::AbstractTriangular, B::AbstractTriangular) = (*)(A, full(B))
Ac_mul_B(A::AbstractTriangular, B::AbstractTriangular) = Ac_mul_B(A, full(B))
At_mul_B(A::AbstractTriangular, B::AbstractTriangular) = At_mul_B(A, full(B))
A_mul_Bc(A::AbstractTriangular, B::AbstractTriangular) = A_mul_Bc(A, full(B))
A_mul_Bt(A::AbstractTriangular, B::AbstractTriangular) = A_mul_Bt(A, full(B))
Ac_mul_Bc(A::AbstractTriangular, B::AbstractTriangular) = Ac_mul_Bc(A, full(B))
At_mul_Bt(A::AbstractTriangular, B::AbstractTriangular) = At_mul_Bt(A, full(B))

(\)(A::AbstractTriangular, B::AbstractTriangular) = (\)(A, full(B))
Ac_ldiv_B(A::AbstractTriangular, B::AbstractTriangular) = Ac_ldiv_B(A, full(B))
At_ldiv_B(A::AbstractTriangular, B::AbstractTriangular) = At_ldiv_B(A, full(B))

A_mul_Bc(A::AbstractTriangular, B::AbstractTriangular) = A_mul_Bc(full(A), B)
A_mul_Bt(A::AbstractTriangular, B::AbstractTriangular) = A_mul_Bt(full(A), B)
Ac_mul_Bc(A::AbstractTriangular, B::AbstractTriangular) = Ac_mul_Bc(full(A), B)
At_mul_Bt(A::AbstractTriangular, B::AbstractTriangular) = At_mul_Bt(full(A), B)

(/)(A::AbstractTriangular, B::AbstractTriangular) = (/)(full(A), B)
A_rdiv_Bc(A::AbstractTriangular, B::AbstractTriangular) = A_rdiv_Bc(full(A), B)
A_rdiv_Bt(A::AbstractTriangular, B::AbstractTriangular) = A_rdiv_Bt(full(A), B)

## The general promotion methods
### Multiplication with triangle to the left and hence rhs cannot be transposed.
for (f, g) in ((:*, :A_mul_B!), (:Ac_mul_B, :Ac_mul_B!), (:At_mul_B, :At_mul_B!))
    @eval begin
        function ($f){TA,TB}(A::AbstractTriangular{TA}, B::StridedVecOrMat{TB})
            TAB = typeof(zero(TA)*zero(TB) + zero(TA)*zero(TB))
            ($g)(TA == TAB ? copy(A) : convert(AbstractArray{TAB}, A), TB == TAB ? copy(B) : convert(AbstractArray{TAB}, B))
        end
    end
end
### Left division with triangle to the left hence rhs cannot be transposed. No quotients.
for (f, g) in ((:\, :A_ldiv_B!), (:Ac_ldiv_B, :Ac_ldiv_B!), (:At_ldiv_B, :At_ldiv_B!))
    @eval begin
        function ($f){TA,TB,S}(A::UnitUpperTriangular{TA,S}, B::StridedVecOrMat{TB})
            TAB = typeof(zero(TA)*zero(TB) + zero(TA)*zero(TB))
            ($g)(TA == TAB ? copy(A) : convert(AbstractArray{TAB}, A), TB == TAB ? copy(B) : convert(AbstractArray{TAB}, B))
        end
    end
end
for (f, g) in ((:\, :A_ldiv_B!), (:Ac_ldiv_B, :Ac_ldiv_B!), (:At_ldiv_B, :At_ldiv_B!))
    @eval begin
        function ($f){TA,TB,S}(A::UnitLowerTriangular{TA,S}, B::StridedVecOrMat{TB})
            TAB = typeof(zero(TA)*zero(TB) + zero(TA)*zero(TB))
            ($g)(TA == TAB ? copy(A) : convert(AbstractArray{TAB}, A), TB == TAB ? copy(B) : convert(AbstractArray{TAB}, B))
        end
    end
end
### Left division with triangle to the left hence rhs cannot be transposed. Quotients.
for (f, g) in ((:\, :A_ldiv_B!), (:Ac_ldiv_B, :Ac_ldiv_B!), (:At_ldiv_B, :At_ldiv_B!))
    @eval begin
        function ($f){TA,TB,S}(A::UpperTriangular{TA,S}, B::StridedVecOrMat{TB})
            TAB = typeof((zero(TA)*zero(TB) + zero(TA)*zero(TB))/one(TA))
            ($g)(TA == TAB ? copy(A) : convert(AbstractArray{TAB}, A), TB == TAB ? copy(B) : convert(AbstractArray{TAB}, B))
        end
    end
end
for (f, g) in ((:\, :A_ldiv_B!), (:Ac_ldiv_B, :Ac_ldiv_B!), (:At_ldiv_B, :At_ldiv_B!))
    @eval begin
        function ($f){TA,TB,S}(A::LowerTriangular{TA,S}, B::StridedVecOrMat{TB})
            TAB = typeof((zero(TA)*zero(TB) + zero(TA)*zero(TB))/one(TA))
            ($g)(TA == TAB ? copy(A) : convert(AbstractArray{TAB}, A), TB == TAB ? copy(B) : convert(AbstractArray{TAB}, B))
        end
    end
end
### Multiplication with triangle to the rigth and hence lhs cannot be transposed.
for (f, g) in ((:*, :A_mul_B!), (:A_mul_Bc, :A_mul_Bc!), (:A_mul_Bt, :A_mul_Bt!))
    @eval begin
        function ($f){TA,TB}(A::StridedVecOrMat{TA}, B::AbstractTriangular{TB})
            TAB = typeof(zero(TA)*zero(TB) + zero(TA)*zero(TB))
            ($g)(TA == TAB ? copy(A) : convert(AbstractArray{TAB}, A), TB == TAB ? copy(B) : convert(AbstractArray{TAB}, B))
        end
    end
end
### Right division with triangle to the right hence lhs cannot be transposed. No quotients.
for (f, g) in ((:/, :A_rdiv_B!), (:A_rdiv_Bc, :A_rdiv_Bc!), (:A_rdiv_Bt, :A_rdiv_Bt!))
    @eval begin
        function ($f){TA,TB,S}(A::StridedVecOrMat{TA}, B::UnitUpperTriangular{TB,S})
            TAB = typeof(zero(TA)*zero(TB) + zero(TA)*zero(TB))
            ($g)(TA == TAB ? copy(A) : convert(AbstractArray{TAB}, A), TB == TAB ? copy(B) : convert(AbstractArray{TAB}, B))
        end
    end
end
for (f, g) in ((:/, :A_rdiv_B!), (:A_rdiv_Bc, :A_rdiv_Bc!), (:A_rdiv_Bt, :A_rdiv_Bt!))
    @eval begin
        function ($f){TA,TB,S}(A::StridedVecOrMat{TA}, B::UnitLowerTriangular{TB,S})
            TAB = typeof(zero(TA)*zero(TB) + zero(TA)*zero(TB))
            ($g)(TA == TAB ? copy(A) : convert(AbstractArray{TAB}, A), TB == TAB ? copy(B) : convert(AbstractArray{TAB}, B))
        end
    end
end
### Right division with triangle to the right hence lhs cannot be transposed. Quotients.
for (f, g) in ((:/, :A_rdiv_B!), (:A_rdiv_Bc, :A_rdiv_Bc!), (:A_rdiv_Bt, :A_rdiv_Bt!))
    @eval begin
        function ($f){TA,TB,S}(A::StridedVecOrMat{TA}, B::UpperTriangular{TB,S})
            TAB = typeof((zero(TA)*zero(TB) + zero(TA)*zero(TB))/one(TA))
            ($g)(TA == TAB ? copy(A) : convert(AbstractArray{TAB}, A), TB == TAB ? copy(B) : convert(AbstractArray{TAB}, B))
        end
    end
end
for (f, g) in ((:/, :A_rdiv_B!), (:A_rdiv_Bc, :A_rdiv_Bc!), (:A_rdiv_Bt, :A_rdiv_Bt!))
    @eval begin
        function ($f){TA,TB,S}(A::StridedVecOrMat{TA}, B::LowerTriangular{TB,S})
            TAB = typeof((zero(TA)*zero(TB) + zero(TA)*zero(TB))/one(TA))
            ($g)(TA == TAB ? copy(A) : convert(AbstractArray{TAB}, A), TB == TAB ? copy(B) : convert(AbstractArray{TAB}, B))
        end
    end
end

function sqrtm{T}(A::UpperTriangular{T})
    TT = typeof(sqrt(zero(T)))
    n = size(A, 1)
    R = zeros(TT, n, n)
    for j = 1:n
        R[j,j] = sqrt(A[j,j])
        for i = j-1:-1:1
            r = A[i,j]
            for k = i+1:j-1
                r -= R[i,k] * R[k,j]
            end
            if r != zero(T)
                R[i,j] = r / (R[i,i] + R[j,j])
            end
        end
    end
    return UpperTriangular(R)
end

function sqrtm{T}(A::UnitUpperTriangular{T})
    TT = typeof(sqrt(zero(T)))
    n = size(A, 1)
    R = zeros(TT, n, n)
    for j = 1:n
        R[j,j] = one(T)
        for i = j-1:-1:1
            r = A[i,j]
            for k = i+1:j-1
                r -= R[i,k] * R[k,j]
            end
            if r != zero(T)
                R[i,j] = r / (R[i,i] + R[j,j])
            end
        end
    end
    return UnitUpperTriangular(R)
end

sqrtm(A::LowerTriangular) = sqrtm(A.').'
sqrtm(A::UnitLowerTriangular) = sqrtm(A.').'

#Generic eigensystems
eigvals(A::AbstractTriangular) = diag(A)

function eigvecs{T}(A::AbstractTriangular{T})
    TT = promote_type(T, Float32)
    if !issubtype(TT, BlasFloat)
        throw(ArgumentError(
            "eigvecs type $(typeof(A)) not supported. Please submit a pull request."))
    end
    return eigvecs(convert(AbstractMatrix{TT}, A))
end

det{T}(A::UnitUpperTriangular{T}) = one(T) * one(T)
det{T}(A::UnitLowerTriangular{T}) = one(T) * one(T)
det{T}(A::UpperTriangular{T}) = prod(diag(A.data))
det{T}(A::LowerTriangular{T}) = prod(diag(A.data))

eigfact(A::AbstractTriangular) = Eigen(eigvals(A), eigvecs(A))

#Generic singular systems
svd(A::AbstractTriangular) = svd(full(A))
svdfact(A::AbstractTriangular) = svdfact(full(A))
svdfact!(A::AbstractTriangular) = svdfact!(full(A))
svdvals(A::AbstractTriangular) = svdvals(full(A))

factorize(A::AbstractTriangular) = A
