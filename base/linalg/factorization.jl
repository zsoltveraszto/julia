## Matrix factorizations and decompositions

eltype{T}(F::Factorization{T}) = T
transpose(F::Factorization) = error("transpose not implemented for $(typeof(F))")
ctranspose(F::Factorization) = error("ctranspose not implemented for $(typeof(F))")

macro assertposdef(A, info)
   :(($info)==0 ? $A : throw(PosDefException($info)))
end

macro assertnonsingular(A, info)
   :(($info)==0 ? $A : throw(SingularException($info)))
end

####################
# QR Factorization #
####################

immutable QR{T,S<:AbstractMatrix} <: Factorization{T}
    factors::S
    τ::Vector{T}

    function QR(factors::AbstractMatrix{T}, τ::Vector{T})
        return new(factors, τ)
    end
end
QR{T}(factors::AbstractMatrix{T}, τ::Vector{T}) =
    QR{T,typeof(factors)}(factors, τ)

# Note. For QRCompactWY factorization without pivoting,
# the WY representation based method introduced in LAPACK 3.4
immutable QRCompactWY{S,M<:AbstractMatrix} <: Factorization{S}
    factors::M
    T::M

    function QRCompactWY(factors::AbstractMatrix{S}, T::AbstractMatrix{S})
        new(factors, T)
    end
end
QRCompactWY{S}(factors::AbstractMatrix{S}, T::AbstractMatrix{S}) =
    QRCompactWY{S,typeof(factors)}(factors, T)

immutable QRPivoted{T,S<:AbstractMatrix} <: Factorization{T}
    factors::S
    τ::Vector{T}
    jpvt::Vector{BlasInt}

    function QRPivoted(factors::AbstractMatrix{T}, τ::Vector{T}, jpvt::Vector{BlasInt})
        return new(factors, τ, jpvt)
    end
end
QRPivoted{T}(factors::AbstractMatrix{T}, τ::Vector{T}, jpvt::Vector{BlasInt}) =
    QRPivoted{T,typeof(factors)}(factors, τ, jpvt)

function qrfact!{T<:BlasFloat}(A::StridedMatrix{T}; pivot::Bool=false)
    if pivot
        return QRPivoted(LAPACK.geqp3!(A)...)
    else
        return QRCompactWY(LAPACK.geqrt!(A, min(minimum(size(A)), 36))...)
    end
end

function qrfact!{T}(A::AbstractMatrix{T}, pivot::TrueOrFalse=Val{false})
    if pivot == Val{true}
        warn("pivoting only implemented for Float32, Float64, Complex64 and Complex128")
    end
    m, n = size(A)
    τ = zeros(T, min(m,n))
    @inbounds begin
        for k = 1:min(m-1 + !(T<:Real), n)
            τk = elementaryLeft!(A, k, k)
            τ[k] = τk
            for j = k+1:n
                vAj = A[k,j]
                for i = k+1:m
                    vAj += conj(A[i,k])*A[i,j]
                end
                vAj = conj(τk)*vAj
                A[k,j] -= vAj
                for i = k+1:m
                    A[i,j] -= A[i,k]*vAj
                end
            end
        end
    end
    return QR(A, τ)
end

function qrfact!{T<:BlasFloat}(A::StridedMatrix{T}, pivot::TrueOrFalse=Val{false})
    if pivot == Val{true}
        return QRPivoted(LAPACK.geqp3!(A)...)
    else
        return QRCompactWY(LAPACK.geqrt!(A, min(minimum(size(A)), 36))...)
    end
end

qrfact{T<:BlasFloat}(A::StridedMatrix{T}, pivot::TrueOrFalse=Val{false}) =
    qrfact!(copy(A), pivot)

copy_oftype{T}(A::StridedMatrix{T}, ::Type{T}) = copy(A)
copy_oftype{T,S}(A::StridedMatrix{T}, ::Type{S}) = convert(AbstractMatrix{S}, A)

qrfact{T}(A::StridedMatrix{T}, pivot::TrueOrFalse=Val{false}) =
    qrfact!(copy_oftype(A, typeof(one(T) / norm(one(T)))), pivot)

qrfact(x::Number) = qrfact(fill(x,1,1))

function qr(A::Union(Number, AbstractMatrix), pivot::TrueOrFalse=Val{false}; thin::Bool=true)
    return _qr(A, pivot, thin=thin)
end

function _qr(A::Union(Number, AbstractMatrix), ::Type{Val{false}}; thin::Bool=true)
    F = qrfact(A, Val{false})
    return full(F[:Q], thin=thin), F[:R]
end

function _qr(A::Union(Number, AbstractMatrix), ::Type{Val{true}}; thin::Bool=true)
    F = qrfact(A, Val{true})
    return full(F[:Q], thin=thin), F[:R], F[:p]
end

convert{T}(::Type{QR{T}},A::QR) =
    QR(convert(AbstractMatrix{T}, A.factors),
       convert(Vector{T}, A.τ))

convert{T}(::Type{Factorization{T}}, A::QR) =
    convert(QR{T}, A)

convert{T}(::Type{QRCompactWY{T}},A::QRCompactWY) =
    QRCompactWY(convert(AbstractMatrix{T}, A.factors),
                convert(AbstractMatrix{T}, A.T))

convert{T}(::Type{Factorization{T}}, A::QRCompactWY) =
    convert(QRCompactWY{T}, A)

convert{T}(::Type{QRPivoted{T}},A::QRPivoted) =
    QRPivoted(convert(AbstractMatrix{T}, A.factors),
              convert(Vector{T}, A.τ), A.jpvt)

convert{T}(::Type{Factorization{T}}, A::QRPivoted) =
    convert(QRPivoted{T}, A)

function getindex(A::QR, d::Symbol)
    m, n = size(A)
    if d == :R
        return triu!(A.factors[1:min(m,n), 1:n])
    elseif d == :Q
        return QRPackedQ(A.factors,A.τ)
    end
    throw(KeyError(d))
end

function getindex(A::QRCompactWY, d::Symbol)
    m, n = size(A)
    if d == :R
        return triu!(A.factors[1:min(m,n), 1:n])
    elseif d == :Q
        return QRCompactWYQ(A.factors,A.T)
    end
    throw(KeyError(d))
end

function getindex{T}(A::QRPivoted{T}, d::Symbol)
    m, n = size(A)
    if d == :R
        return triu!(A.factors[1:min(m,n), 1:n])
    elseif d == :Q
        return QRPackedQ(A.factors,A.τ)
    elseif d == :p
        return A.jpvt
    elseif d == :P
        p = A[:p]
        n = length(p)
        P = zeros(T, n, n)
        for i in 1:n
            P[p[i],i] = one(T)
        end
        return P
    end
    throw(KeyError(d))
end

# Type-stable interface to get Q
getq(A::QRCompactWY) = QRCompactWYQ(A.factors,A.T)
getq(A::QRPivoted) = QRPackedQ(A.factors,A.τ)

immutable QRPackedQ{T,S<:AbstractMatrix} <: AbstractMatrix{T}
    factors::S
    τ::Vector{T}

    function QRPackedQ(factors::AbstractMatrix{T}, τ::Vector{T})
        return new(factors, τ)
    end
end
QRPackedQ{T}(factors::AbstractMatrix{T}, τ::Vector{T}) =
    QRPackedQ{T,typeof(factors)}(factors, τ)

immutable QRPackedWYQ{S,M<:AbstractMatrix} <: AbstractMatrix{S}
    factors::M
    T::Matrix{S}

    function QRPackedWYQ(factors::AbstractMatrix{S}, T::Matrix{S})
        return new(factors, T)
    end
end
QRPackedWYQ{S}(factors::AbstractMatrix{S}, T::Matrix{S}) =
    QRPackedWYQ{S,typeof(factors)}(factors, T)

immutable QRCompactWYQ{S, M<:AbstractMatrix} <: AbstractMatrix{S}
    factors::M
    T::Matrix{S}

    function QRCompactWYQ(factors::AbstractMatrix{S}, T::Matrix{S})
        return new(factors, T)
    end
end
QRCompactWYQ{S}(factors::AbstractMatrix{S}, T::Matrix{S}) =
    QRCompactWYQ{S,typeof(factors)}(factors, T)

convert{T}(::Type{QRPackedQ{T}}, Q::QRPackedQ) =
    QRPackedQ(convert(AbstractMatrix{T}, Q.factors),
              convert(Vector{T}, Q.τ))

convert{T}(::Type{AbstractMatrix{T}}, Q::QRPackedQ) =
    convert(QRPackedQ{T}, Q)

convert{S}(::Type{QRCompactWYQ{S}}, Q::QRCompactWYQ) =
    QRCompactWYQ(convert(AbstractMatrix{S}, Q.factors),
                 convert(AbstractMatrix{S}, Q.T))

convert{S}(::Type{AbstractMatrix{S}}, Q::QRCompactWYQ) =
    convert(QRCompactWYQ{S}, Q)

size(A::Union(QR,QRCompactWY,QRPivoted)) = size(A.factors)
size(A::Union(QR,QRCompactWY,QRPivoted), dim::Integer) = size(A.factors, dim)

size(A::Union(QRPackedQ,QRCompactWYQ)) = (size(A, 1), size(A, 2))

function size(A::Union(QRPackedQ,QRCompactWYQ), dim::Integer)
    if dim < 0
        throw(BoundsError())
    elseif d <= 2
        return size(A.factors, 1)
    else
        return 1
    end
end

function full{T}(A::Union(QRPackedQ{T},QRCompactWYQ{T}); thin::Bool=true)
    if thin
        return A_mul_B!(A, eye(T, size(A.factors,1), minimum(size(A.factors))))
    else
        return A_mul_B!(A, eye(T, size(A.factors,1)))
    end
end

## Multiplication by Q
### QB
A_mul_B!{T<:BlasFloat}(A::QRCompactWYQ{T}, B::StridedVecOrMat{T}) =
    LAPACK.gemqrt!('L','N',A.factors,A.T,B)

A_mul_B!{T<:BlasFloat}(A::QRPackedQ{T}, B::StridedVecOrMat{T}) =
    LAPACK.ormqr!('L','N',A.factors,A.τ,B)

function A_mul_B!{T}(A::QRPackedQ{T}, B::AbstractVecOrMat{T})
    mA, nA = size(A.factors)
    mB, nB = size(B,1), size(B,2)
    if mA != mB
        throw(DimensionMismatch())
    end
    Afactors = A.factors
    @inbounds begin
        for k = min(mA,nA):-1:1
            for j = 1:nB
                vBj = B[k,j]
                for i = k+1:mB
                    vBj += conj(Afactors[i,k])*B[i,j]
                end
                vBj = A.τ[k]*vBj
                B[k,j] -= vBj
                for i = k+1:mB
                    B[i,j] -= Afactors[i,k]*vBj
                end
            end
        end
    end
    return B
end

function (*){TA,Tb}(A::Union(QRPackedQ{TA},QRCompactWYQ{TA}), b::StridedVector{Tb})
    TAb = promote_type(TA, Tb)
    Anew = convert(AbstractMatrix{TAb}, A)
    if size(A.factors, 1) == length
        if Tb == TAb
            bnew = A_nul_B!(Anew, copy(b))
        else
            bnew = A_mul_B!(Anew, convert(Vector{TAb}, b))
        end
    elseif size(A.factors, 2) == length(b)
        bnew = [b, zeros(TAb, size(A.factors,1) - length(b))]
    else
        throw(DimensionMismatch())
    end
    return A_mul_B!(Anew, bnew)
end

function (*){TA,TB}(A::Union(QRPackedQ{TA},QRCompactWYQ{TA}), B::StridedMatrix{TB})
    TAB = promote_type(TA, TB)
    Anew = convert(AbstractMatrix{TAB}, A)
    if size(A.factors,1) == size(B,1)
        if TB == TAB
            Bnew = copy(B)
        else
            Bnew = convert(AbstractMatrix{TAB}, B)
        end
    elseif size(A.factors,2) == size(B,1)
        Bnew = [B; zeros(TAB, size(A.factors,1) - size(B,1), size(B,2))]
    else
        throw(DimensionMismatch())
    end
    return A_mul_B!(Anew, Bnew)
end

### QcB
Ac_mul_B!{T<:BlasReal}(A::QRCompactWYQ{T}, B::StridedVecOrMat{T}) =
    LAPACK.gemqrt!('L','T',A.factors,A.T,B)

Ac_mul_B!{T<:BlasComplex}(A::QRCompactWYQ{T}, B::StridedVecOrMat{T}) =
    LAPACK.gemqrt!('L','C',A.factors,A.T,B)

Ac_mul_B!{T<:BlasReal}(A::QRPackedQ{T}, B::StridedVecOrMat{T}) =
    LAPACK.ormqr!('L','T',A.factors,A.τ,B)

Ac_mul_B!{T<:BlasComplex}(A::QRPackedQ{T}, B::StridedVecOrMat{T}) =
    LAPACK.ormqr!('L','C',A.factors,A.τ,B)

function Ac_mul_B!{T}(A::QRPackedQ{T}, B::AbstractVecOrMat{T})
    mA, nA = size(A.factors)
    mB, nB = size(B,1), size(B,2)
    if mA != mB
        throw(DimensionMismatch())
    end
    Afactors = A.factors
    @inbounds begin
        for k = 1:min(mA,nA)
            for j = 1:nB
                vBj = B[k,j]
                for i = k+1:mB
                    vBj += conj(Afactors[i,k])*B[i,j]
                end
                vBj = conj(A.τ[k])*vBj
                B[k,j] -= vBj
                for i = k+1:mB
                    B[i,j] -= Afactors[i,k]*vBj
                end
            end
        end
    end
    return B
end

function Ac_mul_B{TQ<:Number,TB<:Number,N}(Q::Union(QRPackedQ{TQ},QRCompactWYQ{TQ}),
                                           B::StridedArray{TB,N})
    TQB = promote_type(TQ, TB)
    if TB != TQB
        return Ac_mul_B!(convert(AbstractMatrix{TQB}, Q),
                         convert(AbstractArray{TQB,N}, B))
    else
        return Ac_mul_B!(convert(AbstractMatrix{TQB}, Q),
                         copy(B))
    end
end

### AQ

A_mul_B!{T<:BlasFloat}(A::StridedVecOrMat{T}, B::QRCompactWYQ{T}) =
    LAPACK.gemqrt!('R','N', B.factors, B.T, A)

A_mul_B!(A::StridedVecOrMat, B::QRPackedQ) =
    LAPACK.ormqr!('R', 'N', B.factors, B.τ, A)

function A_mul_B!{T}(A::StridedMatrix{T}, Q::QRPackedQ{T})
    mQ, nQ = size(Q.factors)
    mA, nA = size(A,1), size(A,2)
    if nA != mQ
        throw(DimensionMismatch())
    end
    Qfactors = Q.factors
    @inbounds begin
        for k = 1:min(mQ,nQ)
            for i = 1:mA
                vAi = A[i,k]
                for j = k+1:mQ
                    vAi += A[i,j]*Qfactors[j,k]
                end
                vAi = vAi*Q.τ[k]
                A[i,k] -= vAi
                for j = k+1:nA
                    A[i,j] -= vAi*conj(Qfactors[j,k])
                end
            end
        end
    end
    return A
end

function (*){TA,TQ,N}(A::StridedArray{TA,N},
                      Q::Union(QRPackedQ{TQ},QRCompactWYQ{TQ}))
    TAQ = promote_type(TA, TQ)
    if TA != TAQ
        return A_mul_B!(convert(AbstractArray{TAQ,N}, A),
                        convert(AbstractMatrix{TAQ}, Q))
    else
        return A_mul_B!(copy(A),
                        convert(AbstractMatrix{TAQ}, Q))
    end
end

### AQc
A_mul_Bc!{T<:BlasReal}(A::StridedVecOrMat{T}, B::QRCompactWYQ{T}) =
    LAPACK.gemqrt!('R','T',B.factors,B.T,A)

A_mul_Bc!{T<:BlasComplex}(A::StridedVecOrMat{T}, B::QRCompactWYQ{T}) =
    LAPACK.gemqrt!('R','C',B.factors,B.T,A)

A_mul_Bc!{T<:BlasReal}(A::StridedVecOrMat{T}, B::QRPackedQ{T}) =
    LAPACK.ormqr!('R','T',B.factors,B.τ,A)

A_mul_Bc!{T<:BlasComplex}(A::StridedVecOrMat{T}, B::QRPackedQ{T}) =
    LAPACK.ormqr!('R','C',B.factors,B.τ,A)

function A_mul_Bc!{T}(A::AbstractMatrix{T},Q::QRPackedQ{T})
    mQ, nQ = size(Q.factors)
    mA, nA = size(A,1), size(A,2)
    if nA != mQ
        throw(DimensionMismatch())
    end
    Qfactors = Q.factors
    @inbounds begin
        for k = min(mQ,nQ):-1:1
            for i = 1:mA
                vAi = A[i,k]
                for j = k+1:mQ
                    vAi += A[i,j]*Qfactors[j,k]
                end
                vAi = vAi*conj(Q.τ[k])
                A[i,k] -= vAi
                for j = k+1:nA
                    A[i,j] -= vAi*conj(Qfactors[j,k])
                end
            end
        end
    end
    return A
end

A_mul_Bc(A::AbstractTriangular, B::Union(QRCompactWYQ,QRPackedQ)) =
    A_mul_Bc(full(A), B)

function A_mul_Bc{TA,TB}(A::AbstractArray{TA},
                         B::Union(QRCompactWYQ{TB},QRPackedQ{TB}))
    TAB = promote_type(TA, TB)
    if size(A,2) == size(B.factors,1)
        if TA == TAB
            Anew = copy(A)
        else
            Anew = convert(AbstractMatrix{TAB}, A)
        end
    elseif size(A,2)==size(B.factors,2)
        Anew = [A zeros(TAB, size(A, 1), size(B.factors, 1) - size(B.factors, 2))]
    else
        throw(DimensionMismatch())
    end
    return A_mul_Bc!(Anew, convert(AbstractMatrix{TAB}, B))
end

function A_ldiv_B!{T<:BlasFloat}(A::QRCompactWY{T}, b::StridedVector{T})
    AcmB = Ac_mul_B!(A[:Q], b)
    A_ldiv_B!(UpperTriangular(A[:R]), sub(AcmB, 1:size(A, 2)))
    return b
end

function A_ldiv_B!{T<:BlasFloat}(A::QRCompactWY{T}, B::StridedMatrix{T})
    AcmB = Ac_mul_B!(A[:Q], B)
    A_ldiv_B!(UpperTriangular(A[:R]), sub(AcmB, 1:size(A, 2), 1:size(B, 2)))
    return B
end

# Julia implementation similarly to xgelsy
function A_ldiv_B!{T<:BlasFloat}(A::QRPivoted{T}, B::StridedMatrix{T}, rcond::Real)
    mA, nA = size(A.factors)
    nr = min(mA,nA)
    nrhs = size(B, 2)
    if nr == 0
        return (zeros(T, 0, nrhs), 0)
    end
    ar = abs(A.factors[1])
    if ar == 0
        return (zeros(T, nr, nrhs), 0)
    end
    rnk = 1
    xmin = ones(T, 1)
    xmax = ones(T, 1)
    tmin = tmax = ar
    while rnk < nr
        tmin, smin, cmin = LAPACK.laic1!(2, xmin, tmin,
                                         sub(A.factors, 1:rnk, rnk + 1),
                                         A.factors[rnk + 1, rnk + 1])
        tmax, smax, cmax = LAPACK.laic1!(1, xmax, tmax,
                                         sub(A.factors, 1:rnk, rnk + 1),
                                         A.factors[rnk + 1, rnk + 1])
        if tmax * rcond > tmin
            break
        end
        push!(xmin, cmin)
        push!(xmax, cmax)
        for i = 1:rnk
            xmin[i] *= smin
            xmax[i] = smax*xmin[i]
        end
        rnk += 1
        # if cond(r[1:rnk, 1:rnk])*rcond < 1 break end
    end
    C, τ = LAPACK.tzrzf!(A.factors[1:rnk,:])
    AcmB = Ac_mul_B!(getq(A), sub(B, 1:mA, 1:nrhs))
    A_ldiv_B!(UpperTriangular(C[1:rnk, 1:rnk]), sub(AcmB, 1:rnk, 1:nrhs))
    B[rnk+1:end,:] = zero(T)
    if iseltype(B, Complex)
        LAPACK.ormrz!('L', 'C', C, τ, sub(B, 1:nA, 1:nrhs))
    else
        LAPACK.ormrz!('L', 'T', C, τ, sub(B, 1:nA, 1:nrhs))
    end
    subB = sub(B, 1:nA, :)
    B[1:nA, :] = subB[invperm(A[:p]::Vector{BlasInt}), :]
    return (B, rnk)
end

A_ldiv_B!{T<:BlasFloat}(A::QRPivoted{T}, B::StridedVector{T}) =
    vec(A_ldiv_B!(A, reshape(B, length(B), 1)))

A_ldiv_B!{T<:BlasFloat}(A::QRPivoted{T}, B::StridedVecOrMat{T}) =
    A_ldiv_B!(A, B, maximum(size(A)) * eps(real(float(one(eltype(B))))))[1]

function A_ldiv_B!{T}(A::QR{T},B::StridedMatrix{T})
    m, n = size(A)
    minmn = min(m,n)
    mB, nB = size(B)
    Ac_mul_B!(A[:Q], sub(B,1:m,1:nB)) # Reconsider when arrayviews are merged.
    R = A[:R]
    @inbounds begin
        if n > m # minimum norm solution
            τ = zeros(T,m)
            for k = m:-1:1 # Trapezoid to triangular by elementary operation
                τ[k] = elementaryRightTrapezoid!(R,k)
                for i = 1:k-1
                    vRi = R[i,k]
                    for j = m+1:n
                        vRi += R[i,j]*R[k,j]
                    end
                    vRi *= τ[k]
                    R[i,k] -= vRi
                    for j = m+1:n
                        R[i,j] -= vRi*R[k,j]
                    end
                end
            end
        end
        # solve triangular system.
        # When array views are implemented, consider exporting to function.
        for k = 1:nB
            for i = minmn:-1:1
                for j = i+1:minmn
                    B[i,k] -= R[i,j]*B[j,k]
                end
                B[i,k] /= R[i,i]
            end
        end
        if n > m # Apply elementary transformation to solution
            B[m+1:mB,1:nB] = zero(T)
            for j = 1:nB
                for k = 1:m
                    vBj = B[k,j]
                    for i = m+1:n
                        vBj += B[i,j]*conj(R[k,i])
                    end
                    vBj *= τ[k]
                    B[k,j] -= vBj
                    for i = m+1:n
                        B[i,j] -= R[k,i]*vBj
                    end
                end
            end
        end
    end
    return B
end

A_ldiv_B!(A::QR, B::StridedVector) =
    A_ldiv_B!(A, reshape(B, length(B), 1))[:]

function A_ldiv_B!(A::QRPivoted, b::StridedVector)
    A_ldiv_B!(QR(A.factors,A.τ), b)
    b[1:size(A.factors, 2)] = sub(b, 1:size(A.factors, 2))[invperm(A.jpvt)]
    return b
end

function A_ldiv_B!(A::QRPivoted, B::StridedMatrix)
    A_ldiv_B!(QR(A.factors, A.τ), B)
    B[1:size(A.factors, 2),:] = sub(B, 1:size(A.factors, 2), :)[invperm(A.jpvt)]
    return B
end

function (\){TA,Tb}(A::Union(QR{TA},QRCompactWY{TA},QRPivoted{TA}),
                    b::StridedVector{Tb})
    S = promote_type(TA, Tb)
    m, n = size(A)
    if m != length(b)
        throw(DimensionMismatch(
            "left hand side has $m rows, but right hand side has length $(length(b))"))
    end
    if n > m
        X = A_ldiv_B!(convert(Factorization{S},A), [b, zeros(S,n-m)])
    elseif S == Tb
        X = A_ldiv_B!(convert(Factorization{S},A), copy(b))
    else
        X = A_ldiv_B!(convert(Factorization{S},A), convert(AbstractVector{S}, b))
    end
    if length(X) > n
        return X[1:n]
    end
    return X
end

function (\){TA,TB}(A::Union(QR{TA},QRCompactWY{TA},QRPivoted{TA}),
                    B::StridedMatrix{TB})
    S = promote_type(TA,TB)
    m, n = size(A)
    if m != size(B,1)
        throw(DimensionMismatch(
            "left hand side has $m rows, but right hand side has $(size(B,1)) rows"))
    end
    if n > m
        X = A_ldiv_B!(convert(Factorization{S}, A), [B; zeros(S, n-m, size(B,2))])
    elseif S == TB
        X = A_ldiv_B!(convert(Factorization{S}, A), copy(B))
    else
        X = A_ldiv_B!(convert(Factorization{S}, A), convert(AbstractMatrix{S}, B))
    end
    if size(X, 1) > n
        return X[1:n, :]
    end
    return X
end

##TODO:  Add methods for rank(A::QRP{T}) and adjust the (\) method accordingly
##       Add rcond methods for Cholesky, LU, QR and QRP types
## Lower priority: Add LQ, QL and RQ factorizations

# FIXME! Should add balancing option through xgebal
immutable Hessenberg{T,S<:AbstractMatrix} <: Factorization{T}
    factors::S
    τ::Vector{T}

    function Hessenberg(factors::AbstractMatrix{T}, τ::Vector{T})
        return new(factors, τ)
    end
end
Hessenberg{T}(factors::AbstractMatrix{T}, τ::Vector{T}) =
    Hessenberg{T,typeof(factors)}(factors, τ)

Hessenberg(A::StridedMatrix) = Hessenberg(LAPACK.gehrd!(A)...)

hessfact!{T<:BlasFloat}(A::StridedMatrix{T}) = Hessenberg(A)
hessfact{T<:BlasFloat}(A::StridedMatrix{T}) = hessfact!(copy(A))

function hessfact{T}(A::StridedMatrix{T})
    S = promote_type(Float32, typeof(one(T) / norm(one(T))))
    if S != T
        return hessfact!(convert(AbstractMatrix{S}, A))
    else
        return hessfact!(copy(A))
    end
end

immutable HessenbergQ{T,S<:AbstractMatrix} <: AbstractMatrix{T}
    factors::S
    τ::Vector{T}

    HessenbergQ(factors::AbstractMatrix{T}, τ::Vector{T}) = begin
        new(factors, τ)
    end
end
HessenbergQ{T}(factors::AbstractMatrix{T}, τ::Vector{T}) =
    HessenbergQ{T,typeof(factors)}(factors, τ)

HessenbergQ(A::Hessenberg) = HessenbergQ(A.factors, A.τ)

size(A::HessenbergQ, args...) = size(A.factors, args...)

function getindex(A::Hessenberg, d::Symbol)
    if d == :Q
        return HessenbergQ(A)
    elseif d == :H
        return triu(A.factors, -1)
    end
    throw(KeyError(d))
end

full(A::HessenbergQ) =
    LAPACK.orghr!(1, size(A.factors, 1), copy(A.factors), A.τ)

# Also printing of QRQs
function getindex(A::Union(QRPackedQ,QRCompactWYQ,HessenbergQ),
                  i::Integer, j::Integer)
    x = zeros(eltype(A), size(A, 1))
    x[i] = 1
    y = zeros(eltype(A), size(A, 2))
    y[j] = 1
    return dot(x, A*y)
end

#######################
# Eigendecompositions #
#######################

# Eigenvalues
immutable Eigen{T,V,S<:AbstractMatrix,U<:AbstractVector} <: Factorization{T}
    values::U
    vectors::S

    function Eigen(values::AbstractVector{V}, vectors::AbstractMatrix{T})
        return new(values, vectors)
    end
end
Eigen{T,V}(values::AbstractVector{V}, vectors::AbstractMatrix{T}) =
    Eigen{T,V,typeof(vectors),typeof(values)}(values, vectors)

# Generalized eigenvalue problem.
immutable GeneralizedEigen{T,V,S<:AbstractMatrix,U<:AbstractVector} <: Factorization{T}
    values::U
    vectors::S

    function GeneralizedEigen(values::AbstractVector{V}, vectors::AbstractMatrix{T})
        return new(values, vectors)
    end
end
GeneralizedEigen{T,V}(values::AbstractVector{V}, vectors::AbstractMatrix{T}) =
    GeneralizedEigen{T,V,typeof(vectors),typeof(values)}(values, vectors)

function getindex(A::Union(Eigen,GeneralizedEigen), d::Symbol)
    if d == :values
        return A.values
    elseif d == :vectors
        return A.vectors
    end
    throw(KeyError(d))
end

isposdef(A::Union(Eigen,GeneralizedEigen)) = all(A.values .> 0)

function eigfact!{T<:BlasReal}(A::StridedMatrix{T};
                               permute::Bool=true, scale::Bool=true)
    n = size(A, 2)
    if n==0
        return Eigen(zeros(T, 0), zeros(T, 0, 0))
    end
    if issym(A)
        return eigfact!(Symmetric(A))
    end
    if permute && scale
        PS = 'B'
    elseif permute
        PS = 'P'
    elseif scale
        PS = 'S'
    else
        PS = 'N'
    end
    A, WR, WI, VL, VR, _ = LAPACK.geevx!(PS, 'N', 'V', 'N', A)
    if all(WI .== 0.)
        return Eigen(WR, VR)
    end
    evec = zeros(Complex{T}, n, n)
    j = 1
    while j <= n
        if WI[j] == 0.0
            evec[:,j] = VR[:,j]
        else
            evec[:,j]   = VR[:,j] + im*VR[:,j+1]
            evec[:,j+1] = VR[:,j] - im*VR[:,j+1]
            j += 1
        end
        j += 1
    end
    return Eigen(complex(WR, WI), evec)
end

function eigfact!{T<:BlasComplex}(A::StridedMatrix{T};
                                  permute::Bool=true, scale::Bool=true)
    n = size(A, 2)
    if n == 0
        return Eigen(zeros(T, 0), zeros(T, 0, 0))
    end
    if ishermitian(A)
        return eigfact!(Hermitian(A))
    end
    if permute && scale
        PS = 'B'
    elseif permute
        PS = 'P'
    elseif scale
        PS = 'S'
    else
        PS = 'N'
    end
    A, WR, WI, VL, VR, _ = LAPACK.geevx!(PS, 'N', 'V', 'N', A)
    return Eigen(WR, VL)
end

function eigfact{T}(A::StridedMatrix{T};
                    permute::Bool=true, scale::Bool=true)
    S = promote_type(Float32, typeof(one(T) / norm(one(T))))
    if S != T
        return eigfact!(convert(AbstractMatrix{S}, A),
                        permute=permute, scale=scale)
    else
        return eigfact!(copy(A),
                        permute=permute, scale=scale)
    end
end

eigfact(x::Number) = Eigen([x], fill(one(x), 1, 1))

# function eig(A::Union(Number, AbstractMatrix); permute::Bool=true, scale::Bool=true)
#     F = eigfact(A, permute=permute, scale=scale)
#     F[:values], F[:vectors]
# end
function eig(A::Union(Number, AbstractMatrix), args...; kwargs...)
    F = eigfact(A, args...; kwargs...)
    return (F[:values], F[:vectors])
end

#Calculates eigenvectors
eigvecs(A::Union(Number, AbstractMatrix), args...; kwargs...) =
    eigfact(A, args...; kwargs...)[:vectors]

function eigvals!{T<:BlasReal}(A::StridedMatrix{T};
                               permute::Bool=true, scale::Bool=true)
    if issym(A)
        return eigvals!(Symmetric(A))
    end
    if permute && scale
        PS = 'B'
    elseif permute
        PS = 'P'
    elseif scale
        PS = 'S'
    else
        PS = 'N'
    end
    _, valsre, valsim, _ = LAPACK.geevx!(PS, 'N', 'N', 'N', A)
    if all(valsim .== 0)
        return valsre
    end
    return complex(valsre, valsim)
end

function eigvals!{T<:BlasComplex}(A::StridedMatrix{T};
                                  permute::Bool=true, scale::Bool=true)
    if ishermitian(A)
        return eigvals(Hermitian(A))
    end
    if permute && scale
        PS = 'B'
    elseif permute
        PS = 'P'
    elseif scale
        PS = 'S'
    else
        PS = 'N'
    end
    return LAPACK.geevx!(PS, 'N', 'N', 'N', A)[2]
end

function eigvals{T}(A::StridedMatrix{T};
                    permute::Bool=true, scale::Bool=true)
    S = promote_type(Float32, typeof(one(T) / norm(one(T))))
    if S != T
        return eigvals!(convert(AbstractMatrix{S}, A),
                        permute=permute, scale=scale)
    else
        return eigvals!(copy(A),
                        permute=permute, scale=scale)
    end
end

function eigvals{T<:Number}(x::T; kwargs...)
    val = convert(promote_type(Float32, typeof(one(T) / norm(one(T)))), x)
    if imag(val) == 0
        return [real(val)]
    end
    return [val]
end

#Computes maximum and minimum eigenvalue
function eigmax(A::Union(Number, StridedMatrix);
                permute::Bool=true, scale::Bool=true)
    v = eigvals(A, permute=permute, scale=scale)
    if iseltype(v, Complex)
        error("DomainError: complex eigenvalues cannot be ordered")
    end
    return maximum(v)
end

function eigmin(A::Union(Number, StridedMatrix);
                permute::Bool=true, scale::Bool=true)
    v = eigvals(A, permute=permute, scale=scale)
    if iseltype(v, Complex)
        error("DomainError: complex eigenvalues cannot be ordered")
    end
    return minimum(v)
end

inv(A::Eigen) = A.vectors / Diagonal(A.values) * A.vectors'
det(A::Eigen) = prod(A.values)

# Generalized eigenproblem
function eigfact!{T<:BlasReal}(A::StridedMatrix{T}, B::StridedMatrix{T})
    if issym(A) && isposdef(B)
        return eigfact!(Symmetric(A), Symmetric(B))
    end
    n = size(A, 1)
    alphar, alphai, beta, _, vr = LAPACK.ggev!('N', 'V', A, B)
    if all(alphai .== 0)
        return GeneralizedEigen(alphar ./ beta, vr)
    end
    vecs = zeros(Complex{T}, n, n)
    j = 1
    while j <= n
        if alphai[j] == 0.0
            vecs[:,j] = vr[:,j]
        else
            vecs[:,j  ] = vr[:,j] + im*vr[:,j+1]
            vecs[:,j+1] = vr[:,j] - im*vr[:,j+1]
            j += 1
        end
        j += 1
    end
    return GeneralizedEigen((complex(alphar, alphai) ./ beta), vecs)
end

function eigfact!{T<:BlasComplex}(A::StridedMatrix{T}, B::StridedMatrix{T})
    if ishermitian(A) && isposdef(B)
        return eigfact!(Hermitian(A), Hermitian(B))
    end
    alpha, beta, _, vr = LAPACK.ggev!('N', 'V', A, B)
    return GeneralizedEigen((alpha ./ beta), vr)
end

function eigfact{TA,TB}(A::AbstractMatrix{TA}, B::AbstractMatrix{TB})
    S = promote_type(Float32, typeof(one(TA) / norm(one(TA))), TB)
    if S != TA && S != TB
        return eigfact!(convert(AbstractMatrix{S}, A),
                        convert(AbstractMatrix{S}, B))
    elseif S != TA
        return eigfact!(convert(AbstractMatrix{S}, A),
                        copy(B))
    elseif S != TB
        return eigfact!(copy(A),
                        convert(AbstractMatrix{S}, B))
    else
        return eigfact!(copy(A), copy(B))
    end
end

function eigvals!{T<:BlasReal}(A::StridedMatrix{T}, B::StridedMatrix{T})
    if issym(A) && isposdef(B)
        return eigvals!(Symmetric(A), Symmetric(B))
    end
    alphar, alphai, beta, vl, vr = LAPACK.ggev!('N', 'N', A, B)
    if all(alphai .== 0)
        return alphar ./ beta
    end
    return complex(alphar, alphai) ./ beta
end

function eigvals!{T<:BlasComplex}(A::StridedMatrix{T}, B::StridedMatrix{T})
    if ishermitian(A) && isposdef(B)
        return eigvals!(Hermitian(A), Hermitian(B))
    end
    alpha, beta, vl, vr = LAPACK.ggev!('N', 'N', A, B)
    return alpha ./ beta
end

function eigvals{TA,TB}(A::AbstractMatrix{TA}, B::AbstractMatrix{TB})
    S = promote_type(Float32, typeof(one(TA) / norm(one(TA))), TB)
    if S != TA && S != TB
        return eigvals!(convert(AbstractMatrix{S}, A),
                        convert(AbstractMatrix{S}, B))
    elseif S != TA
        return eigvals!(convert(AbstractMatrix{S}, A),
                        copy(B))
    elseif S != TB
        return eigvals!(copy(A),
                        convert(AbstractMatrix{S}, B))
    else
        return eigvals!(copy(A), copy(B))
    end
end

# SVD
immutable SVD{T<:BlasFloat,Tr,M<:AbstractArray} <: Factorization{T}
    U::M
    S::Vector{Tr}
    Vt::M

    function SVD(U::AbstractArray{T}, S::Vector{Tr}, Vt::AbstractArray{T})
        return new(U, S, Vt)
    end
end
SVD{T<:BlasFloat,Tr}(U::AbstractArray{T}, S::Vector{Tr}, Vt::AbstractArray{T}) =
    SVD{T,Tr,typeof(U)}(U, S, Vt)

function svdfact!{T<:BlasFloat}(A::StridedMatrix{T}; thin::Bool=true)
    m, n = size(A)
    if m == 0 || n == 0
        if thin
            u = eye(T, m, n)
        else
            u = eye(T, m, m)
        end
        s = real(zeros(T, 0))
        vt = eye(T, n, n)
    else
        if thin
            u,s,vt = LAPACK.gesdd!('S', A)
        else
            u,s,vt = LAPACK.gesdd!('A', A)
        end
    end
    return SVD(u,s,vt)
end

svdfact{T<:BlasFloat}(A::StridedMatrix{T}; thin::Bool=true) =
    svdfact!(copy(A), thin=thin)

function svdfact{T}(A::StridedVecOrMat{T}; thin::Bool=true)
    S = promote_type(Float32, typeof(one(T) / norm(one(T))))
    if S != T
        return svdfact!(convert(AbstractMatrix{S}, A), thin=thin)
    else
        return svdfact!(copy(A), thin=thin)
    end
end

function svdfact(x::Number; thin::Bool=true)
    if x == 0
        U = fill(one(x), 1, 1)
    else
        U = fill(x / abs(x), 1, 1)
    end
    return SVD(U, [abs(x)], fill(one(x), 1, 1))
end

svdfact(x::Integer; thin::Bool=true) = svdfact(float(x), thin=thin)

function svd(A::Union(Number, AbstractArray); thin::Bool=true)
    F = svdfact(A, thin=thin)
    return (F.U, F.S, F.Vt')
end

function getindex(F::SVD, d::Symbol)
    if d == :U
        return F.U
    elseif d == :S
        return F.S
    elseif d == :Vt
        return F.Vt
    elseif d == :V
        return F.Vt'
    end
    throw(KeyError(d))
end

function svdvals!{T<:BlasFloat}(A::StridedMatrix{T})
    if any([size(A)...] .== 0)
        return zeros(T, 0)
    end
    return LAPACK.gesdd!('N', A)[2]
end

svdvals{T<:BlasFloat}(A::StridedMatrix{T}) = svdvals!(copy(A))

function svdvals{T}(A::StridedMatrix{T})
    S = promote_type(Float32, typeof(one(T) / norm(one(T))))
    if S != T
        return svdvals!(convert(AbstractMatrix{S}, A))
    else
        return svdvals!(copy(A))
    end
end

svdvals(x::Number) = [abs(x)]

# SVD least squares
function (\){T<:BlasFloat}(A::SVD{T}, B::StridedVecOrMat{T})
    n = length(A.S)
    Sinv = zeros(T, n)
    k = length(find(A.S .> eps(real(float(one(T)))) * maximum(A.S)))
    Sinv[1:k] = one(T) ./ A.S[1:k]
    return A.Vt[1:k,:]' * (Sinv[1:k] .* (A.U[:,1:k]' * B))
end

# Generalized svd
immutable GeneralizedSVD{T,S} <: Factorization{T}
    U::S
    V::S
    Q::S
    a::Vector
    b::Vector
    k::Int
    l::Int
    R::S

    function GeneralizedSVD(U::AbstractMatrix{T},
                            V::AbstractMatrix{T},
                            Q::AbstractMatrix{T},
                            a::Vector,
                            b::Vector,
                            k::Int,
                            l::Int,
                            R::AbstractMatrix{T})
        return new(U, V, Q, a, b, k, l, R)
    end
end

function GeneralizedSVD{T}(U::AbstractMatrix{T},
                           V::AbstractMatrix{T},
                           Q::AbstractMatrix{T},
                           a::Vector,
                           b::Vector,
                           k::Int,
                           l::Int,
                           R::AbstractMatrix{T})
    return GeneralizedSVD{T,typeof(U)}(U, V, Q, a, b, k, l, R)
end

function svdfact!{T<:BlasFloat}(A::StridedMatrix{T}, B::StridedMatrix{T})
    U, V, Q, a, b, k, l, R = LAPACK.ggsvd!('U', 'V', 'Q', A, B)
    return GeneralizedSVD(U, V, Q, a, b, int(k), int(l), R)
end

svdfact{T<:BlasFloat}(A::StridedMatrix{T}, B::StridedMatrix{T}) =
    svdfact!(copy(A), copy(B))

function svdfact{TA,TB}(A::StridedMatrix{TA}, B::StridedMatrix{TB})
    S = promote_type(Float32, typeof(one(TA) / norm(one(TA))), TB)
    if S != TA && S != TB
        return svdfact!(convert(AbstractMatrix{S}, A),
                        convert(AbstractMatrix{S}, B))
    elseif S != TA
        return svdfact!(convert(AbstractMatrix{S},A),
                        copy(B))
    elseif S != TB
        return svdfact!(copy(A),
                        convert(AbstractMatrix{S}, B))
    else
        return svdfact!(copy(A), copy(B))
    end
end

function svd(A::AbstractMatrix, B::AbstractMatrix)
    F = svdfact(A, B)
    return (F[:U], F[:V], F[:Q], F[:D1], F[:D2], F[:R0])
end

function getindex{T}(obj::GeneralizedSVD{T}, d::Symbol)
    if d == :U
        return obj.U
    elseif d == :V
        return obj.V
    elseif d == :Q
        return obj.Q
    elseif (d == :alpha || d == :a)
        return obj.a
    elseif (d == :beta || d == :b)
        return obj.b
    elseif (d == :vals || d == :S)
        return obj.a[1:obj.k + obj.l] ./ obj.b[1:obj.k + obj.l]
    elseif d == :D1
        m = size(obj.U, 1)
        if m - obj.k - obj.l >= 0
            return [eye(T, obj.k) zeros(T, obj.k, obj.l);
                    zeros(T, obj.l, obj.k) diagm(obj.a[obj.k + 1:obj.k + obj.l]);
                    zeros(T, m - obj.k - obj.l, obj.k + obj.l)]
        else
            return [eye(T, m, obj.k)
                    [zeros(T, obj.k, m - obj.k); diagm(obj.a[obj.k + 1:m])]
                    zeros(T, m, obj.k + obj.l - m)]
        end
    elseif d == :D2
        m = size(obj.U, 1)
        p = size(obj.V, 1)
        if m - obj.k - obj.l >= 0
            return [zeros(T, obj.l, obj.k) diagm(obj.b[obj.k + 1:obj.k + obj.l]);
                    zeros(T, p - obj.l, obj.k + obj.l)]
        else
            return [zeros(T, p, obj.k)
                    [diagm(obj.b[obj.k + 1:m]);
                     zeros(T, obj.k + p - m, m - obj.k)]
                    [zeros(T, m - obj.k, obj.k + obj.l - m);
                     eye(T, obj.k + p - m, obj.k + obj.l - m)]]
        end
    elseif d == :R
        return obj.R
    elseif d == :R0
        n = size(obj.Q, 1)
        return [zeros(T, obj.k + obj.l, n - obj.k - obj.l) obj.R]
    end
    throw(KeyError(d))
end

function svdvals!{T<:BlasFloat}(A::StridedMatrix{T}, B::StridedMatrix{T})
    _, _, _, a, b, k, l, _ = LAPACK.ggsvd!('N', 'N', 'N', A, B)
    return a[1:k + l] ./ b[1:k + l]
end

svdvals{T<:BlasFloat}(A::StridedMatrix{T},B::StridedMatrix{T}) =
    svdvals!(copy(A), copy(B))

function svdvals{TA,TB}(A::StridedMatrix{TA}, B::StridedMatrix{TB})
    S = promote_type(Float32, typeof(one(T) / norm(one(TA))), TB)
    if S != TA && S != TB
        return svdvals!(convert(AbstractMatrix{S}, A),
                        convert(AbstractMatrix{S}, B))
    elseif S != TA
        return svdvals!(convert(AbstractMatrix{S}, A),
                        copy(B))
    elseif S != TB
        return svdvals!(copy(A),
                        convert(AbstractMatrix{S}, B))
    else
        return svdvals!(copy(A), copy(B))
    end
end

immutable Schur{Ty<:BlasFloat, S<:AbstractMatrix} <: Factorization{Ty}
    T::S
    Z::S
    values::Vector

    function Schur(T::AbstractMatrix{Ty}, Z::AbstractMatrix{Ty}, values::Vector)
        return new(T, Z, values)
    end
end

Schur{Ty}(T::AbstractMatrix{Ty}, Z::AbstractMatrix{Ty}, values::Vector) =
    Schur{Ty, typeof(T)}(T, Z, values)

schurfact!{T<:BlasFloat}(A::StridedMatrix{T}) =
    Schur(LinAlg.LAPACK.gees!('V', A)...)

schurfact{T<:BlasFloat}(A::StridedMatrix{T}) =
    schurfact!(copy(A))

function schurfact{T}(A::StridedMatrix{T})
    S = promote_type(Float32, typeof(one(T) / norm(one(T))))
    if S != T
        return schurfact!(convert(AbstractMatrix{S}, A))
    else
        return schurfact!(copy(A))
    end
end

function getindex(F::Schur, d::Symbol)
    if (d == :T || d == :Schur)
        return F.T
    elseif (d == :Z || d == :vectors)
        return F.Z
    elseif d == :values
        return F.values
    end
    throw(KeyError(d))
end

function schur(A::StridedMatrix)
    SchurF = schurfact(A)
    return (SchurF[:T], SchurF[:Z], SchurF[:values])
end

ordschur!{Ty<:BlasFloat}(Q::StridedMatrix{Ty}, T::StridedMatrix{Ty}, select::Array{Int}) =
    Schur(LinAlg.LAPACK.trsen!(select, T , Q)...)

function ordschur!{Ty<:BlasFloat}(schur::Schur{Ty}, select::Array{Int})
    res = ordschur!(schur.Z, schur.T, select)
    schur[:values][:] = res[:values]
    return res
end

ordschur{Ty<:BlasFloat}(Q::StridedMatrix{Ty}, T::StridedMatrix{Ty}, select::Array{Int}) =
    ordschur!(copy(Q), copy(T), select)

ordschur{Ty<:BlasFloat}(schur::Schur{Ty}, select::Array{Int}) =
    ordschur(schur.Z, schur.T, select)

immutable GeneralizedSchur{Ty<:BlasFloat, M<:AbstractMatrix} <: Factorization{Ty}
    S::M
    T::M
    alpha::Vector
    beta::Vector{Ty}
    Q::M
    Z::M

    function GeneralizedSchur(S::AbstractMatrix{Ty},
                              T::AbstractMatrix{Ty},
                              alpha::Vector,
                              beta::Vector{Ty},
                              Q::AbstractMatrix{Ty},
                              Z::AbstractMatrix{Ty})
        return new(S, T, alpha, beta, Q, Z)
    end
end
function GeneralizedSchur{Ty}(S::AbstractMatrix{Ty},
                              T::AbstractMatrix{Ty},
                              alpha::Vector,
                              beta::Vector{Ty},
                              Q::AbstractMatrix{Ty},
                              Z::AbstractMatrix{Ty})
    return GeneralizedSchur{Ty,typeof(S)}(S, T, alpha, beta, Q, Z)
end

schurfact!{T<:BlasFloat}(A::StridedMatrix{T}, B::StridedMatrix{T}) =
    GeneralizedSchur(LinAlg.LAPACK.gges!('V', 'V', A, B)...)

schurfact{T<:BlasFloat}(A::StridedMatrix{T},B::StridedMatrix{T}) =
    schurfact!(copy(A), copy(B))

function schurfact{TA,TB}(A::StridedMatrix{TA}, B::StridedMatrix{TB})
    S = promote_type(Float32, typeof(one(TA) / norm(one(TA))), TB)
    if S != TA && S != TB
        return schurfact!(convert(AbstractMatrix{S}, A),
                          convert(AbstractMatrix{S}, B))
    elseif S != TA
        return schurfact!(convert(AbstractMatrix{S}, A),
                          copy(B))
    elseif S != TB
        return schurfact!(copy(A),
                          convert(AbstractMatrix{S}, B))
    else
        return schurfact!(copy(A), copy(B))
    end
end

function ordschur!{Ty<:BlasFloat}(S::StridedMatrix{Ty},
                                  T::StridedMatrix{Ty},
                                  Q::StridedMatrix{Ty},
                                  Z::StridedMatrix{Ty},
                                  select::Array{Int})
    return GeneralizedSchur(
        LinAlg.LAPACK.tgsen!(select, S, T, Q, Z)...)
end

function ordschur{Ty<:BlasFloat}(S::StridedMatrix{Ty},
                                 T::StridedMatrix{Ty},
                                 Q::StridedMatrix{Ty},
                                 Z::StridedMatrix{Ty},
                                 select::Array{Int})
    return ordschur!(copy(S), copy(T), copy(Q), copy(Z), select)
end

function ordschur!{Ty<:BlasFloat}(gschur::GeneralizedSchur{Ty},
                                  select::Array{Int})
    res = ordschur!(gschur.S, gschur.T, gschur.Q, gschur.Z, select)
    gschur[:alpha][:] = res[:alpha]
    gschur[:beta][:] = res[:beta]
    return res
end

ordschur{Ty<:BlasFloat}(gschur::GeneralizedSchur{Ty}, select::Array{Int}) =
    ordschur(gschur.S, gschur.T, gschur.Q, gschur.Z, select)

function getindex(F::GeneralizedSchur, d::Symbol)
    if d == :S
        return F.S
    elseif d == :T
        return F.T
    elseif d == :alpha
        return F.alpha
    elseif d == :beta
        return F.beta
    elseif d == :values
        return F.alpha ./ F.beta
    elseif (d == :Q || d == :left)
        return F.Q
    elseif (d == :Z || d == :right)
        return F.Z
    else
        throw(KeyError(d))
    end
end

function schur(A::StridedMatrix, B::StridedMatrix)
    SchurF = schurfact(A, B)
    return (SchurF[:S], SchurF[:T], SchurF[:Q], SchurF[:Z])
end

### General promotion rules ###
convert{T}(::Type{Factorization{T}}, F::Factorization{T}) = F
inv{T}(F::Factorization{T}) = A_ldiv_B!(F, eye(T, size(F,1)))

function (\){TF<:Number,TB<:Number,N}(F::Factorization{TF},
                                      B::AbstractArray{TB,N})
    TFB = typeof(one(TF) / one(TB))
    fact = convert(Factorization{TFB}, F)
    if TB == TFB
        return A_ldiv_B!(fact, copy(B))
    else
        return A_ldiv_B!(fact, convert(AbstractArray{TFB,N}, B))
    end
end

function Ac_ldiv_B{TF<:Number,TB<:Number,N}(F::Factorization{TF},
                                            B::AbstractArray{TB,N})
    TFB = typeof(one(TF) / one(TB))
    fact = convert(Factorization{TFB}, F)
    if TB == TFB
        return Ac_ldiv_B!(fact, copy(B))
    else
        return Ac_ldiv_B!(fact, convert(AbstractArray{TFB,N}, B))
    end
end

function At_ldiv_B{TF<:Number,TB<:Number,N}(F::Factorization{TF},
                                            B::AbstractArray{TB,N})
    TFB = typeof(one(TF) / one(TB))
    fact = convert(Factorization{TFB}, F)
    if TB == TFB
        return At_ldiv_B!(fact, copy(B))
    else
        return At_ldiv_B!(fact, convert(AbstractArray{TFB,N}, B))
    end
end
