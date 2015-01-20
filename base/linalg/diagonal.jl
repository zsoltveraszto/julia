## Diagonal matrices

immutable Diagonal{T} <: AbstractMatrix{T}
    diag::Vector{T}
end
Diagonal(A::AbstractMatrix) = Diagonal(diag(A))

convert{T}(::Type{Diagonal{T}}, D::Diagonal{T}) = D
convert{T}(::Type{Diagonal{T}}, D::Diagonal) = Diagonal{T}(convert(Vector{T}, D.diag))
convert{T}(::Type{AbstractMatrix{T}}, D::Diagonal) = convert(Diagonal{T}, D)
convert{T}(::Type{UpperTriangular}, A::Diagonal{T}) = UpperTriangular(A)
convert{T}(::Type{LowerTriangular}, A::Diagonal{T}) = LowerTriangular(A)

function similar{T}(D::Diagonal, ::Type{T}, d::(Int,Int))
    if d[1] != d[2]
        throw(ArgumentError("Diagonal matrix must be square"))
    end
    return Diagonal{T}(Array(T,d[1]))
end

copy!(D1::Diagonal, D2::Diagonal) = (copy!(D1.diag, D2.diag); D1)

function size(D::Diagonal,d::Integer)
    if d < 1
        throw(ArgumentError("dimension must be ≥ 1, got $d"))
    elseif d <= 2
        return length(D.diag)
    else
        return 1
    end
end
size(D::Diagonal) = (length(D.diag), length(D.diag))

fill!(D::Diagonal, x) = (fill!(D.diag, x); D)

full(D::Diagonal) = diagm(D.diag)

#TODO: Getindex is broken
function getindex{T}(D::Diagonal{T}, i::Integer, j::Integer)
    if i == j
        return D.diag[i]
    end
    return zero(T)
end

function getindex(D::Diagonal, i::Integer)
    n = length(D.diag)
    id = div(i-1, n)
    if id + id * n == i-1
        return D.diag[id+1]
    end
    return zero(eltype(D.diag))
end

ishermitian{T<:Real}(D::Diagonal{T}) = true
ishermitian(D::Diagonal) = all(D.diag .== real(D.diag))
issym(D::Diagonal) = true
isposdef{T}(D::Diagonal{T}) = all(D.diag .> zero(T))

factorize(D::Diagonal) = D

#TODO: bounds checking
tril(D::Diagonal) = tril!(full(D), 0)
tril(D::Diagonal, i::Integer) = tril!(full(D), i)

tril!(D::Diagonal) = D
function tril!(D::Diagonal, i::Integer)
    if i == 0
        return D
    end
    return tril!(full(D), i)
end

triu(d::Diagonal) = triu!(full(D), 0)
triu(d::Diagonal, i::Integer) = triu!(full(D), i)

triu!(D::Diagonal) = D
function triu!(D::Diagonal, i::Integer)
    if i == 0
        return D
    end
    return triu!(full(D), i)
end

(==)(Da::Diagonal, Db::Diagonal) = Da.diag == Db.diag

(+)(Da::Diagonal, Db::Diagonal) = Diagonal(Da.diag + Db.diag)
(-)(Da::Diagonal, Db::Diagonal) = Diagonal(Da.diag - Db.diag)

(*){T<:Number}(x::T, D::Diagonal) = Diagonal(x * D.diag)
(*){T<:Number}(D::Diagonal, x::T) = Diagonal(D.diag * x)
(/){T<:Number}(D::Diagonal, x::T) = Diagonal(D.diag / x)
(*)(Da::Diagonal, Db::Diagonal) = Diagonal(Da.diag .* Db.diag)
(*)(D::Diagonal, V::Vector) = D.diag .* V
(*)(A::Matrix, D::Diagonal) = scale(A,D.diag)
(*)(D::Diagonal, A::Matrix) = scale(D.diag,A)

A_mul_B!(A::Diagonal,B::AbstractMatrix) = scale!(A.diag,B)
At_mul_B!(A::Diagonal,B::AbstractMatrix) = scale!(A.diag,B)
Ac_mul_B!(A::Diagonal,B::AbstractMatrix) = scale!(conj(A.diag),B)

(/)(Da::Diagonal, Db::Diagonal) = Diagonal(Da.diag ./ Db.diag )

function A_ldiv_B!{T}(D::Diagonal{T}, V::AbstractVector{T})
    if length(v) != length(D.diag)
        throw(DimensionMismatch())
    end
    for i=1:length(D.diag)
        d = D.diag[i]
        if d == zero(T)
            throw(SingularException(i))
        end
        V[i] *= inv(d)
    end
    return V
end

function A_ldiv_B!{T}(D::Diagonal{T}, V::AbstractMatrix{T})
    if size(V,1) != length(D.diag)
        throw(DimensionMismatch())
    end
    for i=1:length(D.diag)
        d = D.diag[i]
        if d == zero(T)
            throw(SingularException(i))
        end
        V[i,:] *= inv(d)
    end
    return V
end

conj(D::Diagonal) = Diagonal(conj(D.diag))
transpose(D::Diagonal) = D
ctranspose(D::Diagonal) = conj(D)

diag(D::Diagonal) = D.diag
trace(D::Diagonal) = sum(D.diag)
det(D::Diagonal) = prod(D.diag)

logdet{T<:Real}(D::Diagonal{T}) = sum(log(D.diag))

function logdet{T<:Complex}(D::Diagonal{T})
    #Make sure branch cut is correct
    x = sum(log(D.diag))
    if -pi < imag(x) < pi
        return x
    end
    return real(x) + (mod2pi(imag(x) + pi) - pi) * im
end

# identity matrices via eye(Diagonal{type},n)
eye{T}(::Type{Diagonal{T}}, n::Int) = Diagonal(ones(T,n))

expm(D::Diagonal) = Diagonal(exp(D.diag))
sqrtm(D::Diagonal) = Diagonal(sqrt(D.diag))

#Linear solver
function (\){TD<:Number,TA<:Number}(D::Diagonal{TD}, A::AbstractArray{TA,1})
    if size(A,2) == 1
        m,n = (size(A,1), 1)
    else
        m,n = size(A)
    end
    if m != length(D.diag)
        throw(DimensionMismatch())
    end
    if m == 0 || n == 0
        return A
    end
    C = Array(typeof(one(TD)/one(TA)),size(A))
    for j = 1:n
        for i = 1:m
            di = D.diag[i]
            if di == 0
                throw(SingularException(i))
            end
            C[i,j] = A[i,j] / di
        end
    end
    return C
end
(\)(Da::Diagonal, Db::Diagonal) = Diagonal(Db.diag ./ Da.diag)

function inv{T}(D::Diagonal{T})
    Di = similar(D.diag)
    for i = 1:length(D.diag)
        if D.diag[i] == zero(T)
            throw(SingularException(i))
        end
        Di[i] = inv(D.diag[i])
    end
    return Diagonal(Di)
end

function pinv{T}(D::Diagonal{T})
    Di = similar(D.diag)
    for i = 1:length(D.diag)
        if isfinite(inv(D.diag[i]))
            Di[i] = inv(D.diag[i])
        else
            Di[i] = zero(T)
        end
    end
    return Diagonal(Di)
end

function pinv{T}(D::Diagonal{T}, tol::Real)
    Di = similar(D.diag)
    if length(D.diag) != zero(T)
        maxabsD = maximum(abs(D.diag))
    end
    for i = 1:length(D.diag)
        if abs(D.diag[i]) > tol*maxabsD && isfinite(inv(D.diag[i]))
            Di[i] = inv(D.diag[i])
        else
            Di[i] = zero(T)
        end
    end
    return Diagonal(Di)
end

#Eigensystem
eigvals{T<:Number}(D::Diagonal{T}) = D.diag
eigvals(D::Diagonal) = [eigvals(x) for x in D.diag] #For block matrices, etc.
eigvecs(D::Diagonal) = eye(D)
eigfact(D::Diagonal) = Eigen(eigvals(D), eigvecs(D))

#Singular system
svdvals(D::Diagonal) = sort(D.diag, rev=true)

function svdfact(D::Diagonal, thin=true)
    S = abs(D.diag)
    piv = sortperm(S, rev=true)
    U = full(Diagonal(D.diag./S))
    Up = hcat([U[:,i] for i=1:length(D.diag)][piv]...)
    V = eye(D)
    Vp = hcat([V[:,i] for i=1:length(D.diag)][piv]...)
    return SVD(Up, S[piv], Vp')
end
