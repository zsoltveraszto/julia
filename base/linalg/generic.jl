## linalg.jl: Some generic Linear Algebra definitions

scale(X::AbstractArray, s::Number) = scale!(copy(X), s)
scale(s::Number, X::AbstractArray) = scale!(copy(X), s)

function scale{R<:Real,S<:Complex}(X::AbstractArray{R}, s::S)
    Y = Array(promote_type(R,S), size(X))
    copy!(Y, X)
    return scale!(Y, s)
end

scale{R<:Real}(s::Complex, X::AbstractArray{R}) = scale(X, s)

# For better performance when input and output are the same array
# See https://github.com/JuliaLang/julia/issues/8415#issuecomment-56608729
function generic_scale!(X::AbstractArray, s::Number)
    for i = 1:length(X)
        @inbounds X[i] *= s
    end
    return X
end

function generic_scale!(C::AbstractArray, X::AbstractArray, s::Number)
    if length(C) != length(X)
        throw(DimensionMismatch(
            string("first array argument must be the same length ",
                   "as the second array argument")))
    end
    for i = 1:length(X)
        @inbounds C[i] = X[i]*s
    end
    return C
end

scale!(C::AbstractArray, s::Number, X::AbstractArray) = generic_scale!(C, X, s)
scale!(C::AbstractArray, X::AbstractArray, s::Number) = generic_scale!(C, X, s)
scale!(X::AbstractArray, s::Number) = generic_scale!(X, s)
scale!(s::Number, X::AbstractArray) = generic_scale!(X, s)

cross(A::AbstractVector, B::AbstractVector) =
    [A[2]*B[3]-A[3]*B[2], A[3]*B[1]-A[1]*B[3], A[1]*B[2]-A[2]*B[1]]

triu(M::AbstractMatrix) = triu!(copy(M))
tril(M::AbstractMatrix) = tril!(copy(M))
triu!(M::AbstractMatrix) = triu!(M,0)
tril!(M::AbstractMatrix) = tril!(M,0)

diff(a::AbstractMatrix) = diff(a, 1)
diff(a::AbstractVector) = [a[i+1] - a[i] for i=1:length(a)-1]

function diff(A::AbstractMatrix, dim::Integer)
    if dim == 1
        return [A[i+1,j] - A[i,j] for i=1:size(A,1)-1, j=1:size(A,2)]
    elseif dim == 2
        return [A[i,j+1] - A[i,j] for i=1:size(A,1), j=1:size(A,2)-1]
    else
        throw(ArgumentError("dimension dim must be 1 or 2, got $dim"))
    end
end


gradient(F::AbstractVector) = gradient(F, [1:length(F)])
gradient(F::AbstractVector, h::Real) = gradient(F, [h*(1:length(F))])

diag(A::AbstractVector) =
    error("use diagm instead of diag to construct a diagonal matrix")

#diagm{T}(v::AbstractVecOrMat{T})

# special cases of vecnorm; note that they don't need to handle isempty(x)
function generic_vecnormMinusInf(x)
    s = start(x)
    (v, s) = next(x, s)
    minabs = abs(v)
    while !done(x, s)
        (v, s) = next(x, s)
        minabs = Base.scalarmin(minabs, abs(v))
    end
    return float(minabs)
end

function generic_vecnormInf(x)
    s = start(x)
    (v, s) = next(x, s)
    maxabs = abs(v)
    while !done(x, s)
        (v, s) = next(x, s)
        maxabs = Base.scalarmax(maxabs, abs(v))
    end
    return float(maxabs)
end

function generic_vecnorm1(x)
    s = start(x)
    (v, s) = next(x, s)
    av = float(abs(v))
    T = typeof(av)
    sum::promote_type(Float64, T) = av
    while !done(x, s)
        (v, s) = next(x, s)
        sum += abs(v)
    end
    return convert(T, sum)
end

function generic_vecnorm2(x)
    maxabs = vecnormInf(x)
    if maxabs == 0
        return maxabs
    end
    s = start(x)
    (v, s) = next(x, s)
    T = typeof(maxabs)
    scale::promote_type(Float64, T) = 1 / maxabs
    y = abs(v)*scale
    sum::promote_type(Float64, T) = y * y
    while !done(x, s)
        (v, s) = next(x, s)
        y = abs(v)*scale
        sum += y*y
    end
    return convert(T, maxabs * sqrt(sum))
end

function generic_vecnormp(x, p)
    if p > 1 || p < 0 # need to rescale to avoid overflow/underflow
        maxabs = vecnormInf(x)
        if maxabs == 0
            return maxabs
        end
        s = start(x)
        (v, s) = next(x, s)
        T = typeof(maxabs)
        spp::promote_type(Float64, T) = p
        scale::promote_type(Float64, T) = 1/maxabs
        ssum::promote_type(Float64, T) = (abs(v) * scale)^spp
        while !done(x, s)
            (v, s) = next(x, s)
            ssum += (abs(v) * scale)^spp
        end
        return convert(T, maxabs * ssum^inv(spp))
    else
        # 0 < p < 1, no need for rescaling (but technically not a true norm)
        s = start(x)
        (v, s) = next(x, s)
        av = float(abs(v))
        T = typeof(av)
        pp::promote_type(Float64, T) = p
        sum::promote_type(Float64, T) = av^pp
        while !done(x, s)
            (v, s) = next(x, s)
            sum += abs(v)^pp
        end
        return convert(T, sum^inv(pp))
    end
end

vecnormMinusInf(x) = generic_vecnormMinusInf(x)
vecnormInf(x) = generic_vecnormInf(x)
vecnorm1(x) = generic_vecnorm1(x)
vecnorm2(x) = generic_vecnorm2(x)
vecnormp(x, p) = generic_vecnormp(x, p)

function vecnorm(itr, p::Real=2)
    if isempty(itr)
        return float(real(zero(eltype(itr))))
    elseif p == 2
        return vecnorm2(itr)
    elseif p == 1
        return vecnorm1(itr)
    elseif p == Inf
        return vecnormInf(itr)
    elseif p == 0
        return convert(typeof(float(real(zero(eltype(itr))))), countnz(itr))
    elseif p == -Inf
        return vecnormMinusInf(itr)
    else
        return vecnormp(itr,p)
    end
end

function vecnorm(x::Number, p::Real=2)
    if p == 0 && x == 0
        return real(zero(x))
    elseif p == 0
        return real(one(x))
    else
        return abs(x)
    end
end

norm(x::AbstractVector, p::Real=2) = vecnorm(x, p)

function norm1{T}(A::AbstractMatrix{T})
    m, n = size(A)
    Tnorm = typeof(float(real(zero(T))))
    Tsum = promote_type(Float64,Tnorm)
    nrm::Tsum = 0
    @inbounds begin
        for j = 1:n
            nrmj::Tsum = 0
            for i = 1:m
                nrmj += abs(A[i,j])
            end
            nrm = max(nrm,nrmj)
        end
    end
    return convert(Tnorm, nrm)
end

function norm2{T}(A::AbstractMatrix{T})
    m,n = size(A)
    if m == 1 || n == 1
        return vecnorm2(A)
    end
    Tnorm = typeof(float(real(zero(T))))
    if m == 0 || n == 0
        return zero(Tnorm)
    end
    return convert(Tnorm, svdvals(A)[1])
end

function normInf{T}(A::AbstractMatrix{T})
    m,n = size(A)
    Tnorm = typeof(float(real(zero(T))))
    Tsum = promote_type(Float64,Tnorm)
    nrm::Tsum = 0
    @inbounds begin
        for i = 1:m
            nrmi::Tsum = 0
            for j = 1:n
                nrmi += abs(A[i,j])
            end
            nrm = max(nrm,nrmi)
        end
    end
    return convert(Tnorm, nrm)
end

function norm{T}(A::AbstractMatrix{T}, p::Real=2)
    if p == 2
        return norm2(A)
    elseif p == 1
        return return norm1(A)
    elseif p == Inf
        return normInf(A)
    end
    throw(ArgumentError("invalid p-norm p=$p. Valid: 1, 2, Inf"))
end

function norm(x::Number, p::Real=2)
    if p == 0
        return convert(typeof(real(x)), ifelse(x != 0, 1, 0))
    else
        return abs(x)
    end
end

function rank(A::AbstractMatrix)
    m,n = size(A)
    if m == 0 || n == 0
        return 0
    end
    sv = svdvals(A)
    return sum(sv .> maximum(size(A)) * eps(sv[1]))
end
rank(A::AbstractMatrix, tol::Real) = sum(svdvals(A) .> tol)
function rank(x::Number)
    if x == 0
        return 0
    else
        return 1
    end
end

function trace(A::AbstractMatrix)
    chksquare(A)
    return sum(diag(A))
end
trace(x::Number) = x

#kron(a::AbstractVector, b::AbstractVector)
#kron{T,S}(a::AbstractMatrix{T}, b::AbstractMatrix{S})

#det(a::AbstractMatrix)

function inv{T}(A::AbstractMatrix{T})
    S = typeof(zero(T) / one(T))
    return A_ldiv_B!(factorize(convert(AbstractMatrix{S}, A)),
                     eye(S, chksquare(A)))
end
inv(a::StridedMatrix) = error("argument must be a square matrix")

function (\){TA,TB}(A::AbstractMatrix{TA}, B::AbstractVecOrMat{TB})
    TC = typeof(one(TA) / one(TB))
    if size(A,1) != size(B,1)
        throw(DimensionMismatch(
            string("LHS and RHS should have the same number of rows. ",
                    "LHS has $(size(A,1)) rows, but RHS has $(size(B,1)) rows.")))
    end
    if TA == TC && TB == TC
        return A_ldiv_B!(factorize(copy(A)), copy(B))
    elseif TA == TC
        return A_ldiv_B!(factorize(copy(A)),
                         convert(AbstractArray{TC}, B))
    elseif TB == TC
        return A_ldiv_B!(factorize(convert(AbstractMatrix{TC}, A)),
                         copy(B))
    else
        return A_ldiv_B!(factorize(convert(AbstractMatrix{TC}, A)),
                         convert(AbstractArray{TC}, B))
    end
end
(\)(a::AbstractVector, b::AbstractArray) = reshape(a, length(a), 1) \ b
(/)(A::AbstractVecOrMat, B::AbstractVecOrMat) = (B' \ A')'

# Should be added at some point when the old elementwise
# version has been deprecated long enough
# \(A::StridedMatrix,x::Number) = inv(A)*x
# /(x::Number,A::StridedMatrix) = x*inv(A)

function cond(x::Number)
    if x == 0
        return Inf
    else
        return 1.0
    end
end
cond(x::Number, p) = cond(x)

#Skeel condition numbers
condskeel(A::AbstractMatrix, p::Real=Inf) =
    norm(abs(inv(A)) * abs(A), p)

condskeel{T<:Integer}(A::AbstractMatrix{T}, p::Real=Inf) =
    norm(abs(inv(float(A))) * abs(A), p)

condskeel(A::AbstractMatrix, x::AbstractVector, p::Real=Inf) =
    norm(abs(inv(A)) * abs(A) * abs(x), p)

condskeel{T<:Integer}(A::AbstractMatrix{T}, x::AbstractVector, p::Real=Inf) =
    norm(abs(inv(float(A))) * abs(A) * abs(x), p)

function issym(A::AbstractMatrix)
    m, n = size(A)
    if m != n
        return false
    end
    for i = 1:(n-1), j = (i+1):n
        if A[i,j] != transpose(A[j,i])
            return false
        end
    end
    return true
end

issym(x::Number) = true

function ishermitian(A::AbstractMatrix)
    m, n = size(A)
    if m != n
        return false
    end
    for i = 1:n, j = i:n
        if A[i,j] != ctranspose(A[j,i])
            return false
        end
    end
    return true
end

ishermitian(x::Number) = (x == conj(x))

function istriu(A::AbstractMatrix)
    m, n = size(A)
    for j = 1:min(n,m-1), i = j+1:m
        if A[i,j] != 0
            return false
        end
    end
    return true
end

function istril(A::AbstractMatrix)
    m, n = size(A)
    for j = 2:n, i = 1:min(j-1,m)
        if A[i,j] != 0
            return false
        end
    end
    return true
end

istriu(x::Number) = true
istril(x::Number) = true

linreg{T<:Number}(X::StridedVecOrMat{T}, y::Vector{T}) =
    hcat(ones(T, size(X,1)), X) \ y

# weighted least squares
function linreg(x::AbstractVector,
                y::AbstractVector,
                w::AbstractVector)
    sw = sqrt(w)
    return hcat(sw, sw .* x) \ (sw.*y)
end

# multiply by diagonal matrix as vector
#diagmm!(C::AbstractMatrix, A::AbstractMatrix, b::AbstractVector)

#diagmm!(C::AbstractMatrix, b::AbstractVector, A::AbstractMatrix)

scale!(A::AbstractMatrix, b::AbstractVector) = scale!(A,A,b)
scale!(b::AbstractVector, A::AbstractMatrix) = scale!(A,b,A)

#diagmm(A::AbstractMatrix, b::AbstractVector)
#diagmm(b::AbstractVector, A::AbstractMatrix)

#^(A::AbstractMatrix, p::Number)

#findmax(a::AbstractArray)
#findmin(a::AbstractArray)

function peakflops(n::Integer=2000; parallel::Bool=false)
    a = rand(100,100)
    t = @elapsed a*a
    a = rand(n,n)
    t = @elapsed a*a
    if parallel
        return sum(pmap(peakflops, [ n for i in 1:nworkers()]))
    else
        return (2*float64(n)^3/t)
    end
end

# BLAS-like in-place y=alpha*x+y function (see also the version in blas.jl
#                                          for BlasFloat Arrays)
function axpy!(alpha, X::AbstractArray, Y::AbstractArray)
    n = length(X)
    if n != length(Y)
        throw(DimensionMismatch())
    end
    for i = 1:n
        @inbounds Y[i] += alpha * X[i]
    end
    return Y
end

function axpy!{Ti<:Integer,Tj<:Integer}(alpha,
                                        X::AbstractArray,
                                        rx::AbstractArray{Ti},
                                        Y::AbstractArray,
                                        ry::AbstractArray{Tj})
    if length(X) != length(Y)
        throw(DimensionMismatch())
    end
    if (minimum(rx) < 1 || maximum(rx) > length(X) ||
        minimum(ry) < 1 || maximum(ry) > length(Y) ||
        length(rx) != length(ry))
        throw(BoundsError())
    end
    for i = 1:length(rx)
        @inbounds Y[ry[i]] += alpha * X[rx[i]]
    end
    return Y
end

# Elementary reflection similar to LAPACK.
# The reflector is not Hermitian but ensures that
# tridiagonalization of Hermitian matrices become real. See lawn72
function elementaryLeft!(A::AbstractMatrix, row::Integer, col::Integer)
    m, n = size(A)
    if !(1 <= row <= m)
        throw(BoundsError(
            "row cannot be less than one or larger than $(size(A,1))"))
    end
    if !(1 <= col <= n)
        throw(BoundsError(
            "col cannot be less than one or larger than $(size(A,2))"))
    end
    @inbounds begin
        ξ1 = A[row,col]
        normu = abs2(ξ1)
        for i = row+1:m
            normu += abs2(A[i,col])
        end
        normu = sqrt(normu)
        ν = copysign(normu,real(ξ1))
        A[row,col] += ν
        ξ1 += ν
        A[row,col] = -ν
        for i = row+1:m
            A[i,col] /= ξ1
        end
    end
    return ξ1/ν
end

function elementaryRight!(A::AbstractMatrix, row::Integer, col::Integer)
    m, n = size(A)
    if !(1 <= row <= m)
        throw(BoundsError(
            "row cannot be less than one or larger than $(size(A,1))"))
    end
    if !(1 <= col <= n)
        throw(BoundsError(
            "col cannot be less than one or larger than $(size(A,2))"))
    end
    if col > row
        throw(BoundsError(
            "col cannot be larger than row"))
    end
    @inbounds begin
        ξ1 = A[row,col]
        normu = abs2(ξ1)
        for i = col+1:n
            normu += abs2(A[row,i])
        end
        normu = sqrt(normu)
        ν = copysign(normu,real(ξ1))
        A[row,col] += ν
        ξ1 += ν
        A[row,col] = -ν
        for i = col+1:n
            A[row,i] /= ξ1
        end
    end
    return conj(ξ1/ν)
end

function elementaryRightTrapezoid!(A::AbstractMatrix, row::Integer)
    m, n = size(A)
    if !(1 <= row <= m)
        throw(BoundsError(
            "row cannot be less than one or larger than $(size(A,1))"))
    end
    @inbounds begin
        ξ1 = A[row,row]
        normu = abs2(A[row,row])
        for i = m+1:n
            normu += abs2(A[row,i])
        end
        normu = sqrt(normu)
        ν = copysign(normu,real(ξ1))
        A[row,row] += ν
        ξ1 += ν
        A[row,row] = -ν
        for i = m+1:n
            A[row,i] /= ξ1
        end
    end
    return conj(ξ1/ν)
end

function det(A::AbstractMatrix)
    if istriu(A) || istril(A)
        return det(UpperTriangular(A))
    end
    return det(lufact(A))
end
det(x::Number) = x

logdet(A::AbstractMatrix) = logdet(lufact(A))
