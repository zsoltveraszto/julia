# This file is a part of Julia. License is MIT: http://julialang.org/license

## hashing a single value ##

hash(x::Any) = hash(x, zero(UInt))
hash(w::WeakRef, h::UInt) = hash(w.value, h)

## hashing general objects ##

hash(x::ANY, h::UInt) = 3*object_id(x) - h

## core data hashing functions ##

function hash_64_64(n::UInt64)
    local a::UInt64 = n
    a = ~a + a << 21
    a =  a $ a >> 24
    a =  a + a << 3 + a << 8
    a =  a $ a >> 14
    a =  a + a << 2 + a << 4
    a =  a $ a >> 28
    a =  a + a << 31
    return a
end

function hash_64_32(n::UInt64)
    local a::UInt64 = n
    a = ~a + a << 18
    a =  a $ a >> 31
    a =  a * 21
    a =  a $ a >> 11
    a =  a + a << 6
    a =  a $ a >> 22
    return a % UInt32
end

function hash_32_32(n::UInt32)
    local a::UInt32 = n
    a = a + 0x7ed55d16 + a << 12
    a = a $ 0xc761c23c $ a >> 19
    a = a + 0x165667b1 + a << 5
    a = a + 0xd3a2646c $ a << 9
    a = a + 0xfd7046c5 + a << 3
    a = a $ 0xb55a4f09 $ a >> 16
    return a
end

if UInt === UInt64
    hash_uint64(x::UInt64) = hash_64_64(x)
    hash_uint(x::UInt)     = hash_64_64(x)
else
    hash_uint64(x::UInt64) = hash_64_32(x)
    hash_uint(x::UInt)     = hash_32_32(x)
end

## symbol & expression hashing ##

hash(x::Symbol, h::UInt) = 3*object_id(x) - h
if UInt === UInt64
    hash(x::Expr, h::UInt) = hash(x.args, hash(x.head, h + 0x83c7900696d26dc6))
else
    hash(x::Expr, h::UInt) = hash(x.args, hash(x.head, h + 0x96d26dc6))
end

hash(x::QuoteNode, h::UInt) = hash(x.value, hash(QuoteNode, h))

## hashing collections ##

const hashaa_seed = UInt === UInt64 ? 0x7f53e68ceb575e76 : 0xeb575e76
const hashrle_seed = UInt == UInt64 ? 0x2aab8909bfea414c : 0xbfea414c
function hash{T}(a::AbstractArray{T}, h::UInt)
    if isleaftype(T)
        if method_exists(-, (T, T))
            val = (x1, x2) -> x2 - x1
        else
            val = (x1, x2) -> x2
        end
    else
        val = (x1, x2) -> applicable(-, x2, x1) ? x2 - x1 : x2
    end

    _hash(a, h, val)
end

function _hash{T}(a::AbstractArray{T}, h::UInt, val::Function)
    h += hashaa_seed
    h += hash(size(a))

    state = start(a)
    done(a, state) && return h
    x1, state = next(a, state)
    # Always hash the first element
    h = hash(x1, h)
    done(a, state) && return h

    # Then hash the difference between two subsequent elements when - is supported,
    # or the elements themselves when not
    x2, state = next(a, state)
    v2 = val(x1, x2)
    done(a, state) && return hash(v2, h)

    v1 = v2
    while !done(a, state)
        x1 = x2
        x2, state = next(a, state)
        v1 = v2
        v2 = val(x1, x2)
        if isequal(v2, v1)
            # For repeated elements, use run length encoding
            # This allows efficient hashing of sparse arrays
            runlength = 2
            while !done(a, state)
                x1 = x2
                x2, state = next(a, state)
                v2 = val(x1, x2)
                isequal(v1, v2) || break
                runlength += 1
            end
            h += hashrle_seed
            h = hash(runlength, h)
        end
        h = hash(v1, h)
    end
    !isequal(v2, v1) && (h = hash(v2, h))
    return h
end

# hashaa_seed and hashrle_seed are defined in abstractarray.jl
function hash(r::Range, h::UInt)
    h += hashaa_seed
    h += hash(size(r))

    length(r) == 0 && return h

    h = hash(first(r), h)
    length(r) == 1 && return h
    length(r) == 2 && return hash(step(r), h)

    h += hashrle_seed
    h = hash(length(r)-1, h)
    hash(step(r), h)
end
