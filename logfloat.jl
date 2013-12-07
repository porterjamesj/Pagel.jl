immutable LogFloat <: FloatingPoint
    exp::Float64
end

# printing
import Base.show
Base.show(io::IO, lf::LogFloat) = print(io, "e^$(lf.exp)")

# converting from Float64
import Base.convert
convert(::Type{LogFloat},float::Float64) = LogFloat(float)

# arithmetic
importall Base.Operators
function (*)(a::LogFloat,b::LogFloat)
    return LogFloat(a.exp + b.exp)
end

function (/)(a::LogFloat,b::LogFloat)
    return LogFloat(a.exp - b.exp)
end

function (+)(a::LogFloat, b::LogFloat)
    af = a.exp
    bf = b.exp
    # ensure a is the maximum
    if (bf > af)
        tmp = af
        af = bf
        bf = tmp
    end

    diff = bf-af
    if af == -Inf || diff < -20
        return a
    else
        return LogFloat(af + log(1.0+exp(diff)))
    end
end
