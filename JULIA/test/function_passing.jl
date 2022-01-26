using Random

function square(x::Float64)
    return x * x
end

function cube(x::Float64)
    return x * x * x
end

function which(g::Function, square::Function, cube::Function, y::Bool)
    if y
        return square
    else
        return cube
    end
end

function main()
    x = 3.0

    function g(x::Float64) end
    y = rand(MersenneTwister())
    #g(any_function::Function) = which(y)
    g = which(g, square, cube, y < 0.5)

    println("x = $x")
    z = g(x)
    println("z = $z")

end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
