using Random

function square(x::Float64)
    return x * x
end

function cube(x::Float64)
    return x * x * x
end

function which(first::Function, second::Function, y::Bool)
    y ? (return first) : (return second)
end

function main()
    x = 3.0

    function g() end
    #g = which(square, cube, rand() < 0.5)
    (rand() < 0.5) ? (g = square) : (g = cube)

    println("x = $x")
    println("The function is $(nameof(var"g"))")
    println("z = $(g(x))")

end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
