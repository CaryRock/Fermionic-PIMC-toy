using Random
using LoopVectorization

function test_simd_turbo()
    n = 5
    sqn = n^2
    r = collect(1:n)
    foo1 = circshift(r', (0, 1))
    foo2 = circshift(r', (0, -1))
    foo3 = circshift(r, 1)
    foo4 = circshift(r, -1)

    # Initialize array1 with random entries
    array1 = randn(n, n);

    # Initialize array2 with zeros
    array2 = zeros(n, n);

    # Initialize temporary arrays used in loops
    bar = zeros(n, n)
    adding = zeros(n, n)

    loop = 0
    while minimum(array1) < -0.8 && loop < 100000
        loop += 1
        fill!(bar, 0)

        # Use simd?
        @turbo for i = 1:sqn
            bar[i] = max(0, array1[i] - 0.2)
        end

        @turbo @. array2 += bar

        fill!(adding, 0)
        # Use simd?
        @turbo for i = 1:n
            for j = 1:n
                s = bar[i, foo1[j]] + bar[i, foo2[j]] + bar[foo3[i], j] + bar[foo4[i], j]
                adding[i, j] = 0.25 * s
            end
        end
        @turbo @. array1 = array1 - bar + adding
    end
    return array1, array2, loop
end
