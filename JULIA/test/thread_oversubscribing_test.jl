using .Threads

function main()
    println("The current number of threads in use: $(nthreads())")

    N = 2^32
    println("N = $N")
    @time begin
        a = zeros(N)
    end

    @time begin
        @threads for i = 1:N
        a[i] = Threads.threadid()
        end
    end
    
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
