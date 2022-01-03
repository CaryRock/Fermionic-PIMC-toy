using ArgParse

function parse_commandline()
    set = ArgParseSettings()


end

function main()
#    println(PROGRAM_FILE)
#    println(ARGS)
    numArgs = length(ARGS)

#    commandString = "julia
    print("julia $(PROGRAM_FILE) ")
    for i = 1:numArgs
        print(ARGS[i] * " ")
    end
    println()

#    println("Julia command: $(Base.julia_cmd())")

#    mycmd = replace(read(joinpath("/", "proc", string(getpid()), "cmdline"), String), "\x00"=>" ")
#    println("myCmd: $mycmd")

end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
