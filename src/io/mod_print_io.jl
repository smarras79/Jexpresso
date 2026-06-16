function println_rank(args...; msg_rank::Int = 0, suppress = false)
    if suppress == true
        return
    end
    if msg_rank == 0
        println(args...)
        flush(stdout)
    end
end

function print_rank(args...; msg_rank::Int = 0, suppress = false)
    if suppress == true
        return
    end
    if msg_rank == 0
        print(args...)
        flush(stdout)
    end
end