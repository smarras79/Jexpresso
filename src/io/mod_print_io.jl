function println_rank(args...; msg_rank::Int = 0)
    if msg_rank == 0
        println(args...)
    end
end

function print_rank(args...; msg_rank::Int = 0)
    if msg_rank == 0
        print(args...)
    end
end