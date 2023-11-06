
"""
    @missinginput "Error message"
Macro used to raise an error, when something is not implemented.
"""
macro missinginput(message=" # Missing input variable in user_input")
    quote
        error($(esc(message)))
    end
end


macro mystop(message=" MY STOP HERE")
    quote
        error($(esc(message)))
    end
end
