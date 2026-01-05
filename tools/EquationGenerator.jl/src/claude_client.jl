"""
Claude API client for interacting with Anthropic's API.
"""

using HTTP
using JSON3

const CLAUDE_API_URL = "https://api.anthropic.com/v1/messages"
const CLAUDE_MODEL = "claude-sonnet-4-20250514"
const CLAUDE_API_VERSION = "2023-06-01"

"""
    call_claude(prompt::String, api_key::Union{String,Nothing}=nothing; max_tokens::Int=4096) -> String

Call Claude API with a prompt and return the response.

# Arguments
- `prompt::String`: The prompt to send to Claude
- `api_key::Union{String,Nothing}`: API key (defaults to ENV["ANTHROPIC_API_KEY"])
- `max_tokens::Int`: Maximum tokens in response (default: 4096)

# Returns
- String response from Claude

# Example
```julia
response = call_claude("Explain quantum physics", api_key="sk-...")
```
"""
function call_claude(
    prompt::String,
    api_key::Union{String,Nothing}=nothing;
    max_tokens::Int=4096
)
    # Get API key
    key = something(api_key, get(ENV, "ANTHROPIC_API_KEY", nothing))
    if isnothing(key) || isempty(key)
        error("API key required. Set ANTHROPIC_API_KEY environment variable or pass api_key parameter.")
    end

    # Prepare request
    headers = [
        "x-api-key" => key,
        "anthropic-version" => CLAUDE_API_VERSION,
        "content-type" => "application/json"
    ]

    body = Dict(
        "model" => CLAUDE_MODEL,
        "max_tokens" => max_tokens,
        "messages" => [
            Dict("role" => "user", "content" => prompt)
        ]
    )

    # Make request
    try
        response = HTTP.post(
            CLAUDE_API_URL,
            headers,
            JSON3.write(body)
        )

        # Parse response
        result = JSON3.read(String(response.body))

        if haskey(result, :content) && length(result.content) > 0
            return result.content[1].text
        else
            error("Unexpected response format from Claude API")
        end

    catch e
        if e isa HTTP.ExceptionRequest.StatusError
            error("Claude API error: $(e.status) - $(String(e.response.body))")
        else
            rethrow(e)
        end
    end
end

"""
    test_api_key(api_key::Union{String,Nothing}=nothing) -> Bool

Test if the API key is valid.

# Arguments
- `api_key::Union{String,Nothing}`: API key to test

# Returns
- true if valid, false otherwise

# Example
```julia
if test_api_key()
    println("API key is valid!")
end
```
"""
function test_api_key(api_key::Union{String,Nothing}=nothing)
    try
        response = call_claude("Hello", api_key; max_tokens=10)
        return !isempty(response)
    catch
        return false
    end
end
