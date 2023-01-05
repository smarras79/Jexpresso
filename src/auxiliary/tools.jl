using BenchmarkTools

d1 = 10000
d2 = 10000

data = randn(d1, d2)

function sum_row(data)
    answer = eltype(data)(0)
    for cnt1 = 1:size(data,1)
        for cnt2 = 1:size(data, 2)
            answer += data[cnt1, cnt2]
        end
    end
    answer
end

function sum_col(data)
    answer = eltype(data)(0)
    for cnt2 = 1:size(data,2)
        for cnt1 = 1:size(data, 1)
            answer += data[cnt1, cnt2]
        end
    end
    answer
end

#
# sum_col is WAY faster than row_col
#
@time sum_row(data)
@time sum_col(data)
