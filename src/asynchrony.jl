using Statistics: cov, std, mean

# sol => couldd be subset to obtain the desired co
function asynchrony(sol, wind = 500)
    ws1 = floor(Int, wind / 2)
    ws2 = wind - ws1
    out = zeros(size(sol)[2] - wind)
    for i in (ws1 + 1):(size(sol)[2] - ws2)
        out[i - ws1] = cov(sol[1, (i - ws1) : (i + ws2 - 1)], sol[2, (i - ws1) : (i + ws2 - 1)])
    end
    return out
end

# Coeficient of variantion with moving window
function cv(val, wind = 500)
    ws1 = floor(Int, wind / 2)
    ws2 = wind - ws1
    out = zeros(length(val) - wind)
    for i in (ws1 + 1):(length(val) - ws2)
        tmp = val[(i - ws1):(i + ws2 - 1)]
        out[i - ws1] = std(tmp) / mean(tmp)
    end
    return out
end