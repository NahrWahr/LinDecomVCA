function SNRcalc(Y::Matrix,Rm,x)::Real
    L,N, = size(Y)
    p,N, = size(x)

    Py = sum(Y .^ 2)/N
    Px = sum(x .^ 2)/N + sum(Rm .^ 2)

    SNRest = 10*log10(Complex((Px - p/L*Py)/(Py - Px)))
    SNRest = norm(SNRest)
    #println("SNRest=" ,SNRest)
    return SNRest
end

function VCA(Y::Matrix{pType},
             R::Int64,
             SNRin::Float64)::Tuple{Matrix{pType},
                                    Vector{Int64},
                                    Matrix{pType}} where pType <: Union{Float32,Float64}

    L,N, = size(Y)
    ym = mean(Y, dims=2)
    Yo = Y .- ym
    Ud = svd((Yo*Yo')/N).U[:,begin:R]
    xp = Ud'*Yo

    #SNR = SNRin==0.0 ? SNRcalc(Y, ym, xp) : SNRin
    SNR = SNRin
    SNRth = 15 + 10*log10(R)

    if SNR < SNRth # Assume Y = 20x1000 & R = 3
        d = R-1 # d = 2

        Ud = Ud[:,begin:d] # Ud = 20x2
        x = Ud'*Yo # xp = 2x1000
        Yp = Ud*x .+ ym # Yp = 20x1000
        c = maximum(sum(x.^2))^0.5 # c = 1x1
        y = vcat(x, c*ones(1,N)) # y = 3x1

    else
        d = R # d = 3

        Ud = Ud[:,begin:d] # Ud = 20x3
        x = Ud'*Yo # x = 3x1000
        Yp = Ud*x .+ ym # Yp = 20x1000
        u = mean(x, dims=2) # u = 3x1
        y = x./(u'*x) # y = 3x1

    end

    indice = zeros(Int,R)
    A = zeros(R,R)
    A[end,begin] = 1

    for i=1:R
        w = rand(R,1)
        f = w - (A* (pinv(A) * w))
        f = f/norm(f)
        v = f' * y
        indice[i] = argmax(abs.(v))[2]
        A[:,i] = y[:, indice[i]]
    end

    Ae = Yp[:,indice]
    return Ae, indice, Yp
end
