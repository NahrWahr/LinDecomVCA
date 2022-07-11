using Statistics
using LinearAlgebra

function SNRcalc(Y::Matrix,Rm,x)::Real
    L,N, = size(Y)
    p,N, = size(x)

    Py = sum(Y .^ 2)/N
    Px = sum(x .^ 2)/N + sum(Rm .^ 2)

    SNRest = 10*log10(Complex((Px - p/L*Py)/(Py - Px)))
    return norm(SNRest)
end

function VCA(Y::Matrix{pType},
             R::Int64,
             SNRin::Float64)::Tuple{Matrix{Float64},
                                    Vector{Int64},
                                    Matrix{Float64}} where pType <: Real

    L,N, = size(Y)

    if SNRin == 0
        ym = mean(Y, dims=2)
        Yo = Y .- ym
        Ud = svd((Yo*Yo')/N).U[:,begin:R]
        xp = Ud'*Yo

        SNR = SNRcalc(Y, ym, xp)
    else
        SNR = SNRin
    end
    SNRth = 15 + 10*log10(R)

    if SNR < SNRth
        d = R-1

        if SNRin == 0
            Ud = Ud[:,begin:d]

        else
            ym = mean(Y, dims=2)
            Yo = Y .- ym
            Ud = svd((Yo*Yo')/N).U[:,begin:d]
            xp = Ud'*Yo
        end

        Yp = Ud*xp[begin:d,:] .+ ym
        x = xp[begin:d,:]
        c = maximum(sum(x.^2))^0.5
        y = vcat(x, c*ones(1,N))

    else
        d = R
        Ud = svd((Yo*Yo')/N).U[:,begin:d]
        xp = Ud'*Y
        Yp = Ud*xp[begin:d,:]
        x = Ud'* Y
        u = mean(x, dims=2)
        y = x/(u'*x)

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
