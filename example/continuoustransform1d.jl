using Plots
using Revise
using Wavelets
WT.paul10
supertype(typeof(WT.paul10))
wavelet(WT.morl)
WT.morl
wavelet(WT.morl, s=8)
plotlyjs()
J=11; n = 2^J
x = testfunction(n,"Bumps")
EmphasizeFrequencyInfo = WT.Morlet(20)
y = cwt(x, wavelet(WT.morl))
heatmap(abs.(y)); plot!(x+90,label="")
y = cwt(x, wavelet(WT.dog1))
heatmap(abs.(y)); plot!(x+90,label="")
y = cwt(x, wavelet(WT.paul4))
heatmap(abs.(y)); plot!(x+90,label="")

gr()
pyplot()
waveType = WT.Morlet(6)

Ψ = wavelet(waveType, s=1,averagingLength=2)
daughters1, ξ = Wavelets.computeWavelets(n,Ψ)
plot([daughters1 sum(daughters1,dims=2)])

Ψ = wavelet(waveType, s=2)
daughters2, ξ = Wavelets.computeWavelets(n,Ψ)
plot([daughters2 sum(daughters2,dims=2)])

Ψ = wavelet(waveType, s=4)
daughters4, ξ = Wavelets.computeWavelets(n,Ψ)
plot([daughters4 sum(daughters4,dims=2)])

waveType = WT.Morlet(5)
Ψ = wavelet(waveType, s=8.0, decreasing=2.0)
daughters8, ξ = Wavelets.computeWavelets(n,Ψ)
plot([daughters8 sum(daughters8,dims=2)],legend=false)

Ψ = wavelet(waveType, s=8.0, decreasing=1.0)
nOctaves, nWaveletsInOctave, totalWavelets, sRanges, sWidths = WT.getNWavelets(n,Ψ)
plot(sWidths)
sRanges
nWaveletsInOctave
linRange = Array{Array{Float64}}(undef, length(nWaveletsInOctave))
startAt(curOctave,nWavelets) = curOctave > 1 ? 0-1 ./ nWavelets[curOctave-1] : 0
stopAt(curOctave,nWavelets) = curOctave < length(nWavelets) ? 1 +
    1 ./ nWavelets[curOctave+1] : 1
for curOctave = 1:length(nWaveletsInOctave)
    linRange[curOctave] = (range(0,
                                 stop=1, length =
                                nWaveletsInOctave[curOctave]+1) .+ curOctave .+ Ψ.averagingLength .- 1)[2:end]
end
gathered = cat(linRange..., dims=1)

justUse = 1:21#[1,3,5,7, 9:21...]
p = 4; plot(scatter(gathered, 1:length(gathered),xaxis = "log freq", yaxis="index", legend=false), scatter(2 .^ gathered,
                                                                                                           1:length(gathered), legend=false),
            scatter((5 .+(1.0:21)).^(1/p),
                    1.0:21,legend=false,title="proposed, 'linear'"),scatter(2 .^((5 .+(1.0:21)).^(1/p))[justUse],
                                                                            1:length(justUse),legend=false,title="proposed
                                                                     exp"),
            layout=4)

p = 8; plot(scatter(gathered, 1:length(gathered),xaxis = "log freq", yaxis="index", legend=false), scatter(2 .^ gathered,
                                                                                                           1:length(gathered), legend=false),
            scatter(WT.polySpacing(nOctaves,p,Ψ.averagingLength, totalWavelets),
                    1.0:totalWavelets,legend=false,title="proposed, 'linear'"), 
            scatter(2 .^(WT.polySpacing(nOctaves,p,Ψ.averagingLength, totalWavelets)),
                    1:totalWavelets,legend=false,title="proposed exp"),
            layout=4)

nOctaves + Ψ.averagingLength+1 -1
sp = WT.polySpacing(nOctaves,4,Ψ.averagingLength, totalWavelets)
Δ
Δ = diff(2 .^ sp)./2
[Δ[1]; Δ]
2 .^ sp
scatter(2 .^(sp),1:length(sp), xerror=[Δ[1]; Δ])
gathered
scatter(WT.polySpacing(nOctaves,4,Ψ.averagingLength, totalWavelets),legend=false)
nOctaves

spaced(nPoints,)
sPoints(nO,Q,β,Ψ) = Q .* ((1.0*nO .- 
                           (Ψ.averagingLength:(nO-1))).^(-β+1))
sPoints(11,8,.9,Ψ)
length(4:10)
Ψ.averagingLength:(11-1)


function getScales(n1::Integer, c::CFW)
    nOctaves = log2(max(n1, 2)) - c.averagingLength
    nWaveletsInOctave = reverse([max(1, round(Int, c.scalingFactor /
                                              x^(c.decreasing))) for
                                 x=1:round(Int, nOctaves)])
    sj = zeros(sum(nWaveletsInOctave))
    for curOctave = 1:round(Int, nOctaves)
        sRange = (2 .^ (range(0, stop=1, length = nWaveletsInOctave[curOctave]+1) .+
                        curOctave .+ c.averagingLength .- 1))[1:end-1]
        ind = (sum(nWaveletsInOctave[1:curOctave-1])+1):sum(nWaveletsInOctave[1:curOctave])
        @debug "ind = $ind"
        sj[ind] = sRange
    end
    return sj
end
