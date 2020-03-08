using Plots
using Revise
using Wavelets
WT.paul10
supertype(typeof(WT.paul10))
wavelet(WT.morl)
WT.morl
wavelet(WT.morl, s=8)
plotlyjs()
J=11; n = 2^J+395
x = testfunction(n,"Bumps")
EmphasizeFrequencyInfo = WT.Morlet(4.7)
y = cwt(x, wavelet(WT.morl))
heatmap(abs.(y)); plot!(x+90,label="")
y = cwt(x, wavelet(WT.dog1))
heatmap(abs.(y)); plot!(x+90,label="")
y = cwt(x, wavelet(WT.paul4))
heatmap(abs.(y)); plot!(x+90,label="")

gr()
pyplot()
waveType = WT.Morlet(6)

Ψ = wavelet(WT.Morlet(5))
daughters, ξ = Wavelets.computeWavelets(n,Ψ)
plot(daughters,legend=false)

Ψ = wavelet(WT.Morlet(4.4))
daughters, ξ = Wavelets.computeWavelets(n,Ψ)
plot([daughters sum(daughters,dims=2)],legend=false)

Ψ = wavelet(waveType, s=1,averagingLength=2)
daughters1, ξ = Wavelets.computeWavelets(n,Ψ)
plot([daughters1 sum(daughters1,dims=2)])

Ψ = wavelet(waveType, s=2)
daughters2, ξ = Wavelets.computeWavelets(n,Ψ)
plot([daughters2 sum(daughters2,dims=2)])

Ψ = wavelet(waveType, s=4)
daughters4, ξ = Wavelets.computeWavelets(n,Ψ)
plot([daughters4 sum(daughters4,dims=2)])

waveType = WT.dog6
waveType = WT.paul1

waveType = WT.Morlet(5)
Ψ = wavelet(waveType, s=4.0, decreasing=4.0,averagingLength=4)
nOctaves, nWaveletsInOctave, totalWavelets, sRanges, sWidths =
    WT.getNWavelets(n,Ψ)
nOctaves
sRanges
nWaveletsInOctave
totalWavelets
daughters8, ω = Wavelets.computeWavelets(n,Ψ)
plt1= plot(abs.([daughters8 sum(daughters8,dims=2)]),legend=false);
vline!(2 .^ (Ψ.averagingLength:Ψ.averagingLength + nOctaves));
plt2 = plot(abs.([daughters8 sum(daughters8,dims=2)][1:512,:]),legend=false);
vline!(min.(2 .^ (Ψ.averagingLength:Ψ.averagingLength + nOctaves), 512));
plot(plt1,plt2, layout=(2,1))
totalWavelets

log2.(sRanges)
6*2
Ψ.σ[2]
k=2; sum((1:2049) .* daughters8[:,k]/sum(daughters8[:,k]))/Ψ.σ[1]
[argmax(daughters8[:,k])/Ψ.σ[1] for k=1:11]
Ψ.averagingLength
6.415

ss = 2^11; sσ = 1.0; (argmax(WT.Mother(Ψ, ss,sσ,ω)), Ψ.σ[1]*ss,ss)
plot(WT.Mother(Ψ, ss,sσ,ω))

Ψ.σ[1] * 2^(6.415+2)
effectiveQualityFactors(sRanges)'
plotSpacing(nOctaves, p, Ψ, totalWavelets)


scatter!(2 .^ (Ψ.averagingLength:Ψ.averagingLength + nOctaves), fill(.025, Int(Ψ.averagingLength + nOctaves)))

sRanges
function effectiveQualityFactors(sRanges)    
    effOct = floor.(Int,log2.(sRanges) .- .0001)
    return [count(effOct .==ii) for ii=minimum(effOct):maximum(effOct)]
end
Ψ.σ[1]
Ψ.σ[2]
Ψ.σ[3]
sum((1:size(daughters8,1)) .*daughters8[:,end]./sum(daughters8[:,end]))
sRanges
4π*sRanges[end]/(Ψ.σ[1] + sqrt(2+Ψ.σ[1]^2))

Ψ = wavelet(waveType, s=8.0, decreasing=1.0)
nOctaves, nWaveletsInOctave, totalWavelets, sRanges, sWidths = WT.getNWavelets(n,Ψ)
minimum(sWidths) ./ sWidths
sWidths
totalWavelets
WT.varianceAdjust(Ψ,totalWavelets)
p =8.0
for p=1.0:.25:5.0
    waveType = WT.Morlet(4.5)
    Ψ = wavelet(waveType, s=8.0, decreasing=p)
    nOctaves, nWaveletsInOctave, totalWavelets, sRanges, sWidths =
        WT.getNWavelets(n,Ψ)
    vAdj = WT.varianceAdjust(Ψ,totalWavelets)
    @test vAdj[1] == p
    @test vAdj[end] == 1
    @test length(vAdj) == totalWavelets
end
plot(1 .+a .*(totalWavelets .- (1:21)).^p)
plot(minimum(sWidths) ./ sWidths)
sWidths
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
                    1.0:totalWavelets,legend=false,xaxis="log freq"), 
            scatter(2 .^(WT.polySpacing(nOctaves,p,Ψ.averagingLength, totalWavelets)),
                    1:totalWavelets,legend=false,xaxis="freq"),
            layout=4)
Ψ.averagingLength
function plotSpacing(nOctaves,Ψ,totalWavelets)
    xVals = WT.polySpacing(nOctaves,Ψ.decreasing,Ψ.averagingLength, totalWavelets)
    plt1 = scatter(xVals, 1.0:totalWavelets,legend=false,xaxis="log freq")
    vline!(Array( Ψ.averagingLength: (Ψ.averagingLength+ nOctaves)))
    plt2 = scatter(2 .^(xVals), 1:totalWavelets, legend=false, xaxis="freq")
    #plot(2.^(1:totalWavelets), 1:totalWavelets)
    #vline!(2 .^ (Ψ.averagingLength:Ψ.averagingLength+totalWavelets))
    plot(plt1, plt2,layout=2)
end

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
