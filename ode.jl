using OrdinaryDiffEq, Plots

let

h = 360
b = .986
c = .014
g(x) = h * exp(-((x - b) ^ 2 / (2 * c ^ 2))) - .5
q(t) = g(t % 1)

S(ζ) = .6
T(ζ) = .001 * ( -1 + exp(8 + 2 * ζ/15))
κ = .001

#t = 0:.001:10
#plot(t, q.(t))
ζ=-30:.01:60
plot(ζ, T.(ζ))
#=
tspan = (0, 100)

ζ₀ = -5

dζ(ζ, p, t) = q(t)/S(ζ) - κ*T(ζ)/S(ζ)

p = ODEProblem(dζ, ζ₀, tspan)
s = solve(p, Tsit5(); dtmax=.01)
plot(s)
=#
end
