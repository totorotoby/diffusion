using LinearAlgebra
using OrdinaryDiffEq
using DifferentialEquations
using Plots
using Printf

function d(x,y)
  return 4*x
end

function dn(u)
  return u*u
end

function gaussian(x, y; x0=0.5, y0=0.5, σx=0.1, σy=0.1, A=1.0)
    return A * exp(-((x - x0)^2 / (2 * σx^2) + (y - y0)^2 / (2 * σy^2)))
end

function init_func(x,y)
   return cos(x)
end

function forcing(x, y, t)

end

# Define a DiscreteCallback to print at each accepted time step
function condition(u, t, integrator)
    return true  # Trigger at every step
end

function affect!(integrator)
    p = integrator.p
    u = integrator.u
    X = p.X
    Y = p.Y
    display(surface(X,Y,u, xlims=(0 , 1), ylims=(0,1), zlims=(0, 1)))
    return nothing
end

function tstepper(du, u, params ,t)

  n = params.n
  p = params.p
  D = params.D
  M = params.M
  X = params.X
  Y = params.Y
  h = params.h

  # bottom
  u[1, 1:n] .= 0
  #top
  u[n, 1:n] .= 0
  #u[1:n, 1] .= 0
  #u[1:n, n] .= 0
  # left - neumann
  du[1:n, 1] .= (4*u[1:n, 2] - u[1:n, 3])/3
  # y derivate
  for i = p:n-(p-1) 
    for sign = -1:2:1  
      for m = 0:(p-1)
        for n = 0:(p-1)
           du[i, 1] += M[m+1, n+1] * D[i + (sign * m), 1] * u[i + (sign * n), 1]
        end  
      end
    end
  end
  # right



  # loop through interior nodes
  for i = p:n-(p-1)
    for j = 1:n
      # compute derivates
      dx = 0
      dy = 0
      for sign = -1:2:1  
        for m = 0:(p-1)
          for n = 0:(p-1)
            if j == 1 || j == n
              du[1:n, n] .= (4*u[1:n, n-1] - u[1:n, n-2])/3
              dy += M[m+1, n+1] * dn(u[i + (sign * m), j]) * u[i + (sign * n), j]
            else
              dy += M[m+1, n+1] * dn(u[i + (sign * m), j]) * u[i + (sign * n), j]
              dx +=  M[m+1, n+1] * dn(u[i, j + (sign * m)]) * u[i, j + (sign * n)]
            end
          end
        end
      end
      du[i,j] = dx + dy
    end
  end

end

let

  p = 2
  n = 150
  L = 1
  h = L/(n - 1)

  x = 0:h:L
  y = 0:h:L

  X = [xi for xi in x, yi in y]
  Y = [yi for xi in x, yi in y]

  # solution length = n^2
  U = gaussian.(X, Y)
  # coefficent length n^2
  D = d.(X, Y)

  # interior stencil
  if p == 4
    M = [1/8 -1/6 1/24 0;
         -1/6 -3/8 2/3 -1/8;
         1/8 -2/3 3/8 1/6;
         0 -1/24 1/6 -1/8]
  elseif p == 2
    M = [ -1/2 1/2;
          -1/2 1/2 ] 
  end


  params = (n = n,
            p = p,
            D = D,
            M = M,
            X = X,
            Y = Y,
            h = h,)

  tspan = (0.0, 100)
  prob = ODEProblem(tstepper, U, tspan, params)
  integrator = init(prob, Tsit5(), dt=0.1)
  cb = DiscreteCallback(condition, affect!)  
  sol = solve(prob, Tsit5(), callback=cb)

  nothing

end
