using Optim
include("kim.jl")

y = readcsv("data.txt")

T = size(y)[1]; n = size(y)[2]
J = 1; s = 2

Q = eye(J)
P0 = eye(J)*1.0
x0 = [1.0]

l = [ kron(0.2, ones(1, J*J*s)) kron(-2.5, ones(1,n*J*s)) kron(0.4, ones(1,n)) kron(0.9, ones(1,s)) ]'
u = [ kron(0.997, ones(1, J*J*s)) kron(2.5, ones(1,n*J*s)) kron(1.0, ones(1,n)) kron(0.999, ones(1,s)) ]'

dimA = J*J*s; dimAF = dimA + n*J*s; dimAFR = dimAF + n

function KimFilterOptim(pars)
  A = reshape(pars[1:dimA], J, J, s)
  F = reshape(pars[(dimA+1):dimAF], n, J, s)
  R = diagm(pars[(dimAF+1):dimAFR])
  p = diagm(pars[(dimAFR+1):(dimAFR+2)])
  p[1,2] = 1 - p[1,1]; p[2,1] = 1 - p[2,2]
  lik = KimFilter(y,F,x0,P0,A,R,Q,p)
  return -lik
end


init = [reshape((rand(J*J*s) +0.3)/1.3, J*J*s, 1), reshape(rand(n*J*s), n*J*s, 1), kron(0.5, ones(1,n))', 0.95,0.95]
# opt = optimize(KimFilterOptim, init[:,1] , method=:l_bfgs, show_trace=true, iterations=5000)

opt = fminbox(DifferentiableFunction(KimFilterOptim), init[:,1], l[:,1], u[:,1], show_trace=true)
