function KimFilter(y,F,x0,P0,A,R,Q,p)
  T = size(y)[1]
  n = size(F)[1]
  J = size(x0)[1]
  s = size(p)[1]

  x = fill(0.0, T,J,s,s)
  xU = fill(0.0, T,J,s,s)
  P = fill(0.0, T,J,J,s,s)
  Pu = fill(0.0, T,J,J,s,s)
  eta = fill(0.0, T,n,s,s)
  H = fill(0.0, T,n,n,s,s)
  K = fill(0.0, T,J,n,s,s)
  lik = fill(0.0, T,s,s)
  loglik = fill(0.0, T,s,s)
  jointP_fut = fill(0.0, T,s,s)
  jointP_cur = fill(0.0, (T+1),s,s)
  stateP_fut = fill(0.0, T,s)
  stateP = fill(0.0, T,s)
  xA = fill(0.0, T,J,s)
  Pa = fill(0.0, T,J,J,s)
  result = fill(0.0 , T,1)

  for i in 1:s
    xA[1,:,i] = x0
  end
  for (i in 1:s)
    Pa[1,:,:,i] = P0
  end
  jointP_cur[1,:,:] = [0.25 0.25; 0.25 0.25]

  for t in 2:T
    for j in 1:s
      for i in 1:s
        x[t,:,i,j] = A[:,:,j] * xA[(t-1),:,i]
        P[t,:,:,i,j] = A[:,:,j] * squeeze(Pa[(t-1),:,:,i],1) * A[:,:,j]' + Q
        eta[t,:,i,j] = y[t,:]' - F[:,:,j] * x[t,:,i,j]
        H[t,:,:,i,j] = F[:,:,j] * squeeze(P[t,:,:,i,j],1) * F[:,:,j]' + R
        K[t,:,:,i,j] = squeeze(P[t,:,:,i,j],1) * F[:,:,j]' * inv(squeeze(H[t,:,:,i,j],1))
        xU[t,:,i,j] = x[t,:,i,j] + squeeze(K[t,:,:,i,j],1) * eta[t,:,i,j]'
        Pu[t,:,:,i,j] = (eye(J) - squeeze(K[t,:,:,i,j],1) * F[:,:,j]) * squeeze(P[t,:,:,i,j],1)
        jointP_fut[t,i,j] = p[i,j]*sum(jointP_cur[(t-1),:,i])
        lik[t,i,j] = abs((2*pi)^(-n/2) * det(squeeze(H[t,:,:,i,j],1))^(-1/2) * jointP_fut[t,i,j] *
                     exp(-1/2*eta[t,:,i,j] * inv(squeeze(H[t,:,:,i,j],1)) * eta[t,:,i,j]')[1,1])
        loglik[t,i,j] = log(lik[t,i,j])
        jointP_cur[t,i,j] = lik[t,i,j]
      end
      stateP[t,j] = sum(jointP_cur[t,:,j])
      stateP_fut[t,j] = sum(jointP_fut[t,:,j])
      xA[t,:,j] = squeeze(xU[t,:,:,j], 1) * squeeze(jointP_cur[t,:,j],1) / stateP[t,j]
      for (i in 1:s)
        Pa[t,:,:,j] = Pa[t,:,:,j] +
                      (Pu[t,:,:,i,j] + (xA[t,:,j] - xU[t,:,i,j]) * transpose(xA[t,:,j] - xU[t,:,i,j])) * 
                      exp(log(jointP_cur[t,i,j]) - log(stateP[t,j]))
      end
    end
    jointP_cur[t,:,:] = exp(log(jointP_cur[t,:,:]) - log(sum(lik[t,:,:])))
    stateP[t,:] = exp(log(stateP[t,:]) - log(sum(lik[t,:,:])))
    result[t,1] = log(sum(lik[t,:,:]))
  end
return sum(result)
end
