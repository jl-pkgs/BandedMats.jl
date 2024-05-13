"""
R version Whittaker Cross validation

Whittaker smoothing with second order differences
Computation of the hat diagonal (Hutchinson and de Hoog, 1986)

- `In` : data vector (y), weigths (w), smoothing parameter (λ)
- `Out`: list with smooth vector (z), hat diagonal (dhat)

#author: 
Gianluca Frasso and Paul HC Eilers, 2015

# References
1. Gianluca Frasso and Paul HC Eilers, L- and V-curves for optimal smoothing, 2015
"""
function whit2_Frasso2015(y::AbstractVector{T1}, w::AbstractVector{T2};
  λ=2.0, include_cve=true) where {T1,T2}

  T = promote_type(T1, T2)
  # w = y*0 .+ 1.0
  n = length(y)
  g0 = ones(T, n) * 6 #rep(6, n)

  g0[1] = g0[n] = 1
  g0[2] = g0[n-1] = 5
  g1 = ones(T, n) * -4
  g1[1] = g1[n-1] = -2
  g1[n] = 0
  g2 = ones(T, n)
  g2[n-1] = 0
  g2[n] = 0

  # Store matrix G = W + λ * D’ * D in vectors
  g0 = g0 * λ .+ w
  g1 = g1 * λ
  g2 = g2 * λ
  # Compute U’VU decomposition (upper triangular U, diagonal V)
  # print(g0)
  v = g0
  u1 = zeros(T, n)
  u2 = zeros(T, n)

  @inbounds for i = 1:n
    if (i > 1)
      v[i] -= v[i-1] * u1[i-1]^2
    end
    if (i > 2)
      v[i] -= v[i-2] * u2[i-2]^2
    end

    if (i < n)
      u = g1[i]
      if (i > 1)
        u = u - v[i-1] * u1[i-1] * u2[i-1]
      end
      u1[i] = u / v[i]
    end
    if (i < n - 1)
      u2[i] = g2[i] / v[i]
    end
  end
  # g0, g1, and g2 can be clear now

  # Solve for smooth vector
  z = 0 * y
  @inbounds for i = 1:n
    z[i] = y[i] * w[i]
    if (i > 1)
      z[i] -= u1[i-1] * z[i-1]
    end
    if (i > 2)
      z[i] -= u2[i-2] * z[i-2]
    end
  end
  z = z ./ v

  @inbounds for i = n:-1:1
    if (i < n)
      z[i] -= u1[i] * z[i+1]
    end
    if (i < n - 1)
      z[i] -= u2[i] * z[i+2]
    end
  end

  cve = T(-999.0)
  if include_cve
    s0 = zeros(T, n)
    s1 = zeros(T, n)
    s2 = zeros(T, n)

    # Compute diagonal of inverse
    # params: v, u1, u2, s0, s1, s2
    @inbounds for i = n:-1:1
      i1 = i + 1
      i2 = i + 2
      s0[i] = 1 / v[i]
      if (i < n)
        s1[i] = -u1[i] * s0[i1]
        s0[i] = 1 / v[i] - u1[i] * s1[i]
      end
      if (i < n - 1)
        s1[i] = -u1[i] * s0[i1] - u2[i] * s1[i1]
        s2[i] = -u1[i] * s1[i1] - u2[i] * s0[i2]
        s0[i] = 1 / v[i] - u1[i] * s1[i] - u2[i] * s2[i]
      end
    end
    
    tol = T(0.0)
    wtol = T(0.0)
    @inbounds for i = 1:n
      s0[i] *= w[i]
      r::T = (y[i] - z[i]) / (1 - s0[i])
      tol += sum(r * r * w[i])
      wtol += w[i]
    end
    cve::T = sqrt(tol / wtol) # sqrt(sum(r .* r) / n)
    # r = @. (y - z) / (1 - s0)
    # cve = sqrt(sum(r .* r / n))
  end
  # return(list(z = z, dhat = s0, cve))
  z, cve
end


export whit2_Frasso2015
