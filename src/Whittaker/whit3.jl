function whit3(y::AbstractVector{T1}, w::AbstractVector{T2};
  λ::Real, include_cve=true) where {T1<:Real,T2<:Real}

  interm = interm_whit{promote_type(T1, T2)}(; n=length(y))

  whit3!(y, w, λ, interm; include_cve)
  # interm.z, cve
end

function whit3!(y::AbstractVector{<:Real}, w::AbstractVector{<:Real}, λ::Real, interm::interm_whit{FT};
  include_cve=true) where {FT<:Real}
  λ = FT(λ)
  @unpack z, c, d, e, f = interm

  d[1] = w[1] + λ # d是分母
  c[1] = -3λ / d[1]
  e[1] = 3λ / d[1]
  f[1] = -λ / d[1]
  z[1] = w[1] * y[1]

  d[2] = w[2] + 10λ - c[1]^2 * d[1]
  c[2] = (-12λ - c[1] * e[1] * d[1]) / d[2]
  e[2] = (6λ + c[1] * λ) / d[2]
  f[2] = -λ / d[2]
  z[2] = w[2] * y[2] - c[1] * z[1]

  d[3] = w[3] + 19λ - c[2]^2 * d[2] - e[1]^2 * d[1]
  c[3] = (-15λ + e[1] * λ - c[2] * e[2] * d[2]) / d[3]
  e[3] = (6λ + c[2] * λ) / d[3]
  f[3] = -λ / d[3]
  z[3] = w[3] * y[3] - e[1] * z[1] - c[2] * z[2]

  n = length(y)
  @fastmath for i = 4:n-3
    d[i] = w[i] + 20λ - c[i-1]^2 * d[i-1] - e[i-2]^2 * d[i-2] - f[i-3]^2 * d[i-3]
    c[i] = (-15λ + e[i-2] * λ - c[i-1] * e[i-1] * d[i-1]) / d[i]
    e[i] = (6λ + c[i-1] * λ) / d[i]
    f[i] = -λ / d[i]
    z[i] = w[i] * y[i] - c[i-1] * z[i-1] - e[i-2] * z[i-2] - f[i-3] * z[i-3]
  end
  i = n - 2
  d[i] = w[i] + 19λ - c[i-1]^2 * d[i-1] - e[i-2]^2 * d[i-2] - f[i-3]^2 * d[i-3]
  c[i] = (-12λ + e[i-2] * λ - c[i-1] * e[i-1] * d[i-1]) / d[i]
  e[i] = (3λ + c[i-1] * λ) / d[i]
  z[i] = w[i] * y[i] - c[i-1] * z[i-1] - e[i-2] * z[i-2] - f[i-3] * z[i-3]

  i = n - 1
  d[i] = w[i] + 10λ - c[i-1]^2 * d[i-1] - e[i-2]^2 * d[i-2] - f[i-3]^2 * d[i-3]
  c[i] = (-3λ + e[i-2] * λ - c[i-1] * e[i-1] * d[i-1]) / d[i]
  z[i] = w[i] * y[i] - c[i-1] * z[i-1] - e[i-2] * z[i-2] - f[i-3] * z[i-3]

  i = n
  d[i] = w[i] + λ - c[i-1]^2 * d[i-1] - e[i-2]^2 * d[i-2] - f[i-3]^2 * d[i-3]
  z[i] = w[i] * y[i] - c[i-1] * z[i-1] - e[i-2] * z[i-2] - f[i-3] * z[i-3]

  z[n] /= d[n]
  z[n-1] = z[n-1] / d[n-1] - c[n-1] * z[n]
  z[n-2] = z[n-2] / d[n-2] - c[n-2] * z[n-1] - e[n-2] * z[n]
  @fastmath for i = (n-3):-1:1
    z[i] = z[i] / d[i] - c[i] * z[i+1] - e[i] * z[i+2] - f[i] * z[i+3]
  end

  cve::FT = include_cve ? whit3_hat(y, w, interm) : FT(-999.0)
  z, cve
end


function cal_cve(U2::AbstractMatrix, d::AbstractVector,
  z::AbstractVector{FT}, y::AbstractVector{FT}, w::AbstractVector; p=3) where {FT<:Real}
  n = length(y)
  # s = inv_diag(U2, d; m=p)
  s0 = cal_diag(U2, d; m=p)
  tol = FT(0.0)
  wsum = FT(0.0)
  @inbounds for i = 1:n
    s0[i] = s0[i] * w[i] # H = (L D L')⁻¹ W, consider W
    r::FT = (y[i] - z[i]) / (1 - s0[i]) # w[i]
    tol += sum(r * r * w[i])
    wsum += w[i]
  end
  cve::FT = sqrt(tol / wsum) # sqrt(sum(r .* r) / n)
  return cve
end

# according to hat and return the generalized cross validation
function whit3_hat(y::AbstractVector{<:Real}, w::AbstractVector{<:Real}, interm::interm_whit{FT}) where {FT<:Real}
  @unpack z, e, c, f, d = interm
  U2 = hcat(c, e, f)
  cal_cve(U2, d, z, y, w; p=3)
end


export whit3_hat
