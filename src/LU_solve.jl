function solve_U(U::AbstractMatrix, b::AbstractVector)
  n = length(b)
  x = similar(b)
  x[n] = b[n] / U[n, n]

  for i = n-1:-1:1
    x[i] = b[i] - sum(U[i, i+1:n] .* x[i+1:n])
    x[i] /= U[i, i]
  end
  x
end

function solve_L(L::AbstractMatrix, b::AbstractVector)
  n = length(b)
  x = similar(b)
  x[1] = b[1] / L[1, 1]

  for i = 2:n
    x[i] = b[i] - sum(L[i, 1:i-1] .* x[1:i-1])
    x[i] /= L[i, i]
  end
  x
end



# L: [p, n]
# U: [n, q+1]
## 带状矩阵的版本
function solve_L(L::BandMatrix, b::AbstractVector; p=2)
  (; p=L)
  n = length(b)
  x = similar(b)
  x[1] = b[1]
  # L, [i, j] -> [i-j, j]
  for i = 2:n
    x[i] = b[i]
    for k = max(i - p, 1):i-1
      x[i] -= L.A[i-k, k] * x[k]
    end
  end
  x
end

function solve_U(U::BandMatrix, b::AbstractVector; q=1)
  n = length(b)
  x = similar(b)
  (; q=U)
  # U, [i, j] -> [i, j-i+1]
  x[n] = b[n] / U.A[n, 1]
  for i = n-1:-1:1
    for k = i+1:min(i + q, n)
      x[i] = b[i] - sum(U.A[i, k-i+1] .* x[k])
    end
    # x[i] = b[i] - sum(U[i, i+1:n] .* x[i+1:n])
    x[i] /= U.A[i, 1]
  end
  x
end
