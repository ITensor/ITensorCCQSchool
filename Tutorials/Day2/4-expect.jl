using ITensorMPS: random_mps, op, siteinds, linkinds
using ITensors: ITensor, sim, prime

function main()
    N = 10
    s = siteinds("S=1/2", N)
    psi = random_mps(s; linkdims = 2)
    # Make a new version of `psi` with randomized link indices.
    dpsi = sim(linkinds, psi)

    O = "X"
    Ls = Vector{ITensor}(undef, N)
    L = ITensor(1.0)
    for j in 1:N
        L = L * psi[j] * dpsi[j]
        Ls[j] = L
    end

    # Make a similar vector of right environments `Rs`
    # by contracting right to left.

    # Use `Ls` and `Rs` as environments to compute the expectation value
    # of operator `O` on every site. As a reminder, you can construct
    # the operator on site `j` with `op(O, s[j])`.

    return (;)
end