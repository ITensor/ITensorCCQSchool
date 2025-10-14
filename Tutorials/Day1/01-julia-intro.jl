# Approximate `pi` using the Leibniz formula.
# See: https://en.wikipedia.org/wiki/Leibniz_formula_for_%CF%80
function approx_pi(n)
    pi_by_4 = 1.0
    for k in 1:n
        pi_by_4 += (-1)^k / (2k + 1)
    end
    return 4 * pi_by_4
end

# We use a `let` block to create a local scope, which
# is good practice in Julia to avoid introducing global
# variables which can negatively impact performance.
function main(; n = 10^5)
    # Print the number of terms to use to approximate `pi`.
    @show n

    # Approximate `pi` by calling the function defined above.
    pi_approx = approx_pi(n)
    @show pi_approx

    # Compute the error in the approximation.
    # `pi` is a built-in constant in Julia that represents
    # the value of pi to arbitrary precision.
    error = abs(pi_approx - pi)
    @show error

    # Return the result and error to the REPL.
    return (; n, pi_approx, error)
end
