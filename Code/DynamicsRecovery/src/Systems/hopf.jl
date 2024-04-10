function hopf_normal!(du, u, p, t)
    # μ: p[1], ω: p[2], A: p[3]
    
    du[1] = p[1]*u[1] - p[2]*u[2] - p[3]*u[1]*(u[1]^2 + u[2]^2)
    du[2] = p[2]*u[1] + p[1]*u[2] - p[3]*u[2]*(u[1]^2 + u[2]^2)
end # function hopf_normal!

