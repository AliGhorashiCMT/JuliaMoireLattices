#=
We provide functionality to create a Moire structure with an additional translation
=#
using ArgParse
using LinearAlgebra
a=1.44*sqrt(3)/0.529
#Original lattice
a1 = a*[1, 0, 0]
a2 = a*[1/2, sqrt(3)/2, 0]
function main(args)
    q, p = parse.(Int, ARGS[1:2])
    dist = parse(Float64, ARGS[3])
    slidex, slidey = parse.(Float64, ARGS[4:5])
    
    moiresize = gcd(p, 3)/(gcd(3*q+p, 3*q-p)^2)*(3*q^2+p^2)
    println("Moire size is: ", moiresize)
    num = 3*q^2-p^2
    denum = 3*q^2+p^2
    println("Running for q = $(q), p = $(p)")
    theta = acos(num/denum)
    println("Theta:  $(theta*180/pi)")
    #Rotated lattice
    b1 = a*[cos(theta), sin(theta), 0]
    b2 = a*[cos(pi/3+theta), sin(pi/3+theta), 0]
    ValidTuples = []
    ntries = 20
    for (ii, jj, kk, ll) in Tuple.(CartesianIndices(rand(ntries, ntries, ntries, ntries)))
            i = ii-ntries/2
            j = jj-ntries/2
            k = kk-ntries/2
            l = ll-ntries/2
            diff = a1*i+a2*j-b1*k-b2*l
            isapprox(diff[1]^2+diff[2]^2, 0, atol=1e-10) || continue
            ((-i, -j, -k, -l) in ValidTuples) && continue
            isapprox(moiresize*a*a, dot(a1*i+a2*j, a1*i+a2*j), atol=1e-7) || continue
            ((a1*i+a2*j)[2] < 0) && continue
            push!(ValidTuples, (i, j, k, l))
            println(i, " ", j, " ", k, " ", l)
    end
    latvec1 = ValidTuples[1][1]*a1+ValidTuples[1][2]*a2
    latvec2 = ValidTuples[3][1]*a1+ValidTuples[3][2]*a2

    println(ValidTuples)
    println(latvec1)
    println(latvec2)
    println(dot(latvec1, latvec2)/sqrt(dot(latvec1, latvec1)*dot(latvec2, latvec2)))
    t1 = latvec1;
    t2 = latvec2

    atot_transform = [ValidTuples[1][1] ValidTuples[3][1]; ValidTuples[1][2] ValidTuples[3][2]]
    btot_transform = [ValidTuples[1][3] ValidTuples[3][3]; ValidTuples[1][4] ValidTuples[3][4]]
    ttoa_transform = inv(atot_transform)
    ttob_transform = inv(btot_transform)

    open("Moire.lattice", "w") do io
            write(io, "lattice\\ \n")
            write(io, "$(t1[1]) $(t2[1]) 0 \\ \n")
            write(io, "$(t1[2]) $(t2[2]) 0 \\ \n")
            write(io, "0 0 40\n")
    end
    nionpos = 10
    open("Moire.ionpos", "w") do io
        for (niter, miter) in Tuple.(CartesianIndices(rand(2*nionpos+1, 2*nionpos+1)))
            n, m = niter-nionpos-1, miter-nionpos-1
            nprime, mprime = ttoa_transform*[n, m]
                ((0 <= nprime < 0.999999999) && (0<= mprime < 0.999999999)) || continue 
                write(io, "ion B $(nprime) $(mprime) 0 1 \n")
        end
        write(io, "\n\n")
        for (niter, miter) in Tuple.(CartesianIndices(rand(2*nionpos+1, 2*nionpos+1)))
            n, m = niter-nionpos-1, miter-nionpos-1
            nprime, mprime = ttoa_transform*[n-1/3, m+2/3]
            (0 <= nprime < 0.999999999 && 0<= mprime < 0.999999999) || continue 
            write(io, "ion N $(nprime) $(mprime) 0 1 \n")
        end
        write(io, "\n\n")
        for (niter, miter) in Tuple.(CartesianIndices(rand(2*nionpos+1, 2*nionpos+1)))
            n, m = niter-nionpos-1, miter-nionpos-1
            nprime, mprime = ttob_transform*[n, m]
            (0 <= nprime < 0.999999999 && 0<= mprime < 0.999999999) || continue 
            write(io, "ion B $(nprime+1/slidex) $(mprime+1/slidey) $(dist/40) 1 \n")
        end
        write(io, "\n\n")
        for (niter, miter) in Tuple.(CartesianIndices(rand(2*nionpos+1, 2*nionpos+1)))
            n, m = niter-nionpos-1, miter-nionpos-1
            nprime, mprime = ttob_transform*[n-1/3, m+2/3]
            (0 <= nprime < 0.999999999 && 0<= mprime < 0.999999999) || continue 
            write(io, "ion N $(nprime+1/slidex) $(mprime+1/slidey) $(dist/40) 1 \n")
        end
    end
end
main(ARGS)
