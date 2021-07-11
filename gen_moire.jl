using ArgParse
using LinearAlgebra
a=1.44*sqrt(3)/0.529
#Original lattice
a1 = a*[1, 0, 0]
a2 = a*[1/2, sqrt(3)/2, 0]
function main(args)
        q, p = parse.(Int, ARGS)
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
        open("Moire.lattice", "w") do io
                write(io, "lattice\\ \n")
                write(io, "$(t1[1]) $(t2[1]) 0 \\ \n")
                write(io, "$(t1[2]) $(t2[2]) 0 \\ \n")
                write(io, "0 0 20\n")
        end
        nionpos = 10
        open("Moire.ionpos", "w") do io
                write(io, "ion B 0 0 0 1\n")
                write(io, "ion B 0 0 6 1\n")
                for n in -nionpos:nionpos
                        for m in -nionpos:nionpos
                                pos = n*a1+m*a2
                                thetapos = acos(dot(pos, t1)/(sqrt(dot(pos, pos))*sqrt(dot(t1, t1))))
                                thetapos < pi/3 || continue
                                thetapos > 0  || continue
                        
                                thetapos3 = acos(dot(pos, t2)/(sqrt(dot(pos, pos))*sqrt(dot(t1, t1))))
                                thetapos3 < pi/3 || continue
                                thetapos3 > 0  || continue

                                thetapos2 = acos(dot(pos-t1-t2, -t1)/(sqrt(dot(pos-t1-t2, pos-t1-t2))*sqrt(dot(t1, t1)))) 
                                thetapos2 < pi/3 || continue
                                thetapos2 > 0  || continue
                        
                                thetapos4 = acos(dot(pos-t1-t2, -t2)/(sqrt(dot(pos-t1-t2, pos-t1-t2))*sqrt(dot(t1, t1))))
                                thetapos4 < pi/3 || continue
                                thetapos4 > 0  || continue
                                write(io, "ion B $(pos[1]) $(pos[2]) 0 1 \n")
                        end
                end
                write(io, "\n\n")
                for n in -nionpos:nionpos
                        for m in -nionpos:nionpos
                                pos = [0, 1.44/0.529, 0] +n*a1+m*a2
                                thetapos = acos(dot(pos, t1)/(sqrt(dot(pos, pos))*sqrt(dot(t1, t1))))#acos(dot(pos, [1, 0, 0]))
                                thetapos < pi/3 || continue
                                thetapos > 0  || continue

                                thetapos3 = acos(dot(pos, t2)/(sqrt(dot(pos, pos))*sqrt(dot(t1, t1))))#acos(dot(pos, [1, 0, 0]))
                                thetapos3 < pi/3 || continue
                                thetapos3 > 0  || continue

                                thetapos2 = acos(dot(pos-t1-t2, -t1)/(sqrt(dot(pos-t1-t2, pos-t1-t2))*sqrt(dot(t1, t1)))) #pos-t1-t2
                                thetapos2 < pi/3 || continue
                                thetapos2 > 0  || continue

                                thetapos4 = acos(dot(pos-t1-t2, -t2)/(sqrt(dot(pos-t1-t2, pos-t1-t2))*sqrt(dot(t1, t1)))) #pos-t1-t2
                                thetapos4 < pi/3 || continue
                                thetapos4 > 0  || continue

                                write(io, "ion N $(pos[1]) $(pos[2]) 0 1 \n")
                        end
                end
                write(io, "\n\n")
                for n in -nionpos:nionpos
                        for m in -nionpos:nionpos
                                pos = n*b1+m*b2
                                thetapos = acos(dot(pos, t1)/(sqrt(dot(pos, pos))*sqrt(dot(t1, t1))))#acos(dot(pos, [1, 0, 0]))
                                thetapos < pi/3 || continue
                                thetapos > 0  || continue

                                thetapos3 = acos(dot(pos, t2)/(sqrt(dot(pos, pos))*sqrt(dot(t1, t1))))#acos(dot(pos, [1, 0, 0]))
                                thetapos3 < pi/3 || continue
                                thetapos3 > 0  || continue

                                thetapos2 = acos(dot(pos-t1-t2, -t1)/(sqrt(dot(pos-t1-t2, pos-t1-t2))*sqrt(dot(t1, t1)))) #pos-t1-t2
                                thetapos2 < pi/3 || continue
                                thetapos2 > 0  || continue

                                thetapos4 = acos(dot(pos-t1-t2, -t2)/(sqrt(dot(pos-t1-t2, pos-t1-t2))*sqrt(dot(t1, t1)))) #pos-t1-t2
                                thetapos4 < pi/3 || continue
                                thetapos4 > 0  || continue

                                write(io, "ion B $(pos[1]) $(pos[2]) 6 1 \n")
                        end
                end
                write(io, "\n\n")
                for n in -nionpos:nionpos
                        for m in -nionpos:nionpos
                                pos = 1.44/0.529*[-sin(theta), cos(theta), 0]+n*b1+m*b2
                                thetapos = acos(dot(pos, t1)/(sqrt(dot(pos, pos))*sqrt(dot(t1, t1))))#acos(dot(pos, [1, 0, 0]))
                                thetapos < pi/3 || continue
                                thetapos > 0  || continue

                                thetapos3 = acos(dot(pos, t2)/(sqrt(dot(pos, pos))*sqrt(dot(t1, t1))))#acos(dot(pos, [1, 0, 0]))
                                thetapos3 < pi/3 || continue
                                thetapos3 > 0  || continue

                                thetapos2 = acos(dot(pos-t1-t2, -t1)/(sqrt(dot(pos-t1-t2, pos-t1-t2))*sqrt(dot(t1, t1)))) #pos-t1-t2
                                thetapos2 < pi/3 || continue
                                thetapos2 > 0  || continue

                                thetapos4 = acos(dot(pos-t1-t2, -t2)/(sqrt(dot(pos-t1-t2, pos-t1-t2))*sqrt(dot(t1, t1)))) #pos-t1-t2
                                thetapos4 < pi/3 || continue
                                thetapos4 > 0  || continue

                                write(io, "ion N $(pos[1]) $(pos[2]) 6 1 \n")
                        end
                end
        end
end
main(ARGS)
