using LinearAlgebra
a = 1.44*sqrt(3)/0.529
a1 = a*[1, 0, 0];
a2 = a*[1/2, sqrt(3)/2, 0];
theta = acos(13/14)

b1 = a*[cos(theta), sin(theta), 0 ]
b2 = a*[cos(pi/3+theta), sin(pi/3+theta), 0]


t1 = 2*a2+a1;
t2 = 3*a2-2*a1

open("Moire.lattice", "w") do io
	write(io, "lattice\\ \n")
	write(io, "$(t1[1]) $(t2[1]) 0 \\ \n")
        write(io, "$(t1[2]) $(t2[2]) 0 \\ \n")
	write(io, "0 0 20\n")
end


open("Moire.ionpos", "w") do io
	write(io, "ion B 0 0 0 1\n")
        write(io, "ion B 0 0 6 1\n")
        for n in -5:5
                for m in -5:5
                        pos = n*a1+m*a2
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

                        write(io, "ion B $(pos[1]) $(pos[2]) 0 1 \n")

		end
	end

	write(io, "\n\n")

        for n in -5:5
                for m in -5:5
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

        for n in -10:10
                for m in -10:10
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

        for n in -10:10
                for m in -10:10
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



#=
open("Moire.ionpos", "w") do io

	for n in -5:5
		for m in -5:5
	
			pos = n*a1+m*a2
			#println(pos)
			apos = sqrt(dot(pos, pos))
			at1 = sqrt(dot(t1, t1))
			at2 = sqrt(dot(t2, t2))
			#println("at2= $(at2)")
                        #println("at2= $(at2)")
			#println(dot(pos, t1)/at1)
			#println(dot(pos, t2)/at2)
			#println(dot(pos, t1))
			0 <= dot(pos, t1)/(at1^2) < 1 || continue
			0 <= dot(pos, t2)/(at2^2) < 1 || continue

			#(pos[1]*t1[1]+pos[2]*t1[2])/(7*sqrt(n^2+m^2+n*m)*a*a) < 1 || continue 
                        #(pos[1]*t1[1]+pos[2]*t1[2])/(7*sqrt(n^2+m^2+n*m)*a*a) >= 0  || continue
                        #(pos[1]*t2[1]+pos[2]*t2[2])/(7*sqrt(n^2+m^2+n*m)*a*a) < 1 || continue
                        #(pos[1]*t2[1]+pos[2]*t2[2])/(7*sqrt(n^2+m^2+n*m)*a*a) >= 0  || continue
			write(io, "ion B $(pos[1]) $(pos[2]) 0 1 \n")
		end
	end

end
=#
