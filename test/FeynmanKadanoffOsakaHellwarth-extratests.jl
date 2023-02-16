# Hellwarth 1999 PRB - Part IV; T-dep of the Feynman variation parameter
# A Friday afternoon of hacking to try and implement the T-dep electron-phonon coupling from the above PRB
# Which was unusually successful! And more or less reproduced Table III

"Print out F(alpha,beta) for a specific v,w; as a test"
function test_fns()
    @printf("\t\t")
    for α in 1:5
        @printf("α=%d\t\t",α)
    end
    @printf("\n")

    for β in 1:0.25:3.0
        v=w=4
        print("β: $β  \t||")
        for α in 1:5
            @printf("%f\t",F(v, w, α, 1, β))
        end
        println()
    end
end

test_fns()  # OK - very primitive!

"
These are 1D traces along the solution for Alpha=Beta=1 in Helwarth PRB TABLE III,
this was used to correct a transcription error in the above typed-in equations
It was also good to see what F(v,w) looked like as a function of v and w near an optimal solution"
function test_trace()
	v=7.20
	w=6.5
	α=1.0
	β=1.0

	for v=6:0.1:8
    	@printf("%f %f\n",v,F(v, w, α, 1, β))
	end

	@printf("\n")
	v=7.20
	for w=6:0.1:7
    	@printf("%f %f\n",w,F(v, w, α, 1, β))
	end
end
test_trace()

# Angle for the ringside seats, when the fall, don't blame me, Bring on the Major Leagues
Fopt(x) = F(x[1],x[2],1,1)

function test_Fopt()
    show(Fopt([7.2,6.5]))
	# OK! It looks like I can bury the alpha, beta parameters (which we don't optimise), by wrapping our function in a function definition.
    initial=[7.2,6.5]

    show(optimize(Fopt,  initial, LBFGS()))
    show(optimize(Fopt, initial, BFGS(), Optim.Options(autodiff=true)))
end

test_Fopt()

function test_Optim()
	# After a bit of fiddling, I figured out how to add bounds, to stop that 'DomainError', 
	# which occurs where the you are evaluating log(-ve Real), i.e. w<0.0 or v<0.0

	initial=[7.2,6.5]

	lower=[0.0,0.0]
	upper=[10.0,10.0]

	@printf("\t\t")
	for α in 1:5
		@printf("α=%d\t\t",α)
	end
	@printf("\n")

	for β in 1:0.25:3.0
		print("β: $β  \t||")
		for α in 1:5
			myf(x) = F(x[1], x[2], α, 1, β)
			solution=optimize(DifferentiableFunction(myf), initial, lower, upper, Fminbox(); optimizer = ConjugateGradient, optimizer_o=Optim.Options(autodiff=true))
			minimum=Optim.minimizer(solution)

			v=minimum[1]
			w=minimum[2]
			#print(solution,"\t")
			@printf("%.2f %.2f\t",v,w)


		end
		println()
	end 
end

test_Optim()

function test_Optimisers()
	# So that looks really good! I was super stoked to see how close these values are to TABLE III in Hellwarth
	# However, the solutions all start on (7.20,6.50) so that top-left data point is cheating, whereas the 
	# others have some disagreement / noise associated with them
	# I was wondering whether it might be a function of the optimiser, so thought I'd try them all

	initial=[7.1,6.5]
	# Main use of these bounds is stopping v or w going negative, at which you get a NaN error as you are evaluating log(-ve Real)
	lower=[1.0,1.0]
	upper=[10.0,10.0]

	for optimizer in [BFGS, LBFGS, ConjugateGradient] # Newton, GradientDescent, NelderMead - steps outside box & log(-ve)->NaN error
		@printf("\n\t\t##### NOW TRIALING: %s #####\n\n",optimizer)

		@printf("\t\t")
		for α in 1:5
			@printf("α=%d\t\t",α)
		end
		@printf("\n")

		for β in 1:0.25:3.0
			print("β: $β  \t||")
			for α in 1:5
				myf(x) = F(x[1], x[2], α, 1, β)
				res=optimize(DifferentiableFunction(myf), initial, lower, upper, Fminbox(); optimizer = optimizer, optimizer_o=Optim.Options(autodiff=true))
				minimum=Optim.minimizer(res)
				#show(Optim.converged(res)) # All came out as 'true'

				#print(solution,"\t")
				@printf("%.2f %.2f\t",minimum[1],minimum[2])
			end
			println()
		end
	end
end

test_Optimisers()

