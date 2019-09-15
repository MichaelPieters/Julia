## AUTHOR: Pieters Michaël
## DATE: 15 September 2019


using Plots

###############################################################
##########CONSTANTS############################################
###############################################################

const global NDIM = 1
const global gamma = 1.4



###############################################################
##########STRUCT for PARTICLE##################################
###############################################################

mutable struct Particle
	x
	vx
	num
	mass
	smooth
	rho
	pressure
	utherm
	numngb
	numpairs
	soundspeed
	acceleration
	denergy
end
	
###############################################################
##########FUNCTIONS############################################
###############################################################

function norm(h)
	if NDIM == 3
		return 1.0/(pi*h*h*h)
	elseif NDIM == 2
		return 10.0/(7.0*pi*h*h)
	else
		return 2.0/(3.0*h)
	end
end
		
function smoothing_kernel(r, h)
    u = r/h
    if (u >= 0.0) && (u < 1.0)
        return norm(h)*(1.0 - (3.0/2.0)*u*u + (3.0/4.0)*u*u*u)
    elseif (u >= 1.0) && (u < 2.0)
        return norm(h)*(1.0/4.0)*(2.0 - u)^3
    else
        return 0.0
	end
end

function dsmoothing_kernel(r, h)
    u = r/h
    if (u >= 0.0) && (u < 1.0)
        return (norm(h)/h)*3.0*u*((3.0/4.0)*u - 1.0)
    elseif (u >= 1.0) && (u < 2.0)
        return -(norm(h)/h)*(3.0/4.0)*(2.0 - u)^2
    else
        return 0.0
	end
end

function artificial_viscosity(a, b)
    alfa = 1.0
    beta = 1.0

    dotproduct = (particles[a].x - particles[b].x)*(particles[a].vx - particles[b].vx)
    innerproduct = (particles[a].x - particles[b].x)*(particles[a].x - particles[b].x)
    
    monaghansoundspeed = (particles[a].soundspeed + particles[b].soundspeed)/2.0
    monaghandensity = (particles[a].rho + particles[b].rho)/2.0
    monaghansmooth = (particles[a].smooth + particles[b].smooth)/2.0
    phi = 0.1*monaghansmooth
    omega = monaghansmooth*dotproduct/(innerproduct + phi*phi)

    if dotproduct < 0.0
        return omega*(beta*omega - alfa*monaghansoundspeed)/monaghandensity
    else
        return 0.0
	end
end

###############################################################
##########DATA EXTRACTING######################################
###############################################################


function extract_data(filename)
	ROW = []
	
	open(filename) do file
		readline(file)
		readline(file)
		
		while !eof(file)
			line = readline(file)
			words = split(line)
		
			x = parse(Float64, words[1])
			vx = parse(Float64, words[2])
			num = parse(Int64, words[3])
			mass = parse(Float64, words[4])
			smooth = parse(Float64, words[5])
			rho = parse(Float64, words[6])
			pressure = parse(Float64, words[7])
			utherm = parse(Float64, words[8])
			numngb = parse(Float64, words[9])
			
			push!(ROW, Particle(x, vx, num, mass, smooth, rho, pressure, utherm, numngb, 0, 0, 0, 0))
		end
	end
	
	return ROW
end

particles = extract_data("shock-tube.dat")
number = length(particles)
println(number)

#for j in 1:number
#    println(particles[j].num, particles[j].numngb)
#end

###############################################################
##########WORK HORSE###########################################
###############################################################

function get_pairs()
    pairnumbers = []

    # resetting the neighbour count
    for k in 1:number
        particles[k].numngb = 0.0
	end
        
    for i in 1:number-1
        for j in i+1:number
            xi = particles[i].x
            xj = particles[j].x
            hi = particles[i].smooth
            hj = particles[j].smooth
            if (abs(xj - xi) < 2*max(hi, hj))
				push!(pairnumbers, [particles[i].num+1, particles[j].num+1])
                particles[i].numngb += 1
                particles[j].numngb += 1
			end
		end
	end
	
    return pairnumbers
end

function get_densities(PAIRS)

    for k in 1:number
        particles[k].rho = particles[k].mass*smoothing_kernel(0.0, particles[k].smooth)
	end
        
    for l in 1:length(PAIRS)
        particles[PAIRS[l][1]].rho += particles[PAIRS[l][2]].mass*smoothing_kernel(abs(particles[PAIRS[l][2]].x - particles[PAIRS[l][1]].x), particles[PAIRS[l][1]].smooth)
        particles[PAIRS[l][2]].rho += particles[PAIRS[l][1]].mass*smoothing_kernel(abs(particles[PAIRS[l][1]].x - particles[PAIRS[l][2]].x), particles[PAIRS[l][2]].smooth)
	end

    for k in 1:number
        particles[k].pressure = (gamma - 1.0)*(particles[k].rho)*(particles[k].utherm)
        particles[k].soundspeed = sqrt((gamma - 1.0)*(particles[k].utherm))
	end

end
	
function get_energy()

    ENERGY = 0
    for k in 1:number
        ENERGY += particles[k].mass*(particles[k].utherm + particles[k].vx^2/2.0)
	end
	
	return ENERGY
end

function compute_accelerations(PAIRS)

    for k in 1:number
        particles[k].acceleration = 0.0
	end

    for l in 1:length(PAIRS)
        particles[PAIRS[l][1]].acceleration += -particles[PAIRS[l][2]].mass*(particles[PAIRS[l][1]].pressure/(particles[PAIRS[l][1]].rho)^2 +
                                                                             particles[PAIRS[l][2]].pressure/(particles[PAIRS[l][2]].rho)^2 +
                                                                             artificial_viscosity(PAIRS[l][1], PAIRS[l][2]))*
																			 (particles[PAIRS[l][1]].x - particles[PAIRS[l][2]].x)/(abs(particles[PAIRS[l][1]].x - particles[PAIRS[l][2]].x))*
																			 dsmoothing_kernel(abs(particles[PAIRS[l][1]].x - particles[PAIRS[l][2]].x), particles[PAIRS[l][1]].smooth)
        particles[PAIRS[l][2]].acceleration += -particles[PAIRS[l][1]].mass*(particles[PAIRS[l][1]].pressure/(particles[PAIRS[l][1]].rho)^2 + 
                                                                             particles[PAIRS[l][2]].pressure/(particles[PAIRS[l][2]].rho)^2 + 
                                                                             artificial_viscosity(PAIRS[l][2], PAIRS[l][1]))*
																			 (particles[PAIRS[l][2]].x - particles[PAIRS[l][1]].x)/(abs(particles[PAIRS[l][2]].x - particles[PAIRS[l][1]].x))*
																			 dsmoothing_kernel(abs(particles[PAIRS[l][2]].x - particles[PAIRS[l][1]].x), particles[PAIRS[l][2]].smooth)
	end
end


function compute_denergy(PAIRS)

    for k in 1:number
        particles[k].denergy = 0.0
	end

    for l in 1:length(PAIRS)
        particles[PAIRS[l][1]].denergy += (1.0/2.0)*particles[PAIRS[l][2]].mass*(particles[PAIRS[l][1]].pressure/(particles[PAIRS[l][1]].rho^2) +
                                                                             particles[PAIRS[l][2]].pressure/(particles[PAIRS[l][2]].rho^2) +
                                                                             artificial_viscosity(PAIRS[l][1], PAIRS[l][2]))*
																			 (particles[PAIRS[l][1]].vx - particles[PAIRS[l][2]].vx)*
																			 (particles[PAIRS[l][1]].x - particles[PAIRS[l][2]].x)/(abs(particles[PAIRS[l][1]].x - particles[PAIRS[l][2]].x))*
																			 dsmoothing_kernel(abs(particles[PAIRS[l][2]].x - particles[PAIRS[l][1]].x), particles[PAIRS[l][1]].smooth)
        particles[PAIRS[l][2]].denergy += (1.0/2.0)*particles[PAIRS[l][1]].mass*(particles[PAIRS[l][1]].pressure/(particles[PAIRS[l][1]].rho^2) +
                                                                             particles[PAIRS[l][2]].pressure/(particles[PAIRS[l][2]].rho^2) +
                                                                             artificial_viscosity(PAIRS[l][2], PAIRS[l][1]))*
																			 (particles[PAIRS[l][2]].vx - particles[PAIRS[l][1]].vx)*
																			 (particles[PAIRS[l][2]].x - particles[PAIRS[l][1]].x)/(abs(particles[PAIRS[l][2]].x - particles[PAIRS[l][1]].x))*
																			 dsmoothing_kernel(abs(particles[PAIRS[l][1]].x - particles[PAIRS[l][2]].x), particles[PAIRS[l][2]].smooth)
	end
end

																		 
###############################################################
##########ALGORITHM############################################
###############################################################

# 1) Determine the Pairs
# 2) Calculate the initial densities
# 3) Get the accelerations
# 4) Update v to halfway the timestep
# 5) Get the changes in energy
# 6) Update the energy
# 7) Update v to the end of the time step
# 8) Now the velocities are updated, update the positions

dt = 0.005

# first inital step
PAIRS = get_pairs()
#println(PAIRS)
println(PAIRS[1][1])
get_densities(PAIRS)
compute_accelerations(PAIRS)



for k in 1:number
    particles[k].x += particles[k].vx*(dt/2.0) + particles[k].acceleration*(dt*dt/8.0)
end


function update()
    PAIRS = get_pairs()
    get_densities(PAIRS)
    compute_accelerations(PAIRS)
	
    for k in 1:number
        particles[k].vx += particles[k].acceleration*(dt/2.0)
	end
	
    compute_denergy(PAIRS)
	
    for k in 1:number
        particles[k].utherm += particles[k].denergy*dt
	end
	
    for k in 1:number
        particles[k].vx += particles[k].acceleration*(dt/2.0)
        particles[k].x += particles[k].vx*dt
	end
end
    
    
for countdown in 1:40
    update()
    println(countdown)
end

###############################################################
##########PLOTTING#############################################
###############################################################

X = []
V = []
D = []
P = []

for k in 1:number
	push!(X, particles[k].x)
	push!(V, particles[k].vx)
	push!(D, particles[k].rho)
	push!(P, particles[k].pressure)
end

plot(X, V, xlim=(-0.4, 0.4), ylim=(0.0, 0.8), xlabel="x (m)", title="Sod Shock tube Test", label="Velocity (m/s)")
plot!(X, D, xlim=(-0.4, 0.4), ylim=(0.0, 1.2), xlabel="x (m)", title="Sod Shock tube Test", label="Density (kg/m^3)")
plot!(X, P, xlim=(-0.4, 0.4), ylim=(0.0, 1.2), xlabel="x (m)", title="Sod Shock tube Test", label="Pressure (N/m^2)")