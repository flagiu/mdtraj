# time is in units of dt

# base of the logarithm
alpha: float = 1.25
# number of points per cycle (logarithmic scale)
npc: int =       36
# number of cycles (linear scale)
ncycles: int =   10
# spacing (in timesteps) btw. consecutive cycles (Better if multiple of alpha!)
delta: float =  5.0

file="schedule.times"

#---------------------------------------------------------------#

t0=int(delta/alpha)
period = t0*(alpha**npc - alpha)

print()
print("<---------------- SUMMARY -------------------->")
print("Log base: %.2f"%alpha)
print("Points per cycle: %d"%npc)
print("Spacing btw consecutive cycles: %.2f"%delta)
print("--> Log resolution: %.2f ~ %.0f"%(t0*(alpha-1), round(t0*(alpha-1))) )
print("--> Log capacity  : %.2f ~ %.0f"%(t0*alpha**(npc-1), round(t0*alpha**(npc-1))) )
print("--> Lin resolution (timesteps): %.2f ~ %.0f"%(period,round(period)) )
print("Number of cycles: %d"%ncycles)
print("--> Lin capacity   (timesteps): %.2f ~ %.0f"%( ncycles*period, round(ncycles*period)) )
print("<--------------------------------------------->")

f = open(file,"w")
for i in range(ncycles):
	for j in range(npc):
		#t = round(t0*(alpha**npc)*i + t0*(alpha**(j+1)) ) # wrong!
		t = round(t0 * alpha**npc * i) + round(t0 * alpha**(j+1))
		f.write("%d\n"%t)
print("Maximum step for LAMMPS: %d"%t)
print("<--------------------------------------------->")
print("Times saved into %s"%file)
print("<--------------------------------------------->")
print()
# add an extra point for safety in LAMMPS
i=ncycles
j=0
t = t0 * (alpha**npc * i + alpha**(j+1) )
f.write("%d\n"%round(t))
f.close()
