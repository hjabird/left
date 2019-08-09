---
layout: post
title: Redistributing vortex particles
category: posts
draft: true
---

	julia> using LiftingLineTheory
	[ Info: Recompiling stale cache file C:\Users\hjabi\.julia\compiled\v1.1\LiftingLineTheory.ji for LiftingLineTheory [top-level]

	julia> using CVortex

	julia> prob = LAUTAT(;kinematics=RigidKinematics2D(x->0, x->deg2rad(25) * sin(x), 0.0), lesp_critical=0.11, num_fourier_terms=16);

	julia> for i = 1 : 250
			   println(i)
			   advance_one_step(prob)
			   if i % 10 == 0
				   ppos = vcat(prob.te_particles.positions[1:end-1,:], prob.le_particles.positions[1:end-1,:])
				   pvort = vcat(prob.te_particles.vorts[1:end-1], prob.le_particles.vorts[1:end-1])
				   ppos, pvort, ~ = CVortex.redistribute_particles_on_grid(ppos, pvort, lambda3_redistribution(), 0.025)
				   prob.te_particles.positions = vcat(ppos, prob.te_particles.positions[end,:]')
				   prob.te_particles.vorts = vcat(pvort, prob.te_particles.vorts[end])
				   if length(prob.le_particles.vorts) > 0
					   prob.le_particles.positions = prob.le_particles.positions[end,:]'
					   prob.le_particles.vorts = [prob.le_particles.vorts[end]]
				   end
			   end
			   to_vtk(prob, "hhabb_"*string(i))
		   end
