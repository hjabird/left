---
layout: post
title: ! 'New Conference Paper: Leading edge vortex formation on finite wings using vortex particles'
category: posts
draft: false
---

![A visualisation of the VoFFLE method]({{ site.url }}/images/2021-01-12-New_conference_paper_aiaa_scitech_2021/LMPEVLM_iso.png "A visualisation of the VoFFLE method")

# Can I have a copy?

A pre-print is available <a href="https://www.hjab.co.uk/pdf_files/Bird2021_VoFFLE_preprint.pdf">here</a>, and the presentation for the AIAA SciTech 2021 conference is available <a href="https://youtu.be/K4W09EfRaVw">here.</a>. Related code is available at <a href="https://github.com/hjabird/LiftingLineTheory.jl/commit/1efa38267ab84675a6e157f5dca1d74205db398b">on Github.</a>

# What is the paper about?

Leading edge separation can occur in unsteady flow conditions. This leads to a leading edge vortex structure that can create significant forces on a wing. As a consequence, it is important to be able to be able to simulate this LEV structure. The problem is, unsteady finite volume method CFD can be computationally expensive to run. Plus meshes at many timesteps can be combersome to work with.

This paper provides an alternative: the Vortex Formation on Finite Leading Edge (VoFFLE) method. The VoFFLE method modelled pitch up-down ramp problem in just over a minute on a desktop, compared with CFD which required 2200 core hours on HPC resources.

# So how does VoFFLE work?

VoFFLE is a low order method - it discards important physics from a problem, to make it easier to solve. VoFFLE uses a potential flow model, discarding compressibility and viscosity. To model the wing, a vortex lattice is used. The both leading edge wake (ie. LEV) and trailing edge wakes are modelled using vortex lattices that are transformed into vortex particles after a few rows of lattice. See the picture.

The separation at the leading edge is a viscous phenomena. Since a potential flow model is being used, the separation must be modelled explicity. This is done using the critical leading edge suction model from 2D (turns out that it can be applied in 3D too). We can relate the wing lattice's leading edge filament strength to the singular leading edge vorticity distribution from 2D. With  a critical filament strength value and some linear algebra a leading edge filmant strength can be obtained. 

As for the LDVM model of Ramesh et al., the vorticity distribution and newly shed wake strengths are computed using singular distributions. For convection, everything is regularised. This helps with stuff blowing up.

# But does it work?

Kind of.

It works in some respects. It produces  a sensible looking wake, that convects nicely. An LEV of vaguely the right shape forms. It can even be shed from the wing. So, new ground has been broken in that an LEV can form and shed in a low order model. Wonderful.

However, the LEV isn't as strong as it should be, and doesn't roll up as completely as would be desired. And the simulation blows up soon after the LEV is shed (depending on the case setup). 

# What next?

Further development centres on making things work better. How to improve the LEV shedding itself isn't clear - in 2D it works, in 3D it doesn't. I don't think the code has errors, but it wouldn't be the first time I've thought something doesn't have errors. Secondly, the vortex particle method must be made more stable for longer running simulations. How to do this isn't completely clear - to some extent producing long term stable vortex particle methods seem to be a dark art known only to Winckelmans, Leonard, Cottet and others. Viscoscity might be the way forward, and can be added for only a small additional computational cost using the particle strength exchange method. Currently, my PSE method implementation fails to match the results of validation cases (it isn't nearly viscous enough). Another question is whether some kind of subgrid viscosity model is needed as for LES models. Its unclear whether using a non-uniform background grid provides benefits.


# FAQ

### What hardware is the time comparisons for?
For a test case, I claim 72s on a desktop, and 2200 core hours on HPC resources. The desktop features an AMD FX 8320E CPU and an AMD RX Vega 56 GPU. The code is far from optimal, and the n-body problem is solved naively. The 2200 cores hours is for a OpenFOAM case run on the CIRRUS HPC service, which uses 4 nodes each featuring two 18 core Broadwell CPUs. The solver parameters are described in the paper.

### Can you model tip separation?
I tried that. It seemed to have an unfortunate impact on simulation stability.

### What about non-rectangular wings?
I haven't tried this yet. Rounded wing tips are important, but there may be complications in how to relate the leading edge filament strength to the LESP criterion.


