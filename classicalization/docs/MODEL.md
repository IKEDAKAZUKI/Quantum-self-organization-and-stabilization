# Autonomous environmental conversion and persistence of spatial order

This study extends the bosonic dimer lattice of **Quantum entanglement as a resource for self-organization and stabilization**. It introduces explicit environmental degrees of freedom.

## Fixed enlarged-system generator

Each cell has an environmental two-state register, with projectors H_x = |1><1| and 1-H_x. State 1 is an irreversible converted state. The registers begin in state 0. No parameter is changed during a run.

Use the physical occupation contrast u_x, intermode coherence operator C_x, and original finite-range field drift. Replace C_x in the transport term by

    Q_x = (1 - H_x) C_x + chi H_x.

The transport contribution on edge (x,y) is j_xy (Q_x + Q_y)(u_x-u_y)/2. The original, weak coherence-independent coefficient remains Dcl = 0.001 throughout. The additional classical transport appears when H changes, rather than through an external schedule for Dcl. The coherence-generating jump becomes

    sqrt(r) (1 - H_x) R_x.

The original number dephasing, leakage, repair, background noise, nonlinear saturation, and stabilizing fourth-order spatial term are retained. Field jumps are T_x^± sqrt(q_x^±) with the same positive-rate formula as the original model, using the enlarged drift operator.

## Environmental response

Define the local occupation activity on a site and its four nearest neighbors,

    A_x = [u_x^2 + sum_{y nearest x} u_y^2] / 5.

The activation operator is

    Lambda_x = Gamma A_x^4 / (A_x^4 + A_*^4)
               + (beta/4) sum_{y nearest x} H_y.

For each physical number state a of dimer x, add the Lindblad jump

    K_{x,a} = |1><0|_{env,x} P_{x,a} sqrt(Lambda_x).

Here P_{x,a} is the physical number-state projector, not a projector onto an entangled relative-phase basis state. The rate is therefore proportional to the local state operator, not to an externally evaluated expectation value. The enlarged master equation is linear in its density operator, completely positive, and time homogeneous.

The activation channel preserves the instantaneous occupation distribution. It records which physical number state was present and removes that cell's intermode coherence. Its environment register then disables further local coherence pumping and enables coherence-independent transport. Neighbor-assisted conversion allows low-contrast cells to convert as well. There is no prescribed conversion time, stored pattern template, global pattern classifier, or postselection in the generator.

The construction is an effective model of an active, irreversible environment. It is not derived from a particular passive thermal bath. The source, activated transport, and irreversible environmental record require a nonequilibrium physical implementation. A particular microscopic experimental implementation and its energetic costs are not established by these simulations.

## Parameters

All original parameters are unchanged: nx=48, ny=56, eps=0.01, U=1, delta=0.3, D4=1, b=1, c=0.1, r=20, eta=1, leak=0.0002, repair=2, Db=0.00001, smooth=0.002, Dcl=0.001, and disorder=0.05. The three original families use (D2,nu,zeta) = (1.6,0,0.8), (1.4,0.55,0), and (1.8,0,0).

Added baseline parameters:

| Parameter | Value | Meaning |
|---|---:|---|
| Gamma | 0.002 | Local activity-induced conversion rate scale |
| beta | 0.08 | Neighbor-assisted conversion rate scale |
| A_* | 0.10 | Squared-contrast activity scale |
| chi | 0.91 | Gain of the converted environment |
| Hill exponent | 4 | Cooperative dependence on local squared contrast |

All coefficients are constant. In the fully converted, homogeneous limit, the effective classical gain is Dcl + chi D2: 1.457, 1.275, and 1.639. These exceed 2 sqrt(delta D4) = 1.095445..., the continuous-wavevector instability threshold. Static edge disorder is retained.

## Exact stochastic representation

The invariant density operators are mixtures in a hybrid basis. Unconverted nonleaked cells use |n,+> and |n,->. Converted cells use mixtures of physical number states |n,L> and |n,R>. Leakage uses |n,0> in both sectors. The field jumps are diagonal-rate shifts in this basis. Coherence pumping and dephasing generate the original phase-label transitions only on unconverted cells. Environmental conversion maps each |n,+> or |n,-> population to an equal incoherent mixture of |n,L> and |n,R>. Their labels can be lumped because subsequent rates depend on total occupation and leakage, not on L versus R.

The Gillespie simulation is an exact continuous-time sampling of this invariant population process. It does not use a deterministic mean-field closure, time-step discretization, or continuous external projection in the primary runs.

## Entanglement and classicality

With spatially averaged unconditional expectations e = M^{-1} sum_x <C_x> and ell = M^{-1} sum_x <P0,x>, the original witness gives a lower bound on the mean local negativity,

    max(e - 2 ell, 0)/4 <= mean_x N(rho_x).

Convexity of negativity in the hybrid-basis mixture gives an independent upper bound,

    mean_x N(rho_x) <= [sum_x Pr(H_x=0, nonleaked)]/(2 M).

The two bounds are evaluated from ensemble means. They are not inferred from entanglement of individual sampled |n,+> states alone.

Once every environmental register is converted, the physical state in that environmental sector is exactly number-diagonal and remains so. This is an invariant separable sector, established by the jump operators. Its classicality does not depend on the failure of an entanglement witness.

Finite stochastic samples all entering this sector do not prove that the exact unconditioned density operator has exactly zero negativity at a fixed finite time. Rare unconverted branches can have nonzero weight. The upper bound tends to zero as their weight vanishes. On a finite connected lattice with positive background noise and the stated irreversible conversion, complete conversion occurs almost surely. This supports asymptotically exact classicalization and finite-time entry of individual environmental records into the classical sector. It is not a demonstration of a thermodynamic quantum phase transition, critical exponents, or a nonanalytic stationary-state boundary.

## Controls

- `no_quantum`: remove coherence pumping and the unconverted coherence-dependent field channel; retain the environmental response and converted classical transport. This is a channel-removal control, not the exact projected generator.
- `projected_control`: impose physical-number projection for comparison only. Nonlinear field rates are averaged over all unconverted relative-phase signs. It is an independent, exact number-projected control and retains the corresponding noise.
- `no_classical_gain`: keep the environmental dephasing and pump shutdown, but set chi=0. This tests whether loss of quantum coherence alone preserves order.
- `no_environment`: set Gamma=beta=0. Quantum formation and maintenance remain active.
- `initial_coherent`: initialize the local relative-phase mixture with e(0)=0.91 instead of starting with e(0)=0. All subsequent dynamics remains autonomous.

The baseline main runs start with the original weak, number-diagonal field. The quantum source develops local entanglement before environmental conversion and macroscopic ordering. The initial-coherence runs additionally test a literally entangled initial stage.

Quantum-free controls not activating during the observation window do not prove that activation is impossible at all times. Rare noise-induced nucleation remains possible. The comparison is over the stated finite interval.

## Measurements

The Fourier band, morphology classifier, and periodic spatial-extent criteria are inherited from the supplied study and are held fixed. Morphology is checked at a specified set of saved times; spectral power is computed at every saved time. All realizations, including mixed morphologies and changes of morphology, enter the statistics. Correlation is the zero-shift Pearson correlation with the stated reference configuration. It is not optimized over translations or rotations.

The persistence reference t_* is the first common saved time at which all 48 baseline autonomous environmental records are fully converted. The latest individual conversion event occurs at t=340.76007873663474, so the shared observation grid gives t_*=345. The generator is unchanged at this time. Final/reference ratios average per-run ratios, not a ratio of ensemble means. Reported uncertainty is one standard error across independent initial fields and disorder realizations within each family. Corresponding control runs share their initial fields and disorder and are paired comparisons.

## Relation to the original separable bounds

The enlarged model has additional environmental channels and a converted transport sector. The original fixed-channel bounds cannot be applied unchanged to this enlarged system. The result does not contradict the original paper's conclusion that entanglement is necessary for the original fixed-parameter dynamics.

## Context

Time-independent local dissipation as an engineering resource is discussed by Verstraete, Wolf, and Cirac, Nature Physics 5, 633-636 (2009), DOI 10.1038/nphys1342. Absorbing-state open quantum dynamics is studied, for example, by Carollo et al., Physical Review Letters 123, 100604 (2019), DOI 10.1103/PhysRevLett.123.100604. These works motivate the general construction method; they do not establish this particular environmental model or its simulation results.
