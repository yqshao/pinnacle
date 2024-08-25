---
title: 'PiNNAcLe: Machine-Learning Potential Workflow Featuring the
Adaptive Learning-On-The-Fly Algorithm'

tags:
  - machine-learning potential
  - molecular dynamics
  - workflow

authors:
  - name: Yunqi Shao
    orcid: 0000-0002-5769-5558
    corresponding: true
    affiliation: "1, 2"
  - name: Chao Zhang
    orcid: 0000-0002-7167-0840
    affiliation: "1"

affiliations:
  - name: Department of Chemistry-Ångström Laboratory, Uppsala University, Lägerhyddsvägen 1, BOX 538, 75121 Uppsala, Sweden
    index: 1
  - name: Chalmers eCommons, Chalmers University of Technology, 41296 Gothenburg , Sweden
    index: 2
date: 2 September 2024
bibliography: paper.bib

---

# Summary

PiNNAcLe is an implementation of the so-called "activated learning"
algorithm for running machine-learning potential (MLP)-based
molecular dynamics (MD) simulations -- an emerging approach to
simulate the large-scale and long-time dynamics of systems where
empirical forms of the PES are difficult to obtain.

The algorithm aims to solve the challenge of parameterizing MLPs for
large-time-scale MD simulations, by validating simulation results at
adaptive time intervals. This approach eliminates the need of
uncertainty quantification methods for labelling new data, and thus
avoids the additional computational cost and arbitrariness thereof.

The algorithm is implemented in the NextFlow workflow language
[@2017_DiTommasoChatzouEtAl]. Components such as MD simulation and MLP
engines are designed in a modular fashion, and the workflows are
agnostic to the implementation of such modules. This makes it easy to
apply the same algorithm to different references, as well as scaling
the workflow to a variety of computational resources.

The code is published under BSD 3-Clause License, the source code and
documentation are hosted on Github. It currently supports MLP
generation with PiNN [@2020_ShaoHellstroemEtAl], reference
calculations with CP2K [@2020_KuhneIannuzziEtAl] and DFTB+
[@2020_HourahineAradiEtAl], and MD simulation with ASE
[@2017_LarsenMortensenEtAl].

# Statement of need

Recent development of MLPs see great progress in the understanding
complete descriptions of atomic configurations,
[@2021_MusilGrisafiEtAl] underpinning the feasibility of
parameterizing MLPs at high accuracy and low computational cost.  That
said, generation of MLPs is still challenging since reference data are
scarce and their distribution is unknown to the primary application of
MLPs -- to sample statistic distributions where labelling is
expensive.

A full understanding of MLPs thus calls for understanding the
parameterization procedure, the following aspects:

1. The coupling between label generation and sampling with MLPs,
   i.e. how limited initial data limit the extrapolation ability of
   MLPs, which in turn affects the data generation;
2. The most representative distribution of data given a specific
   system or dataset, i.e. how procedures like query-by-committee
   affects the performance of MLPs; [@2018_SmithNebgenEtAl;
   @2021_Zaverkin]
3. How do hyperparameters of the MLPs, such as the architecture or
   training algorithm, influence the performance of MLPs
   [@2021_ShaoDietrichEtAl]

While 2 and 3 are method-specific, 1 is a general problem for the
development of all MLPs and it precedes 2 and 3 in most cases.
Classical learn-on-the-fly (LOTF) and active learning (AL) are two
available paradigms to tackle the data-generation challenge mentioned
above.

In classical LOTF, labelling is carried out at a fixed time interval
and MLPs were used to extrapolate for 5 to 30 molecular dynamics (MD)
steps in between [@2004_CsanyiAlbaretEtAl;
@2015_LiKermodeDeVita]. This has the clear limitation that reference
calculations must be continuously carried out during the simulation.

In AL, the performance of MLPs during the extrapolative sampling is
gauged with an uncertainty estimation, through e.g. Bayesian inference
[@2019_JinnouchiKarsaiEtAl; @2020_VandermauseTorrisiEtAl],
D-optimality [@2017_PodrybinkinEvgenyEtAl], or the
query-by-committee strategy [@2017_GasteggerBehlerMarquetand;
@2020_SchranBrezinaEtAl; @2020_ZhangWangEtAl]. One can in principle
minimize the frequency of (and eventually skip) labelling by focusing
on the extrapolated data points estimated with large uncertainty.

However, the assumed correlation between the test error and the
uncertainty estimation is only verified for very short MD runs
[@2015_Behler; @2019_JinnouchiKarsaiEtAl] and not guaranteed for
either the long timescale dynamics or the out-of-distribution samples
[@2021_Zaverkin; @2018_UtevaGrahamEtAl;
@2022_VazquezSalazarLuisEtAl]. This calls for an alternative
parameterization procedure by introducing a direct feedback mechanism
to the classical LOTF without involving the uncertainty estimation,
which we name as *the adaptive LOTF*.

# Adaptive LOFT algorithm

The adaptive LOTF algorithm can be illustrated in terms of different
data and the transformation thereof, as shown in \autoref{fig1}.

![Flowchart of the adaptive LOTF algorithm, Green boxes denote
  processes, red boxes denote data which are updated at each
  iteration, and blue boxes denote loops and
  decisions.\label{fig1}](figs/flowchart.png){width=8cm}

The iterative sampling and training of NNP consists of a training
process ("PiNN Train"), a sampling process ("PiNN MD"), and a
labelling process ("CP2K Label"). Two auxiliary processes ("Converge"
and "Mixing") checks the error of the trajectories and generates a new
training set. We denote such an iteration as a "generation". In each
generation, the model, the dataset, and the training step are updated
for "PiNN Train"; the initial geometry for the extrapolative sampling
and the corresponding sampling time are updated for "PiNN MD"; the
snapshots taken from the extrapolative sampling are updated for "CP2K
Label".

Those processes and the iterative workflow are common in applications
of MLPs. What makes the adaptive LOTF protocol distinct is how the
timescale of the sampling processes are adjusted for generation
according to the convergence of test set error in the previous
generation, as highlighted in \autoref{fig2}.

![Schematics of the adaptive LOTF scheme. Compared to the typical
  feedback loop in either classical LOTF or AL, the labelling process
  in adaptive LOTF directly affects the sampling
  timescale.\label{fig2}](figs/schematics.png){width=4cm}

Specifically, the convergence of each sampled trajectory is validated
against references. If all sampled trajectories are deemed as
converged, the "PiNN MD" sampling time ("time") is multiplied by a
factor for the next generation, starting from a given small timescale
to generation to a maximal timescale where the sampling is deemed as
sufficient. Otherwise, "PiNN Train" in the next generation to
incorporate new reference data.

# Acknowledgements

This work has been supported by the Swedish Research Council (VR),
grant 2019-05012. The authors are thankful for the funding from the
Swedish National Strategic e-Science program eSSENCE, STandUP for
Energy and BASE (Batteries Sweden). The simulations were performed
on the resources provided by the Swedish National Infrastructure for
Computing (SNIC) at C3SE and PDC.

We thank A. Laio for reading the manuscript and for providing helpful
feedback. Y.S. also thanks A. Laio for hosting a research visit at the
International School for Advanced Studies (SISSA) and D. Doimo for
many useful discussions.

# References
