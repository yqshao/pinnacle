from ase.calculators.dftb import Dftb

calc = Dftb(
    Hamiltonian_ = "xTB",
    Hamiltonian_Method = "GFN2-xTB",
    Hamiltonian_MaxAngularMomentum_='',
    Hamiltonian_MaxAngularMomentum_H='s',
    kpts=(1,1,1),
    Parallel_ = "",
    Parallel_UseOmpThreads = "Yes",
)
