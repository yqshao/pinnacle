from ase.calculators.dftb import Dftb

calc = Dftb(
    Hamiltonian_ = "DFTB",
    Hamiltonian_Scc = "Yes",
    Hamiltonian_ReadInitialCharges = "Yes",
    Hamiltonian_ShellResolvedSCC = "Yes",
    Hamiltonian_MaxAngularMomentum_='',
    Hamiltonian_MaxAngularMomentum_Au='d',
    Hamiltonian_MaxSCCIterations='1000',
    Hamiltonian_SCCTolerance='1e-06',
    Hamiltonian_SlaterKosterFiles_='Type2FileNames',
    Hamiltonian_SlaterKosterFiles_Type2FileNames_Prefix='sktable/',
    Hamiltonian_SlaterKosterFiles_Type2FileNames_Separator='-',
    Hamiltonian_SlaterKosterFiles_Type2FileNames_Suffix='.skf',
    Hamiltonian_Filling_='MethfesselPaxton',
    Hamiltonian_Filling_Order='2',
    Hamiltonian_Filling_Temperature='300.0',
    kpts=(1,1,1),
    Parallel_ = "",
    Parallel_UseOmpThreads = "Yes",
)
