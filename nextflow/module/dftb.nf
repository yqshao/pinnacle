// An Alias for ASE-based calculators for DFTB+

nextflow.enable.dsl=2

params.publish = 'dftb'
params.dftb_aux = '/dev/null'

include { sp; md } from './ase.nf' addParams (
  ase_label: 'dftb',
  ase_aux: params.dftb_aux
)
