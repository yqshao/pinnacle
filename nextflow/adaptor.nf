#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.trainer  = 'pinn'
params.sampler  = 'pinn'
params.labeller = 'lammps'
params.filter   = 'tips'


// All known implementations
switch (params.filter) {
    case {it.startsWith('./')}:
        include {filter} from params.filter;
        break;
    case 'tips':
        include {tipsFilter as filter} from "$moduleDir/tips";
        break;
    default:
        throw new Exception("Unkown filter $params.filter.");
}

switch (params.trainer) {
    case {it.startsWith('./')}:
        include {trainer} from params.trainer;
        break;
    case 'pinn':
        include {pinnTrain as trainer} from "$moduleDir/pinn";
        break;
    default:
        throw new Exception("Unkown trainer $params.trainer.");
}

switch (params.sampler) {
    case {it.startsWith('./')}:
        include {sampler} from params.sampler;
        break;
    case 'pinn':
        include {pinnSample as sampler} from "$moduleDir/pinn";
        break;
    case 'lammps':
        include {lammpsSample as sampler} from "$moduleDir/lammps";
        break;
    case 'gromacs':
        include {gromacsSample as sampler} from "$moduleDir/gromacs";
        break;
    default:
        throw new Exception("Unkown sampler $params.sampler.");
}

switch (params.labeller) {
    case {it.startsWith('./')}:
        include {labeller} from params.labeller;
        break;
    case 'pinn':
        include {pinnLabel as labeller} from "$moduleDir/pinn";
        break;
    case 'lammps':
        include {lammpsLabel as labeller} from "$moduleDir/lammps";
        break;
    case 'gromacs':
        include {gromacsLabel as labeller} from "$moduleDir/gromacs";
        break;
    default:
        throw new Exception("Unkown labeller $params.labeller.");
}
