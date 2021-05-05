#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// hard coded
initInp = 'csvr.lmp'
initGeo = 'init.xyz'
trainInp = 'pinet.yml'
labelInp = 'label.lmp'

// params, will be parsed to all sub-workflows
params.publishDir = 'explore'
params.augFilter  = '-vmax "e:-20" -amax "f:10"'
params.resFilter  = '-vmax "e:-42" -amax "f:4"'
params.qbcFilter  = '-a qbc'
params.trainSeeds = '2'
params.sampleInit = 'init.xyz'
// lammps specific params
params.lmpEmap    = '1:1,2:8,3:11,4:17'
params.lmpInit    = 'nacl.init'
params.lmpData    = 'nacl.data'
params.lmpSetting = 'nacl.setting'
// pinn specific params
params.pinnBatch  = '1'
params.repeat = 3

tipsDir = '../../nextflow'
include {sampler} from "$tipsDir/adaptor" addParams(sampler:'lammps')
include {explore} from "$tipsDir/explore" addParams(trainInp:trainInp,
                                                    labelInp:labelInp)

workflow {
    initDs = sampler([null, [inp:initInp, init:initGeo, subDir:'init']])
    inputs = initDs.flatMap{(1..params.repeat)
                            .collect{run->[it[0], [subDir:"explore$run",initDs:it[1]]]}}
    explore(inputs)
}
