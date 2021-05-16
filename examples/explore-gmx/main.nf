#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// hard coded
initInp = 'water_nve.mdp'
initGeo = 'water.xyz'
trainInp = 'pinet.yml'
labelInp = 'water_label.mdp'

// params, will be parsed to all sub-workflows
params.publishDir = 'explore'
params.augFilter  = ''
params.resFilter  = ''
params.qbcFilter  = '-a qbc'
params.trainSeeds = '2'
params.sampleInit = 'water.xyz'
params.gmxTop     = 'water.top'
params.pinnBatch  = '1'

tipsDir = '../../nextflow'
include {sampler} from "$tipsDir/adaptor" addParams(sampler:'gromacs')
include {explore} from "$tipsDir/explore" addParams(labeller:'gromacs',
                                                    trainInp:trainInp,
                                                    labelInp:labelInp)

workflow {
    initDs = sampler([null, [inp:initInp, init:initGeo, subDir:'init']])
    inputs = initDs.map{[it[0], [initDs:it[1]]]}
    explore(inputs)
}
