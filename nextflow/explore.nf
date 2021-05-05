#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// parameters trainer/sampler/labeller can be replaced with processes
params.publishDir     = 'explore'
params.trainer        = 'pinn'
params.sampler        = 'pinn'
params.labeller       = 'lammps'
params.maxIter        = '20'
params.resFilter      = '-amax "f:4"'
params.augFilter      = '-amax "f:8"'
params.qbcFilter      = '-a qbc'
// other parameters can be are adjusted with the inputs channel, defaults listed below
defaults = [:]
defaults.subDir       = '.'
defaults.initDs       = null
defaults.trainInp     = null
defaults.trainSeeds   = '3'
defaults.trainSteps   = '200000'
defaults.retrainSteps = '200000'
defaults.sampleInit   = null
defaults.sampleParams = [:]
defaults.labelInp     = 'label.lmp'

include {iterUntil; setNext; getParams} from "$moduleDir/utils"
include {trainer; sampler; labeller} from "$moduleDir/adaptor"
include {labeller as qbcLabel} from "$moduleDir/adaptor" addParams(labeller:params.trainer)
include {filter as qbcFilter} from "$moduleDir/adaptor" addParams(inp:params.qbcFilter)
include {filter as augFilter} from "$moduleDir/adaptor" addParams(inp:params.augFilter)
include {filter as resFilter} from "$moduleDir/adaptor" addParams(inp:params.resFilter)

workflow explore {
    take:
    inputs

    main:
    condition = {it[0].iter>=params.maxIter.toInteger()}
    defaults = getParams(defaults, params)
    setup = inputs.map{[it[1]+[meta:it[0]], getParams(defaults, it[1])]}
        .map{[it[0]+[seeds:it[1].trainSeeds.toInteger(), subDir:it[1].subDir], it[1]]}
        .flatMap{(1..it[0].seeds).collect(seed->[it[0]+[seed:seed], it[1]])}

    // initalize the first iteration, only sample and label the first seed
    initDs      = setup.map{[it[0], [ds: it[1].initDs]]}
    initCkpt    = setup.map{[it[0], [inp: it[1].trainInp]]}
    initSteps   = setup.map{[it[0], [maxSteps: it[1].trainSteps.toInteger(), retrainSteps: it[1].retrainSteps.toInteger()]]}
    initSample  = setup.filter{it[0].seed==1}.map{[it[0], [init: it[1].sampleInit]+it[1].sampleParams]}
    initLabel   = setup.filter{it[0].seed==1}.map{[it[0], [inp: it[1].labelInp]]}

    // // create the iterating channels
    (trainDs,      nextDs    ) = iterUntil(initDs,     condition)
    (trainCkpt,    nextCkpt  ) = iterUntil(initCkpt,   condition)
    (trainSteps,   nextSteps ) = iterUntil(initSteps,  condition)
    (sampleParams, nextSample) = iterUntil(initSample, condition)
    (labelInp,     nextLabel ) = iterUntil(initLabel,  condition)

    // // the actual work for each iteraction
    trainInp   = trainDs.join(trainCkpt).join(trainSteps).map{[it[0], [ds:it[1].ds]+it[2]+[genDress:it[0].iter==1, maxSteps:it[3].maxSteps]]}
    models     = trainer(trainInp.map{[it[0],it[1]+[subDir:"${it[0].subDir}/models/iter${it[0].iter}/seed${it[0].seed}"]]})
    sampleInp  = sampleParams.join(models).filter{it[0].seed==1}.map{[it[0], it[1]+[inp:it[2]]]}
    traj       = sampler(sampleInp.map{[it[0],it[1]+[subDir:"${it[0].subDir}/trajs/iter${it[0].iter}"]]})
    qbcInps    = traj.flatMap{(1..it[0].seeds).collect(seed->[it[0]+[seed:seed], it[1]])}
        .join(models.filter{it[0].seed!=1}).map{[it[0], [ds:it[1], inp:it[2], subDir:"${it[0].subDir}/trajs/iter${it[0].iter}/qbclabel${it[0].seed}"]]}
    qbcLabels  = traj.mix(qbcLabel(qbcInps)).map{meta,inp -> [groupKey(meta+[seed:1],meta.seeds), inp]}.groupTuple()
    toLabel    = qbcFilter(qbcLabels.map{[it[0], [ds:it[1], subDir:"${it[0].subDir}/trajs/iter${it[0].iter}/qbcfilter"]]})
    labelInp2  = labelInp.join(toLabel).map{[it[0], it[1]+[ds:it[2], subDir:"${it[0].subDir}/trajs/iter${it[0].iter}/label"]]}
    labels     = labeller(labelInp2).map{[it[0], [ds:it[1]]]}
    augDs      = augFilter(labels.map{[it[0],it[1]+[subDir:"${it[0].subDir}/trajs/iter${it[0].iter}/filter"]]})
    restart    = resFilter(labels.map{[it[0],it[1]+[subDir:"${it[0].subDir}/trajs/iter${it[0].iter}/restart"]]})

    // prepare for the next iteration
    setNext(nextDs,     trainDs.join(augDs).flatMap{(1..it[0].seeds).collect{seed->[it[0]+[seed:seed], it[1]+[ds:([]<<it[1].ds<<it[2]).flatten()]]}})
    setNext(nextCkpt,   models.map{[it[0], [inp:it[1]]]})
    setNext(nextSteps,  trainSteps.map{[it[0], it[1]+[maxSteps: it[1].maxSteps+it[1].retrainSteps]]})
    setNext(nextSample, sampleParams.join(restart).map{[it[0],it[1]+[init:it[2]]]})
    setNext(nextLabel,  labelInp)

    emit:
    models.filter{it[0].seed==1}
        .map{[it[0].findAll{k,v->k!='iter'&k!='seed'}, it[1]]}
        .groupTuple().map{[it[0], it[1][-1]]}
}

workflow {
    explore(Channel.of([:]))
}
