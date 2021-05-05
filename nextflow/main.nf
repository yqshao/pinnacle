#!/usr/bin/env nextflow
nextflow.enable.dsl=2
params.mode = 'explore'
include {trainer; labeller; filter; sampler} from "$moduleDir/adaptor"
include {explore} from "$moduleDir/explore"

workflow {
    input = Channel.of([null, [:]],)
    switch (params.mode){
        case 'explore':
            explore(input);
            break;
        case 'train':
            trainer(input);
            break;
        case 'label':
            labeller(input);
            break;
        case 'sample':
            sampler(input);
            break;
        case 'filter':
            filter(input);
            break;
        default:
            throw new Exception("Unkown mode $params.mode.");
    }
}
