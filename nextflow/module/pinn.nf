nextflow.enable.dsl=2

params.publish = 'pinn'

process pinnTrain {
  tag "$name"
  label 'pinn'
  publishDir "$params.publish/$name"

  input:
    tuple val(name), path(dataset), path(input, stageAs:'input'), val(flags)

  output:
    tuple val(name), path('model', type:'dir'), emit: model
    tuple val(name), path('pinn.log'), emit: log

  script:
    convert_flag = "${(flags =~ /--seed[\s,\=]\d+/)[0]}"
    train_flags = "${flags.replaceAll(/--seed[\s,\=]\d+/, '')}"
    dataset = (dataset instanceof Path) ? dataset : dataset[0].baseName+'.yml'
    """
    #!/bin/bash

    pinn convert $dataset -o 'train:9,eval:1' $convert_flag

    if [ ! -f $input/params.yml ];  then
        mkdir -p model; cp $input model/params.yml
    else
        cp -rL $input model
    fi
    pinn train model/params.yml --model-dir='model'\
        --train-ds='train.yml' --eval-ds='eval.yml'\
        $train_flags
    pinn log model/eval > pinn.log
    """
}
