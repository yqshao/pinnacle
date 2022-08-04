nextflow.enable.dsl=2

process pinnTrain {
    label 'pinn'
    publishDir "$publish", mode: 'link'

    input:
      val tag
      path dataset
      path input
      val flags
      val publish

    output:
      tuple val(tag), path('model/', type:'dir')

    script:
      """
      convert_flag="${(flags =~ /--seed[\s,\=]\d+/)[0]}"
      train_flags="${flags.replaceAll(/--seed[\s,\=]\d+/, '')}"
      pinn convert $dataset -o 'train:8,eval:2' \$convert_flag
      if [ ! -f $input/params.yml ];  then
          mkdir -p model; cp $input model/params.yml
      else
          cp -rL $input model
      fi
      pinn train model/params.yml --model-dir='model'\
          --train-ds='train.yml' --eval-ds='eval.yml'\
          \$train_flags
      """
}
