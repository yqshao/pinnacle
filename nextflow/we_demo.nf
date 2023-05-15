nextflow.enable.dsl=2

// tuples of names and SMILES tags
inputs = [
  ['w32',       'O,32', '10.0'],
  ['w8e8', 'O,8;OCC,8', '10.0'],
  ['e10',     'OCC,10', '10.0'],
]

// --8<-- [start:params]
params.proj = 'we_demo'
params.init_flags = '--t 0.10 --log-every 1' // 200 steps each
params.init_seeds = 2
params.init_steps = 200000
params.init_model = 'input/pinn/pinet-hco-adam.yml'
params.init_time = 0.5
// --8<-- [end:params]

dftb_inp = file(params.dftb_calc)
geo_size = inputs.size * params.init_seeds.toInteger()

// --8<-- [start:include]
include { mol2box } from "./module/molutils.nf" addParams (publish: "$params.proj/init/pack")
include { md } from "./module/dftb.nf" addParams (publish:"$params.proj/init/md/")
include { convert as mkgeo } from "./module/tips.nf" addParams (publish:"$params.proj/init/ds")
include { convert as mkds } from "./module/tips.nf" addParams (publish:"$params.proj/init/geo")
include { acle } from "./acle.nf" addParams (
  publish: "$params.proj/acle",
  geo_size: geo_size)
// --8<-- [end:include]

// --8<-- [start:wf]
workflow entry {
  ch_inputs = Channel.fromList(inputs) \
    | combine(Channel.of(1..params.init_seeds)) \
    | map {name, tag, box, seed -> ["$name-$seed", tag, box.toFloat(), seed]} \
    | mol2box
// --8<-- [end:wf]

  mol2box.out.geo \
    | map {name, geo -> [name, dftb_inp, geo, params.init_flags]} \
    | md \

  md.out.traj \
    | map {name, traj -> traj} \
    | collect \
    | map {it -> ['.', it, '-f asetraj -o initds -of pinn']} \
    | mkds

  md.out.traj \
    | map {name, traj -> ['.', traj, "-f asetraj -o $name -of extxyz --subsample tail --nsample 1"] } \
    | mkgeo | map {it -> it[1]} \
    | collect | map {it -> [it]} \
    | combine(mkds.out) \
    | map {it -> ['1', it[0], it[2], [file(params.init_model)],
                  params.init_steps, params.init_time, false]} \
    | set { ch_init }

  acle (ch_init)
}

workflow {entry()}
