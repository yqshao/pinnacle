#!/usr/bin/env nextflow

// The activated learning workflow  ======================================================
//
// The '--proj' parameter controls the output directory. See the parameters
// sections below for other parameters that can be tuned for the workflow.
//
//                                         written by  Yunqi Shao, first ver.: 2022.Aug.29
//                                                 adapted as PiNNAcLe recipe: 2023.Apr.24
//========================================================================================

nextflow.enable.dsl = 2
nextflow.preview.recursion = true

def logger (msg) {
  logfile = file("$params.publish/pinnacle.log")
  if (!logfile.getParent().exists()) {logfile.getParent().mkdirs()}
  logfile.append("$msg \n")
}

// entrypoint parameters ==================================================================
params.restart_from  = false
params.restart_conv  = false
params.init_geo      = 'input/geo/*.xyz'
params.init_model    = 'input/pinn/pinet-adam.yml'
params.init_ds       = 'input/dataset/init-ds.{yml,tfr}'
params.init_time     = 0.5
params.init_steps    = 200000
params.ens_size      = 1
params.geo_size      = 6
params.sp_points     = 10
//========================================================================================

// acle parameter ========================================================================
params.ref           = 'dftb' // reference (module name)
params.mpl           = 'pinn' // machine learning potential (module name)
params.frmsetol      = 0.150
params.ermsetol      = 0.005
params.fmaxtol       = 2.000
params.emaxtol       = 0.020
params.filters       = "" // ""--filter 'peratom(energy)<-125.1' --filter 'abs(force)<100.0'"
params.retrain_step  = 100000
params.collect_flags = '-f asetraj --subsample uniform --nsample 10 -of idx.xyz -o ds'
params.merge_flags   = '-f asetraj'
params.old_flag      = '--nsample 240' // 20% new data
params.new_flag      = '--psample 100' // 6*10 = 60 pts per iter
params.acc_fac       = 4.0
params.brake_fac     = 1.0 // do not slow down
params.min_time      = 0.5
params.max_time      = 1000.0 // 0.5 2  8 32 128 512 | stop at 1000 ps
params.max_gen       = 40     // 0   1  2  3   4   5
params.t_start       = 340
params.t_end         = 340
params.t_step        = 10
params.ref_inp       = null
params.train_flags   = '--log-every 10000 --ckpt-every 100000 --batch 1 --max-ckpts 1 --shuffle-buffer 3000'
params.train_init    = '--init'
params.md_flags      = '--ensemble nvt --dt 0.5 --log-every 100' // every 0.05 ps
params.exit_at_max_time = false
//========================================================================================

// Imports (publish directories are set here) ============================================
include { convert} from './module/tips.nf' addParams(publish: "$params.publish/collect")
include { dsmix } from './module/tips.nf' addParams(publish: "$params.publish/dsmix")
include { merge } from './module/tips.nf' addParams(publish: "$params.publish/merge")
include { check } from './module/tips.nf' addParams(publish: "$params.publish/check")
include { train } from "./module/${params.mpl}.nf" addParams(publish: "$params.publish/models")
include { md } from "./module/${params.mpl}.nf" addParams(publish: "$params.publish/md")
include { sp } from "./module/${params.ref}.nf" addParams(publish: "$params.publish/label")
//========================================================================================

// Main Iteration and Loops ==============================================================
workflow acle {
  take:
    ch_init

  main:
  loop.recurse(ch_init.first())
    .until{ it[0].toInteger()>params.max_gen || (it[5]>=params.max_time.toFloat() && params.exit_at_max_time) }
}

// Loop for each iteration =================================================================
workflow loop {
  take: ch_inp

  main:
  // retrain or keep the model ============================================================
  ch_inp \
    | branch {gen, geo, ds, models, step, time, converge -> \
              keep: converge
                return [gen, models]
              retrain: !converge
                return [gen, models, ds, (1..params.ens_size).toList(), step]} \
    | set {ch_model}


  ch_model.retrain.transpose(by:[1,3]) \
    | map {gen, model, ds, seed, steps -> \
           ["gen$gen/model$seed", ds, model,
            params.train_flags+
            " --seed $seed --train-steps $steps"+
            (gen.toInteger()==1?" $params.train_init":'')]}\
    | train

  train.out.model \
    | map {name, model -> (name=~/gen(\d+)\/model(\d+)/)[0][1,2]+[model]} \
    | map {gen, seed, model -> [gen, model]} \
    | mix (ch_model.keep.transpose()) \
    | groupTuple(size:params.ens_size) \
    | set {nx_models}
  //=======================================================================================

  // sampling with ensable NN =============================================================
  t_start = params.t_start.toFloat()
  t_end = params.t_end.toFloat()
  t_step = params.t_step.toFloat()
  ch_inp | map {[it[0], it[1], it[5]]} | transpose | set {ch_init_t} // init and time

  nx_models \
    | combine (ch_init_t, by:0)  \
    | map {gen, models, init, t -> \
           ["gen$gen/$init.baseName", models, init, params.md_flags+" --t $t"+" --T ${Math.max(t_end, t_start-gen.toInteger()*t_step)}"]} \
    | md
  md.out.traj.set {ch_trajs}
  //=======================================================================================

  // relabel with reference ===============================================================
  ref_inp = file(params.ref_inp)
  ch_trajs \
    | map {name, traj -> [name, traj, params.collect_flags]} \
    | convert \
    | flatMap {name, inps -> inps.collect {["$name/$it.baseName", it]}} \
    | map {name, geo -> [name, ref_inp, geo]} \
    | sp

  sp.out \
    | map {name, logs -> (name=~/(gen\d+\/.+)\/(\d+)/)[0][1,2]+[logs]} \
    | map {name, idx, logs -> [name, idx.toInteger(), logs]} \
    | groupTuple(size:params.sp_points) \
    | map {name, idx, logs -> [name, idx, logs, params.merge_flags]} \
    | merge \
    | set {ch_new_ds}
  //=======================================================================================

  // check convergence ====================================================================
  ch_new_ds \
    | join(ch_trajs) \
    | check \

  check.out \
    | map{name,geo,msg-> \
          [(name=~/gen(\d+)\/.+/)[0][1], geo, msg.contains('Converged')]} \
    | groupTuple(size:params.geo_size.toInteger()) \
    | map {gen, geo, conv -> [gen, geo, conv.every()]}
    | set {nx_geo_converge}

  //=======================================================================================

  // mix the new dataset ==================================================================
  ch_inp.map {[it[0], it[2]]}.set{ ch_old_ds }
  ch_new_ds \
    | map {name, idx, ds -> [(name=~/gen(\d+)\/.+/)[0][1], ds]} \
    | groupTuple(size:params.geo_size.toInteger()) \
    | join(ch_old_ds) \
    | map {it+[params.new_flag, params.old_flag]} \
    | dsmix \
    | set {nx_ds}
  //=======================================================================================

  // combine everything for new inputs ====================================================
  ch_inp.map{[it[0], it[4]]}.set {nx_step}
  ch_inp.map{[it[0], it[5]]}.set {nx_time}

  acc_fac = params.acc_fac.toFloat()
  brake_fac = params.brake_fac.toFloat()
  min_time = params.min_time.toFloat()
  max_time = params.max_time.toFloat()
  retrain_step = params.retrain_step.toInteger()

  nx_geo_converge | join(nx_models) | join(nx_ds) | join(nx_time) | join (nx_step) \
    | map {gen, geo, converge, models, ds, time, step -> \
           [(gen.toInteger()+1).toString(),
            geo, ds, models, \
            converge ? step : step+retrain_step, \
            converge ? Math.min(time*acc_fac, max_time) : Math.max(time*brake_fac, min_time), \
            converge]} \
    | set {nx_inp}
  //=======================================================================================

  ch_inp.subscribe {logger("[gen${it[0]}] ${it[-1]? 'not training': 'training'} the models.")}
  check.out.subscribe {name, geo, msg -> logger("[$name] ${msg.trim()}")}
  nx_inp.subscribe {logger('-'*80+'\n'+"[gen${it[0]}] next time scale ${it[5]} ps, ${it[6] ? 'no training planned' : 'next training step '+it[4]}.") } \

  emit:
  nx_inp
}
