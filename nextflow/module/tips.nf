nextflow.enable.dsl=2

params.publish = "."

def space_sep(in) {(in instanceof Path) ?in :in.join(' ')}

process convert {
  label 'tips'
  publishDir "$params.publish/$name"

  input:
    tuple val(name), path(in, stageAs:'.in*/*'), val(flags)

  output:
    tuple val(name), path('*')

  script:
    """
    tips convert ${space_sep(in)} $flags
    """
}

process dsmix {
  label 'tips'
  publishDir "$params.publish/$name"
  input: tuple val(name), path(newDS, stageAs:'*.traj'), path(oldDS, stageAs:'old/*'), val(newFlag), val(oldFlag)
  output: tuple val(name), path('mix-ds.{tfr,yml}')

  script:
  """
  tips convert old/${oldDS[0].baseName}.yml -f pinn -o old-ds -of asetraj $oldFlag
  tips convert ${space_sep(newDS)} -f asetraj -o tmp.traj -of asetraj
  tips convert tmp.traj -f asetraj -o new-ds -of asetraj $newFlag
  tips convert new-ds.traj old-ds.traj -f asetraj -o mix-ds -of pinn --shuffle $params.filters
  rm {new-ds,old-ds,tmp}.*
  """
}

process merge {
  label 'tips'
  publishDir "$params.publish/$name"
  input: tuple val(name), val(idx), path(in, stageAs:'.in*/*'), val(flags)
  output: tuple val(name), path('merged.idx'), path('merged.traj')

  script:
  """
  printf "${idx.join('\\n')}" > merged.idx
  tips convert ${space_sep(in)} -o merged -of asetraj $flags
  """
}

process check {
  label 'tips'
  publishDir "$params.publish/$name"

  input:
  tuple val(name), path(idx), path(logs), path(traj)

  output:
  tuple val(name), path('*.xyz'), stdout

  script:
  fmaxtol = params.fmaxtol
  emaxtol = params.emaxtol
  frmsetol = params.frmsetol
  ermsetol = params.ermsetol
  sp_points = params.sp_points
  """
  #!/usr/bin/env python
  import numpy as np
  from ase import Atoms
  from ase.io import read, write
  from tips.io import load_ds
  from tips.io.filter import filters2fn

  filters = "$params.filters".replace("'", '').split(' ')[1::2]
  filter_fn = filters2fn(filters) # ^ a crude extractor

  idx = [int(i) for i in np.loadtxt("$idx")]
  logs = load_ds("$logs", fmt='asetraj')
  traj = load_ds("$traj", fmt='asetraj')

  idx, logs = tuple(zip(*(
      (i, datum) for (i, datum) in zip(idx, logs) if filter_fn(datum))))

  e_label = np.array([datum['energy']/len(datum['elem']) for datum in logs])
  f_label = np.array([datum['force'] for datum in logs])
  e_pred = np.array([traj[i]['energy']/len(traj[i]['elem']) for i in idx])
  f_pred = np.array([traj[i]['force'] for i in idx])

  ecnt = np.sum(np.abs(e_pred-e_label)>$emaxtol)
  fcnt = np.sum(np.any(np.abs(f_pred-f_label)>$fmaxtol,axis=(1,2)))
  emax = np.max(np.abs(e_pred-e_label))
  fmax = np.max(np.abs(f_pred-f_label))
  ermse = np.sqrt(np.mean((e_pred-e_label)**2))
  frmse = np.sqrt(np.mean((f_pred-f_label)**2))
  converged = (emax<$emaxtol) and (fmax<$fmaxtol) and (ermse<$ermsetol) and (frmse<$frmsetol) and (len(idx)==$sp_points)

  geoname = "$name".split('/')[1]
  if converged:
      msg = f'Converged; will restart from latest frame.'
      new_geo = logs[np.argmax(idx)]
  else:
      msg = f'energy: {ecnt}/{len(idx)} failed, max={emax:.2f} rmse={ermse:.2f}; '\
            f'force: {fcnt}/{len(idx)} failed, max={fmax:.2f} rmse={frmse:.2f}.'
      new_geo = logs[np.argmin(idx)]
  atoms = Atoms(new_geo['elem'], positions=new_geo['coord'], cell=new_geo['cell'],
                pbc=True)
  write(f'{geoname}.xyz', atoms)
  print(msg)
  """
}
