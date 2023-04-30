nextflow.enable.dsl=2

params.publish = "."

process mol2box {
  label 'molutils'
  publishDir "$params.publish/molutils/$name"

  input:
    tuple val(name), val(tags), val(box), val(seed)

  output:
    tuple val(name), path("${name}.xyz"), emit: geo
    tuple val(name), path("packmol.{in,log}"), emit: logs

  shell:
  '''
  rm -rf packmol.in

  cat << EOF >> packmol.in
  tolerance 2.0
  filetype xyz
  output !{name}.xyz
  seed !{seed}
  EOF

  mols=$(echo "!{tags}" | tr ";" "\\n")
  for mol in $mols
  do
    IFS=, read smiles number <<< "$mol"
    obabel -h --gen3d -oxyz -:"$smiles" -O $smiles.xyz
    cat << EOF >> packmol.in
  structure $smiles.xyz
    number $number
    inside cube 0. 0. 0. !{box-2.0}
  end structure
  EOF
  done

  packmol < packmol.in > packmol.log
  sed -i '2s/.*/Lattice="!{box} 0.0 0.0 0.0 !{box} 0.0 0.0 0.0 !{box}"/' !{name}.xyz
  '''
}
