nextflow.enable.dsl=2

params.publish = 'analysis'

process rdf {
  tag "$name"
  label 'analysis'
  publishDir "$param.publish/$name"

  input:
    tuple val(name), path(traj), val(flags)

  output:
    tuple val(tag), path('rdf*.dat')

  script:
    """
    #!/usr/bin/env python3
    import re
    import itertools
    import numpy as np
    from ase.io import iread

    setup = {
      'tags': '',
      'rc': 5,
      'bins': 50,
      'start': 0
    }
    flags = {
      k: v for k,v in
        re.findall('--(.*?)[\\s,\\=]([^\\s]*)', "$flags")
    }
    setup.update(flags)

    def rdf_compute(iterator, tag):
        # this computes only O-O rdf here,
        assert (tag[0]=='O') and (tag[1]=='O')
        bins=np.linspace(0,float(setup['rc']), int(setup['bins']))
        r_mid = (bins[1:] + bins[:-1])/2
        bin_vol = (bins[1:]**3 - bins[:-1]**3)*4*np.pi/3
        hist = 0
        count = 0
        for data in iterator:
            count += 1
            self_pair = lambda x: np.array(list(itertools.combinations(x,2)))
            v = self_pair([idx for idx in range(len(data)) if data[idx].symbol == 'O'])
            diff = data.positions[v[:,0]] - data.positions[v[:,1]]
            diff = diff-np.rint(diff/data.cell[0][0])*data.cell[0][0]
            dist = np.linalg.norm(diff, axis=1)
            h, edges = np.histogram(dist, bins)
            hist += h/bin_vol/len(v)*np.prod(np.diag(data.cell))
        rdf = hist/count
        return r_mid, rdf

    tags = [tag.split('-') for tag in setup['tags'].split(',')]
    for tag in tags:
        traj = iread('$traj', index=f'{setup["start"]}:')
        try:
            r_mid, rdf = rdf_compute(traj, tag)
            np.savetxt(f"rdf_{''.join(tag)}.dat", np.stack([r_mid, rdf], axis=1))
        except:
            np.savetxt(f"rdf_{''.join(tag)}.dat", [np.nan])
    """
}


process mdlog {
  tag "$name"
  label 'analysis'
  publishDir "$param.publish/$name"

  input:
    tuple val(tag), path(traj), val(flags)

  output:
    tuple val(tag), path('*.dat')

  script:
    """
    #!/usr/bin/env python3
    import re
    import itertools
    import numpy as np
    from ase.io import iread

    setup = {
      'tags': 'density',
    }
    flags = {
      k: v for k,v in
        re.findall('--(.*?)[\\s,\\=]([^\\s]*)', "$flags")
    }
    setup.update(flags)

    # for now only density log supported, to be extended
    assert setup['tags']=='density'

    d_log = []
    traj = iread('$traj', index=':')
    for atoms in traj:
        density = atoms.get_masses().sum() / 6.022e23 / atoms.get_volume() * 1e30 /1e6
        d_log.append(density)
    np.savetxt('density.dat', d_log)
    """
}
