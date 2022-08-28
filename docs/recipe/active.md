# Active workflows

The active learnign recipe provided by TIPS runs active learning loops. The
workflow is controlled by the several subworkflows, the trianing, sampling, and
labelling process. A typical workflow is shown below.

=== "Flowchart"

    ```mermaid
    flowchart TD
    %%{init:{'flowchart':{'nodeSpacing': 30, 'rankSpacing': 30, 'curve':'linear'}}}%%

    subgraph init [ ]
    initgeo
    initds
    initmodel
    end
    initgeo["gen${i}/geo"] ----> pinnMD([pinnMD])
    initmodel["gen${i}/model"] --> pinnTrain
    initmodel --> model
    initds["gen${i}/ds"] --> pinnTrain([pinnTrain])
    subgraph iteration [ ]
    pinnTrain --> model["model"]
    model --> pinnMD
    pinnMD --> traj[trajectory]
    traj --> tipsDiff([tipsDiff])
    traj --> cp2kLabel([cp2kLabel])
    cp2kLabel --> tipsDiff
    cp2kLabel --> tipsMix([tipsMix])
    end
    subgraph nx [ ]
    tipsDiff --> nxgeo["gen${i+1}/geo"]
    model ------> nxmodel["gen${i+1}/model"]
    initds --> tipsMix
    tipsMix --> nxds["gen${i+1}/ds"]
    end
    
    classDef hide fill:none, stroke:none;
    classDef data fill:#fff,stroke:#000,stroke-width:2px;
    classDef process stroke-width:2px;
    class init,nx hide;
    class pinnTrain,pinnMD,cp2kLabel,tipsMix,tipsDiff process;
    class initgeo,initds,initmodel,model,traj,nxmodel,nxgeo,nxds data;
    style iteration fill:none, stroke-dasharray: 6 2,stroke:#333
    linkStyle 1 stroke:#900,color:black;
    linkStyle 2 stroke:#090,color:black;
    ```

## Parameters

The workflow supports a number of parameters to adjust a runtime, most of those
parameters can be chose on your request with `tips wizard`.

### Initalization options

| Parameter      | Default    | Description                                           |
| -------------- | ---------- | ----------------------------------------------------- |
| `proj `        | user input | folder for storing results                            |
| `init_geo `    | user input | inital geometries for sampling                        |
| `init_model `  | user input | initial model or model parameters                     |
| `init_ds `     | user input | initial dataset                                       |
| `restart_from` | `false`    | restart from a given generation                       |
| `init_time `   | `1.0`      | sampling time scale in ps                             |
| `init_steps `  | `100000`   | training steps for initial model                      |
| `ens_size `    | `5`        | size of the model ensemble (for ensemble MD sampling) |
| `geo_size `    | `6`        | size of the starting geometries for sampling          |

### Iteration options

| Parameter      | Default            | Description                                         |
|----------------|--------------------|-----------------------------------------------------|
| `retrain_step` | `50000 `           | number of retrain steps per generation              |
| `label_flags ` | `'--nsample 50' `  | selection rule for the data to label                |
| `old_flag `    | `'--nsample 2700'` | selection rule for the old dataset                  |
| `new_flag `    | `'--nsample 300'`  | selection rule for the new dataset                  |
| `sp_points `   | `50`               | number of single points for each sampled trajectory |
| `acc_fac `     | `2.0 `             | factor to acceralate/slow down the sampling         |
| `min_time `    | `1.0 `             | minimal timescale for sampling                      |
| `emaxtol `     | `0.020 `           | toleranace for max error error                      |
| `ermsetol `    | `0.005 `           | toleranace for energy RMSE                          |
| `fmaxtol `     | `0.800 `           | toleranace for max force component error            |
| `frmsetol `    | `0.200 `           | toleranace for force RMSE                           |

See Also:

- biased dynamics with: `pinnMD`
- biased subsampling with: `subsample`
