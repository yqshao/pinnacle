# Extending TIPS

## Including a new model

Suppose you have a model that takes `.xyz` formated dataset and outputs the
weights you may define the training and evaluation process and include it in 
the workflow as shown below.

Here we assumed that your model runs in `bash` and outputs the parameters as
`weights.dat` and errors as `error.log`, and your model can be used as a
calculator in ASE. In other cases you might need to define also the process to
evaluate the model on a dataset, to sample a trajectory with the model, etc.

``` groovy
process myTrainer {
    input:
        ds dataset
    
    output:
        path "weights.dat"
        path "error.log"
        
    """
    mymodel train ${dataset.asXyz}
    """
}

model2ase = {"""
from mylib import myModel 
calc = myModel($it)
"""}

myModel = new Model(trianer: myTrainer, ase: model2ase)

workflow = activeLearn(myModel, ...)
```

