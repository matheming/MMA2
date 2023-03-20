This is a code review abput MMA, notice `rnorm` and `mvrnorm` is different, and Simple Averaging sometimes performs better than MMA.

In fact, from the soruce code of `mvrnorm` and file `sourcecode_mvrnorm`, if we set `mu=(0,0,0)` and `Sigma=diag()`, then the only difference between `rnorm` and `mvrnorm` is the column order.