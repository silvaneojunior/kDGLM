# Auxiliary function for forecasting

Auxiliary function for forecasting

## Usage

``` r
# S3 method for class 'fitted_dlm'
forecast(
  object,
  t = 1,
  plot = ifelse(requireNamespace("plotly", quietly = TRUE), "plotly",
    ifelse(requireNamespace("ggplot2", quietly = TRUE), "ggplot2", "base")),
  pred.cred = 0.95,
  safe.mode = TRUE,
  ...
)
```

## Arguments

- object:

  fitted_dlm object: The fitted model to be use for predictions.

- t:

  numeric: Time window for prediction.

- plot:

  boolean or character: A flag indicating if a plot should be produced.
  Should be one of FALSE, TRUE, 'base', 'ggplot2' or 'plotly'.

- pred.cred:

  numeric: The credibility level for the C.I.

- safe.mode:

  boolean: A flag indicating if consistency check should be performed at
  each time step. Recommended to be left on, but if you know what you
  are doing (i.e., you tested the model and it is safe) and need to fit
  it several times, you can disable the checks to save some time.

- ...:

  Extra variables necessary for prediction (covariates, etc.).

## Value

A list containing:

- data data.frame: A data frame contain the mean, variance and
  credibility intervals for the outcomes, including both the observed
  data and the predictions for future observations.

- forecast data.frame: Same as data, but restricted to predictions for
  future observations.

- outcomes list: A named list containing predictions for each outcome.
  Each element of this list is a list containing predictions (mean,
  variance and credibility intervals), the distribution of the linear
  predictor for the parameter of the observational model and the
  parameters of the predictive distribution (if available).

- theta.mean matrix: A matrix with the values for the latent states at
  each time. Dimensions are n x t, where n is the number of latent
  states

- theta.cov array: A 3D-array with the covariance of the latent states
  at each time. Dimensions are n x n x t, where n is the number of
  latent predictors.

- lambda.mean matrix: A matrix with the values for the linear predictors
  at each time. Dimensions are k x t, where k is the number of linear
  predictors

- lambda.cov array: A 3D-array with the covariance of the linear
  predictors at each time. Dimensions are k x k x t, where k is the
  number of linear predictors.

- plot (if so chosen): A plotly or ggplot object.

A list containing:

- data data.frame: A table with the model evaluated at each observed
  time, plus the forecasted period.

- forecast data.frame: A table with the model evaluated at the
  forecasted period.

- outcomes list: A list containing the parameters of the predictive
  distribution for each outcome at the forecasted period.

- theta.mean matrix: The mean of the latent states at each forecasted
  time. Dimensions are n x t.forecast, where t.forecast is the size of
  the forecast windows and n is the number of latent states.

- theta.cov array: A 3D-array containing the covariance matrix of the
  latent states at each forecasted time. Dimensions are n x n x
  t.forecast, where t.forecast is the size of the forecast windows and n
  is the number of latent states.

- lambda.mean matrix: The mean of the linear predictor at each
  forecasted time. Dimensions are k x t.forecast, where t.forecast is
  the size of the forecast windows and k is the number of linear
  predictors.

- lambda.cov array: A 3D-array containing the covariance matrix for the
  linear predictor at each forecasted time. Dimensions are k x k x
  t.forecast, where t.forecast is the size of the forecast windows and k
  is the number of linear predictors.

## Details

If an a covariate is necessary for forecasting, it should be passed as a
named argument. Its name must follow this structure: \<block
name\>.Covariate\<.index\>. If there is only one covariate in the
associated block the index is omitted. If an a pulse is necessary for
forecasting, it should be passed as a named argument. Its name must
follow this structure: \<block name\>.Pulse\<.index\>. If there is only
one pulse in the associated block the index is omitted. The user may
pass the observed values at the prediction windows (optional). See
example. As an special case, if the model has an Multinomial outcome,
the user may pass the N parameter instead of the observations. If an
offset is necessary for forecasting, it should be passed with the same
syntax as the observed data. See example.

## See also

Other auxiliary functions for fitted_dlm objects:
[`coef.fitted_dlm()`](coef.fitted_dlm.md),
[`eval_dlm_norm_const()`](eval_dlm_norm_const.md),
[`fit_model()`](fit_model.md), [`kdglm()`](kdglm.md),
[`simulate.fitted_dlm()`](simulate.fitted_dlm.md),
[`smoothing()`](smoothing.md),
[`update.fitted_dlm()`](update.fitted_dlm.md)

## Examples

``` r
structure <-
  polynomial_block(p = 1, order = 2, D = 0.95) +
  harmonic_block(p = 1, period = 12, D = 0.975) +
  noise_block(p = 1, R1 = 0.1) +
  regression_block(
    p = chickenPox$date >= as.Date("2013-09-1"),
    # Vaccine was introduced in September of 2013
    name = "Vaccine"
  )

outcome <- Multinom(p = c("p.1", "p.2"), data = chickenPox[, c(2, 3, 5)])
fitted.data <- fit_model(structure * 2,
  chickenPox = outcome
)

forecast(fitted.data, 24,
  chickenPox = list(Total = rep(175, 24)), # Optional
  Vaccine.1.Covariate = rep(TRUE, 24),
  Vaccine.2.Covariate = rep(TRUE, 24)
)
#> Scale for y is already present.
#> Adding another scale for y, which will replace the existing scale.
#> $data
#>     Time                   Serie Observation Prediction   Variance C.I.lower
#> 1      1     chickenPox.< 5 year         218  213.58336  1970.0224       127
#> 2      2     chickenPox.< 5 year         215  213.08166  1826.3294       129
#> 3      3     chickenPox.< 5 year         256  248.92116  2230.5153       156
#> 4      4     chickenPox.< 5 year         271  262.87385  2058.2753       172
#> 5      5     chickenPox.< 5 year         348  333.57497  2708.1150       228
#> 6      6     chickenPox.< 5 year         333  320.15775  2177.4260       225
#> 7      7     chickenPox.< 5 year         397  385.62449  2918.9654       275
#> 8      8     chickenPox.< 5 year         459  478.25619  4110.7692       346
#> 9      9     chickenPox.< 5 year         811  769.97956  8030.5584       583
#> 10    10     chickenPox.< 5 year         932  889.36173 10562.4534       675
#> 11    11     chickenPox.< 5 year         730  694.17741  7446.8652       516
#> 12    12     chickenPox.< 5 year         460  453.79876  4141.3474       324
#> 13    13     chickenPox.< 5 year         214  227.84081  1349.9293       155
#> 14    14     chickenPox.< 5 year         133  148.90345   670.2821        98
#> 15    15     chickenPox.< 5 year         162  155.79359   724.2655       103
#> 16    16     chickenPox.< 5 year         227  196.43298  1054.3681       133
#> 17    17     chickenPox.< 5 year         231  212.43181  1125.4755       146
#> 18    18     chickenPox.< 5 year         179  196.43252   863.3400       138
#> 19    19     chickenPox.< 5 year         256  262.10029  1264.7064       191
#> 20    20     chickenPox.< 5 year         343  332.63869  1693.4840       249
#> 21    21     chickenPox.< 5 year         403  427.57466  2512.0559       325
#> 22    22     chickenPox.< 5 year         544  514.78942  3454.2190       394
#> 23    23     chickenPox.< 5 year         518  482.05868  3224.7552       366
#> 24    24     chickenPox.< 5 year         446  422.45053  2826.6894       315
#> 25    25     chickenPox.< 5 year         248  253.12430  1234.4500       183
#> 26    26     chickenPox.< 5 year         126  134.75607   434.3968        94
#> 27    27     chickenPox.< 5 year         113  143.53711   524.4366        99
#> 28    28     chickenPox.< 5 year         152  140.86369   490.6935        98
#> 29    29     chickenPox.< 5 year         157  154.50156   525.3190       109
#> 30    30     chickenPox.< 5 year         181  177.78759   571.2457       130
#> 31    31     chickenPox.< 5 year         222  226.91655   735.3450       172
#> 32    32     chickenPox.< 5 year         292  284.88427   946.1142       222
#> 33    33     chickenPox.< 5 year         521  489.48026  2425.7620       389
#> 34    34     chickenPox.< 5 year         587  574.50167  3318.6183       457
#> 35    35     chickenPox.< 5 year         495  456.24162  2311.7290       358
#> 36    36     chickenPox.< 5 year         363  323.04310  1382.0081       248
#> 37    37     chickenPox.< 5 year         157  185.44184   580.5863       137
#> 38    38     chickenPox.< 5 year         102  118.11201   294.5789        84
#> 39    39     chickenPox.< 5 year         132  136.94759   412.6463        97
#> 40    40     chickenPox.< 5 year         155  144.29627   449.6885       103
#> 41    41     chickenPox.< 5 year         217  211.47770   852.7351       154
#> 42    42     chickenPox.< 5 year         364  326.65769  1709.7904       245
#> 43    43     chickenPox.< 5 year         333  321.96834  1459.9463       246
#> 44    44     chickenPox.< 5 year         308  318.37241  1329.8531       245
#> 45    45     chickenPox.< 5 year         535  534.98008  2980.6070       423
#> 46    46     chickenPox.< 5 year         705  683.24551  4778.4721       542
#> 47    47     chickenPox.< 5 year         623  570.26250  3575.8874       449
#> 48    48     chickenPox.< 5 year         542  465.57851  2736.3861       360
#> 49    49     chickenPox.< 5 year         295  277.52804  1216.1768       208
#> 50    50     chickenPox.< 5 year         108  124.96477   321.2862        90
#> 51    51     chickenPox.< 5 year          85  108.54534   278.8636        76
#> 52    52     chickenPox.< 5 year         122  125.04087   368.6349        88
#> 53    53     chickenPox.< 5 year         145  159.28629   545.7949       114
#> 54    54     chickenPox.< 5 year         208  187.57580   657.4364       137
#> 55    55     chickenPox.< 5 year         191  202.02059   680.1540       150
#> 56    56     chickenPox.< 5 year         191  195.93931   597.7855       147
#> 57    57     chickenPox.< 5 year         289  291.02267  1213.6264       221
#> 58    58     chickenPox.< 5 year         374  358.47281  1846.4432       273
#> 59    59     chickenPox.< 5 year         305  285.28703  1277.6387       214
#> 60    60     chickenPox.< 5 year         297  262.95845  1206.1113       194
#> 61    61     chickenPox.< 5 year         156  156.87922   517.1706       112
#> 62    62     chickenPox.< 5 year          64   77.43775   159.3300        53
#> 63    63     chickenPox.< 5 year          97   98.33762   257.1332        67
#> 64    64     chickenPox.< 5 year         119   98.02634   255.3345        67
#> 65    65     chickenPox.< 5 year          96  101.10096   268.8404        70
#> 66    66     chickenPox.< 5 year         108  113.88600   318.7548        79
#> 67    67     chickenPox.< 5 year         127  141.19799   446.6542       100
#> 68    68     chickenPox.< 5 year         143  146.91046   459.8712       105
#> 69    69     chickenPox.< 5 year         191  184.75125   698.0615       133
#> 70    70     chickenPox.< 5 year         195  191.66894   770.6620       137
#> 71    71     chickenPox.< 5 year         153  146.94228   494.2250       104
#> 72    72     chickenPox.< 5 year         134  128.31627   409.6046        89
#> 73    73     chickenPox.< 5 year          79   92.57443   242.7741        63
#> 74    74     chickenPox.< 5 year          76   72.26792   162.2861        48
#> 75    75     chickenPox.< 5 year          88   81.43944   203.9909        54
#> 76    76     chickenPox.< 5 year          76   79.88421   195.5635        53
#> 77    77     chickenPox.< 5 year          84   91.69628   239.5374        62
#> 78    78     chickenPox.< 5 year          94   89.70525   217.3170        61
#> 79    79     chickenPox.< 5 year         126  124.00926   375.4526        87
#> 80    80     chickenPox.< 5 year         135  134.15747   432.9562        94
#> 81    81     chickenPox.< 5 year         152  148.87920   536.6530       104
#> 82    82     chickenPox.< 5 year         193  192.49156   903.0650       135
#> 83    83     chickenPox.< 5 year         151  146.47017   562.1278       101
#> 84    84     chickenPox.< 5 year         146  133.51518   487.5810        91
#> 85    85     chickenPox.< 5 year         105  108.33834   341.0566        73
#> 86    86     chickenPox.< 5 year          72   77.06402   187.9047        51
#> 87    87     chickenPox.< 5 year          90   94.09285   265.7415        63
#> 88    88     chickenPox.< 5 year          77   73.92722   169.1858        49
#> 89    89     chickenPox.< 5 year          90   86.80048   218.2579        59
#> 90    90     chickenPox.< 5 year          87   94.51491   246.3712        65
#> 91    91     chickenPox.< 5 year         110  106.02820   297.5535        73
#> 92    92     chickenPox.< 5 year         134  129.65846   436.4795        90
#> 93    93     chickenPox.< 5 year         132  133.79274   481.6768        92
#> 94    94     chickenPox.< 5 year         151  149.27467   616.5594       102
#> 95    95     chickenPox.< 5 year         138  127.81807   477.3924        86
#> 96    96     chickenPox.< 5 year         104  103.70684   331.9335        69
#> 97    97     chickenPox.< 5 year          67   91.03763   263.6324        60
#> 98    98     chickenPox.< 5 year          58   64.29270   140.5495        42
#> 99    99     chickenPox.< 5 year          94   85.36832   219.0392        57
#> 100  100     chickenPox.< 5 year          79   72.95833   160.3592        49
#> 101  101     chickenPox.< 5 year          93   93.88886   239.7493        64
#> 102  102     chickenPox.< 5 year          84   76.40348   162.0571        52
#> 103  103     chickenPox.< 5 year          96  100.03248   260.1230        69
#> 104  104     chickenPox.< 5 year          82   83.73816   194.3599        57
#> 105  105     chickenPox.< 5 year          80   88.56045   222.3970        60
#> 106  106     chickenPox.< 5 year          92   90.81995   242.2150        61
#> 107  107     chickenPox.< 5 year          77   80.35493   199.4477        54
#> 108  108     chickenPox.< 5 year          79   81.13588   202.1706        54
#> 109  109     chickenPox.< 5 year          73   72.87132   163.6432        49
#> 110  110     chickenPox.< 5 year          54   64.19351   127.1832        43
#> 111  111     chickenPox.< 5 year         133  101.96153   266.8406        71
#> 112  112     chickenPox.< 5 year         162  104.46116   276.8494        73
#> 113  113     chickenPox.< 5 year         102   86.42489   207.4812        59
#> 114  114     chickenPox.< 5 year          70   75.20901   169.4770        50
#> 115  115     chickenPox.< 5 year          79   80.53083   202.3640        54
#> 116  116     chickenPox.< 5 year          59   74.58019   193.4642        48
#> 117  117     chickenPox.< 5 year          49   67.78072   178.1302        43
#> 118  118     chickenPox.< 5 year          63   63.84081   171.8494        39
#> 119  119     chickenPox.< 5 year          63   66.32271   195.3912        40
#> 120  120     chickenPox.< 5 year          48   51.15106   132.2350        30
#> 121    1 chickenPox.5 to 9 years          64   65.59179  1970.0224        16
#> 122    2 chickenPox.5 to 9 years          43   53.88936  1826.3294        11
#> 123    3 chickenPox.5 to 9 years          60   61.31019  2230.5153        13
#> 124    4 chickenPox.5 to 9 years          47   55.70363  2058.2753        11
#> 125    5 chickenPox.5 to 9 years          66   69.43924  2708.1150        16
#>     C.I.upper type
#> 1         299  Fit
#> 2         295  Fit
#> 3         339  Fit
#> 4         348  Fit
#> 5         430  Fit
#> 6         406  Fit
#> 7         484  Fit
#> 8         595  Fit
#> 9         931  Fit
#> 10       1075  Fit
#> 11        852  Fit
#> 12        574  Fit
#> 13        298  Fit
#> 14        199  Fit
#> 15        208  Fit
#> 16        259  Fit
#> 17        277  Fit
#> 18        252  Fit
#> 19        329  Fit
#> 20        409  Fit
#> 21        520  Fit
#> 22        623  Fit
#> 23        588  Fit
#> 24        522  Fit
#> 25        320  Fit
#> 26        175  Fit
#> 27        188  Fit
#> 28        184  Fit
#> 29        199  Fit
#> 30        223  Fit
#> 31        278  Fit
#> 32        342  Fit
#> 33        581  Fit
#> 34        681  Fit
#> 35        546  Fit
#> 36        393  Fit
#> 37        231  Fit
#> 38        151  Fit
#> 39        177  Fit
#> 40        186  Fit
#> 41        268  Fit
#> 42        406  Fit
#> 43        395  Fit
#> 44        387  Fit
#> 45        636  Fit
#> 46        812  Fit
#> 47        682  Fit
#> 48        564  Fit
#> 49        344  Fit
#> 50        160  Fit
#> 51        141  Fit
#> 52        163  Fit
#> 53        205  Fit
#> 54        237  Fit
#> 55        252  Fit
#> 56        243  Fit
#> 57        357  Fit
#> 58        440  Fit
#> 59        354  Fit
#> 60        330  Fit
#> 61        201  Fit
#> 62        102  Fit
#> 63        130  Fit
#> 64        130  Fit
#> 65        134  Fit
#> 66        149  Fit
#> 67        183  Fit
#> 68        189  Fit
#> 69        236  Fit
#> 70        246  Fit
#> 71        190  Fit
#> 72        168  Fit
#> 73        123  Fit
#> 74         98  Fit
#> 75        110  Fit
#> 76        108  Fit
#> 77        123  Fit
#> 78        119  Fit
#> 79        162  Fit
#> 80        175  Fit
#> 81        195  Fit
#> 82        252  Fit
#> 83        193  Fit
#> 84        177  Fit
#> 85        145  Fit
#> 86        105  Fit
#> 87        127  Fit
#> 88        100  Fit
#> 89        116  Fit
#> 90        126  Fit
#> 91        140  Fit
#> 92        171  Fit
#> 93        178  Fit
#> 94        199  Fit
#> 95        172  Fit
#> 96        140  Fit
#> 97        124  Fit
#> 98         88  Fit
#> 99        115  Fit
#> 100        98  Fit
#> 101       125  Fit
#> 102       102  Fit
#> 103       132  Fit
#> 104       112  Fit
#> 105       118  Fit
#> 106       122  Fit
#> 107       109  Fit
#> 108       110  Fit
#> 109        99  Fit
#> 110        87  Fit
#> 111       135  Fit
#> 112       138  Fit
#> 113       115  Fit
#> 114       101  Fit
#> 115       109  Fit
#> 116       103  Fit
#> 117        95  Fit
#> 118        91  Fit
#> 119        95  Fit
#> 120        75  Fit
#> 121       141  Fit
#> 122       122  Fit
#> 123       136  Fit
#> 124       126  Fit
#> 125       152  Fit
#>  [ reached 'max' / getOption("max.print") -- omitted 307 rows ]
#> 
#> $forecast
#>    Time                     Serie Observation Variance Prediction C.I.lower
#> 1   121       chickenPox.< 5 year          NA 135.5736   49.71787        28
#> 2   122       chickenPox.< 5 year          NA 135.5736   48.69890        27
#> 3   123       chickenPox.< 5 year          NA 135.5736   48.20482        26
#> 4   124       chickenPox.< 5 year          NA 135.5736   48.26928        26
#> 5   125       chickenPox.< 5 year          NA 135.5736   48.74705        26
#> 6   126       chickenPox.< 5 year          NA 135.5736   49.30286        25
#> 7   127       chickenPox.< 5 year          NA 135.5736   49.51656        25
#> 8   128       chickenPox.< 5 year          NA 135.5736   49.07372        23
#> 9   129       chickenPox.< 5 year          NA 135.5736   47.93041        22
#> 10  130       chickenPox.< 5 year          NA 135.5736   46.30816        20
#> 11  131       chickenPox.< 5 year          NA 135.5736   44.51784        18
#> 12  132       chickenPox.< 5 year          NA 135.5736   42.80074        17
#> 13  133       chickenPox.< 5 year          NA 135.5736   41.32312        15
#> 14  134       chickenPox.< 5 year          NA 135.5736   40.24110        14
#> 15  135       chickenPox.< 5 year          NA 135.5736   39.69342        14
#> 16  136       chickenPox.< 5 year          NA 135.5736   39.72013        13
#> 17  137       chickenPox.< 5 year          NA 135.5736   40.19498        13
#> 18  138       chickenPox.< 5 year          NA 135.5736   40.82722        12
#> 19  139       chickenPox.< 5 year          NA 135.5736   41.24545        12
#> 20  140       chickenPox.< 5 year          NA 135.5736   41.14590        11
#> 21  141       chickenPox.< 5 year          NA 135.5736   40.42889        10
#> 22  142       chickenPox.< 5 year          NA 135.5736   39.21549         9
#> 23  143       chickenPox.< 5 year          NA 135.5736   37.73478         8
#> 24  144       chickenPox.< 5 year          NA 135.5736   36.21121         7
#> 25  121   chickenPox.5 to 9 years          NA 135.5736   21.65348         8
#> 26  122   chickenPox.5 to 9 years          NA 135.5736   18.60570         5
#> 27  123   chickenPox.5 to 9 years          NA 135.5736   16.39643         4
#> 28  124   chickenPox.5 to 9 years          NA 135.5736   15.30112         3
#> 29  125   chickenPox.5 to 9 years          NA 135.5736   15.34896         3
#> 30  126   chickenPox.5 to 9 years          NA 135.5736   16.46339         3
#> 31  127   chickenPox.5 to 9 years          NA 135.5736   18.45239         4
#> 32  128   chickenPox.5 to 9 years          NA 135.5736   20.88650         5
#> 33  129   chickenPox.5 to 9 years          NA 135.5736   23.03546         6
#> 34  130   chickenPox.5 to 9 years          NA 135.5736   24.07000         6
#> 35  131   chickenPox.5 to 9 years          NA 135.5736   23.50624         5
#> 36  132   chickenPox.5 to 9 years          NA 135.5736   21.53524         4
#> 37  133   chickenPox.5 to 9 years          NA 135.5736   18.90456         3
#> 38  134   chickenPox.5 to 9 years          NA 135.5736   16.45364         2
#> 39  135   chickenPox.5 to 9 years          NA 135.5736   14.73308         1
#> 40  136   chickenPox.5 to 9 years          NA 135.5736   13.95315         1
#> 41  137   chickenPox.5 to 9 years          NA 135.5736   14.12701         1
#> 42  138   chickenPox.5 to 9 years          NA 135.5736   15.18255         1
#> 43  139   chickenPox.5 to 9 years          NA 135.5736   16.94978         1
#> 44  140   chickenPox.5 to 9 years          NA 135.5736   19.06053         1
#> 45  141   chickenPox.5 to 9 years          NA 135.5736   20.90129         2
#> 46  142   chickenPox.5 to 9 years          NA 135.5736   21.78778         2
#> 47  143   chickenPox.5 to 9 years          NA 135.5736   21.33572         1
#> 48  144   chickenPox.5 to 9 years          NA 135.5736   19.72452         1
#> 49  121 chickenPox.15 to 49 years          NA 135.5736   96.62864        72
#> 50  122 chickenPox.15 to 49 years          NA 135.5736  100.69540        75
#> 51  123 chickenPox.15 to 49 years          NA 135.5736  103.39875        77
#> 52  124 chickenPox.15 to 49 years          NA 135.5736  104.42960        78
#> 53  125 chickenPox.15 to 49 years          NA 135.5736  103.90399        76
#> 54  126 chickenPox.15 to 49 years          NA 135.5736  102.23375        73
#> 55  127 chickenPox.15 to 49 years          NA 135.5736  100.03105        70
#> 56  128 chickenPox.15 to 49 years          NA 135.5736   98.03978        67
#> 57  129 chickenPox.15 to 49 years          NA 135.5736   97.03413        65
#> 58  130 chickenPox.15 to 49 years          NA 135.5736   97.62184        64
#> 59  131 chickenPox.15 to 49 years          NA 135.5736   99.97592        66
#> 60  132 chickenPox.15 to 49 years          NA 135.5736  103.66402        69
#> 61  133 chickenPox.15 to 49 years          NA 135.5736  107.77232        73
#> 62  134 chickenPox.15 to 49 years          NA 135.5736  111.30526        76
#> 63  135 chickenPox.15 to 49 years          NA 135.5736  113.57350        78
#> 64  136 chickenPox.15 to 49 years          NA 135.5736  114.32672        78
#> 65  137 chickenPox.15 to 49 years          NA 135.5736  113.67801        76
#> 66  138 chickenPox.15 to 49 years          NA 135.5736  111.99023        72
#> 67  139 chickenPox.15 to 49 years          NA 135.5736  109.80477        68
#> 68  140 chickenPox.15 to 49 years          NA 135.5736  107.79357        65
#> 69  141 chickenPox.15 to 49 years          NA 135.5736  106.66982        62
#> 70  142 chickenPox.15 to 49 years          NA 135.5736  106.99673        61
#> 71  143 chickenPox.15 to 49 years          NA 135.5736  108.92951        62
#> 72  144 chickenPox.15 to 49 years          NA 135.5736  112.06428        65
#>    C.I.upper
#> 1         74
#> 2         73
#> 3         73
#> 4         74
#> 5         75
#> 6         77
#> 7         78
#> 8         79
#> 9         79
#> 10        78
#> 11        77
#> 12        76
#> 13        74
#> 14        74
#> 15        74
#> 16        75
#> 17        76
#> 18        79
#> 19        81
#> 20        82
#> 21        83
#> 22        83
#> 23        82
#> 24        80
#> 25        41
#> 26        37
#> 27        35
#> 28        34
#> 29        34
#> 30        37
#> 31        40
#> 32        45
#> 33        49
#> 34        51
#> 35        51
#> 36        49
#> 37        46
#> 38        42
#> 39        40
#> 40        40
#> 41        41
#> 42        44
#> 43        48
#> 44        53
#> 45        57
#> 46        60
#> 47        60
#> 48        58
#> 49       121
#> 50       125
#> 51       128
#> 52       129
#> 53       130
#> 54       129
#> 55       128
#> 56       127
#> 57       127
#> 58       129
#> 59       132
#> 60       135
#> 61       139
#> 62       142
#> 63       144
#> 64       145
#> 65       146
#> 66       146
#> 67       145
#> 68       145
#> 69       145
#> 70       146
#> 71       148
#> 72       150
#> 
#> $outcomes
#> $outcomes$chickenPox
#> $outcomes$chickenPox$conj.param
#>       alpha_1  alpha_2   alpha_3
#> 121 16.905993 7.363019 32.857462
#> 122 15.633856 5.973007 32.326348
#> 123 14.538583 4.945167 31.185083
#> 124 13.558861 4.298091 29.334317
#> 125 12.629936 3.976782 26.920622
#> 126 11.687022 3.902573 24.234050
#> 127 10.690872 3.983962 21.597200
#> 128  9.659035 4.111028 19.296883
#> 129  8.662543 4.163236 17.537140
#> 130  7.781518 4.044669 16.404155
#> 131  7.061989 3.728861 15.859457
#> 132  6.503896 3.272443 15.752533
#> 133  6.072196 2.777916 15.836524
#> 134  5.717802 2.337875 15.815208
#> 135  5.399034 2.003969 15.448081
#> 136  5.090677 1.788287 14.652534
#> 137  4.777907 1.679252 13.512707
#> 138  4.449060 1.654486 12.203897
#> 139  4.099407 1.684648 10.913553
#> 140  3.740901 1.732945  9.800373
#> 141  3.402080 1.758838  8.976237
#> 142  3.113179 1.729655  8.494092
#> 143  2.890148 1.634126  8.343031
#> 144  2.729502 1.486780  8.447100
#> 
#> $outcomes$chickenPox$ft
#>           [,1]       [,2]       [,3]       [,4]       [,5]       [,6]
#> [1,] -0.679083 -0.7432195 -0.7818002 -0.7919054 -0.7782452 -0.7518979
#> [2,] -1.549856 -1.7591244 -1.9299130 -2.0242632 -2.0246982 -1.9389059
#>            [,7]       [,8]       [,9]      [,10]      [,11]      [,12]
#> [1,] -0.7273409 -0.7185724 -0.7353596 -0.7806223 -0.8496504 -0.9313656
#> [2,] -1.7976785 -1.6466622 -1.5341261 -1.4980284 -1.5558459 -1.6998908
#>          [,13]     [,14]     [,15]     [,16]     [,17]     [,18]     [,19]
#> [1,] -1.011290 -1.075427 -1.114008 -1.124113 -1.110453 -1.084105 -1.059548
#> [2,] -1.899371 -2.108639 -2.279428 -2.373778 -2.374213 -2.288421 -2.147194
#>          [,20]     [,21]     [,22]     [,23]     [,24]
#> [1,] -1.050780 -1.067567 -1.112830 -1.181858 -1.263573
#> [2,] -1.996177 -1.883641 -1.847543 -1.905361 -2.049406
#> 
#> $outcomes$chickenPox$Qt
#> , , 1
#> 
#>             [,1]        [,2]
#> [1,] 0.111031031 0.002952275
#> [2,] 0.002952275 0.114880251
#> 
#> , , 2
#> 
#>             [,1]        [,2]
#> [1,] 0.123672883 0.003275323
#> [2,] 0.003275323 0.130104017
#> 
#> , , 3
#> 
#>             [,1]        [,2]
#> [1,] 0.136830380 0.003541592
#> [2,] 0.003541592 0.145910042
#> 
#> , , 4
#> 
#>             [,1]        [,2]
#> [1,] 0.149593768 0.003738337
#> [2,] 0.003738337 0.161078684
#> 
#> , , 5
#> 
#>             [,1]        [,2]
#> [1,] 0.161467867 0.003876132
#> [2,] 0.003876132 0.174921870
#> 
#> , , 6
#> 
#>             [,1]        [,2]
#> [1,] 0.172493799 0.003984471
#> [2,] 0.003984471 0.187492253
#> 
#> , , 7
#> 
#>             [,1]        [,2]
#> [1,] 0.183192092 0.004103871
#> [2,] 0.004103871 0.199515810
#> 
#> , , 8
#> 
#>             [,1]        [,2]
#> [1,] 0.194434525 0.004277022
#> [2,] 0.004277022 0.212174787
#> 
#> , , 9
#> 
#>             [,1]        [,2]
#> [1,] 0.207308809 0.004539212
#> [2,] 0.004539212 0.226871186
#> 
#> , , 10
#> 
#>             [,1]        [,2]
#> [1,] 0.222942264 0.004908079
#> [2,] 0.004908079 0.244987991
#> 
#> , , 11
#> 
#>             [,1]        [,2]
#> [1,] 0.242227680 0.005375258
#> [2,] 0.005375258 0.267588593
#> 
#> , , 12
#> 
#>             [,1]        [,2]
#> [1,] 0.265491167 0.005904707
#> [2,] 0.005904707 0.295057987
#> 
#> , , 13
#> 
#>             [,1]        [,2]
#> [1,] 0.292266282 0.006441238
#> [2,] 0.006441238 0.326834106
#> 
#> , , 14
#> 
#>             [,1]        [,2]
#> [1,] 0.321351876 0.006927658
#> [2,] 0.006927658 0.361440228
#> 
#> , , 15
#> 
#>             [,1]        [,2]
#> [1,] 0.351188944 0.007323712
#> [2,] 0.007323712 0.396907427
#> 
#> , , 16
#> 
#>             [,1]        [,2]
#> [1,] 0.380397001 0.007618777
#> [2,] 0.007618777 0.431442027
#> 
#> , , 17
#> 
#>             [,1]        [,2]
#> [1,] 0.408221257 0.007833981
#> [2,] 0.007833981 0.464034915
#> 
#> , , 18
#> 
#>             [,1]        [,2]
#> [1,] 0.434723899 0.008014972
#> [2,] 0.008014972 0.494756726
#> 
#> , , 19
#> 
#>             [,1]        [,2]
#> [1,] 0.460721555 0.008219696
#> [2,] 0.008219696 0.524685622
#> 
#> , , 20
#> 
#>             [,1]        [,2]
#> [1,] 0.487577799 0.008504881
#> [2,] 0.008504881 0.555595865
#> 
#> , , 21
#> 
#>            [,1]       [,2]
#> [1,] 0.51693605 0.00891269
#> [2,] 0.00891269 0.58956268
#> 
#> , , 22
#> 
#>            [,1]       [,2]
#> [1,] 0.55039436 0.00945864
#> [2,] 0.00945864 0.62854309
#> 
#> , , 23
#> 
#>            [,1]       [,2]
#> [1,] 0.58910514 0.01012381
#> [2,] 0.01012381 0.67392152
#> 
#> , , 24
#> 
#>            [,1]       [,2]
#> [1,] 0.63337341 0.01085602
#> [2,] 0.01085602 0.72606499
#> 
#> 
#> $outcomes$chickenPox$data
#>       [,1] [,2] [,3]
#>  [1,]   48   26   94
#>  [2,]   48   26   94
#>  [3,]   48   26   94
#>  [4,]   48   26   94
#>  [5,]   48   26   94
#>  [6,]   48   26   94
#>  [7,]   48   26   94
#>  [8,]   48   26   94
#>  [9,]   48   26   94
#> [10,]   48   26   94
#> [11,]   48   26   94
#> [12,]   48   26   94
#> [13,]   48   26   94
#> [14,]   48   26   94
#> [15,]   48   26   94
#> [16,]   48   26   94
#> [17,]   48   26   94
#> [18,]   48   26   94
#> [19,]   48   26   94
#> [20,]   48   26   94
#> [21,]   48   26   94
#> [22,]   48   26   94
#> [23,]   48   26   94
#> [24,]   48   26   94
#> 
#> $outcomes$chickenPox$offset
#>       [,1] [,2] [,3]
#>  [1,]    1    1    1
#>  [2,]    1    1    1
#>  [3,]    1    1    1
#>  [4,]    1    1    1
#>  [5,]    1    1    1
#>  [6,]    1    1    1
#>  [7,]    1    1    1
#>  [8,]    1    1    1
#>  [9,]    1    1    1
#> [10,]    1    1    1
#> [11,]    1    1    1
#> [12,]    1    1    1
#> [13,]    1    1    1
#> [14,]    1    1    1
#> [15,]    1    1    1
#> [16,]    1    1    1
#> [17,]    1    1    1
#> [18,]    1    1    1
#> [19,]    1    1    1
#> [20,]    1    1    1
#> [21,]    1    1    1
#> [22,]    1    1    1
#> [23,]    1    1    1
#> [24,]    1    1    1
#> 
#> $outcomes$chickenPox$pred
#>          [,1]     [,2]      [,3]      [,4]      [,5]      [,6]      [,7]
#> [1,] 49.71787  48.6989  48.20482  48.26928  48.74705  49.30286  49.51656
#> [2,] 21.65348  18.6057  16.39643  15.30112  15.34896  16.46339  18.45239
#> [3,] 96.62864 100.6954 103.39875 104.42960 103.90399 102.23375 100.03105
#>          [,8]     [,9]    [,10]    [,11]     [,12]     [,13]     [,14]
#> [1,] 49.07372 47.93041 46.30816 44.51784  42.80074  41.32312  40.24110
#> [2,] 20.88650 23.03546 24.07000 23.50624  21.53524  18.90456  16.45364
#> [3,] 98.03978 97.03413 97.62184 99.97592 103.66402 107.77232 111.30526
#>          [,15]     [,16]     [,17]     [,18]     [,19]     [,20]     [,21]
#> [1,]  39.69342  39.72013  40.19498  40.82722  41.24545  41.14590  40.42889
#> [2,]  14.73308  13.95315  14.12701  15.18255  16.94978  19.06053  20.90129
#> [3,] 113.57350 114.32672 113.67801 111.99023 109.80477 107.79357 106.66982
#>          [,22]     [,23]     [,24]
#> [1,]  39.21549  37.73478  36.21121
#> [2,]  21.78778  21.33572  19.72452
#> [3,] 106.99673 108.92951 112.06428
#> 
#> $outcomes$chickenPox$var.pred
#> , , 1
#> 
#>            [,1]      [,2]       [,3]
#> [1,]  135.57355 -24.81896 -110.75459
#> [2,]  -24.81896  73.05560  -48.23663
#> [3,] -110.75459 -48.23663  158.99122
#> 
#> , , 2
#> 
#>            [,1]      [,2]       [,3]
#> [1,]  139.71454 -21.78929 -117.92524
#> [2,]  -21.78929  66.84333  -45.05404
#> [3,] -117.92524 -45.05404  162.97928
#> 
#> , , 3
#> 
#>           [,1]      [,2]       [,3]
#> [1,]  145.4718 -19.91080 -125.56098
#> [2,]  -19.9108  62.61923  -42.70843
#> [3,] -125.5610 -42.70843  168.26941
#> 
#> , , 4
#> 
#>            [,1]      [,2]       [,3]
#> [1,]  153.61139 -19.63094 -133.98045
#> [2,]  -19.63094  62.10207  -42.47113
#> [3,] -133.98045 -42.47113  176.45159
#> 
#> , , 5
#> 
#>            [,1]      [,2]       [,3]
#> [1,]  164.37960 -21.15718 -143.22242
#> [2,]  -21.15718  66.25356  -45.09638
#> [3,] -143.22242 -45.09638  188.31880
#> 
#> , , 6
#> 
#>            [,1]      [,2]       [,3]
#> [1,]  177.33166 -24.59605 -152.73561
#> [2,]  -24.59605  75.59808  -51.00203
#> [3,] -152.73561 -51.00203  203.73764
#> 
#> , , 7
#> 
#>            [,1]      [,2]       [,3]
#> [1,]  191.39239 -29.80709 -161.58530
#> [2,]  -29.80709  90.02198  -60.21489
#> [3,] -161.58530 -60.21489  221.80020
#> 
#> , , 8
#> 
#>            [,1]      [,2]      [,3]
#> [1,]  205.03358 -36.00914 -169.0244
#> [2,]  -36.00914 107.94844  -71.9393
#> [3,] -169.02444 -71.93930  240.9637
#> 
#> , , 9
#> 
#>            [,1]      [,2]       [,3]
#> [1,]  216.66009 -41.56643 -175.09366
#> [2,]  -41.56643 125.71681  -84.15038
#> [3,] -175.09366 -84.15038  259.24404
#> 
#> , , 10
#> 
#>            [,1]      [,2]       [,3]
#> [1,]  225.18629 -44.54065 -180.64563
#> [2,]  -44.54065 138.43644  -93.89579
#> [3,] -180.64563 -93.89579  274.54142
#> 
#> , , 11
#> 
#>            [,1]      [,2]       [,3]
#> [1,]  230.34786 -43.84935 -186.49850
#> [2,]  -43.84935 142.32401  -98.47466
#> [3,] -186.49850 -98.47466  284.97316
#> 
#> , , 12
#> 
#>            [,1]      [,2]       [,3]
#> [1,]  232.68624 -40.02383 -192.66241
#> [2,]  -40.02383 136.96216  -96.93832
#> [3,] -192.66241 -96.93832  289.60073
#> 
#> , , 13
#> 
#>            [,1]      [,2]       [,3]
#> [1,]  233.73598 -34.88147 -198.85451
#> [2,]  -34.88147 125.85370  -90.97223
#> [3,] -198.85451 -90.97223  289.82674
#> 
#> , , 14
#> 
#>            [,1]      [,2]       [,3]
#> [1,]  236.08564 -30.40468 -205.68096
#> [2,]  -30.40468 114.50279  -84.09812
#> [3,] -205.68096 -84.09812  289.77907
#> 
#> , , 15
#> 
#>            [,1]      [,2]       [,3]
#> [1,]  242.57423 -27.85411 -214.72012
#> [2,]  -27.85411 107.55218  -79.69807
#> [3,] -214.72012 -79.69807  294.41818
#> 
#> , , 16
#> 
#>            [,1]      [,2]       [,3]
#> [1,]  255.12392 -27.75012 -227.37380
#> [2,]  -27.75012 107.62351  -79.87338
#> [3,] -227.37380 -79.87338  307.24718
#> 
#> , , 17
#> 
#>            [,1]      [,2]      [,3]
#> [1,]  274.09621 -30.29739 -243.7988
#> [2,]  -30.29739 115.98339  -85.6860
#> [3,] -243.79882 -85.68600  329.4848
#> 
#> , , 18
#> 
#>            [,1]      [,2]       [,3]
#> [1,]  298.22233 -35.60334 -262.61899
#> [2,]  -35.60334 133.26429  -97.66096
#> [3,] -262.61899 -97.66096  360.27995
#> 
#> , , 19
#> 
#>            [,1]       [,2]      [,3]
#> [1,]  324.77085  -43.42878 -281.3421
#> [2,]  -43.42878  159.04606 -115.6173
#> [3,] -281.34207 -115.61728  396.9593
#> 
#> , , 20
#> 
#>            [,1]       [,2]      [,3]
#> [1,]  349.88314  -52.57187 -297.3113
#> [2,]  -52.57187  190.29908 -137.7272
#> [3,] -297.31126 -137.72721  435.0385
#> 
#> , , 21
#> 
#>            [,1]       [,2]      [,3]
#> [1,]  369.39340  -60.52154 -308.8719
#> [2,]  -60.52154  220.20491 -159.6834
#> [3,] -308.87186 -159.68337  468.5552
#> 
#> , , 22
#> 
#>            [,1]       [,2]      [,3]
#> [1,]  380.22632  -64.32675 -315.8996
#> [2,]  -64.32675  239.83782 -175.5111
#> [3,] -315.89957 -175.51107  491.4106
#> 
#> , , 23
#> 
#>            [,1]       [,2]      [,3]
#> [1,]  381.61813  -62.50399 -319.1141
#> [2,]  -62.50399  242.93513 -180.4311
#> [3,] -319.11414 -180.43114  499.5453
#> 
#> , , 24
#> 
#>            [,1]       [,2]      [,3]
#> [1,]  375.59874  -56.21497 -319.3838
#> [2,]  -56.21497  230.18569 -173.9707
#> [3,] -319.38377 -173.97072  493.3545
#> 
#> 
#> $outcomes$chickenPox$icl.pred
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14]
#> [1,]   28   27   26   26   26   25   25   23   22    20    18    17    15    14
#> [2,]    8    5    4    3    3    3    4    5    6     6     5     4     3     2
#> [3,]   72   75   77   78   76   73   70   67   65    64    66    69    73    76
#>      [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24]
#> [1,]    14    13    13    12    12    11    10     9     8     7
#> [2,]     1     1     1     1     1     1     2     2     1     1
#> [3,]    78    78    76    72    68    65    62    61    62    65
#> 
#> $outcomes$chickenPox$icu.pred
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14]
#> [1,]   74   73   73   74   75   77   78   79   79    78    77    76    74    74
#> [2,]   41   37   35   34   34   37   40   45   49    51    51    49    46    42
#> [3,]  121  125  128  129  130  129  128  127  127   129   132   135   139   142
#>      [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24]
#> [1,]    74    75    76    79    81    82    83    83    82    80
#> [2,]    40    40    41    44    48    53    57    60    60    58
#> [3,]   144   145   146   146   145   145   145   146   148   150
#> 
#> 
#> $outcomes$show
#> [1] FALSE
#> 
#> 
#> $theta.mean
#>              [,1]        [,2]         [,3]        [,4]        [,5]         [,6]
#>  [1,] -1.01019636 -1.03788032 -1.065564268 -1.09324822 -1.12093217 -1.148616125
#>  [2,] -0.02768395 -0.02768395 -0.027683953 -0.02768395 -0.02768395 -0.027683953
#>  [3,] -0.05892289 -0.09537541 -0.106272170 -0.08869338 -0.04734928  0.006682028
#>  [4,] -0.08869338 -0.04734928  0.006682028  0.05892289  0.09537541  0.106272170
#>  [5,]  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000  0.000000000
#>  [6,]  0.39003624  0.39003624  0.390036237  0.39003624  0.39003624  0.390036237
#>  [7,] -2.05310860 -2.08223486 -2.111361115 -2.14048737 -2.16961363 -2.198739882
#>  [8,] -0.02912626 -0.02912626 -0.029126256 -0.02912626 -0.02912626 -0.029126256
#>  [9,]  0.03653268 -0.14360984 -0.285272227 -0.35049615 -0.32180491 -0.206886307
#> [10,] -0.35049615 -0.32180491 -0.206886307 -0.03653268  0.14360984  0.285272227
#> [11,]  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000  0.000000000
#> [12,]  0.46672032  0.46672032  0.466720323  0.46672032  0.46672032  0.466720323
#>              [,7]        [,8]         [,9]       [,10]       [,11]        [,12]
#>  [1,] -1.17630008 -1.20398403 -1.231667983 -1.25935194 -1.28703589 -1.314719841
#>  [2,] -0.02768395 -0.02768395 -0.027683953 -0.02768395 -0.02768395 -0.027683953
#>  [3,]  0.05892289  0.09537541  0.106272170  0.08869338  0.04734928 -0.006682028
#>  [4,]  0.08869338  0.04734928 -0.006682028 -0.05892289 -0.09537541 -0.106272170
#>  [5,]  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000  0.000000000
#>  [6,]  0.39003624  0.39003624  0.390036237  0.39003624  0.39003624  0.390036237
#>  [7,] -2.22786614 -2.25699239 -2.286118649 -2.31524490 -2.34437116 -2.373497416
#>  [8,] -0.02912626 -0.02912626 -0.029126256 -0.02912626 -0.02912626 -0.029126256
#>  [9,] -0.03653268  0.14360984  0.285272227  0.35049615  0.32180491  0.206886307
#> [10,]  0.35049615  0.32180491  0.206886307  0.03653268 -0.14360984 -0.285272227
#> [11,]  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000  0.000000000
#> [12,]  0.46672032  0.46672032  0.466720323  0.46672032  0.46672032  0.466720323
#>             [,13]       [,14]        [,15]       [,16]       [,17]        [,18]
#>  [1,] -1.34240379 -1.37008775 -1.397771698 -1.42545565 -1.45313960 -1.480823556
#>  [2,] -0.02768395 -0.02768395 -0.027683953 -0.02768395 -0.02768395 -0.027683953
#>  [3,] -0.05892289 -0.09537541 -0.106272170 -0.08869338 -0.04734928  0.006682028
#>  [4,] -0.08869338 -0.04734928  0.006682028  0.05892289  0.09537541  0.106272170
#>  [5,]  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000  0.000000000
#>  [6,]  0.39003624  0.39003624  0.390036237  0.39003624  0.39003624  0.390036237
#>  [7,] -2.40262367 -2.43174993 -2.460876183 -2.49000244 -2.51912869 -2.548254950
#>  [8,] -0.02912626 -0.02912626 -0.029126256 -0.02912626 -0.02912626 -0.029126256
#>  [9,]  0.03653268 -0.14360984 -0.285272227 -0.35049615 -0.32180491 -0.206886307
#> [10,] -0.35049615 -0.32180491 -0.206886307 -0.03653268  0.14360984  0.285272227
#> [11,]  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000  0.000000000
#> [12,]  0.46672032  0.46672032  0.466720323  0.46672032  0.46672032  0.466720323
#>             [,19]       [,20]        [,21]       [,22]       [,23]        [,24]
#>  [1,] -1.50850751 -1.53619146 -1.563875413 -1.59155937 -1.61924332 -1.646927271
#>  [2,] -0.02768395 -0.02768395 -0.027683953 -0.02768395 -0.02768395 -0.027683953
#>  [3,]  0.05892289  0.09537541  0.106272170  0.08869338  0.04734928 -0.006682028
#>  [4,]  0.08869338  0.04734928 -0.006682028 -0.05892289 -0.09537541 -0.106272170
#>  [5,]  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000  0.000000000
#>  [6,]  0.39003624  0.39003624  0.390036237  0.39003624  0.39003624  0.390036237
#>  [7,] -2.57738121 -2.60650746 -2.635633716 -2.66475997 -2.69388623 -2.723012483
#>  [8,] -0.02912626 -0.02912626 -0.029126256 -0.02912626 -0.02912626 -0.029126256
#>  [9,] -0.03653268  0.14360984  0.285272227  0.35049615  0.32180491  0.206886307
#> [10,]  0.35049615  0.32180491  0.206886307  0.03653268 -0.14360984 -0.285272227
#> [11,]  0.00000000  0.00000000  0.000000000  0.00000000  0.00000000  0.000000000
#> [12,]  0.46672032  0.46672032  0.466720323  0.46672032  0.46672032  0.466720323
#> 
#> $theta.cov
#> , , 1
#> 
#>                [,1]          [,2]          [,3]          [,4]      [,5]
#>  [1,]  9.559869e-02  1.868326e-03  2.177845e-04  2.520295e-03 0.0000000
#>  [2,]  1.868326e-03  4.604931e-04 -3.709447e-05  1.628894e-04 0.0000000
#>  [3,]  2.177845e-04 -3.709447e-05  6.451331e-03  1.343320e-04 0.0000000
#>  [4,]  2.520295e-03  1.628894e-04  1.343320e-04  6.562663e-03 0.0000000
#>  [5,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 0.0716381
#>  [6,] -6.199776e-02 -8.659647e-07 -6.660636e-04  6.198757e-04 0.0000000
#>  [7,]  5.740983e-03  1.057553e-04  1.793818e-04 -4.734687e-06 0.0000000
#>  [8,]  1.053915e-04  5.737582e-06  4.118280e-06  3.150242e-06 0.0000000
#>  [9,]  1.844940e-04  4.187968e-06  4.309714e-04  2.915415e-06 0.0000000
#> [10,]  2.845902e-06  3.353556e-06  2.289083e-06  4.250154e-04 0.0000000
#> [11,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 0.0000000
#> [12,] -3.395756e-03  1.197936e-06 -1.080394e-04  1.108729e-04 0.0000000
#>                [,6]          [,7]          [,8]          [,9]        [,10]
#>  [1,] -6.199776e-02  5.740983e-03  1.053915e-04  1.844940e-04 2.845902e-06
#>  [2,] -8.659647e-07  1.057553e-04  5.737582e-06  4.187968e-06 3.353556e-06
#>  [3,] -6.660636e-04  1.793818e-04  4.118280e-06  4.309714e-04 2.289083e-06
#>  [4,]  6.198757e-04 -4.734687e-06  3.150242e-06  2.915415e-06 4.250154e-04
#>  [5,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 0.000000e+00
#>  [6,]  6.223498e-02 -3.393720e-03  1.244983e-06 -1.125074e-04 1.144439e-04
#>  [7,] -3.393720e-03  1.092513e-01  2.228886e-03  2.307867e-04 3.328256e-03
#>  [8,]  1.244983e-06  2.228886e-03  5.406644e-04 -4.712594e-05 1.967476e-04
#>  [9,] -1.125074e-04  2.307867e-04 -4.712594e-05  8.001389e-03 1.299627e-04
#> [10,]  1.144439e-04  3.328256e-03  1.967476e-04  1.299627e-04 8.234927e-03
#> [11,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 0.000000e+00
#> [12,]  3.426467e-03 -6.900170e-02  1.662161e-07 -8.172003e-04 7.844873e-04
#>            [,11]         [,12]
#>  [1,] 0.00000000 -3.395756e-03
#>  [2,] 0.00000000  1.197936e-06
#>  [3,] 0.00000000 -1.080394e-04
#>  [4,] 0.00000000  1.108729e-04
#>  [5,] 0.00000000  0.000000e+00
#>  [6,] 0.00000000  3.426467e-03
#>  [7,] 0.00000000 -6.900170e-02
#>  [8,] 0.00000000  1.662161e-07
#>  [9,] 0.00000000 -8.172003e-04
#> [10,] 0.00000000  7.844873e-04
#> [11,] 0.06747194  0.000000e+00
#> [12,] 0.00000000  6.933184e-02
#> 
#> , , 2
#> 
#>                [,1]          [,2]          [,3]          [,4]      [,5]
#>  [1,]  1.045778e-01  2.328819e-03  1.498074e-03  2.233361e-03 0.0000000
#>  [2,]  2.328819e-03  4.835345e-04  4.931995e-05  1.596136e-04 0.0000000
#>  [3,]  1.498074e-03  4.931995e-05  6.756898e-03  1.153742e-04 0.0000000
#>  [4,]  2.233361e-03  1.596136e-04  1.153742e-04  6.582627e-03 0.0000000
#>  [5,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 0.0716381
#>  [6,] -6.199862e-02 -8.659647e-07 -2.668902e-04  8.698599e-04 0.0000000
#>  [7,]  5.957867e-03  1.114929e-04  1.581235e-04 -9.312221e-05 0.0000000
#>  [8,]  1.111291e-04  5.737582e-06  5.141656e-06  6.690499e-07 0.0000000
#>  [9,]  1.665031e-04  5.303664e-06  4.317360e-04 -9.647155e-07 0.0000000
#> [10,] -8.897209e-05  8.102811e-07 -1.591048e-06  4.242508e-04 0.0000000
#> [11,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 0.0000000
#> [12,] -3.394558e-03  1.197936e-06 -3.812842e-05  1.500384e-04 0.0000000
#>                [,6]          [,7]         [,8]          [,9]         [,10]
#>  [1,] -6.199862e-02  5.957867e-03 1.111291e-04  1.665031e-04 -8.897209e-05
#>  [2,] -8.659647e-07  1.114929e-04 5.737582e-06  5.303664e-06  8.102811e-07
#>  [3,] -2.668902e-04  1.581235e-04 5.141656e-06  4.317360e-04 -1.591048e-06
#>  [4,]  8.698599e-04 -9.312221e-05 6.690499e-07 -9.647155e-07  4.242508e-04
#>  [5,]  0.000000e+00  0.000000e+00 0.000000e+00  0.000000e+00  0.000000e+00
#>  [6,]  6.223498e-02 -3.392475e-03 1.244983e-06 -4.021231e-05  1.553650e-04
#>  [7,] -3.392475e-03  1.197175e-01 2.769550e-03  1.921557e-03  2.960913e-03
#>  [8,]  1.244983e-06  2.769550e-03 5.677054e-04  5.756152e-05  1.939514e-04
#>  [9,] -4.021231e-05  1.921557e-03 5.756152e-05  8.373676e-03  1.661064e-04
#> [10,]  1.553650e-04  2.960913e-03 1.939514e-04  1.661064e-04  8.268741e-03
#> [11,]  0.000000e+00  0.000000e+00 0.000000e+00  0.000000e+00  0.000000e+00
#> [12,]  3.426467e-03 -6.900154e-02 1.662161e-07 -3.154726e-04  1.087986e-03
#>            [,11]         [,12]
#>  [1,] 0.00000000 -3.394558e-03
#>  [2,] 0.00000000  1.197936e-06
#>  [3,] 0.00000000 -3.812842e-05
#>  [4,] 0.00000000  1.500384e-04
#>  [5,] 0.00000000  0.000000e+00
#>  [6,] 0.00000000  3.426467e-03
#>  [7,] 0.00000000 -6.900154e-02
#>  [8,] 0.00000000  1.662161e-07
#>  [9,] 0.00000000 -3.154726e-04
#> [10,] 0.00000000  1.087986e-03
#> [11,] 0.06747194  0.000000e+00
#> [12,] 0.00000000  6.933184e-02
#> 
#> , , 3
#> 
#>                [,1]          [,2]          [,3]          [,4]      [,5]
#>  [1,]  0.1145008962  2.812354e-03  2.536570e-03  1.298679e-03 0.0000000
#>  [2,]  0.0028123539  5.065759e-04  1.225191e-04  1.135695e-04 0.0000000
#>  [3,]  0.0025365700  1.225191e-04  6.974647e-03 -1.777445e-05 0.0000000
#>  [4,]  0.0012986794  1.135695e-04 -1.777445e-05  6.690410e-03 0.0000000
#>  [5,]  0.0000000000  0.000000e+00  0.000000e+00  0.000000e+00 0.0716381
#>  [6,] -0.0619994881 -8.659647e-07  2.037963e-04  8.867658e-04 0.0000000
#>  [7,]  0.0061862268  1.172305e-04  9.516520e-05 -1.616994e-04 0.0000000
#>  [8,]  0.0001168667  5.737582e-06  4.787330e-06 -1.991414e-06 0.0000000
#>  [9,]  0.0001047081  4.998249e-06  4.287580e-04 -3.566964e-06 0.0000000
#> [10,] -0.0001622538 -1.950108e-06 -4.193297e-06  4.272288e-04 0.0000000
#> [11,]  0.0000000000  0.000000e+00  0.000000e+00  0.000000e+00 0.0000000
#> [12,] -0.0033933606  1.197936e-06  4.199903e-05  1.490013e-04 0.0000000
#>                [,6]          [,7]          [,8]          [,9]         [,10]
#>  [1,] -6.199949e-02  0.0061862268  1.168667e-04  1.047081e-04 -1.622538e-04
#>  [2,] -8.659647e-07  0.0001172305  5.737582e-06  4.998249e-06 -1.950108e-06
#>  [3,]  2.037963e-04  0.0000951652  4.787330e-06  4.287580e-04 -4.193297e-06
#>  [4,]  8.867658e-04 -0.0001616994 -1.991414e-06 -3.566964e-06  4.272288e-04
#>  [5,]  0.000000e+00  0.0000000000  0.000000e+00  0.000000e+00  0.000000e+00
#>  [6,]  6.223498e-02 -0.0033912296  1.244983e-06  4.285764e-05  1.546562e-04
#>  [7,] -3.391230e-03  0.1312919874  3.337256e-03  3.291399e-03  1.742633e-03
#>  [8,]  1.244983e-06  0.0033372555  5.947465e-04  1.468254e-04  1.391860e-04
#>  [9,]  4.285764e-05  0.0032913988  1.468254e-04  8.692646e-03  3.761480e-05
#> [10,]  1.546562e-04  0.0017426332  1.391860e-04  3.761480e-05  8.355871e-03
#> [11,]  0.000000e+00  0.0000000000  0.000000e+00  0.000000e+00  0.000000e+00
#> [12,]  3.426467e-03 -0.0690013708  1.662161e-07  2.707857e-04  1.099960e-03
#>            [,11]         [,12]
#>  [1,] 0.00000000 -3.393361e-03
#>  [2,] 0.00000000  1.197936e-06
#>  [3,] 0.00000000  4.199903e-05
#>  [4,] 0.00000000  1.490013e-04
#>  [5,] 0.00000000  0.000000e+00
#>  [6,] 0.00000000  3.426467e-03
#>  [7,] 0.00000000 -6.900137e-02
#>  [8,] 0.00000000  1.662161e-07
#>  [9,] 0.00000000  2.707857e-04
#> [10,] 0.00000000  1.099960e-03
#> [11,] 0.06747194  0.000000e+00
#> [12,] 0.00000000  6.933184e-02
#> 
#> , , 4
#> 
#>                [,1]          [,2]          [,3]          [,4]      [,5]
#>  [1,]  1.254141e-01  3.318930e-03  3.008963e-03 -1.065011e-04 0.0000000
#>  [2,]  3.318930e-03  5.296173e-04  1.628894e-04  3.709447e-05 0.0000000
#>  [3,]  3.008963e-03  1.628894e-04  7.049594e-03 -1.319653e-04 0.0000000
#>  [4,] -1.065011e-04  3.709447e-05 -1.319653e-04  6.940995e-03 0.0000000
#>  [5,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 0.0716381
#>  [6,] -6.200035e-02 -8.659647e-07  6.198757e-04  6.660636e-04 0.0000000
#>  [7,]  6.426062e-03  1.229681e-04  4.716040e-06 -1.917366e-04 0.0000000
#>  [8,]  1.226043e-04  5.737582e-06  3.150242e-06 -4.118280e-06 0.0000000
#>  [9,]  1.290657e-05  3.353556e-06  4.250154e-04 -2.289083e-06 0.0000000
#> [10,] -1.970579e-04 -4.187968e-06 -2.915415e-06  4.309714e-04 0.0000000
#> [11,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 0.0000000
#> [12,] -3.392163e-03  1.197936e-06  1.108729e-04  1.080394e-04 0.0000000
#>                [,6]          [,7]          [,8]          [,9]         [,10]
#>  [1,] -6.200035e-02  6.426062e-03  1.226043e-04  1.290657e-05 -1.970579e-04
#>  [2,] -8.659647e-07  1.229681e-04  5.737582e-06  3.353556e-06 -4.187968e-06
#>  [3,]  6.198757e-04  4.716040e-06  3.150242e-06  4.250154e-04 -2.915415e-06
#>  [4,]  6.660636e-04 -1.917366e-04 -4.118280e-06 -2.289083e-06  4.309714e-04
#>  [5,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
#>  [6,]  6.223498e-02 -3.389985e-03  1.244983e-06  1.144439e-04  1.125074e-04
#>  [7,] -3.389985e-03  1.440290e-01  3.932002e-03  3.918499e-03 -8.940886e-05
#>  [8,]  1.244983e-06  3.932002e-03  6.217875e-04  1.967476e-04  4.712594e-05
#>  [9,]  1.144439e-04  3.918499e-03  1.967476e-04  8.842380e-03 -1.270205e-04
#> [10,]  1.125074e-04 -8.940886e-05  4.712594e-05 -1.270205e-04  8.612239e-03
#> [11,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
#> [12,]  3.426467e-03 -6.900120e-02  1.662161e-07  7.844873e-04  8.172003e-04
#>            [,11]         [,12]
#>  [1,] 0.00000000 -3.392163e-03
#>  [2,] 0.00000000  1.197936e-06
#>  [3,] 0.00000000  1.108729e-04
#>  [4,] 0.00000000  1.080394e-04
#>  [5,] 0.00000000  0.000000e+00
#>  [6,] 0.00000000  3.426467e-03
#>  [7,] 0.00000000 -6.900120e-02
#>  [8,] 0.00000000  1.662161e-07
#>  [9,] 0.00000000  7.844873e-04
#> [10,] 0.00000000  8.172003e-04
#> [11,] 0.06747194  0.000000e+00
#> [12,] 0.00000000  6.933184e-02
#> 
#> , , 5
#> 
#>                [,1]          [,2]          [,3]          [,4]      [,5]
#>  [1,]  1.373635e-01  3.848547e-03  2.712202e-03 -1.646034e-03 0.0000000
#>  [2,]  3.848547e-03  5.526587e-04  1.596136e-04 -4.931995e-05 0.0000000
#>  [3,]  2.712202e-03  1.596136e-04  7.069558e-03 -1.130075e-04 0.0000000
#>  [4,] -1.646034e-03 -4.931995e-05 -1.130075e-04  7.246562e-03 0.0000000
#>  [5,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 0.0716381
#>  [6,] -6.200122e-02 -8.659647e-07  8.698599e-04  2.668902e-04 0.0000000
#>  [7,]  6.677371e-03  1.287056e-04 -9.111506e-05 -1.735485e-04 0.0000000
#>  [8,]  1.283419e-04  5.737582e-06  6.690499e-07 -5.141656e-06 0.0000000
#>  [9,] -8.654125e-05  8.102811e-07  4.242508e-04  1.591048e-06 0.0000000
#> [10,] -1.824141e-04 -5.303664e-06  9.647155e-07  4.317360e-04 0.0000000
#> [11,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 0.0000000
#> [12,] -3.390965e-03  1.197936e-06  1.500384e-04  3.812842e-05 0.0000000
#>                [,6]          [,7]          [,8]          [,9]         [,10]
#>  [1,] -6.200122e-02  6.677371e-03  1.283419e-04 -8.654125e-05 -1.824141e-04
#>  [2,] -8.659647e-07  1.287056e-04  5.737582e-06  8.102811e-07 -5.303664e-06
#>  [3,]  8.698599e-04 -9.111506e-05  6.690499e-07  4.242508e-04  9.647155e-07
#>  [4,]  2.668902e-04 -1.735485e-04 -5.141656e-06  1.591048e-06  4.317360e-04
#>  [5,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
#>  [6,]  6.223498e-02 -3.388740e-03  1.244983e-06  1.553650e-04  4.021231e-05
#>  [7,] -3.388740e-03  1.579825e-01  4.553790e-03  3.542767e-03 -2.094241e-03
#>  [8,]  1.244983e-06  4.553790e-03  6.488286e-04  1.939514e-04 -5.756152e-05
#>  [9,]  1.553650e-04  3.542767e-03  1.939514e-04  8.876193e-03 -1.631642e-04
#> [10,]  4.021231e-05 -2.094241e-03 -5.756152e-05 -1.631642e-04  8.984526e-03
#> [11,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
#> [12,]  3.426467e-03 -6.900104e-02  1.662161e-07  1.087986e-03  3.154726e-04
#>            [,11]         [,12]
#>  [1,] 0.00000000 -3.390965e-03
#>  [2,] 0.00000000  1.197936e-06
#>  [3,] 0.00000000  1.500384e-04
#>  [4,] 0.00000000  3.812842e-05
#>  [5,] 0.00000000  0.000000e+00
#>  [6,] 0.00000000  3.426467e-03
#>  [7,] 0.00000000 -6.900104e-02
#>  [8,] 0.00000000  1.662161e-07
#>  [9,] 0.00000000  1.087986e-03
#> [10,] 0.00000000  3.154726e-04
#> [11,] 0.06747194  0.000000e+00
#> [12,] 0.00000000  6.933184e-02
#> 
#> , , 6
#> 
#>                [,1]          [,2]          [,3]          [,4]      [,5]
#>  [1,]  0.1503952412  4.401206e-03  1.639388e-03 -2.904127e-03 0.0000000
#>  [2,]  0.0044012059  5.757001e-04  1.135695e-04 -1.225191e-04 0.0000000
#>  [3,]  0.0016393878  1.135695e-04  7.177341e-03  2.014110e-05 0.0000000
#>  [4,] -0.0029041274 -1.225191e-04  2.014110e-05  7.464311e-03 0.0000000
#>  [5,]  0.0000000000  0.000000e+00  0.000000e+00  0.000000e+00 0.0716381
#>  [6,] -0.0620020860 -8.659647e-07  8.867658e-04 -2.037963e-04 0.0000000
#>  [7,]  0.0069401565  1.344432e-04 -1.676736e-04 -1.095272e-04 0.0000000
#>  [8,]  0.0001340794  5.737582e-06 -1.991414e-06 -4.787330e-06 0.0000000
#>  [9,] -0.0001681041 -1.950108e-06  4.272288e-04  4.193297e-06 0.0000000
#> [10,] -0.0001197029 -4.998249e-06  3.566964e-06  4.287580e-04 0.0000000
#> [11,]  0.0000000000  0.000000e+00  0.000000e+00  0.000000e+00 0.0000000
#> [12,] -0.0033897667  1.197936e-06  1.490013e-04 -4.199903e-05 0.0000000
#>                [,6]          [,7]          [,8]          [,9]         [,10]
#>  [1,] -6.200209e-02  0.0069401565  1.340794e-04 -1.681041e-04 -1.197029e-04
#>  [2,] -8.659647e-07  0.0001344432  5.737582e-06 -1.950108e-06 -4.998249e-06
#>  [3,]  8.867658e-04 -0.0001676736 -1.991414e-06  4.272288e-04  3.566964e-06
#>  [4,] -2.037963e-04 -0.0001095272 -4.787330e-06  4.193297e-06  4.287580e-04
#>  [5,]  0.000000e+00  0.0000000000  0.000000e+00  0.000000e+00  0.000000e+00
#>  [6,]  6.223498e-02 -0.0033874947  1.244983e-06  1.546562e-04 -4.285764e-05
#>  [7,] -3.387495e-03  0.1732065904  5.202618e-03  2.160191e-03 -3.731875e-03
#>  [8,]  1.244983e-06  0.0052026181  6.758696e-04  1.391860e-04 -1.468254e-04
#>  [9,]  1.546562e-04  0.0021601913  1.391860e-04  8.963324e-03 -3.467260e-05
#> [10,] -4.285764e-05 -0.0037318751 -1.468254e-04 -3.467260e-05  9.303497e-03
#> [11,]  0.000000e+00  0.0000000000  0.000000e+00  0.000000e+00  0.000000e+00
#> [12,]  3.426467e-03 -0.0690008722  1.662161e-07  1.099960e-03 -2.707857e-04
#>            [,11]         [,12]
#>  [1,] 0.00000000 -3.389767e-03
#>  [2,] 0.00000000  1.197936e-06
#>  [3,] 0.00000000  1.490013e-04
#>  [4,] 0.00000000 -4.199903e-05
#>  [5,] 0.00000000  0.000000e+00
#>  [6,] 0.00000000  3.426467e-03
#>  [7,] 0.00000000 -6.900087e-02
#>  [8,] 0.00000000  1.662161e-07
#>  [9,] 0.00000000  1.099960e-03
#> [10,] 0.00000000 -2.707857e-04
#> [11,] 0.06747194  0.000000e+00
#> [12,] 0.00000000  6.933184e-02
#> 
#> , , 7
#> 
#>                [,1]          [,2]          [,3]          [,4]      [,5]
#>  [1,]  1.645553e-01  4.976906e-03  4.782306e-06 -3.497631e-03 0.0000000
#>  [2,]  4.976906e-03  5.987415e-04  3.709447e-05 -1.628894e-04 0.0000000
#>  [3,]  4.782306e-06  3.709447e-05  7.427926e-03  1.343320e-04 0.0000000
#>  [4,] -3.497631e-03 -1.628894e-04  1.343320e-04  7.539258e-03 0.0000000
#>  [5,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 0.0716381
#>  [6,] -6.200295e-02 -8.659647e-07  6.660636e-04 -6.198757e-04 0.0000000
#>  [7,]  7.214417e-03  1.401808e-04 -2.040915e-04 -1.416677e-05 0.0000000
#>  [8,]  1.398170e-04  5.737582e-06 -4.118280e-06 -3.150242e-06 0.0000000
#>  [9,] -2.096218e-04 -4.187968e-06  4.309714e-04  2.915415e-06 0.0000000
#> [10,] -2.296724e-05 -3.353556e-06  2.289083e-06  4.250154e-04 0.0000000
#> [11,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 0.0000000
#>                [,6]          [,7]          [,8]          [,9]         [,10]
#>  [1,] -6.200295e-02  7.214417e-03  1.398170e-04 -2.096218e-04 -2.296724e-05
#>  [2,] -8.659647e-07  1.401808e-04  5.737582e-06 -4.187968e-06 -3.353556e-06
#>  [3,]  6.660636e-04 -2.040915e-04 -4.118280e-06  4.309714e-04  2.289083e-06
#>  [4,] -6.198757e-04 -1.416677e-05 -3.150242e-06  2.915415e-06  4.250154e-04
#>  [5,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
#>  [6,]  6.223498e-02 -3.386250e-03  1.244983e-06  1.125074e-04 -1.144439e-04
#>  [7,] -3.386250e-03  1.897554e-01  5.878488e-03  5.196896e-05 -4.508742e-03
#>  [8,]  1.244983e-06  5.878488e-03  7.029107e-04  4.712594e-05 -1.967476e-04
#>  [9,]  1.125074e-04  5.196896e-05  4.712594e-05  9.219692e-03  1.299627e-04
#> [10,] -1.144439e-04 -4.508742e-03 -1.967476e-04  1.299627e-04  9.453230e-03
#> [11,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
#>            [,11]         [,12]
#>  [1,] 0.00000000 -3.388569e-03
#>  [2,] 0.00000000  1.197936e-06
#>  [3,] 0.00000000  1.080394e-04
#>  [4,] 0.00000000 -1.108729e-04
#>  [5,] 0.00000000  0.000000e+00
#>  [6,] 0.00000000  3.426467e-03
#>  [7,] 0.00000000 -6.900071e-02
#>  [8,] 0.00000000  1.662161e-07
#>  [9,] 0.00000000  8.172003e-04
#> [10,] 0.00000000 -7.844873e-04
#> [11,] 0.06747194  0.000000e+00
#> 
#>  [ reached 'max' / getOption("max.print") -- omitted 17 slices ] 
#> 
#> $lambda.mean
#>           [,1]       [,2]       [,3]       [,4]       [,5]       [,6]
#> [1,] -0.679083 -0.7432195 -0.7818002 -0.7919054 -0.7782452 -0.7518979
#> [2,] -1.549856 -1.7591244 -1.9299130 -2.0242632 -2.0246982 -1.9389059
#>            [,7]       [,8]       [,9]      [,10]      [,11]      [,12]
#> [1,] -0.7273409 -0.7185724 -0.7353596 -0.7806223 -0.8496504 -0.9313656
#> [2,] -1.7976785 -1.6466622 -1.5341261 -1.4980284 -1.5558459 -1.6998908
#>          [,13]     [,14]     [,15]     [,16]     [,17]     [,18]     [,19]
#> [1,] -1.011290 -1.075427 -1.114008 -1.124113 -1.110453 -1.084105 -1.059548
#> [2,] -1.899371 -2.108639 -2.279428 -2.373778 -2.374213 -2.288421 -2.147194
#>          [,20]     [,21]     [,22]     [,23]     [,24]
#> [1,] -1.050780 -1.067567 -1.112830 -1.181858 -1.263573
#> [2,] -1.996177 -1.883641 -1.847543 -1.905361 -2.049406
#> 
#> $lambda.cov
#> , , 1
#> 
#>             [,1]        [,2]
#> [1,] 0.111031031 0.002952275
#> [2,] 0.002952275 0.114880251
#> 
#> , , 2
#> 
#>             [,1]        [,2]
#> [1,] 0.123672883 0.003275323
#> [2,] 0.003275323 0.130104017
#> 
#> , , 3
#> 
#>             [,1]        [,2]
#> [1,] 0.136830380 0.003541592
#> [2,] 0.003541592 0.145910042
#> 
#> , , 4
#> 
#>             [,1]        [,2]
#> [1,] 0.149593768 0.003738337
#> [2,] 0.003738337 0.161078684
#> 
#> , , 5
#> 
#>             [,1]        [,2]
#> [1,] 0.161467867 0.003876132
#> [2,] 0.003876132 0.174921870
#> 
#> , , 6
#> 
#>             [,1]        [,2]
#> [1,] 0.172493799 0.003984471
#> [2,] 0.003984471 0.187492253
#> 
#> , , 7
#> 
#>             [,1]        [,2]
#> [1,] 0.183192092 0.004103871
#> [2,] 0.004103871 0.199515810
#> 
#> , , 8
#> 
#>             [,1]        [,2]
#> [1,] 0.194434525 0.004277022
#> [2,] 0.004277022 0.212174787
#> 
#> , , 9
#> 
#>             [,1]        [,2]
#> [1,] 0.207308809 0.004539212
#> [2,] 0.004539212 0.226871186
#> 
#> , , 10
#> 
#>             [,1]        [,2]
#> [1,] 0.222942264 0.004908079
#> [2,] 0.004908079 0.244987991
#> 
#> , , 11
#> 
#>             [,1]        [,2]
#> [1,] 0.242227680 0.005375258
#> [2,] 0.005375258 0.267588593
#> 
#> , , 12
#> 
#>             [,1]        [,2]
#> [1,] 0.265491167 0.005904707
#> [2,] 0.005904707 0.295057987
#> 
#> , , 13
#> 
#>             [,1]        [,2]
#> [1,] 0.292266282 0.006441238
#> [2,] 0.006441238 0.326834106
#> 
#> , , 14
#> 
#>             [,1]        [,2]
#> [1,] 0.321351876 0.006927658
#> [2,] 0.006927658 0.361440228
#> 
#> , , 15
#> 
#>             [,1]        [,2]
#> [1,] 0.351188944 0.007323712
#> [2,] 0.007323712 0.396907427
#> 
#> , , 16
#> 
#>             [,1]        [,2]
#> [1,] 0.380397001 0.007618777
#> [2,] 0.007618777 0.431442027
#> 
#> , , 17
#> 
#>             [,1]        [,2]
#> [1,] 0.408221257 0.007833981
#> [2,] 0.007833981 0.464034915
#> 
#> , , 18
#> 
#>             [,1]        [,2]
#> [1,] 0.434723899 0.008014972
#> [2,] 0.008014972 0.494756726
#> 
#> , , 19
#> 
#>             [,1]        [,2]
#> [1,] 0.460721555 0.008219696
#> [2,] 0.008219696 0.524685622
#> 
#> , , 20
#> 
#>             [,1]        [,2]
#> [1,] 0.487577799 0.008504881
#> [2,] 0.008504881 0.555595865
#> 
#> , , 21
#> 
#>            [,1]       [,2]
#> [1,] 0.51693605 0.00891269
#> [2,] 0.00891269 0.58956268
#> 
#> , , 22
#> 
#>            [,1]       [,2]
#> [1,] 0.55039436 0.00945864
#> [2,] 0.00945864 0.62854309
#> 
#> , , 23
#> 
#>            [,1]       [,2]
#> [1,] 0.58910514 0.01012381
#> [2,] 0.01012381 0.67392152
#> 
#> , , 24
#> 
#>            [,1]       [,2]
#> [1,] 0.63337341 0.01085602
#> [2,] 0.01085602 0.72606499
#> 
#> 
#> $plot
#> 
```
