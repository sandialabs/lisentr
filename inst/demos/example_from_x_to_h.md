Understanding the Conversion from `x` to `h` within \`listenr
================
<br> Updated on March 21, 2024

Below is an example using `listenr` to fit an echo-state network (ESN)
model to the MERRA-2 data included in the package (`merra2`). The model
is trained to forecast stratospheric temperatures given lagged
stratospheric temperatures and lagged AOD values on monthly data from
1986 to 1993.

Load R packages:

``` r
library(dplyr)
library(listenr)
library(tidyr)
```

## Model Training

Prepare training data for ESN (only AOD and stratospheric temperature):

``` r
merra2_train = merra2 %>% filter(data == "train")

train_mat_aod <- 
  merra2_train %>%
  select(date, lon, lat, aod_stdzd) %>%
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = date,
    values_from = aod_stdzd
  ) %>%
  select(-lon,-lat) %>%
  t()

train_mat_temp_strat <- 
  merra2_train %>%
  select(date, lon, lat, temp_strat_stdzd) %>%
  pivot_wider(
    id_cols = c(lon, lat),
    names_from = date,
    values_from = temp_strat_stdzd
  ) %>%
  select(-lon,-lat) %>%
  t()
```

Compute EOFs:

``` r
eofs_aod = compute_eofs(Ztrain = train_mat_aod, n_eofs = 10)
eofs_temp_strat = compute_eofs(Ztrain = train_mat_temp_strat, n_eofs = 10)
```

Specify model inputs/outputs (AOD and stratospheric temperature are used
to forecast stratospheric temperature):

``` r
x = cbind(eofs_aod$train, eofs_temp_strat$train)
y = eofs_temp_strat$train
```

Extract the times associated with the training data (note that the times
for `x` and `y` are the same):

``` r
t = rownames(x)
```

Train the ESN:

``` r
esn <-
  fit_esn(
    x = x,
    y = y,
    t = as.character(t),
    tau = 1,
    m = 5,
    tau_emb = 1,
    nh = 50,
    add_quad = TRUE, 
    internal_scaling = "joint",
    seed = 20220915
  )
```

## From `x` to `h`

Explaining how the input data changes as the ESN is trained:

1.  Start with matrix of 20 columns containing specified input variables
    (in this case PCs 1-10 for AOD and stratospheric temperature) and 96
    rows containing the number of times for training:

    ``` r
    # Dimensions
    dim(esn$data_input$x)
    ## [1] 96 20

    # First 6 rows and columns 
    head(esn$data_input$x)[, 1:6]
    ##                  [,1]       [,2]       [,3]       [,4]       [,5]      [,6]
    ## 1986-01-01  -9.813371 -2.7221423 -2.9670116  1.1736340  4.4132166 0.1832662
    ## 1986-02-01 -11.231271 -3.9372459 -1.8107796  0.6145001  2.8378613 1.4372247
    ## 1986-03-01 -10.569711 -1.0830847 -2.9901863 -0.9537097  2.1107133 1.8808492
    ## 1986-04-01 -10.389466 -0.5819202 -2.3548835 -2.7844826  0.1056463 2.9105547
    ## 1986-05-01 -10.679590 -0.9002462 -2.4194503 -3.6795298 -1.6297185 2.0407384
    ## 1986-06-01  -9.612050 -2.1341306  0.2575888 -4.1924677 -2.3984229 0.3436124

    # Last 5 rows and columns 
    tail(esn$data_input$x)[, 1:6]
    ##                 [,1]        [,2]       [,3]      [,4]       [,5]        [,6]
    ## 1993-07-01 34.360066 -12.4796463  3.5104701 -3.614514  0.2479648 -1.03801381
    ## 1993-08-01 28.204443 -10.4497726  4.9715430  1.160067 -1.0031331 -0.87301575
    ## 1993-09-01 21.432044  -8.2172735  5.7101922  6.553072 -2.3650500 -1.17995373
    ## 1993-10-01 13.114014  -3.4531962  3.4714111  8.003122 -1.9213596 -1.17449153
    ## 1993-11-01  4.557627   0.0524969 -0.9472726  4.504729 -0.2122131  0.48301545
    ## 1993-12-01  2.122742   1.5795372 -3.8945770  2.482843  3.1886899 -0.08174707
    ```

2.  Then some rows at the end of the input data matrix are removed based
    on the forecast lag specified (note that one row has been removed
    here since `tau` is 1:

    ``` r
    # Dimensions
    dim(esn$data_train$x_train)
    ## [1] 95 20

    # First 6 rows and columns 
    head(esn$data_train$x_train)[,1:6]
    ##                  [,1]       [,2]       [,3]       [,4]       [,5]      [,6]
    ## 1986-01-01  -9.813371 -2.7221423 -2.9670116  1.1736340  4.4132166 0.1832662
    ## 1986-02-01 -11.231271 -3.9372459 -1.8107796  0.6145001  2.8378613 1.4372247
    ## 1986-03-01 -10.569711 -1.0830847 -2.9901863 -0.9537097  2.1107133 1.8808492
    ## 1986-04-01 -10.389466 -0.5819202 -2.3548835 -2.7844826  0.1056463 2.9105547
    ## 1986-05-01 -10.679590 -0.9002462 -2.4194503 -3.6795298 -1.6297185 2.0407384
    ## 1986-06-01  -9.612050 -2.1341306  0.2575888 -4.1924677 -2.3984229 0.3436124

    # Last 6 rows and columns 
    tail(esn$data_train$x_train)[,1:6]
    ##                 [,1]        [,2]       [,3]      [,4]       [,5]       [,6]
    ## 1993-06-01 23.606502  -7.7439838  0.7749936 -2.106507 -0.7395371  0.8087662
    ## 1993-07-01 34.360066 -12.4796463  3.5104701 -3.614514  0.2479648 -1.0380138
    ## 1993-08-01 28.204443 -10.4497726  4.9715430  1.160067 -1.0031331 -0.8730158
    ## 1993-09-01 21.432044  -8.2172735  5.7101922  6.553072 -2.3650500 -1.1799537
    ## 1993-10-01 13.114014  -3.4531962  3.4714111  8.003122 -1.9213596 -1.1744915
    ## 1993-11-01  4.557627   0.0524969 -0.9472726  4.504729 -0.2122131  0.4830154
    ```

3.  Next, the embedding vectors are created and converted into a “design
    matrix”. The first column contains all 1s to represent the
    intercept, and the other columns contain various lags of the input
    variables. In this case, the times included in the embedding vector
    are $t$, $t-1$, $t-2$, $t-3$, $t-4$, and $t-5$). Note that these
    times will be included in the design matrix for each input variable.
    In this case, the design matrix has 121 columns: 1 for the intercept
    and 120 for the 20 principal components times the 6 time points
    included (1+(20\*6)). The number of rows (times) is reduced further
    based on the embedding vector length (`m`) and embedding vector lag
    (`tau_emb`). Here, the PCs have also been centered and scaled
    internally before creating the design matrix:

    ``` r
    # Dimensions
    dim(esn$data_train$design_matrix)
    ## [1]  90 121

    # First 6 rows and columns
    head(esn$data_train$design_matrix)[,1:6]
    ##      [,1]      [,2]      [,3]      [,4]      [,5]      [,6]
    ## [1,]    1 -1.934558 -2.214181 -2.083715 -2.048169 -2.105384
    ## [2,]    1 -2.214181 -2.083715 -2.048169 -2.105384 -1.894856
    ## [3,]    1 -2.083715 -2.048169 -2.105384 -1.894856 -1.634303
    ## [4,]    1 -2.048169 -2.105384 -1.894856 -1.634303 -1.663443
    ## [5,]    1 -2.105384 -1.894856 -1.634303 -1.663443 -1.768453
    ## [6,]    1 -1.894856 -1.634303 -1.663443 -1.768453 -1.818069

    # Last 6 rows and columns
    tail(esn$data_train$design_matrix)[,1:6]
    ##       [,1]     [,2]     [,3]     [,4]     [,5]     [,6]
    ## [85,]    1 2.246864 1.793144 1.403657 1.195938 2.067344
    ## [86,]    1 1.793144 1.403657 1.195938 2.067344 4.656129
    ## [87,]    1 1.403657 1.195938 2.067344 4.656129 6.776825
    ## [88,]    1 1.195938 2.067344 4.656129 6.776825 5.562883
    ## [89,]    1 2.067344 4.656129 6.776825 5.562883 4.227307
    ## [90,]    1 4.656129 6.776825 5.562883 4.227307 2.586920

    # Note: The times associated with all values in the design matrix are included in the output

    # Times in first 6 rows and columns
    head(esn$data_train$design_matrix_times)[,1:6]
    ##      [,1]        [,2]         [,3]         [,4]         [,5]        
    ## [1,] "Intercept" "1986-01-01" "1986-02-01" "1986-03-01" "1986-04-01"
    ## [2,] "Intercept" "1986-02-01" "1986-03-01" "1986-04-01" "1986-05-01"
    ## [3,] "Intercept" "1986-03-01" "1986-04-01" "1986-05-01" "1986-06-01"
    ## [4,] "Intercept" "1986-04-01" "1986-05-01" "1986-06-01" "1986-07-01"
    ## [5,] "Intercept" "1986-05-01" "1986-06-01" "1986-07-01" "1986-08-01"
    ## [6,] "Intercept" "1986-06-01" "1986-07-01" "1986-08-01" "1986-09-01"
    ##      [,6]        
    ## [1,] "1986-05-01"
    ## [2,] "1986-06-01"
    ## [3,] "1986-07-01"
    ## [4,] "1986-08-01"
    ## [5,] "1986-09-01"
    ## [6,] "1986-10-01"

    # Times in last 6 rows and columns
    tail(esn$data_train$design_matrix_times)[,1:6]
    ##       [,1]        [,2]         [,3]         [,4]         [,5]        
    ## [85,] "Intercept" "1993-01-01" "1993-02-01" "1993-03-01" "1993-04-01"
    ## [86,] "Intercept" "1993-02-01" "1993-03-01" "1993-04-01" "1993-05-01"
    ## [87,] "Intercept" "1993-03-01" "1993-04-01" "1993-05-01" "1993-06-01"
    ## [88,] "Intercept" "1993-04-01" "1993-05-01" "1993-06-01" "1993-07-01"
    ## [89,] "Intercept" "1993-05-01" "1993-06-01" "1993-07-01" "1993-08-01"
    ## [90,] "Intercept" "1993-06-01" "1993-07-01" "1993-08-01" "1993-09-01"
    ##       [,6]        
    ## [85,] "1993-05-01"
    ## [86,] "1993-06-01"
    ## [87,] "1993-07-01"
    ## [88,] "1993-08-01"
    ## [89,] "1993-09-01"
    ## [90,] "1993-10-01"
    ```

4.  Finally, the $h$ matrix is created, which is used to fit the ridge
    regression. The number of columns corresponds to the number of rows
    in the design matrix. The number of rows is dependent on the number
    of hidden units specified by the user (`nh`) and whether the user
    requested a quadratic term to be included (`add_quad`). In this
    case, the $h$ matrix has 100 columns since `nh = 50` and
    `add_quad = TRUE` (50x2). If `add_quad = FALSE`, then $h$ would only
    have 50 columns. When `add_quad = TRUE`, the first `nh` rows are
    squared to create the second `nh` rows.

    ``` r
    # Dimensions
    dim(esn$h$h)
    ## [1] 100  90

    # First 6 rows and columns
    head(esn$h$h)[,1:6]
    ##              [,1]        [,2]         [,3]        [,4]        [,5]        [,6]
    ## [1,]  0.133280903  0.18066854  0.175721162  0.06982338 -0.01100350 -0.08876392
    ## [2,] -0.076400413 -0.07412147 -0.078557967 -0.10868062 -0.12843916 -0.13373244
    ## [3,] -0.104669051 -0.15203790 -0.092392503 -0.08773511 -0.13661907 -0.22923288
    ## [4,] -0.079521519 -0.02449933  0.001073625  0.04479090  0.06859767  0.09087651
    ## [5,] -0.094811081 -0.11723325 -0.095202810 -0.10298391 -0.09835920 -0.10990292
    ## [6,] -0.001034056  0.03535584  0.041506622  0.05515962  0.04629742  0.05533932

    # Last 6 rows and columns
    tail(esn$h$h)[,1:6]
    ##                [,1]         [,2]         [,3]         [,4]        [,5]
    ##  [95,] 3.961213e-03 0.0026136274 0.0040397847 0.0116685279 0.018794131
    ##  [96,] 4.583901e-02 0.0141732440 0.0086814409 0.0057350103 0.002376154
    ##  [97,] 6.307066e-02 0.0625247583 0.0511151577 0.0557471874 0.069238994
    ##  [98,] 2.346862e-06 0.0008018631 0.0006530502 0.0050765032 0.003700211
    ##  [99,] 1.699889e-02 0.0131812142 0.0008951453 0.0001158859 0.001069244
    ## [100,] 2.675877e-03 0.0076341740 0.0001911897 0.0033208986 0.016364641
    ##               [,6]
    ##  [95,] 0.012315880
    ##  [96,] 0.001469464
    ##  [97,] 0.075245538
    ##  [98,] 0.003430946
    ##  [99,] 0.002400701
    ## [100,] 0.024223858

    # Rows 51 to 100 are the squared values of rows 1 to 50
    identical((esn$h$h[1:50,])^2, esn$h$h[51:100,])
    ## [1] TRUE
    ```
