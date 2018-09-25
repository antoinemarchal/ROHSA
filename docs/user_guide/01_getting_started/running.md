title: Running ROHSA

#Running ROHSA

![ROHSA](|media|/LogoMakr_0dTJ9B.png)
{: style="text-align: center" }

---

## User Parameters

### Input file 

Once the source code is compiled, you can run ROHSA using the `parameters.txt` file. 
```bash
&user_parameters
   filename = ''           !! name of the data file
   ,fileout = ''           !! name of the output result
   ,filename_noise = ''    !! name of the file with STD map (if noise .eq. true)
   ,n_gauss = 8            !! number of gaussian to fit
   ,lambda_amp = 1.d0      !! lambda for amplitude parameter
   ,lambda_mu = 1.d0       !! lamnda for mean position parameter
   ,lambda_sig = 1.d0      !! lambda for dispersion parameter
   ,lambda_var_amp = 0.d0  !! lambda for dispersion amplitude parameter
   ,lambda_var_mu = 0.d0   !! lambda for dispersion mean position parameter
   ,lambda_var_sig = 1.d0  !! lambda for dispersion variance parameter
   ,amp_fact_init = 0.66d0 !! times max amplitude of additional Gaussian
   ,sig_init = 5.d0        !! dispersion of additional Gaussian
   ,init_option = "mean"   !! init ROHSA with the mean or the std spectrum
   ,maxiter_init = 15000   !! max iteration for L-BFGS-B alogorithm (init mean spectrum)
   ,maxiter = 400          !! max iteration for L-BFGS-B alogorithm
   ,m = 10                 !! number of corrections used in the limited memory matrix by LBFGS-B
   ,noise = .false.        !! if false --> STD map computed by ROHSA with lstd and ustd (if true given by the user)
   ,regul = .true.         !! if true --> activate regulation
   ,descent = .true.       !! if true --> activate hierarchical descent to initiate the optimization
   ,lstd = 1               !! lower bound to compute the standard deviation map of the cube (if noise .eq. false)
   ,ustd = 20              !! upper bound to compute the standrad deviation map of the cube (if noise .eq. false)
   ,iprint = 1             !! print option
   ,iprint_init = -1       !! print option init
   /
```
`parameters.txt` contains all the free parameters of ROHSA. 

@Warning If one of the parameters if missing, ROHSA wil take the default value encoded in the `main.f90` source file. 

### Print option (L-BFGS-B argument)

iprint is an INTEGER variable that must be set by the user. It controls the frequency and type of output generated:

```bash
iprint<0    no output is generated;
iprint=0    print only one line at the last iteration;
0<iprint<99 print also f and |proj g| every iprint iterations;
iprint=99   print details of every iteration except n-vectors;
iprint=100  print also the changes of active set and final x;
iprint>100  print details of every iteration including x and g;
When iprint > 0, the file iterate.dat will be created to summarize the iteration.
```

@Note For any question about the parameters, please contact us at the following address : antoine.marchal@ias.u-psud.fr

### Data format
ROHSA does not support FITS file in version {!../version!}. It wiil be added in a next release. 
You can find a python code : `fits2dat.py` in `src/` directory to convert your `.fits` file into a `.dat` file. 

## Running ROHSA
You can now run your Gaussian decomposition with ROHSA. Enjoy !

```bash
./ROHSA parameters.txt 
```  

```bash
(marchalenv) dapmcw133:src antoinemarchal$ ./ROHSA parameters.txt 
 
 opening file and reading data
 -------------------------------------------------------------------------
25 September 2018   5:09:11.592 PM
 
   ____     ___    _   _   ____       _    
  |  _ \   / _ \  | | | | / ___|     / \   
  | |_) | | | | | | |_| | \___ \    / _ \  
  |  _ <  | |_| | |  _  |  ___) |  / ___ \ 
  |_| \_\  \___/  |_| |_| |____/  /_/   \_\ 
 
  Version 1.0.0
  ROHSA is released as open source code
  Check out the documentation: https://antoinemarchal.github.io/ROHSA/
 
 run: ./ROHSA parameters.txt
 -------------------------------------------------------------------------
 fileout = 'result.dat'

 ______Parameters_____

 n_gauss =            6
 lambda_amp =    1.0000000000000000     
 lambda_mu =    1.0000000000000000     
 lambda_sig =    1.0000000000000000     
 lambda_var_amp =    0.0000000000000000     
 lambda_var_mu =    0.0000000000000000     
 lambda_var_sig =    1.0000000000000000     
 amp_fact_init =   0.66000000000000003     
 sig_init =    5.0000000000000000     
 init_option = mean    
 maxiter_init =        15000
 maxiter =           20
 lstd =            1
 ustd =           20
 noise =  F
 regul =  T
 descent =  T

 dim_v, dim_y, dim_x =          101          16          16
 
 number of los =          256
 nside =            4

 Reshape cube, new dimensions :
 dim_v, dim_y, dim_x =          101          16          16

 Compute mean and std spectrum
                     Start iteration

 Start hierarchical descent
 Init mean spectrum
 
 Update parameters level            0 >           1
 
 Update parameters level            1 >           2
 
 Update parameters level            2 >           4
 
 Update parameters level            3 >           8
 
 Update parameters level            4 >          16

 Reshape cube, restore initial dimensions :
 dim_v, dim_y, dim_x =          101          16          16

 Update last level ...


 _____ Write output file _____

 ##################################################################
25 September 2018   5:09:11.654 PM
  Terminate
                         ROHSA ALGORITHM
 
 ##################################################################
```