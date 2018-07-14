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
   filename = ''         !! name of the data file
   ,fileout = ''         !! name of the output result
   ,filename_noise = ''  !! name of the file with STD map (if noise .eq. true)
   ,n_gauss = 8          !! number of gaussian to fit
   ,lambda_amp = 1.d0    !! lambda for amplitude parameter
   ,lambda_mu = 1.d0     !! lamnda for mean position parameter
   ,lambda_sig = 1.d0    !! lambda for dispersion parameter
   ,maxiter_init = 15000 !! max iteration for L-BFGS-B alogorithm (init mean spectrum)
   ,maxiter = 800        !! max iteration for L-BFGS-B alogorithm
   ,m = 10               !! number of corrections used in the limited memory matrix by LBFGS-B
   ,noise = .false.      !! if false --> STD map computed by ROHSA with lstd and ustd (if true given by the user)
   ,regul = .true.       !! if true --> activate regulation
   ,lstd = 160           !! lower bound to compute the standard deviation map of the cube (if noise .eq. false)
   ,ustd = 198           !! upper bound to compute the standrad deviation map of the cube (if noise .eq. false)
   ,iprint = -1          !! print option
   ,iprint_init = -1     !! print option init
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
(Clumpix) amarchal@glx-calcul3:~/ROHSA/src$ ./ROHSA parameters.txt 
 -------------------------------------------------------------------------
17 April 2018   9:45:45.824 AM
 
   ____     ___    _   _   ____       _    
  |  _ \   / _ \  | | | | / ___|     / \   
  | |_) | | | | | | |_| | \___ \    / _ \  
  |  _ <  | |_| | |  _  |  ___) |  / ___ \ 
  |_| \_\  \___/  |_| |_| |____/  /_/   \_\ 
 
  Version {!../version!}
  ROHSA is released as open source code
  Check out the documentation: https://antoinemarchal.github.io/ROHSA/
 
 run: ./ROHSA parameters.txt
 -------------------------------------------------------------------------
 
 filename = './data.dat'
 fileout = './result.dat'

 ______Parameters_____

 n_gauss =            8
 lambda_amp =    1.0000000000000000     
 lambda_mu =    1.0000000000000000     
 lambda_sig =    1.0000000000000000     
 maxiter_itit =        15000
 maxiter =          800
 lstd =          160
 ustd =          198
 noise =  F
 regul =  T

 
 opening file and reading data
 dim_v, dim_y, dim_x =          395          94          67
 
 number of los =         6298
 nside =            7

 Reshape cube, new dimensions :
 dim_v, dim_y, dim_x =          395         128         128

                     Start iteration

 Init mean spectrum
 
 Update parameters level            0 >           1
RUNNING THE L-BFGS-B CODE

           * * *

Machine precision = 2.220D-16
 N =           24     M =           10

At X0         3 variables are exactly at the bounds

At iterate    0    f=  1.29989D-02    |proj g|=  4.42658D-05

At iterate    1    f=  1.29989D-02    |proj g|=  6.62203D-05

...

```