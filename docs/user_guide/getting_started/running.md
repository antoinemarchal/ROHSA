title: Running ROHSA

#Running ROHSA

---
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

```bash
(marchalenv) Antoines-MacBook-Pro:src antoinemarchal$ ./ROHSA parameters.txt 
 -------------------------------------------------------------------------
16 April 2018  11:23:06.413 AM
 
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
 
 filename = ''
 fileout = ''

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

STOP opening file error
```