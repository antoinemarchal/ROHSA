title: Installation

#Installing ROHSA

![ROHSA](|media|/LogoMakr_0dTJ9B.png)
{: style="text-align: center" }

---
##Dependencies
ROHSA has the following dependencies:

* [Make](https://www.gnu.org/software/make/) for building ROHSA. (Required)
* [Fortran](https://gcc.gnu.org/) (Required)
* [git](https://git-scm.com/) for version control (Optional)
* [Ford](https://github.com/cmacmackin/ford) for creating HTML documentation (Optional)

##Download the code source
The last version of ROHSA ({!../version!}) can be downloaded from [here](https://github.com/antoinemarchal/ROHSA/releases/)
Once you have downloaded the `.tar.gz` or `.zip` file unpack it using the following commands.
```bash
tar -zxf ROHSA-{!../VERSION!}.tar.gz
```
or if you downloaded the `.zip` file
```bash
unzip ROHSA-{!../VERSION!}.zip
```

##Get the git repository
You can also use git if you want to follow the developement of the code
```bash
git clone https://github.com/antoinemarchal/ROHSA.git ROHSA
```

##Building ROHSA
Once you are in the src/ directory you are now able to build the code running the following
```bash
cd ROHSA/src/
```
```bash
make
```

You can now check if ROHSA is bluit correctly
```bash
(marchalenv) Antoines-MacBook-Pro:src antoinemarchal$ ./ROHSA 
 -------------------------------------------------------------------------
16 April 2018  11:08:50.786 AM
 
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
 
STOP opening file error
```

