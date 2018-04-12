---
project: ROHSA
Date: April 8, 2018
version: {!../version!}
output_dir: ./Doc
src_dir: /Users/antoinemarchal/Desktop/ROHSA/src/
page_dir: /Users/antoinemarchal/Desktop/ROHSA/docs/user_guide/
project_github: https://github.com/antoinemarchal/
project_download: https://github.com/antoinemarchal/ROHSA/releases/
summary: ![ROHSA](|media|/LogoMakr_0dTJ9B.png)
         {: style="text-align: center" }
author: Antoine Marchal
author_description: PhD Student in Astrophysics at IAS/CEA
github: https://github.com/antoinemarchal/
email: antoine.marchal@ias.u-psud.fr
fpp_extensions: fpp
predocmark: >
media_dir: ./media
docmark_alt: #
predocmark_alt: <
display: public
         protected
         private
source: true
graph: true
search: true
macro: TEST
       LOGIC=.true.
extra_mods: json_module: http://jacobwilliams.github.io/json-fortran/
            futility: http://cmacmackin.github.io
license: by-nc
extra_filetypes: sh #
---

ROHSA (Regularized Optimization for Hyper-Spectral Analysis) was developped in Paris-Saclay University
(IAS/CEA) to study the statistical properties of interstellar gas through atomic and molecular lines.  
This code is a "Gaussian Decomposition Algorithm" designed to decompose any kind of hyper-spectral 
observations into a sum of coherent Gaussian. It is written in Fortran 90 on CPU. An extended version 
with GPU implementation is currently in development. 

![ROHSA](|media|/screenshot.png)
{: style="text-align: center" }

Exemple of a three componant Gaussian decomposition with ROHSA.
{: style="text-align: center" }

This work is part of the [Hyperstars](http://hyperstars.lmpa.eu) project supported by 
[MASTODONS](http://www.cnrs.fr/mi/spip.php?article53&lang=fr).







