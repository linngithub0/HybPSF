# HybPSF
This is PSF reconstruction code for JWST NIRCam image

This is an introduction for compiling and preparing for this code
first you need zunzip the fftw-3.3.10.zip and cfitsio-4.1.10.zip respectively
and then you can just run: ./setup.sh, under the path in the command line
and then see the example.py for the simple usage of the code

the configure paramaters used to generate the catalogue are needed, which
named "config.py", and these variables should be given by user.

After this code are prepared, and the configure are ready, the the user can run:python gen_cat.py
if you want to run gen_cat.py, you should configure the default.sex file with some parameters:
PARAMETERS_NAME FILTER_NAME STARNNW_NAME respectively

to generate the catalogue files, after the catalogue files are generated,then run: python example.py


one more thing, on MacOS, you may encounter the error report such as :**** image not found, 
run:install_name_tool -add_rpath <jwst_psf_path>/lib <the application> , and then run the
python scripy may pass. 
for example: on my Mac, it reports: 
"dyld: Library not loaded: @rpath/libcfitsio.9.dylib
  Referenced from: /Users/linn/Documents/code/JWST/fits/jwst_psf/source/staridf/cut
  Reason: image not found".
  
And then run: "install_name_tool -add_rpath /Users/linn/Documents/code/JWST/fits/jwst_psf/lib /Users/linn/Documents/code/JWST/fits/jwst_psf/cut " in terminal
would fix this problem
May this could help you!


Attentions: sextractor need to be installed on you platform, and the commmand of run the sextractor mya be different in different platform, you should change the 
code in gen_cat.py of line 63 and 69

any problem are welcome, you can throw the proplem in the Wechat group of "JWST revised"
or send email to linn@shao.ac.cn
