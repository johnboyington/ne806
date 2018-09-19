          Code Name & Version = MCNP6, 1.0
  
     _/      _/        _/_/_/       _/      _/       _/_/_/         _/_/_/
    _/_/  _/_/      _/             _/_/    _/       _/    _/     _/       
   _/  _/  _/      _/             _/  _/  _/       _/_/_/       _/_/_/    
  _/      _/      _/             _/    _/_/       _/           _/    _/   
 _/      _/        _/_/_/       _/      _/       _/             _/_/      
  
  +---------------------------------------------------------------------+
  | Copyright 2008. Los Alamos National Security, LLC.  All rights      |
  | reserved.                                                           |
  |   This material was produced under U.S. Government contract         |
  | DE-AC52-06NA25396 for Los Alamos National Laboratory, which is      |
  | operated by Los Alamos National Security, LLC, for the U.S.         |
  | Department of Energy. The Government is granted for itself and      |
  | others acting on its behalf a paid-up, nonexclusive, irrevocable    |
  | worldwide license in this material to reproduce, prepare derivative |
  | works, and perform publicly and display publicly. Beginning five    |
  | (5) years after 2008, subject to additional five-year worldwide     |
  | renewals, the Government is granted for itself and others acting on |
  | its behalf a paid-up, nonexclusive, irrevocable worldwide license   |
  | in this material to reproduce, prepare derivative works, distribute |
  | copies to the public, perform publicly and display publicly, and to |
  | permit others to do so. NEITHER THE UNITED STATES NOR THE UNITED    |
  | STATES DEPARTMENT OF ENERGY, NOR LOS ALAMOS NATIONAL SECURITY, LLC, |
  | NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, |
  | OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY,  |
  | COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT, |
  | OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE |
  | PRIVATELY OWNED RIGHTS.                                             |
  +---------------------------------------------------------------------+
  
1mcnp     version 6     ld=05/08/13                     09/19/18 10:58:59 
 *************************************************************************                 probid =  09/19/18 10:58:59 
 r name=solid.i tasks 26                                                         

 
  warning.  Physics models disabled.
         1-       Bonner Sphere Template                                                          
         2-       c                                                                               
         3-       c Updated  7/12/17 by John Boyington                                            
         4-       c                                                                               
         5-       c ******************************************************************************
         6-       c                               CELL CARDS                                      
         7-       c ******************************************************************************
         8-       c                                                                               
         9-       c                          -----A Sphere*------                                 
        10-       1 0         1 -2 -3      IMP:N=1                                                
        11-       2 0         #1   -3      IMP:N=1                                                
        12-       3 0         3            IMP:N=0                                                
        13-                                                                                       
        14-       c ******************************************************************************
        15-       c                               SURFACE CARDS                                   
        16-       c ******************************************************************************
        17-       1 PX   100.0                                                                    
        18-       2 C/X    0.0  20.0  4.5135166683820502                                          
        19-       3 SO   200.0                                                                    
        20-                                                                                       
        21-       c ******************************************************************************
        22-       c                               DATA CARDS                                      
        23-       c ******************************************************************************
        24-       SDEF   POS=0 0 0                                                                
        25-       c                                                                               
        26-       NPS 1E7                                                                         
        27-       c                                                                               
        28-       c  -----------------------------------------------------------------------------
        29-       c                                                          MATERIAL CARDS       
        30-       c  -----------------------------------------------------------------------------
        31-       c  -----------------------------------------------------------------------------
        32-       c                                                                               
        33-       c                                                                               
        34-       c                                                                               
        35-       c  -----------------------------------------------------------------------------
        36-       c                                                             TALLY CARDS       
        37-       c  -----------------------------------------------------------------------------
        38-       c                                                                               
        39-       c  -----------------------------------------------------------------------------
        40-       c  TALLY 2:         That disk.                                                  
        41-       c  -----------------------------------------------------------------------------
        42-       F2:N 1                                                                          
        43-       FS2  2 T                                                                        
        44-       FM2   12.566370614359172                                                        
        45-       SD2 1 1 1                                                                       
        46-       c                                                                               
        47-       c                                                                               
        48-       c ******************************************************************************
        49-       c                             END OF INPUT FILE                                 
        50-       c ******************************************************************************
        51-                                                                                       
 
  comment.  total nubar used if fissionable isotopes are present.
 
  warning.  no cross-section tables are called for in this problem.
1cells                                                                                                  print table 60

                               atom        gram                                            neutron                                     
              cell      mat   density     density     volume       mass            pieces importance                                   

        1        1        0  0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00           0  1.0000E+00                                   
        2        2        0  0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00           0  1.0000E+00                                   
        3        3        0  0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00           0  0.0000E+00                                   

 total                                               0.00000E+00 0.00000E+00

    minimum source weight = 1.0000E+00    maximum source weight = 1.0000E+00

 ***************************************************
 * Random Number Generator  =                    1 *
 * Random Number Seed       =       19073486328125 *
 * Random Number Multiplier =       19073486328125 *
 * Random Number Adder      =                    0 *
 * Random Number Bits Used  =                   48 *
 * Random Number Stride     =               152917 *
 ***************************************************

 
  comment.  threading will be used when possible in portions of mcnp6.
 
  comment.  threading will be used for n/p/e table physics.
 
  comment.  threading will generally not be used for model physics.

         2 warning messages so far.

 ***********************************************************************************************************************

 dump no.    1 on file solid.ir     nps =           0     coll =              0     ctm =        0.00   nrn =           
      0

         2 warning messages so far.
1problem summary                                                                                                           

      run terminated when    10000000  particle histories were done.
+                                                                                                    09/19/18 11:01:09 
      Bonner Sphere Template                                                               probid =  09/19/18 10:58:59 

 neutron creation    tracks      weight        energy            neutron loss        tracks      weight        energy
                                 (per source particle)                                           (per source particle)

 source            10000000    1.0000E+00    1.4000E+01          escape            10000000    1.0000E+00    1.4000E+01
 nucl. interaction        0    0.            0.                  energy cutoff            0    0.            0.        
 particle decay           0    0.            0.                  time cutoff              0    0.            0.        
 weight window            0    0.            0.                  weight window            0    0.            0.        
 cell importance          0    0.            0.                  cell importance          0    0.            0.        
 weight cutoff            0    0.            0.                  weight cutoff            0    0.            0.        
 e or t importance        0    0.            0.                  e or t importance        0    0.            0.        
 dxtran                   0    0.            0.                  dxtran                   0    0.            0.        
 forced collisions        0    0.            0.                  forced collisions        0    0.            0.        
 exp. transform           0    0.            0.                  exp. transform           0    0.            0.        
 upscattering             0    0.            0.                  downscattering           0    0.            0.        
 photonuclear             0    0.            0.                  capture                  0    0.            0.        
 (n,xn)                   0    0.            0.                  loss to (n,xn)           0    0.            0.        
 prompt fission           0    0.            0.                  loss to fission          0    0.            0.        
 delayed fission          0    0.            0.                  nucl. interaction        0    0.            0.        
 prompt photofis          0    0.            0.                  particle decay           0    0.            0.        
 tabular boundary         0    0.            0.                  tabular boundary         0    0.            0.        
 tabular sampling         0    0.            0.                  elastic scatter          0    0.            0.        
     total         10000000    1.0000E+00    1.4000E+01              total         10000000    1.0000E+00    1.4000E+01

   number of neutrons banked                       0        average time of (shakes)              cutoffs
   neutron tracks per source particle     1.0000E+00          escape            3.9076E+00          tco   1.0000E+33
   neutron collisions per source particle 0.0000E+00          capture           0.0000E+00          eco   0.0000E+00
   total neutron collisions                        0          capture or escape 3.9076E+00          wc1  -5.0000E-01
   net multiplication              1.0000E+00 0.0000          any termination   3.9076E+00          wc2  -2.5000E-01

 computer time so far in this run  1452.67 minutes            maximum number ever in bank         0
 computer time in mcrun            1451.82 minutes            bank overflows to backup file       0
 source particles per minute            6.8879E+03
 random numbers generated                 25467146            most random numbers used was          22 in history      347006

 range of sampled source weights = 1.0000E+00 to 1.0000E+00

 number of histories processed by each thread
      384845      384632      384637      384607      384612      384621      384591      384533      384609      384591
      384625      384624      384563      384632      384607      384581      384647      384641      384610      384561
      384640      384573      384595      384617      384606      384600
1neutron  activity in each cell                                                                         print table 126

                       tracks     population   collisions   collisions     number        flux        average      average
              cell    entering                               * weight     weighted     weighted   track weight   track mfp
                                                          (per history)    energy       energy     (relative)      (cm)

        1        1        8244         8244            0    0.0000E+00   1.4000E+01   1.4000E+01   1.0000E+00   0.0000E+00
        2        2    10006971     10000000            0    0.0000E+00   1.4000E+01   1.4000E+01   1.0000E+00   0.0000E+00

           total      10015215     10008244            0    0.0000E+00

1tally        2        nps =    10000000
           tally type 2    particle flux averaged over a surface.                              
           particle(s): neutrons 

           this tally is all multiplied by  1.25664E+01

           areas   
                surface:       1                                                                                   
               segment
                 1       1.00000E+00
                 2       1.00000E+00
              whole      1.00000E+00
 
 surface  1                                                                                                                            
 segment:             2                                                                                                                
                 0.00000E+00 0.0000
 
 surface  1                                                                                                                            
 segment:            -2                                                                                                                
                 6.09189E-03 0.0145
 
 surface  1                                                                                                                            
 segment:   whole surface                                                                                                              
                 6.09189E-03 0.0145


 ===================================================================================================================================

           results of 10 statistical checks for the estimated answer for the tally fluctuation chart (tfc) bin of tally        2

 tfc bin     --mean--      ---------relative error---------      ----variance of the variance----      --figure of merit--     -pdf-
 behavior    behavior      value   decrease   decrease rate      value   decrease   decrease rate       value     behavior     slope

 desired      random       <0.10      yes      1/sqrt(nps)       <0.10      yes        1/nps           constant    random      >3.00
 observed     random        0.01      yes          yes            0.00      yes         yes            constant    random      10.00
 passed?        yes          yes      yes          yes             yes      yes         yes               yes        yes         yes

 ===================================================================================================================================


 this tally meets the statistical criteria used to form confidence intervals: check the tally fluctuation chart to verify.
 the results in other bins associated with this tally may not meet these statistical criteria.

 ----- estimated confidence intervals:  -----

 estimated asymmetric confidence interval(1,2,3 sigma): 6.0042E-03 to 6.1809E-03; 5.9158E-03 to 6.2692E-03; 5.8274E-03 to 6.3576E-03
 estimated  symmetric confidence interval(1,2,3 sigma): 6.0035E-03 to 6.1802E-03; 5.9152E-03 to 6.2686E-03; 5.8268E-03 to 6.3569E-03

1analysis of the results in the tally fluctuation chart bin (tfc) for tally        2 with nps =    10000000  print table 160


 normed average tally per history  = 6.09189E-03          unnormed average tally per history  = 6.09189E-03
 estimated tally relative error    = 0.0145               estimated variance of the variance  = 0.0002
 relative error from zero tallies  = 0.0145               relative error from nonzero scores  = 0.0001

 number of nonzero history tallies =        4752          efficiency for the nonzero tallies  = 0.0005
 history number of largest  tally  =      568490          largest  unnormalized history tally = 1.29376E+01
 (largest  tally)/(average tally)  = 2.12375E+03          (largest  tally)/(avg nonzero tally)= 1.00921E+00

 (confidence interval shift)/mean  = 0.0001               shifted confidence interval center  = 6.09253E-03


 if the largest  history score sampled so far were to occur on the next history, the tfc bin quantities would change as follows:

      estimated quantities           value at nps           value at nps+1           value(nps+1)/value(nps)-1.

      mean                            6.09189E-03             6.09318E-03                     0.000212
      relative error                  1.45032E-02             1.45016E-02                    -0.000105
      variance of the variance        2.10154E-04             2.10110E-04                    -0.000211
      shifted center                  6.09253E-03             6.09253E-03                     0.000000
      figure of merit                 3.27463E+00             3.27532E+00                     0.000211

 the 100 largest  history tallies appear to have a  maximum value of about 1.29376E+01
 the large score tail of the empirical history score probability density function appears to have no unsampled regions.

 fom = (histories/minute)*(f(x) signal-to-noise ratio)**2 = (6.888E+03)*( 2.180E-02)**2 = (6.888E+03)*(4.754E-04) = 3.275E+00

1status of the statistical checks used to form confidence intervals for the mean for each tally bin


 tally   result of statistical checks for the tfc bin (the first check not passed is listed) and error magnitude check for all bins

        2   passed the 10 statistical checks for the tally fluctuation chart bin result               
         passed all bin error check:     3 tally bins had     1 bins with zeros and     0 bins with relative errors exceeding 0.10


 the 10 statistical checks are only for the tally fluctuation chart bin and do not apply to other tally bins.

 the tally bins with zeros may or may not be correct: compare the source, cutoffs, multipliers, et cetera with the tally bins.

1tally fluctuation charts                              

                            tally        2
          nps      mean     error   vov  slope    fom
       512000   6.0588E-03 0.0643 0.0041  0.0 3.2E+00
      1024000   6.0095E-03 0.0456 0.0021  0.0 3.2E+00
      1536000   6.0348E-03 0.0372 0.0014 10.0 3.3E+00
      2048000   6.0597E-03 0.0321 0.0010 10.0 3.1E+00
      2560000   6.0748E-03 0.0287 0.0008 10.0 3.2E+00
      3072000   6.0051E-03 0.0264 0.0007 10.0 3.1E+00
      3584000   6.0094E-03 0.0244 0.0006 10.0 3.2E+00
      4096000   6.0313E-03 0.0228 0.0005 10.0 3.2E+00
      4608000   5.9231E-03 0.0217 0.0005 10.0 3.1E+00
      5120000   5.9941E-03 0.0204 0.0004 10.0 3.2E+00
      5632000   6.0114E-03 0.0195 0.0004 10.0 3.2E+00
      6144000   5.9883E-03 0.0187 0.0003 10.0 3.2E+00
      6656000   6.0516E-03 0.0178 0.0003 10.0 3.2E+00
      7168000   6.0218E-03 0.0172 0.0003 10.0 3.2E+00
      7680000   6.0877E-03 0.0166 0.0003 10.0 3.2E+00
      8192000   6.0703E-03 0.0161 0.0003 10.0 3.2E+00
      8704000   6.0505E-03 0.0156 0.0002 10.0 3.2E+00
      9216000   6.0732E-03 0.0151 0.0002 10.0 3.3E+00
      9728000   6.0895E-03 0.0147 0.0002 10.0 3.3E+00
     10000000   6.0919E-03 0.0145 0.0002 10.0 3.3E+00

 ***********************************************************************************************************************

 dump no.    2 on file solid.ir     nps =    10000000     coll =              0     ctm =     1451.82   nrn =         
 25467146

         2 warning messages so far.


 run terminated when    10000000  particle histories were done.

 computer time = 1452.67 minutes

 mcnp     version 6     05/08/13                     09/19/18 11:01:09                     probid =  09/19/18 10:58:59 
