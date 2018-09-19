Bonner Sphere Template
c
c Updated  7/12/17 by John Boyington 
c 
c ****************************************************************************** 
c                               CELL CARDS 
c ******************************************************************************
c
c                          -----A Sphere*------
1 0         1 -2 -3      IMP:N=1
2 0         #1   -3      IMP:N=1
3 0         3            IMP:N=0

c ****************************************************************************** 
c                               SURFACE CARDS 
c ****************************************************************************** 
1 PX   100.0
2 C/X    0.0  20.0  4.5135166683820502
3 SO   200.0

c ****************************************************************************** 
c                               DATA CARDS 
c ******************************************************************************
SDEF   POS=0 0 0
c
NPS 1E9
c
c  -----------------------------------------------------------------------------
c                                                          MATERIAL CARDS
c  -----------------------------------------------------------------------------
c  -----------------------------------------------------------------------------
c
c
c
c  -----------------------------------------------------------------------------
c                                                             TALLY CARDS       
c  -----------------------------------------------------------------------------
c
c  -----------------------------------------------------------------------------
c  TALLY 2:         That disk.
c  -----------------------------------------------------------------------------
F2:N 1
FS2  2 T
FM2   12.566370614359172
SD2 1 1 1
c
c
c ****************************************************************************** 
c                             END OF INPUT FILE
c ******************************************************************************
