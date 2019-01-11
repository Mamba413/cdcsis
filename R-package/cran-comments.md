## Maintainer comments     
The following NOTEs result from:    
* the change of Maintainer    
* (possibly) invalid doi 10.5705/ss.202014.0117. (Actually, this doi is definitely right. You can search it through google.com or search.crossref.org.)      
The ERROR is arose from the package dependencies of ks package.

## Test environments
* Local PC, Windows, R 3.5.1    
* Debian Linux, R-devel, GCC ASAN/UBSAN
* Windows Server 2008 R2 SP1, R-devel
* Ubuntu Linux 16.04 LTS, R-release, GCC
* Fedora Linux, R-devel, clang, gfortran

## R CMD check results
* Local PC, Windows, R 3.5.1 and Debian Linux, R-devel, GCC ASAN/UBSAN      
Status: OK

* Ubuntu Linux 16.04 LTS, R-release, GCC and Windows Server 2008 R2 SP1, R-devel      
NOTE:    

```
  Maintainer: 'Jin Zhu <zhuj37@mail2.sysu.edu.cn>'
  
  New maintainer:
    Jin Zhu <zhuj37@mail2.sysu.edu.cn>
  Old maintainer(s):
    Canhong Wen <wencanhong@gmail.com>
  
  Found the following (possibly) invalid DOIs:
    DOI: 10.5705/ss.202014.0117
      From: DESCRIPTION
      Status: Not Found
      Message: 404
```

* Fedora Linux, R-devel, clang, gfortran      
ERROR:      
```
** byte-compile and prepare package for lazy loading
Error in loadNamespace(j <- i[[1L]], c(lib.loc, .libPaths()), versionCheck = vI[[j]]) :
  there is no package called ‘Matrix’
```

NOTE:
```
  Maintainer: 'Jin Zhu <zhuj37@mail2.sysu.edu.cn>'
  
  New maintainer:
    Jin Zhu <zhuj37@mail2.sysu.edu.cn>
  Old maintainer(s):
    Canhong Wen <wencanhong@gmail.com>
  
  Found the following (possibly) invalid DOIs:
    DOI: 10.5705/ss.202014.0117
      From: DESCRIPTION
      Status: Not Found
      Message: 404
```

## Reverse dependencies
There are no reverse dependencies.
