## Maintainer comments     
Response to Uwe Ligges:       
* We have fixed the **keyword** and **concept** entries errors which appeared in the last submission.       
* We have cross checked the doi:10.5705/ss.202014.0117 in DESCRIPTION file. And it is right. For the DOI, doi.org redierects to
http://www3.stat.sinica.edu.tw/statistica/J28N1/J28N114/J27N114.html. However, the website can not be visited directly. Alternative, we have to first visit http://www3.stat.sinica.edu.tw/statistica/J28N1/28-1.html, then jump to 
http://www3.stat.sinica.edu.tw/statistica/J28N1/J28N114/J28N114.html.

The following NOTE occurs when checking:      
* *(possibly) invalid doi 10.5705/ss.202014.0117.*    
Actually, this doi is definitely right. R users can search it through google.com or search.crossref.org.    

## Test environments
* Local PC, Windows, R 3.6.0    
* Oracle Solaris 10, x86, 32 bit, R-patched
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit
* Windows Server 2008 R2 SP1, R-release, 32/64 bit

## R CMD check results
* **Local PC, Windows, R 3.6.0** and **Oracle Solaris 10, x86, 32 bit, R-patched**      
Status: OK

* **Windows Server 2008 R2 SP1, R-devel, 32/64 bit** and **Windows Server 2008 R2 SP1, R-release, 32/64 bit**       
NOTE:    
  
```
Maintainer: 'Jin Zhu <zhuj37@mail2.sysu.edu.cn>'

Found the following (possibly) invalid DOIs:
  DOI: 10.5705/ss.202014.0117
From: DESCRIPTION
Status: Not Found
Message: 404
```

## Reverse dependencies
There are no reverse dependencies.
