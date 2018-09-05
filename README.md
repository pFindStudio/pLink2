# pLink 2

pLink® is a software dedicated for the analysis of chemically cross-linked proteins or protein complexes using mass spectrometry.

pLink 2 is developed as an upgrade of pLink 1. Compared with pLink 1, pLink 2 provides a graphical user interface, and is ~40 times faster with a newly designed index structure. There are also some improvements in the precision.

![](http://pfind.ict.ac.cn/software/pLink/pLink2.png)
![](http://pfind.ict.ac.cn/software/pLink/pLabel.png)


## Cite us
1. [Identification of cross-linked peptides from complex samples.](http://pfind.ict.ac.cn/paper/2012_yang.pdf)
Bing Yang, Yan-Jie Wu, Ming Zhu, Sheng-Bo Fan, Jin-Zhong Lin, Kun Zhang, Shuang Li, Hao Chi, Yu-Xin Li, Hai-Feng Chen, Shu-Kun Luo, Yue-He Ding, Le-Heng Wang, Zhi-Qi Hao, Li-Yun Xiu, She Chen, Ke-Qiong Ye, Si-Min He and Meng-Qiu Dong.
Nature Methods. July 08, 2012. [[abstract]](https://www.nature.com/articles/nmeth.2099)

2. [Mapping native disulfide bonds at a proteome scale.](http://pfind.ict.ac.cn/paper/2015_fan.pdf)
Shan Lu, Sheng-Bo Fan, Bing Yang, Yu-Xin Li, Jia-Ming Meng, Long Wu, Pin Li, Kun Zhang, Mei-Jun Zhang, Yan Fu, Jin-Cai Luo, Rui-Xiang Sun, Si-Min He, Meng-Qiu Dong.
Nature Methods. Feb. 9, 2015. [[abstract]](https://www.nature.com/articles/nmeth.3283)

## Downloads
Please download and read [user_guide.pdf](http://pfind.ict.ac.cn/software/pLink/pLink2%20User%20Guide.pdf) before download and use pLink 2.

pLink 2 is currently free to use. **[Download pLink 2.3](http://pfind.ict.ac.cn/software/pLink/index.html#Downloads)**.

If you have any questions about it, please contact [pLink@ict.ac.cn.](mailto:pLink@ict.ac.cn)

Online discussion: [https://github.com/pFindStudio/pLink2/issues](https://github.com/pFindStudio/pLink2/issues), see [github.pdf](http://pfind.ict.ac.cn/file/github.pdf) for usage.

## pLink Release Notes

#### Version 2.3.4 - September 5 2018
* Fixed a bug when searching super long peptides. Please note that the min length of peptides cannot be shorter than 4aa and the max length of peptides cannot be longer than 120aa.
* Fixed a bug when the unit of fragment tolerance is Da.
* Fixed a bug when activating software in Windows 7.
* Merged the SS flow and the SS_0 flow in disulfide bonds identification.
* Provided the global FDR estimation for intra-protein and inter-protein cross-links.
* Supported lower case amino acids in fasta file.
* Added a warning when the path of pf2 file is invalid.
* Improved GUI usability.
* Updated expiration time to August 1 2019.
* A new license is required since this version and will be valid until August 1 2019.
#### Version 2.3.3 - May 30 2018
* Fixed a bug when labeling the -NH3/-H2O peaks in pLabel.
* Fixed a bug when drawing the FDR curve.
* Fixed a bug when running on Windows 10 Build 1803.
* Improved the spectra preprocessing algorithm.
* License since version 2.3.0 is still valid for this version.
#### Version 2.3.2 - April 4 2018
* Fixed a bug when the m/z of a peak is negative.
* Fixed a bug when the result is empty.
* Fixed a bug when computing E-value with SS_0 linker.
* Fixed a bug when clicking HELP button in pLabel.
* Corrected the mono mass of SS_0 linker.
* Found solution to pLabel's annotation problem on non-Chinese Edition of Windows.
* License since version 2.3.0 is still valid for this version.
#### Version 2.3.1 - February 9 2018
* Fixed a bug on French Edition of Windows.
* License for version 2.3.0 is still valid for this version.
#### Version 2.3.0 - January 31 2018
* Fixed the decimal point bug on French Edition of Windows.
* Fixed a bug on some Windows 10 Pro.
* Fixed bugs when configuring Enzymes and Quantifications.
* Fixed pLabel annotation problem with MGF extracted by PD.
* Disabled the requirement of administrative privileges.
* Improved software activation methods.
* Improved robustness when searching against *.fasta file containing non-alphabet characters.
* Improved GUI usability.
* Added validity check when configuring meta data.
* Extended the max missed cleavages to 5.
* Extended quantification support to BS3 labeling.
* Changed the default variable modification when linker SS used.
* Removed pBuild folder in results.
* Updated the user guide.
#### Version 2.2.1649 - December 30 2017
* First public beta version.

