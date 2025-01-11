# pLink 2

## Notification
**pLink® 3 is now available. [Github Release: pLink3](https://github.com/pFindStudio/pLink3/releases).**

* pLink3 has replaced pLink2.

* If you need to extend the validity of a license for a previous version, please contact us at support@pfind.org. However, we recommend upgrading to pLink3, as maintenance for pLink2 will be discontinued in the future.

* Technical support or bug reporting: https://github.com/pFindStudio/pLink3/issues
	* We encourage to use the GitHub issues platform for communication, as it can help more people efficiently. 
	* Of course, you are also welcome to reach out via email.

* License self-service: http://i.pfind.org (Alternative: http://i.pfind.net )

* Contact: support@pfind.org. (Alternative: plink@ict.ac.cn )

* Website: http://pfind.org (Alternative: http://pfind.net )

---
pLink® is a software dedicated for the analysis of chemically cross-linked proteins or protein complexes using mass spectrometry.

pLink 2 is developed as an upgrade of pLink 1. Compared with pLink 1, pLink 2 provides a graphical user interface, and is ~40 times faster with a newly designed index structure. There are also some improvements in the precision.

Our paper [A high-speed search engine pLink 2 with systematic evaluation for proteome-scale identification of cross-linked peptides](https://www.nature.com/articles/s41467-019-11337-z) has been accepted for publication by *Nature Communications*. Congratulations and please [cite this paper](https://github.com/pFindStudio/pLink2#cite-us) if you used pLink 2.

![](http://pfind.org/software/pLink/pLink_2_work%EF%AC%82ow.png)
![](http://pfind.org/software/pLink/pLink2.png)
![](http://pfind.org/software/pLink/pLabel.png)


## Cite us
[A high-speed search engine pLink 2 with systematic evaluation for proteome-scale identification of cross-linked peptides.](https://www.nature.com/articles/s41467-019-11337-z.pdf)
Zhen-Lin Chen, Jia-Ming Meng, Yong Cao, Ji-Li Yin, Run-Qian Fang, Sheng-Bo Fan, Chao Liu, Wen-Feng Zeng, Yue-He Ding, Dan Tan, Long Wu, Wen-Jing Zhou, Hao Chi, Rui-Xiang Sun, Meng-Qiu Dong & Si-Min He.
*Nature Communications.* July 30, 2019. [[abstract]](https://www.nature.com/articles/s41467-019-11337-z)

## Downloads
Please download and read [user_guide.pdf](http://pfind.org/software/pLink/pLink2%20User%20Guide.pdf) before download and use pLink 2.

pLink 2 is currently free to use. **[Download pLink 2.3](https://github.com/pFindStudio/pLink2/raw/master/installer/pLink2.3.11.exe)**.

If you have any questions about it, please contact [pLink@ict.ac.cn.](mailto:pLink@ict.ac.cn)

Online discussion: [https://github.com/pFindStudio/pLink2/issues](https://github.com/pFindStudio/pLink2/issues), see [github.pdf](http://pfind.org/file/github.pdf) for usage.

## Notice
* The expiration date of **pLink 2.3.11** is set on **Jan. 10, 2025**.
* If RAW file read required, [MSFileReader](https://github.com/pFindStudio/pLink2/wiki/FAQ#how-to-install-msfilereader), 3.0 SP2 or a **lower** version, both 32 bit and 64 bit, need be installed first.
    * *If you have problem to run pLink 2 on Windows 7, please see this [issue](https://github.com/pFindStudio/pLink2/issues/68). Re-install MSFileReader using "Run as administrator".*
* [.NET framework 4.5](https://dotnet.microsoft.com/en-us/download/dotnet-framework/thank-you/net452-offline-installer) or a higher version is required.
* [Java](https://www.java.com/en/download/) is also required for quantitation.

## pLink Release Notes

#### Version 2.3.11 - March 28 2022
* Change the website to [http://pfind.org](http://pfind.org) from http://pfind.ict.ac.cn.
* Add the self-service registration website [http://i.pfind.org](http://i.pfind.org).
* The license for version 2.3.10 is still valid for version 2.3.11. The same expiration date is Jan. 10, 2025.
* If you have problem to run pLink 2 on Windows 7, please see this [issue](https://github.com/pFindStudio/pLink2/issues/68). Re-install MSFileReader using "Run as administrator".

#### Version 2.3.10 - January 12 2022
* Fixed a major bug in calculating E-value that causes software to crash, thank you Olexandr Dybkov and Liu Lab who reported.
* Changed the format of 'Precursor_Mass' and 'Peptide_Mass' in the report to decimal, thank you [@IwanParf](https://github.com/pFindStudio/pLink2/issues/70).
* Fixed several other bugs.
* All users are recommended to update to this version if required.
* A new license is required for this version, but we have extended the license validity to three years. This means you can use pLink 2.3.10 until 2025.1.10 without worrying about expiration.
* If you have problem to run pLink 2 on Windows 7, please see this [issue](https://github.com/pFindStudio/pLink2/issues/68). Re-install MSFileReader using "Run as administrator".
#### Version 2.3.9 - January 7 2020
* Fixed a pParse bug, thank you lili.
* Fixed a fasta reading bug, thank you [@vladtheimpalerSr](https://github.com/pFindStudio/pLink2/issues/66).
* Updated the user guide, please read it before installation.
* All users using version 2.3.8 are recommended to update to this version.
* License for version 2.3.8 is still valid for this version.
* If you have problem to run pLink 2.3.9 on Windows 7, please see this [issue](https://github.com/pFindStudio/pLink2/issues/68).
#### Version 2.3.8 - December 27 2019
* Happy New Year 2020!
* Fixed a quantitation bug on French Edition of Windows.
* Updated the configuration for cross-linker.
* Updated pLabel for MS-cleavable cross-linked spectra annotation.
* Updated the user guide, please read it before installation.
* A new license is required for this version, but we have extended the license validity to three years. This means you can use pLink 2.3.8 until 2023 without worrying about expiration.
#### Version 2.3.7 - November 5 2019
* Fixed a bug when open pLabel.
* License for version 2.3.5 is still valid for this version.
#### Version 2.3.6 - September 30 2019
* Congratulations! Our [paper of pLink 2](https://www.nature.com/articles/s41467-019-11337-z) is finally published, and please cite this paper if you used pLink 2.
* Fixed a bug when register the software.
* Fixed a bug when showing the symbol of β in pLabel.
* Disabled the SUMO identification flow temporarily, as we are optimizing and testing it. This feature will be available as soon as possible.
* Updated the lists of modifications and enzymes.
* Updated the GUI.
* License for version 2.3.5 is still valid for this version.
#### Version 2.3.5 - December 29 2018
* [Happy New Year 2019, and please see greetings from your pLink developers!](http://pfind.org/news.html#pLink_Greetings_2019)
* [Fixed a bug when infer proteins from peptides](https://github.com/pFindStudio/pLink2/issues/39), thank you @zrpeak.
* Fixed a bug when the database is small.
* Fixed a bug when save a disulfide bond search task.
* Fixed a bug in modification tab of pConfig.
* Add a warning when combineSS.exe failed.
* Updated expiration time to January 1 2020.
* A new license is required since this version and will be valid until January 1 2020.
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

