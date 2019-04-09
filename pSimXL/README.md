# pSimXL

pSimXL is developed to generate simulated MS/MS spectra, including regular, mono-linked, loop-linked, and cross-linked spectra.


# Usage

#### Compile pSimXL

1.	Download `MinGW` and `MSYS`
2.	Set the paths of `MINGW_HOME` and `MSYS_HOME`
3.	Goto the folder of `pSimXL`, and type `make.bat kernel` in cmd
4.	Once the compile is done, the executable program is located in `installer/pSimXL`


#### Generate simulated spectra

1.	Goto the `installer/pSimXL` folder
2.	Customize the parameter file in `example` folder
3.	Goto the `bin` folder, and type `pSimXL.exe bs3.cfg`, where the `bs3.cfg` is the parameter file used to generate bs3 cross-linked spectra, you change it to `ss.cfg` in `example` folder
4.	Find the simulated dataset in `result_output_path` which was set in the `bs3.cfg`



# License
GNU General Public License v3.0

See [LICENSE](https://github.com/pFindStudio/pLink2/blob/master/pSimXL/LICENSE) to see the full text.
