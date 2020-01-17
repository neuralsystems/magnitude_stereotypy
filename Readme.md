# Description 

The repository contains all the source code and data files to reproduce the results

## Getting Started

Download the source files, extract to a new folder, and add this folder and all its subfolders to your MATLAB path

### Prerequisites

* MATLAB (this code has been tested in version r2016a)

## Running the code

* Please run the ```main.m``` file present in the folder. It will run all simulations, run analyses on these simulations and create the figure files (NOTE that running the entire code might need atleast 32 gb RAM)

## Dataset details

* The folder named ```data``` contains all the dataset files used in the code. The ```preparealldatasets``` function generates the ```.mat``` files used for plotting the datasets. The details of each dataset is as follows:
** Fly GCaMP3 data: The stereotypy values for each lobe for each odor are specified in the file ```fly_mbon_gcamp3_stereotypy.xlsx```. 'wt' and 'APL>TNT' data is marked accordingly. The raw data images are available on request from the authors.
** Locust bLN1 data: The file ```bln1_stereotypy.mat``` contains the spike rates extracted from the raster plots in Supplementary Fig. 1a in the variable ```meanfiringNorm``` in the format num_individual x num_odors. Raw traces are available on request from the authors.
** Murthy et al 2008 data: The file ```murthy_2008.xlsx``` contains the KC response data that was manually extracted from Fig. 3A of the paper. '0', '1' and '-' denote 'response', 'no response' and 'odor not tested' respectively.
** Shimizu and Stopfer 2017 data: The file ```shimizu_2017.xlsx``` describes the data files (present in the ```raw``` folder) containing spike time data for each PN and odor. Raw traces are available on request from the authors.
** Schaffer et al. 2018: The raw data folder contains the results of the 6 different runs of the simulations as described in the paper. The files ```mc_pop4thOrderRespSmall.mat``` and ```wmc_pop4thOrderRespSmall.mat``` contain the results of running the script ```calculate4thOrderPopResp_v2small.m``` to calculate the SNR for the different odors with and without weight normalization, respectively.

## Built With

* MATLAB (version r2016a)
* modified version of the **gramm** plotting package available at [piermorel/gramm](https://github.com/piermorel/gramm). I modified the ```stat_violin``` function to suit our plotting needs.
* **export_fig** toolbox available at [altmany/export_fig](https://github.com/altmany/export_fig)

## Authors

* **Aarush Mohit Mittal** - *All files except gramm and export_fig packages*

## License

This project is licensed under the MIT License - see the [LICENSE.md](License.md) file for details
