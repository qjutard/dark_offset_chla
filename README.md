# DARK : computation of dark offsets for CHLA delayed mode
DARK allows the visualization of the time series of minima of an argo float. It is meant to help the Delayed Mode operator decide on a method for the correction of bias or drift by introducing *DARK_OFFSET*. 

By default DARK plots the minima and their median, as well as the offset that was applied in the automatic computation of *CHLA_ADJUSTED*. But it also allows the user to plot offsets computed by a running median on minima, by [DMMC](https://github.com/qjutard/chl_bbp_ttt), or by a Kalman filter on minima (experimental). All computed offsets can be extracted in text files that can directly be used in [DMMC](https://github.com/qjutard/chl_bbp_ttt).

## Recommended "installation" process :
* Clone the repository where you want it ( `$ git clone https://github.com/qjutard/dark_offset_chla` ) or download the release you want to use
* You will need to have R installed as well as the libraries called in *main.R*
* Get the latest profile index and greylist if necessary
* Adapt pathway definitions in *pathways.R* and *main.R*
* Create an alias in your *.bashrc* or *.bash_aliases* ( `alias DARK="~/path/to/repository/DARK.sh"` )
* Go to the working directory (where you want outputs to be written)
* READ THE HELP ( `$ DARK -h` )
* Updates can then be obtained with `$ git pull` if you cloned the repository
* Feel free to contact me directly for help and to open issues for bugs or desired features