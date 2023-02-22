
# model.energy/efuels: online optimisation of hydrogen-based e-fuel delivery

This is the code for the online optimisation of hydrogen-based e-fuel
delivery from wind and solar.

It uses only free software and open data, including [Python for Power
System Analysis (PyPSA)](https://github.com/PyPSA/PyPSA) for the
optimisation framework, the European Centre for Medium-Range Weather
Forecasts (ECMWF) [ERA5
dataset](https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels)
for the open weather data, the [atlite
library](https://github.com/FRESNA/atlite) for converting weather data
to generation profiles, [Clp](https://projects.coin-or.org/Clp) for
the solver, [D3.js](https://d3js.org/) for graphics,
[Mapbox](https://www.mapbox.com/), [Leaflet](http://leafletjs.com/)
and [Natural Earth](https://www.naturalearthdata.com/) for maps, and
free software for the server infrastructure (GNU/Linux, nginx, Flask,
gunicorn, Redis).

You can find a live version at:

<https://model.energy/efuels/>

The tool is a simplified version of the model behind the research
paper [Import options for chemical energy carriers from renewable
sources to Germany](https://doi.org/10.1371/journal.pone.0281380) by
Johannes Hampp, Michael Düren and Tom Brown in PLOS ONE, 2023.



## Requirements

### Software

This software has only been tested on the Ubuntu distribution of GNU/Linux.

Ubuntu packages:

`sudo apt install coinor-clp coinor-cbc redis-server`

To install, we recommend using [miniconda](https://docs.conda.io/en/latest/miniconda.html) in combination with [mamba](https://github.com/QuantStack/mamba).

	conda install -c conda-forge mamba
	mamba env create -f environment.yaml

For (optional) server deployment:

	sudo apt install nginx
	mamba install gunicorn


### Data

This uses the same data as <https://model.energy/>, see the [WHOBS-server repository](https://github.com/PyPSA/whobs-server).


### Regenerate default assumptions

The basic assumptions in `defaults-initial.csv` are added to by the script `generate_defaults.py`, which adds assumptions from the [trace repository](https://github.com/euronion/trace) and [technology-data repository](https://github.com/PyPSA/technology-data) (which is in turn largely based on the [Danish Energy Agency's Technology Data](https://ens.dk/en/our-services/projections-and-models/technology-data)). The CSV `defaults.csv` is an output of this script and should not be modified. All modifications should be make in `defaults-initial.csv`.

To run the script do:

`python generate_defaults.py`


## Run server locally on your own computer

To run locally you need to start the Python Flask server in one terminal, and redis in another:

Start the Flask server in one terminal with:

`python server.py`

This will serve to local address:

http://127.0.0.1:5002/

In the second terminal start Redis:

`rq worker efuels`

where `efuels` is the name of the queue. No jobs will be solved until
this is run. You can run multiple workers to process jobs in parallel.


## Deploy on a publicly-accessible server

Use nginx, gunicorn for the Python server, rq, and manage with supervisor.


## License

Copyright 2022-3 Tom Brown <https://nworbmot.org/>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation; either [version 3 of the
License](LICENSE.txt), or (at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the [GNU
Affero General Public License](LICENSE.txt) for more details.
