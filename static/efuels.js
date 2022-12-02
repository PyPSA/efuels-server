// Copyright 2022 Tom Brown

// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU Affero General Public License as
// published by the Free Software Foundation; either version 3 of the
// License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.

// License and more information at:
// https://github.com/PyPSA/efuels-server


// Centered on Frankfurt
var map = L.map('mapid').setView([20,0], 2);

L.tileLayer('https://api.mapbox.com/styles/v1/{id}/tiles/{z}/{x}/{y}?access_token={accessToken}', {
    attribution: 'Map data &copy; <a href="https://www.openstreetmap.org/">OpenStreetMap</a> contributors, <a href="https://creativecommons.org/licenses/by-sa/2.0/">CC-BY-SA</a>, Imagery © <a href="https://www.mapbox.com/">Mapbox</a>',
    tileSize: 512,
    maxZoom: 18,
    zoomOffset: -1,
    id: 'mapbox/streets-v11',
    accessToken: 'pk.eyJ1IjoibndvcmJtb3QiLCJhIjoiY2prbWxibTUyMjZsMDNwcGp2bHR3OWZsaSJ9.MgSprgR6BEbBLXl5rPvXvQ'
}).addTo(map);



assumptions = {"year": 2011,
	       "cf_exponent": 2,
	       "destination" : {"lat": 50.11, "lng" : 8.68},
	       "source" : {"lat": -30.8, "lng" : 18.5},
	    "frequency" : 3,
	    "efuels_load" : 100,
	    "wind_min" : 0,
	    "solar_min" : 0,
	    "wind_max" : 1e7,
	    "solar_max" : 1e7,
	    "wind_discount" : 5,
	    "solar_discount" : 5,
	    "battery_energy_discount" : 5,
	    "battery_power_discount" : 5,
	    "hydrogen_energy_discount" : 5,
	    "hydrogen_electrolyser_discount" : 5,
	    "hydrogen_turbine_discount" : 5,
	    "dispatchable1_discount" : 10,
	       "dispatchable2_discount" : 10,
      	    "wind" : true,
	    "solar" : true,
	    "battery" : true,
	    "hydrogen" : true,
	    "dispatchable1" : false,
	    "dispatchable2" : false,
	      };




let tech_assumptions = {"2020" : {"wind_cost" : 1120,
				  "wind_fom" : 3,
				  "wind_lifetime" : 25,
				  "solar_cost" : 420,
				  "solar_fom" : 3,
				  "solar_lifetime" : 25,
				  "battery_energy_cost" : 232,
				  "battery_energy_fom" : 3,
				  "battery_energy_lifetime" : 15,
				  "battery_power_cost" : 270,
				  "battery_power_fom" : 3,
				  "battery_power_lifetime" : 15,
				  "battery_power_efficiency_charging" : 95,
				  "battery_power_efficiency_discharging" : 95,
				  "hydrogen_energy_cost" : 0.7,
				  "hydrogen_energy_fom" : 14,
				  "hydrogen_energy_lifetime" : 25,
				  "hydrogen_electrolyser_cost" : 1100,
				  "hydrogen_electrolyser_efficiency" : 58,
				  "hydrogen_electrolyser_fom" : 3,
				  "hydrogen_electrolyser_lifetime" : 20,
				  "hydrogen_turbine_cost" : 880,
				  "hydrogen_turbine_efficiency" : 56,
				  "hydrogen_turbine_fom" : 3,
				  "hydrogen_turbine_lifetime" : 25,
				  "dispatchable1_cost" : 400,
				  "dispatchable1_marginal_cost" : 50,
				  "dispatchable1_emissions" : 500,
				  "dispatchable1_fom" : 3,
				  "dispatchable1_lifetime" : 25,
				  "dispatchable2_cost" : 6000,
				  "dispatchable2_marginal_cost" : 10,
				  "dispatchable2_emissions" : 0,
				  "dispatchable2_fom" : 3,
				  "dispatchable2_lifetime" : 25,
				 },
			"2030" : {"wind_cost" : 1040,
				  "wind_fom" : 3,
				  "wind_lifetime" : 25,
				  "solar_cost" : 300,
				  "solar_fom" : 3,
				  "solar_lifetime" : 25,
				  "battery_energy_cost" : 142,
				  "battery_energy_fom" : 3,
				  "battery_energy_lifetime" : 15,
				  "battery_power_cost" : 160,
				  "battery_power_fom" : 3,
				  "battery_power_lifetime" : 15,
				  "battery_power_efficiency_charging" : 95,
				  "battery_power_efficiency_discharging" : 95,
				  "hydrogen_energy_cost" : 0.7,
				  "hydrogen_energy_fom" : 14,
				  "hydrogen_energy_lifetime" : 25,
				  "hydrogen_electrolyser_cost" : 600,
				  "hydrogen_electrolyser_efficiency" : 62,
				  "hydrogen_electrolyser_fom" : 3,
				  "hydrogen_electrolyser_lifetime" : 20,
				  "hydrogen_turbine_cost" : 830,
				  "hydrogen_turbine_efficiency" : 58,
				  "hydrogen_turbine_fom" : 3,
				  "hydrogen_turbine_lifetime" : 25,
				  "dispatchable1_cost" : 400,
				  "dispatchable1_marginal_cost" : 50,
				  "dispatchable1_emissions" : 500,
				  "dispatchable1_fom" : 3,
				  "dispatchable1_lifetime" : 25,
				  "dispatchable2_cost" : 6000,
				  "dispatchable2_marginal_cost" : 10,
				  "dispatchable2_emissions" : 0,
				  "dispatchable2_fom" : 3,
				  "dispatchable2_lifetime" : 25,
				 },
			"2050" : {"wind_cost" : 960,
				  "wind_fom" : 3,
				  "wind_lifetime" : 25,
				  "solar_cost" : 240,
				  "solar_fom" : 3,
				  "solar_lifetime" : 25,
				  "battery_energy_cost" : 75,
				  "battery_energy_fom" : 3,
				  "battery_energy_lifetime" : 15,
				  "battery_power_cost" : 60,
				  "battery_power_fom" : 3,
				  "battery_power_lifetime" : 15,
				  "battery_power_efficiency_charging" : 95,
				  "battery_power_efficiency_discharging" : 95,
				  "hydrogen_energy_cost" : 0.7,
				  "hydrogen_energy_fom" : 14,
				  "hydrogen_energy_lifetime" : 25,
				  "hydrogen_electrolyser_cost" : 400,
				  "hydrogen_electrolyser_efficiency" : 67,
				  "hydrogen_electrolyser_fom" : 3,
				  "hydrogen_electrolyser_lifetime" : 20,
				  "hydrogen_turbine_cost" : 800,
				  "hydrogen_turbine_efficiency" : 60,
				  "hydrogen_turbine_fom" : 3,
				  "hydrogen_turbine_lifetime" : 25,
				  "dispatchable1_cost" : 400,
				  "dispatchable1_marginal_cost" : 50,
				  "dispatchable1_emissions" : 500,
				  "dispatchable1_fom" : 3,
				  "dispatchable1_lifetime" : 25,
				  "dispatchable2_cost" : 6000,
				  "dispatchable2_marginal_cost" : 10,
				  "dispatchable2_emissions" : 0,
				  "dispatchable2_fom" : 3,
				  "dispatchable2_lifetime" : 25,
				 }
		       };


let default_tech_scenario = "2030";

for (let i = 0; i < Object.keys(tech_assumptions[default_tech_scenario]).length; i++){
    let key = Object.keys(tech_assumptions[default_tech_scenario])[i];
    if(!(key in assumptions)){
	assumptions[key] = tech_assumptions[default_tech_scenario][key];
    };
};


d3.select("#tech_scenario").on("change", function(){
    let scenario = this.value;
    console.log("tech scenario change to ",scenario);
    for (let i = 0; i < Object.keys(tech_assumptions[scenario]).length; i++){
	let key = Object.keys(tech_assumptions[scenario])[i];
	let value = tech_assumptions[scenario][key];
	assumptions[key] = value;
	document.getElementsByName(key)[0].value = value;
    };
});



let locations = ["destination", "source"];

for (let i=0; i<locations.length; i++){
    let location = locations[i];
    let lat = assumptions[location]["lat"];
    let lng = assumptions[location]["lng"];
    let marker = L.marker([lat,lng], {draggable:'true'});
    map.addLayer(marker);
    marker.bindPopup('efuel ' + location);
    document.getElementsByName(location+"_location")[0].value = printLocation(lat,lng);
    marker.on('dragend', function(){
	lat = marker.getLatLng().lat;
	lng = marker.getLatLng().lng;
	assumptions[location]["lat"] = lat;
	assumptions[location]["lng"] = lng;
	document.getElementsByName(location+"_location")[0].value = printLocation(lat,lng);
	console.log(location, 'changed location to: lat:', lat.toFixed(1), ", lon:", lng.toFixed(1));
    });
};



function printLocation(lat,lng){
    return "latitude: " + lat.toFixed(1) + ", longitude: " + lng.toFixed(1);
};




for (let i = 0; i < Object.keys(assumptions).length; i++){
    let key = Object.keys(assumptions)[i];
    let value = assumptions[key];
    console.log(key,value);
    if(typeof value === "boolean"){
	document.getElementsByName(key)[0].checked = value;
	d3.selectAll("input[name='" + key + "']").on("change", function(){
	    assumptions[key] = this.checked;
	    console.log(key,"changed to",assumptions[key]);
	});
    }
    else if(["job_type","location","version","jobid","timestamp","queue_length","weather_hex","results_hex"].includes(key)||locations.includes(key)){
    }
    else{
	document.getElementsByName(key)[0].value = value;
	d3.selectAll("input[name='" + key + "']").on("change", function(){
	    assumptions[key] = this.value;
	    console.log(key,"changed to",assumptions[key]);
	});
    }
};
