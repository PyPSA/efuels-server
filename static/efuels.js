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
// https://github.com/PyPSA/whobs-server


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



assumptions = {"year": 2012,
	       "cf_exponent": 2,
	       "destination" : {"lat": 50.11, "lng" : 8.68},
	       "source" : {"lat": -30.8, "lng" : 18.5}
	      };

var locations = ["destination", "source"];

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



for (let i=0; i< Object.keys(assumptions).length; i++){
    let key = Object.keys(assumptions)[i];
    let value = assumptions[key];
    console.log(key,value);
    if(!locations.includes(key)){
	document.getElementsByName(key)[0].value = value;
    };
};
