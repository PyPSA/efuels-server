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


var licenceText = '(Licence: <a href="https://creativecommons.org/licenses/by/4.0/">CC BY 4.0</a>, Attribution: <a href="https://model.energy">model.energy</a> & <a href="https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5">ECMWF ERA5</a> via <a href="https://cds.climate.copernicus.eu/">Copernicus Climate Change Service</a> (see <a href="https://apps.ecmwf.int/datasets/licences/copernicus/">weather data licence</a>))';

var parseDate = d3.timeParse("%Y-%m-%d %H:%M:00");

var formatDate = d3.timeFormat("%b %d %H:%M");



var colors = {"wind":"#3B6182",
              "solar" :"#FFFF00",
              "battery" : "#999999",
              "battery_power" : "#999999",
              "battery_energy" : "#666666",
              "hydrogen_turbine" : "red",
              "hydrogen_electrolyser" : "cyan",
              "hydrogen_energy" : "magenta",
	      "dispatchable1" : "orange",
	      "dispatchable2" : "lime",
	      "hydrogen_submarine_pipeline" : "green",
	      "load" : "black",
	      "efuels_load" : "purple",
	      "hydrogen_load" : "purple",
             };


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


let assumptions = {};

for (let key in defaults){
    assumptions[key] = defaults[key]["value"];
}

let tech_assumptions = {"2020" : {},
			"2030" : {},
			"2050" : {},
		       };

let default_tech_scenario = "2030";

for (let key in defaults_t[default_tech_scenario]){
    assumptions[key] = defaults_t[default_tech_scenario][key]["value"];
    for (let year in tech_assumptions){
	tech_assumptions[year][key] = defaults_t[year][key]["value"];
    }
}


let assets = ["solar","wind","battery_power",
	  "battery_energy","hydrogen_electrolyser",
	  "hydrogen_turbine","hydrogen_energy",
	      "dispatchable1","dispatchable2","hydrogen_submarine_pipeline"];

let vre = ["solar","wind"];

let electricity = ["solar","wind","dispatchable1","dispatchable2"];

let storage = ["battery","hydrogen"];


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
    let lat = assumptions[location +"_lat"];
    let lng = assumptions[location +"_lng"];
    let marker = L.marker([lat,lng], {draggable:'true'});
    map.addLayer(marker);
    marker.bindPopup('efuel ' + location);
    document.getElementsByName(location+"_location")[0].value = printLocation(lat,lng);
    marker.on('dragend', function(){
	lat = marker.getLatLng().lat;
	lng = marker.getLatLng().lng;
	assumptions[location +"_lat"] = lat;
	assumptions[location +"_lng"] = lng;
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
    else if(["job_type","source_lat","source_lng","destination_lat","destination_lng","version","jobid","timestamp","queue_length","weather_hex","results_hex"].includes(key)){
    }
    else{
	document.getElementsByName(key)[0].value = value;
	d3.selectAll("input[name='" + key + "']").on("change", function(){
	    assumptions[key] = this.value;
	    console.log(key,"changed to",assumptions[key]);
	});
    }
}


var solveButton = d3.select("#solve-button");

var solveButtonText = {"before" : "Solve",
		       "after" : "Solving"};

var jobid = "";

var timer;
var timeout;
var timerStart;
var timerExpected = 15;


// time between status polling in milliseconds
var poll_interval = 500;

// time out for polling if it doesn't finish after 10 minutes
// Shouldn't be divisible by poll_interval
var poll_timeout = 10*60*1000 + poll_interval/2;



function solve() {
    if (solveButton.text() == solveButtonText["before"]) {
	clear_results();
	var send_job = new XMLHttpRequest();
	send_job.open('POST', './jobs', true);
	send_job.setRequestHeader("Content-Type", "application/json");
	send_job.onload = function () {
	    var data = JSON.parse(this.response);
	    if("jobid" in data){
		jobid = data["jobid"];
		console.log("Jobid:",jobid);
		timer = setInterval(poll_result, poll_interval);
		timerStart = new Date().getTime();
		document.getElementById("countdown").innerHTML="Ready in around " + timerExpected + " seconds";
		console.log("timer",timer,"polling every",poll_interval,"milliseconds");
		timeout = setTimeout(poll_kill, poll_timeout);
		solveButton.text(solveButtonText["after"]);
		solveButton.attr("disabled","");
		document.getElementById("status").innerHTML="Sending job to solver";
	    } else if("status" in data && data["status"] == "Error") {
		console.log("results:", data);
		document.getElementById("status").innerHTML = data["status"] + ": " + data["error"];
	    } else {
		console.log("results:", data);
		document.getElementById("status").innerHTML = "Found previously calculated version";
		results = data;
		display_results();
	    };
	};
	assumptions["job_type"] = "solve";
	send_job.send(JSON.stringify(assumptions));
    };
};


solveButton.on("click", solve);


function poll_result() {

    var poll = new XMLHttpRequest();

    poll.open('GET', './jobs/' + jobid, true);

    poll.onload = function () {
	results = JSON.parse(this.response);
	status = results["status"];
	document.getElementById("status").innerHTML=status;
	console.log("status is",status);

	document.getElementById("countdown").innerHTML = "Ready in around " + Math.round(timerExpected - (new Date().getTime() - timerStart)/1000.) + " seconds";

	if(status == "Error"){
	    clearInterval(timer);
	    clearTimeout(timeout);
	    document.getElementById("countdown").innerHTML = "Ready in around " + timerExpected + " seconds";
	    console.log("results:",results);
	    document.getElementById("status").innerHTML=status + ": " + results["error"];
	    solveButton.text(solveButtonText["before"]);
	    $('#solve-button').removeAttr("disabled");
	};
	if(status == "Finished"){
	    clearInterval(timer);
	    clearTimeout(timeout);
	    document.getElementById("countdown").innerHTML = "Solved in " + Math.round((new Date().getTime() - timerStart)/1000.) + " seconds";
	    console.log("results:",results);
	    solveButton.text(solveButtonText["before"]);
	    $('#solve-button').removeAttr("disabled");
	    display_results();
	};
    };
    poll.send();
};


function poll_kill() {
    clearInterval(timer);
    solveButton.text(solveButtonText["before"]);
    $('#solve-button').removeAttr("disabled");
    document.getElementById("status").innerHTML="Error: Timed out";
};




function clear_results(){
    document.getElementById("results_assumptions").innerHTML="";
    document.getElementById("average_cost").innerHTML="";
    document.getElementById("distance").innerHTML="";
    for (let i = 0; i < assets.length; i++){
	document.getElementById(assets[i] + "_capacity").innerHTML="";
	document.getElementById(assets[i] + "_cf_used").innerHTML="";
	if(!assets[i].includes("energy")){
	    document.getElementById(assets[i] + "_rmv").innerHTML="";
	};
    };
    document.getElementById("battery_discharge_rmv").innerHTML="";
    for (let i = 0; i < vre.length; i++){
	document.getElementById(vre[i] + "_cf_available").innerHTML="";
	document.getElementById(vre[i] + "_curtailment").innerHTML="";
    };
    for (let i = 0; i < electricity.length; i++){
	document.getElementById(electricity[i] + "_lcoe").innerHTML="";
    };
    for (let i = 0; i < storage.length; i++){
	document.getElementById(storage[i] + "_lcos").innerHTML="";
    };
    d3.select("#power").selectAll("g").remove();
    d3.select("#average_cost_graph").selectAll("g").remove();
    d3.select("#power_capacity_bar").selectAll("g").remove();
    d3.select("#energy_capacity_bar").selectAll("g").remove();

    document.getElementById("results-overview-download").innerHTML = '';
    document.getElementById("results-series-download").innerHTML = '';
    document.getElementById("results-link").innerHTML = '';

};




function display_results(){

    $('#collapseResults').addClass("show");

    document.getElementById("results_assumptions").innerHTML=" for weather year " + results["assumptions"]["year"];
    document.getElementById("average_cost").innerHTML="Average cost [EUR/MWh]: " + results["average_cost"].toFixed(1) + ", [EUR/kg]: " + (results["average_cost"]*0.033).toFixed(2);
    document.getElementById("distance").innerHTML="Distance as crow flies [km]: " + results["distance"].toFixed(0);

    for (let i = 0; i < assets.length; i++){
	document.getElementById(assets[i] + "_capacity").innerHTML=Math.abs(results[assets[i] + "_capacity"].toFixed(1));
	document.getElementById(assets[i] + "_cf_used").innerHTML=Math.abs((results[assets[i] + "_cf_used"]*100)).toFixed(1);
	if(!assets[i].includes("energy")){
	    document.getElementById(assets[i] + "_rmv").innerHTML=Math.abs((results[assets[i] + "_rmv"]*100)).toFixed(1);
	};
    };
    document.getElementById("battery_discharge_rmv").innerHTML=Math.abs((results["battery_discharge_rmv"]*100)).toFixed(1);
    for (let i = 0; i < vre.length; i++){
	document.getElementById(vre[i] + "_cf_available").innerHTML=Math.abs((results[vre[i] + "_cf_available"]*100)).toFixed(1);
	document.getElementById(vre[i] + "_curtailment").innerHTML=Math.abs((results[vre[i] + "_curtailment"]*100)).toFixed(1);
    };

    for (let i = 0; i < electricity.length; i++){
	if(results[electricity[i] + "_capacity"] > 0.1){
	    var cost = results[electricity[i]+"_cost"];
	    if(results.hasOwnProperty(electricity[i]+"_marginal_cost")){
		cost += results[electricity[i]+"_marginal_cost"];
	    };
	    document.getElementById(electricity[i] + "_lcoe").innerHTML=Math.abs(cost/(results[electricity[i] + "_capacity"]*results[electricity[i] + "_cf_used"])).toFixed(1);
	};
    };

    // all per hour
    if (results["battery_power_capacity"] > 0.1){
	let battery_fedin = results["battery_power_capacity"]*results["battery_power_cf_used"];
	let battery_charged = battery_fedin/0.95**2;
	document.getElementById("battery_lcos").innerHTML=Math.abs((results["battery_power_cost"]+results["battery_energy_cost"]+battery_charged*results["battery_power_rmv"]*results["average_price"])/battery_fedin).toFixed(1);
    };
    if (results["hydrogen_turbine_capacity"] > 0.1){
	let hydrogen_fedin = results["hydrogen_turbine_capacity"]*results["hydrogen_turbine_cf_used"];
	let hydrogen_charged = hydrogen_fedin/(results["assumptions"]["hydrogen_turbine_efficiency"]/100)/(results["assumptions"]["hydrogen_electrolyser_efficiency"]/100);
	document.getElementById("hydrogen_lcos").innerHTML=Math.abs((results["hydrogen_turbine_cost"]+results["hydrogen_energy_cost"]+results["hydrogen_electrolyser_cost"]+hydrogen_charged*results["hydrogen_electrolyser_rmv"]*results["average_price"])/hydrogen_fedin).toFixed(1);
    };

    for(var j=0; j < results.snapshots.length; j++) {
	results.snapshots[j] = parseDate(results.snapshots[j]);
    };

    draw_power_graph();
    draw_power_capacity_bar();
    draw_energy_capacity_bar();
    draw_cost_stack();

    document.getElementById("results-overview-download").innerHTML = '<a href="data/results-overview-' + results.assumptions.results_hex + '.csv">Download Comma-Separated-Variable (CSV) file of results overview</a> ' + licenceText;
    document.getElementById("results-series-download").innerHTML = '<a href="data/results-series-' + results.assumptions.results_hex + '.csv">Download Comma-Separated-Variable (CSV) file of results time series</a> ' + licenceText;
    document.getElementById("results-link").innerHTML = '<a href="https://model.energy/?results=' + results.assumptions.results_hex + '#solve">Link to these results</a>';
};




function draw_power_graph(){

    let snapshots = results["snapshots"];
    let selection = [...Array(results["snapshots"].length).keys()];

    // Inspired by https://bl.ocks.org/mbostock/3885211

    var svgGraph = d3.select("#power"),
	margin = {top: 20, right: 20, bottom: 110, left: 50},
	marginContext = {top: 430, right: 20, bottom: 30, left: 50},
	width = svgGraph.attr("width") - margin.left - margin.right,
	height = svgGraph.attr("height") - margin.top - margin.bottom,
	heightContext = svgGraph.attr("height") - marginContext.top - marginContext.bottom;

    // remove existing
    svgGraph.selectAll("g").remove();

    var x = d3.scaleTime().range([0, width]).domain(d3.extent(snapshots));
    var y = d3.scaleLinear().range([height, 0]);
    var xContext = d3.scaleTime().range([0, width]).domain(d3.extent(snapshots));
    var yContext = d3.scaleLinear().range([heightContext, 0]);

    var xAxis = d3.axisBottom(x),
	xAxisContext = d3.axisBottom(xContext),
	yAxis = d3.axisLeft(y);


    var brush = d3.brushX()
        .extent([[0, 0], [width, heightContext]])
        .on("start brush end", brushed);


    var zoom = d3.zoom()
        .scaleExtent([1, Infinity])
        .translateExtent([[0, 0], [width, height]])
        .extent([[0, 0], [width, height]])
        .on("zoom", zoomed);


    var data = [];

    // Custom version of d3.stack

    var previous = new Array(selection.length).fill(0);

    for (var j = 0; j < results["positive"].columns.length; j++){
	var item = [];
	for (var k = 0; k < selection.length; k++){
	    item.push([previous[k], previous[k] + results["positive"]["data"][selection[k]][j]]);
	    previous[k] = previous[k] + results["positive"]["data"][selection[k]][j];
	    }
	data.push(item);
    }
    var previous = new Array(selection.length).fill(0);

    for (var j = 0; j < results["negative"].columns.length; j++){
	var item = [];
	for (var k = 0; k < selection.length; k++){
	    item.push([-previous[k] - results["negative"]["data"][selection[k]][j],-previous[k]]);
	    previous[k] = previous[k] + results["negative"]["data"][selection[k]][j];
	    }
	data.push(item);
    }

    var ymin = 0, ymax = 0;
    for (var k = 0; k < selection.length; k++){
	if(data[results["positive"].columns.length-1][k][1] > ymax){ ymax = data[results["positive"].columns.length-1][k][1];};
	if(data[results["positive"].columns.length+results["negative"].columns.length-1][k][0] < ymin){ ymin = data[results["positive"].columns.length+results["negative"].columns.length-1][k][0];};
    };

    y.domain([ymin,ymax]);
    yContext.domain([ymin,ymax]);

    var area = d3.area()
        .curve(d3.curveMonotoneX)
        .x(function(d,i) { return x(snapshots[i]); })
        .y0(function(d) { return y(d[0]); })
        .y1(function(d) { return y(d[1]); });

    var areaContext = d3.area()
        .curve(d3.curveMonotoneX)
        .x(function(d,i) { return xContext(snapshots[i]); })
        .y0(function(d) { return yContext(d[0]); })
        .y1(function(d) { return yContext(d[1]); });


    svgGraph.append("defs").append("clipPath")
        .attr("id", "clip")
	.append("rect")
        .attr("width", width)
        .attr("height", height);

    var focus = svgGraph.append("g")
        .attr("class", "focus")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    var context = svgGraph.append("g")
        .attr("class", "context")
        .attr("transform", "translate(" + marginContext.left + "," + marginContext.top + ")");

    var layer = focus.selectAll(".layer")
        .data(data)
        .enter().append("g")
        .attr("class", "layer");

    layer.append("path")
        .attr("class", "area")
        .style("fill", function(d, i) {if(i < results["positive"].color.length){ return results["positive"].color[i];} else{return results["negative"].color[i-results["positive"].color.length];} })
        .attr("d", area);

    // add demand curve

    var lineFunction = d3.line()
	.x(function(d) { return x(d[0]) })
	.y(function(d) { return y(d[1]) })
	.curve(d3.curveLinear);

    var demand = focus.append("g");

    //demand.append("path")
     //   .attr("d", lineFunction([[snapshots[0],assumptions["load"]],[snapshots[snapshots.length-1],assumptions["load"]]]))
     //   .attr("id", "indicator")
     //   .attr("stroke", "#000000")
     //   .attr("stroke-width", 3);


    focus.append("g")
        .attr("class", "axis axis--x")
        .attr("transform", "translate(0," + height + ")")
        .call(xAxis);

    focus.append("g")
        .attr("class", "axis axis--y")
        .call(yAxis);


    var label = svgGraph.append("g").attr("class", "y-label");

    // text label for the y axis
    label.append("text")
        .attr("transform", "rotate(-90)")
        .attr("y", 0)
        .attr("x",0 - (height / 2))
        .attr("dy", "1em")
        .style("text-anchor", "middle")
        .text("Power [MW]");


    var layerContext = context.selectAll(".layerContext")
        .data(data)
        .enter().append("g")
        .attr("class", "layerContext");

    layerContext.append("path")
        .attr("class", "area")
        .style("fill", function(d, i) {if(i < results["positive"].color.length){ return results["positive"].color[i];} else{return results["negative"].color[i-results["positive"].color.length];} })
        .attr("d", areaContext);

    context.append("g")
        .attr("class", "axis axis--x")
        .attr("transform", "translate(0," + heightContext + ")")
        .call(xAxisContext);

    var gBrush = context.append("g")
        .attr("class", "brush")
        .call(brush);

    // brush handle follows
    // https://bl.ocks.org/Fil/2d43867ba1f36a05459c7113c7f6f98a

    // following handle looks nicer
    // https://bl.ocks.org/robyngit/89327a78e22d138cff19c6de7288c1cf

    var brushResizePath = function(d) {
	var e = +(d.type == "e"),
	    x = e ? 1 : -1,
	    y = heightContext / 2;
	return "M" + (.5 * x) + "," + y + "A6,6 0 0 " + e + " " + (6.5 * x) + "," + (y + 6) + "V" + (2 * y - 6) + "A6,6 0 0 " + e + " " + (.5 * x) + "," + (2 * y) + "Z" + "M" + (2.5 * x) + "," + (y + 8) + "V" + (2 * y - 8) + "M" + (4.5 * x) + "," + (y + 8) + "V" + (2 * y - 8);
    }

    var handle = gBrush.selectAll(".handle--custom")
	.data([{type: "w"}, {type: "e"}])
	.enter().append("path")
        .attr("class", "handle--custom")
        .attr("stroke", "#000")
        .attr("cursor", "ew-resize")
        .attr("d", brushResizePath);

    gBrush.call(brush.move, x.range()); //this sets initial position of brush


    svgGraph.append("rect")
        .attr("class", "zoom")
        .attr("width", width)
        .attr("height", height)
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")")
        .call(zoom);

    function brushed() {
	if (d3.event.sourceEvent && d3.event.sourceEvent.type === "zoom") return; // ignore brush-by-zoom
	var s = d3.event.selection || xContext.range();
	x.domain(s.map(xContext.invert, xContext));
	layer.attr("d", area);
	focus.select(".axis--x").call(xAxis);
	svgGraph.select(".zoom").call(zoom.transform, d3.zoomIdentity
				      .scale(width / (s[1] - s[0]))
				      .translate(-s[0], 0));
	handle.attr("transform", function(d, i) { return "translate(" + [ s[i], - heightContext / 4] + ")"; });
    }

    function zoomed() {
	if (d3.event.sourceEvent && d3.event.sourceEvent.type === "brush") return; // ignore zoom-by-brush
	var t = d3.event.transform;
	x.domain(t.rescaleX(xContext).domain());
	layer.select(".area").attr("d", area);
	focus.select(".axis--x").call(xAxis);
	var newRange = x.range().map(t.invertX, t);
	context.select(".brush").call(brush.move, newRange);
	handle.attr("transform", function(d, i) { return "translate(" + [ newRange[i], - heightContext / 4] + ")"; });
    }


};




function draw_cost_stack(){

    let data = [];
    let color = [];
    let labels = [];

    for(let i=0; i < assets.length; i++){
	let cost = results[assets[i]+"_cost"];
	if(results.hasOwnProperty(assets[i]+"_marginal_cost")){
	    cost += results[assets[i]+"_marginal_cost"];
	};
	data.push(cost/results["assumptions"]["efuels_load"]);
	color.push(colors[assets[i]]);
	labels.push(assets[i].replaceAll("_"," "));
    };

    draw_stack(data, labels, color, "Breakdown of avg. sys. cost [EUR/MWh]", "#average_cost_graph", " EUR/MWh");
};


function draw_power_capacity_stack(){

    let data = [];
    let color = [];
    let labels = [];

    for(let i=0; i < assets.length; i++){
	if(!assets[i].includes("energy")){
	    data.push(results[assets[i]+"_capacity"]);
	    color.push(colors[assets[i]]);
	    labels.push(assets[i].replaceAll("_"," "));
	};
    };

    draw_stack(data, labels, color, "Power capacity [MW]", "#power_capacity_graph", " MW");
};


function draw_power_capacity_bar(){

    let data = [];
    let color = [];
    let labels = [];

    data.push(results["assumptions"]["efuels_load"]);
    color.push(colors["efuels_load"]);
    labels.push("efuels demand");

    for(let i=0; i < assets.length; i++){
	if(!assets[i].includes("energy") && results[assets[i]+"_capacity"] >= 0.1){
	    data.push(results[assets[i]+"_capacity"]);
	    color.push(colors[assets[i]]);
	    labels.push(assets[i].replaceAll("_"," ").replace("hydrogen","H2"));
	};
    };

    draw_bar(data, labels, color, "Power capacity [MW]", "#power_capacity_bar", " MW");
};



function draw_energy_capacity_stack(){

    let data = [];
    let color = [];
    let labels = [];

    for(let i=0; i < assets.length; i++){
	if(assets[i].includes("energy")){
	    data.push(results[assets[i]+"_capacity"]/1000.);
	    color.push(colors[assets[i]]);
	    labels.push(assets[i].replaceAll("_"," "));
	};
    };

    draw_stack(data, labels, color, "Energy storage capacity [GWh]", "#energy_capacity_graph", " GWh");
};



function draw_energy_capacity_bar(){

    let data = [];
    let color = [];
    let labels = [];

    data.push(results["assumptions"]["efuels_load"]*24/1000.);
    color.push(colors["efuels_load"]);
    labels.push("24h efuels demand");

    for(let i=0; i < assets.length; i++){
	if(assets[i].includes("energy") && results[assets[i]+"_capacity"] >= 10){
	    data.push(results[assets[i]+"_capacity"]/1000.);
	    color.push(colors[assets[i]]);
	    labels.push(assets[i].replace("energy","storage").replaceAll("_"," ").replace("hydrogen","H2"));
	};
    };

    draw_bar(data, labels, color, "Energy capacity [GWh]", "#energy_capacity_bar", " GWh");
};



function draw_energy_stack(){

    let data = [];
    let color = [];
    let labels = [];

    for(let i=0; i < assets.length; i++){
	if(!assets[i].includes("energy")){
	    data.push(results[assets[i]+"_used"]);
	    color.push(colors[assets[i]]);
	    labels.push(assets[i].replaceAll("_"," "));
	};
    };

    draw_stack(data, labels, color, "Average power dispatch [MW]", "#energy_graph", " MW");
};



function draw_stack(data, labels, color, ylabel, svgName, suffix){

    // Inspired by https://bl.ocks.org/mbostock/3885211 and
    // https://bl.ocks.org/mbostock/1134768

    let totals = [0.];

    for(let i=0; i < data.length; i++){
	totals.push(totals[i] + data[i]);
    };


    let svgGraph = d3.select(svgName),
	margin = {top: 20, right: 20, bottom: 30, left: 50},
	width = svgGraph.attr("width") - margin.left - margin.right,
	height = svgGraph.attr("height") - margin.top - margin.bottom;

    // remove existing
    svgGraph.selectAll("g").remove();
    let x = d3.scaleLinear().range([0, width]);
    let y = d3.scaleLinear().range([height, 0]);

    x.domain([0,1]);
    y.domain([0,totals[totals.length-1]]).nice();

    var g = svgGraph.append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    let tip = d3.tip()
	.attr('class', 'd3-tip')
	.offset([-8, 0])
	.html(function(d,i) {
	    return labels[i] + ": " + Math.abs(data[i]).toFixed(1) + suffix;
	});
    svgGraph.call(tip);

    var layer = g.selectAll("rect")
        .data(data)
        .enter().append("rect")
	.attr("x", x(0.1))
        .attr("y", function(d,i) { return y(totals[i+1]);})
        // following abs avoids rect with negative height e.g. -1e10
	.attr("height", function(d,i) { return Math.abs((y(totals[i]) - y(totals[i+1])).toFixed(2)); })
    	.attr("width", x(0.8))
        .style("fill", function(d, i) { return color[i];})
        .on('mouseover', tip.show)
        .on('mouseout', tip.hide);

    //g.append("g")
    //    .attr("class", "axis axis--x")
    //    .attr("transform", "translate(0," + height + ")")
    //    .call(d3.axisBottom(x));

    g.append("g")
        .attr("class", "axis axis--y")
        .call(d3.axisLeft(y));

    var label = svgGraph.append("g").attr("class", "y-label");

    // text label for the y axis
    label.append("text")
        .attr("transform", "rotate(-90)")
        .attr("y", 0)
        .attr("x",0 - (height / 2))
        .attr("dy", "1em")
        .style("text-anchor", "middle")
        .text(ylabel);


    var label = svgGraph.append("g").attr("class", "column-total")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");


    // text label for the y axis
    label.append("text")
        .attr("y", y(totals[totals.length-1])-25)
        .attr("x",x(0.5))
        .attr("dy", "1em")
        .style("text-anchor", "middle")
        .text(totals[totals.length-1].toFixed(1));


};





function draw_bar(data, labels, color, ylabel, svgName, suffix){

    let svgGraph = d3.select(svgName),
	margin = {top: 20, right: 20, bottom: 30, left: 50},
	width = svgGraph.attr("width") - margin.left - margin.right,
	height = svgGraph.attr("height") - margin.top - margin.bottom;

    // remove existing
    svgGraph.selectAll("g").remove();
    let x = d3.scaleBand().range([0, width]).padding(0.1);
    let y = d3.scaleLinear().range([height,0]);

    x.domain(labels);
    y.domain([0,d3.max(data)]).nice();

    var g = svgGraph.append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    let tip = d3.tip()
	.attr('class', 'd3-tip')
	.offset([-8, 0])
	.html(function(d,i) {
	    return labels[i] + ": " + Math.abs(data[i]).toFixed(1) + suffix;
	});
    svgGraph.call(tip);

    var layer = g.selectAll("rect")
        .data(data)
        .enter().append("rect")
	.attr("x", function(d,i) {return x(labels[i]);})
        .attr("y", function(d,i) { return Math.abs(y(d).toFixed(2));})
        // following abs avoids rect with negative height e.g. -1e10
	.attr("height", function(d,i) { return Math.abs((height -y(d)).toFixed(2)); })
    	.attr("width", x.bandwidth())
        .style("fill", function(d, i) { return color[i];})
        .on('mouseover', tip.show)
        .on('mouseout', tip.hide);

    //g.append("g")
    //    .attr("class", "axis axis--x")
    //    .attr("transform", "translate(0," + height + ")")
    //    .call(d3.axisBottom(x));

    g.append("g")
        .attr("class", "axis axis--y")
        .call(d3.axisLeft(y));

    var label = svgGraph.append("g").attr("class", "y-label");

    // text label for the y axis
    label.append("text")
        .attr("transform", "rotate(-90)")
        .attr("y", 0)
        .attr("x",0 - (height / 2))
        .attr("dy", "1em")
        .style("text-anchor", "middle")
        .text(ylabel);

    // add the x Axis
    svgGraph.append("g")
        .attr("transform", "translate(" + margin.left + "," + (height + margin.top) + ")")
        .call(d3.axisBottom(x));

};
