var EventEmitter = require('events').EventEmitter;
var d3 = require('d3');
var Globals = require('./globals.js')
var ANIMATION_DURATION = Globals.ANIMATION_DURATION;

var RootLhPlot = {};

RootLhPlot.lh = [];

RootLhPlot.padding = 30;
RootLhPlot.left_padding = 50;

RootLhPlot.create = function(el, props, state){
    console.log("CREATING lh PLOT")
    var svg = d3.select(el).append('svg')
      .attr('class', 'd3')
      .attr('width', props.width)
      .attr('height', props.height);
    
    svg.append('g')
      .attr('class', 'd3_lh_axes')

    svg.append('g')
      .attr('class', 'd3_lh_points')

    this.lh = props.lh
    
    this._update_points();
    var dispatcher = new EventEmitter();
    this.update(el, this.lh, state, dispatcher);
    return dispatcher;

};

RootLhPlot._update_points = function(){
    this.points = this.lh;
};

RootLhPlot.update = function(el, lh, state, dispatcher){
    
    console.log("UPDATING lh");
    
    this.lh = lh;
    this._update_points();

    var scales = this._scales(el);
    this._draw_axis(el, scales)
    this._draw_points(el, scales, dispatcher)
    
    // selected node
};

RootLhPlot._scales = function(el){

  
  var width = el.offsetWidth;
  var height = el.offsetHeight;
  var xs = this.points.map(function(d){return d.x});
  var ys = this.points.map(function(d){return d.y});
  var x = d3.scale.linear()
    .domain([d3.min(xs) , d3.max(xs)])
    .range([this.left_padding, width-this.padding]);

  var y = d3.scale.linear()
      .domain([d3.max(ys), d3.min(ys)])
      .range([this.padding,height-this.padding])
  return {x: x, y: y};

};

RootLhPlot._tipRadius = function(d){
    return d.selected ? 8.0 : 6.0;
};

RootLhPlot._draw_axis = function(el, scales){
    
    var width = el.offsetWidth;
    var height = el.offsetHeight;

    var xAxis = d3.svg.axis()
        .scale(scales.x)
        .orient("bottom")
        .ticks(5)
    
    var yAxis = d3.svg.axis()
        .scale(scales.y)
        .orient("left")
        .ticks(10)

    // function for the y grid lines
    function make_x_axis() {
    return d3.svg.axis()
      .scale(scales.x)
      .orient("bottom")
      .ticks(5)
    }

    function make_y_axis() {
    return d3.svg.axis()
      .scale(scales.y)
      .orient("left")
      .ticks(10)
    }

    var svg = d3.select(el).append("svg")
        .attr("class", "d3_lh_axis") 
        .attr("width", width)
        .attr("height", height)

    svg.append("g")
        .attr("class", "d3_lh_axis")
        .attr("transform", "translate(0," + (height - this.padding) + ")")
        .call(xAxis)

    svg.append("g")
        .attr("class", "d3_lh_grid")
        .attr("transform", "translate(0," + (this.padding) + ")")
        .call(make_x_axis()
            .tickSize(height-2*this.padding, 0, 0)
            .tickFormat("")
            )
    
    svg.append("g")
        .attr("transform", "translate(" + (this.left_padding) + ",0)")
        .call(yAxis);
    
    svg.append("g")
        .attr("class", "d3_lh_grid")
        .attr("transform", "translate(" + (width - this.padding) + ",0)")
        .call(make_y_axis()
            .tickSize(width-this.padding-this.left_padding, 0, 0)
            .tickFormat("")
            )
};

RootLhPlot._draw_points = function(el, scales, dispatcher){

    console.log("DRAW POINTS")
    var g = d3.select(el).selectAll('.d3_lh_points');
    var tip = g.selectAll('.d3-lh-point')
          .data(this.points);

    tip.enter()
      .append("circle")
      .attr("class", "d3-lh-point")
      .attr("id", function(d){return d.x})
    
    tip
      .attr("cx", function(d){return scales.x(d.x)})
      .attr("cy", function(d){return scales.y(d.y)})
      .attr("r", 3)
      .style("fill", "#4D92BF")
      .on('mouseover', function(d) {
          dispatcher.emit('point:point_mouseover', d);
      })
      .on('mouseout', function(d) {
          dispatcher.emit('point:point_mouseout', d);
      })

    var line = d3.svg.line()
      .x(function(d) { return scales.x(d.x); })
      .y(function(d) { return scales.y(d.y); });
    
    var svg = d3.select(el).append("svg")
      .append("g")

    svg.append("path")
      .datum(this.points)
      .attr("class", "d3_lh_line")
      .attr("d", line);
};

export default RootLhPlot;