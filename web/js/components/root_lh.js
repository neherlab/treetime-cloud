var EventEmitter = require('events').EventEmitter;
var d3 = require('d3');
var Globals = require('./globals.js')
var ANIMATION_DURATION = Globals.ANIMATION_DURATION;

var RootLhPlot = {};

RootLhPlot.lh = [];

RootLhPlot.padding_bottom = 80;

RootLhPlot.padding_text = 20;

RootLhPlot.padding_top = 10;

RootLhPlot.padding_right = 10;

RootLhPlot.padding_left = 120;

RootLhPlot.create = function(el, props, state){

    this.width = el.offsetWidth
    this.height = el.offsetHeight

    //console.log("CREATING lh PLOT")
    this.svg = d3.select(el).append('svg')
      .attr('class', 'd3_lh')
      .attr('width', el.offsetWidth)
      .attr('height', el.offsetHeight);

    ////console.log(svg)
    this.svg.append('svg')
      .attr('class', 'd3_lh_axis')

    this.svg.append('g')
      .attr('class', 'd3_lh_points')

    var dispatcher = new EventEmitter();
    this.update(el, this.lh, state, dispatcher);
    return dispatcher;

};

RootLhPlot._update_points = function(){
    this.points = this.lh;
};

RootLhPlot.update = function(el, lh, state, dispatcher){

    if (this.points != lh){

      // update the whole plot
      this.lh = lh;
      this._update_points();

      var scales = this._scales(el);
      this._draw_axis(el, scales)
      this._draw_points(el, scales, dispatcher)

    }

    if (this.width != el.offsetHeight || this.height != el.offsetHeight){

      var g = d3.select(el).select('.d3_lh_axis').selectAll("*");
      g.remove();
      g = d3.select(el).select('.d3_lh_points').selectAll("*");
      g.remove();

      this.svg
        .attr("width", el.offsetWidth)
        .attr("height", el.offsetHeight);

      var scales = this._scales(el);
      this._draw_axis(el, scales)
      this._draw_points(el, scales, dispatcher)
      this.width = el.offsetWidth;
      this.height = el.offsetHeight;

    }


    // selected node
};

RootLhPlot._scales = function(el){


  var width = el.offsetWidth;
  var height = el.offsetHeight;
  var xs = this.points.map(function(d){return d.x});
  var ys = this.points.map(function(d){return d.y});

  var x = d3.scale.linear()
    .domain([d3.min(xs) , d3.max(xs)])
    .range([this.padding_left, width-this.padding_right]);

  var y = d3.scale.linear()
      .domain([d3.max(ys), d3.min(ys)])
      .range([this.padding_top,height-this.padding_bottom])
  return {x: x, y: y};

};

RootLhPlot._tipRadius = function(d){
    return d.selected ? 4.0 : 4.0;
};

RootLhPlot._draw_axis = function(el, scales){

    var width = el.offsetWidth;
    var height = el.offsetHeight;
    var h = -height + 30

    var xAxis = d3.svg.axis()
        .scale(scales.x)
        .orient("bottom")
        .tickValues(scales.x.ticks(5).concat(scales.x.domain()[1]))
        .tickFormat(d3.format("d"))
        .tickSize(h, 0, 0)
    
    var yAxis = d3.svg.axis()
        .scale(scales.y)
        .orient("left")
        .tickValues(scales.y.ticks(10).concat(scales.y.domain()[1]))
        .tickFormat(d3.format(".1f"))
        .tickSize(-width, 0, 0)

    
    var svg = d3.select(el).select('.d3_lh_axis')


    svg.append("g")
        .attr("class", "d3_lh_x_axis")
        .attr("transform", "translate(0," + (height -  this.padding_bottom) + ")")
        .call(xAxis)

    var axis_start = scales.x.range()[0];
    var axis_width = scales.x.range()[1]-scales.x.range()[0];
    svg.append("text")      // text label for the x axis
        .attr("x", axis_start + axis_width / 2 )
        .attr("y", height - this.padding_text )
        .style("text-anchor", "middle")
        .text("Inferred root date");


    svg.append("g")
        .attr("class", "d3_lh_y_axis")
        .attr("transform", "translate(" + (this.padding_left) + ",0)")
        .call(yAxis);

    svg.append("text")
        .attr("transform", "rotate(-90)")
        .attr("y", this.padding_text)
        .attr("x",  - (height - this.padding_bottom)/2)
        .attr("dy", "1em")
        .style("text-anchor", "middle")
        .text("Normalized likelihood");

};

RootLhPlot._draw_points = function(el, scales, dispatcher){

    //console.log("DRAW POINTS")
    var g = d3.select(el).selectAll('.d3_lh_points');
    var tip = g.selectAll('.d3_lh-point')
          .data(this.points);

    tip.enter()
      .append("circle")
      .attr("class", "d3_lh-point")
      .attr("id", function(d){return d.x})

    tip
      .attr("cx", function(d){return scales.x(d.x)})
      .attr("cy", function(d){return scales.y(d.y)})
      .attr("r", this._tipRadius)
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

    g.append("path")
      .datum(this.points)
      .attr("class", "d3_lh_line")
      .attr("d", line);
};

export default RootLhPlot;