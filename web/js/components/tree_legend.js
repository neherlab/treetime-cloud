var EventEmitter = require('events').EventEmitter;
var d3 = require('d3');

var Globals = require('./globals.js')
var colors = Globals.colors;
var ANIMATION_DURATION = Globals.ANIMATION_DURATION;


var TreeLegend = {};

TreeLegend.create = function(el, props, state){
    
    console.log("CREATING legend");
    var width = el.offsetWidth;
    var height = el.offsetHeight;

    var svg = d3.select(el).append('svg')
      .attr('class', 'd3-legend')
      .attr('width',  width)
      .attr('height', height);
    console.log(svg);

    
    svg.append('g')
        .attr('class', 'd3-legend-grid')

    var dispatcher = new EventEmitter();
    
//    this.update(el, state);
    this.update(el, state, dispatcher);
    return dispatcher;
};

TreeLegend._style={};

TreeLegend.state = {};

TreeLegend.update = function (el, state, dispatcher){

    console.log("UPDATING legend");
    if (this.state.cscale != state.cscale){
        this.state.cscale = state.cscale;
        this._draw_legend(el, dispatcher);        
    }    
};

TreeLegend._draw_legend = function(el, dispatcher){
    
    console.log ("Redraw legend...")
    var width = el.offsetWidth;
    var height = el.offsetHeight;
    var cScale = this.state.cscale;
    if (!cScale) return;

    var m = cScale.get_cmap();
    var size = this._style.rect_size ? this._style.rect_size : 20;
    var legendSpacing = 6
    var nRows = d3.round((m.length + 0.5)/2);
    var rowHeight = size + legendSpacing;

    var grid_data = m.map(function (d, i){
        
        var x = i >= nRows ? d3.min([150, width/2]) : 0;
        var y = i >= nRows ? rowHeight*(i - nRows) : rowHeight * (i);
        
        return ({
            color: d.color,
            value: d.value,
            x: x,
            y: y,
            h: size, 
            w: size
       });
   
    });

    var g = d3.select(el).select('.d3-legend').selectAll('.d3-legend-grid');

    var legend = g.selectAll('.d3-legend-entry')
        .data(grid_data)
        .enter()
        .append('g');

    legend.append('rect')
        .attr('class', 'd3_legend_rect')
        .attr('width', function(d){return d.w})
        .attr('height', function(d){return d.h})
        .attr("x", function(d){return d.x})
        .attr("y", function(d){return d.y})
        .style('fill', function(d) {return d.color;});

    legend.append('text')
        .text(function(d){return d.value;})
        .attr('x', function(d){return d.x + d.w + legendSpacing })
        .attr('y', function(d) {return d.y + d.h / 2 + 4})

};

export default TreeLegend;
