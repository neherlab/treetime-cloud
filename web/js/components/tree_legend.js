var EventEmitter = require('events').EventEmitter;
var d3 = require('d3');

var Globals = require('./globals.js')
var colors = Globals.colors;
var ANIMATION_DURATION = Globals.ANIMATION_DURATION;


        // var tmp_leg = canvas.selectAll(".legend")
        //     .data(cScale.domain())
        //     .enter().append('g')
        //     .attr('class', 'legend')
        //     .attr('transform', function(d, i) {
        //         var stack = 5;
        //         var height = legendRectSize + legendSpacing;
        //         var fromRight = Math.floor(i / stack);
        //         var fromTop = i % stack;
        //         var horz = fromRight * 145 + 5;
        //         var vert = fromTop * height + 5;
        //         return 'translate(' + horz + ',' + vert + ')';
        //      });
        // tmp_leg.append('rect')
        //     .attr('width', legendRectSize)
        //     .attr('height', legendRectSize)
        //     .style('fill', function (d) {
        //         var col = cScale(d);
        //         return d3.rgb(col).brighter([0.35]).toString();
        //      })
        //     .style('stroke', function (d) {
        //         var col = cScale(d);
        //         return d3.rgb(col).toString();
        //     })
        //     .on('mouseover', mouseover_func)
        //     .on('mouseout', mouseout_func)
        //     .on('click', click_func);

        // tmp_leg.append('text')
        //     .attr('x', legendRectSize + legendSpacing + 5)
        //     .attr('y', legendRectSize - legendSpacing)
        //     .text(label_fmt_func)
        //     .on('mouseover', mouseover_func)
        //     .on('mouseout', mouseout_func)
        //     .on('click', click_func);
        // }


var TreeLegend = {};

TreeLegend.create = function(el, props, state){
    
    console.log("CREATING legend");

    var svg = d3.select(el).append('svg')
      .attr('class', 'd3')
      .attr('width',  props.width)
      .attr('height', props.height);
    
    svg.append('g')
        .attr('class', 'd3-legend-grid')

    var dispatcher = new EventEmitter();
    
    this.update(el, state);
//    this.update(el, state, dispatcher);
    return dispatcher;
};

TreeLegend.update = function (el, state){

    console.log("UPDATING legend");
    var width = el.offsetWidth;
    var height = el.offsetHeight;
    var cScale = state.cscale;
    var m = cScale.get_cmap();
    var legendSpacing = 6
    var nRows = d3.round((m.length + 0.5)/2);
    var rowHeight = d3.round(height / nRows);

    var grid_data = m.map(function (d, i){
        
        var x = i >= nRows ? d3.min([100, width/2]) : 0;
        var y = i >= nRows ? rowHeight*(i - nRows) : rowHeight * (i);
        
        return ({
            color: d.color,
            value: d.value,
            x: x,
            y: y,
            h: rowHeight - legendSpacing,
            w: rowHeight - legendSpacing
        });
    });

    var g = d3.select(el).selectAll('.d3-legend-grid');

    var legend = g.selectAll('.d3-legend-entry')
        .data(grid_data)
        .enter()
        .append('g');

    legend.append('rect')
        .attr('width', function(d){return d.w})
        .attr('height', function(d){return d.h})
        .attr("x", function(d){return d.x})
        .attr("y", function(d){return d.y})
        .style('fill', function(d) {return d.color;});

    legend.append('text')
        .text(function(d){return d.value;})
        .attr('x', function(d){return d.x + rowHeight + legendSpacing })
        .attr('y', function(d) {return d.y  + rowHeight - legendSpacing})

};

export default TreeLegend;
