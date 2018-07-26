(function(factory) {
    if (typeof module === 'object' && module.exports) {
        module.exports = factory;
    } else {
        factory(Highcharts);
    }
}(function(Highcharts) {
    (function(Highcharts) {
        'use strict';
        Highcharts.theme = {
            colors: [
                "#3366cc", "#dc3912", "#ff9900", "#109618", "#990099",
                "#0099c6", "#dd4477", "#66aa00", "#b82e2e", "#316395"
            ],
               chart: {
                   backgroundColor: {
                       radialGradient: {cx: 0, cy: 1, r: 1},
                       stops: [
                           [0, '#ffffff'],
                           [1, '#f2f2ff']
                       ]
                   },
                   style: {
                       fontFamily: 'arial, sans-serif',
                       color: '#333'
                   }
               },
               title: {
                   style: {
                       color: '#222',
                       fontSize: '21px',
                       fontWeight: 'bold'
                   }
               },
               subtitle: {
                   style: {
                       fontSize: '16px',
                       fontWeight: 'bold'
                   }
               },
               xAxis: {
                   lineWidth: 1,
                   lineColor: '#cccccc',
                   tickWidth: 1,
                   tickColor: '#cccccc',
                   labels: {
                       style: {
                           fontSize: '12px'
                       }
                   }
               },
               yAxis: {
                   gridLineWidth: 1,
                   gridLineColor: '#d9d9d9',
                   labels: {
                      style: {
                           fontSize: '12px'
                       }
                   }
               },
               legend: {
                   itemStyle: {
                       color: '#666',
                       fontSize: '9px'
                   },
                   itemHoverStyle:{
                       color: '#222'
                   }
               }
        };

        // Apply the theme
        Highcharts.setOptions(Highcharts.theme);

    }(Highcharts));
}));
