var step3_ageDetection = function(){
  
      var maskClouds = function(image){ 
        var cfmask = image.select("pixel_qa");    
        return image.updateMask(cfmask.eq(66));    
      };
      
      var addBSI = function(image) {
        return image
              .addBands(image.expression(
                            '((b6+b4)-(b5+b2))/((b6+b4)+(b5+b2))*1.0', {
                            'b2': image.select('B1'),
                            'b4': image.select('B3'),
                            'b5': image.select('B4'),
                            'b6': image.select('B5')}).rename("bsi"));
                };
      
      var l5 = ee.ImageCollection("LANDSAT/LT05/C01/T1_SR"),
          l7 = ee.ImageCollection("LANDSAT/LE07/C01/T1_SR"),
          danylo2 = ee.ImageCollection("users/OlhaDanylo/palmoil/danylo_v2"),
          grid = ee.FeatureCollection("users/OlhaDanylo/global/tropics/tropics_grid"),
          dem = ee.Image("USGS/SRTMGL1_003");
          
      
      var landsat = l7.merge(l5)
                      .filterDate('1984-01-01', '2018-01-01')
                      .map(maskClouds)
                      .map(addBSI);

      var bsi_min = landsat.select("bsi")
                           .filterDate('2017-01-01', '2018-01-01')
                           .map(function(im){
                                              return im.updateMask(im.select("bsi").gte(-55.42));
                                             })
                           .median();
      
      var im_selected = landsat
                          .select("bsi")
                          .map(function(im){
                                            var diff = im.date().difference(ee.Date('2000-01-01'), 'year');
                                            return ee.Image(im).addBands(ee.Image(diff).float());
                              });
                          
 
      var years = ee.List.sequence(1985, 2018);
      var month = ee.List.sequence(1, 12);
      var dates = years.map(function(year){
                              var cur_dates = month.map(function(month_cur){
                                return ee.Date.fromYMD(year, month_cur, 1);
                              });
                              return cur_dates;
                            }).flatten();
      
      dates = dates.slice(6, dates.length().subtract(6));
      
      //moving yearly average
      var im_selected_year_mean = dates.map(function(date){
                                                    var start = ee.Date(date).advance(-6, "month");
                                                    var end = ee.Date(date).advance(6, "month");
                                                    var reducer = ee.Reducer.median().combine({reducer2:ee.Reducer.count(), sharedInputs:true});
                                                    return ee.Image(landsat.filterDate(start, end)
                                                                            .select("bsi")
                                                                            .reduce(reducer))
                                                                    .rename(["bsi", "bsi_count"])
                                                    .set("system:time_start", start.millis());
                                                  });
      
      im_selected_year_mean = ee.ImageCollection(im_selected_year_mean)
                                                      .map(function(im){
                                                                          return im.updateMask(im.select("bsi_count").gte(3))
                                                                                   .select("bsi");
                                                       });
                              
      var bsi_thres = -0.2031250;
      
      var im_disturb_yearly = ee.ImageCollection(im_selected_year_mean)
                                      .map(function(im){return  im.select("bsi")
                                                                  .updateMask(im.select("bsi").gte(bsi_thres))
                                                                  .addBands(ee.Image(im.date().difference(ee.Date('1980-01-01'), 'year')).float()
                                                                                            .updateMask(im.select("bsi").gte(bsi_thres))
                                                                  .rename("constant"));
                                            })
                                      .qualityMosaic("constant").select(["constant"]) ;
          
      var tableOPregions = ee.FeatureCollection("users/OlhaDanylo/global/tropics/GAUL_OPregions_validation")
                                          .map(function(ft){return ft.buffer(3000)});
      
      var reg_list = ee.Dictionary(tableOPregions.aggregate_histogram("OP_reg")).keys();
      
      var OPregions_image = tableOPregions.remap(reg_list,
                                                  ee.List.sequence(1,7), "OP_reg")
                                                  .reduceToImage(["OP_reg"], ee.Reducer.first()).rename("op_regions");
      
      var regions  = grid.filterBounds(ee.FeatureCollection([kalimantan, sulawesi, sumatra, thailand]));
      
      for (var i = 0; i<ee.List(regions.aggregate_array("FID")).length().getInfo(); ++i){
      
              var roi = regions.filterMetadata("FID", "equals", ee.List(regions.aggregate_array("FID")).get(i));
              var input = im_disturb_yearly.mask(OPregions_image.mask()).multiply(100).round().int16().clip(roi);
              Export.image.toAsset({image:input, 
                    assetId:ee.String("users/OlhaDanylo/palmoil/danylo_po_age/").cat(ee.String(ee.Number(ee.List(regions.aggregate_array("FID")).get(i)))).getInfo(), 
                    region:roi, 
                    maxPixels:3608235476110, 
                    scale:30, 
                    description: ee.String("po_age_").cat(ee.String(ee.Number(i))).getInfo()});
      }                            
};
